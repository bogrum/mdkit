#!/usr/bin/env python
"""mdkit cizim katmani.

Cizim ANALIZ ADINA gore degil, veri SEKLINE gore dallanir: timeseries ve
profile. Yeni bir timeseries analizi (gyrate, sasa) eklendiginde bu dosyaya
dokunmak gerekmez.

Birim donusumu burada yapilir: .xvg'ler ps ve nm cinsindendir (spec 2.5).
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

MDKIT = Path(__file__).resolve().parent

NM_TO_ANGSTROM = 10.0
PS_TO_NS = 1e-3

# dataviz skill'inin dogrulanmis paleti (categorical, slot 1/2/3 ve 1/8).
# Sadece bu sabitler degistirilir; kodun geri kalani ayni kalir.
REP_COLORS = {"rep1": "#2a78d6", "rep2": "#eb6834", "rep3": "#1baf7a"}
# Kompleks gruplarinin renkleri SIRAYLA atanir: hangi grubun var oldugu
# config.sh'deki COMPLEX_GROUPS'a baglidir, bu dosyaya degil.
GROUP_PALETTE = ["#2a78d6", "#e34948"]
# Gruplar tanimli ama kompleks hicbirine uymuyorsa: notr gri. Bir kategori
# rengi vermek onu yanlis gruba aitmis gibi gosterirdi.
UNGROUPED_COLOR = "#8a8f98"
# Seriler cizgi TIPIYLE ayrilir; renk replikaya ayrilmis durumda.
SERIES_LINESTYLES = ["-", "--", ":", "-."]

DEFAULT_COMPARE_OUTPUT = "rmsd_pep_on_mhc.xvg"

# Isi haritasi icin algisal olarak duzgun, TEK YONLU bir skala. Kategorik
# REP_COLORS paleti burada kullanilmaz: o palet ust uste binen replika
# CIZGILERINI ayirmak icin secilmisti, sirali bir buyuklugu kodlamak icin degil.
MATRIX_CMAP = "viridis"

# Profil konum kutulari. Peptidler 8-11 residue arasinda degistigi icin ham
# residue numarasi kompleksler arasi karsilastirilamaz: bir 8-mer'in 5.
# residue'su ortada, 11-mer'in 5.'si degil. Kutular iki UCTAN sayilarak
# tanimlanir, boylece ayni etiket ayni ROLE denk gelir.
POSITION_BINS = ["P1", "P2", "orta", "PO-1", "PO"]

# Denge egrisinin yuvarlanan ortalama penceresi, ns. SABIT bir SUREDIR;
# serinin oranı (n/20) DEGILDIR. Gerekce: yuvarlanan ortalama bir alcak
# geciren filtredir ve kesme frekansi ~1/pencere; dolayisiyla pencere,
# ayrilmak istenen iki FIZIKSEL zaman olceginin arasina konur. Oransal bir
# kural filtreyi kosu uzunluguna baglar -- ayni sistem 20 ns yerine 100 ns
# kosuldugunda ayni surec farkli duzlestirilir.
#
# 1.0 ns secildi cunku gercek veride olculen otokorelasyon sureleri
# (tau_int) 0.5-7.8 ns araliginda, medyan ~3.9 ns; 1 ns bunlarin hemen
# hepsinin ALTINDA kalir, yani korelasyonlu gercek yapiyi bozmadan hizli
# titresimi alir. Onceki n/20 kurali 100 ns'lik kosuda 4.6 ns veriyor ve
# bir konformasyonel cikisin genliginin %69'unu yutuyordu.
DEFAULT_MATRIX_SMOOTH_NS = 1.0


def unit_label(unit):
    """xmgrace bicim kodlarini matplotlib mathtext'e cevirir.

    gmx birimleri xmgrace bicimiyle yazar: `nm\\S2\\N` = nm^2
    (`\\S` ust simge baslat, `\\s` alt simge, `\\N` normale don). Ham
    hali eksen etiketinde `nm\\S2\\N` olarak gorunur.

    Bu YALNIZCA GOSTERIM icindir. Birim tanima ve cevrim bundan
    etkilenmez -- nm^2'yi nm gibi 10 ile carpmak sessizce YANLIS olurdu.
    """
    if not unit:
        return unit
    s = re.sub(r"\\S(.*?)\\N", r"$^{\1}$", unit)
    s = re.sub(r"\\s(.*?)\\N", r"$_{\1}$", s)
    return re.sub(r"\\[A-Za-z]", "", s)


def scale_and_label(unit):
    """(carpan, eksen etiketi). Yalnizca nm cevrilir; bilinmeyen birim
    cevrilmeden, kendi etiketiyle cizilir -- bilmedigimiz bir birimi
    cevirmis gibi yapmiyoruz. Etiket xmgrace bicim kodlarindan arindirilir
    ama bu KARSILASTIRMAYI etkilemez: unit ham haliyle kiyaslanir."""
    if unit == "nm":
        return NM_TO_ANGSTROM, "Å"
    return 1.0, unit_label(unit) or "birimsiz"


def _warn_unknown_unit(output, unit):
    print(f"uyari: {output} icin bilinmeyen birim '{unit}' - "
          "donusum uygulanmadi, oldugu gibi cizildi", file=sys.stderr)


def load(results_dir, required=True):
    ts_path = results_dir / "timeseries_long.csv"
    pr_path = results_dir / "profile_long.csv"
    if not ts_path.exists() and not pr_path.exists():
        # --matrix tek basina istendiginde CSV'ler olmayabilir: matrisler
        # ayri bir urun ve cross_rmsd tek basina kosulmus olabilir.
        if required:
            sys.exit(f"sonuc CSV'leri bulunamadi: {results_dir}")
        return pd.DataFrame(), pd.DataFrame()
    # unit str okunur: birimsiz ciktilarda (gmx hbond: y ekseni "Hbonds",
    # parantezli birim yok) collect bos yazar ve pandas bunu NaN yapar --
    # NaN etiketi re.sub'i patlatiyor, onceki surumde eksene "nan" yaziyordu.
    kw = dict(converters={"unit": str})
    ts = pd.read_csv(ts_path, **kw) if ts_path.exists() else pd.DataFrame()
    pr = pd.read_csv(pr_path, **kw) if pr_path.exists() else pd.DataFrame()
    return ts, pr


# --- kompleks gruplari: projeye ozel, config.sh'den okunur ---------------

def parse_complex_groups(raw):
    """'onek:etiket' ogelerini [(onek, etiket)] olarak cozer.

    'top'/'last' ayrimi BU PROJENIN kurgusudur, aracin degil: sabit kodlu
    oldugu surece baska bir veri setinde her kompleks 'last' rengini alir
    ve efsane 'top*/last*' yazar."""
    groups = []
    for item in raw.split():
        prefix, sep, label = item.partition(":")
        if not sep or not prefix:
            print(f"uyari: COMPLEX_GROUPS ogesi 'onek:etiket' bicimine "
                  f"uymuyor, atlandi: {item!r}", file=sys.stderr)
            continue
        groups.append((prefix, label or prefix))
    return groups


def read_complex_groups(config=None):
    """config.sh'deki COMPLEX_GROUPS; okunamazsa gruplama YOK.

    Config'i bash tarafindaki dogrulayiciyla okur (collect_results.py ile
    ayni tek kaynak). Sonuc CSV'leri baska bir makineye tasinip config'siz
    cizdirilebildigi icin basarisizlik olumcul degildir."""
    arg = f'"{config}"' if config else ""
    r = subprocess.run(
        ["bash", "-c",
         f'source "{MDKIT}/analysis/lib.sh" && mdkit_load_config {arg} '
         f'>/dev/null 2>&1 && printf "%s" "${{COMPLEX_GROUPS[*]:-}}"'],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        return []
    return parse_complex_groups(r.stdout)


def group_color(complex_name, groups):
    """Kompleksin rengi. groups bossa tek renk: gruplama yapilandirilmamis."""
    if not groups:
        return GROUP_PALETTE[0]
    for i, (prefix, _label) in enumerate(groups):
        if complex_name.startswith(prefix):
            return GROUP_PALETTE[i % len(GROUP_PALETTE)]
    return UNGROUPED_COLOR


def series_note(g):
    """Tek serili ciktilarda seri adi; cok serilide bos.

    gmx rms, NEYIN NEYE FIT edildigini .xvg'nin subtitle'ina yazar
    ("LIGAND after lsq fit to RECEPTOR_BB") ve collect bunu `series`
    kolonuna alir. Ama tek seri varken cizgi etiketi replika adidir
    (`rep1`) ve seri adi DUSERDI -- bilgi veride var, figurde yok.
    Baslikta bir kez gosterilir.

    Cok serili ciktilarda (or. sasa: Total + LIGAND) legend zaten
    gosterir, baslikta tekrarlanmaz.
    """
    s = g["series"].dropna().unique() if "series" in g else []
    return str(s[0]) if len(s) == 1 else ""


def _panels(ts, pr):
    """(df, x kolonu, x etiketi, x carpani) uclusu uretir."""
    if not ts.empty:
        yield ts, "time_ps", "Zaman (ns)", PS_TO_NS
    if not pr.empty:
        yield pr, "residue", "Residue", 1.0


# --- seri ayrimi ---------------------------------------------------------

def line_specs(g, xcol):
    """Bir (kompleks, cikti) panelinde cizilecek cizgiler.

    Her (replika, SERI) icin bir cizgi: renk replikayi, cizgi tipi seriyi
    kodlar. Seriyi yok sayip yalnizca replikaya gore gruplamak cok serili
    bir ciktida (gmx gyrate: Rg/RgX/RgY/RgZ) tum serileri tek bir zikzak
    cizgide birlestirirdi. Tek serili ciktida davranis degismez: duz cizgi,
    etiket = replika adi."""
    series_names = sorted(g["series"].astype(str).unique())
    style = {s: SERIES_LINESTYLES[i % len(SERIES_LINESTYLES)]
             for i, s in enumerate(series_names)}
    specs = []
    for (rep, series), gr in g.groupby(["replica", "series"], sort=True):
        gr = gr.sort_values(xcol)
        specs.append({
            "replica": str(rep),
            "series": str(series),
            "label": str(rep) if len(series_names) == 1
                     else f"{rep} · {series}",
            "color": REP_COLORS.get(rep),
            "linestyle": style[str(series)],
            "x": gr[xcol].to_numpy(),
            "y": gr["value"].to_numpy(),
        })
    return specs


def series_stats(g, xcol):
    """(seri, mean/std tablosu) ciftleri -- SD her serinin ICINDE hesaplanir.

    Serileri havuzlamak, farkli fiziksel buyukluklerin (Rg ile RgX) arasindaki
    yayilimi replika degiskenligi gibi gosteren, yayin gorunumlu ama anlamsiz
    bir band uretir."""
    for series, gs in g.groupby("series", sort=True):
        stats = (gs.groupby(xcol)["value"]
                   .agg(["mean", "std"])
                   .sort_index()
                   .fillna(0.0))
        yield str(series), stats


def plot_per_complex(ts, pr, out_dir, groups=()):
    """Kompleks basina bir figure; rep1/rep2/rep3 ust uste. Yakinsama denetimi."""
    target = out_dir / "per_complex"
    target.mkdir(parents=True, exist_ok=True)
    for df, xcol, xlabel, xconv in _panels(ts, pr):
        for (cx, output), g in df.groupby(["complex", "output"], sort=True):
            fig = None
            try:
                unit = g["unit"].iloc[0]
                yconv, ylabel = scale_and_label(unit)
                if unit != "nm":
                    _warn_unknown_unit(output, unit)
                fig, ax = plt.subplots(figsize=(7, 4))
                for spec in line_specs(g, xcol):
                    ax.plot(spec["x"] * xconv, spec["y"] * yconv,
                            label=spec["label"], color=spec["color"],
                            linestyle=spec["linestyle"], linewidth=1.0)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                not_ = series_note(g)
                ax.set_title(f"{cx} — {Path(output).stem}"
                             + (f"\n{not_}" if not_ else ""), fontsize=10.5)
                ax.legend(frameon=False)
                ax.spines[["top", "right"]].set_visible(False)
                fig.tight_layout()
                fig.savefig(target / f"{cx}_{Path(output).stem}.png", dpi=150)
            except Exception as exc:
                print(f"per-complex cizimi basarisiz ({cx}, {output}): {exc}",
                      file=sys.stderr)
            finally:
                if fig is not None:
                    plt.close(fig)


def plot_mean_sd(ts, pr, out_dir, groups=()):
    """Replika ortalamasi + ±SD seridi. Yayina/teze giden temiz figure."""
    target = out_dir / "mean_sd"
    target.mkdir(parents=True, exist_ok=True)
    for df, xcol, xlabel, xconv in _panels(ts, pr):
        for (cx, output), g in df.groupby(["complex", "output"], sort=True):
            fig = None
            try:
                unit = g["unit"].iloc[0]
                yconv, ylabel = scale_and_label(unit)
                if unit != "nm":
                    _warn_unknown_unit(output, unit)
                color = group_color(cx, groups)

                fig, ax = plt.subplots(figsize=(7, 4))
                bands = list(series_stats(g, xcol))
                for i, (series, stats) in enumerate(bands):
                    x = stats.index.to_numpy() * xconv
                    mean = stats["mean"].to_numpy() * yconv
                    sd = stats["std"].to_numpy() * yconv
                    ls = SERIES_LINESTYLES[i % len(SERIES_LINESTYLES)]
                    ax.fill_between(x, mean - sd, mean + sd, color=color,
                                    alpha=0.25, linewidth=0)
                    ax.plot(x, mean, color=color, linewidth=1.4, linestyle=ls,
                            label=series)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                not_ = series_note(g)
                ax.set_title(
                    f"{cx} — {Path(output).stem} "
                    f"(n={g['replica'].nunique()} replika)"
                    + (f"\n{not_}" if not_ else ""), fontsize=10.5)
                if len(bands) > 1:
                    ax.legend(frameon=False)
                ax.spines[["top", "right"]].set_visible(False)
                fig.tight_layout()
                fig.savefig(target / f"{cx}_{Path(output).stem}.png", dpi=150)
            except Exception as exc:
                print(f"mean-sd cizimi basarisiz ({cx}, {output}): {exc}",
                      file=sys.stderr)
            finally:
                if fig is not None:
                    plt.close(fig)


def compare_order(per_rep):
    """Kompleksleri medyan degere gore artan sirada dondurur (saf fonksiyon)."""
    return (per_rep.groupby("complex")["value"]
                   .median()
                   .sort_values()
                   .index
                   .tolist())


def plot_compare(ts, out_dir, output=DEFAULT_COMPARE_OUTPUT, groups=()):
    """Tum kompleksler tek panelde, medyana gore sirali boxplot."""
    if ts.empty:
        return
    df = ts[ts["output"] == output]
    if df.empty:
        return
    unit = df["unit"].iloc[0]
    yconv, ylabel = scale_and_label(unit)
    if unit != "nm":
        _warn_unknown_unit(output, unit)
    per_rep = (df.groupby(["complex", "replica"])["value"].mean()
                 .mul(yconv)
                 .reset_index())
    order = compare_order(per_rep)
    data = [per_rep.loc[per_rep["complex"] == c, "value"].to_numpy() for c in order]

    fig, ax = plt.subplots(figsize=(max(8.0, len(order) * 0.35), 4.5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.6)
    for patch, cx in zip(bp["boxes"], order):
        patch.set_facecolor(group_color(cx, groups))
        patch.set_alpha(0.75)
        patch.set_edgecolor("#333333")
    for median in bp["medians"]:
        median.set_color("#222222")

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel(f"Ortalama {Path(output).stem} ({ylabel})")
    ax.set_title("Kompleksler arasi karsilastirma (replika basina ortalama)")
    ax.spines[["top", "right"]].set_visible(False)
    # Gruplama yapilandirilmamissa efsane de YOK: uydurma bir 'top*/last*'
    # efsanesi baska bir veri setinde duz yalan olurdu.
    if groups:
        handles = [
            plt.Line2D([], [], linewidth=6, alpha=0.75,
                       color=GROUP_PALETTE[i % len(GROUP_PALETTE)])
            for i in range(len(groups))
        ]
        ax.legend(handles, [label for _p, label in groups],
                  frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / f"compare_{Path(output).stem}.png", dpi=150)
    plt.close(fig)


def load_matrices(results_dir):
    """results/matrices/*.npz -> {(kompleks, analiz): [kayit]}.

    Bozuk tek bir dosya digerlerini dusurmez; adiyla stderr'e yazilir."""
    mdir = results_dir / "matrices"
    groups = {}
    if not mdir.is_dir():
        return groups
    for f in sorted(mdir.glob("*.npz")):
        try:
            with np.load(f, allow_pickle=False) as d:
                rec = {
                    "values": d["values"],
                    "x_ps": d["x_ps"], "y_ps": d["y_ps"],
                    "unit": str(d["unit"]),
                    "complex": str(d["complex"]),
                    "replica_i": str(d["replica_i"]),
                    "replica_j": str(d["replica_j"]),
                    "analysis": str(d["analysis"]),
                    # Bu alan sonradan eklendi; onceden yazilmis .npz'ler
                    # onu icermez ve yeniden cizim patlamamali.
                    "title": str(d["title"]) if "title" in d.files else "",
                    "subtitle": (str(d["subtitle"])
                                 if "subtitle" in d.files else ""),
                }
        except Exception as exc:
            print(f"matris okunamadi ({f.name}): {exc}", file=sys.stderr)
            continue
        groups.setdefault((rec["complex"], rec["analysis"]), []).append(rec)
    return groups


def equilibrium_curves(recs):
    """rep_i -> (x_ps, egri). Egri: rep_i'nin her frame'inin TUM eslerin
    TUM frame'lerine ortalama uzakligi.

    Matris duzeni [y=es frame, x=kendi frame] oldugu icin esler axis=0'da
    yigilir ve sonuc kendi frame sayisi kadar uzunluktadir. Bu, uc replikayi
    uc uca ekleyip 3Nx3N matriste satir ortalamasi almanin karsiligidir --
    ama replikalar burada ORTAK bir zaman eksenine oturur, uc uca eklenmez,
    yani dogrudan karsilastirilabilirler.

    Egri duzlestiginde replika yeni konformasyonel bolge bulmayi birakmis
    demektir; denge argumani icin isi haritasindan daha okunakli.
    """
    by_i = {}
    for r in recs:
        by_i.setdefault(r["replica_i"], []).append(r)

    curves = {}
    for rep_i, rs in by_i.items():
        widths = {int(r["values"].shape[1]) for r in rs}
        if len(widths) != 1:
            # Yigmak sessizce yanlis olurdu: ayni replikanin matrisleri ayni
            # kendi-frame sayisina sahip OLMALI.
            print(f"denge egrisi atlandi ({rep_i}): kendi frame sayisi "
                  f"matrisler arasinda tutmuyor {sorted(widths)}",
                  file=sys.stderr)
            continue
        stack = np.concatenate([r["values"] for r in rs], axis=0)
        curves[rep_i] = (rs[0]["x_ps"], stack.mean(axis=0))
    return curves


def matrix_caption(rec):
    """Figurde gosterilecek aciklama: NE olculdu, NEYE fit edildi.

    subtitle tercih edilir cunku ikisini birden tasir
    ("LIGAND_BB after lsq fit to RECEPTOR_BB"); title yalnizca olculeni
    soyler. Ikisi de yoksa bos -- uydurulmaz."""
    return rec.get("subtitle") or rec.get("title") or ""


def equilibrium_window(x_ps, smooth_ns):
    """Verilen SURE penceresinin kac frame ettigini dondurur; 0 = duzlestirme yok.

    Kayit araligi verinin kendisinden okunur: 50 frame tek basina bir sey
    ifade etmez, 200 ps araliktaki 5 frame ile 10 ps araliktaki 5 frame
    farkli surelerdir. 3 frame'in altinda duzlestirme anlamsizdir."""
    if smooth_ns <= 0 or len(x_ps) < 2:
        return 0
    dt_ns = float(x_ps[1] - x_ps[0]) * PS_TO_NS
    if dt_ns <= 0:
        return 0
    w = int(round(smooth_ns / dt_ns))
    return w if 3 <= w <= len(x_ps) else 0


def window_label(window, x_ps):
    """Pencereyi hem frame hem SURE olarak yazar.

    Yalnizca frame sayisi yazmak yetmez: ayni 5 frame, -dt 200 ile 1 ns,
    -dt 1000 ile 5 ns eder. dt cikarilamiyorsa uydurulmaz."""
    if len(x_ps) > 1:
        dt_ns = float(x_ps[1] - x_ps[0]) * PS_TO_NS
        return f"pencere={window} frame ({window * dt_ns:g} ns)"
    return f"pencere={window} frame"


def rolling_mean(y, window):
    """(baslangic_indeksi, ortalamalar). Pencere ORTALANIR.

    Uc noktalari pencerenin SONUNA baglamak (x[window-1:]) egriyi yarim
    pencere kadar saga kaydirir ve "ne zaman dengelendi" sorusunu sistematik
    olarak GEC cevaplar. mode='same' ise kenarlari sifirla doldurup uclari
    asagi ceker; bu yuzden 'valid' + ortalanmis indeks.
    """
    if window < 1 or len(y) < window:
        return 0, np.asarray([], dtype=float)
    vals = np.convolve(y, np.ones(window) / window, mode="valid")
    return (window - 1) // 2, vals


def global_spans(groups):
    """(analiz, birim) -> (vmin, vmax), TUM kompleksler uzerinden.

    Birim anahtarin PARCASI: nm ile nm^2 matrisleri ayni renk skalasina
    konamaz, cevrim carpanlari farklidir."""
    spans = {}
    for (_cx, analysis), recs in groups.items():
        unit = recs[0]["unit"]
        conv, _ulabel = scale_and_label(unit)
        lo = min(float(r["values"].min()) for r in recs) * conv
        hi = max(float(r["values"].max()) for r in recs) * conv
        key = (analysis, unit)
        if key in spans:
            spans[key] = (min(spans[key][0], lo), max(spans[key][1], hi))
        else:
            spans[key] = (lo, hi)
    return spans


def _draw_matrix_grid(cx, analysis, recs, target, conv, ulabel, vlim=None):
    """Kompleks basina N x N isi haritasi izgarasi, ORTAK renk skalasiyla.

    Ortak skala sart: panel basina ayri skala, farkli replika ciftlerini
    gorsel olarak karsilastirilamaz kilardi -- bu figurun tek amaci o
    karsilastirma."""
    rows = sorted({r["replica_j"] for r in recs})
    cols = sorted({r["replica_i"] for r in recs})
    by_cell = {(r["replica_i"], r["replica_j"]): r for r in recs}

    if vlim is None:
        vmin = min(float(r["values"].min()) for r in recs) * conv
        vmax = max(float(r["values"].max()) for r in recs) * conv
    else:
        vmin, vmax = vlim

    fig, axes = plt.subplots(
        len(rows), len(cols), squeeze=False,
        figsize=(2.6 * len(cols) + 1.6, 2.6 * len(rows)),
        sharex=True, sharey=True,
    )
    try:
        im = None
        for ri, rep_j in enumerate(rows):
            for ci, rep_i in enumerate(cols):
                ax = axes[ri][ci]
                rec = by_cell.get((rep_i, rep_j))
                if rec is None:
                    ax.set_axis_off()
                    continue
                x, y = rec["x_ps"], rec["y_ps"]
                im = ax.imshow(
                    rec["values"] * conv, origin="lower", aspect="auto",
                    vmin=vmin, vmax=vmax, cmap=MATRIX_CMAP,
                    extent=[x[0] * PS_TO_NS, x[-1] * PS_TO_NS,
                            y[0] * PS_TO_NS, y[-1] * PS_TO_NS],
                )
                if ri == len(rows) - 1:
                    ax.set_xlabel(f"{rep_i} (ns)")
                if ci == 0:
                    ax.set_ylabel(f"{rep_j} (ns)")
        if im is not None:
            fig.colorbar(im, ax=axes, label=f"RMSD ({ulabel})",
                         fraction=0.046, pad=0.02)
        # gmx'in kendi basligi: NEYIN olculdugu. Olmadan okuyucu grafikten
        # hangi buyuklugun cizildigini anlayamaz.
        caption = matrix_caption(recs[0])
        # Araligi HER IKI sette de basliga yaz: kompleks basina skalada
        # "bu yesil kac angstrom?" sorusunun cevabi figurde olmazsa okuyucu
        # renkleri kompleksler arasi kiyaslamaya kalkar -- ve yanilir.
        kind = "ortak renk skalasi" if vlim is not None else "renk skalasi"
        scale_note = f"  [{kind}: {vmin:.1f}-{vmax:.1f} {ulabel}]"
        fig.suptitle(f"{cx} \u2014 {analysis}{scale_note}"
                     + (f"\n{caption}" if caption else ""))
        fig.savefig(target / f"{cx}_{analysis}.png", dpi=150,
                    bbox_inches="tight")
    finally:
        plt.close(fig)


def _draw_equilibrium(cx, analysis, recs, target, conv, ulabel,
                      smooth_ns=DEFAULT_MATRIX_SMOOTH_NS):
    """Replika basina denge egrisi, ortak zaman ekseninde."""
    curves = equilibrium_curves(recs)
    if not curves:
        return
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    try:
        for rep_i, (x_ps, curve) in sorted(curves.items()):
            color = REP_COLORS.get(rep_i, UNGROUPED_COLOR)
            x_ns = x_ps * PS_TO_NS
            y = curve * conv
            ax.plot(x_ns, y, color=color, alpha=0.35, linewidth=1)
            # Pencere frame SAYISIYLA degil, serinin oranıyla secilir:
            # sabit bir pencere farkli -dt degerlerinde farkli sureye denk
            # gelir ve egriler kiyaslanamaz hale gelir.
            # Pencere AYNI figurdeki butun egrilerde ayni: her egriye kendi
            # otokorelasyon suresinden pencere vermek onlari
            # karsilastirilamaz kilardi.
            window = equilibrium_window(x_ps, smooth_ns)
            start, roll = rolling_mean(y, window)
            if roll.size:
                ax.plot(x_ns[start:start + roll.size], roll, color=color,
                        linewidth=2,
                        label=f"{rep_i} \u00b7 {window_label(window, x_ps)}")
            else:
                ax.plot([], [], color=color, linewidth=2, label=rep_i)
        ax.set_xlabel("Zaman (ns)")
        ax.set_ylabel(f"Ortalama RMSD ({ulabel})")
        caption = matrix_caption(recs[0])
        ax.set_title(f"{cx} \u2014 {analysis}: denge analizi"
                     + (f"\n{caption}" if caption else "")
                     + "\n(her frame'in tum es frame'lere ortalama uzakligi)")
        ax.legend()
        fig.savefig(target / f"{cx}_{analysis}_equilibrium.png", dpi=150,
                    bbox_inches="tight")
    finally:
        plt.close(fig)


def position_bin(residue, residue_from_end):
    """Residue -> konum kutusu. N-ucu ONCELIKLI.

    Oncelik sart: 3 residue'luk bir peptidde 2. residue hem `P2` hem
    `PO-1` olurdu ve ayni olcum iki kutuya birden girerdi.
    """
    if residue == 1:
        return "P1"
    if residue == 2:
        return "P2"
    if residue_from_end == 0:
        return "PO"
    if residue_from_end == -1:
        return "PO-1"
    return "orta"


def group_of(complex_name, groups):
    """Kompleksin grup etiketi; hicbir onege uymuyorsa None."""
    for prefix, label in groups:
        if complex_name.startswith(prefix):
            return label
    return None


def profile_length_summary(g):
    """Grup basina profil uzunlugu: [(etiket, medyan, min, max, n_kompleks)].

    Uzunluk kolonlardan turetilir: `residue_from_end = residue - n_res`
    oldugu icin `n_res = residue - residue_from_end`. Ayri bir groupby
    gerekmez.
    """
    t = g.assign(_n=g["residue"] - g["residue_from_end"])
    per_cx = t.groupby(["grup", "complex"])["_n"].first()
    out = []
    for lbl, s in per_cx.groupby(level=0):
        out.append([lbl, float(s.median()), int(s.min()), int(s.max()),
                    int(s.size)])
    return out


def profile_length_note(g):
    """Figure yazilacak uzunluk dagilimi satiri.

    HER ZAMAN gosterilir, bir esige BAGLANMAZ. Konuma dayali
    karsilastirmalarda uzunluk en sik karisan degiskendir ve uclardan
    sayilan konumlar (PO, PO-1) ozellikle duyarlidir: "son residue" tanim
    geregi sistemin BITTIGI yerdir, ayni residue degil. Bir esik koymak
    ("medyanlar X'ten fazla farkliysa uyar") tam esikte olan bir vakanin
    sessizce gecmesi demek olurdu; karar okuyucunun.
    """
    parcalar = [f"{lbl}: medyan {med:g} ({lo}-{hi}, n={n})"
                for lbl, med, lo, hi, n in profile_length_summary(g)]
    return "profil uzunlugu -- " + "  |  ".join(parcalar)


def profile_complex_means(g, conv=1.0):
    """(kutu, grup, kompleks) basina tek deger -- replikalar ORTALANIR.

    Sart: aksi halde ayni kompleksin uc replikasi uc BAGIMSIZ gozlem
    sayilir, orneklem yapay olarak ucer katina cikar ve p degeri
    oldugundan kucuk cikar. Replikalar ayni sistemin tekrarlaridir,
    bagimsiz ornekler degil."""
    return (g.groupby(["kutu", "grup", "complex"])["value"]
            .mean().mul(conv).reset_index())


def plot_profile_compare(pr, out_dir, groups=()):
    """Profil ciktilarini KONUMA gore gruplar arasi karsilastirir.

    Her `profile` ciktisi icin bir figur: konum kutulari x ekseninde, her
    kutuda grup basina bir boxplot, ustunde testin p degeri.

    Gruplama `config.sh`'deki COMPLEX_GROUPS'tan gelir; projeye ozel bilgi
    (hangi onek hangi etiket) koda GIRMEZ. Grup tanimli degilse ya da
    veride tek grup varsa mod sessizce gecilir -- karsilastirilacak sey
    yoksa uydurma p uretmek yanlis olur.

    DIKKAT: `rmsf_mhc` gibi uzun profillerde "orta" kutusu yuzlerce residue
    icerir; figur teknik olarak dogru ama az bilgilendiricidir. Sihirli bir
    uzunluk esigiyle gizlenmiyor: sessizce kaybolan cikti, zayif ciktidan
    kotudur.
    """
    if pr is None or len(pr) == 0 or not groups:
        return
    if "residue_from_end" not in pr.columns:
        # Bu kolon sonradan eklendi. Eski bir profile_long.csv ile cizim
        # yapiliyorsa mod atlanir -- ama SESSIZCE degil: kullanici neden
        # figur uretilmedigini bilmeli ve cozumu tek komut.
        print("uyari: profile_long.csv'de 'residue_from_end' kolonu yok "
              "(eski sema); --profile-compare atlandi. "
              "collect_results.py'yi yeniden calistirin.", file=sys.stderr)
        return
    df = pr.copy()
    df["grup"] = df["complex"].map(lambda c: group_of(c, groups))
    df = df[df["grup"].notna()]
    if df.empty or df["grup"].nunique() < 2:
        return
    df["kutu"] = [position_bin(r, e) for r, e
                  in zip(df["residue"], df["residue_from_end"])]

    target = out_dir / "profile_compare"
    target.mkdir(parents=True, exist_ok=True)
    etiketler = [lbl for _p, lbl in groups]

    for output, g in df.groupby("output", sort=True):
        fig = None
        try:
            unit = g["unit"].iloc[0]
            conv, ylabel = scale_and_label(unit)
            if unit != "nm":
                _warn_unknown_unit(output, unit)
            per_cx = profile_complex_means(g, conv)
            mevcut = [k for k in POSITION_BINS
                      if k in set(per_cx["kutu"])]
            gruplar = [lbl for lbl in etiketler
                       if lbl in set(per_cx["grup"])]
            n = len(gruplar)
            w = 0.8 / n

            fig, ax = plt.subplots(figsize=(max(7.0, 1.9 * len(mevcut)), 4.8))
            ust = per_cx["value"].max()
            for i, kutu in enumerate(mevcut):
                orn = []
                for j, lbl in enumerate(gruplar):
                    v = per_cx[(per_cx.kutu == kutu)
                               & (per_cx.grup == lbl)]["value"].to_numpy()
                    orn.append(v)
                    if v.size == 0:
                        continue
                    bp = ax.boxplot(
                        [v], positions=[i + (j - (n - 1) / 2) * w],
                        widths=w * 0.85, patch_artist=True,
                        medianprops=dict(color="black", lw=1.3),
                        flierprops=dict(ms=3, mfc="0.5", mec="0.5"))
                    bp["boxes"][0].set_facecolor(
                        GROUP_PALETTE[j % len(GROUP_PALETTE)])
                    bp["boxes"][0].set_alpha(0.75)
                gecerli = [v for v in orn if v.size > 0]
                if len(gecerli) < 2:
                    continue
                if len(gecerli) == 2:
                    p = stats.mannwhitneyu(*gecerli,
                                           alternative="two-sided").pvalue
                else:
                    p = stats.kruskal(*gecerli).pvalue
                ax.text(i, ust * 1.08, f"p={p:.3f}", ha="center", fontsize=9,
                        color="black" if p < 0.05 else "0.45",
                        fontweight="bold" if p < 0.05 else "normal")

            ax.set_xticks(range(len(mevcut)))
            ax.set_xticklabels(mevcut)
            ax.set_xlabel("Konum (zincirin uclarindan sayilarak)")
            ax.set_ylabel(ylabel)
            ax.set_ylim(top=ust * 1.16)
            ax.spines[["top", "right"]].set_visible(False)
            test_adi = "Mann-Whitney" if n == 2 else "Kruskal-Wallis"
            # Coklu test uyarisi BASLIKTA: figuru tek basina goren biri
            # p<0.05'i anlamli sanmasin.
            ax.set_title(
                f"{Path(output).stem} \u2014 konuma gore gruplar\n"
                f"{test_adi}; {len(mevcut)} konum test edildi, "
                f"Bonferroni esigi 0.05/{len(mevcut)} = "
                f"{0.05 / len(mevcut):.3f}\n"
                f"{profile_length_note(g)}", fontsize=9.5)
            ax.legend(handles=[
                mpatches.Patch(fc=GROUP_PALETTE[j % len(GROUP_PALETTE)],
                               alpha=0.75, label=lbl)
                for j, lbl in enumerate(gruplar)],
                frameon=False, fontsize=9)
            fig.tight_layout()
            fig.savefig(target / f"{Path(output).stem}.png", dpi=150)
        except Exception as exc:
            print(f"profil karsilastirma cizimi basarisiz ({output}): {exc}",
                  file=sys.stderr)
        finally:
            if fig is not None:
                plt.close(fig)


def plot_matrix(results_dir, out_dir, smooth_ns=DEFAULT_MATRIX_SMOOTH_NS):
    """--matrix modu: kompleks basina isi haritasi izgarasi + denge egrisi."""
    groups = load_matrices(results_dir)
    if not groups:
        return
    target = out_dir / "matrix"
    target.mkdir(parents=True, exist_ok=True)

    # Kompleks basina skala, kompleks ICI kontrasti korur ama kompleksler
    # ARASI renk okumayi imkansiz kilar (olculdu: top1'de "yesil" 2.9 A,
    # last10'da 8.9 A). Ikisi de gerekli, o yuzden ikinci bir set ortak
    # skalayla yazilir. Tek kompleksli bir (analiz, birim) grubunda ortak
    # skala kendi skalasiyla ayni oldugundan yazilmaz.
    spans = global_spans(groups)
    per_key = {}
    for (_cx, analysis), recs in groups.items():
        per_key[(analysis, recs[0]["unit"])] = \
            per_key.get((analysis, recs[0]["unit"]), 0) + 1
    global_target = out_dir / "matrix_global"

    for (cx, analysis), recs in sorted(groups.items()):
        # Birim kurali .xvg ile birebir ayni: yalnizca nm cevrilir,
        # taninmayan birim cevrilmeden kendi etiketiyle cizilir.
        unit = recs[0]["unit"]
        conv, ulabel = scale_and_label(unit)
        if unit != "nm":
            _warn_unknown_unit(f"{cx} ({analysis})", unit)

        # Izgara ve denge egrisi AYRI yalitilir: biri patlarsa digeri yine
        # de diske yazilmis olur.
        try:
            _draw_matrix_grid(cx, analysis, recs, target, conv, ulabel)
        except Exception as exc:
            print(f"matris cizimi basarisiz ({cx}, {analysis}): {exc}",
                  file=sys.stderr)

        key = (analysis, unit)
        if per_key.get(key, 0) > 1:
            try:
                global_target.mkdir(parents=True, exist_ok=True)
                _draw_matrix_grid(cx, analysis, recs, global_target, conv,
                                  ulabel, vlim=spans[key])
            except Exception as exc:
                print(f"ortak skalali matris cizimi basarisiz "
                      f"({cx}, {analysis}): {exc}", file=sys.stderr)
        try:
            _draw_equilibrium(cx, analysis, recs, target, conv, ulabel,
                              smooth_ns)
        except Exception as exc:
            print(f"denge cizimi basarisiz ({cx}, {analysis}): {exc}",
                  file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="mdkit cizim katmani")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("-c", "--config", type=Path, default=None,
                    help="config.sh (COMPLEX_GROUPS icin; varsayilan: mdkit/config.sh)")
    ap.add_argument("--per-complex", action="store_true")
    ap.add_argument("--mean-sd", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--profile-compare", action="store_true")
    ap.add_argument("--matrix-smooth-ns", type=float,
                    default=DEFAULT_MATRIX_SMOOTH_NS,
                    help="denge egrisinin yuvarlanan ortalama penceresi, ns "
                         f"(varsayilan {DEFAULT_MATRIX_SMOOTH_NS}; 0 kapatir)")
    ap.add_argument("--compare-output", default=DEFAULT_COMPARE_OUTPUT)
    args = ap.parse_args()

    run_all = not (args.per_complex or args.mean_sd or args.compare
                   or args.matrix or args.profile_compare)
    needs_csv = (args.per_complex or args.mean_sd or args.compare
                 or args.profile_compare or run_all)
    ts, pr = load(args.results_dir, required=needs_csv)
    groups = read_complex_groups(args.config)
    out_dir = args.results_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.per_complex or run_all:
        plot_per_complex(ts, pr, out_dir, groups)
    if args.mean_sd or run_all:
        plot_mean_sd(ts, pr, out_dir, groups)
    if args.compare or run_all:
        plot_compare(ts, out_dir, args.compare_output, groups)
    if args.profile_compare or run_all:
        plot_profile_compare(pr, out_dir, groups)
    if args.matrix or run_all:
        plot_matrix(args.results_dir, out_dir, args.matrix_smooth_ns)

    print(f"grafikler -> {out_dir}")


if __name__ == "__main__":
    main()
