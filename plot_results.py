#!/usr/bin/env python
"""mdkit cizim katmani.

Cizim ANALIZ ADINA gore degil, veri SEKLINE gore dallanir: timeseries ve
profile. Yeni bir timeseries analizi (gyrate, sasa) eklendiginde bu dosyaya
dokunmak gerekmez.

Birim donusumu burada yapilir: .xvg'ler ps ve nm cinsindendir (spec 2.5).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

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


def scale_and_label(unit):
    """(carpan, eksen etiketi). Yalnizca nm cevrilir; bilinmeyen birim
    cevrilmeden, kendi etiketiyle cizilir -- bilmedigimiz bir birimi
    cevirmis gibi yapmiyoruz."""
    if unit == "nm":
        return NM_TO_ANGSTROM, "Å"
    return 1.0, unit or "birimsiz"


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
    ts = pd.read_csv(ts_path) if ts_path.exists() else pd.DataFrame()
    pr = pd.read_csv(pr_path) if pr_path.exists() else pd.DataFrame()
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
                ax.set_title(f"{cx} — {Path(output).stem}")
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
                ax.set_title(f"{cx} — {Path(output).stem} (n={g['replica'].nunique()} replika)")
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
                }
        except Exception as exc:
            print(f"matris okunamadi ({f.name}): {exc}", file=sys.stderr)
            continue
        groups.setdefault((rec["complex"], rec["analysis"]), []).append(rec)
    return groups


def plot_matrix(results_dir, out_dir):
    """Kompleks basina N x N isi haritasi izgarasi, ORTAK renk skalasiyla.

    Ortak skala sart: panel basina ayri skala, farkli replika ciftlerini
    gorsel olarak karsilastirilamaz kilardi -- bu figurun tek amaci o
    karsilastirma."""
    groups = load_matrices(results_dir)
    if not groups:
        return
    target = out_dir / "matrix"
    target.mkdir(parents=True, exist_ok=True)

    for (cx, analysis), recs in sorted(groups.items()):
        fig = None
        try:
            rows = sorted({r["replica_j"] for r in recs})
            cols = sorted({r["replica_i"] for r in recs})
            by_cell = {(r["replica_i"], r["replica_j"]): r for r in recs}

            unit = recs[0]["unit"]
            conv, ulabel = scale_and_label(unit)
            if unit != "nm":
                _warn_unknown_unit(f"{cx} ({analysis})", unit)
            vmin = min(float(r["values"].min()) for r in recs) * conv
            vmax = max(float(r["values"].max()) for r in recs) * conv

            fig, axes = plt.subplots(
                len(rows), len(cols), squeeze=False,
                figsize=(2.6 * len(cols) + 1.6, 2.6 * len(rows)),
                sharex=True, sharey=True,
            )
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
            fig.suptitle(f"{cx} \u2014 {analysis}")
            fig.savefig(target / f"{cx}_{analysis}.png", dpi=150,
                        bbox_inches="tight")
        except Exception as exc:
            print(f"matris cizimi basarisiz ({cx}, {analysis}): {exc}",
                  file=sys.stderr)
        finally:
            if fig is not None:
                plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="mdkit cizim katmani")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("-c", "--config", type=Path, default=None,
                    help="config.sh (COMPLEX_GROUPS icin; varsayilan: mdkit/config.sh)")
    ap.add_argument("--per-complex", action="store_true")
    ap.add_argument("--mean-sd", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--compare-output", default=DEFAULT_COMPARE_OUTPUT)
    args = ap.parse_args()

    run_all = not (args.per_complex or args.mean_sd or args.compare
                   or args.matrix)
    needs_csv = args.per_complex or args.mean_sd or args.compare or run_all
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
    if args.matrix or run_all:
        plot_matrix(args.results_dir, out_dir)

    print(f"grafikler -> {out_dir}")


if __name__ == "__main__":
    main()
