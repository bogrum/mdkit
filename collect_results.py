#!/usr/bin/env python
"""mdkit: rep*/analysis/*.{xvg,xpm} -> results/ CSV'leri + matrisler

    *.xvg -> results/{timeseries,profile}_long.csv
    *.xpm -> results/matrices/*.npz + results/matrix_summary.csv

Analiz manifestosu `run_analysis.sh --list` ciktisindan okunur; boylece hangi
ciktinin hangi analize ve hangi KIND'a ait oldugu tek kaynakta (analysis/*.sh)
kalir ve burada kopyalanmaz.

Degerler .xvg'deki ham haliyle yazilir (zaman ps, mesafe nm) ve birim `unit`
kolonunda belirtilir. Donusum cizim katmaninda yapilir: korlemesine nm->A
carpmak ileride gmx sasa (nm^2) eklendiginde sessizce yanlis olurdu.
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

MDKIT = Path(__file__).resolve().parent

XVG_LEGEND = re.compile(r'@\s+s(\d+)\s+legend\s+"(.*)"')
XVG_LABEL = re.compile(r'@\s+(title|subtitle)\s+"(.*)"')
XVG_AXIS = re.compile(r'@\s+(xaxis|yaxis)\s+label\s+"(.*)"')
UNIT_IN_LABEL = re.compile(r"\(([^)]*)\)")

XPM_TITLE = re.compile(r'/\*\s*title:\s*"(.*)"\s*\*/')
XPM_SUBTITLE = re.compile(r'/\*\s*subtitle:\s*"(.*)"\s*\*/')
XPM_LEGEND = re.compile(r'/\*\s*legend:\s*"(.*)"\s*\*/')
XPM_AXIS = re.compile(r'/\*\s*([xy])-axis:\s*(.*?)\s*\*/')
XPM_HEADER = re.compile(r'^"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"')
XPM_VALUE = re.compile(r'/\*\s*"(.*?)"\s*\*/')

TS_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "time_ps", "value", "unit"]
# residue_from_end: C-ucundan sayan konum (son residue 0, oncesi -1, ...).
# Peptidler 8-11 residue arasinda degistigi icin HAM residue numarasi
# kompleksler arasi karsilastirilamaz: bir 8-mer'in 5. residue'su ortada,
# 11-mer'in 5.'si degil. Iki kolon birlikte her turlu konum normalizasyonunu
# turetilebilir kilar. ETIKET (P1/PO gibi) burada TUTULMAZ -- o bir yorumdur
# ve cizim katmanina aittir.
PR_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "residue", "residue_from_end", "value", "unit"]
MX_FIELDS = ["complex", "replica_i", "replica_j", "analysis", "output",
             "n_x", "n_y", "dt_ps", "min", "mean", "max", "mean_vs_self",
             "unit"]


def parse_xvg(path):
    """xvg dosyasini (meta, satirlar) olarak dondurur."""
    meta = {"legends": {}}
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            m = XVG_LEGEND.match(line)
            if m:
                meta["legends"][int(m.group(1))] = m.group(2)
                continue
            m = XVG_LABEL.match(line)
            if m:
                meta[m.group(1)] = m.group(2)
                continue
            m = XVG_AXIS.match(line)
            if m:
                meta[m.group(1)] = m.group(2)
            continue
        try:
            rows.append([float(p) for p in line.split()])
        except ValueError:
            continue
    return meta, rows


def parse_xpm(path):
    """gmx .xpm matrisini (meta, values, x_ps, y_ps) olarak dondurur.

    UC TUZAK, ucu de GROMACS 2025.4 uzerinde olculerek saptandi (spec 2.4):

    1. Uzun eksenler BIRDEN FAZLA `/* x-axis: */` yorumuna bolunur. Yalnizca
       ilkini okuyan bir ayristirici ekseni sessizce keser.
    2. Piksel satirlari y ekseninin TERSI sirada yazilir. Cevrilmezse matris
       yatayda aynalanir -- self-matris disinda gozle fark edilmez.
    3. Karakter alani CPP genisliginde SABITTIR ve bosluk da gecerli bir
       karakterdir. split() ile ayristirmak CPP>1'de anahtari bozar.
    """
    colors = {}
    axes = {"x": [], "y": []}
    pixel_rows = []
    meta = {}
    width = height = cpp = None
    colors_left = 0

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("/*"):
            # gmx buraya NEYIN olculdugunu yazar (or. "LIGAND_BB RMSD
            # matrix"). Figurde gosterilir; okuyucu aksi halde grafikten
            # ne oldugunu anlayamaz.
            # 'subtitle' ONCE denenir: gmx'in yazdigi satir degil, eklentinin
            # ekledigi satirdir ve HEM olculen HEM fit grubunu tasir.
            m = XPM_SUBTITLE.search(line)
            if m:
                meta["subtitle"] = m.group(1)
                continue
            m = XPM_TITLE.search(line)
            if m:
                meta["title"] = m.group(1)
                continue
            m = XPM_LEGEND.search(line)
            if m:
                meta["legend"] = m.group(1)
                continue
            m = XPM_AXIS.search(line)
            if m:
                # TUZAK 1: extend, assign DEGIL.
                axes[m.group(1)].extend(float(v) for v in m.group(2).split())
            continue
        if not line.startswith('"'):
            continue
        if width is None:
            m = XPM_HEADER.match(line)
            if m:
                width, height, ncolors, cpp = (int(g) for g in m.groups())
                colors_left = ncolors
            continue
        if colors_left > 0:
            # TUZAK 3: konumsal dilim.
            key = line[1:1 + cpp]
            vm = XPM_VALUE.search(line)
            if vm is None:
                raise ValueError(f"{path.name}: renk satiri cozulemedi: {line}")
            colors[key] = float(vm.group(1))
            colors_left -= 1
            continue
        pixel_rows.append(line[1:1 + width * cpp])

    if width is None:
        raise ValueError(f"{path.name}: .xpm basligi bulunamadi")
    if len(pixel_rows) != height:
        raise ValueError(
            f"{path.name}: basligi {height} satir diyor, {len(pixel_rows)} bulundu"
        )

    values = np.empty((height, width), dtype=np.float32)
    for r, row in enumerate(pixel_rows):
        for c in range(width):
            key = row[c * cpp:(c + 1) * cpp]
            try:
                values[r, c] = colors[key]
            except KeyError:
                raise ValueError(
                    f"{path.name}: renk tablosunda olmayan karakter {key!r}"
                ) from None
    # TUZAK 2.
    values = values[::-1].copy()

    x_ps = np.asarray(axes["x"], dtype=np.float64)
    y_ps = np.asarray(axes["y"], dtype=np.float64)
    if values.shape != (len(y_ps), len(x_ps)):
        raise ValueError(
            f"{path.name}: matris {values.shape}, eksenler "
            f"({len(y_ps)}, {len(x_ps)}) -- uyusmuyor"
        )

    m = UNIT_IN_LABEL.search(meta.get("legend", ""))
    meta["unit"] = m.group(1) if m else ""
    meta.setdefault("title", "")
    meta.setdefault("subtitle", "")
    return meta, values, x_ps, y_ps


def _bash(snippet):
    r = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr.strip() or "bash cagrisi basarisiz")
    return r.stdout


def read_config(config):
    """Config'i bash tarafindaki dogrulayiciyla okur (tek kaynak)."""
    arg = f'"{config}"' if config else ""
    out = _bash(
        f'source "{MDKIT}/analysis/lib.sh" && mdkit_load_config {arg} && '
        'printf "%s\\n%s\\n%s\\n%s\\n" '
        '"$DATA_ROOT" "$RESULTS_DIR" "$COMPLEX_GLOB" "${REPS[*]}"'
    )
    data_root, results_dir, complex_glob, reps = out.splitlines()[:4]
    return Path(data_root), Path(results_dir), complex_glob, reps.split()


def read_manifest(config):
    """cikti_dosya_adi -> (analiz_adi, KIND)"""
    cmd = ["bash", str(MDKIT / "run_analysis.sh"), "--list"]
    if config:
        cmd += ["--config", str(config)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr.strip() or "analiz manifestosu okunamadi")
    manifest = {}
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        # Ilk dort alan sabit; 5. alan (opsiyonel ciktilar) sonradan eklendi
        # ve zaten 4. alanin bir alt kumesidir -- eski bir --list ciktisiyla
        # da calissin diye varligi zorunlu tutulmaz.
        parts = line.split("\t")
        name, kind, _desc, outputs = parts[:4]
        optional = parts[4] if len(parts) > 4 else ""
        for out in outputs.split(",") + optional.split(","):
            if out.strip():
                manifest[out.strip()] = (name, kind)
    return manifest


def parse_bin(path, n_x, n_y):
    """gmx rms -bin ham dump'ini (n_y, n_x) float32 matrise cevirir.

    Duzen GROMACS 2025.4'te OLCULDU: dosya BASLIKSIZDIR, float32'dir ve
    [x][y] sirasinda yazilir -- yani .xpm'in TRANSPOZESI. .xpm'den farkli
    olarak satir cevirme (rows[::-1]) GEREKMEZ; cevirmek matrisi bozar.

    Baslik olmadigi icin boyutlar disaridan gelir (.xpm'den) ve tek
    dogrulama dosya boyutudur: yanlis sekillendirilen bir matris sessizce
    cop olurdu.
    """
    raw = np.fromfile(path, dtype=np.float32)
    if raw.size != n_x * n_y:
        raise ValueError(
            f"{path.name}: {raw.size} deger var, {n_x}x{n_y}={n_x * n_y} bekleniyordu"
        )
    return raw.reshape(n_x, n_y).T.copy()


def peer_replica(output_name, reps):
    """'cross_rmsd_rep2.xpm' -> 'rep2'; eslesme yoksa None.

    Dosya adi kalibini burada yeniden tanimlamiyoruz: config'den gelen
    BILINEN replika adlariyla eslestiriyoruz. Boylece `matrix` destegi
    replikalar-arasi analizlere kilitlenmez -- eslesmeyen bir matris
    (or. residue x zaman) da toplanir, yalnizca replica_j'si bos olur."""
    stem = Path(output_name).stem
    for rep in reps:
        if stem.endswith("_" + rep):
            return rep
    return None


def collect(data_root, complex_glob, reps, manifest, matrix_dir=None):
    timeseries, profile, matrices = [], [], []
    warned = set()
    warned_unlisted = set()
    for cx in sorted(data_root.glob(complex_glob)):
        if not cx.is_dir():
            continue
        cname = cx.name.split("_")[0]
        for rep in reps:
            adir = cx / rep / "analysis"
            if not adir.is_dir():
                continue

            if matrix_dir is not None:
                for xpm in sorted(adir.glob("*.xpm")):
                    entry = manifest.get(xpm.name)
                    if entry is None:
                        if xpm.name not in warned_unlisted:
                            warned_unlisted.add(xpm.name)
                            print(
                                f"mdkit: manifestoda olmayan .xpm atlandi: "
                                f"{xpm.name} (ilk gorulen: {adir})",
                                file=sys.stderr,
                            )
                        continue
                    analysis, kind = entry
                    if kind != "matrix":
                        key = (xpm.name, kind)
                        if key not in warned:
                            warned.add(key)
                            print(
                                f"mdkit: .xpm ciktisi {kind!r} kind'i ile ilan "
                                f"edilmis ({xpm.name}) -- toplanmadan atlandi",
                                file=sys.stderr,
                            )
                        continue
                    # Tek bir bozuk matris 105 dizinlik toplamayi dusurmemeli
                    # (bash tarafindaki hata yalitimiyla ayni gerekce).
                    try:
                        rec = collect_matrix(
                            xpm, cname, rep, reps, analysis, matrix_dir)
                    except Exception as exc:
                        print(f"mdkit: {xpm.name} okunamadi ({adir}): {exc}",
                              file=sys.stderr)
                        continue
                    matrices.append(rec)

            for xvg in sorted(adir.glob("*.xvg")):
                entry = manifest.get(xvg.name)
                if entry is None:
                    # Sessiz atlama, manifesto ile disk arasindaki her
                    # uyusmazligi gorunmez kiliyordu: bir eklenti ciktisini
                    # kosullu ilan ettiginde saatler suren hesaplar hicbir
                    # CSV'ye girmeden kayboluyordu. Dosya adi basina bir satir.
                    if xvg.name not in warned_unlisted:
                        warned_unlisted.add(xvg.name)
                        print(
                            f"mdkit: manifestoda olmayan .xvg atlandi: "
                            f"{xvg.name} (ilk gorulen: {adir})",
                            file=sys.stderr,
                        )
                    continue
                analysis, kind = entry
                if kind not in ("timeseries", "profile"):
                    key = (xvg.name, kind)
                    if key not in warned:
                        warned.add(key)
                        print(
                            f"mdkit: taninmayan ANALYSIS_KIND {kind!r} "
                            f"({xvg.name}) -- toplanmadan atlandi",
                            file=sys.stderr,
                        )
                    continue
                meta, rows = parse_xvg(xvg)
                m = UNIT_IN_LABEL.search(meta.get("yaxis", ""))
                unit = m.group(1) if m else ""
                # Uzunluk DOSYA basina okunur: ayni kosuda farkli uzunlukta
                # profiller (peptid ~10, MHC ~275) bir arada bulunur.
                n_res = int(max(r[0] for r in rows)) if rows else 0
                for row in rows:
                    x, ys = row[0], row[1:]
                    for i, y in enumerate(ys):
                        series = (meta["legends"].get(i)
                                  or meta.get("subtitle")
                                  or f"y{i}")
                        rec = {
                            "complex": cname, "replica": rep,
                            "analysis": analysis, "output": xvg.name,
                            "series": series, "value": y, "unit": unit,
                        }
                        if kind == "timeseries":
                            rec["time_ps"] = x
                            timeseries.append(rec)
                        elif kind == "profile":
                            rec["residue"] = int(x)
                            rec["residue_from_end"] = int(x) - n_res
                            profile.append(rec)

    # Oran ancak TUM matrisler toplandiktan sonra hesaplanabilir: bir
    # capraz ciftin paydasi, iki farkli replikanin self matrislerinden gelir.
    add_self_ratio(matrices)
    return timeseries, profile, matrices


def add_self_ratio(matrices):
    """Her matrise `mean_vs_self` ekler: ortalamanin, ilgili SELF
    matrislerin ortalamasina orani.

        mean_vs_self = mean / ((self_i + self_j) / 2)

    "Replikalar ayni konformasyonel alani mi ornekliyor?" sorusunun NICEL
    cevabi. 1'e yakin = evet (capraz uzaklik, replika ICI uzakliktan farkli
    degil); buyudukce = hayir (replikalar farkli havzalarda). Isi
    haritasina bakip goz karariyla soylemek yerine sayiyla soylenebilir --
    ve renk skalasi kompleksten komplekse degistigi icin goz karari zaten
    guvenilmez.

    Self satirlarda tanim geregi 1.0. Referans self matrisi yoksa (or.
    replika adiyla eslesmeyen bir matris, replica_j bos) oran UYDURULMAZ,
    bos birakilir.

    DIKKAT: self ortalamalari kosegen HARIC hesaplanir (bkz. collect_matrix),
    capraz ortalamalar ise tum hucreleri icerir. Payda bu yuzden bir miktar
    BUYUK, yani oran biraz muhafazakar -- heterojenligi abartmaz.
    """
    selfs = {}
    for r in matrices:
        if r["replica_i"] == r["replica_j"] and r["replica_i"]:
            selfs[(r["complex"], r["analysis"], r["replica_i"])] = r["mean"]

    for r in matrices:
        a = selfs.get((r["complex"], r["analysis"], r["replica_i"]))
        b = selfs.get((r["complex"], r["analysis"], r["replica_j"]))
        if a is None or b is None or (a + b) <= 0:
            r["mean_vs_self"] = ""
        else:
            r["mean_vs_self"] = r["mean"] / ((a + b) / 2)
    return matrices


def collect_matrix(xpm, cname, rep, reps, analysis, matrix_dir):
    """Bir .xpm'i .npz olarak yazar ve ozet kaydini dondurur.

    Matrisler uzun-format CSV'ye GIRMEZ: 451x451'lik 315 matris ~64 milyon
    satir ederdi. 1D veri icin dogru olan bicim 2B icin degil."""
    meta, values, x_ps, y_ps = parse_xpm(xpm)

    # .xpm degerleri 80 renk seviyesine yuvarlanmistir. Kardes .dat ayni
    # matrisin ham float32 halidir; varsa DEGERLER ondan alinir. Eksenler,
    # birim ve sekil yine .xpm'den gelir -- .dat'ta baslik bile yoktur.
    dat = xpm.with_suffix(".dat")
    if dat.is_file():
        try:
            values = parse_bin(dat, len(x_ps), len(y_ps))
        except ValueError as exc:
            # .xpm dogru bir yedektir, yalnizca daha kaba. Toplamayi
            # dusurmek 105 dizinlik kosuyu tek bozuk dosyaya feda ederdi.
            print(f"mdkit: {dat.name} kullanilamadi, .xpm degerlerine "
                  f"donuldu ({exc})", file=sys.stderr)

    # dt katmanlar arasinda tasinmaz, eksen araligindan geri okunur.
    dt_ps = float(x_ps[1] - x_ps[0]) if len(x_ps) > 1 else ""
    peer = peer_replica(xpm.name, reps)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    out = matrix_dir / f"{cname}_{rep}_{xpm.stem}.npz"
    np.savez_compressed(
        out,
        values=values, x_ps=x_ps, y_ps=y_ps,
        unit=meta["unit"], complex=cname, replica_i=rep,
        replica_j=peer or "", analysis=analysis, output=xpm.name,
        dt_ps=dt_ps, title=meta["title"], subtitle=meta["subtitle"],
    )
    # Self-matriste kosegen TANIM GEREGI sifirdir ve min'i anlamsiz kilar;
    # tezde kullanilacak sayi capraz ciftin minimumudur.
    if peer == rep and values.shape[0] == values.shape[1]:
        sample = values[~np.eye(values.shape[0], dtype=bool)]
    else:
        sample = values.ravel()

    return {
        "complex": cname, "replica_i": rep, "replica_j": peer or "",
        "analysis": analysis, "output": xpm.name,
        "n_x": int(values.shape[1]), "n_y": int(values.shape[0]),
        "dt_ps": dt_ps,
        "min": float(sample.min()), "mean": float(sample.mean()),
        "max": float(sample.max()), "unit": meta["unit"],
    }


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description="mdkit .xvg/.xpm toplayici")
    ap.add_argument("-c", "--config", type=Path, default=None,
                    help="config.sh (varsayilan: mdkit/config.sh)")
    ap.add_argument("-o", "--out-dir", type=Path, default=None,
                    help="cikti dizini (varsayilan: config'deki RESULTS_DIR)")
    args = ap.parse_args()

    data_root, results_dir, complex_glob, reps = read_config(args.config)
    out_dir = args.out_dir or results_dir
    manifest = read_manifest(args.config)

    matrix_dir = out_dir / "matrices"
    timeseries, profile, matrices = collect(
        data_root, complex_glob, reps, manifest, matrix_dir=matrix_dir)
    write_csv(out_dir / "timeseries_long.csv", TS_FIELDS, timeseries)
    write_csv(out_dir / "profile_long.csv", PR_FIELDS, profile)
    write_csv(out_dir / "matrix_summary.csv", MX_FIELDS, matrices)

    print(f"timeseries: {len(timeseries)} satir -> {out_dir / 'timeseries_long.csv'}")
    print(f"profile   : {len(profile)} satir -> {out_dir / 'profile_long.csv'}")
    print(f"matris    : {len(matrices)} dosya -> {matrix_dir} "
          f"(+ {out_dir / 'matrix_summary.csv'})")


if __name__ == "__main__":
    main()
