#!/usr/bin/env python
"""mdkit: rep*/analysis/*.xvg -> results/{timeseries,profile}_long.csv

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

XPM_LEGEND = re.compile(r'/\*\s*legend:\s*"(.*)"\s*\*/')
XPM_AXIS = re.compile(r'/\*\s*([xy])-axis:\s*(.*?)\s*\*/')
XPM_HEADER = re.compile(r'^"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"')
XPM_VALUE = re.compile(r'/\*\s*"(.*?)"\s*\*/')

TS_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "time_ps", "value", "unit"]
PR_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "residue", "value", "unit"]


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


def collect(data_root, complex_glob, reps, manifest):
    timeseries, profile = [], []
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
                            profile.append(rec)
    return timeseries, profile


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description="mdkit .xvg toplayici")
    ap.add_argument("-c", "--config", type=Path, default=None,
                    help="config.sh (varsayilan: mdkit/config.sh)")
    ap.add_argument("-o", "--out-dir", type=Path, default=None,
                    help="cikti dizini (varsayilan: config'deki RESULTS_DIR)")
    args = ap.parse_args()

    data_root, results_dir, complex_glob, reps = read_config(args.config)
    out_dir = args.out_dir or results_dir
    manifest = read_manifest(args.config)

    timeseries, profile = collect(data_root, complex_glob, reps, manifest)
    write_csv(out_dir / "timeseries_long.csv", TS_FIELDS, timeseries)
    write_csv(out_dir / "profile_long.csv", PR_FIELDS, profile)

    print(f"timeseries: {len(timeseries)} satir -> {out_dir / 'timeseries_long.csv'}")
    print(f"profile   : {len(profile)} satir -> {out_dir / 'profile_long.csv'}")


if __name__ == "__main__":
    main()
