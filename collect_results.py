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

MDKIT = Path(__file__).resolve().parent

XVG_LEGEND = re.compile(r'@\s+s(\d+)\s+legend\s+"(.*)"')
XVG_LABEL = re.compile(r'@\s+(title|subtitle)\s+"(.*)"')
XVG_AXIS = re.compile(r'@\s+(xaxis|yaxis)\s+label\s+"(.*)"')
UNIT_IN_LABEL = re.compile(r"\(([^)]*)\)")

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
        name, kind, _desc, outputs = line.split("\t")
        for out in outputs.split(","):
            manifest[out.strip()] = (name, kind)
    return manifest


def collect(data_root, complex_glob, reps, manifest):
    timeseries, profile = [], []
    warned = set()
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
