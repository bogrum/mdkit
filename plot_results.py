#!/usr/bin/env python
"""mdkit cizim katmani.

Cizim ANALIZ ADINA gore degil, veri SEKLINE gore dallanir: timeseries ve
profile. Yeni bir timeseries analizi (gyrate, sasa) eklendiginde bu dosyaya
dokunmak gerekmez.

Birim donusumu burada yapilir: .xvg'ler ps ve nm cinsindendir (spec 2.5).
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

NM_TO_ANGSTROM = 10.0
PS_TO_NS = 1e-3

# dataviz skill'inin dogrulanmis paleti (categorical, slot 1/2/3 ve 1/8).
# Sadece bu iki sabit degistirilir; kodun geri kalani ayni kalir.
REP_COLORS = {"rep1": "#2a78d6", "rep2": "#eb6834", "rep3": "#1baf7a"}
GROUP_COLORS = {"top": "#2a78d6", "last": "#e34948"}

DEFAULT_COMPARE_OUTPUT = "rmsd_pep_on_mhc.xvg"


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


def load(results_dir):
    ts_path = results_dir / "timeseries_long.csv"
    pr_path = results_dir / "profile_long.csv"
    if not ts_path.exists() and not pr_path.exists():
        sys.exit(f"sonuc CSV'leri bulunamadi: {results_dir}")
    ts = pd.read_csv(ts_path) if ts_path.exists() else pd.DataFrame()
    pr = pd.read_csv(pr_path) if pr_path.exists() else pd.DataFrame()
    return ts, pr


def _group_color(complex_name):
    return GROUP_COLORS["top" if complex_name.startswith("top") else "last"]


def _panels(ts, pr):
    """(df, x kolonu, x etiketi, x carpani) uclusu uretir."""
    if not ts.empty:
        yield ts, "time_ps", "Zaman (ns)", PS_TO_NS
    if not pr.empty:
        yield pr, "residue", "Residue", 1.0


def plot_per_complex(ts, pr, out_dir):
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
                for rep, gr in g.groupby("replica", sort=True):
                    gr = gr.sort_values(xcol)
                    ax.plot(gr[xcol] * xconv, gr["value"] * yconv,
                            label=rep, color=REP_COLORS.get(rep), linewidth=1.0)
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


def plot_mean_sd(ts, pr, out_dir):
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
                stats = (g.groupby(xcol)["value"]
                           .agg(["mean", "std"])
                           .sort_index()
                           .fillna(0.0))
                x = stats.index.to_numpy() * xconv
                mean = stats["mean"].to_numpy() * yconv
                sd = stats["std"].to_numpy() * yconv
                color = _group_color(cx)

                fig, ax = plt.subplots(figsize=(7, 4))
                ax.fill_between(x, mean - sd, mean + sd, color=color, alpha=0.25,
                                linewidth=0)
                ax.plot(x, mean, color=color, linewidth=1.4)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                ax.set_title(f"{cx} — {Path(output).stem} (n={g['replica'].nunique()} replika)")
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


def plot_compare(ts, out_dir, output=DEFAULT_COMPARE_OUTPUT):
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
        patch.set_facecolor(_group_color(cx))
        patch.set_alpha(0.75)
        patch.set_edgecolor("#333333")
    for median in bp["medians"]:
        median.set_color("#222222")

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel(f"Ortalama {Path(output).stem} ({ylabel})")
    ax.set_title("Kompleksler arasi karsilastirma (replika basina ortalama)")
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], color=c, linewidth=6, alpha=0.75)
               for c in (GROUP_COLORS["top"], GROUP_COLORS["last"])]
    ax.legend(handles, ["top*", "last*"], frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / f"compare_{Path(output).stem}.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="mdkit cizim katmani")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--per-complex", action="store_true")
    ap.add_argument("--mean-sd", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--compare-output", default=DEFAULT_COMPARE_OUTPUT)
    args = ap.parse_args()

    ts, pr = load(args.results_dir)
    out_dir = args.results_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_all = not (args.per_complex or args.mean_sd or args.compare)
    if args.per_complex or run_all:
        plot_per_complex(ts, pr, out_dir)
    if args.mean_sd or run_all:
        plot_mean_sd(ts, pr, out_dir)
    if args.compare or run_all:
        plot_compare(ts, out_dir, args.compare_output)

    print(f"grafikler -> {out_dir}")


if __name__ == "__main__":
    main()
