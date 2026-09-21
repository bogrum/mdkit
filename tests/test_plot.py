import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_results  # noqa: E402

pd = pytest.importorskip("pandas")

TS_HEADER = "complex,replica,analysis,output,series,time_ps,value,unit\n"
PR_HEADER = "complex,replica,analysis,output,series,residue,value,unit\n"


@pytest.fixture
def results_dir(tmp_path):
    d = tmp_path / "results"
    d.mkdir()
    ts_rows, pr_rows = [], []
    for cx, base in [("last1", 0.20), ("top1", 0.10), ("last2", 0.30)]:
        for i, rep in enumerate(["rep1", "rep2", "rep3"]):
            for t in range(0, 50, 10):
                ts_rows.append(
                    f"{cx},{rep},rmsd,rmsd_pep_on_mhc.xvg,LIGAND,{t}.0,"
                    f"{base + 0.01 * i + 0.001 * t},nm"
                )
            for res in range(1, 6):
                pr_rows.append(
                    f"{cx},{rep},rmsf,rmsf_pep_self.xvg,LIGAND,{res},"
                    f"{base + 0.02 * res},nm"
                )
    (d / "timeseries_long.csv").write_text(TS_HEADER + "\n".join(ts_rows) + "\n")
    (d / "profile_long.csv").write_text(PR_HEADER + "\n".join(pr_rows) + "\n")
    return d


def run_plot(mdkit, results_dir, *args):
    return subprocess.run(
        [sys.executable, str(mdkit / "plot_results.py"),
         "--results-dir", str(results_dir), *args],
        capture_output=True, text=True,
    )


def test_compare_order_medyana_gore_siralar():
    df = pd.DataFrame({
        "complex": ["a", "a", "b", "b", "c", "c"],
        "value": [3.0, 3.2, 1.0, 1.1, 2.0, 2.2],
    })
    assert plot_results.compare_order(df) == ["b", "c", "a"]


def test_per_complex_kompleks_basina_figure_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--per-complex")
    assert r.returncode == 0, r.stderr
    out = results_dir / "plots" / "per_complex"
    names = sorted(p.name for p in out.glob("*.png"))
    assert "last1_rmsd_pep_on_mhc.png" in names
    assert "last1_rmsf_pep_self.png" in names
    assert len(names) == 3 * 2          # 3 kompleks x 2 cikti
    assert all((out / n).stat().st_size > 1000 for n in names)


def test_mean_sd_figure_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--mean-sd")
    assert r.returncode == 0, r.stderr
    out = results_dir / "plots" / "mean_sd"
    assert len(list(out.glob("*.png"))) == 3 * 2


def test_compare_tek_panel_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--compare")
    assert r.returncode == 0, r.stderr
    p = results_dir / "plots" / "compare_rmsd_pep_on_mhc.png"
    assert p.exists() and p.stat().st_size > 1000


def test_mod_verilmezse_hepsi_calisir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir)
    assert r.returncode == 0, r.stderr
    assert (results_dir / "plots" / "per_complex").is_dir()
    assert (results_dir / "plots" / "mean_sd").is_dir()
    assert (results_dir / "plots" / "compare_rmsd_pep_on_mhc.png").exists()


def test_eksik_csv_anlasilir_hata(mdkit, tmp_path):
    r = run_plot(mdkit, tmp_path / "yok")
    assert r.returncode != 0
    assert "bulunamadi" in (r.stderr + r.stdout)


# --- Fix 1: unit-aware scaling -------------------------------------------

def test_scale_and_label_nm_angstroma_cevirir():
    assert plot_results.scale_and_label("nm") == (10.0, "Å")


def test_scale_and_label_bilinmeyen_birimi_cevirmeden_birakir():
    assert plot_results.scale_and_label("nm^2") == (1.0, "nm^2")


def test_scale_and_label_bos_birim_birimsiz_doner():
    assert plot_results.scale_and_label("") == (1.0, "birimsiz")


@pytest.fixture
def results_dir_odd_unit(tmp_path):
    """timeseries'i nm disinda bir birimle (nm^2, sasa benzeri) tasiyan set."""
    d = tmp_path / "results"
    d.mkdir()
    rows = []
    for rep in ["rep1", "rep2"]:
        for t in range(0, 30, 10):
            rows.append(
                f"only1,{rep},sasa,sasa_total.xvg,LIGAND,{t}.0,"
                f"{1.5 + 0.01 * t},nm^2"
            )
    (d / "timeseries_long.csv").write_text(TS_HEADER + "\n".join(rows) + "\n")
    (d / "profile_long.csv").write_text(PR_HEADER)
    return d


def test_bilinmeyen_birim_cevrilmeden_cizilir_ve_stderrde_ciktiyi_adlandirir(
    mdkit, results_dir_odd_unit
):
    r = run_plot(mdkit, results_dir_odd_unit, "--per-complex")
    assert r.returncode == 0, r.stderr
    p = results_dir_odd_unit / "plots" / "per_complex" / "only1_sasa_total.png"
    assert p.exists() and p.stat().st_size > 1000
    assert "sasa_total.xvg" in r.stderr


# --- Fix 2: error isolation in the render loops --------------------------

def _rows(cx, output, unit="nm", n=3, base=0.1):
    return [
        {
            "complex": cx, "replica": "rep1", "analysis": "rmsd",
            "output": output, "series": "LIGAND",
            "time_ps": float(t * 10), "value": base + 0.001 * t, "unit": unit,
        }
        for t in range(n)
    ]


def test_plot_per_complex_bir_kompleksin_hatasi_digerlerini_durdurmaz(tmp_path):
    output = "rmsd_pep_on_mhc.xvg"
    rows = (_rows("good1", output, base=0.10)
            + _rows("bad/name", output, base=0.15)
            + _rows("good2", output, base=0.20))
    ts = pd.DataFrame(rows)
    pr = pd.DataFrame(columns=[
        "complex", "replica", "analysis", "output", "series",
        "residue", "value", "unit",
    ])
    out_dir = tmp_path / "plots"

    plot_results.plot_per_complex(ts, pr, out_dir)  # not raise

    target = out_dir / "per_complex"
    assert (target / f"good1_{Path(output).stem}.png").exists()
    assert (target / f"good2_{Path(output).stem}.png").exists()


def test_plot_mean_sd_tek_replika_std_nan_ile_bile_calisir(tmp_path):
    output = "rmsd_pep_on_mhc.xvg"
    ts = pd.DataFrame(_rows("solo", output, base=0.15))
    pr = pd.DataFrame(columns=[
        "complex", "replica", "analysis", "output", "series",
        "residue", "value", "unit",
    ])
    out_dir = tmp_path / "plots"

    plot_results.plot_mean_sd(ts, pr, out_dir)  # not raise

    p = out_dir / "mean_sd" / f"solo_{Path(output).stem}.png"
    assert p.exists() and p.stat().st_size > 1000
