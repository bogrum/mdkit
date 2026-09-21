import csv
import subprocess
import sys
import textwrap

import pytest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import collect_results  # noqa: E402

TIMESERIES_XVG = textwrap.dedent("""\
    # gmx rms tarafindan uretildi
    @    title "RMSD"
    @    xaxis  label "Time (ps)"
    @    yaxis  label "RMSD (nm)"
    @ subtitle "LIGAND after lsq fit to RECEPTOR_BB"
    @TYPE xy
       0.0000000    0.0000005
      10.0000000    0.1457492
      20.0000000    0.1278536
    """)

PROFILE_XVG = textwrap.dedent("""\
    @    title "RMS fluctuation"
    @    xaxis  label "Residue"
    @    yaxis  label "(nm)"
    @TYPE xy
       1   0.1823
       2   0.1490
    """)

MULTI_SERIES_XVG = textwrap.dedent("""\
    @    xaxis  label "Time (ps)"
    @    yaxis  label "RMSD (nm)"
    @ s0 legend "birinci"
    @ s1 legend "ikinci"
       0.0   0.10   0.20
      10.0   0.11   0.21
    """)


def write_outputs(fake_dataset, mapping):
    """mapping: {dosya_adi: icerik} -> her kompleks/replikaya yazar."""
    for cx in fake_dataset.iterdir():
        if not cx.name.endswith("_pandora"):
            continue
        for rep in ["rep1", "rep2", "rep3"]:
            adir = cx / rep / "analysis"
            adir.mkdir(parents=True, exist_ok=True)
            for name, content in mapping.items():
                (adir / name).write_text(content)


def test_parse_xvg_yorum_ve_basliklari_atlar(tmp_path):
    f = tmp_path / "a.xvg"
    f.write_text(TIMESERIES_XVG)
    meta, rows = collect_results.parse_xvg(f)
    assert len(rows) == 3
    assert rows[1] == [10.0, 0.1457492]
    assert meta["yaxis"] == "RMSD (nm)"
    assert meta["subtitle"] == "LIGAND after lsq fit to RECEPTOR_BB"


def test_parse_xvg_legend_okur(tmp_path):
    f = tmp_path / "b.xvg"
    f.write_text(MULTI_SERIES_XVG)
    meta, rows = collect_results.parse_xvg(f)
    assert meta["legends"] == {0: "birinci", 1: "ikinci"}
    assert rows[0] == [0.0, 0.10, 0.20]


def test_manifest_runner_listesinden_okunur(mdkit, fake_config, monkeypatch):
    monkeypatch.setattr(collect_results, "MDKIT", mdkit)
    manifest = collect_results.read_manifest(fake_config)
    assert manifest["rmsd_pep_on_mhc.xvg"] == ("rmsd", "timeseries")
    assert manifest["rmsf_pep_self.xvg"] == ("rmsf", "profile")


def run_collect(mdkit, config):
    return subprocess.run(
        [sys.executable, str(mdkit / "collect_results.py"), "-c", str(config)],
        capture_output=True, text=True,
    )


def test_timeseries_ve_profile_ayri_dosyalara_gider(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "rmsf_pep_self.xvg": PROFILE_XVG,
    })
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr

    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    pr = list(csv.DictReader((tmp_path / "results" / "profile_long.csv").open()))
    assert len(ts) == 2 * 3 * 3   # 2 kompleks x 3 replika x 3 satir
    assert len(pr) == 2 * 3 * 2
    assert ts[0]["analysis"] == "rmsd"
    assert ts[0]["time_ps"] == "0.0"
    assert ts[0]["unit"] == "nm"
    assert pr[0]["residue"] == "1"
    assert pr[0]["analysis"] == "rmsf"


def test_ham_deger_korunur_donusturulmez(mdkit, fake_config, fake_dataset, tmp_path):
    """nm->A donusumu burada YAPILMAZ; ileride sasa (nm^2) eklendiginde
    korlemesine carpma sessizce yanlis olurdu."""
    write_outputs(fake_dataset, {"rmsd_pep_on_mhc.xvg": TIMESERIES_XVG})
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    assert float(ts[1]["value"]) == pytest.approx(0.1457492)


def test_manifestte_olmayan_xvg_yok_sayilir(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "elle_yazilmis_baska.xvg": TIMESERIES_XVG,
    })
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    assert {row["output"] for row in ts} == {"rmsd_pep_on_mhc.xvg"}


def test_cok_serili_xvg_her_seriyi_ayri_satira_yazar(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {"rmsd_pep_on_mhc.xvg": MULTI_SERIES_XVG})
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    series = {row["series"] for row in ts}
    assert series == {"birinci", "ikinci"}
    assert len(ts) == 2 * 3 * 2 * 2   # kompleks x replika x satir x seri


def test_bos_veri_seti_sadece_baslik_yazar(mdkit, fake_config, tmp_path):
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr
    lines = (tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("complex,replica,analysis,output,series,time_ps")
