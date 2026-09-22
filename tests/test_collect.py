import csv
import subprocess
import sys
import textwrap

import numpy as np
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

    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
    pr = list(csv.DictReader((tmp_path / "results" / "profile_long.csv").read_text().splitlines()))
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
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
    assert float(ts[1]["value"]) == pytest.approx(0.1457492)


def test_manifestte_olmayan_xvg_yok_sayilir(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "elle_yazilmis_baska.xvg": TIMESERIES_XVG,
    })
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
    assert {row["output"] for row in ts} == {"rmsd_pep_on_mhc.xvg"}


def test_cok_serili_xvg_her_seriyi_ayri_satira_yazar(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {"rmsd_pep_on_mhc.xvg": MULTI_SERIES_XVG})
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
    series = {row["series"] for row in ts}
    assert series == {"birinci", "ikinci"}
    assert len(ts) == 2 * 3 * 2 * 2   # kompleks x replika x satir x seri


def test_bilinmeyen_kind_uyarilir_ve_csvye_girmez(mdkit, fake_config, fake_dataset, tmp_path):
    """collect(): ANALYSIS_KIND ne 'timeseries' ne 'profile' ise CSV'lere
    sessizce girmemeli; stderr'e cikti adini ve kind'i belirten EN FAZLA BIR
    satir yazilmali -- fake_dataset 2 kompleks x 3 replikadir, hepsi ayni
    (output, kind) cifti icin uyari uretir, dedup bunu teke indirmeli.

    Gecici bir 'matrix' eklentisi kullanilir ki test gercek --list yolunu
    (run_analysis.sh --list -> collect_results.read_manifest) egzersiz
    etsin, sahte/stub bir manifest degil."""
    plugin = mdkit / "analysis" / "zz_test_matrix_kind.sh"
    plugin.write_text(textwrap.dedent("""\
        ANALYSIS_NAME="zzmatrix"
        ANALYSIS_DESC="test icin (matrix kind, henuz desteklenmiyor)"
        ANALYSIS_KIND="matrix"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(matrix_test.xvg)
        analysis_run() { echo "@ dummy" > "$2/matrix_test.xvg"; }
        """))
    try:
        write_outputs(fake_dataset, {"matrix_test.xvg": PROFILE_XVG})
        r = run_collect(mdkit, fake_config)
        assert r.returncode == 0, r.stderr

        warnings = [ln for ln in r.stderr.splitlines() if "matrix_test.xvg" in ln]
        assert len(warnings) == 1, r.stderr
        assert "matrix" in warnings[0]

        ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
        pr = list(csv.DictReader((tmp_path / "results" / "profile_long.csv").read_text().splitlines()))
        assert "matrix_test.xvg" not in {row["output"] for row in ts}
        assert "matrix_test.xvg" not in {row["output"] for row in pr}
    finally:
        plugin.unlink()


def test_bos_veri_seti_sadece_baslik_yazar(mdkit, fake_config, tmp_path):
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr
    lines = (tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("complex,replica,analysis,output,series,time_ps")


# --- Fix 1: kosullu ilan edilen ciktilar ve sessiz atlamalar --------------

def test_groove_fit_ciktisi_manifestoda_ve_csvde_yer_alir(
    mdkit, fake_config, fake_dataset, tmp_path
):
    """rmsf_pep_groovefit.xvg diskte varsa CSV'ye GIRMELI.

    Manifesto `run_analysis.sh --list` ile TAZE bir kabukta okunur; o kabuk
    --groove-fit gormez. Ucuncu cikti kosullu ILAN edildigi surece manifestoya
    hic girmiyor, collect() de `entry is None` dalinda dosyayi sessizce
    atliyordu: saatler suren ve ~24 GB ara dosya ureten bir kosunun 105
    .xvg'sinin tamami hicbir CSV'ye ve hicbir figure'e girmiyordu."""
    write_outputs(fake_dataset, {"rmsf_pep_groovefit.xvg": PROFILE_XVG})
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr

    pr = list(csv.DictReader((tmp_path / "results" / "profile_long.csv").read_text().splitlines()))
    groove = [row for row in pr if row["output"] == "rmsf_pep_groovefit.xvg"]
    assert len(groove) == 2 * 3 * 2          # kompleks x replika x residue
    assert {row["analysis"] for row in groove} == {"rmsf"}
    assert "rmsf_pep_groovefit.xvg" not in r.stderr   # manifestoda, uyari yok


def test_manifestte_olmayan_xvg_tek_satirlik_uyari_uretir(
    mdkit, fake_config, fake_dataset, tmp_path
):
    """Manifestoda olmayan bir .xvg atlanmaya devam eder ama artik SESSIZ
    degil: dosya adi basina TEK bir stderr satiri. fake_dataset 2 kompleks x
    3 replikadir, yani dosya 6 kez gorulur; dedup bunu teke indirmeli.
    Bu uyari, groove-fit hatasini kosmadan once yakalayacak olan uyaridir."""
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "elle_yazilmis_baska.xvg": TIMESERIES_XVG,
    })
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr

    warnings = [ln for ln in r.stderr.splitlines()
                if "elle_yazilmis_baska.xvg" in ln]
    assert len(warnings) == 1, r.stderr

    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()))
    assert {row["output"] for row in ts} == {"rmsd_pep_on_mhc.xvg"}


MATRIX_XPM = textwrap.dedent("""\
    /* XPM */
    /* legend:  "RMSD (nm)" */
    static char *gromacs_xpm[] = {
    "3 3   3 1",
    "A  c #FFFFFF " /* "0" */,
    "B  c #808080 " /* "0.5" */,
    "C  c #000000 " /* "1" */,
    /* x-axis:  0 100 200 */
    /* y-axis:  0 100 200 */
    "CBA",
    "BAB",
    "ABC"
    """)


def test_peer_replica_bilinen_adlarla_eslesir():
    reps = ["rep1", "rep2", "rep3"]
    assert collect_results.peer_replica("cross_rmsd_rep2.xpm", reps) == "rep2"
    assert collect_results.peer_replica("dssp.xpm", reps) is None


def test_matris_npz_olarak_yazilir(fake_dataset, fake_config, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep2.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep2.xpm": ("cross_rmsd", "matrix")}
    mdir = tmp_path / "matrices"
    ts, pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=mdir,
    )
    assert ts == [] and pr == []
    # 2 kompleks x 3 replika
    assert len(mx) == 6
    f = mdir / "last1_rep1_cross_rmsd_rep2.npz"
    assert f.exists(), sorted(p.name for p in mdir.glob("*"))
    d = np.load(f, allow_pickle=False)
    assert d["values"].shape == (3, 3)
    assert np.allclose(np.diag(d["values"]), 0.0)
    assert np.allclose(d["x_ps"], [0.0, 100.0, 200.0])
    assert str(d["unit"]) == "nm"
    assert str(d["replica_i"]) == "rep1"
    assert str(d["replica_j"]) == "rep2"


def test_matris_kind_artik_uyari_uretmiyor(fake_dataset, capsys, tmp_path):
    """Once .xpm hic gorulmuyordu (*.xvg glob'u); simdi toplaniyor."""
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m",
    )
    err = capsys.readouterr().err
    assert "taninmayan ANALYSIS_KIND" not in err


def test_manifestoda_olmayan_xpm_uyarir(fake_dataset, capsys, tmp_path):
    write_outputs(fake_dataset, {"baskabir.xpm": MATRIX_XPM})
    collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], {},
        matrix_dir=tmp_path / "m",
    )
    err = capsys.readouterr().err
    assert "baskabir.xpm" in err


def test_bozuk_xpm_kosuyu_durdurmaz(fake_dataset, capsys, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": "bozuk icerik\n"})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m",
    )
    assert mx == []
    assert "cross_rmsd_rep1.xpm" in capsys.readouterr().err


def test_matrix_dir_yoksa_matris_atlanir(fake_dataset, tmp_path):
    """collect(), matrix_dir verilmediginde .xpm'leri hic islemez."""
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest)
    assert mx == []
