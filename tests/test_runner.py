import csv
import subprocess
import textwrap

from conftest import needs_gmx


def run_cli(mdkit, *args, **kw):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True, **kw,
    )


def _stub_analysis(mdkit, name, body):
    """Gecici bir analiz eklentisi yazar; test sonunda silinir."""
    p = mdkit / "analysis" / f"{name}.sh"
    p.write_text(body)
    return p


def test_eksik_onkosul_skip_missing_loglanir(mdkit, fake_config, fake_dataset, tmp_path):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "traj_compact_center_dry.xtc").unlink()
    run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "rmsd",
            str(fake_dataset / "last1_AAA_A0201_pandora"))
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    skipped = [r for r in rows if r["status"] == "SKIP_MISSING"]
    assert len(skipped) == 1
    assert skipped[0]["replica"] == "rep1"
    assert "NO_TRAJ" in skipped[0]["error"]


def test_bir_replikanin_hatasi_digerlerini_durdurmaz(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() {
            if [[ "$1" == *rep2* ]]; then echo "kasitli hata" >&2; return 1; fi
            echo "veri" > "$2/zz.xvg"
        }
        """))
    try:
        r = run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "zzstub",
                    str(fake_dataset / "last1_AAA_A0201_pandora"))
        assert r.returncode == 0, r.stderr
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        by_rep = {r["replica"]: r["status"] for r in rows}
        assert by_rep == {"rep1": "OK", "rep2": "HATA", "rep3": "OK"}
        assert "kasitli hata" in [r["error"] for r in rows if r["status"] == "HATA"][0]
    finally:
        stub.unlink()


def test_idempotency_ikinci_kosuda_skip_done(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() { echo "veri" > "$2/zz.xvg"; }
        """))
    try:
        args = ["-c", str(fake_config), "-y", "-a", "zzstub",
                str(fake_dataset / "last1_AAA_A0201_pandora")]
        run_cli(mdkit, *args)
        run_cli(mdkit, *args)
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert [r["status"] for r in rows[:3]] == ["OK", "OK", "OK"]
        assert [r["status"] for r in rows[3:]] == ["SKIP_DONE"] * 3
    finally:
        stub.unlink()


def test_force_idempotencyi_ezer(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() { echo "veri" > "$2/zz.xvg"; }
        """))
    try:
        args = ["-c", str(fake_config), "-y", "-a", "zzstub",
                str(fake_dataset / "last1_AAA_A0201_pandora")]
        run_cli(mdkit, *args)
        run_cli(mdkit, *args, "--force")
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert [r["status"] for r in rows[3:]] == ["OK"] * 3
    finally:
        stub.unlink()


def test_analiz_scripti_ortam_degiskenlerini_gorur(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzenv", textwrap.dedent("""\
        ANALYSIS_NAME="zzenv"
        ANALYSIS_DESC="ortam testi"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(env.xvg)
        analysis_run() {
            printf 'REF=%s\\nXTC=%s\\nB=%s\\nGF=%s\\n' \\
                "$REF_PDB" "$XTC" "$B_PS" "$GROOVE_FIT" > "$2/env.xvg"
        }
        """))
    try:
        run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "zzenv", "-r", "rep1",
                "-b", "777", "--groove-fit",
                str(fake_dataset / "last1_AAA_A0201_pandora"))
        out = (fake_dataset / "last1_AAA_A0201_pandora" / "rep1" / "analysis" / "env.xvg").read_text()
        assert "check_ref.pdb" in out
        assert "traj_compact_center_dry.xtc" in out
        assert "B=777" in out
        assert "GF=1" in out
    finally:
        stub.unlink()


def test_eksik_ref_onay_reddedilirse_uretilmez(mdkit, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_cli(mdkit, "-c", str(fake_config), "-a", "rmsd",
                str(fake_dataset / "last1_AAA_A0201_pandora"), input="h\n")
    assert "1 replikada" in r.stdout
    assert not (rep / "check_ref.pdb").exists()


@needs_gmx
def test_eksik_ref_onaylanirsa_uretilir_ve_analiz_devam_eder(mdkit, real_config, tmp_path):
    """Self-heal + devam: rmsd.sh henuz Task 7 stub'u oldugu icin gercek analiz
    yerine yerel bir stub eklenti kullanilir -- Task 6'nin sahip oldugu sey
    self-heal onayi ve dongunun devami, rmsd'nin analiz govdesi degil.
    mdkit_make_ref gercek gmx cagirir, bu yuzden @needs_gmx/real_config kalir.
    """
    stub = _stub_analysis(mdkit, "zzheal", textwrap.dedent("""\
        ANALYSIS_NAME="zzheal"
        ANALYSIS_DESC="self-heal sonrasi devam testi"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(heal.xvg)
        analysis_run() { echo "veri" > "$2/heal.xvg"; }
        """))
    try:
        rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
        (rep / "check_ref.pdb").unlink()
        r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "zzheal",
                    str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
        assert (rep / "check_ref.pdb").exists(), r.stderr
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert any(row["status"] == "OK" for row in rows), rows
    finally:
        stub.unlink()


@needs_gmx
def test_force_indeksi_de_yeniden_kurar(mdkit, real_config, tmp_path):
    """--force sadece analiz ciktilarini degil, index.ndx'i de yeniden kurmali.

    Aksi halde config.sh'de zincir haritasi degistirilip yeniden kosulunca
    (ör. CHAIN_LIGAND), --force verilse bile eski/gecersiz index sessizce
    kullanilmaya devam eder.
    """
    target = str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora")
    r1 = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd", target)
    ndx = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis" / "index.ndx"
    assert ndx.exists(), r1.stderr

    ndx.write_text("bozuk ve gecersiz icerik\n")

    r2 = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd", "--force", target)
    content = ndx.read_text()
    assert "bozuk ve gecersiz icerik" not in content, r2.stderr
    for g in ["RECEPTOR", "AUX", "LIGAND", "RECEPTOR_BB", "LIGAND_BB"]:
        assert f"[ {g} ]" in content
