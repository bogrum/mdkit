import subprocess

from conftest import needs_gmx

OUTPUTS = ["rmsd_pep_on_mhc.xvg", "rmsd_pep_internal.xvg", "rmsd_mhc_bb.xvg"]


def run_cli(mdkit, *args):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True,
    )


def data_rows(path):
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


def header(path, key):
    for ln in path.read_text().splitlines():
        if ln.startswith(f"@") and key in ln:
            return ln
    return ""


@needs_gmx
def test_uc_cikti_uretilir(mdkit, real_config, tmp_path):
    r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
                str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    assert r.returncode == 0, r.stderr
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    for name in OUTPUTS:
        assert (out / name).exists(), f"{name} uretilmedi. stderr={r.stderr[-800:]}"


@needs_gmx
def test_frame_sayisi_dogru(mdkit, real_config, tmp_path):
    """small_rep 0-200 ps, dt=10 ps -> 21 frame."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    for name in OUTPUTS:
        assert len(data_rows(out / name)) == 21, name


@needs_gmx
def test_zaman_ps_cinsinden(mdkit, real_config, tmp_path):
    """-tu kullanilmadigi icin x ekseni ps olmali: son satir 200 civari."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    f = out / "rmsd_pep_on_mhc.xvg"
    last_x = float(data_rows(f)[-1].split()[0])
    assert abs(last_x - 200.0) < 1e-6
    assert "(ps)" in header(f, "xaxis")


@needs_gmx
def test_subtitle_gruplari_belgeler(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert "LIGAND after lsq fit to RECEPTOR_BB" in header(
        out / "rmsd_pep_on_mhc.xvg", "subtitle")
    assert "LIGAND after lsq fit to LIGAND_BB" in header(
        out / "rmsd_pep_internal.xvg", "subtitle")


@needs_gmx
def test_ilk_frame_sifira_yakin(mdkit, real_config, tmp_path):
    """Referans frame 0 oldugu icin ilk RMSD ~0 olmali."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    first = float(data_rows(out / "rmsd_pep_on_mhc.xvg")[0].split()[1])
    assert first < 0.01


@needs_gmx
def test_begin_secenegi_frameleri_kirpar(mdkit, real_config, tmp_path):
    """-b 100 (ps) -> 100..200 ps arasi 11 frame."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd", "-b", "100",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    rows = data_rows(out / "rmsd_pep_on_mhc.xvg")
    assert len(rows) == 11
    assert abs(float(rows[0].split()[0]) - 100.0) < 1e-6


@needs_gmx
def test_ic_ve_toplam_rmsd_farkli_fit_grubundan_gelir(mdkit, real_config, tmp_path):
    """rmsd_pep_on_mhc (fit RECEPTOR_BB) ve rmsd_pep_internal (fit LIGAND_BB)
    farkli fit gruplari kullanir. Ayni sayida frame uretmeliler, ama
    egrileri birebir ayni OLMAMALI -- ayniysa fit grubu hesaplamayi
    etkilemiyor demektir (ör. iki gmx rms cagrisina yanlislikla ayni fit
    grubu verilmis olabilir)."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    total = [float(ln.split()[1]) for ln in data_rows(out / "rmsd_pep_on_mhc.xvg")]
    internal = [float(ln.split()[1]) for ln in data_rows(out / "rmsd_pep_internal.xvg")]
    assert len(total) == len(internal)
    assert total != internal
