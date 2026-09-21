import subprocess

from conftest import needs_gmx


def run_cli(mdkit, *args):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True,
    )


def data_rows(path):
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


@needs_gmx
def test_iki_profil_uretilir(mdkit, real_config, tmp_path):
    r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
                str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    assert r.returncode == 0, r.stderr
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert (out / "rmsf_pep_self.xvg").exists(), r.stderr[-800:]
    assert (out / "rmsf_mhc.xvg").exists()


@needs_gmx
def test_peptid_profili_residue_sayisi_kadar(mdkit, real_config, tmp_path):
    """last10 peptidi IMGQQPAPQV -> 10 residue."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert len(data_rows(out / "rmsf_pep_self.xvg")) == 10


@needs_gmx
def test_mhc_profili_275_residue(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert len(data_rows(out / "rmsf_mhc.xvg")) == 275


@needs_gmx
def test_groove_fit_ucuncu_cikti_ve_fitli_traj_uretir(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    assert (rep / "analysis" / "rmsf_pep_groovefit.xvg").exists()
    assert (rep / "traj_fit_mhc.xtc").exists()
    assert len(data_rows(rep / "analysis" / "rmsf_pep_groovefit.xvg")) == 10


@needs_gmx
def test_groove_ve_self_farkli_hesaplamadan_gelir(mdkit, real_config, tmp_path):
    """rmsf_pep_self (LIGAND uzerinde fit+hesap) ve rmsf_pep_groovefit
    (RECEPTOR_BB'ye fitli trajektoride -nofit ile LIGAND) farkli fit
    yontemleri kullanir. Ayni sayida satir uretmeliler, ama profiller
    birebir ayni OLMAMALI -- ayniysa oluk-cercevesi yolu yanlislikla
    self-fit hesabina yonlenmis demektir (ör. -nofit unutulmus ya da
    fitli trajektori yerine orijinal $XTC kullanilmis olabilir).

    Kayitli gozlem (last10 fixture, 21 frame, -b 0): self toplam=0.9583,
    groove toplam=1.1121 (~%16 fark). Bu, oluk cercevesindeki RMSF'nin
    kayma/sallanmayi da icerdigi icin daha buyuk olmasi beklentisiyle
    tutarli -- ama residue bazinda bu iliski tersine donebiliyor
    (residue 1: self=0.087 > groove=0.080), yani 'groove >= self'
    toplamda gozlenen ama ispatlanmamis bir alan beklentisi; kod
    seviyesinde garanti edilebilecek tek sey iki profilin gercekten
    FARKLI bir hesaplamadan geldigidir.
    """
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    self_v = [float(ln.split()[1]) for ln in data_rows(out / "rmsf_pep_self.xvg")]
    groove = [float(ln.split()[1]) for ln in data_rows(out / "rmsf_pep_groovefit.xvg")]
    assert len(self_v) == len(groove)
    assert self_v != groove


@needs_gmx
def test_groove_fit_idempotency_ucuncu_ciktiyi_da_sayar(mdkit, real_config, tmp_path):
    """--groove-fit ile ANALYSIS_OUTPUTS uc elemanli olmali; aksi halde
    ucuncu cikti hic uretilmeden SKIP_DONE verilirdi."""
    import csv
    target = str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora")
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0", target)
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit", target)
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert (out / "rmsf_pep_groovefit.xvg").exists()
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert rows[-1]["status"] == "OK"
