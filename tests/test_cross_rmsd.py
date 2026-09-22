import subprocess

from conftest import run_bash


def _list(mdkit, config):
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(config), "--list"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        if parts and parts[0] == "cross_rmsd":
            return parts
    raise AssertionError(f"--list ciktisinda cross_rmsd yok:\n{r.stdout}")


def test_manifestoda_matrix_kind(mdkit, fake_config):
    parts = _list(mdkit, fake_config)
    assert parts[1] == "matrix"


def test_ciktilar_repsten_turetilir(mdkit, fake_config):
    """fake_config REPS=(rep1 rep2 rep3)."""
    parts = _list(mdkit, fake_config)
    assert parts[3] == ("cross_rmsd_rep1.xpm,cross_rmsd_rep2.xpm,"
                        "cross_rmsd_rep3.xpm")


def test_farkli_reps_farkli_manifesto(mdkit, fake_config, tmp_path):
    """Cikti kumesi config'e bagli; dosyada sabit kodlu degil."""
    alt = tmp_path / "alt_config.sh"
    alt.write_text(
        fake_config.read_text().replace("REPS=(rep1 rep2 rep3)",
                                        "REPS=(repA repB)")
    )
    parts = _list(mdkit, alt)
    assert parts[3] == "cross_rmsd_repA.xpm,cross_rmsd_repB.xpm"


def test_opsiyonel_cikti_yok(mdkit, fake_config):
    """Her replika her matrisi URETIR; kosullu cikti yok, yoksa idempotency
    kontrolu eksik dosya arar ve analiz her kosuda yeniden kosar."""
    parts = _list(mdkit, fake_config)
    assert parts[4].strip() == ""


def test_sozlesme_dogrulamasindan_gecer(mdkit, lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'mdkit_clear_plugin && source "{mdkit}/analysis/cross_rmsd.sh" && '
        f'mdkit_validate_plugin "{mdkit}/analysis/cross_rmsd.sh"'
    )
    assert r.returncode == 0, r.stderr


def test_yardimci_degisken_sizmiyor(mdkit, lib, fake_config):
    """Eklenti ANA KABUGA source edilir; dongu degiskeni orada kalmamali."""
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'source "{mdkit}/analysis/cross_rmsd.sh" && '
        f'echo "SIZDI=${{_cross_rmsd_rep:-yok}}"'
    )
    assert "SIZDI=yok" in r.stdout, r.stdout
