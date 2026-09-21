from conftest import run_bash


def test_gecerli_config_yuklenir(lib, fake_config, fake_dataset):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && echo "$DATA_ROOT"'
    )
    assert r.returncode == 0, r.stderr
    assert str(fake_dataset) in r.stdout


def test_eksik_degisken_adiyla_raporlanir(lib, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text('DATA_ROOT="/x"\nREPS=(rep1)\n')
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "TRAJ_NAME" in r.stderr
    assert "REF_NAME" in r.stderr


def test_eksik_python_hata_verir(lib, tmp_path):
    """PYTHON zorunlu: PATH fallback'i yok, cunku PATH'teki python3 bilimsel
    yigina sahip olmayabilir. Erken yakalanmali."""
    bad = tmp_path / "bad.sh"
    bad.write_text(
        'DATA_ROOT=/x\nCOMPLEX_GLOB="*"\nTRAJ_NAME=a\nTPR_NAME=b\nREF_NAME=c\n'
        'RESULTS_DIR=/y\nCHAIN_RECEPTOR=A\nCHAIN_AUX=B\nCHAIN_LIGAND=C\nREPS=(rep1)\n'
    )
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "PYTHON" in r.stderr


def test_bos_reps_hata_verir(lib, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text(
        'DATA_ROOT=/x\nCOMPLEX_GLOB="*"\nTRAJ_NAME=a\nTPR_NAME=b\nREF_NAME=c\n'
        'RESULTS_DIR=/y\nPYTHON=/usr/bin/python3\n'
        'CHAIN_RECEPTOR=A\nCHAIN_AUX=B\nCHAIN_LIGAND=C\nREPS=()\n'
    )
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "REPS" in r.stderr


def test_olmayan_config_anlasilir_hata(lib, tmp_path):
    r = run_bash(f'source "{lib}" && mdkit_load_config "{tmp_path}/yok.sh"')
    assert r.returncode != 0
    assert "bulunamadi" in r.stderr


def test_varsayilan_config_script_yaninda_aranir(lib, mdkit):
    r = run_bash(f'source "{lib}" && mdkit_load_config && echo "$TRAJ_NAME"')
    assert r.returncode == 0, r.stderr
    assert "traj_compact_center_dry.xtc" in r.stdout
