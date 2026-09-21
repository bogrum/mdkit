from conftest import run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && '


def test_all_modu_kompleksleri_bulur(lib, fake_config, fake_dataset):
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_find_complexes "{fake_dataset}" 1'
    )
    assert r.returncode == 0, r.stderr
    lines = r.stdout.split()
    assert len(lines) == 2
    assert lines[0].endswith("last1_AAA_A0201_pandora")
    assert lines[1].endswith("top1_BBB_A0201_pandora")
    assert not any("not_a_complex" in ln for ln in lines)


def test_tek_kompleks_modu_hedefi_aynen_dondurur(lib, fake_config, fake_dataset):
    target = fake_dataset / "last1_AAA_A0201_pandora"
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_find_complexes "{target}" 0'
    )
    assert r.stdout.strip() == str(target)


def test_rep_status_tam_replikada_ok(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "OK"


def test_rep_status_eksik_ref_bildirir(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "NO_REF"


def test_rep_status_eksik_traj_bildirir(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep2"
    (rep / "traj_compact_center_dry.xtc").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "NO_TRAJ"


def test_complex_name_kisaltir(lib, fake_config, fake_dataset):
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_complex_name "{cx}"')
    assert r.stdout.strip() == "last1"
