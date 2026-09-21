import csv

from conftest import run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && '


def test_log_init_baslik_yazar(lib, fake_config, tmp_path):
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && echo "$MDKIT_LOG"')
    assert r.returncode == 0, r.stderr
    log = tmp_path / "results" / "run_log.csv"
    assert log.exists()
    assert log.read_text().splitlines()[0] == (
        "timestamp,complex,replica,analysis,status,seconds,error"
    )


def test_log_row_satir_ekler(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd OK 12 && '
          'mdkit_log_row last1 rep2 rmsd HATA 3 "gmx patladi"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert len(rows) == 2
    assert rows[0]["status"] == "OK"
    assert rows[1]["error"] == "gmx patladi"


def test_log_row_virgul_ve_yenisatiri_bozmaz(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 '
          '"$(printf \'a, b\\nc\')"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert len(rows) == 1
    assert "a, b" in rows[0]["error"]
    assert "c" in rows[0]["error"]


def test_log_row_cift_tirnagi_standart_kacirir(lib, fake_config, tmp_path):
    """CSV standardi tirnagi ikileyerek kacirir; tek tirnaga cevirmek
    hata mesajini degistirir ve yaniltici olur."""
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 '
          "'gmx: \"can not find group\"'"
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert rows[0]["error"] == 'gmx: "can not find group"'


def test_log_init_mevcut_dosyayi_ezmez(lib, fake_config, tmp_path):
    snippet = PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && mdkit_log_row a b c OK 1'
    run_bash(snippet)
    run_bash(snippet)
    lines = (tmp_path / "results" / "run_log.csv").read_text().splitlines()
    assert lines.count("timestamp,complex,replica,analysis,status,seconds,error") == 1
    assert len(lines) == 3


def test_outputs_present_hepsi_varsa_dogru(lib, fake_config, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.xvg").write_text("veri")
    (out / "b.xvg").write_text("veri")
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_outputs_present "{out}" a.xvg b.xvg'
    )
    assert r.returncode == 0


def test_outputs_present_bos_dosyayi_yok_sayar(lib, fake_config, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.xvg").write_text("veri")
    (out / "b.xvg").write_text("")
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_outputs_present "{out}" a.xvg b.xvg'
    )
    assert r.returncode != 0

    (out / "b.xvg").write_text("veri")
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_outputs_present "{out}" a.xvg b.xvg'
    )
    assert r.returncode == 0


def test_run_isolated_hatayi_yakalar_ve_devam_eder(lib, fake_config, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text('analysis_run() { echo "birinci satir"; echo "patladim" >&2; return 1; }\n')
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_run_isolated "{bad}" /tmp /tmp; echo "rc=$?"; echo "err=$MDKIT_LAST_ERROR"'
    )
    assert "rc=1\n" in r.stdout
    assert "patladim" in r.stdout


def test_run_isolated_basarida_sifir_doner(lib, fake_config, tmp_path):
    good = tmp_path / "good.sh"
    good.write_text('analysis_run() { echo calisti; return 0; }\n')
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_run_isolated "{good}" /tmp /tmp; echo "rc=$?"'
    )
    assert "rc=0\n" in r.stdout
