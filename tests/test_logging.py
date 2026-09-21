from pathlib import Path

import csv

from conftest import run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && '


def test_log_init_baslik_yazar(lib, fake_config, tmp_path):
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && echo "$MDKIT_LOG"')
    assert r.returncode == 0, r.stderr
    log = tmp_path / "results" / "run_log.csv"
    assert log.exists()
    assert log.read_text().splitlines()[0] == (
        "timestamp,complex,replica,analysis,status,seconds,begin_ps,error"
    )


def test_log_row_satir_ekler(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd OK 12 10000 && '
          'mdkit_log_row last1 rep2 rmsd HATA 3 0 "gmx patladi"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").read_text().splitlines()))
    assert len(rows) == 2
    assert rows[0]["status"] == "OK"
    assert rows[0]["begin_ps"] == "10000"
    assert rows[1]["begin_ps"] == "0"
    assert rows[1]["error"] == "gmx patladi"


def test_log_row_virgul_ve_yenisatiri_bozmaz(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 0 '
          '"$(printf \'a, b\\nc\')"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").read_text().splitlines()))
    assert len(rows) == 1
    assert "a, b" in rows[0]["error"]
    assert "c" in rows[0]["error"]


def test_log_row_cift_tirnagi_standart_kacirir(lib, fake_config, tmp_path):
    """CSV standardi tirnagi ikileyerek kacirir; tek tirnaga cevirmek
    hata mesajini degistirir ve yaniltici olur."""
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 0 '
          "'gmx: \"can not find group\"'"
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").read_text().splitlines()))
    assert rows[0]["error"] == 'gmx: "can not find group"'


def test_log_init_mevcut_dosyayi_ezmez(lib, fake_config, tmp_path):
    snippet = PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && mdkit_log_row a b c OK 1 0'
    run_bash(snippet)
    run_bash(snippet)
    lines = (tmp_path / "results" / "run_log.csv").read_text().splitlines()
    assert lines.count(
        "timestamp,complex,replica,analysis,status,seconds,begin_ps,error") == 1
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


def test_run_isolated_sozdizimi_hatasini_da_yakalar(lib, fake_config, tmp_path):
    """2>&1 yalnizca analysis_run'a baglandiginda, sozdizimi hatali bir
    eklentinin mesaji `source`tan terminale siziyor ve HATA satirinin error
    kolonu BOS kaliyordu -- teshis edilemez bir basarisizlik. Yonlendirme
    bilesik komuta baglanmali."""
    bad = tmp_path / "bozuk.sh"
    bad.write_text('analysis_run() {\n  echo "kapanis parantezi yok"\n')
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_run_isolated "{bad}" /tmp /tmp; echo "rc=$?"; '
          'echo "err=[$MDKIT_LAST_ERROR]"'
    )
    assert "rc=1\n" in r.stdout
    assert "err=[]" not in r.stdout, r.stdout
    assert r.stderr == "", r.stderr


def test_analysis_scripts_alt_cizgili_dosyalari_atlar(lib, fake_config, mdkit):
    """Yarida kesilmis bir testin geride biraktigi stub, varsayilan (tum
    analizler) kosusunda gercek bir analiz gibi calistirilmamali."""
    stub = mdkit / "analysis" / "_zztest_gizli.sh"
    stub.write_text('ANALYSIS_NAME="gizli"\n')
    try:
        r = run_bash(
            PRE.format(lib=lib, cfg=fake_config) + "mdkit_analysis_scripts"
        )
        assert r.returncode == 0, r.stderr
        listed = [Path(ln).name for ln in r.stdout.split()]
        assert "_zztest_gizli.sh" not in listed, listed
        assert "rmsd.sh" in listed and "rmsf.sh" in listed
        assert "lib.sh" not in listed
    finally:
        stub.unlink()


def test_mandatory_outputs_opsiyonelleri_duser(lib, fake_config):
    """ANALYSIS_OUTPUTS eksi ANALYSIS_OPTIONAL_OUTPUTS. Opsiyonel cikti
    zorunlu sayilirsa, onu hic uretmeyen bir analiz her kosuda bastan
    hesaplanir; hic dusulmezse (eski hal) manifestodan da duser."""
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'ANALYSIS_OUTPUTS=(a.xvg b.xvg c.xvg); '
          'ANALYSIS_OPTIONAL_OUTPUTS=(c.xvg); mdkit_mandatory_outputs'
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["a.xvg", "b.xvg"]


def test_mandatory_outputs_opsiyonel_tanimsizken_hepsini_dondurur(lib, fake_config):
    """ANALYSIS_OPTIONAL_OUTPUTS tanimlamayan eski/basit eklentiler
    `set -u` altinda patlamamali."""
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'ANALYSIS_OUTPUTS=(a.xvg b.xvg); mdkit_mandatory_outputs'
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["a.xvg", "b.xvg"]
