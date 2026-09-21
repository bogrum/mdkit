import subprocess

import pytest


def run_cli(mdkit, *args, **kw):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True, **kw,
    )


def test_help_sifirla_doner(mdkit):
    r = run_cli(mdkit, "--help")
    assert r.returncode == 0
    assert "kullanim:" in r.stdout


def test_bilinmeyen_secenek_hata(mdkit):
    r = run_cli(mdkit, "--boyle-bir-sey-yok")
    assert r.returncode == 2
    assert "bilinmeyen secenek" in r.stderr


def test_hedefsiz_cagri_hata(mdkit, fake_config):
    r = run_cli(mdkit, "--config", str(fake_config))
    assert r.returncode == 2
    assert "hedef" in r.stderr.lower()


def test_list_tsv_uretir(mdkit, fake_config):
    r = run_cli(mdkit, "--config", str(fake_config), "--list")
    assert r.returncode == 0, r.stderr
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert len(lines) >= 2
    names = set()
    for ln in lines:
        parts = ln.split("\t")
        assert len(parts) == 4, f"4 alan bekleniyordu: {ln!r}"
        name, kind, desc, outputs = parts
        assert kind in {"timeseries", "profile", "matrix"}
        assert desc.strip()
        assert outputs.strip()
        names.add(name)
    assert {"rmsd", "rmsf"} <= names


def test_bilinmeyen_analiz_adi_hata(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "-a", "yokboylebiranaliz",
        str(fake_dataset / "last1_AAA_A0201_pandora"),
    )
    assert r.returncode == 2
    assert "bilinmeyen analiz" in r.stderr


def test_dry_run_hicbir_sey_yazmaz(mdkit, fake_config, fake_dataset, tmp_path):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "--all", str(fake_dataset)
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.count("DRY-RUN") == 2 * 3 * 2   # 2 kompleks x 3 replika x 2 analiz
    assert "last1/rep1/rmsd" in r.stdout
    assert not (tmp_path / "results" / "run_log.csv").exists()
    assert not (fake_dataset / "last1_AAA_A0201_pandora" / "rep1" / "analysis").exists()


def test_dry_run_tek_kompleks(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == 3 * 2
    assert "last1" not in r.stdout


def test_reps_secenegi_daraltir(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-r", "rep1",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == 2
    assert "rep2" not in r.stdout


def test_begin_secenegi_dry_run_ciktisinda_gorunur(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-b", "5000", "-a", "rmsf",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=5000)" in r.stdout


def test_begin_verilmezse_analiz_varsayilani_kullanilir(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-a", "rmsf",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=10000)" in r.stdout
    r2 = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-a", "rmsd",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=0)" in r2.stdout
