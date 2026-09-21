import subprocess


def test_calisma_kodunda_mutlak_yol_yok(mdkit):
    """config.sh disinda hicbir calisma dosyasinda makineye ozel yol olmamali.
    tests/ haric tutulur: testler bu makinenin gercek verisine bakar ve
    yollari MDKIT_TEST_* ortam degiskenleriyle ezilebilir (bkz. conftest.py)."""
    r = subprocess.run(
        ["grep", "-rEn", "/mnt/|/usr/local/|/home/", str(mdkit),
         "--include=*.sh", "--include=*.py",
         "--exclude=config.sh", "--exclude-dir=tests",
         "--exclude-dir=__pycache__"],
        capture_output=True, text=True,
    )
    assert r.stdout == "", f"mutlak yol bulundu:\n{r.stdout}"


def test_proje_adi_kodda_gecmiyor(mdkit):
    """'pandora', 'TUSEB', 'A0201' gibi projeye ozel isimler config'de kalmali."""
    r = subprocess.run(
        ["grep", "-rEni", "pandora|TUSEB|A0201", str(mdkit),
         "--include=*.sh", "--include=*.py",
         "--exclude=config.sh", "--exclude-dir=tests",
         "--exclude-dir=__pycache__"],
        capture_output=True, text=True,
    )
    assert r.stdout == "", f"projeye ozel isim bulundu:\n{r.stdout}"


def test_kopya_configle_calisir(mdkit, fake_config, fake_dataset):
    """Aracin config disinda hicbir seye bagli olmadiginin uctan uca kaniti."""
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(fake_config),
         "--dry-run", "--all", str(fake_dataset)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "DRY-RUN" in r.stdout
