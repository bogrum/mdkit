import subprocess

import pytest


def run_cli(mdkit, *args, **kw):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True, **kw,
    )


def n_analyses(mdkit, config):
    """--list'teki analiz sayisi. Sabit kodlanmaz: analysis/ altina bir
    eklenti eklemek dry-run SAYIM testlerini dusurmemeli -- bu testler
    kompleks x replika yayilimini olcuyor, eklenti sayisini degil."""
    r = run_cli(mdkit, "-c", str(config), "--list")
    assert r.returncode == 0, r.stderr
    return len([ln for ln in r.stdout.splitlines() if ln.strip()])


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
    """Ilk DORT alan konumunu ve anlamini korur (makine okunur manifesto);
    5. alan opsiyonel ciktilari tasir ve 4. alanin alt kumesidir."""
    r = run_cli(mdkit, "--config", str(fake_config), "--list")
    assert r.returncode == 0, r.stderr
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert len(lines) >= 2
    names = set()
    for ln in lines:
        parts = ln.split("\t")
        assert len(parts) == 5, f"5 alan bekleniyordu: {ln!r}"
        name, kind, desc, outputs, optional = parts
        assert kind in {"timeseries", "profile", "matrix"}
        assert desc.strip()
        assert outputs.strip()
        opt = {o for o in optional.split(",") if o}
        assert opt <= {o for o in outputs.split(",") if o}, ln
        names.add(name)
    assert {"rmsd", "rmsf"} <= names


def test_list_opsiyonel_ciktilari_bes_alanda_bildirir(mdkit, fake_config):
    """--list TAZE bir kabukta kosar, yani --groove-fit gormez. Ucuncu rmsf
    ciktisi yine de 4. alanda ILAN edilmeli (yoksa collect_results.py diskteki
    dosyayi manifestoda bulamaz ve sessizce atlar) ve 5. alanda opsiyonel
    olarak isaretlenmeli (yoksa idempotency kontrolu onu arar)."""
    r = run_cli(mdkit, "--config", str(fake_config), "--list")
    assert r.returncode == 0, r.stderr
    row = [ln for ln in r.stdout.splitlines() if ln.startswith("rmsf\t")][0]
    _name, _kind, _desc, outputs, optional = row.split("\t")
    assert "rmsf_pep_groovefit.xvg" in outputs.split(",")
    assert optional.split(",") == ["rmsf_pep_groovefit.xvg"]


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
    # 2 kompleks x 3 replika x analiz sayisi
    assert r.stdout.count("DRY-RUN") == 2 * 3 * n_analyses(mdkit, fake_config)
    assert "last1/rep1/rmsd" in r.stdout
    assert not (tmp_path / "results" / "run_log.csv").exists()
    assert not (fake_dataset / "last1_AAA_A0201_pandora" / "rep1" / "analysis").exists()


def test_dry_run_tek_kompleks(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == 3 * n_analyses(mdkit, fake_config)
    assert "last1" not in r.stdout


def test_reps_secenegi_daraltir(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-r", "rep1",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == n_analyses(mdkit, fake_config)
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


@pytest.mark.parametrize("flag", ["--config", "-b"])
def test_deger_gerektiren_secenek_degersiz_hata(mdkit, flag):
    # Regresyon: -c/-a/-r/-b degersiz birakilirsa shift 2 basarisiz olur ve
    # set -e olmadigi icin dongu ayni $1'i sonsuza kadar isler. timeout=10
    # olmadan bu test, gelecekte bu koruma bozulursa, TAKILIR (fail degil,
    # hic bitmez) -- suit'in geri kalanini da kilitler. timeout burada
    # "guzel olsun" degil, bu testin kendisinin asla asilamamasi icin sart.
    r = run_cli(mdkit, flag, timeout=10)
    assert r.returncode == 2
    assert "deger gerektirir" in r.stderr


def test_iki_pozisyonel_hedef_hata(mdkit, fake_config, fake_dataset, tmp_path):
    other = tmp_path / "baska_dizin"
    other.mkdir()
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run",
        str(fake_dataset / "top1_BBB_A0201_pandora"), str(other),
    )
    assert r.returncode == 2
    assert "hedef" in r.stderr.lower()


def test_bozuk_plugin_list_basarisiz_olur(mdkit, fake_config):
    # Alt cizgi ONEKI YOK: mdkit_analysis_scripts artik '_' ile baslayan
    # dosyalari kesiften duser, oysa bu testin amaci kesfedilen bozuk bir
    # eklentinin --list'i dusurmesini dogrulamak.
    plugin = mdkit / "analysis" / "zz_test_bozuk_plugin.sh"
    plugin.write_text(
        "#!/usr/bin/env bash\n"
        'ANALYSIS_NAME="bozuk"\n'
        'ANALYSIS_DESC="ANALYSIS_KIND kasten eksik"\n'
        "ANALYSIS_NEEDS_INDEX=1\n"
        "ANALYSIS_DEFAULT_BEGIN=0\n"
        "ANALYSIS_OUTPUTS=(x.xvg)\n"
        'analysis_run() { return 1; }\n'
    )
    try:
        r = run_cli(mdkit, "--config", str(fake_config), "--list", timeout=10)
        assert r.returncode != 0
        assert "zz_test_bozuk_plugin.sh" in r.stderr
    finally:
        plugin.unlink()


def test_full_postmd_postmd_script_yoksa_hicbir_is_yapmadan_durur(
    mdkit, fake_config, fake_dataset, tmp_path
):
    """POSTMD_SCRIPT dogrulamasi eksik-referans dongusunun ICINDE duruyordu:
    hata ancak log dosyasi olusturulduktan ve ilk replikaya gelindikten sonra
    veriliyordu. Dogrulama tum isten ONCE yapilmali."""
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    (cx / "rep1" / "check_ref.pdb").unlink()     # missing_refs bos olmasin
    r = run_cli(mdkit, "-c", str(fake_config), "-y", "--full-postmd",
                "-a", "rmsd", str(cx), timeout=30)
    assert r.returncode == 2
    assert "POSTMD_SCRIPT" in r.stderr
    assert not (tmp_path / "results" / "run_log.csv").exists(), r.stdout
