import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from conftest import needs_gmx, run_bash


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
    """fake_config REPS=(rep1 rep2 rep3).

    Es basina IKI cikti: .xpm (eksenler, birim, legend) ve .dat (tam
    hassasiyetli degerler). Ikisi de zorunlu -- .dat'in metadata'si,
    .xpm'in de tam degerleri yoktur, yani biri eksikken matris eksiktir."""
    parts = _list(mdkit, fake_config)
    assert parts[3] == (
        "cross_rmsd_rep1.xpm,cross_rmsd_rep1.dat,"
        "cross_rmsd_rep2.xpm,cross_rmsd_rep2.dat,"
        "cross_rmsd_rep3.xpm,cross_rmsd_rep3.dat")


def test_farkli_reps_farkli_manifesto(mdkit, fake_config, tmp_path):
    """Cikti kumesi config'e bagli; dosyada sabit kodlu degil."""
    alt = tmp_path / "alt_config.sh"
    alt.write_text(
        fake_config.read_text().replace("REPS=(rep1 rep2 rep3)",
                                        "REPS=(repA repB)")
    )
    parts = _list(mdkit, alt)
    assert parts[3] == ("cross_rmsd_repA.xpm,cross_rmsd_repA.dat,"
                        "cross_rmsd_repB.xpm,cross_rmsd_repB.dat")


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


CONTRACT = """
source "{lib}"
mdkit_load_config "{config}"
GMX="{gmx}"
REF_PDB="{rep}/check_ref.pdb"
XTC="{rep}/traj_compact_center_dry.xtc"
NDX="{rep}/analysis/index.ndx"
B_PS=0
FORCE=0
mdkit_clear_plugin
source "{mdkit}/analysis/cross_rmsd.sh"
analysis_run "{rep}" "{out}"
"""


def _code_lines(path):
    """Yorumlari atilmis kaynak.

    Bayrak testleri KODU olcmeli. Eklentinin icinde neden -skip/-tu
    kullanilmadigini anlatan aciklama satirlari var; ham metinde arama
    yapmak bu aciklamalarin kendisini ihlal sayardi.
    """
    return "\n".join(
        line for line in path.read_text().splitlines()
        if not line.strip().startswith("#")
    )


def test_uyusmayan_es_sistemi_reddedilir(mdkit, lib, fake_config,
                                         fake_dataset, tmp_path):
    """Capraz rms TEK index'i IKI trajektoriye birden uygular; esin sistemi
    farkliysa matris sessizce cop olur. gmx yolu bilerek gecersiz: koruma
    gmx cagrilmadan ONCE donmeli."""
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    (cx / "rep1" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    (cx / "rep2" / "check_ref.pdb").write_text("ATOM      1  N\n" * 4)
    (cx / "rep3" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    out = tmp_path / "out"
    out.mkdir()
    r = run_bash(CONTRACT.format(
        lib=lib, config=fake_config, gmx="/nonexistent/gmx",
        rep=cx / "rep1", mdkit=mdkit, out=out,
    ))
    assert r.returncode != 0
    assert "uyusmuyor" in (r.stdout + r.stderr).lower()


def test_eksik_es_trajektorisi_reddedilir(mdkit, lib, fake_config,
                                          fake_dataset, tmp_path):
    cx = fake_dataset / "top1_BBB_A0201_pandora"
    (cx / "rep2" / "traj_compact_center_dry.xtc").unlink()
    out = tmp_path / "out2"
    out.mkdir()
    r = run_bash(CONTRACT.format(
        lib=lib, config=fake_config, gmx="/nonexistent/gmx",
        rep=cx / "rep1", mdkit=mdkit, out=out,
    ))
    assert r.returncode != 0
    assert "trajektori" in (r.stdout + r.stderr).lower()


def test_koruma_gmx_cagrilmadan_once_doner(mdkit, lib, fake_config,
                                           fake_dataset, tmp_path):
    """Dogrulama hesap dongusunun ICINDE olsaydi self-matris (ilk es) once
    uretilir, bozuk ikinci es ancak ondan SONRA fark edilirdi: out_dir'de
    yarim bir cikti kumesi kalirdi."""
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    (cx / "rep1" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    (cx / "rep2" / "check_ref.pdb").write_text("ATOM      1  N\n" * 4)
    (cx / "rep3" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    out = tmp_path / "out3"
    out.mkdir()
    r = run_bash(CONTRACT.format(
        lib=lib, config=fake_config, gmx="/nonexistent/gmx",
        rep=cx / "rep1", mdkit=mdkit, out=out,
    ))
    assert r.returncode != 0
    assert list(out.iterdir()) == [], "hicbir matris yazilmamaliydi"
    assert "gmx" not in r.stderr.lower() or "uyusmuyor" in r.stderr.lower()


def test_skip_ve_tu_kullanilmiyor(mdkit):
    """spec 2.2: -skip .xpm eksen zamanlarini bozar, -tu -b/-e'yi cevirir."""
    src = _code_lines(mdkit / "analysis" / "cross_rmsd.sh")
    assert " -skip" not in src
    assert " -tu" not in src
    assert " -dt " in src


@needs_gmx
def test_gercek_capraz_matris(mdkit, real_config):
    """real_config tek replikali (REPS=(rep1)); self-matris uretilmeli.

    small_rep 0-200 ps'lik 21 frame'dir, bu yuzden -b 0 ve CROSS_RMSD_DT=50
    verilir: eklentinin varsayilani (ANALYSIS_DEFAULT_BEGIN=10000) bu kisa
    trajektoride bos cikti uretirdi."""
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "cross_rmsd", "-b", "0", "--all", root],
        capture_output=True, text=True,
        env={**os.environ, "CROSS_RMSD_DT": "50"},
    )
    assert r.returncode == 0, r.stderr
    xpm = next(Path(root).glob("*/rep1/analysis/cross_rmsd_rep1.xpm"), None)
    assert xpm is not None, r.stdout

    sys.path.insert(0, str(mdkit))
    import collect_results
    meta, values, x, y = collect_results.parse_xpm(xpm)
    assert values.shape[0] == values.shape[1] == len(x) == len(y)
    assert meta["unit"] == "nm"
    # Self-matriste kosegen TANIM GEREGI sifir. Ayristiricinin satir cevirmesi
    # ile gmx'in GERCEK yazma sirasinin uyustugunun kaniti -- fixture'lar bunu
    # gosteremez, cunku fixture'i da biz yaziyoruz.
    assert np.allclose(np.diag(values), 0.0, atol=1e-3)


@needs_gmx
def test_gercek_capraz_dal_transpoze_tutarli(mdkit, pair_config):
    """-f2 dalinin GERCEK regresyonu: rep1'in rep2-matrisi, rep2'nin
    rep1-matrisinin TRANSPOZESI olmalidir.

    Bu ozdeslik ancak (a) piksel satirlari y ekseninin tersi sirada okunursa
    VE (b) genislik=x=-f, yukseklik=y=-f2 yorumu dogruysa saglanir. Ikisinden
    biri yanlis olsa da self-matris hala simetrik ve kosegeni sifir cikardi --
    yani self-matris bu iki hatayi GOREMEZ.
    """
    cfg, root, reps = pair_config
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(cfg),
         "-a", "cross_rmsd", "-b", "0", "--all", str(root)],
        capture_output=True, text=True,
        env={**os.environ, "CROSS_RMSD_DT": "50"},
    )
    assert r.returncode == 0, r.stderr

    sys.path.insert(0, str(mdkit))
    import collect_results

    cx = next(Path(root).glob("*_pandora"))
    a, b = reps[0], reps[1]
    _m1, c_ab, _x, _y = collect_results.parse_xpm(
        cx / a / "analysis" / f"cross_rmsd_{b}.xpm")   # -f a, -f2 b
    _m2, c_ba, _x, _y = collect_results.parse_xpm(
        cx / b / "analysis" / f"cross_rmsd_{a}.xpm")   # -f b, -f2 a
    assert np.allclose(c_ab, c_ba.T, atol=1e-6), (
        f"max fark {np.abs(c_ab - c_ba.T).max()}")

    # Capraz matrisin kosegeni self-matrisin aksine SIFIR DEGILDIR: ayni anda
    # farkli replikalar farkli konformasyonlardadir.
    _m3, self_a, _x, _y = collect_results.parse_xpm(
        cx / a / "analysis" / f"cross_rmsd_{a}.xpm")
    assert np.allclose(np.diag(self_a), 0.0, atol=1e-3)
    assert not np.allclose(np.diag(c_ab), 0.0, atol=1e-3)


@needs_gmx
def test_gercek_dat_tam_hassasiyetli(mdkit, real_config):
    """.dat gercekten .xpm'den daha hassas mi? Iddia edilmiyor, olculuyor.

    Ayrica self-matrisin kosegeni .dat'ta TAM sifir olmalidir: .xpm'de
    sifir, 80 seviyenin en dusugune denk geldigi icin zaten sifir cikar --
    yani kosegen tek basina hassasiyeti kanitlamaz, farkli deger SAYISI
    kanitlar.
    """
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "cross_rmsd", "-b", "0", "--all", root],
        capture_output=True, text=True,
        env={**os.environ, "CROSS_RMSD_DT": "50"},
    )
    assert r.returncode == 0, r.stderr

    adir = next(Path(root).glob("*/rep1/analysis"))
    xpm, dat = adir / "cross_rmsd_rep1.xpm", adir / "cross_rmsd_rep1.dat"
    assert dat.is_file(), sorted(f.name for f in adir.iterdir())

    sys.path.insert(0, str(mdkit))
    import collect_results
    _meta, v_xpm, x, y = collect_results.parse_xpm(xpm)
    v_dat = collect_results.parse_bin(dat, len(x), len(y))

    assert v_dat.shape == v_xpm.shape
    assert np.allclose(v_dat, v_xpm, atol=0.05), "ayni matris olmali"
    assert len(np.unique(v_dat)) > len(np.unique(v_xpm)), (
        f".dat {len(np.unique(v_dat))} farkli deger, "
        f".xpm {len(np.unique(v_xpm))} -- hassasiyet kazanci yok")
    assert np.all(np.diag(v_dat) == 0.0)
