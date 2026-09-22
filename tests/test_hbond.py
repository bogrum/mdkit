import os
import subprocess
import sys
from pathlib import Path

from conftest import needs_gmx, run_bash


def _list(mdkit, config, ad):
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(config), "--list"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        if parts and parts[0] == ad:
            return parts
    raise AssertionError(f"--list ciktisinda {ad} yok:\n{r.stdout}")


def _kod(path):
    """Yorumlari atilmis kaynak: aciklama satirlari bayrak testini
    dusurmemeli."""
    return "\n".join(l for l in path.read_text().splitlines()
                     if not l.strip().startswith("#"))


def test_manifesto(mdkit, fake_config):
    parts = _list(mdkit, fake_config, "hbond")
    assert parts[1] == "timeseries"
    assert parts[3] == "hbond_pep_mhc.xvg"
    assert parts[4].strip() == ""


def test_sozlesme_dogrulamasindan_gecer(mdkit, lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'mdkit_clear_plugin && source "{mdkit}/analysis/hbond.sh" && '
        f'mdkit_validate_plugin "{mdkit}/analysis/hbond.sh"'
    )
    assert r.returncode == 0, r.stderr


def test_topoloji_kullaniliyor_pdb_degil(mdkit):
    """OLCULDU: gmx hbond -s check_ref.pdb ile SESSIZCE 0 H-bagi dondurur
    (PDB'de bag bilgisi yok, donor/akseptor cikarilamiyor). Ayni komut
    -s md_0_10.tpr ile 10.96 buluyor. Hata da vermiyor -- yani yanlis
    dosyayi vermek sessiz veri kaybidir."""
    src = _kod(mdkit / "analysis" / "hbond.sh")
    assert "TPR_NAME" in src, "topoloji kullanilmali"
    assert "REF_PDB" not in src, "-s icin PDB KULLANILMAMALI"


def test_hedef_receptor(mdkit):
    """OLCULDU: bes komplekste de b2m (AUX) peptide TAM SIFIR H-bagi
    katkisi veriyor. Hedef agir zincir."""
    src = _kod(mdkit / "analysis" / "hbond.sh")
    assert 'group "RECEPTOR"' in src
    assert "AUX" not in src


@needs_gmx
def test_gercek_kosu_sifir_dondurmez(mdkit, real_config):
    """Asil regresyon: PDB tuzagina dusulurse bu test sifir gorur."""
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "hbond", "-b", "0", "--all", root],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    f = next(Path(root).glob("*/rep1/analysis/hbond_pep_mhc.xvg"))
    sys.path.insert(0, str(mdkit))
    import collect_results
    _meta, rows = collect_results.parse_xvg(f)
    assert rows, "bos cikti"
    ort = sum(r[1] for r in rows) / len(rows)
    assert ort > 0, f"ortalama {ort} -- PDB tuzagina dusulmus olabilir"


def test_o_acikca_veriliyor(mdkit):
    """gmx hbond, -o'suz cagrildiginda H-bagi index dosyasini CALISMA
    DIZININE yazar ve her kosuda #hbond.ndx.N# yedegi birakir. GROMACS
    99. yedekte durur:

        Will not make more than 99 backups

    105 replikalik gercek bir kosuda bu, SON 21 REPLIKANIN hata vermesi
    demekti (olculdu). Ayni tuzak spec'te gmx rms icin belgelenmisti.
    """
    src = _kod(mdkit / "analysis" / "hbond.sh")
    assert " -o " in src, "-o acikca verilmeli"
    assert "mktemp" in src, "-o gecici dizine yonlendirilmeli"


@needs_gmx
def test_calisma_dizinine_dosya_sizmaz(mdkit, real_config, tmp_path):
    """Kosu, cagrildigi dizine hicbir sey birakmamali."""
    bos = tmp_path / "bos"
    bos.mkdir()
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "hbond", "-b", "0", "--force", "--all", root],
        capture_output=True, text=True, cwd=str(bos),
    )
    assert r.returncode == 0, r.stderr
    kalan = sorted(p.name for p in bos.iterdir())
    assert kalan == [], f"calisma dizinine sizan dosyalar: {kalan}"
