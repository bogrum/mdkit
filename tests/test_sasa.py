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
    return "\n".join(l for l in path.read_text().splitlines()
                     if not l.strip().startswith("#"))


def test_manifesto(mdkit, fake_config):
    parts = _list(mdkit, fake_config, "sasa")
    assert parts[1] == "timeseries"
    assert parts[3] == "sasa_pep_in_complex.xvg,sasa_pep_alone.xvg"
    assert parts[4].strip() == ""


def test_sozlesme_dogrulamasindan_gecer(mdkit, lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'mdkit_clear_plugin && source "{mdkit}/analysis/sasa.sh" && '
        f'mdkit_validate_plugin "{mdkit}/analysis/sasa.sh"'
    )
    assert r.returncode == 0, r.stderr


def test_seyreltme_ortam_degiskeniyle_ezilir(mdkit):
    """SASA pahali: tam trajektoride replika basina 3 dk 10 sn, -dt 100 ile
    19.7 sn (10 kat) ve ortalama farki %0.03. Varsayilan seyreltme sart,
    ama ezilebilir olmali."""
    src = _kod(mdkit / "analysis" / "sasa.sh")
    assert "SASA_DT" in src
    assert " -dt " in src


def test_residue_bazli_cikti_yok(mdkit):
    """-or KULLANILMAZ: OLCULDU, iki blogu (yuzey + -output secimi)
    birlestirip 386 satir ve 100 TEKRARLI residue numarasi uretiyor,
    numaralar zincir basina sifirlaniyor. profile kind'i residue'yu
    benzersiz varsayar; bu dosya onu sessizce bozar."""
    src = _kod(mdkit / "analysis" / "sasa.sh")
    assert " -or " not in src
    assert " -oa " not in src


@needs_gmx
def test_gercek_kosu_kompleks_icinde_daha_az_acik(mdkit, real_config):
    """Fiziksel dogrulama: peptid olukta kismen gomulu oldugu icin
    kompleks icindeki acik yuzeyi, tek basinakinden KUCUK olmali. Iki
    ciktinin karistirilmasi da bu testle yakalanir."""
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "sasa", "-b", "0", "--all", root],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    sys.path.insert(0, str(mdkit))
    import collect_results

    def ort(ad):
        f = next(Path(root).glob(f"*/rep1/analysis/{ad}"))
        _m, rows = collect_results.parse_xvg(f)
        assert rows, f"{ad} bos"
        return sum(r[-1] for r in rows) / len(rows)

    icinde = ort("sasa_pep_in_complex.xvg")
    tek = ort("sasa_pep_alone.xvg")
    assert 0 < icinde < tek, f"kompleks icinde {icinde:.2f}, tek basina {tek:.2f}"
