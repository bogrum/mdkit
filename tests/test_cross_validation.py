"""mdkit'in gmx cagri katmaninin capraz dogrulamasi.

DIKKAT: bu dosyadaki testlerin NEYI kanitladigi birbirinden farklidir ve
karistirilmamalidir -- biri mdkit'i, digeri yalnizca bu makinenin
GROMACS'ini dogrular. Docstring'ler bunu acikca soyler.
"""
import csv
import os
import re
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

from conftest import GMX_BIN, REAL_REP, REAL_ROOT, needs_gmx

KNOWN_CSV = REAL_ROOT / "rmsf_per_position.csv"

pytestmark = pytest.mark.slow


def _parse_xvg_values(path):
    return [float(ln.split()[1]) for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


def _ndx_group_number(ndx, name):
    """index.ndx icindeki grubun 0 tabanli SIRA numarasi.

    Testin mdkit'ten BAGIMSIZ oldugu yer burasi: mdkit gruplari kanonik
    ADIYLA secer, bu test ayni grubu kendi ayristirdigi NUMARAYLA secer.
    Ikisi ayni sayilari veriyorsa 'LIGAND' adi gercekten bekledigimiz
    atom kumesine cozulmus demektir."""
    names = re.findall(r"^\[ (.*) \]$", ndx.read_text(), re.M)
    return names.index(name)


def _ndx_group_size(ndx, name):
    size, cur = 0, None
    for line in ndx.read_text().splitlines():
        m = re.match(r"^\[ (.*) \]$", line)
        if m:
            cur = m.group(1)
        elif cur == name:
            size += len(line.split())
    return size


@pytest.fixture
def full_traj_config(tmp_path):
    """TAM trajektoriyi tek replikali bir veri agacinda gosteren config.

    Gercek veri koku YAZILMAZ: girdiler symlink, analysis/ ve results/
    tmp_path altindadir."""
    if not (GMX_BIN.exists() and REAL_REP.exists()):
        pytest.skip("gmx veya gercek veri yok")
    root = tmp_path / "data"
    rep = root / "xval_PEPTIDE_A0201_pandora" / "rep1"
    rep.mkdir(parents=True)
    for name in ("traj_compact_center_dry.xtc", "check_ref.pdb", "md_0_10.tpr"):
        os.symlink(REAL_REP / name, rep / name)
    cfg = tmp_path / "config.sh"
    cfg.write_text(textwrap.dedent(f"""\
        DATA_ROOT="{root}"
        GMX="{GMX_BIN}"
        COMPLEX_GLOB="*_pandora"
        REPS=(rep1)
        TRAJ_NAME="traj_compact_center_dry.xtc"
        TPR_NAME="md_0_10.tpr"
        REF_NAME="check_ref.pdb"
        RESULTS_DIR="{tmp_path}/results"
        PYTHON="{sys.executable}"
        CHAIN_RECEPTOR="A"
        CHAIN_AUX="B"
        CHAIN_LIGAND="C"
        """))
    return cfg


def _standalone_rmsf(ndx, rep, out, group_number, begin):
    subprocess.run(
        [str(GMX_BIN), "rmsf",
         "-s", str(rep / "check_ref.pdb"),
         "-f", str(rep / "traj_compact_center_dry.xtc"),
         "-n", str(ndx), "-o", str(out), "-res", "-b", str(begin)],
        input=f"{group_number}\n", text=True, capture_output=True, check=True,
    )
    return _parse_xvg_values(out)


@needs_gmx
def test_mdkit_rmsf_ciktisi_bagimsiz_gmx_cagrisiyla_eslesir(
    mdkit, full_traj_config, tmp_path
):
    """mdkit'in URETTIGI rmsf_pep_self.xvg, aynisini elle hesaplayan bagimsiz
    bir gmx cagrisiyla eslesmeli.

    Bu test mdkit'i GERCEKTEN kosar: run_analysis.sh, lib.sh, index kurucusu,
    kanonik grup adlari ve rmsf.sh eklentisi devrededir. Karsilastirma
    grubunu mdkit'in yazdigi index.ndx'ten NUMARAYLA secer, yani 'LIGAND'
    adinin dogru atom kumesine cozuldugunu de dogrular. Son olarak -b'nin
    gercekten gmx'e ulastigini, -b 0 ile hesaplanan profilin FARKLI
    cikmasiyla gosterir.

    (Onceki surumu bu dosyanin tamamini mdkit'siz kosuyordu: mdkit/ silinse
    bile gecerdi ve yalnizca GROMACS'in deterministik oldugunu olcuyordu.)"""
    target = tmp_path / "data" / "xval_PEPTIDE_A0201_pandora"
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(full_traj_config),
         "-y", "-a", "rmsf", str(target)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rep = target / "rep1"
    out_dir = rep / "analysis"
    ndx = out_dir / "index.ndx"
    assert ndx.exists(), r.stderr[-800:]

    mdkit_profile = _parse_xvg_values(out_dir / "rmsf_pep_self.xvg")
    # Index kurulumu: LIGAND gercekten peptiddir (10 residue / 151 atom).
    assert len(mdkit_profile) == 10
    assert _ndx_group_size(ndx, "LIGAND") == 151

    ligand_no = _ndx_group_number(ndx, "LIGAND")
    # ANALYSIS_DEFAULT_BEGIN=10000 (ps) -- rmsf.sh'nin varsayilani.
    same = _standalone_rmsf(ndx, rep, tmp_path / "xval_b10000.xvg",
                            ligand_no, 10000)
    assert len(same) == len(mdkit_profile)
    for got, expected in zip(mdkit_profile, same):
        assert abs(got - expected) < 1e-6, (mdkit_profile, same)

    # -b gercekten uygulandi mi: bastan baslayan ayni hesap farkli olmali.
    from_zero = _standalone_rmsf(ndx, rep, tmp_path / "xval_b0.xvg",
                                 ligand_no, 0)
    assert from_zero != pytest.approx(mdkit_profile, abs=1e-6)


@needs_gmx
@pytest.mark.skipif(not KNOWN_CSV.exists(), reason="rmsf_per_position.csv yok")
def test_bu_makinenin_gromacsi_tarihsel_rmsf_sonucunu_yeniden_uretir():
    """Bu test mdkit'i KOSMAZ ve mdkit'i dogrulamaz.

    Kanitladigi tek sey sudur: bu makinedeki GROMACS, eski `rmsf_analysis.py`
    kosusunun `rmsf_per_position.csv`'ye yazdigi sayilari (Protein grubuna
    fit, -b 10000) bugun de yeniden uretiyor. Yani tarihsel sonuclarla yeni
    sonuclar ayni GROMACS davranisina dayaniyor. mdkit'in gmx cagri
    katmanini dogrulayan test yukaridaki
    test_mdkit_rmsf_ciktisi_bagimsiz_gmx_cagrisiyla_eslesir'dir."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "rmsf_protein.xvg"
        subprocess.run(
            [str(GMX_BIN), "rmsf",
             "-s", str(REAL_REP / "md_0_10.tpr"),
             "-f", str(REAL_REP / "traj_compact_center_dry.xtc"),
             "-o", str(out), "-res", "-b", "10000"],
            input="1\n", text=True, capture_output=True, check=True,
        )
        computed = [v * 10.0 for v in _parse_xvg_values(out)][-10:]

    known = [
        float(row["rmsf_angstrom"])
        for row in csv.DictReader(KNOWN_CSV.read_text().splitlines())
        if row["complex"] == "last10" and row["replica"] == "rep1"
    ]
    assert len(known) == 10
    for got, expected in zip(computed, known):
        assert abs(got - expected) < 0.002, (computed, known)


@needs_gmx
def test_mdkit_rmsf_farkli_bir_buyukluk_olcer(mdkit, real_config, tmp_path):
    """mdkit'in rmsf_pep_self'i LIGAND grubuna fit edilir; rmsf_analysis.py
    ise Protein grubuna fit eder. Bunlar TANIM GEREGI farkli sayilardir ve
    birbirinin yerine kullanilamaz -- bu testi gecmesi o ayrimi kilitler."""
    target = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"
    subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-y", "-a", "rmsf", "-b", "0", str(target)],
        capture_output=True, text=True, check=False,
    )
    out = target / "rep1" / "analysis"
    pep_self = _parse_xvg_values(out / "rmsf_pep_self.xvg")

    out_protein = tmp_path / "rmsf_protein_short.xvg"
    subprocess.run(
        [str(GMX_BIN), "rmsf",
         "-s", str(target / "rep1" / "md_0_10.tpr"),
         "-f", str(target / "rep1" / "traj_compact_center_dry.xtc"),
         "-o", str(out_protein), "-res", "-b", "0"],
        input="1\n", text=True, capture_output=True, check=True,
    )
    protein_fit_peptide = _parse_xvg_values(out_protein)[-10:]

    assert len(pep_self) == 10
    assert pep_self != pytest.approx(protein_fit_peptide, abs=1e-6)
