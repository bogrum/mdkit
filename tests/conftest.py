import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

MDKIT = Path(__file__).resolve().parents[1]

# Testler bu makinenin gercek verisine bakar; baska bir makinede kosmak icin
# asagidaki ortam degiskenleri ayarlanir. (mdkit'in KENDISI hicbir mutlak yol
# icermez -- bkz. test_portability.py; bu istisna yalnizca test verisi icindir.)
REAL_ROOT = Path(os.environ.get("MDKIT_TEST_DATA_ROOT",
                                "/mnt/data/scratch-simulations-TUSEB"))
REAL_REP_NAME = os.environ.get("MDKIT_TEST_REP",
                               "last10_IMGQQPAPQV_A0201_pandora/rep1")
REAL_REP = REAL_ROOT / REAL_REP_NAME
GMX_BIN = Path(os.environ.get("MDKIT_TEST_GMX",
                              "/usr/local/gromacs-2025.4-cuda/bin/gmx"))

needs_gmx = pytest.mark.skipif(
    not (GMX_BIN.exists() and REAL_REP.exists()),
    reason="gmx veya gercek veri yok",
)


def run_bash(snippet, **kw):
    """Bash parcacigini calistir, CompletedProcess dondur."""
    return subprocess.run(["bash", "-c", snippet], capture_output=True, text=True, **kw)


@pytest.fixture
def mdkit():
    return MDKIT


@pytest.fixture
def lib(mdkit):
    return mdkit / "analysis" / "lib.sh"


@pytest.fixture
def fake_dataset(tmp_path):
    """gmx gerektirmeyen sahte veri agaci: 2 kompleks x 3 replika."""
    root = tmp_path / "data"
    for cx in ["last1_AAA_A0201_pandora", "top1_BBB_A0201_pandora"]:
        for rep in ["rep1", "rep2", "rep3"]:
            d = root / cx / rep
            d.mkdir(parents=True)
            (d / "traj_compact_center_dry.xtc").write_text("x")
            (d / "md_0_10.tpr").write_text("x")
            (d / "check_ref.pdb").write_text("x")
    (root / "not_a_complex").mkdir()
    return root


@pytest.fixture
def fake_config(tmp_path, fake_dataset):
    """fake_dataset'i isaret eden gecerli bir config.sh."""
    cfg = tmp_path / "config.sh"
    cfg.write_text(
        textwrap.dedent(
            f"""\
            DATA_ROOT="{fake_dataset}"
            GMX=""
            COMPLEX_GLOB="*_pandora"
            REPS=(rep1 rep2 rep3)
            TRAJ_NAME="traj_compact_center_dry.xtc"
            TPR_NAME="md_0_10.tpr"
            REF_NAME="check_ref.pdb"
            RESULTS_DIR="{tmp_path}/results"
            PYTHON="{sys.executable}"
            CHAIN_RECEPTOR="A"
            CHAIN_AUX="B"
            CHAIN_LIGAND="C"
            """
        )
    )
    return cfg


@pytest.fixture(scope="session")
def small_rep(tmp_path_factory):
    """Gercek veriden 21 frame'lik kucuk replika kopyasi (oturumda bir kez)."""
    if not (GMX_BIN.exists() and REAL_REP.exists()):
        pytest.skip("gmx veya gercek veri yok")
    d = tmp_path_factory.mktemp("small_rep")
    subprocess.run(
        [
            str(GMX_BIN), "trjconv",
            "-s", str(REAL_REP / "md_0_10.tpr"),
            "-f", str(REAL_REP / "traj_compact_center_dry.xtc"),
            "-o", str(d / "traj_compact_center_dry.xtc"),
            "-b", "0", "-e", "200",
        ],
        input="Protein\n", text=True, capture_output=True, check=True,
    )
    shutil.copy(REAL_REP / "check_ref.pdb", d / "check_ref.pdb")
    os.symlink(REAL_REP / "md_0_10.tpr", d / "md_0_10.tpr")
    return d


@pytest.fixture(scope="session")
def small_pair(tmp_path_factory):
    """Gercek veriden IKI replikali kucuk kopya (oturumda bir kez).

    real_config tek replikalidir (REPS=(rep1)), yani cross_rmsd'nin -f2 dali
    orada HIC kosmaz: yalnizca self-matris uretilir. Capraz dalin kendi
    fixture'i olmadan, arac aslinda yapmak icin yazildigi seyi test etmemis
    olurdu."""
    if not (GMX_BIN.exists() and REAL_REP.exists()):
        pytest.skip("gmx veya gercek veri yok")
    src_cx = REAL_REP.parent
    reps = [r for r in ["rep1", "rep2"] if (src_cx / r / "check_ref.pdb").exists()]
    if len(reps) < 2:
        pytest.skip("ikinci gercek replika yok")

    root = tmp_path_factory.mktemp("small_pair") / "data"
    for rep in reps:
        d = root / "test1_PEPTIDE_A0201_pandora" / rep
        d.mkdir(parents=True)
        subprocess.run(
            [
                str(GMX_BIN), "trjconv",
                "-s", str(src_cx / rep / "md_0_10.tpr"),
                "-f", str(src_cx / rep / "traj_compact_center_dry.xtc"),
                "-o", str(d / "traj_compact_center_dry.xtc"),
                "-b", "0", "-e", "200",
            ],
            input="Protein\n", text=True, capture_output=True, check=True,
        )
        shutil.copy(src_cx / rep / "check_ref.pdb", d / "check_ref.pdb")
        os.symlink(src_cx / rep / "md_0_10.tpr", d / "md_0_10.tpr")
    return root, reps


@pytest.fixture
def pair_config(tmp_path, small_pair):
    """small_pair'i iki replikali tek kompleks gibi gosteren config."""
    root, reps = small_pair
    cfg = tmp_path / "pair_config.sh"
    cfg.write_text(
        textwrap.dedent(
            f"""\
            DATA_ROOT="{root}"
            GMX="{GMX_BIN}"
            COMPLEX_GLOB="*_pandora"
            REPS=({" ".join(reps)})
            TRAJ_NAME="traj_compact_center_dry.xtc"
            TPR_NAME="md_0_10.tpr"
            REF_NAME="check_ref.pdb"
            RESULTS_DIR="{tmp_path}/pair_results"
            PYTHON="{sys.executable}"
            CHAIN_RECEPTOR="A"
            CHAIN_AUX="B"
            CHAIN_LIGAND="C"
            """
        )
    )
    return cfg, root, reps


@pytest.fixture
def real_config(tmp_path, small_rep):
    """small_rep'i tek replikali tek kompleks gibi gosteren config + veri agaci.

    DIKKAT: icerik GERCEKTEN last10_IMGQQPAPQV_A0201_pandora/rep1 verisidir,
    yalnizca dizin adi 'test1_PEPTIDE_A0201_pandora' olarak maskelenmistir
    (testler kesif mantigini gercek isimden bagimsiz dogrulasin diye).
    Bu yuzden testlerdeki "10 residue", "275 residue", "151 atom" gibi
    beklenen degerler last10'un gercek sayilaridir.
    """
    root = tmp_path / "data"
    rep = root / "test1_PEPTIDE_A0201_pandora" / "rep1"
    rep.mkdir(parents=True)
    for name in ["traj_compact_center_dry.xtc", "check_ref.pdb", "md_0_10.tpr"]:
        os.symlink(small_rep / name, rep / name)
    cfg = tmp_path / "config.sh"
    cfg.write_text(
        textwrap.dedent(
            f"""\
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
            """
        )
    )
    return cfg


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: tam trajektori gerektiren, dakikalar suren testler"
    )
