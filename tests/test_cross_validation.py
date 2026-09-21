import csv
import subprocess

import pytest

from conftest import GMX_BIN, REAL_REP, REAL_ROOT, needs_gmx

KNOWN_CSV = REAL_ROOT / "rmsf_per_position.csv"

pytestmark = pytest.mark.slow


def _parse_xvg_values(path):
    return [float(ln.split()[1]) for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


@needs_gmx
@pytest.mark.skipif(not KNOWN_CSV.exists(), reason="rmsf_per_position.csv yok")
def test_gmx_cagri_katmani_bilinen_ciktiyi_yeniden_uretir(tmp_path):
    """rmsf_analysis.py'nin urettigi rmsf_per_position.csv, 'Protein' grubuna
    fit edilmis RMSF'dir (-b 10000). Ayni komutu tekrarlayip son 10 residue'yu
    karsilastirmak, gmx cagri katmanini bilinen-dogru bir ciktiya karsi
    dogrular. TAM trajektori gerektirir, birkac dakika surer."""
    out = tmp_path / "rmsf_protein.xvg"
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
        for row in csv.DictReader(KNOWN_CSV.open())
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
