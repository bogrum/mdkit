import re

import pytest

from conftest import needs_gmx, run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && mdkit_resolve_gmx && '


def test_gmx_bulunamazsa_anlasilir_hata(lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        'PATH=/nonexistent mdkit_resolve_gmx'
    )
    assert r.returncode != 0
    assert "gmx bulunamadi" in r.stderr


@needs_gmx
def test_make_ref_zincirleri_koruyarak_uretir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=real_config) + f'mdkit_make_ref "{rep}"')
    assert r.returncode == 0, r.stderr
    text = (rep / "check_ref.pdb").read_text()
    chains = {ln[21] for ln in text.splitlines() if ln.startswith("ATOM")}
    assert chains == {"A", "B", "C"}


@needs_gmx
def test_build_index_kanonik_gruplari_uretir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    r = run_bash(
        PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"'
    )
    assert r.returncode == 0, r.stderr
    names = re.findall(r"^\[ (.*) \]$", (out / "index.ndx").read_text(), re.M)
    for g in ["RECEPTOR", "AUX", "LIGAND", "RECEPTOR_BB", "LIGAND_BB"]:
        assert g in names
    assert not any(n.startswith("ch") for n in names)


@needs_gmx
def test_build_index_grup_boyutlari_dogru(lib, real_config, tmp_path):
    """LIGAND_BB, peptid residue sayisinin 3 kati atom icermeli (N, CA, C)."""
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    r = run_bash(PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"')
    assert r.returncode == 0, r.stderr

    counts, cur = {}, None
    for line in (out / "index.ndx").read_text().splitlines():
        m = re.match(r"^\[ (.*) \]$", line)
        if m:
            cur = m.group(1)
            counts[cur] = 0
        elif cur:
            counts[cur] += len(line.split())

    assert counts["LIGAND"] == 151        # last10 peptidi, 10 residue
    assert counts["LIGAND_BB"] == 30      # 10 residue x 3 backbone atomu
    assert counts["RECEPTOR_BB"] == 825   # 275 residue x 3
    assert counts["LIGAND_BB"] % 3 == 0


@needs_gmx
def test_build_index_bos_zincirde_hata_verir(lib, real_config, tmp_path):
    """Var olmayan bir zincir harfi verilirse index kurulmamali (bos grup
    sessizce kabul edilirse sonraki analizler anlamsiz sayi uretir)."""
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    cfg = tmp_path / "config_badchain.sh"
    cfg.write_text(
        real_config.read_text().replace('CHAIN_LIGAND="C"', 'CHAIN_LIGAND="Z"')
    )
    r = run_bash(PRE.format(lib=lib, cfg=cfg) + f'mdkit_build_index "{rep}" "{out}"')
    assert r.returncode != 0
    assert not (out / "index.ndx").exists()
    assert "mdkit:" in r.stderr, r.stderr
    # gmx make_ndx silently skips creating a group for a chain letter that
    # matches 0 atoms ("Group is empty"), so LIGAND (and LIGAND_BB, which
    # depends on it) never gets created at all -- this is caught by the
    # canonical-name presence loop, not the later size check.
    assert "index grubu eksik: LIGAND" in r.stderr, r.stderr


@needs_gmx
def test_build_index_bozuk_referansta_hata_verir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    (rep / "check_ref.pdb").unlink()
    (rep / "check_ref.pdb").write_text("bu bir pdb degil\n")
    r = run_bash(
        PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"'
    )
    assert r.returncode != 0
    assert "mdkit:" in r.stderr, r.stderr
