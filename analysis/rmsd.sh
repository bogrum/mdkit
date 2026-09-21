#!/usr/bin/env bash
# RMSD eklentisi. run_analysis.sh tarafindan source edilir.
ANALYSIS_NAME="rmsd"
ANALYSIS_DESC="Peptid ve MHC RMSD zaman serileri (oluk-uzeri, ic, global)"
ANALYSIS_KIND="timeseries"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=0
ANALYSIS_OUTPUTS=(rmsd_pep_on_mhc.xvg rmsd_pep_internal.xvg rmsd_complex_bb.xvg)

analysis_run() {
    local rep_dir="$1" out_dir="$2"

    # -tu KULLANILMAZ: -tu ns, -b/-e degerlerini de ns'e cevirir (spec 2.5).
    # Araç icinde zaman daima ps, mesafe daima nm.

    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_pep_on_mhc.xvg" -b "$B_PS" <<< $'RECEPTOR_BB\nLIGAND' \
        || { echo "gmx rms basarisiz: rmsd_pep_on_mhc" >&2; return 1; }

    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_pep_internal.xvg" -b "$B_PS" <<< $'LIGAND_BB\nLIGAND' \
        || { echo "gmx rms basarisiz: rmsd_pep_internal" >&2; return 1; }

    # GLOBAL (tum kompleks) backbone RMSD: Backbone grubu agir zincir + b2m +
    # PEPTIDI birlikte kapsar (gercek indekste 1155 atom / 385 residue, oysa
    # RECEPTOR_BB 825 / 275'tir). Equilibration tespiti icin dogru buyukluk
    # budur; dosya adi bu yuzden "mhc" degil "complex" der.
    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_complex_bb.xvg" -b "$B_PS" <<< $'Backbone\nBackbone' \
        || { echo "gmx rms basarisiz: rmsd_complex_bb" >&2; return 1; }

    return 0
}
