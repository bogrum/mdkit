#!/usr/bin/env bash
# RMSF eklentisi. run_analysis.sh tarafindan source edilir.
#
# DIKKAT (spec 2.4): gmx rmsf varsayilan olarak fit yapar ve fit'i SECILEN
# GRUBUN uzerinde yapar. Bu yuzden:
#   LIGAND + -fit            -> peptidin IC esnekligi (kayma/sallanma cikarilmis)
#   RECEPTOR_BB'ye fitli traj + LIGAND + -nofit -> OLUK CERCEVESINDEKI esneklik
ANALYSIS_NAME="rmsf"
ANALYSIS_DESC="Residue bazli RMSF profilleri (peptid ic, MHC, opsiyonel oluk-cercevesi)"
ANALYSIS_KIND="profile"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
# ILAN kosulsuzdur, URETIM kosulludur. Manifesto (`--list`) --groove-fit
# bilmeyen TAZE bir kabukta okunur; ucuncu ciktiyi kosullu ilan etmek onu
# manifestodan dusuruyor, collect_results.py da diskteki dosyayi SESSIZCE
# atliyordu -- saatlerce hesaplanan 105 dosya hicbir CSV'ye girmiyordu.
ANALYSIS_OUTPUTS=(rmsf_pep_self.xvg rmsf_mhc.xvg rmsf_pep_groovefit.xvg)

# --groove-fit verilmediginde ucuncu cikti mesru sekilde YOKTUR; idempotency
# kontrolu onu aramamalidir. Verildiginde zorunlu olur, yoksa daha once
# --groove-fit'siz kosmus bir replika SKIP_DONE alip ucuncu ciktiyi hic
# uretmezdi.
if [[ "${GROOVE_FIT:-0}" == "1" ]]; then
    ANALYSIS_OPTIONAL_OUTPUTS=()
else
    ANALYSIS_OPTIONAL_OUTPUTS=(rmsf_pep_groovefit.xvg)
fi

analysis_run() {
    local rep_dir="$1" out_dir="$2"

    "$GMX" rmsf -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsf_pep_self.xvg" -res -b "$B_PS" <<< 'LIGAND' \
        || { echo "gmx rmsf basarisiz: rmsf_pep_self" >&2; return 1; }

    "$GMX" rmsf -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsf_mhc.xvg" -res -b "$B_PS" <<< 'RECEPTOR' \
        || { echo "gmx rmsf basarisiz: rmsf_mhc" >&2; return 1; }

    if [[ "${GROOVE_FIT:-0}" == "1" ]]; then
        local fitted="$rep_dir/traj_fit_mhc.xtc"
        # Fitli trajektori saklanir: ileride PCA/DCCM ayni dosyayi isteyecek.
        # --force ile bu eklentinin sahip oldugu artefakt da yeniden kurulur
        # (index.ndx'in --force ile yeniden kurulmasiyla simetrik).
        [[ "${FORCE:-0}" == "1" ]] && rm -f "$fitted"
        if [[ ! -s "$fitted" ]]; then
            "$GMX" trjconv -s "$REF_PDB" -f "$XTC" -n "$NDX" \
                -fit rot+trans -o "$fitted" <<< $'RECEPTOR_BB\nSystem' \
                || { echo "gmx trjconv -fit basarisiz" >&2; rm -f "$fitted"; return 1; }
        fi
        "$GMX" rmsf -s "$REF_PDB" -f "$fitted" -n "$NDX" \
            -o "$out_dir/rmsf_pep_groovefit.xvg" -res -nofit -b "$B_PS" <<< 'LIGAND' \
            || { echo "gmx rmsf basarisiz: rmsf_pep_groovefit" >&2; return 1; }
    fi

    return 0
}
