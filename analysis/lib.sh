#!/usr/bin/env bash
# mdkit ortak altyapi. run_analysis.sh ve analiz scriptleri bunu source eder.
# Bu dosya analiz mantigi icermez.

MDKIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Index kurulumunun uretecegi kanonik grup adlari.
MDKIT_GROUPS=(RECEPTOR AUX LIGAND RECEPTOR_BB LIGAND_BB)

mdkit_load_config() {
    local cfg="${1:-$MDKIT_DIR/config.sh}"
    if [[ ! -f "$cfg" ]]; then
        echo "mdkit: config bulunamadi: $cfg" >&2
        return 1
    fi
    # shellcheck source=/dev/null
    source "$cfg" || { echo "mdkit: config okunamadi: $cfg" >&2; return 1; }

    # PYTHON zorunludur ve GMX'in aksine PATH fallback'i YOKTUR: PATH'teki
    # python3 bu makinede bilimsel yigina sahip olmayan bir venv'dir, yani
    # fallback sessizce yanlis yorumlayiciyi secip collect/plot adiminda patlardi.
    local missing=() v
    for v in DATA_ROOT COMPLEX_GLOB TRAJ_NAME TPR_NAME REF_NAME RESULTS_DIR \
             PYTHON CHAIN_RECEPTOR CHAIN_AUX CHAIN_LIGAND; do
        [[ -n "${!v:-}" ]] || missing+=("$v")
    done
    if [[ -z "${REPS+x}" ]] || [[ ${#REPS[@]} -eq 0 ]]; then
        missing+=("REPS")
    fi
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "mdkit: config eksik degisken: ${missing[*]}  ($cfg)" >&2
        return 1
    fi
    return 0
}
