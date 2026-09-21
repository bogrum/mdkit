#!/usr/bin/env bash
# RMSF eklentisi. run_analysis.sh tarafindan source edilir.
ANALYSIS_NAME="rmsf"
ANALYSIS_DESC="Residue bazli RMSF profilleri (peptid ic, MHC, opsiyonel oluk-cercevesi)"
ANALYSIS_KIND="profile"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
ANALYSIS_OUTPUTS=(rmsf_pep_self.xvg rmsf_mhc.xvg)

analysis_run() {
    echo "rmsf: Task 8'de doldurulacak" >&2
    return 1
}
