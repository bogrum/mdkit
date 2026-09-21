#!/usr/bin/env bash
# RMSD eklentisi. run_analysis.sh tarafindan source edilir.
ANALYSIS_NAME="rmsd"
ANALYSIS_DESC="Peptid ve MHC RMSD zaman serileri (oluk-uzeri, ic, global)"
ANALYSIS_KIND="timeseries"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=0
ANALYSIS_OUTPUTS=(rmsd_pep_on_mhc.xvg rmsd_pep_internal.xvg rmsd_mhc_bb.xvg)

analysis_run() {
    echo "rmsd: Task 7'de doldurulacak" >&2
    return 1
}
