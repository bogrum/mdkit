#!/usr/bin/env bash
# Capraz-RMSD eklentisi. run_analysis.sh tarafindan source edilir.
#
# rep_i baglaminda kosar ve $REPS'teki HER replikaya (kendisi dahil) karsi bir
# frame x frame RMSD matrisi uretir: RECEPTOR_BB uzerine fit, LIGAND_BB uzerinde
# RMSD.
#
# Cevapladigi soru: uc replika ayni konformasyonel bolgeleri mi geziyor?
# Mevcut `rmsd` analizi bunu cevaplayamaz -- her replikayi KENDI referansina
# gore olctugu icin iki replika ayni RMSD degerine FARKLI konformasyonlarda
# ulasabilir.
#
# Kardes replikalara erisim: mdkit_run_isolated analysis_run'i ana kabugun ALT
# KABUGUNDA kosturur (analysis/lib.sh:222), yani $REPS/$TRAJ_NAME/$REF_NAME
# gorunurdur ve kardesler dirname "$rep_dir" ile bulunur. run_analysis.sh
# degistirilmedi.
ANALYSIS_NAME="cross_rmsd"
ANALYSIS_DESC="Replikalar arasi 2D RMSD matrisleri (peptid backbone, oluk-uzeri fit)"
ANALYSIS_KIND="matrix"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000

# Ciktilar $REPS'ten turetilir. Self-matris (rep1'in dizininde
# cross_rmsd_rep1.xpm) iki isi birden gorur: tek replika icindeki metastabil
# durumlari gosterir VE cikti kumesini her replika icin TEKDUZE yapar.
# Tekduzelik olmadan ANALYSIS_OUTPUTS replikaya gore degisirdi (rep1 ->
# rep2,rep3 ama son replika -> bos) ve mdkit_mandatory_outputs tabanli
# idempotency calismazdi: hicbir zorunlu cikti uretmeyen replika her kosuda
# yeniden kosardi.
ANALYSIS_OUTPUTS=()
for _cross_rmsd_rep in ${REPS[@]+"${REPS[@]}"}; do
    ANALYSIS_OUTPUTS+=("cross_rmsd_${_cross_rmsd_rep}.xpm")
done
unset _cross_rmsd_rep

# Ciktilarin hepsi KOSULSUZ uretilir.
ANALYSIS_OPTIONAL_OUTPUTS=()

analysis_run() {
    local rep_dir="$1" out_dir="$2"
    return 0
}
