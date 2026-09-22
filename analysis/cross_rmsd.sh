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

# Es basina IKI cikti: .xpm ve .dat.
#   .xpm -> eksen zamanlari, birim, legend (ama degerler 80 renk seviyesine
#           yuvarlanmistir)
#   .dat -> gmx rms -bin ham dump'i; tam float32 degerler, ama HICBIR
#           metadata yok (baslik bile yok)
# Ikisi de zorunlu: tek baslarina eksiktirler. matrix_summary.csv'deki
# min/mean/max tez tablosuna girecegi icin yuvarlanmis degerlerle yetinilmiyor.
#
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
    ANALYSIS_OUTPUTS+=("cross_rmsd_${_cross_rmsd_rep}.dat")
done
unset _cross_rmsd_rep

# Ciktilarin hepsi KOSULSUZ uretilir.
ANALYSIS_OPTIONAL_OUTPUTS=()

_cross_rmsd_atom_count() {
    # $1 = PDB. Mutlak dogruluk degil, iki replika ARASINDA esitlik onemli.
    local n
    n="$(grep -c -E '^(ATOM|HETATM)' "$1" 2>/dev/null)" || n=0
    printf '%s' "${n:-0}"
}

_cross_rmsd_check_peers() {
    # $1 = kompleks dizini. TUM esler, HICBIR gmx cagrisi yapilmadan ONCE
    # dogrulanir.
    #
    # Sira onemli: dogrulama hesap dongusunun ICINDE olsaydi self-matris (ilk
    # es) once uretilir, bozuk bir ikinci es ancak dakikalar suren bir gmx
    # kosusundan SONRA fark edilirdi -- ve out_dir'de yarim bir cikti kumesi
    # kalirdi. Yarim kume, mdkit_mandatory_outputs acisindan "eksik" demektir,
    # yani analiz her kosuda bastan kosardi.
    local cx_dir="$1" peer peer_dir own_atoms peer_atoms
    own_atoms="$(_cross_rmsd_atom_count "$REF_PDB")"

    for peer in ${REPS[@]+"${REPS[@]}"}; do
        peer_dir="$cx_dir/$peer"

        if [[ ! -s "$peer_dir/$TRAJ_NAME" ]]; then
            echo "es replikanin trajektorisi yok: $peer_dir/$TRAJ_NAME" >&2
            return 1
        fi

        # Capraz rms TEK bir index dosyasini IKI trajektoriye birden uygular.
        # Esin sistemi farkliysa (atom sayisi/sirasi) matris SESSIZCE cop olur.
        peer_atoms="$(_cross_rmsd_atom_count "$peer_dir/$REF_NAME")"
        if [[ "$peer_atoms" != "$own_atoms" ]]; then
            echo "es replikanin sistemi uyusmuyor ($peer: $peer_atoms atom, kendi: $own_atoms)" >&2
            return 1
        fi
    done
    return 0
}

_cross_rmsd_pairs() {
    # $1 = kompleks dizini, $2 = out_dir, $3 = gecici .xvg yolu.
    # Esler _cross_rmsd_check_peers tarafindan ZATEN dogrulandi.
    local cx_dir="$1" out_dir="$2" tmp_xvg="$3"
    local peer peer_xtc out out_bin dt
    dt="${CROSS_RMSD_DT:-200}"

    for peer in ${REPS[@]+"${REPS[@]}"}; do
        peer_xtc="$cx_dir/$peer/$TRAJ_NAME"
        out="$out_dir/cross_rmsd_${peer}.xpm"
        out_bin="$out_dir/cross_rmsd_${peer}.dat"

        # -skip KULLANILMAZ: .xpm eksen zamanlarini bozar -- ilk zaman dogru,
        #   gerisi 0 yazilir (GROMACS 2025.4'te olculdu). Seyreltme -dt ile.
        # -tu KULLANILMAZ: -b/-e degerlerini de cevirir.
        # -b HER IKI trajektoriye de uygulanir (olculdu), es icin ayrica
        #   kirpma gerekmez.
        # -o acikca verilir: gmx rms -o'suz cagrildiginda rmsd.xvg'yi CALISMA
        #   DIZININE yazar. out_dir'e yazmak da olmaz -- manifestoda ilan
        #   edilmeyen bir .xvg her toplamada uyari uretirdi.
        if [[ "$peer_xtc" == "$XTC" ]]; then
            # Self-matris: -f2 verilmez. Olculen ve kosegeni sifir dogrulanan
            # bicim budur.
            "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
                -m "$out" -bin "$out_bin" -o "$tmp_xvg" -b "$B_PS" -dt "$dt" \
                <<< $'RECEPTOR_BB\nLIGAND_BB' \
                || { echo "gmx rms basarisiz: cross_rmsd_${peer}" >&2; return 1; }
        else
            "$GMX" rms -s "$REF_PDB" -f "$XTC" -f2 "$peer_xtc" -n "$NDX" \
                -m "$out" -bin "$out_bin" -o "$tmp_xvg" -b "$B_PS" -dt "$dt" \
                <<< $'RECEPTOR_BB\nLIGAND_BB' \
                || { echo "gmx rms basarisiz: cross_rmsd_${peer}" >&2; return 1; }
        fi

        [[ -s "$out" && -s "$out_bin" ]] || {
            echo "matris uretilmedi: cross_rmsd_${peer} (.xpm ve .dat birlikte gerekli)" >&2
            return 1
        }
    done
    return 0
}

analysis_run() {
    local rep_dir="$1" out_dir="$2" cx_dir tmp_dir rc
    cx_dir="$(dirname "$rep_dir")"

    _cross_rmsd_check_peers "$cx_dir" || return 1

    tmp_dir="$(mktemp -d)" || { echo "gecici dizin acilamadi" >&2; return 1; }
    _cross_rmsd_pairs "$cx_dir" "$out_dir" "$tmp_dir/rms.xvg"
    rc=$?
    rm -rf "$tmp_dir"
    return "$rc"
}
