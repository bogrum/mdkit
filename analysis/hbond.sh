#!/usr/bin/env bash
# Peptid-MHC hidrojen bagi eklentisi. run_analysis.sh tarafindan source edilir.
#
# Cevapladigi soru: peptid oluga kac hidrojen bagiyla tutunuyor ve bu sayi
# zaman icinde nasil dalgalaniyor? TCR olmadan immunojeniteye bakilirken
# baglanma sikiligi dogrudan ilgili bir buyukluk.
ANALYSIS_NAME="hbond"
ANALYSIS_DESC="Peptid-MHC hidrojen bagi sayisi (zaman serisi)"
ANALYSIS_KIND="timeseries"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
ANALYSIS_OUTPUTS=(hbond_pep_mhc.xvg)
ANALYSIS_OPTIONAL_OUTPUTS=()

analysis_run() {
    local rep_dir="$1" out_dir="$2" tmp_dir rc
    tmp_dir="$(mktemp -d)" || { echo "gecici dizin acilamadi" >&2; return 1; }
    _hbond_run "$rep_dir" "$out_dir" "$tmp_dir"
    rc=$?
    rm -rf "$tmp_dir"
    return "$rc"
}

_hbond_run() {
    local rep_dir="$1" out_dir="$2" tmp_dir="$3"

    # -s TOPOLOJI olmali, referans PDB DEGIL.
    #
    # OLCULDU (GROMACS 2025.4): ayni komut
    #   -s check_ref.pdb  -> ortalama 0.00 H-bagi
    #   -s md_0_10.tpr    -> ortalama 10.96 H-bagi
    # PDB'de bag bilgisi yoktur, bu yuzden donor/akseptor cikarilamaz ve
    # arac HATA VERMEDEN bos sonuc dondurur. Yani yanlis dosyayi vermek
    # sessiz veri kaybidir; testi bu yuzden var.
    #
    # Topoloji eklentiye run_analysis.sh tarafindan degisken olarak
    # verilmez, ama $TPR_NAME ve $rep_dir gorunurdur (mdkit_run_isolated
    # analysis_run'i ana kabugun alt kabugunda kosturur) -- cross_rmsd'nin
    # kardes replikalara eristigi mekanizmanin aynisi. run_analysis.sh
    # degismedi.
    local tpr="$rep_dir/$TPR_NAME"
    [[ -s "$tpr" ]] || {
        echo "topoloji yok: $tpr" >&2
        return 1
    }

    # NE HESAPLANIYOR: her frame'de peptid ile MHC arasindaki hidrojen
    # bagi SAYISI. Geometrik olcut: donor-akseptor mesafesi ve
    # donor-H-akseptor acisi esikleri (gmx varsayilanlari). -num zaman
    # serisini yazar; -r referans, -t hedef secimdir.
    #
    # Yuksek ve kararli sayi = peptid oluga siki tutunuyor. Dalgalanma =
    # baglar surekli kopup yeniden kuruluyor.
    #
    # Hedef RECEPTOR (agir zincir), AUX (b2m) DEGIL.
    # OLCULDU: bes farkli komplekste de b2m'nin peptide H-bagi katkisi TAM
    # SIFIR; beklenen, cunku b2m olugun karsi tarafindadir. RECEPTOR ile
    # RECEPTOR+AUX ayni sonucu veriyor, o yuzden dar olani seciliyor.
    #
    # Seyreltme YOK: tam trajektoride replika basina ~2 saniye suruyor,
    # seyreltmenin kazanci yok.
    # -o ACIKCA VERILMELIDIR. gmx hbond, -o'suz cagrildiginda H-bagi index
    # dosyasini (hbond.ndx) CALISMA DIZININE yazar -- yani kullanicinin
    # kabuk dizinine. Dahasi GROMACS her kosuda eskisini #hbond.ndx.N#
    # olarak yedekler ve 99. yedekte DURUR:
    #
    #   "Will not make more than 99 backups"
    #
    # 105 replikalik bir kosuda bu, son 21 replikanin HATA vermesi demekti
    # (olculdu). Ayni tuzak spec'te gmx rms -o icin belgelenmisti; buraya
    # da uygulanmasi gerekiyordu.
    #
    # Dosya out_dir'e de yazilamaz: manifestoda ilan edilmeyen bir .ndx her
    # toplamada "manifestoda olmayan dosya" uyarisi uretirdi. Bu yuzden
    # gecici dizine yazilir ve analysis_run cikista siler.
    "$GMX" hbond -s "$tpr" -f "$XTC" -n "$NDX" \
        -r 'group "LIGAND"' -t 'group "RECEPTOR"' \
        -num "$out_dir/hbond_pep_mhc.xvg" -o "$tmp_dir/hbond.ndx" -b "$B_PS" \
        || {
            echo "gmx hbond basarisiz: hbond_pep_mhc" >&2
            return 1
        }

    return 0
}
