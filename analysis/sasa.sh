#!/usr/bin/env bash
# Cozucuye acik yuzey (SASA) eklentisi. run_analysis.sh tarafindan source edilir.
#
# Cevapladigi soru: peptidin ne kadari olugun disinda, yani TCR'in gorebilecegi
# yerde? RMSF "ne kadar oynuyor" der; SASA "ne kadar gorunuyor" der. TCR
# olmadan immunojeniteye bakilirken ikincisi soruya daha yakindir.
ANALYSIS_NAME="sasa"
ANALYSIS_DESC="Peptidin acik yuzeyi: kompleks icinde ve tek basina"
ANALYSIS_KIND="timeseries"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
ANALYSIS_OUTPUTS=(sasa_pep_in_complex.xvg sasa_pep_alone.xvg)
ANALYSIS_OPTIONAL_OUTPUTS=()

# GOMULU YUZEY BURADA SAKLANMAZ. alone - in_complex farkiyla turetilir ve
# ucuncu bir dosya olarak yazmak iki saklanan buyuklugun farkini veriyi
# cogaltarak tekrarlamak olurdu.

analysis_run() {
    local rep_dir="$1" out_dir="$2" dt
    dt="${SASA_DT:-100}"

    # SEYRELTME SART. Olculdu (GROMACS 2025.4, 9001 frame):
    #   tam (10 ps)  190 s   ortalama 5.5857 nm^2
    #   -dt 100       19.7 s ortalama 5.5841 nm^2   (10 kat hizli, %0.03 fark)
    #   -dt 200       10.3 s ortalama 5.5737 nm^2
    # 105 replika x 2 kosu: tam cozunurlukte ~11 saat, -dt 100 ile ~70 dakika.
    # Ortalamaya etkisi olcum hatasinin altinda.
    #
    # -or / -oa KULLANILMAZ. Olculdu: -or, yuzey secimi ile -output secimini
    # AYNI dosyada arka arkaya yazar (386 satir) ve residue numaralari zincir
    # basina sifirlandigi icin 100 TEKRARLI numara icerir. collect_results'in
    # profile kind'i residue'yu benzersiz varsayar; bu dosya onu sessizce
    # bozardi. Residue bazli SASA icin blok-farkinda bir ayristirici gerekir.

    # 1) PEPTIDIN KOMPLEKS ICINDEKI ACIK YUZEYI
    # Cozucuye acik yuzey, bir prob kuresi (varsayilan 0.14 nm, su
    # yaricapi) molekulun uzerinde yuvarlanarak hesaplanir; degdigi alan
    # "acik" sayilir.
    #
    # -surface TUM protein verilir, boylece MHC peptidi FIZIKSEL OLARAK
    # ORTER; -output ile yalnizca LIGAND'in payi raporlanir. Cikti iki
    # kolonlu: 1) tum kompleksin yuzeyi, 2) bunun icinde peptide dusen pay.
    # Ikincisi TCR'in gorebilecegi yuzeydir.
    "$GMX" sasa -s "$REF_PDB" -f "$XTC" -n "$NDX" -b "$B_PS" -dt "$dt" \
        -surface 'group "Protein"' -output 'group "LIGAND"' \
        -o "$out_dir/sasa_pep_in_complex.xvg" \
        || {
            echo "gmx sasa basarisiz: sasa_pep_in_complex" >&2
            return 1
        }

    # 2) PEPTIDIN TEK BASINA ACIK YUZEYI
    # -surface yalnizca LIGAND: MHC hesaba hic girmez, yani orten bir sey
    # yok. Peptidin "tamamen aciktaki" yuzeyi budur ve gomulu yuzeyin
    # referansidir: gomulu = (2) - (1). Oran olarak okununca "peptidin
    # yuzde kaci olugun icinde saklaniyor" sorusunu cevaplar.
    "$GMX" sasa -s "$REF_PDB" -f "$XTC" -n "$NDX" -b "$B_PS" -dt "$dt" \
        -surface 'group "LIGAND"' \
        -o "$out_dir/sasa_pep_alone.xvg" \
        || {
            echo "gmx sasa basarisiz: sasa_pep_alone" >&2
            return 1
        }

    return 0
}
