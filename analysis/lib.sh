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

mdkit_find_complexes() {
    # $1 = hedef dizin, $2 = 1 ise hedefin altindaki kompleksleri tara
    local target="${1%/}" all="${2:-0}"
    if [[ "$all" == "1" ]]; then
        find "$target" -mindepth 1 -maxdepth 1 -type d -name "$COMPLEX_GLOB" | sort
    else
        printf '%s\n' "$target"
    fi
}

mdkit_rep_status() {
    # $1 = replika dizini -> OK | NO_TRAJ | NO_TPR | NO_REF
    local rd="$1"
    if [[ ! -s "$rd/$TRAJ_NAME" ]]; then echo "NO_TRAJ"; return 0; fi
    if [[ ! -s "$rd/$TPR_NAME"  ]]; then echo "NO_TPR";  return 0; fi
    if [[ ! -s "$rd/$REF_NAME"  ]]; then echo "NO_REF";  return 0; fi
    echo "OK"
}

mdkit_complex_name() {
    # /yol/<kompleks>_<ek>_<ek> -> <kompleks>
    local base
    base="$(basename "${1%/}")"
    printf '%s\n' "${base%%_*}"
}

mdkit_resolve_gmx() {
    if [[ -n "${GMX:-}" && -x "${GMX}" ]]; then
        return 0
    fi
    local found
    found="$(command -v gmx 2>/dev/null || true)"
    if [[ -n "$found" ]]; then
        GMX="$found"
        return 0
    fi
    echo "mdkit: gmx bulunamadi. config.sh icindeki GMX degerini ayarlayin veya PATH'e ekleyin." >&2
    return 1
}

mdkit_make_ref() {
    # $1 = replika dizini. $TPR_NAME zincir ID'lerini tasiyan topolojidir (spec 2.2).
    local rd="$1" err
    rm -f "$rd/$REF_NAME"
    err="$("$GMX" trjconv -s "$rd/$TPR_NAME" -f "$rd/$TRAJ_NAME" \
        -o "$rd/$REF_NAME" -dump 0 2>&1 >/dev/null <<< "Protein")"
    if [[ ! -s "$rd/$REF_NAME" ]]; then
        echo "mdkit: $REF_NAME uretilemedi: $rd -- $(printf '%s' "$err" | tail -3 | tr '\n' ' ')" >&2
        return 1
    fi
    return 0
}

mdkit_build_index() {
    # $1 = replika dizini, $2 = cikti dizini. Grup NUMARASI aritmetigi yok:
    # make_ndx'e isimle kesisim verilir, sonra basliklar kanonik adlara cevrilir.
    #
    # SURUM BAGIMLILIGI: asagidaki sed, "chain A" komutunun "chA", isimle
    # kesisimin "chA_&_Backbone" adini uretmesine dayanir. GROMACS 2025.4'te
    # dogrulandi. Surum degisirse isimler tutmaz ve fonksiyon SESSIZCE degil,
    # asagidaki kanonik-ad denetiminde GURULTULU sekilde basarisiz olur.
    # (Alternatif olarak "son eklenen 5 grup"u konuma gore yeniden adlandirmak
    # dusunuldu; reddedildi, cunku bos bir zincir secimi grup sayisini degistirip
    # konumsal eslemeyi sessizce KAYDIRIR -- isimle eslesme yanlis olamaz, sadece
    # bulunamaz.)
    local rep_dir="$1" out_dir="$2"
    local ref="$rep_dir/$REF_NAME" ndx="$out_dir/index.ndx"
    local err

    rm -f "$ndx"
    err="$("$GMX" make_ndx -f "$ref" -o "$ndx" 2>&1 >/dev/null <<EOF
chain $CHAIN_RECEPTOR
chain $CHAIN_AUX
chain $CHAIN_LIGAND
"ch$CHAIN_RECEPTOR" & "Backbone"
"ch$CHAIN_LIGAND" & "Backbone"
q
EOF
)"

    if [[ ! -s "$ndx" ]]; then
        echo "mdkit: make_ndx basarisiz: $rep_dir -- $(printf '%s' "$err" | tail -3 | tr '\n' ' ')" >&2
        return 1
    fi

    sed -i \
        -e "s/^\[ ch$CHAIN_RECEPTOR \]/[ RECEPTOR ]/" \
        -e "s/^\[ ch$CHAIN_AUX \]/[ AUX ]/" \
        -e "s/^\[ ch$CHAIN_LIGAND \]/[ LIGAND ]/" \
        -e "s/^\[ ch${CHAIN_RECEPTOR}_&_Backbone \]/[ RECEPTOR_BB ]/" \
        -e "s/^\[ ch${CHAIN_LIGAND}_&_Backbone \]/[ LIGAND_BB ]/" \
        "$ndx"

    local g
    for g in "${MDKIT_GROUPS[@]}"; do
        if ! grep -q "^\[ $g \]" "$ndx"; then
            echo "mdkit: index grubu eksik: $g ($rep_dir)" >&2
            rm -f "$ndx"
            return 1
        fi
    done

    # Boyut denetimi: bos ya da yanlis eslesmis grubu yakalar. Backbone
    # kesisimi residue basina tam 3 atom (N, CA, C) icermelidir.
    local n_lig n_ligbb
    n_lig="$(mdkit_group_size "$ndx" LIGAND)"
    n_ligbb="$(mdkit_group_size "$ndx" LIGAND_BB)"
    if [[ "$n_lig" -eq 0 || "$n_ligbb" -eq 0 || $((n_ligbb % 3)) -ne 0 ]]; then
        echo "mdkit: index grup boyutlari tutarsiz (LIGAND=$n_lig LIGAND_BB=$n_ligbb): $rep_dir" >&2
        rm -f "$ndx"
        return 1
    fi
    return 0
}

mdkit_group_size() {
    # $1 = index dosyasi, $2 = grup adi -> atom sayisi
    awk -v want="$2" '
        /^\[ .* \]$/ { cur = $2; next }
        cur == want   { n += NF }
        END           { print n + 0 }
    ' "$1"
}

mdkit_log_init() {
    mkdir -p "$RESULTS_DIR" || return 1
    MDKIT_LOG="$RESULTS_DIR/run_log.csv"
    if [[ ! -s "$MDKIT_LOG" ]]; then
        echo "timestamp,complex,replica,analysis,status,seconds,error" > "$MDKIT_LOG"
    fi
    return 0
}

mdkit_log_row() {
    # $1=complex $2=replica $3=analysis $4=status $5=seconds $6=error(opsiyonel)
    local err="${6:-}"
    err="${err//\"/\"\"}"      # standart CSV kacisi: tirnagi ikile
    err="${err//$'\n'/ }"      # yeni satiri bosluga cevir
    printf '%s,%s,%s,%s,%s,%s,"%s"\n' \
        "$(date -Iseconds)" "$1" "$2" "$3" "$4" "$5" "$err" >> "$MDKIT_LOG"
}

mdkit_outputs_present() {
    # $1 = cikti dizini, kalani beklenen dosya adlari
    local out_dir="$1"; shift
    local f
    for f in "$@"; do
        [[ -s "$out_dir/$f" ]] || return 1
    done
    return 0
}

mdkit_run_isolated() {
    # $1 = analiz scripti, $2 = rep_dir, $3 = out_dir
    # analysis_run alt kabukta kosar; hata ana donguyu dusurmez.
    local out
    MDKIT_LAST_ERROR=""
    if out="$( source "$1" && analysis_run "$2" "$3" 2>&1 )"; then
        return 0
    fi
    MDKIT_LAST_ERROR="$(printf '%s' "$out" | tail -3 | tr '\n' ' ')"
    return 1
}

mdkit_analysis_scripts() {
    local f
    for f in "$MDKIT_DIR/analysis"/*.sh; do
        [[ -f "$f" ]] || continue
        [[ "$(basename "$f")" == "lib.sh" ]] && continue
        printf '%s\n' "$f"
    done
}

mdkit_analysis_meta() {
    # $1 = analiz scripti -> name<TAB>kind<TAB>desc<TAB>output1,output2
    # Alt kabukta source edilir; degiskenler ana kabuga sizmaz.
    (
        source "$1" || exit 1
        local IFS=,
        printf '%s\t%s\t%s\t%s\n' \
            "$ANALYSIS_NAME" "$ANALYSIS_KIND" "$ANALYSIS_DESC" "${ANALYSIS_OUTPUTS[*]}"
    )
}
