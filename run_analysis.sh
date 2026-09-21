#!/usr/bin/env bash
# mdkit ust script. DIKKAT: set -e YOK -- tek bir replikanin hatasi 105 elemanli
# batch'i dusurmemeli. Hata yalitimi mdkit_run_isolated ile yapilir.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=analysis/lib.sh
source "$SCRIPT_DIR/analysis/lib.sh"

usage() {
    cat <<'USAGE'
kullanim: run_analysis.sh [secenekler] <hedef>

<hedef>   Bir kompleks dizini (icinde rep1/rep2/rep3),
          veya --all ile birlikte kok dizin.

Secenekler:
  -c, --config FILE     Config dosyasi (varsayilan: scriptin yanindaki config.sh)
  -a, --analysis LIST   Virgulle ayrilmis analiz adlari (varsayilan: tumu)
      --all             Hedefi kok kabul et, altindaki tum kompleksleri gez
  -r, --reps LIST       Replika listesi (varsayilan: config'deki REPS)
  -b, --begin PS        Equilibration cutoff (ps); analiz varsayilanini ezer
      --groove-fit      MHC-cerceveli fitlenmis trajektori uret ve kullan
      --force           Cikti varsa bile yeniden hesapla
      --full-postmd     Eksik referans PDB icin post_md_script.sh'in tamamini kos
  -y, --yes             Onay sorularini otomatik onayla
  -l, --list            Analizleri TSV olarak listele ve cik
  -n, --dry-run         Ne yapilacagini yazdir, calistirma
  -h, --help            Bu yardim
USAGE
}

CONFIG=""; ANALYSES=""; TARGET=""; REPS_OPT=""; BEGIN_OPT=""
ALL=0; GROOVE_FIT=0; FORCE=0; FULL_POSTMD=0; ASSUME_YES=0; DRY_RUN=0; LIST=0

require_value() {
    [[ $# -ge 2 ]] || { echo "secenek deger gerektirir: $1" >&2; exit 2; }
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)   require_value "$@"; CONFIG="$2";    shift 2 ;;
        -a|--analysis) require_value "$@"; ANALYSES="$2";  shift 2 ;;
        -r|--reps)     require_value "$@"; REPS_OPT="$2";  shift 2 ;;
        -b|--begin)    require_value "$@"; BEGIN_OPT="$2"; shift 2 ;;
        --all)         ALL=1;          shift ;;
        --groove-fit)  GROOVE_FIT=1;   shift ;;
        --force)       FORCE=1;        shift ;;
        --full-postmd) FULL_POSTMD=1;  shift ;;
        -y|--yes)      ASSUME_YES=1;   shift ;;
        -l|--list)     LIST=1;         shift ;;
        -n|--dry-run)  DRY_RUN=1;      shift ;;
        -h|--help)     usage; exit 0 ;;
        -*)            echo "bilinmeyen secenek: $1" >&2; usage >&2; exit 2 ;;
        *)
            if [[ -n "$TARGET" ]]; then
                echo "birden fazla hedef verildi: $TARGET, $1" >&2
                exit 2
            fi
            TARGET="$1"; shift ;;
    esac
done

mdkit_load_config "$CONFIG" || exit 1

# --full-postmd dogrulamasi HERHANGI bir is baslamadan once yapilir: bu kontrol
# eksik-referans dongusunun icinde dururken, onlarca replika islendikten sonra
# exit 2 ile dusuyordu.
if [[ "$FULL_POSTMD" == 1 ]]; then
    if [[ -z "${POSTMD_SCRIPT:-}" || ! -f "${POSTMD_SCRIPT:-}" ]]; then
        echo "--full-postmd icin config'de POSTMD_SCRIPT tanimli olmali" >&2
        exit 2
    fi
fi

# --- analiz secimi ---
selected=()
if [[ -n "$ANALYSES" ]]; then
    IFS=',' read -r -a wanted <<< "$ANALYSES"
    for w in "${wanted[@]}"; do
        f="$MDKIT_DIR/analysis/$w.sh"
        if [[ ! -f "$f" ]]; then
            echo "bilinmeyen analiz: $w" >&2
            exit 2
        fi
        selected+=("$f")
    done
else
    mapfile -t selected < <(mdkit_analysis_scripts)
fi

if [[ "$LIST" == 1 ]]; then
    list_status=0
    for f in "${selected[@]}"; do
        mdkit_analysis_meta "$f" || { echo "mdkit: metadata okunamadi: $f" >&2; list_status=1; }
    done
    exit "$list_status"
fi

if [[ -z "$TARGET" ]]; then
    echo "hedef dizin verilmedi" >&2
    usage >&2
    exit 2
fi

if [[ ! -d "$TARGET" ]]; then
    echo "dizin degil: $TARGET" >&2
    exit 2
fi

if [[ -n "$REPS_OPT" ]]; then
    IFS=',' read -r -a REPS <<< "$REPS_OPT"
fi

mapfile -t complexes < <(mdkit_find_complexes "$TARGET" "$ALL")
if [[ ${#complexes[@]} -eq 0 ]]; then
    echo "hedef altinda '$COMPLEX_GLOB' eslesen dizin yok: $TARGET" >&2
    exit 1
fi

# --- dry-run: ne yapilacagini yazdir, hicbir sey olusturma ---
if [[ "$DRY_RUN" == 1 ]]; then
    for cx in "${complexes[@]}"; do
        cname="$(mdkit_complex_name "$cx")"
        for rep in "${REPS[@]}"; do
            [[ -d "$cx/$rep" ]] || continue
            st="$(mdkit_rep_status "$cx/$rep")"
            if [[ "$st" != "OK" ]]; then
                echo "DRY-RUN $cname/$rep: $st -- atlanacak"
                continue
            fi
            for f in "${selected[@]}"; do
                IFS=$'\t' read -r a_name a_kind a_desc a_outs a_opt \
                    < <(mdkit_analysis_meta "$f")
                a_begin="$BEGIN_OPT"
                if [[ -z "$a_begin" ]]; then
                    a_begin="$( source "$f"; printf '%s' "$ANALYSIS_DEFAULT_BEGIN" )"
                fi
                echo "DRY-RUN $cname/$rep/$a_name -> ${a_outs//,/ } (b=$a_begin)"
            done
        done
    done
    exit 0
fi

mdkit_resolve_gmx || exit 1
# mdkit_log_init her mdkit_log_row cagrisindan ONCE calismalidir -- MDKIT_LOG'u
# o ayarlar; aksi halde ilk log_row cagrisi "ambiguous redirect" ile patlar.
mdkit_log_init || exit 1

echo "mdkit: ${#complexes[@]} kompleks, ${#REPS[@]} replika, ${#selected[@]} analiz"

# --- eksik referans PDB'leri: tek seferlik onay, sonra kaldigi yerden devam ---
missing_refs=()
for cx in "${complexes[@]}"; do
    for rep in "${REPS[@]}"; do
        [[ -d "$cx/$rep" ]] || continue
        if [[ "$(mdkit_rep_status "$cx/$rep")" == "NO_REF" ]]; then
            missing_refs+=("$cx/$rep")
        fi
    done
done

if [[ ${#missing_refs[@]} -gt 0 ]]; then
    echo "${#missing_refs[@]} replikada $REF_NAME eksik."
    approve=0
    if [[ "$ASSUME_YES" == 1 ]]; then
        approve=1
    else
        read -r -p "Mevcut kuru trajektoriden uretilsin mi? [e/H] " answer
        [[ "$answer" == [eE] ]] && approve=1
    fi
    if [[ "$approve" == 1 ]]; then
        for rd in "${missing_refs[@]}"; do
            if [[ "$FULL_POSTMD" == 1 ]]; then
                ( cd "$rd" && bash "$POSTMD_SCRIPT" ) >/dev/null 2>&1 \
                    || echo "post_md_script basarisiz: $rd" >&2
            else
                mdkit_make_ref "$rd" || echo "$REF_NAME uretilemedi: $rd" >&2
            fi
        done
    else
        echo "uretilmedi; bu replikalar SKIP_MISSING olarak loglanacak."
    fi
fi

# --- ana dongu ---
for cx in "${complexes[@]}"; do
    cname="$(mdkit_complex_name "$cx")"
    for rep in "${REPS[@]}"; do
        rep_dir="$cx/$rep"
        [[ -d "$rep_dir" ]] || continue

        status="$(mdkit_rep_status "$rep_dir")"
        if [[ "$status" != "OK" ]]; then
            echo "$cname/$rep: $status -- atlaniyor"
            mdkit_log_row "$cname" "$rep" "-" SKIP_MISSING 0 "${BEGIN_OPT:--}" "$status"
            continue
        fi

        out_dir="$rep_dir/analysis"
        mkdir -p "$out_dir"
        [[ "$FORCE" == 1 ]] && rm -f "$out_dir/index.ndx"

        for script in "${selected[@]}"; do
            # Onceki eklentinin sozlesmesi TEMIZLENIR. Aksi halde eksik bir
            # metadata ya set -u ile tum batch'i dusurur ya da bir onceki
            # eklentinin degerini (hatta analysis_run govdesini) miras alir.
            mdkit_clear_plugin
            # shellcheck source=/dev/null
            source "$script" || true
            if ! mdkit_validate_plugin "$script"; then
                echo "$cname/$rep: HATA -- $MDKIT_LAST_ERROR" >&2
                mdkit_log_row "$cname" "$rep" "$(basename "$script" .sh)" HATA 0 \
                    "${BEGIN_OPT:--}" "$MDKIT_LAST_ERROR"
                continue
            fi
            begin_ps="${BEGIN_OPT:-$ANALYSIS_DEFAULT_BEGIN}"

            # Idempotency yalnizca ZORUNLU ciktilara bakar; opsiyonel ciktilar
            # (ör. --groove-fit) manifestoda ilan edilir ama uretilmeyebilir.
            mapfile -t mandatory < <(mdkit_mandatory_outputs)
            if [[ "$FORCE" != 1 ]] && [[ ${#mandatory[@]} -gt 0 ]] && \
               mdkit_outputs_present "$out_dir" "${mandatory[@]}"; then
                echo "$cname/$rep/$ANALYSIS_NAME: mevcut -- atlaniyor"
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" SKIP_DONE 0 "$begin_ps"
                continue
            fi

            if [[ "${ANALYSIS_NEEDS_INDEX:-0}" == 1 && ! -s "$out_dir/index.ndx" ]]; then
                if ! mdkit_build_index "$rep_dir" "$out_dir"; then
                    mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" HATA 0 \
                        "$begin_ps" "index kurulamadi"
                    continue
                fi
            fi

            # Analiz scriptlerinin gordugu sozlesme degiskenleri
            REF_PDB="$rep_dir/$REF_NAME"
            XTC="$rep_dir/$TRAJ_NAME"
            NDX="$out_dir/index.ndx"
            B_PS="$begin_ps"
            FORCE="$FORCE"

            t0=$SECONDS
            if mdkit_run_isolated "$script" "$rep_dir" "$out_dir"; then
                echo "$cname/$rep/$ANALYSIS_NAME: tamam ($((SECONDS - t0))s)"
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" OK \
                    "$((SECONDS - t0))" "$begin_ps"
            else
                echo "$cname/$rep/$ANALYSIS_NAME: HATA -- $MDKIT_LAST_ERROR" >&2
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" HATA \
                    "$((SECONDS - t0))" "$begin_ps" "$MDKIT_LAST_ERROR"
            fi
        done
    done
done

echo "mdkit: bitti. log -> $MDKIT_LOG"
