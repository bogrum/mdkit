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

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)   CONFIG="${2:-}";    shift 2 ;;
        -a|--analysis) ANALYSES="${2:-}";  shift 2 ;;
        -r|--reps)     REPS_OPT="${2:-}";  shift 2 ;;
        -b|--begin)    BEGIN_OPT="${2:-}"; shift 2 ;;
        --all)         ALL=1;          shift ;;
        --groove-fit)  GROOVE_FIT=1;   shift ;;
        --force)       FORCE=1;        shift ;;
        --full-postmd) FULL_POSTMD=1;  shift ;;
        -y|--yes)      ASSUME_YES=1;   shift ;;
        -l|--list)     LIST=1;         shift ;;
        -n|--dry-run)  DRY_RUN=1;      shift ;;
        -h|--help)     usage; exit 0 ;;
        -*)            echo "bilinmeyen secenek: $1" >&2; usage >&2; exit 2 ;;
        *)             TARGET="$1";    shift ;;
    esac
done

mdkit_load_config "$CONFIG" || exit 1

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
    for f in "${selected[@]}"; do
        mdkit_analysis_meta "$f"
    done
    exit 0
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
                IFS=$'\t' read -r a_name a_kind a_desc a_outs < <(mdkit_analysis_meta "$f")
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

echo "mdkit: ${#complexes[@]} kompleks, ${#REPS[@]} replika, ${#selected[@]} analiz"
echo "mdkit: yurutme dongusu Task 6'da eklenecek"
