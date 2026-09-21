# mdkit — post-MD analiz araci

GROMACS trajektorileri uzerinde, kompleks x replika duzeninde tekrarlanabilir
analiz kosar ve sonuclari tek bir birlesik tablodan cizdirir.

## Kurulum

Ayri bir kurulum yok. Gereksinimler:

- bash 5+
- GROMACS (2025.4 ile gelistirildi ve test edildi)
- Python 3.11+ ve numpy, pandas, matplotlib (cizim icin); pytest (testler icin)

## Yapilandirma

Projeye ozel **tum** varsayimlar `config.sh` icindedir. Baska bir veri setiyle
calistirmak icin bu dosyanin bir kopyasini duzenleyip `--config` ile verin;
kodda degisiklik gerekmez.

| Degisken | Anlami |
|---|---|
| `DATA_ROOT` | Kompleks dizinlerinin bulundugu kok |
| `RESULTS_DIR` | Birlesik CSV'ler, log ve grafikler |
| `GMX` | gmx ikili dosyasi (bos birakilirsa PATH'ten bulunur) |
| `PYTHON` | numpy/pandas/matplotlib iceren yorumlayici (**PATH fallback'i yoktur**) |
| `COMPLEX_GLOB` | Kompleks dizinlerini eslestiren desen |
| `REPS` | Replika alt dizinleri |
| `TRAJ_NAME` / `TPR_NAME` / `REF_NAME` | Replika icindeki dosya adlari |
| `CHAIN_RECEPTOR` / `CHAIN_AUX` / `CHAIN_LIGAND` | Zincir harfleri |
| `POSTMD_SCRIPT` | Opsiyonel; yalnizca `--full-postmd` icin |

Zincir harfleri index kurulumunda **kanonik adlara** cevrilir:
`RECEPTOR`, `AUX`, `LIGAND`, `RECEPTOR_BB`, `LIGAND_BB`. Analiz scriptleri
yalnizca bu adlari bilir, zincir harflerini gormez.

## Kullanim

```bash
# Tek kompleks, tum replikalar, tum analizler
./run_analysis.sh /veri/kok/last10_IMGQQPAPQV_A0201_pandora

# Tum kompleksler, yalnizca RMSD
./run_analysis.sh --all -a rmsd /veri/kok

# Once ne yapilacagini gor
./run_analysis.sh --all --dry-run /veri/kok

# Mevcut analizleri listele (TSV: ad, kind, aciklama, ciktilar)
./run_analysis.sh --list

# Sonuclari topla ve cizdir
python collect_results.py
python plot_results.py --results-dir /veri/kok/results
```

Butun secenekler icin `./run_analysis.sh --help`.

## Birimler

`.xvg` ciktilari ve birlesik CSV'ler **ham** birimdedir: zaman ps, mesafe nm.
`unit` kolonu bunu belirtir. ns ve Angstrom'a donusum yalnizca cizim
katmanindadir.

> `gmx`'in `-tu ns` secenegi `-b`/`-e` degerlerini de ns'e cevirir. Bu yuzden
> arac icinde `-tu` **hic kullanilmaz**; aksi halde `-b 10000` "10 ns sonrasi"
> yerine "10 000 ns sonrasi" anlamina gelir ve gmx uyarmadan bos cikti uretir.

## Yeni analiz ekleme

`analysis/` icine bir `.sh` dosyasi koyun. `run_analysis.sh` degismez.

```bash
ANALYSIS_NAME="gyrate"
ANALYSIS_DESC="Donme yaricapi zaman serisi"
ANALYSIS_KIND="timeseries"          # timeseries | profile | matrix
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=0            # ps
ANALYSIS_OUTPUTS=(gyrate.xvg)       # idempotency + toplama manifestosu

analysis_run() {
    local rep_dir="$1" out_dir="$2"
    # Kullanilabilir degiskenler: $GMX $REF_PDB $XTC $NDX $B_PS $GROOVE_FIT
    "$GMX" gyrate -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/gyrate.xvg" -b "$B_PS" <<< 'LIGAND'
}
```

`ANALYSIS_KIND` cizim katmaninin hangi fonksiyonu kullanacagini belirler.
`timeseries` ve `profile` icin yeni cizim kodu gerekmez.

**`matrix` henuz desteklenmiyor:** `ANALYSIS_KIND` sozlesme olarak `matrix`
degerini kabul eder, ama ne `collect_results.py` ne de `plot_results.py` bunu
tuketir. `matrix` ilan eden bir analiz calisir ve `.xvg` ciktisi uretir, fakat
toplama asamasinda hicbir CSV'ye girmez; `collect_results.py` bu durumda
stderr'e bir uyari yazar (cikti dosyasi ve kind adiyla, cikti basina en fazla
bir kez). Bu bilinen bir sinirdir, sessiz veri kaybi degil.

## Durum sozlugu

`results/run_log.csv` icindeki `status`:

| Deger | Anlami |
|---|---|
| `OK` | Analiz kostu, ciktilar uretildi |
| `SKIP_DONE` | Ciktilar zaten vardi, kosulmadi (`--force` ile ezilir) |
| `SKIP_MISSING` | On kosul dosyasi yoktu |
| `HATA` | gmx veya dogrulama hatasi; `error` kolonunda mesaj |

## Testler

```bash
cd tests && python -m pytest -v            # hizli testler
cd tests && python -m pytest -v -m slow    # tam trajektori gerektirenler
```

Gercek veriye bakan testler `MDKIT_TEST_DATA_ROOT`, `MDKIT_TEST_REP` ve
`MDKIT_TEST_GMX` ortam degiskenleriyle baska bir makineye yonlendirilebilir.
Veri veya gmx yoksa ilgili testler atlanir.
