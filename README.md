# mdkit — post-MD analiz araci

GROMACS trajektorileri uzerinde, **kompleks x replika** duzeninde tekrarlanabilir
analiz kosar; ciktilari tek bir birlesik tabloda toplar ve oradan cizdirir.

Tek bir trajektoriyi elle analiz etmek icin degil, ayni analizi 105 trajektoride
gozetimsiz kosup sonra hepsini birlikte gormek icin yazildi.

---

## Calisma prensibi

### Boru hatti

```
config.sh ──► run_analysis.sh ──► analysis/<ad>.sh ──► rep*/analysis/*.{xvg,xpm,dat}
                  │                    (eklenti)              │
                  │                                           ▼
                  └──► results/run_log.csv          collect_results.py
                                                              │
                                                              ▼
                                            results/{timeseries,profile}_long.csv
                                            results/matrices/*.npz + matrix_summary.csv
                                                              │
                                                              ▼
                                                     plot_results.py
                                                              │
                                                              ▼
                                                   results/plots/*.png
```

### Uc katman

| Katman | Dosyalar | Sorumlulugu | Bilmedigi sey |
|---|---|---|---|
| **Yapilandirma** | `config.sh` | Projeye ozel her sey: yollar, dosya adlari, zincir harfleri, replika listesi | Analizlerin ne yaptigi |
| **Kosum (bash)** | `run_analysis.sh`, `analysis/lib.sh`, `analysis/<ad>.sh` | Kesif, index kurulumu, gmx cagrilari, log, idempotency, hata yalitimi | Sonuclarin nasil gosterilecegi |
| **Toplama/cizim (Python)** | `collect_results.py`, `plot_results.py` | `.xvg`/`.xpm` ayristirma, birlesik CSV + matris `.npz`, grafik | gmx'in nasil cagrildigi |

Katmanlar arasindaki tek sozlesme iki sey: **`run_analysis.sh --list`** ciktisi
(hangi cikti hangi analize ve hangi veri sekline ait) ve **cikti dosyalarinin
kendisi** (`.xvg` 1B, `.xpm` 2B). Python tarafi gmx komutlarini hic bilmez; bash tarafi pandas'i hic
bilmez.

### Bir kosuda ne olur

`./run_analysis.sh --all /veri/kok` dediginizde sirayla:

1. **Config yuklenir ve dogrulanir.** Zorunlu bir degisken eksikse ADIYLA
   bildirilir ve cikilir. (`PYTHON`'un PATH fallback'i yoktur: yanlis
   yorumlayici sessizce secilip cizim asamasinda patlamasin diye.)
2. **Hedefler kesfedilir.** `COMPLEX_GLOB` eslesen dizinler, her birinde `REPS`
   alt dizinleri.
3. **On kosullar kontrol edilir.** Her replika icin trajektori, tpr ve referans
   PDB. Eksik olan `SKIP_MISSING` olarak loglanir; kosu durmaz.
4. **Eksik referans PDB'ler icin TEK bir onay sorulur.** Kabul edilirse mevcut
   kuru trajektoriden `-dump 0` ile uretilir (saniyeler), sonra analiz kaldigi
   yerden devam eder. `-y` ile sorulmaz.
5. **Index kurulur** (gereken analizler icin, replika basina bir kez):
   zincir harfleri `chain A/B/C` ile secilir, Backbone kesisimleri alinir, ve
   basliklar **kanonik adlara** cevrilir: `RECEPTOR`, `AUX`, `LIGAND`,
   `RECEPTOR_BB`, `LIGAND_BB`. Bes adin da olustugu ve grup boyutlarinin
   tutarli oldugu dogrulanir.
6. **Her (replika x analiz) icin:** ciktilar zaten varsa `SKIP_DONE`; yoksa
   eklenti **alt kabukta** kosturulur. Cikis kodu, sure ve hata mesaji
   `run_log.csv`'ye yazilir.
7. **Hicbir hata kosuyu durdurmaz.** Bir replikanin gmx hatasi `HATA` satiri
   olur; sonraki replika ve sonraki kompleks normal devam eder.

Sonra ayri bir adimda `collect_results.py` butun `.xvg`'leri iki uzun-format
CSV'ye, butun `.xpm`/`.dat` ciftlerini de `results/matrices/*.npz` + `matrix_summary.csv`
haline indirger; `plot_results.py` de yalnizca `results/` altini okuyarak cizer.

### Diskte ne nereye yazilir

```
$DATA_ROOT/
├── <kompleks>/rep{1,2,3}/
│   ├── traj_compact_center_dry.xtc     (girdi, dokunulmaz)
│   ├── md_0_10.tpr                     (girdi, dokunulmaz)
│   ├── check_ref.pdb                   (yoksa uretilir)
│   ├── traj_fit_mhc.xtc                (yalnizca --groove-fit; saklanir)
│   └── analysis/
│       ├── index.ndx                   (bir kez kurulur)
│       ├── *.xvg                       (ham 1B ciktilar, trajektorinin yaninda)
│       ├── *.xpm                       (ham 2B ciktilar: eksenler + birim)
│       └── *.dat                       (ayni matrisin tam float32 degerleri)
└── results/
    ├── run_log.csv
    ├── timeseries_long.csv
    ├── profile_long.csv
    ├── matrix_summary.csv
    ├── matrices/<kompleks>_<replika>_<cikti>.npz
    └── plots/
        ├── per_complex/<kompleks>_<cikti>.png
        ├── mean_sd/<kompleks>_<cikti>.png
        ├── matrix/<kompleks>_<analiz>.png
        ├── matrix/<kompleks>_<analiz>_equilibrium.png
        └── compare_<cikti>.png
```

Ham `.xvg`/`.xpm` bilerek trajektorinin yaninda durur: dizin tasinirsa sonuc da
birlikte gider ve hangi kosudan geldigi belirsizlesmez. Birlesik CSV ise
cizim katmaninin 105 dizini taramasini onler.

### Dayandigi dort karar

**1. Projeye ozel her sey tek dosyada.** `config.sh` disinda hicbir calisma
dosyasinda mutlak yol ya da proje adi yoktur; bu bir testle zorlanir
(`tests/test_portability.py`). Baska bir veri setiyle calistirmak kod
degistirmeyi degil, config kopyalamayi gerektirir.

**2. Gruplar isimle secilir, numarayla degil.** `make_ndx`'in urettigi
numaralar sisteme gore kayabilir; kanonik adlar kaymaz. Ustelik gmx secimi
`.xvg` basligina yazar — `@ subtitle "LIGAND after lsq fit to RECEPTOR_BB"` —
yani **cikti dosyasi kendi tanimini tasir.**

**3. Cizim, analizin ADINA degil verinin SEKLINE bakar.** `ANALYSIS_KIND`
(`timeseries` / `profile`) hangi cizim fonksiyonunun kullanilacagini belirler.
Bu yuzden yeni bir zaman serisi analizi (Rg, SASA) eklendiginde cizim koduna
hic dokunulmaz.

**4. Hata yalitimi ve idempotency.** `set -e` hic kullanilmaz; her eklenti alt
kabukta kosar ve hatasi loglanip gecilir. Ciktilari zaten var olan is
tekrarlanmaz, boylece yarida kesilen bir batch kaldigi yerden devam eder
(`--force` bunu ezer).

---

## Kurulum

Ayri bir kurulum yok. Gereksinimler:

- bash 5+
- GROMACS (2025.4 ile gelistirildi ve test edildi)
- Python 3.11+ ve numpy, pandas, matplotlib (cizim icin); pytest (testler icin).
  Tamami `requirements.txt` icinde: `pip install -r requirements.txt`.

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
| `COMPLEX_GROUPS` | Opsiyonel; cizimde kompleks gruplari, `"onek:etiket"` listesi |

Zincir harfleri index kurulumunda kanonik adlara cevrilir; analiz scriptleri
yalnizca `RECEPTOR`/`AUX`/`LIGAND`/`RECEPTOR_BB`/`LIGAND_BB` adlarini bilir,
zincir harflerini hic gormez.

## Kullanim

```bash
# Tek kompleks, tum replikalar, tum analizler
./run_analysis.sh /veri/kok/<kompleks>

# Tum kompleksler, yalnizca RMSD
./run_analysis.sh --all -a rmsd /veri/kok

# Once ne yapilacagini gor (hicbir sey yazmaz)
./run_analysis.sh --all --dry-run /veri/kok

# Mevcut analizleri listele
# TSV: ad, kind, aciklama, tum ciktilar, bunlarin opsiyonel olanlari
./run_analysis.sh --list

# Replikalar arasi 2D RMSD matrisleri (seyreltme varsayilani 200 ps)
./run_analysis.sh --all -a cross_rmsd /veri/kok
CROSS_RMSD_DT=500 ./run_analysis.sh --all -a cross_rmsd /veri/kok

# Sonuclari topla, sonra cizdir
python collect_results.py
python plot_results.py --results-dir /veri/kok/results
```

Butun secenekler icin `./run_analysis.sh --help`.

---

## Grafikler

Cizim **iki adimdir**: once toplama, sonra cizme. Ikisini ayirmanin sebebi,
cizim yaparken 105 dizini yeniden taramamak — ve grafik parametreleriyle
oynarken gmx'i yeniden kosmamak.

### Adim 1 — toplama

```bash
python collect_results.py                 # config.sh'teki RESULTS_DIR'e yazar
python collect_results.py -c baska.sh     # baska bir config ile
python collect_results.py -o /tmp/cikti   # baska bir dizine
```

Butun `rep*/analysis/*.xvg` dosyalarini gezer ve iki uzun-format CSV uretir:

| Dosya | Kolonlar |
|---|---|
| `timeseries_long.csv` | `complex, replica, analysis, output, series, time_ps, value, unit` |
| `profile_long.csv` | `complex, replica, analysis, output, series, residue, value, unit` |

Hangi `.xvg`'nin hangi analize ve hangi sekle ait oldugunu `run_analysis.sh
--list` ciktisindan ogrenir — bu eslesme Python'da **tekrarlanmaz**. Manifestoda
olmayan bir `.xvg` atlanir ve stderr'e bir uyari yazilir (dosya basina bir kez).

`series`, `.xvg` icindeki `@ s0 legend` satirlarindan gelir. `gmx rms` legend
yazmadigi icin orada `@ subtitle` kullanilir — yani seri adi
`LIGAND after lsq fit to RECEPTOR_BB` gibi kendini aciklayan bir metin olur.

### Adim 2 — cizme

```bash
python plot_results.py --results-dir /veri/kok/results                  # dordu birden
python plot_results.py --results-dir ... --per-complex                  # yalnizca biri
python plot_results.py --results-dir ... --mean-sd --compare            # secerek
python plot_results.py --results-dir ... --matrix                       # yalnizca matrisler
python plot_results.py --results-dir ... --compare-output rmsf_mhc.xvg  # baska metrik
```

Mod verilmezse **dordu birden** kosar. `-c/--config` ile `COMPLEX_GROUPS` okunur.
`--matrix` tek basina verildiginde uzun-format CSV'ler olmasa da calisir:
matrisler ayri bir urundur ve `cross_rmsd` tek basina kosulmus olabilir.

### Dort mod, dort soru

| Mod | Cikti | Hangi soruya cevap verir |
|---|---|---|
| `--per-complex` | `plots/per_complex/<kompleks>_<cikti>.png` | *Bu simulasyon guvenilir mi?* Kompleks basina tek figur; rep1/rep2/rep3 ayri renkte ust uste. Replikalar ayrisiyorsa yakinsama yok demektir. |
| `--mean-sd` | `plots/mean_sd/<kompleks>_<cikti>.png` | *Yayina/teze ne koyacagim?* Replika ortalamasi + ±SD seridi. Replikalar arasi sacilim belirsizlik bandi olarak gosterilir. |
| `--compare` | `plots/compare_<cikti>.png` | *Hangi peptidler kararli?* Tum kompleksler tek panelde, **medyana gore artan** sirali boxplot. Her kutu bir kompleks, kutu icindeki noktalar replikalarin ortalamalari. |
| `--matrix` | `plots/matrix/<kompleks>_<analiz>.png` ve `..._equilibrium.png` | *Replikalar ayni konformasyonlari mi geziyor?* Kompleks basina N x N isi haritasi izgarasi, ortak renk skalasiyla. Kosegende self-matrisler (tek replika icindeki metastabil durumlar), kosegen disinda capraz ciftler. Koyu bir capraz panel, iki replikanin AYNI bolgeyi ziyaret ettigini soyler. Yaninda **denge egrisi**: her frame'in tum es frame'lere ortalama uzakligi; egri duzlestiginde replika yeni bolge bulmayi birakmistir. Replikalar uc uca eklenmez, ortak zaman ekseninde ust uste cizilir. |

`--compare` varsayilan olarak `rmsd_pep_on_mhc.xvg`'yi kullanir (ana metrik);
`--compare-output` ile baska bir cikti secilebilir.

### Cok serili ciktilar

Bir `.xvg` birden fazla seri tasiyorsa (ör. `gmx gyrate`: Rg, RgX, RgY, RgZ;
`gmx hbond`: hbond sayisi + temas sayisi) her seri ayri cizilir:

- `--per-complex`: **renk replikayi**, **cizgi tipi seriyi** gosterir. Boylece
  uc replika x dort seri tek panelde karisiklik yaratmadan okunur.
- `--mean-sd`: her seri kendi ortalama/SD bandini alir. Farkli fiziksel
  buyuklukler **birbirine karistirilmaz** — aksi halde SD bandi replika
  varyansi degil, seriler arasi fark olurdu ve anlamsiz bir hata cubugu cikardi.

Tek serili ciktilarda gorunum degismez: duz cizgi, efsanede yalnizca replika adi.

### Birimler

CSV'ler **ham** birimdedir (zaman ps, mesafe nm) ve `unit` kolonu gmx'in ne
dedigini tasir. Donusum yalnizca cizim aninda yapilir ve **birime bakarak**:
`nm` ise x10 ve eksen "Å" olur; tanimadigi bir birimde (ör. `nm^2`) deger
**cevrilmez**, eksen ham birimi gosterir ve stderr'e bir uyari yazilir.

Bu bilincli: bilmedigimiz bir birimi cevirmis gibi yapmak, `gmx sasa` (nm²)
eklendigi gun grafikleri sessizce 10 kat yanlis gosterirdi.

### Renkler ve gruplar

Replika renkleri ve grup renkleri, renk korlugu esiklerinden gecen bir paletten
secilmistir (uc replika cizgisi ust uste bindigi icin gerekli).

`COMPLEX_GROUPS` bu projede `("top:top*" "last:last*")`'tir: kompleks adi onekle
basliyorsa karsilastirma panelinde o grubun rengini ve efsane etiketini alir.
**Bos birakilirsa panel tek renkte ve grup efsanesi olmadan cizilir** — baska
bir veri setinde uydurma bir `top*/last*` efsanesi cikmasin diye.

### Bir kompleksin cizimi patlarsa

Cizim dongusu de yalitimlidir: tek bir (kompleks, cikti) ciftinin hatasi
stderr'e adiyla yazilir ve digerleri cizilmeye devam eder. Gozetimsiz bir
kosuda tek bozuk kompleks yuzunden butun grafikleri kaybetmezsiniz.

---

## Yeni analiz ekleme

`analysis/` icine bir `.sh` dosyasi koyun. `run_analysis.sh` degismez,
`collect_results.py` degismez, `plot_results.py` degismez.

```bash
ANALYSIS_NAME="gyrate"
ANALYSIS_DESC="Donme yaricapi zaman serisi"
ANALYSIS_KIND="timeseries"          # timeseries | profile | matrix
ANALYSIS_NEEDS_INDEX=1              # 1 ise index.ndx yoksa once kurulur
ANALYSIS_DEFAULT_BEGIN=0            # ps; -b verilmediginde kullanilir
ANALYSIS_OUTPUTS=(gyrate.xvg)       # idempotency + toplama manifestosu
ANALYSIS_OPTIONAL_OUTPUTS=()        # ilan edilen ama uretilmeyebilen ciktilar

analysis_run() {
    local rep_dir="$1" out_dir="$2"
    # Kullanilabilir degiskenler:
    #   $GMX        gmx ikili dosyasi
    #   $REF_PDB    referans yapi (zincir ID'leri korunmus)
    #   $XTC        trajektori
    #   $NDX        kanonik gruplu index.ndx
    #   $B_PS       etkin equilibration cutoff (ps)
    #   $FORCE      1 ise eklentinin kendi ara artefaktlari da yeniden kurulur
    #   $GROOVE_FIT 1 ise --groove-fit verilmistir (opsiyonel yolu acar)
    "$GMX" gyrate -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/gyrate.xvg" -b "$B_PS" <<< 'LIGAND' \
        || { echo "gmx gyrate basarisiz" >&2; return 1; }
}
```

Bu kadar. `--list` onu otomatik gorur, `run_analysis.sh -a gyrate` calistirir,
`collect_results.py` `timeseries_long.csv`'ye ekler, `plot_results.py` cizer.

### Dikkat edilecek uc sey

**`ANALYSIS_OUTPUTS` KOSULSUZ ilan edilir.** Manifesto (`--list`) `--groove-fit`
gibi secenekleri bilmeyen TAZE bir kabukta okunur; bir ciktiyi kosullu ilan
etmek onu manifestodan dusurur ve `collect_results.py` diskteki dosyayi
atlar. Kosula bagli uretilen ciktilar ayrica `ANALYSIS_OPTIONAL_OUTPUTS`
icinde listelenir: boylece idempotency kontrolu onlari ARAMAZ (yoksa hic
uretilmeyen bir cikti yuzunden analiz her seferinde yeniden kosardi), ama
toplama manifestosunda kalirlar. `rmsf.sh` bunun ornegidir.

**Her gmx cagrisi hata durumunda `return 1` yapmali.** Yoksa runner kismi bir
sonucu basari sayar. Hata mesajlari `>&2`'ye yazilir; runner onlari yakalayip
`run_log.csv`'nin `error` kolonuna koyar.

**Eksik sozlesme sessiz kalmaz.** `ANALYSIS_NAME`, `ANALYSIS_KIND`,
`ANALYSIS_DESC`, `ANALYSIS_DEFAULT_BEGIN`, `ANALYSIS_OUTPUTS` veya
`analysis_run` tanimsizsa o eklenti `HATA` olarak loglanir ve kosu diger
analizlerle devam eder. Dosya adi `_` ile baslayan eklentiler otomatik
kesiften duser (yalnizca `-a _ad` ile acikca secilebilirler).

### `matrix` — 2B ciktilar

`ANALYSIS_KIND="matrix"` ilan eden bir analiz her matris icin **iki** dosya
uretir; `collect_results.py` bunlari
`results/matrices/<kompleks>_<replika>_<cikti>.npz` dosyalarina ve
`results/matrix_summary.csv` ozetine cevirir, `plot_results.py --matrix` de
kompleks basina bir isi haritasi izgarasi + bir denge egrisi cizer. Ilk ornegi
`cross_rmsd`.

| Dosya | Tasidigi | Tasimadigi |
|---|---|---|
| `.xpm` (`gmx -m`) | eksen zamanlari, birim, legend | tam degerler (80 seviyeye yuvarlanmis) |
| `.dat` (`gmx -bin`) | tam `float32` degerler | hicbir metadata -- baslik bile yok |

Bu yuzden **ikisi de zorunlu cikti olarak ilan edilir**: tek baslarina eksiktirler.
Toplama katmani sekli ve eksenleri `.xpm`'den, degerleri `.dat`'tan alir. `.dat`
bozuksa (boyut tutmuyorsa) uyari yazilir ve `.xpm` degerlerine donulur -- daha
kaba, ama dogru bir yedek.

Matrisler uzun-format CSV'ye **girmez**: 451x451'lik 315 matris ~64 milyon satir
ederdi. `timeseries_long.csv`/`profile_long.csv` icin dogru olan bicim 2B veri
icin degil.

> **`-skip` TUZAGI.** `gmx`'in `-skip`/`-skip2` secenekleri `.xpm`'in eksen
> zaman degerlerini bozar: ilk zaman dogru yazilir, geri kalani `0` olur
> (GROMACS 2025.4'te olculdu). Seyreltme **her zaman `-dt`** (ps) ile yapilir.
> `-dt` frame'leri okuma aninda eler, matris zaten seyreltilmis veriden kurulur
> ve eksenler dogru cikar.

#### Bilinen sinirlar

**`min`/`mean`/`max` SEYRELTMEYE baglidir.** Ayni trajektoriden farkli `-dt`
farkli sayilar verir -- daha sik ornekleme daha fazla uc deger yakalar. Bu yuzden
`matrix_summary.csv` kullanilan seyreltmeyi `dt_ps` kolonunda tasir; bu sayilari
bir tabloya koyarken `dt`'yi de yazin.

**Capraz ciftin iki yonu son basamaklarda ayrisir.** `rep1 x rep2` ile
`rep2 x rep1` matematiksel olarak birbirinin transpozesidir, ama `gmx` cifti fit
sirasi ters cevrilmis hesapladigi icin `float32` yuvarlamasi farkli birikir:
gercek veride olculen sapma ~`1e-6` nm (`float32` eps ~`1.2e-7`). Self-matris
boyle bir sapma icermez -- TAM simetriktir ve kosegeni TAM sifirdir.

**`-r/--reps` ile kisitlanmis kosu eksik matris uretir.** Cikti manifestosu
`config.sh`'teki `REPS`'ten turetilir, `-r` ise onu yalnizca kosu icin ezer.
Sessiz bir yanlislik degil: eksik dosya toplanmaz, sonraki tam kosu eksikligi
gorup yeniden uretir.

**DSSP gibi farkli semantikli 2B ciktilar** (residue x zaman) bu izgara
ciziminden faydalanmaz; onlar icin ayri bir cizim bicimi gerekir. Toplama
katmani ise genel: `collect_results.py` replika adiyla eslesmeyen bir matris
ciktisini da toplar, yalnizca `replica_j` kolonu bos kalir.

---

## Birimler ve `-tu` tuzagi

Arac icinde tum zamanlar **ps**, tum mesafeler **nm**'dir. Donusum yalnizca
cizim katmanindadir.

> `gmx`'in `-tu ns` secenegi yalnizca cikti eksenini degil, **`-b`/`-e`
> degerlerini de** ns'e cevirir. Bu yuzden arac icinde `-tu` **hic
> kullanilmaz**; aksi halde `-b 10000` "10 ns sonrasi" yerine "10 000 ns
> sonrasi" anlamina gelir ve gmx uyarmadan bos cikti uretir.

## Durum sozlugu

`results/run_log.csv` icindeki `status`:

| Deger | Anlami |
|---|---|
| `OK` | Analiz kostu, ciktilar uretildi |
| `SKIP_DONE` | Ciktilar zaten vardi, kosulmadi (`--force` ile ezilir) |
| `SKIP_MISSING` | On kosul dosyasi yoktu |
| `HATA` | gmx veya dogrulama hatasi; `error` kolonunda mesaj |

Kolonlar: `timestamp,complex,replica,analysis,status,seconds,begin_ps,error`.

`begin_ps`, o satiri ureten **etkin** equilibration cutoff'udur (`-b` verilmisse
o, yoksa analizin `ANALYSIS_DEFAULT_BEGIN`'i). `SKIP_DONE` ve `SKIP_MISSING`
satirlarinda kullanilacak OLAN degeri tasir. Bu kolon olmadan, farkli `-b` ile
hesaplanmis replikalar birlesik bir CSV'de ayirt edilemezdi.

## Testler

```bash
cd tests && python -m pytest -v            # hizli testler
cd tests && python -m pytest -v -m slow    # tam trajektori gerektirenler
```

Gercek veriye bakan testler `MDKIT_TEST_DATA_ROOT`, `MDKIT_TEST_REP` ve
`MDKIT_TEST_GMX` ortam degiskenleriyle baska bir makineye yonlendirilebilir.
Veri veya gmx yoksa ilgili testler atlanir.

---

## Tasarim belgeleri

`docs/` altinda, aracin neden boyle yazildigini anlatan belgeler var:

- **`docs/design.md`** — tasarim spec'i. GROMACS davranisi hakkinda komutla
  dogrulanmis bulgular (`-tu` tuzagi, hangi tpr'nin zincir ID'lerini tasidigi,
  `gmx rmsf`'in fit davranisi), kanonik grup adlarinin gerekcesi, eklenti
  sozlesmesinin sekli. Araci degistirecek biri once bunu okumali.
- **`docs/implementation-plan.md`** — aracin 11 adimda nasil kuruldugunu
  gosteren uygulama plani. Tarihsel kayit; gunluk kullanim icin gerekmez.
- **`docs/superpowers/specs/2026-09-21-cross-rmsd-matrix-design.md`** — capraz-RMSD
  matrisleri ve `matrix` katman destegi. `.xpm` biciminin uc tuzagi (cok satirli
  eksen yorumlari, ters sirada yazilan piksel satirlari, sabit genislikli
  karakter alani) ve `-skip`/`-b` davranisi komutla dogrulanmis olarak burada.

Bu depo `TUSEB-Bitirme/scratch` deposundan `git subtree split` ile cikarilmistir;
19 commit'lik gecmis korunmustur.
