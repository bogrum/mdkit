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
        ├── matrix/<kompleks>_<analiz>.png            (kompleks basina skala)
        ├── matrix/<kompleks>_<analiz>_equilibrium.png
        ├── profile_compare/<cikti>.png               (gruplar arasi, konuma gore)
        ├── matrix_global/<kompleks>_<analiz>.png     (ortak skala)
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
- Python 3.11+ ve numpy, pandas, matplotlib, scipy (cizim ve istatistik icin);
  pytest (testler icin).
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

# Acik yuzey (SASA) -- pahali, seyreltme varsayilani 100 ps
./run_analysis.sh --all -a sasa /veri/kok
SASA_DT=200 ./run_analysis.sh --all -a sasa /veri/kok

# Peptid-MHC hidrojen baglari
./run_analysis.sh --all -a hbond /veri/kok

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
python plot_results.py --results-dir /veri/kok/results                  # besi birden
python plot_results.py --results-dir ... --per-complex                  # yalnizca biri
python plot_results.py --results-dir ... --mean-sd --compare            # secerek
python plot_results.py --results-dir ... --matrix                       # yalnizca matrisler
python plot_results.py --results-dir ... --profile-compare              # gruplar arasi profil
python plot_results.py --results-dir ... --compare-output rmsf_mhc.xvg  # baska metrik
```

Mod verilmezse **besi birden** kosar. `-c/--config` ile `COMPLEX_GROUPS` okunur.
`--matrix` tek basina verildiginde uzun-format CSV'ler olmasa da calisir:
matrisler ayri bir urundur ve `cross_rmsd` tek basina kosulmus olabilir.

### Bes mod, bes soru

| Mod | Cikti | Hangi soruya cevap verir |
|---|---|---|
| `--per-complex` | `plots/per_complex/<kompleks>_<cikti>.png` | *Bu simulasyon guvenilir mi?* Kompleks basina tek figur; rep1/rep2/rep3 ayri renkte ust uste. Replikalar ayrisiyorsa yakinsama yok demektir. |
| `--mean-sd` | `plots/mean_sd/<kompleks>_<cikti>.png` | *Yayina/teze ne koyacagim?* Replika ortalamasi + ±SD seridi. Replikalar arasi sacilim belirsizlik bandi olarak gosterilir. |
| `--compare` | `plots/compare_<cikti>.png` | *Hangi peptidler kararli?* Tum kompleksler tek panelde, **medyana gore artan** sirali boxplot. Her kutu bir kompleks, kutu icindeki noktalar replikalarin ortalamalari. |
| `--profile-compare` | `plots/profile_compare/<cikti>.png` | *Gruplar birbirinden ayriliyor mu, ve NEREDE?* Her `profile` ciktisi icin konum kutulari (`P1`, `P2`, `orta`, `PO-1`, `PO`) ve her kutuda grup basina boxplot + test p degeri. Gruplar `COMPLEX_GROUPS`'tan gelir. |
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

### Konum normalizasyonu: neden gerekli

Peptidler ayni uzunlukta degildir (bu veri setinde 8-11 residue). Bu yuzden
**ham residue numarasi kompleksler arasi karsilastirilamaz**: bir 8-mer'in 5.
residue'su peptidin ortasindadir, 11-mer'in 5.'si degil. Residue 5'in RMSF'ini
35 komplekste ortalamak, farkli rolleri olan residue'lari ayni kovaya koymaktir.

`profile_long.csv` bu yuzden **iki** konum kolonu tasir:

| Kolon | Sayma yonu | 9-mer'de ornek |
|---|---|---|
| `residue` | N-ucundan, 1'den baslar | 1, 2, 3, … 9 |
| `residue_from_end` | C-ucundan, son residue 0 | -8, -7, … -1, 0 |

Ikisi birlikte her turlu konum normalizasyonunu turetilebilir kilar.
**Etiket (`P2`, `PO` gibi) CSV'ye YAZILMAZ**: o bir yorumdur ve katman
sozlesmesi geregi cizim tarafina aittir. Baska bir semaya (or. yalnizca ankraj
/ ankraj-disi) gecmek isteyen biri CSV'yi yeniden uretmek zorunda kalmaz.

`plot_results.py` bu iki kolondan su kutulari turetir:

    P1  P2  orta  PO-1  PO

N-ucu **oncelikli**: 3 residue'luk bir peptidde 2. residue hem `P2` hem
`PO-1` olurdu ve ayni olcum iki kutuya birden girerdi.

> **Figurde uzunluk dagilimi HER ZAMAN yazar.** Basligin ucuncu satiri her
> grubun profil uzunlugu medyanini, araligini ve kompleks sayisini gosterir.
> Bu bir esige BAGLANMAZ ("medyanlar X'ten fazla farkliysa uyar" gibi):
> tam esikte olan bir vaka sessizce gecerdi. Karar okuyucunun.
>
> Gerekce olculdu. Bu veri setinde konuma dayali iki "anlamli" sonuc
> cikti ve **ikisi de uzunluk farkindan** geliyordu:
>
> | bulgu | ham p | karisan | uzunluga gore tabakali |
> |---|---|---|---|
> | peptid orta bolge | 0.023 | peptid uzunlugu (rho +0.29) | 9-mer'de p=0.955 |
> | MHC son residue | 0.005 | MHC zincir uzunlugu (rho **+0.888**) | 275'te 0.234, 276'da 1.000 |
>
> Uclardan sayilan konumlar (`PO`, `PO-1`) uzunluga ozellikle duyarlidir:
> "son residue" tanim geregi sistemin BITTIGI yerdir, her komplekste ayni
> residue degil.

> **`--profile-compare` coklu test uretir.** Bes konum test edilir; figurun
> basliginda Bonferroni esigi (`0.05/5 = 0.01`) yazar ve esigi gecmeyen p
> degerleri **soluk** cizilir. Tek bir `p < 0.05` gorup anlamli saymayin --
> figur bunu kasten zorlastirir.

> **Uzun profillerde az bilgilendiricidir.** `rmsf_mhc` ~275 residue'dur ve
> "orta" kutusu bunlarin neredeyse tamamini icerir. Figur yine de uretilir:
> sihirli bir uzunluk esigiyle sessizce gizlemek, zayif bir ciktidan kotudur.

Replikalar **once kompleks basina ortalanir**, sonra gruplar karsilastirilir.
Aksi halde ayni kompleksin uc replikasi uc bagimsiz gozlem sayilir ve p
degeri oldugundan kucuk cikar.

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

## Iki gmx tuzagi: `hbond` ve `sasa`

Ikisi de olcumle saptandi ve eklentilerin icinde yorum olarak yazili.

### `gmx hbond` referans PDB ile SESSIZCE sifir dondurur

Ayni komut, yalnizca `-s` degisiyor:

| `-s` | sonuc |
|---|---|
| `check_ref.pdb` | ortalama **0.00** H-bagi |
| `md_0_10.tpr` | ortalama **10.96** H-bagi |

PDB'de bag bilgisi yoktur, bu yuzden donor/akseptor cikarilamaz -- ve arac
**hata vermez**. Yani yanlis dosyayi vermek sessiz veri kaybidir. `hbond.sh`
bu yuzden `$REF_PDB` degil `"$rep_dir/$TPR_NAME"` kullanir, ve
`tests/test_hbond.py` hem kaynagi hem gercek kosuyu denetler.

Topoloji eklentiye degisken olarak verilmez ama `$TPR_NAME` ve `$rep_dir`
gorunurdur (`mdkit_run_isolated` alt kabukta kosturur) -- `cross_rmsd`'nin
kardes replikalara eristigi mekanizmanin aynisi. `run_analysis.sh` degismedi.

Hedef grup `RECEPTOR`'dur, `AUX` (b2m) degil: **bes farkli komplekste de**
b2m'nin peptide H-bagi katkisi tam sifir olctuldu (b2m olugun karsi
tarafindadir). `RECEPTOR` ile `RECEPTOR`+`AUX` ayni sonucu verir.

### `gmx hbond -o` verilmezse CALISMA DIZININE yazar -- ve 99'da durur

`gmx hbond`'un `-num` disinda bir de `-o` (H-bagi index dosyasi) ciktisi
vardir. Verilmezse `hbond.ndx`'i **calisma dizinine** yazar, yani
`run_analysis.sh`'in cagrildigi yere. Dahasi GROMACS her kosuda eskisini
`#hbond.ndx.N#` olarak yedekler ve **99. yedekte durur**:

    Will not make more than 99 backups

105 replikalik gercek bir kosuda bu yasandi: ilk 84 replika calisti, **son
21'i HATA verdi** ve depo 98 adet yedek dosyayla doldu. Ayni tuzak spec'te
`gmx rms -o` icin zaten belgelenmisti; `hbond.sh`'e uygulanmasi atlanmisti.

Cozum `cross_rmsd`'nin kullandigi kalibin aynisi: `-o` bir `mktemp -d`
dizinine yonlendirilir ve `analysis_run` cikista siler. `out_dir`'e yazmak
da olmaz -- manifestoda ilan edilmeyen bir `.ndx` her toplamada uyari
uretirdi.

Iki test birden baglar: kaynakta `-o` ve `mktemp` aranir, ve kosu **bos bir
dizinde** calistirilip dizinin bos kaldigi dogrulanir. `.gitignore`'daki
`#*#` satiri yalnizca ikinci savunma hattidir.

### `gmx sasa -or` `profile` kind'i ile UYUMSUZ

`-or` (residue basina alan), yuzey secimi ile `-output` secimini **ayni
dosyada arka arkaya** yazar ve residue numaralari zincir basina sifirlanir.
Olculen ornek: **386 satir, 100 TEKRARLI residue numarasi**, en buyuk numara
276, son satir residue 10.

`collect_results.py`'nin `profile` kind'i residue'yu benzersiz varsayar; bu
dosya onu **sessizce** bozardi. Bu yuzden `sasa.sh` yalnizca `-o` (zaman
serisi) uretir; residue bazli SASA icin blok-farkinda bir ayristirici gerekir
ve bu henuz yazilmadi.

### `sasa` pahalidir, seyreltme varsayilan olarak aciktir

Olculdu (9001 frame, tek replika):

| ayar | sure | peptid SASA ortalamasi |
|---|---|---|
| tam (10 ps) | 190 s | 5.5857 nm² |
| `-dt 100` (varsayilan) | **19.7 s** | 5.5841 nm² |
| `-dt 200` | 10.3 s | 5.5737 nm² |

105 replika x 2 kosu: tam cozunurlukte ~11 saat, `-dt 100` ile ~70 dakika.
Ortalamaya etkisi %0.03. `SASA_DT` ile ezilir.

`hbond` seyreltilmez: tam trajektoride replika basina ~2 saniye surdugu icin
kazanci yok.

### Gomulu yuzey saklanmaz, turetilir

`sasa` iki cikti verir: peptidin **kompleks icindeki** acik yuzeyi
(`sasa_pep_in_complex.xvg`, `-surface Protein -output LIGAND`) ve **tek
basina** acik yuzeyi (`sasa_pep_alone.xvg`, `-surface LIGAND`). Gomulu yuzey
bu ikisinin farkidir ve ucuncu bir dosya olarak yazilmaz -- saklanan iki
buyuklugun farkini dosyaya dokmek veriyi cogaltmak olurdu.

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

`.xpm`'in kendi `title`'i yalnizca NE olculdugunu soyler (`LIGAND_BB RMSD
matrix`); NEYE FIT EDILDIGINI ise `gmx` sadece `-o` ile yazdigi `.xvg`'nin
`subtitle`'ina koyar (`LIGAND_BB after lsq fit to RECEPTOR_BB`) ve o `.xvg`
atilir. Eklenti bu satiri `.xpm`'e ikinci bir yorum olarak **tasir**; figurler
onu gosterir. Boylece matris dosyasi kendi kendini aciklar ve cizim katmani gmx
grup adlarini bilmek zorunda kalmaz -- yalnizca dosyada yazani gosterir.
Alternatifi, `plot_results.py`'ye `RECEPTOR_BB` string'ini gommekti; bu,
katmanlar arasi sozlesmeyi bozardi.
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

**`matrix_summary.csv`'deki `mean_vs_self` kolonu**, "replikalar ayni
konformasyonel alani mi ornekliyor?" sorusunun nicel cevabidir:

    mean_vs_self = mean / ((self_i + self_j) / 2)

1'e yakin = evet (capraz uzaklik, replika ICI uzakliktan farkli degil);
buyudukce = hayir. Isi haritasina bakip goz karariyla karar vermenin yerini
alir -- renk skalasi kompleksten komplekse degistigi icin goz karari zaten
guvenilmezdir. Oran CIFT basinadir: gercek veride heterojenligin cogunlukla
komplekse yayilmadigi, TEK bir replika ciftinden geldigi gorulur. Self
satirlarda tanim geregi 1.0; referans self matrisi yoksa bos birakilir.
Self ortalamalari kosegen haric hesaplandigi icin payda bir miktar buyuktur,
yani oran muhafazakardir.

**Renk skalasi IKI sette birden uretilir.** `plots/matrix/` her kompleksi
KENDI araligiyla cizer (kompleks ICI kontrast korunur); `plots/matrix_global/`
ayni analizin butun komplekslerini ORTAK bir aralikla cizer (kompleksler ARASI
kiyas mumkun olur). Ikisi de gerekli, cunku takas gercek: `top2`'nin araligi
0-4.4 A, `last10`'unki 0-14.9 A; ortak skalada `top2` neredeyse duz cikar ama
"bu kompleks digerlerine gore dusuk" bir bakista gorunur. Her iki setin de
BASLIGINDA kullanilan aralik yazar -- yazmasaydi okuyucu renkleri kompleksler
arasi kiyaslamaya kalkardi ve yanilirdi (olculdu: `top1`'de "yesil" 2.9 A,
`last10`'da 8.9 A). Tek kompleksli bir (analiz, birim) grubunda ortak skala
kendi skalasiyla ayni oldugu icin ikinci set yazilmaz.

**Denge egrisinin penceresi SABIT bir suredir** (`--matrix-smooth-ns`,
varsayilan 1.0 ns; `0` kapatir). Serinin oranina baglanmis bir pencere
(`n/20` gibi) filtrenin kesme frekansini kosu uzunluguna baglar: ayni sistem
20 ns yerine 100 ns kosuldugunda ayni surec farkli duzlestirilir. 1.0 ns,
gercek veride olculen otokorelasyon surelerinin (`tau_int` 0.5-7.8 ns, medyan
~3.9 ns) hemen hepsinin altinda kalir. Pencere ORTALANIR ve ham egri her zaman
figurde kalir.

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
