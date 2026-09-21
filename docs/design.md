# Post-MD Analiz Çerçevesi — Tasarım

**Tarih:** 2026-09-21
**Durum:** Onay bekliyor
**Yerini alacağı:** `mdsimulations/postmd_analysis/rmsd_rmsf.sh`

---

## 1. Amaç

Tek bir trajektori üzerinde interaktif çalışan `rmsd_rmsf.sh`'ı, **35 kompleks × 3 replika = 105
trajektori** üzerinde tekrarlanabilir şekilde koşan, yeni analizlerin dosya ekleyerek takıldığı bir
çerçeveye dönüştürmek; sonuçları tek bir birleşik tablodan çizdirmek.

Tasarımı yönlendiren kısıt: RMSD ve RMSF ilk iki analiz, son ikisi değil. İleride Rg, SASA, hbond,
DSSP, PCA/kovaryans, temas/mesafe analizleri ve equilibration süresi tespiti eklenecek. Bu yüzden
asıl ürün iki script değil, **aralarındaki sözleşme**.

---

## 2. Doğrulanmış mevcut durum

Bu bölümdeki her madde tasarımdan önce komutla doğrulandı (2026-09-21, GROMACS 2025.4-cuda).

| Bulgu | Değer |
|---|---|
| Veri kökü | `/mnt/data/scratch-simulations-TUSEB/` (repo dışında; 940 G kullanımda, 3.4 T boş) |
| Kompleks sayısı | 35 (`*_pandora`), her birinde `rep1`, `rep2`, `rep3` |
| `traj_compact_center_dry.xtc` | **105/105** mevcut |
| `md_0_10.tpr` | **105/105** mevcut |
| `dry_reference.tpr` | **105/105** mevcut (gmx_MMPBSA ürünü) |
| `check_ref.pdb` | **75/105** — 30 replikada EKSİK |
| Trajektori uzunluğu | 10001 frame, 0–100000 ps, dt = 10 ps (last19/rep1'de ölçüldü) |
| Zincirler | A = MHC ağır zincir, B = β2m, C = peptid |

### 2.1 `check_ref.pdb` eksik olan replikalar

`last6`, `last19`, `last20`, `last26`, `last27`, `last36`, `top6`, `top7`, `top8`, `top10` —
her birinde 3 replika, toplam 30.

Bu replikalarda `traj_compact_center_dry.xtc` **mevcut ve sağlam**. Yani `post_md_script.sh`'ın
yalnızca son adımı (frame 0 dump'ı) eksik. Tam script'i yeniden koşmak, zaten var olan kuru
trajektoriyi üretmek için sulu tam sistem üzerinde üç `trjconv` geçişi daha yapar —
30 replika × ~4 GB ara dosya, saatler. Gereksiz.

### 2.2 Zincir ID'leri hangi topolojiden korunur

`gmx trjconv -dump 0` ile referans PDB üretirken:

- `-s dry_reference.tpr` → **zincir ID'leri kaybolur** (tüm atomlar boş chain kolonu).
  `make_ndx`'te `chain A` seçimi boş grup döner.
- `-s md_0_10.tpr` → **zincir ID'leri korunur** (last19/rep1'de A=4343, B=1640, C=220 atom).

Dolayısıyla eksik `check_ref.pdb` üretimi `md_0_10.tpr` ile yapılmalıdır.

### 2.3 Index grup numaraları

Kuru protein sisteminde `make_ndx` varsayılan grupları her zaman 0–9 (System, Protein, Protein-H,
C-alpha, Backbone, MainChain, MainChain+Cb, MainChain+H, SideChain, SideChain-H). `chain A/B/C`
eklenince 10/11/12 olur. Mevcut script'teki `echo -e "10\n12"` bu yüzden doğru çalışıyor —
kırılgan değil. Yine de **isimle seçim** kullanılacak (`chA`, `chC`): numara ne yaptığını söylemiyor.

### 2.4 `gmx rmsf` fit davranışı

`gmx rmsf -[no]fit` varsayılanı **`yes`** ve fit, *seçilen grubun kendi üzerinde* yapılır.
Bu, aynı komutun üç farklı büyüklük ölçebileceği anlamına gelir:

| Seçim | Ölçülen |
|---|---|
| `chC`, `-fit` | Peptidin **iç** esnekliği (oluk içi rijit-cisim hareketi çıkarılmış) |
| Traj önce `chA_BB`'ye fitlenir, sonra `chC` + `-nofit` | Peptidin **oluk çerçevesindeki** esnekliği (kayma/sallanma dahil) |
| `Protein`, `-fit` | Tüm kompleks çerçevesinde peptid RMSF'i |

Mevcut `rmsf_analysis.py` üçüncüsünü hesaplamış (`rmsf_per_position.csv`, 105/105 OK, `-b 10000`).
Bu dosya geçerlidir; yeni çerçeve onu ezmez.

### 2.5 Zaman birimi tuzağı ve isimle grup seçimi

İkisi de komutla doğrulandı:

- **`gmx rms -tu ns`, `-b`/`-e` değerlerini de ns olarak yorumlar.** `-tu ns -e 1` ile
  (tu'suz) `-e 1000` aynı 101 frame'i verdi. Dolayısıyla `-tu ns -b 10000` "10 ns sonrası"
  değil "10 000 ns sonrası" demektir ve gmx uyarı vermeden boş çıktı üretir.
  **Karar: araç içinde `-tu` hiç kullanılmaz, tüm zamanlar ps'dir; ns'e çevirim yalnızca
  çizim katmanında yapılır.**
- **`make_ndx` isimle grup kesişimi kabul eder** (`"chA" & "Backbone"` → `chA_&_Backbone`),
  ve `gmx rms`/`rmsf` grup sorgularına numara yerine isim yazılabilir. Üstelik gmx seçimi
  çıktının başlığına yazar: `@ subtitle "LIGAND after lsq fit to RECEPTOR_BB"`. Yani
  `.xvg` dosyası kendi tanımını taşır.

Doğrulanan atom sayıları (`last10/rep1`): `RECEPTOR_BB` 825 (275 residue),
`LIGAND_BB` 30 (10 residue = IMGQQPAPQV).

---

## 3. Mimari

### 3.1 Dosya yapısı (repo)

Çerçeve, **kendi kendine yeten bir alt dizine** (`mdkit/`) yazılır. Bu dizin ileride ayrı bir
repo'ya çıkarılacak birimdir (§10); mevcut `postmd_analysis/` klasörü ise araç olmaya aday
dosyalarla projeye özel tek seferlik script'leri (`replika_check.py`, `hbond_analysis_sep12.py`,
`pymbar_reanalysis.py`) ve büyük notebook'ları karışık barındırdığı için olduğu gibi çıkarılamaz.

```
mdsimulations/postmd_analysis/          # mevcut karışık klasör — dokunulmuyor
└── mdkit/                              # ← çıkarılabilir birim; gelecekteki repo kökü
    ├── config.sh                       # projeye özel TÜM varsayımlar burada (§3.3)
    ├── run_analysis.sh                 # ÜST SCRIPT: keşif, setup, dispatch, log
    ├── analysis/
    │   ├── lib.sh                      # ortak altyapı
    │   ├── rmsd.sh
    │   └── rmsf.sh                     # ileride: gyrate.sh, sasa.sh, dssp.sh, covar.sh, ...
    ├── collect_results.py              # .xvg -> birleşik uzun-format CSV
    ├── plot_results.py                 # üç grafik modu
    └── README.md                       # kurulum, config, örnek koşu, analiz ekleme
```

`.gitignore` değişikliği gerekmiyor: `!*.sh`, `!*.py`, `!*.md`, `!*/` zaten whitelist'li.

**Sert kural:** `mdkit/` içindeki hiçbir dosyada mutlak yol veya proje adı geçmez. Tek istisna
`config.sh`. Bu kural, laboratuvardan ikinci bir kullanıcının aracı kendi verisiyle çalıştırmasını
mümkün kılan tek şeydir ve §8'de test edilir.

### 3.2 Çıktı yapısı (veri tarafı)

```
$DATA_ROOT/                              # örn. /mnt/data/scratch-simulations-TUSEB
├── <kompleks>/rep{1,2,3}/analysis/
│   ├── index.ndx                    # bir kez kurulur, tekrar kullanılır
│   ├── rmsd_pep_on_mhc.xvg          # ham çıktılar trajektorinin yanında kalır
│   └── ...
└── results/
    ├── timeseries_long.csv          # collect_results.py ürünü
    ├── profile_long.csv
    ├── run_log.csv                  # complex,replica,analysis,status,error,seconds,timestamp
    └── plots/*.png
```

Gerekçe: ham `.xvg` trajektoriyle birlikte taşınır ve hangi koşudan geldiği belirsizleşmez;
birleşik CSV ise çizim katmanının 105 dizini taramasını önler.

### 3.3 `config.sh` — projeye özel varsayımların tek adresi

Aracı "MD analizi" yapan kısım ile "bu pMHC projesi" yapan kısım arasındaki sınır. Başka bir
veri setiyle çalıştırmak = farklı bir `config.sh` vermek; kodda değişiklik gerekmez.

```bash
DATA_ROOT="/mnt/data/scratch-simulations-TUSEB"
GMX="/usr/local/gromacs-2025.4-cuda/bin/gmx"   # boş bırakılırsa PATH'ten bulunur
COMPLEX_GLOB="*_pandora"
REPS=(rep1 rep2 rep3)
TRAJ_NAME="traj_compact_center_dry.xtc"
TPR_NAME="md_0_10.tpr"
REF_NAME="check_ref.pdb"
RESULTS_DIR="$DATA_ROOT/results"
PYTHON="/home/emre/anaconda3/bin/python"        # numpy+pandas+matplotlib+pytest olan yorumlayıcı

# Zincir eşlemesi — analiz script'leri gruplara bu adlarla erişir
CHAIN_RECEPTOR="A"     # MHC ağır zincir
CHAIN_AUX="B"          # beta-2 mikroglobulin
CHAIN_LIGAND="C"       # peptid
```

`run_analysis.sh --config FILE` ile başka bir config verilebilir; verilmezse script'in yanındaki
`config.sh` kullanılır.

### 3.4 Eklenti sözleşmesi

Her `analysis/<ad>.sh` aşağıdakileri tanımlar. `run_analysis.sh` dosyayı `source` edip çağırır.

```bash
ANALYSIS_NAME="rmsd"
ANALYSIS_DESC="Peptid ve MHC RMSD zaman serileri (oluk-uzeri, ic, global)"
ANALYSIS_KIND="timeseries"       # timeseries | profile | matrix
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=0         # ps; --begin ile ezilir
ANALYSIS_OUTPUTS=(rmsd_pep_on_mhc.xvg rmsd_pep_internal.xvg rmsd_complex_bb.xvg)

analysis_run() {
    # $1 = rep_dir, $2 = out_dir
    # lib.sh'tan gelen değişkenler: $REF_PDB $NDX $XTC $B_PS $GMX
    ...
}
```

- `ANALYSIS_KIND` çizim katmanının hangi fonksiyonu kullanacağını belirler. `gyrate` ve `sasa`
  `timeseries` olacağı için yeni çizim kodu gerektirmezler; `dssp` `matrix` olarak yeni bir
  çizim fonksiyonu ister.
- `ANALYSIS_DESC` yalnızca `run_analysis.sh --list` çıktısı için kullanılır.
- `ANALYSIS_OUTPUTS` iki iş görür: idempotency kontrolü (hepsi varsa atla) ve `collect_results.py`
  için dosya manifestosu.
- Çok adımlı analizler (`covar` + `anaeig`) `analysis_run` içinde serbestçe birden fazla gmx
  çağırabilir; sözleşme buna engel değildir.

Yeni analiz eklemek = `analysis/` içine bir dosya koymak. `run_analysis.sh` değişmez.

---

## 4. `analysis/lib.sh` sorumlulukları

0. **Config yükleme** — `--config` ile verilen dosya, yoksa script'in yanındaki `config.sh`
   source edilir. Zorunlu değişkenlerin (`DATA_ROOT`, `TRAJ_NAME`, `TPR_NAME`, `REF_NAME`)
   tanımlı olduğu doğrulanır; eksikse hangi değişkenin eksik olduğunu söyleyerek çıkar.
1. **GROMACS keşfi** — `config.sh`'teki `GMX`; boşsa `PATH`'teki `gmx`; ikisi de yoksa
   anlaşılır hata. Kodda gömülü yol yok.
2. **Hedef keşfi** — verilen dizin bir kompleks mi (içinde `rep*`) yoksa kök mü (`--all` ile
   altındaki `*_pandora`) ayrımı; replika listesi.
3. **Ön koşul kontrolü** — rep başına `traj_compact_center_dry.xtc` ve `md_0_10.tpr`.
4. **`check_ref.pdb` self-heal** — eksikse, koşunun başında **tek seferlik onay sorusu** sorulur:

   ```
   N replikada check_ref.pdb eksik. Mevcut kuru trajektoriden üretilsin mi? [e/H]
   ```

   Onay verilirse her eksik rep için:
   `gmx trjconv -s md_0_10.tpr -f traj_compact_center_dry.xtc -o check_ref.pdb -dump 0 <<< "Protein"`
   Üretim bitince analiz kaldığı yerden devam eder. Onay verilmezse o replikalar `SKIP` olarak
   loglanır, diğerleri koşmaya devam eder.
   `--full-postmd` verilirse bunun yerine `post_md_script.sh` tamamı çalıştırılır (yavaş yol).
   `-y/--yes` soruyu otomatik onaylar.
5. **Index kurulumu** — `rep*/analysis/index.ndx` yoksa `check_ref.pdb`'den kurulur.
   Grup **numarası aritmetiği kullanılmaz**; make_ndx'e isimle kesişim verilir, sonra
   üretilen dosyadaki başlıklar `sed` ile **kanonik adlara** çevrilir:

   ```
   chain $CHAIN_RECEPTOR         -> chA
   chain $CHAIN_AUX              -> chB
   chain $CHAIN_LIGAND           -> chC
   "chA" & "Backbone"            -> chA_&_Backbone
   "chC" & "Backbone"            -> chC_&_Backbone
   ```
   ```
   chA             -> RECEPTOR
   chB             -> AUX
   chC             -> LIGAND
   chA_&_Backbone  -> RECEPTOR_BB
   chC_&_Backbone  -> LIGAND_BB
   ```

   Kanonik adlar iki iş görür: analiz script'leri zincir harflerini bilmeden çalışır
   (config sınırı korunur), ve gmx seçimi `.xvg` başlığına yazdığı için çıktı kendi
   tanımını taşır. Kurulum sonrası beş kanonik adın da dosyada olduğu doğrulanır.
   Kurulumdan sonra üretilen grup isimleri doğrulanır (`grep '^\[' index.ndx`); beklenen isim
   yoksa o replika hata olarak loglanır.
6. **Idempotency** — `ANALYSIS_OUTPUTS`'un tamamı varsa ve boyutu > 0 ise analiz atlanır;
   `--force` bunu ezer.
7. **Hata yalıtımı** — her (replika × analiz) çağrısı alt kabukta ve hata tuzağıyla koşar.
   Tek bir replikanın hatası batch'i düşürmez; `run_log.csv`'ye `HATA` + mesaj yazılıp devam edilir.
8. **Log** — `results/run_log.csv`'ye satır ekler; ekranda kompakt ilerleme yazar.
   `status` sözlüğü tek anlamlıdır:

   | Değer | Anlamı |
   |---|---|
   | `OK` | Analiz koştu, çıktılar üretildi |
   | `SKIP_DONE` | Çıktılar zaten vardı (idempotency), koşulmadı |
   | `SKIP_MISSING` | Ön koşul dosyası yoktu (ör. onaylanmamış `check_ref.pdb`) |
   | `HATA` | gmx veya doğrulama hatası; `error` kolonunda mesaj |

---

## 5. `run_analysis.sh` arayüzü

```
run_analysis.sh [seçenekler] <hedef>

<hedef>   Bir kompleks dizini (içinde rep1/rep2/rep3),
          veya --all ile birlikte kök dizin.

Seçenekler:
  -c, --config FILE     Config dosyası (varsayılan: script'in yanındaki config.sh)
  -a, --analysis LIST   Virgülle ayrılmış analiz adları (varsayılan: tümü). Örn: rmsd,rmsf
      --all             Hedefi kök kabul et, altındaki tüm *_pandora dizinlerini gez
  -r, --reps LIST       Replika listesi (varsayılan: rep1,rep2,rep3)
  -b, --begin PS        Equilibration cutoff (ps); analiz varsayılanını ezer
      --groove-fit      MHC-çerçeveli fitlenmiş trajektori üret, chC RMSF'ini onunla hesapla
      --force           Çıktı varsa bile yeniden hesapla
      --full-postmd     Eksik check_ref.pdb için post_md_script.sh'ın tamamını koş
  -y, --yes             Onay sorularını otomatik onayla (batch için)
  -l, --list            Mevcut analizleri ve açıklamalarını listele, çık
  -n, --dry-run         Ne yapılacağını yazdır, çalıştırma
  -h, --help
```

**Argümansız çağrı:** stdin bir TTY ise, mevcut script'in menü davranışı korunur — çalışılacak
dizin ve analizler sorulur. TTY değilse kullanım özeti yazdırılıp çıkılır. Böylece interaktif
alışkanlık kaybolmaz ama batch'i bloklamaz.

---

## 6. Analiz tanımları

### 6.1 `analysis/rmsd.sh`

`KIND=timeseries`, `DEFAULT_BEGIN=0` — **tüm trajektori**. Equilibration süresi bu eğrilerden
tespit edileceği için baştan kesilmemeli.

| Çıktı | fit grubu | hesap grubu | Ne söyler |
|---|---|---|---|
| `rmsd_pep_on_mhc.xvg` | `RECEPTOR_BB` | `LIGAND` | Peptidin oluk içindeki toplam hareketi (ana metrik) |
| `rmsd_pep_internal.xvg` | `LIGAND_BB` | `LIGAND` | Peptidin konformasyon değişimi |
| `rmsd_complex_bb.xvg` | `Backbone` | `Backbone` | TÜM KOMPLEKS backbone RMSD'si (ağır zincir + β2m + peptid, 385 residue) — global denge — equilibration tespitinin dayanağı |

> **Düzeltme (final review):** bu çıktının adı §6.1'de önce `rmsd_complex_bb.xvg` idi, ama
> hem fit hem hesap grubu `Backbone`'dur; gerçek indekste `Backbone` = 1155 atom /
> 385 residue (ağır zincir + β2m + **peptid**), `RECEPTOR_BB` ise 825 / 275. Yani dosya
> MHC omurgasını değil, tüm kompleksin omurgasını taşıyordu — tablodaki "global denge"
> açıklaması doğru, ad yanlıştı. Gruplar korunup ad `rmsd_complex_bb.xvg` yapıldı.
> (Belgenin başka yerlerinde geçen eski ad tarihsel kayıt olarak bırakıldı.)

İlk iki eğrinin farkı yaklaşık olarak rijit-cisim kayma/sallanmadır: tek bir RMSD eğrisi
"peptid kaydı mı, büküldü mü" sorusunu ayıramaz; bu ikisi ayırır.

Komut kalıbı (`-tu` **yok**, zaman ps — §2.5):
```bash
"$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" -o "$2/rmsd_pep_on_mhc.xvg" \
    -b "$B_PS" <<< $'RECEPTOR_BB\nLIGAND'
```

### 6.2 `analysis/rmsf.sh`

`KIND=profile`, `DEFAULT_BEGIN=10000` — ilk 10 ns atılır, mevcut `rmsf_analysis.py` ile tutarlı.

| Çıktı | Nasıl | Koşul |
|---|---|---|
| `rmsf_chC_self.xvg` | `LIGAND` + `-fit` + `-res` | Varsayılan |
| `rmsf_chA.xvg` | `RECEPTOR` + `-fit` + `-res` | Varsayılan |
| `rmsf_chC_groovefit.xvg` | Traj `RECEPTOR_BB`'ye fitlenir → `LIGAND` + `-nofit` + `-res` | `--groove-fit` |

`--groove-fit` yolu rep başına ek bir `trjconv -fit rot+trans` geçişi ve ~230 MB ara dosya
(`traj_fit_mhc.xtc`) üretir; 105 replikada ~24 GB ve birkaç saat. Disk müsait (3.4 T boş).
Üretilen fitlenmiş trajektori **saklanır**, çünkü ileride PCA/DCCM aynı dosyayı isteyecek.

### 6.3 `rmsf_analysis.py` ile ilişki

Bu sürümde ikisi ayrı durur. `rmsf_analysis.py` peptid RMSF'inin (Protein-fit tanımı) tek kaynağı
olmaya devam eder; yeni `rmsf.sh` farklı tanımlar üretir ve `.xvg`'leri saklar. Birleştirme kararı
ileriye bırakılmıştır.

---

## 7. Toplama ve çizim

### 7.1 `collect_results.py`

`rep*/analysis/` altındaki `.xvg` dosyalarını gezer, `#` ve `@` başlıklarını atlar,
`@ s<N> legend "..."` satırlarından seri isimlerini okur. Çıktı, uzun format:

```
complex, replica, analysis, output, series, x, y
```

`KIND`'a göre ayrı dosyalara yazılır: `results/timeseries_long.csv` ve `results/profile_long.csv`.
Bu sürümde `matrix` türünde analiz yok; eklendiğinde uzun format ona uymayacağı için kendi
dosya biçimini getirecek (§9).
Birim dönüşümü burada yapılır ve kolon adında belirtilir: RMSD/RMSF nm → Å, zaman ps → ns
(§2.5 gereği `.xvg`'ler ps ve nm cinsindendir).

### 7.2 `plot_results.py`

| Mod | Ne üretir |
|---|---|
| `--per-complex` | Kompleks başına bir figure; rep1/rep2/rep3 farklı renkte üst üste. Yakınsama denetimi. |
| `--mean-sd` | 3 replika ortalaması + ±SD şeridi. Yayına/teze giden temiz figure. |
| `--compare` | Tüm kompleksler tek panelde; ortalama peptid RMSD'sine göre sıralı boxplot, `top*` vs `last*` renk ayrımı. |

Çizim kodu `KIND`'a göre üç fonksiyona dallanır (timeseries / profile / matrix), analiz adına göre
değil. Böylece `gyrate` ve `sasa` eklendiğinde yeni çizim kodu gerekmez.

Renk paleti ve eksen/legend kuralları `dataviz` rehberinden alınır; 35 kategorili boxplot'ta
okunabilirliği belirleyen budur.

---

## 8. Doğrulama planı

1. **`--dry-run --all`** — 35 kompleks × 3 replika keşfi doğru mu, eksik dosyalar doğru raporlanıyor mu.
2. **Tek kompleks smoke testi** (`last10`, 3 replika, `rmsd,rmsf`):
   - `rmsd_*.xvg` satır sayısı = 10001 (dt = 10 ps, 0–100 ns)
   - `rmsf_chC_self.xvg` satır sayısı = peptid residue sayısı (`last10` → 10)
   - `rmsf_chA.xvg` satır sayısı = MHC ağır zincir residue sayısı
3. **Bilinen sonuçla çapraz kontrol** — `last10/rep1` üzerinde `gmx rmsf -b 10000` `Protein`
   grubuyla koşulur, son 10 residue `rmsf_per_position.csv`'deki `last10,rep1` satırlarıyla
   karşılaştırılır. Üç ondalık basamağa kadar eşleşmeli. Bu, çerçevenin gmx çağrı katmanını
   bilinen-doğru bir çıktıya karşı doğrular.
4. **Idempotency** — aynı komut ikinci kez koşulduğunda tüm analizler `SKIP` olmalı, süre ~0.
5. **Hata yalıtımı** — bir replikanın `index.ndx`'i kasten bozulur; batch'in devam ettiği ve
   `run_log.csv`'de tek bir `HATA` satırı olduğu doğrulanır.
6. **`check_ref.pdb` self-heal** — eksik 30 replikadan biri üzerinde onay akışı denenir; üretilen
   PDB'nin zincir dağılımı (A/B/C) ve atom sayısı kontrol edilir.
7. **Taşınabilirlik** — `grep -rE '/mnt/|/usr/local/|/home/' mdkit/ --exclude=config.sh` **boş
   dönmeli**. Ek olarak, `DATA_ROOT`'u tek bir komplekse işaret eden bir kopya config ile
   `--dry-run` koşulur; aracın config dışında hiçbir yere bağlı olmadığı doğrulanır.

---

## 9. Kapsam dışı (bu sürümde yapılmayacak)

- Rg, SASA, hbond, DSSP, PCA/kovaryans, temas/mesafe analizleri — sözleşme bunları kaldıracak
  şekilde tasarlandı, ama bu sürümde yazılmıyor.
- Equilibration süresinin **otomatik** tespiti. Bu sürüm gerekli veriyi üretir
  (`rmsd_complex_bb.xvg`, tüm trajektori); tespit algoritması ayrı bir analiz olarak sonra gelir.
- `rmsf_analysis.py`'nin emekliye ayrılması veya çıktısının yeniden üretilmesi.
- Paralel çalıştırma. 105 replika seri koşulur; gerekirse `--reps` ile bölünüp elle paralelleştirilir.

---

## 10. Ayrı repo'ya çıkarma planı

Karar: **şimdilik ayrılmıyor.** Gerekçe — geliştirme sırasında araç ve analiz eş zamanlı
değişecek, iki repo her değişikliği iki commit'e böler. Ancak laboratuvardan ikinci kullanıcılar
beklendiği için ayrılma *olacak* kabul edilip maliyeti sıfıra indirilmiştir:

- Tüm projeye özel varsayımlar `config.sh`'te (§3.3), kodda mutlak yol yok (§3.1, §8.7).
- `mdkit/` kendi `README.md`'si ile kendi kendine yeten bir birim.

Ayırma zamanı geldiğinde:

```bash
git subtree split -P mdsimulations/postmd_analysis/mdkit -b mdkit
# -> mdkit branch'i yalnızca o dizinin geçmişini taşır; yeni repo'ya push edilir
```

**Ayırmayı tetikleyecek sinyaller:** laboratuvardan ikinci bir kişi/proje aracı kullanmaya
başlarsa, tezde veya makalede yazılım olarak atıf gerekirse (Zenodo DOI + sürüm etiketi),
ya da araç kendi sürüm numarasını hak edecek kadar oturursa.

## 11. Açık notlar

- Replikaların `md_0_10.tpr`'si zincir bilgisini taşıyan tek güvenilir topoloji; `dry_reference.tpr`
  bu iş için kullanılamaz (§2.2). İleride gmx_MMPBSA çıktılarıyla çalışan bir analiz eklenirse
  bu ayrım tekrar önem kazanacak.
- `set -e` ile 105 elemanlı döngü birlikte tehlikelidir: tek bir kötü replika saatlerce süren
  batch'i düşürür. Hata yalıtımı (§4.7) bu yüzden opsiyonel değil.
