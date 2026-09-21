# Çapraz-RMSD Matrisleri ve `matrix` Katman Desteği — Tasarım

**Tarih:** 2026-09-21
**Durum:** Onay bekliyor
**Kapsam:** `analysis/cross_rmsd.sh` (yeni), `collect_results.py`, `plot_results.py`, `README.md`

---

## 1. Amaç

Replikalar arası **2D RMSD matrisi** üretmek ve bunu aracın mevcut boru hattına
uçtan uca bağlamak.

Cevaplanan soru: *aynı kompleksin üç replikası aynı konformasyonel bölgeleri mi
geziyor, yoksa her biri kendi havzasında mı kalıyor?* Mevcut `rmsd` analizi bunu
cevaplayamaz — her replikayı kendi referansına göre ölçtüğü için iki replika aynı
RMSD değerine **farklı** konformasyonlarda ulaşabilir.

Bu, `README.md`'de "bilinen sınır" olarak yazılı olan `ANALYSIS_KIND="matrix"`
boşluğunu da kapatır: `matrix` bugün ilan edilebiliyor ama hiçbir katman
tüketmiyor.

---

## 2. Doğrulanmış bulgular

Bu bölümdeki her madde tasarımdan önce komutla doğrulandı
(2026-09-21, GROMACS 2025.4-cuda, `last10_IMGQQPAPQV_A0201_pandora`).

### 2.1 `-b` ikinci trajektoriye de uygulanır

`gmx rms -h` çıktısında `-b2`/`-e2` yoktur, yalnızca `-skip2` vardır. Bu,
`-b`'nin `-f2` ile verilen trajektoriye uygulanmadığı izlenimini verir; öyle olsaydı
`-b 10000` dengelenmiş frame'leri eşin **dengelenme dönemi dahil** tamamına karşı
kıyaslar, sonuç sessizce anlamsız çıkardı.

Ölçüm bunun tersini gösterdi. Trajektori 10001 frame (0–100000 ps, dt = 10 ps):

```
-b 90000 -skip 200 -skip2 100  ->  "Building RMSD matrix, 6x11 elements"
```

`-b` yalnızca `-f`'e uygulansaydı kolon sayısı 10001/100 ≈ 101 olurdu; 11 çıktı,
yani 1001/100. **`-b` her iki trajektoriyi de kırpar.**

**Sonucu:** eşin trajektorisini `trjconv -b` ile önceden kırpma, bu kırpılmış
dosyayı saklama, kesim değerini dosya adına gömerek önbelleği geçersizleştirme —
bunların hiçbirine gerek yok. Tasarımdan çıkarıldılar.

### 2.2 `-skip`, `.xpm` eksen zamanlarını bozar

`-skip > 1` verildiğinde `.xpm`'in eksen yorumlarında yalnızca ilk zaman doğru
yazılır, geri kalanı sıfır olur:

```
-b 90000 -skip 200        ->  /* x-axis:  90000 0 0 0 0 0 */
-b 99000 (skip yok)       ->  /* x-axis:  99000 99010 99020 ... 100000 */
```

`-f2` masum; `-skip` tek başına da bozuyor. Eksen yorumuna güvenen bir
ayrıştırıcı bütün frame'leri t = 0'da sanırdı.

**Sonucu:** seyreltme `-skip`/`-skip2` ile değil, **`-dt` (ps)** ile yapılır.
`-dt` frame'leri okuma anında eler, matris zaten seyrekleştirilmiş veriden
kurulur ve eksenler doğru çıkar:

```
-b 90000 -dt 1000  ->  "11 11", x-axis: 90000 91000 ... 100000
```

`-tu` ise araç genelinde olduğu gibi burada da kullanılmaz (`README.md`, `-tu`
tuzağı: `-b`/`-e` değerlerini de çevirir).

### 2.3 Maliyet

`-b 10000 -dt 200` ile çapraz bir çift (rep1 × rep2):

| Ölçüm | Değer |
|---|---|
| Matris | 451 × 451 |
| Süre | 2.9 s |
| `.xpm` boyutu | 209 KB |

35 kompleks × 9 matris = 315 çağrı ≈ **16 dakika**, ≈ **66 MB**. Maliyet tasarımı
kısıtlamıyor; `-dt 200` seçimi hassasiyet değil okunabilirlik gerekçesiyle yapıldı.

### 2.4 `.xpm` biçimi

```
/* legend:  "RMSD (nm)" */        <- birim buradan okunur
static char *gromacs_xpm[] = {
"451 451   80 1",                <- genislik yukseklik renk_sayisi karakter/piksel
"A  c #FFFFFF " /* "0.565" */,   <- karakter -> deger tablosu (80 satir)
...
/* x-axis:  90000 91000 ... */   <- UZUNSA BIRDEN FAZLA SATIRA BOLUNUR
/* y-axis:  90000 91000 ... */
"Hk1TeU",                        <- piksel satirlari
```

Üç tuzak, üçü de ölçümle saptandı:

1. **Eksen yorumları bölünür.** 101 elemanlı bir eksen iki ayrı
   `/* x-axis: ... */` satırına yazılır. Ayrıştırıcı aynı etiketli satırları
   **birleştirmelidir**; yalnızca ilkini okuyan bir ayrıştırıcı ekseni sessizce
   keser.

2. **Piksel satırları y ekseninin TERSİ sırada yazılır.** Self-matriste köşegen
   tanım gereği sıfırdır (`"A  c #FFFFFF " /* "0" */`). 9×9 bir self-matriste
   `A`'nın konumu ölçüldü:

   | Dosya satırı | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
   |---|---|---|---|---|---|---|---|---|---|
   | `A` kolonu   | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |

   Sıfır köşegen dosyada **ters köşegen** üzerinde. Yani dosyanın ilk satırı
   y ekseninin **son** değerine karşılık gelir. Ayrıştırıcı satır sırasını
   çevirmelidir (`rows[::-1]`), aksi halde matris yatay eksende aynalanır —
   ve bu, self-matris dışında **gözle fark edilmez**.

3. **Genişlik = x ekseni = `-f`; yükseklik = y ekseni = `-f2`.** `"6 11"`
   başlığı 6 kolon × 11 satır demektir; doğrulandı (`-skip 200`/`-skip2 100`
   ile 1001/200 = 6 ve 1001/100 = 11).

### 2.5 Değerler 80 seviyeye yuvarlanmıştır

`.xpm` sürekli değerleri `-nlevels` (varsayılan 80) renk seviyesine indirger.
`-bin` ile ham binary dump alınabilir. **Karar: `.xpm` yeterli.** Isı haritası
için 80 seviye fazlasıyla yeterli; `matrix_summary.csv`'deki min/mean/max
sayılarının da bu çözünürlükte olduğu belgelenir.

---

## 3. Katman 1 — `analysis/cross_rmsd.sh`

### Sözleşme

```bash
ANALYSIS_NAME="cross_rmsd"
ANALYSIS_KIND="matrix"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
# $REPS'ten turetilir: cross_rmsd_rep1.xpm cross_rmsd_rep2.xpm cross_rmsd_rep3.xpm
ANALYSIS_OUTPUTS=( ... )
ANALYSIS_OPTIONAL_OUTPUTS=()
```

`rep_i` bağlamında koşan eklenti, `$REPS`'teki **her** replikaya karşı
(kendisi dahil) bir matris üretir. Fit `RECEPTOR_BB`, RMSD `LIGAND_BB`.

### Kardeş replikalara erişim

Eklenti tek bir `rep_dir` alır, ama `mdkit_run_isolated` `analysis_run`'ı ana
kabuğun **alt kabuğunda** koşturur (`analysis/lib.sh:222`). Bu yüzden `$REPS`,
`$TRAJ_NAME`, `$REF_NAME` görünürdür ve kardeşler `dirname "$rep_dir"` ile
bulunur. **`run_analysis.sh` değişmez.**

### Self-matrisin gerekçesi

Kendisiyle olan matris (`rep1`'in dizininde `cross_rmsd_rep1.xpm`) iki işi
birden görür: tek replika içindeki metastabil durumları gösterir, **ve** çıktı
kümesini her replika için tekdüze yapar. Tekdüzelik olmadan `ANALYSIS_OUTPUTS`
replikaya göre değişirdi ve `mdkit_mandatory_outputs` tabanlı idempotency
çalışmazdı (bkz. §6.1).

### Seyreltme

`-dt "${CROSS_RMSD_DT:-200}"`. Ortam değişkeniyle ezilir; `config.sh`'e
dokunulmaz, çünkü bu bir grafik okunabilirliği ayarıdır, projeye özel bir
varsayım değil.

### Güvenlik kontrolü

Çapraz `rms` **tek bir index dosyasını iki trajektoriye birden** uygular. Eşin
sistemi farklıysa (atom sayısı/sırası) matris sessizce çöp olur. Bu yüzden her
çift için, eşin `check_ref.pdb` atom sayısı kendi referansınınkiyle
karşılaştırılır; tutmazsa `>&2` mesajı ve `return 1`.

Eşin trajektorisi/referansı hiç yoksa aynı şekilde `return 1` — runner bunu
`run_log.csv`'ye `HATA` olarak yazar, kalan replikalar etkilenmez.

### Çağrı

```bash
"$GMX" rms -s "$REF_PDB" -f "$XTC" -f2 "$peer_xtc" -n "$NDX" \
    -m "$out_dir/cross_rmsd_${peer}.xpm" -o "$tmp_xvg" \
    -b "$B_PS" -dt "$dt" <<< $'RECEPTOR_BB\nLIGAND_BB'
```

`-o` açıkça verilmelidir. `gmx rms` `-o`'suz çağrıldığında `rmsd.xvg`'yi
**çalışma dizinine** yazar — yani `run_analysis.sh`'ın çağrıldığı yere, ki bu
kullanıcının kabuk dizinidir. Çapraz koşuda bu `.xvg`'nin bir değeri de yoktur.

Bu yüzden `-o`, `mktemp` ile `out_dir` **dışında** açılan geçici bir dosyaya
yönlendirilir ve `trap ... RETURN` ile her çıkış yolunda silinir. `out_dir`
içine yazılması yanlış olurdu: manifestoda ilan edilmeyen bir `.xvg`,
`collect_results.py`'nin "manifestoda olmayan .xvg atlandı" uyarısını her
toplamada tetiklerdi — üstelik eklenti yarıda hata verirse dosya kalıcı olurdu.

`-skip`/`-skip2` **kullanılmaz** (§2.2), `-tu` **kullanılmaz** (§2.2).

---

## 4. Katman 2 — `collect_results.py`

### `.xpm` ayrıştırıcısı

`parse_xvg`'nin yanına `parse_xpm(path) -> (meta, values, x_ps, y_ps)`:

1. `/* legend: "RMSD (nm)" */` → birim, parantez içinden.
2. `"W H NC CPP"` başlığı → boyutlar ve karakter/piksel.
3. `NC` adet renk satırı → `{karakter: float}` tablosu. Karakter alanı
   **`CPP` kadar sabit genişliktedir**, boşluk da geçerli bir karakterdir;
   `split()` ile ayrıştırılamaz, konumsal dilim gerekir.
4. `x-axis` / `y-axis` yorumları → aynı etiketli **tüm** satırlar birleştirilip
   float listesine çevrilir (§2.4-1).
5. Piksel satırları → `CPP` genişliğinde dilimlenip tabloyla eşlenir, sonra
   **satır sırası çevrilir** (§2.4-2).

Doğrulama: `len(x_ps) == W`, `len(y_ps) == H`, `values.shape == (H, W)`.
Tutmazsa dosya atlanır ve `stderr`'e adıyla bir uyarı yazılır.

### Keşif

`collect()` bugün yalnızca `adir.glob("*.xvg")` yapıyor
(`collect_results.py:115`) — bir `.xpm` bu yüzden "tanınmayan kind" uyarısını
bile tetiklemez, **hiç görülmez**. Glob `*.xvg` ve `*.xpm` olacak şekilde
genişletilir; manifesto eşlemesi (`--list`) aynen kullanılır.

### Çıktı biçimi

Matrisler uzun-format CSV'ye **girmez**: 451² × 315 ≈ 64 milyon satır olurdu.
Bunun yerine iki ürün:

**a) Matris başına `.npz`** — `results/matrices/<kompleks>_<rep_i>_<rep_j>.npz`:

| Anahtar | İçerik |
|---|---|
| `values` | `float32`, şekil `(len(y_ps), len(x_ps))`; satır = `rep_j` zamanı, kolon = `rep_i` zamanı |
| `x_ps` | `float64`, `rep_i` frame zamanları (ps) |
| `y_ps` | `float64`, `rep_j` frame zamanları (ps) |
| `unit` | gmx'in dediği birim (ör. `nm`), ham |
| `complex`, `replica_i`, `replica_j`, `analysis`, `output` | köken bilgisi |

**b) `results/matrix_summary.csv`** — tezde tabloya girecek sayılar:

```
complex,replica_i,replica_j,analysis,output,n_x,n_y,min,mean,max,unit
```

`rep_i == rep_j` olan self-matrislerde istatistikler **köşegen hariç** hesaplanır
(köşegen tanım gereği sıfırdır ve `min`'i anlamsız kılar). Bu, CSV'de ayrı bir
kolonla değil, `replica_i == replica_j` karşılaştırmasıyla anlaşılır.

Değerler `.xvg` tarafında olduğu gibi **ham** birimde yazılır; dönüşüm yalnızca
çizim katmanındadır.

---

## 5. Katman 3 — `plot_results.py`

Yeni mod `--matrix`, mod verilmediğinde koşan varsayılan kümeye katılır.

**Çıktı:** `plots/matrix/<kompleks>_cross_rmsd.png` — kompleks başına tek figür,
N×N ısı haritası ızgarası. Köşegende self-matrisler, köşegen dışında çapraz
çiftler. Satır/kolon başlıkları replika adları.

**Ortak renk skalası.** `vmin`/`vmax` o kompleksin **tüm** matrisleri üzerinden
hesaplanır ve tek bir colorbar çizilir. Panel başına ayrı skala, farklı çiftleri
görsel olarak karşılaştırılamaz kılardı — bu grafiğin tek amacı o karşılaştırma.

**Birim kuralı aynen geçerli.** `unit == "nm"` ise değerler ×10 ve colorbar
etiketi "Å"; tanınmayan birimde değer **çevrilmez**, ham birim gösterilir ve
`stderr`'e uyarı yazılır (mevcut `.xvg` davranışıyla birebir aynı).

**Eksenler** ns cinsinden gösterilir (ham veri ps), `origin="lower"` ile — `.npz`
zaten y ekseni artan sırada saklanıyor.

**Renk haritası:** algısal olarak düzgün ve tek yönlü bir skala (`viridis`).
Mevcut palet replika **çizgilerini** ayırmak için seçilmişti; ısı haritasının
ihtiyacı farklı, kategorik palet burada kullanılmaz.

**Yalıtım:** mevcut davranış korunur — tek bir kompleksin hatası adıyla
`stderr`'e yazılır, diğerleri çizilmeye devam eder.

**Boş durum:** `results/matrices/` yoksa veya boşsa mod tek satırlık bir bilgi
mesajıyla sessizce geçer, hata vermez (`cross_rmsd` hiç koşulmamış olabilir).

---

## 6. Reddedilen alternatifler

**6.1 Yalnızca çapraz çiftler (self-matris yok).** Kompleks başına 9 yerine 6
matris, %33 tasarruf. Reddedildi: çıktı kümesi replikaya göre değişirdi
(`rep1 -> rep2,rep3` ama `rep3 -> ` boş), `ANALYSIS_OUTPUTS` tekdüze
ilan edilemezdi ve son replika hiç zorunlu çıktı üretmediği için **her koşuda
yeniden koşardı**. §2.3'e göre tasarruf zaten 5 dakikalık.

**6.2 Matrisleri uzun-format CSV'ye koymak.** `timeseries_long.csv` ile simetrik
olurdu ama ≈64 milyon satır üretirdi. 1D veri için doğru olan biçim 2D için
değil.

**6.3 `plot_results.py`'nin `.xpm`'leri doğrudan okuması.** Toplama adımını
atlardı, ama çizim katmanının 105 dizini taramama sözleşmesini bozardı
(`README.md`, "Boru hattı").

**6.4 `-bin` ile ham değerler.** 80 seviyelik yuvarlamayı ortadan kaldırırdı;
ikinci bir dosya biçimi ve boyut bilgisi olmayan bir binary ayrıştırıcı
gerektirir. Isı haritası için kazanç yok (§2.5).

**6.5 Eşin trajektorisini önceden kırpmak.** §2.1 ile gereksiz olduğu ölçüldü.

---

## 7. Bilinen sınırlar

- Değerler 80 seviyeye yuvarlanmıştır (§2.5); `matrix_summary.csv` sayıları da
  bu çözünürlüktedir.
- `-r/--reps` ile kısıtlanmış bir koşuda üretilen matris sayısı `--list`
  manifestosundan azdır. Sessiz bir yanlışlık değil: eksik dosya toplanmaz,
  sonraki tam koşu eksikliği görüp yeniden üretir.
- `matrix` desteği yalnızca `cross_rmsd`'nin ürettiği **kare olmayan da
  olabilen** RMSD matrisleri için yazılır. DSSP gibi farklı semantikli 2B
  çıktılar (residue × zaman) bu ızgara çiziminden faydalanmaz; onlar için ayrı
  bir çizim biçimi gerekir.

---

## 8. Test planı

| Test | Tür | Doğruladığı |
|---|---|---|
| `parse_xpm` sabit fixture üzerinde | hızlı | Boyutlar, renk tablosu, **çok satırlı eksen birleştirme**, **satır sırası çevirme** |
| Self-matris fixture'ında köşegen | hızlı | Çevirme sonrası `values[i, i] == 0` — §2.4-2'nin regresyon testi |
| `CPP > 1` ve boşluk içeren karakter | hızlı | Konumsal dilimleme (`split()` kullanılmadığı) |
| `--list` çıktısında `cross_rmsd` | hızlı | `$REPS`'ten türetilen `ANALYSIS_OUTPUTS` manifestoya doğru giriyor |
| `tests/test_portability.py` | hızlı | Yeni dosyalarda mutlak yol/proje adı yok |
| Gerçek veride tek çift `cross_rmsd` | yavaş | gmx çağrısı, `-dt`, atom sayısı kontrolü |
| Uçtan uca: koş → topla → çiz | yavaş | Üç katmanın sözleşmesi |

Fixture'lar `tests/` altına küçük (≤ 9×9) elle yazılmış `.xpm` dosyaları olarak
konur; gerçek veriye bağlı testler mevcut `MDKIT_TEST_*` değişkenleriyle
yönlendirilebilir kalır.
