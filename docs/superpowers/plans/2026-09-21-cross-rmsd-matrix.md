# Çapraz-RMSD Matrisleri ve `matrix` Katman Desteği — Uygulama Planı

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replikalar arası 2D RMSD matrisi üreten bir `cross_rmsd` eklentisi yazmak ve `ANALYSIS_KIND="matrix"`'i toplama ve çizim katmanlarında uçtan uca tüketilebilir hâle getirmek.

**Architecture:** Üç katmanlı mevcut sözleşme korunur. Eklenti (`analysis/cross_rmsd.sh`) `gmx rms -f -f2 -m` ile `.xpm` matrisleri üretir; `collect_results.py` bunları ayrıştırıp `results/matrices/*.npz` + `results/matrix_summary.csv` hâline getirir; `plot_results.py` yalnızca `results/` altını okuyup N×N ısı haritası ızgarası çizer. `run_analysis.sh` ve `analysis/lib.sh` **değişmez**.

**Tech Stack:** bash 5+, GROMACS 2025.4, Python 3.11+, numpy, pandas, matplotlib, pytest.

**Spec:** [`docs/superpowers/specs/2026-09-21-cross-rmsd-matrix-design.md`](../specs/2026-09-21-cross-rmsd-matrix-design.md)

## Global Constraints

- **`-skip`/`-skip2` asla kullanılmaz.** `.xpm` eksen zaman değerlerini bozar (ilk zaman doğru, gerisi `0`). Seyreltme yalnızca `-dt` (ps) ile. — spec §2.2
- **`-tu` asla kullanılmaz.** `-b`/`-e` değerlerini de çevirir. Araç içinde zaman daima **ps**, mesafe daima **nm**. — `README.md`
- **`-b` her iki trajektoriye de uygulanır.** Eş trajektoriyi önceden kırpmak **yasak** — gereksizdir ve ölçümle doğrulanmıştır. — spec §2.1
- **`.xpm` piksel satırları y ekseninin tersi sırada yazılır.** Ayrıştırıcı satır sırasını çevirmelidir. — spec §2.4-2
- **`.xpm` karakter alanı `CPP` genişliğinde sabittir ve boşluk geçerli bir karakterdir.** `split()` ile ayrıştırılamaz, konumsal dilim şarttır. — spec §2.4-3
- **Değerler ham birimde saklanır** (`nm`, `ps`). Dönüşüm yalnızca `plot_results.py` içinde ve yalnızca `unit == "nm"` için. Tanınmayan birim çevrilmez, `stderr`'e uyarı yazılır.
- **Mutlak yol ve proje adı yasak.** `tests/test_portability.py` `*.sh`/`*.py` dosyalarında `/mnt/`, `/usr/local/`, `/home/`, `pandora`, `TUSEB`, `A0201` arar (`config.sh` ve `tests/` hariç).
- **`set -e` kullanılmaz**; her gmx hatası `>&2` mesajı + `return 1`.

### Spec'ten tek sapma (kasıtlı)

Spec §4 `.npz` adını `<kompleks>_<rep_i>_<rep_j>.npz` diyor. Plan bunun yerine
**`<kompleks>_<rep_i>_<çıktı_kökü>.npz`** kullanır (ör.
`last10_rep1_cross_rmsd_rep2.npz`). Gerekçe: `<rep_j>` kalıbı `matrix` desteğini
replikalar-arası analizlere kilitler; çıktı köküyle adlandırmak spec §7'nin
"`matrix` desteği ileride başka 2B çıktılar için de gerekecek" notunu karşılar ve
aynı komplekste iki farklı matris analizi olduğunda çakışmayı önler.

---

## Dosya Yapısı

| Dosya | Sorumluluk | Durum |
|---|---|---|
| `analysis/cross_rmsd.sh` | Eklenti sözleşmesi + `gmx rms` çağrıları | **Oluştur** (mevcut taslağın üstüne yazılır) |
| `collect_results.py` | `parse_xpm` + matris keşfi + `.npz`/`matrix_summary.csv` yazımı | Değiştir |
| `plot_results.py` | `--matrix` modu, N×N ısı haritası ızgarası | Değiştir |
| `tests/test_xpm.py` | `parse_xpm` birim testleri (fixture'lar dosya içinde) | **Oluştur** |
| `tests/test_cross_rmsd.py` | Eklenti sözleşmesi, atom sayısı koruması, gerçek gmx koşusu | **Oluştur** |
| `tests/test_collect.py` | Matris toplama entegrasyonu | Değiştir |
| `tests/test_plot.py` | `--matrix` çizimi | Değiştir |
| `README.md` | "matrix henüz desteklenmiyor" bölümünün yerine matris boru hattı | Değiştir |

`analysis/lib.sh` ve `run_analysis.sh` **hiç değişmez** — eklentinin kardeş
replikalara erişimi zaten mevcut altyapıdan geliyor (`mdkit_run_isolated`
`analysis_run`'ı ana kabuğun alt kabuğunda koşturur, `analysis/lib.sh:222`).

---

## Task 1: `.xpm` ayrıştırıcısı

**Files:**
- Modify: `collect_results.py` (yeni `parse_xpm`, `parse_xvg`'nin hemen altına)
- Test: `tests/test_xpm.py` (yeni)

**Interfaces:**
- Consumes: yok (ilk task)
- Produces: `collect_results.parse_xpm(path) -> (meta: dict, values: np.ndarray[float32, (H, W)], x_ps: np.ndarray[float64, (W,)], y_ps: np.ndarray[float64, (H,)])`. `meta["unit"]` gmx'in legend'ından çıkarılan ham birim (ör. `"nm"`), `meta["legend"]` ham legend metni. Hatalı dosyada `ValueError` fırlatır.

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_xpm.py` oluştur:

```python
import sys
import textwrap

import numpy as np
import pytest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import collect_results  # noqa: E402

# 3x3 self-matris. Sifir kosegen DOSYADA ters kosegen uzerinde durur
# (gmx piksel satirlarini y ekseninin TERSI sirada yazar, spec 2.4-2).
# x-ekseni bilerek IKI yoruma bolundu: uzun eksenlerde gmx boyle yazar.
SELF_XPM = textwrap.dedent("""\
    /* XPM */
    /* title:   "LIGAND_BB RMSD matrix" */
    /* legend:  "RMSD (nm)" */
    /* x-label: "Time (ps)" */
    /* y-label: "Time (ps)" */
    /* type:    "Continuous" */
    static char *gromacs_xpm[] = {
    "3 3   3 1",
    "A  c #FFFFFF " /* "0" */,
    "B  c #808080 " /* "0.5" */,
    "C  c #000000 " /* "1" */,
    /* x-axis:  10000 10200 */
    /* x-axis:  10400 */
    /* y-axis:  10000 10200 10400 */
    "CBA",
    "BAB",
    "ABC"
    """)

# CPP=2: ilk anahtar 'A' + BOSLUK. split() ile ayristirilirsa anahtar 'A'
# olur ve piksel eslemesi KeyError verir.
CPP2_XPM = textwrap.dedent("""\
    /* XPM */
    /* legend:  "RMSD (nm)" */
    static char *gromacs_xpm[] = {
    "2 2   2 2",
    "A  c #FFFFFF " /* "0" */,
    "BB c #000000 " /* "1" */,
    /* x-axis:  0 100 */
    /* y-axis:  0 100 */
    "A BB",
    "BBA "
    """)

# Genislik != yukseklik: -f ve -f2 farkli uzunlukta olabilir.
RECT_XPM = textwrap.dedent("""\
    /* XPM */
    /* legend:  "RMSD (nm)" */
    static char *gromacs_xpm[] = {
    "2 3   2 1",
    "A  c #FFFFFF " /* "0" */,
    "B  c #000000 " /* "1" */,
    /* x-axis:  0 100 */
    /* y-axis:  0 100 200 */
    "AA",
    "AB",
    "BB"
    """)


def _write(tmp_path, text, name="m.xpm"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_boyutlar_ve_eksenler(tmp_path):
    meta, values, x, y = collect_results.parse_xpm(_write(tmp_path, SELF_XPM))
    assert values.shape == (3, 3)
    assert np.allclose(x, [10000.0, 10200.0, 10400.0])
    assert np.allclose(y, [10000.0, 10200.0, 10400.0])


def test_cok_satirli_eksen_birlestirilir(tmp_path):
    """x-ekseni iki /* x-axis: */ satirina bolunmus; yalnizca ilkini okuyan
    bir ayristirici ekseni 2 elemanda keser ve sekil dogrulamasi patlar."""
    _meta, values, x, _y = collect_results.parse_xpm(_write(tmp_path, SELF_XPM))
    assert len(x) == 3
    assert values.shape[1] == 3


def test_satir_sirasi_cevrilir_kosegen_sifir(tmp_path):
    """spec 2.4-2 regresyonu: cevirme yapilmazsa matris yatayda aynalanir ve
    bu, self-matris disinda GOZLE FARK EDILMEZ."""
    _meta, values, _x, _y = collect_results.parse_xpm(_write(tmp_path, SELF_XPM))
    assert np.allclose(np.diag(values), 0.0)
    assert values[0, 2] == pytest.approx(1.0)
    assert values[2, 0] == pytest.approx(1.0)


def test_cpp2_ve_bosluklu_karakter(tmp_path):
    """spec 2.4-3 regresyonu: karakter alani konumsal dilimle okunmali."""
    _meta, values, _x, _y = collect_results.parse_xpm(_write(tmp_path, CPP2_XPM))
    assert np.allclose(values, [[1.0, 0.0], [0.0, 1.0]])


def test_birim_legendden_okunur(tmp_path):
    meta, _v, _x, _y = collect_results.parse_xpm(_write(tmp_path, SELF_XPM))
    assert meta["unit"] == "nm"


def test_kare_olmayan_matris(tmp_path):
    _meta, values, x, y = collect_results.parse_xpm(_write(tmp_path, RECT_XPM))
    assert values.shape == (3, 2) == (len(y), len(x))
    assert np.allclose(values[0], [1.0, 1.0])


def test_bozuk_dosya_valueerror(tmp_path):
    """Basliktaki yukseklik piksel satiri sayisiyla tutmuyor."""
    bad = SELF_XPM.replace('"3 3   3 1"', '"3 5   3 1"')
    with pytest.raises(ValueError):
        collect_results.parse_xpm(_write(tmp_path, bad))
```

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_xpm.py -v`
Expected: FAIL — `AttributeError: module 'collect_results' has no attribute 'parse_xpm'`

- [ ] **Step 3: `parse_xpm`'i uygula**

`collect_results.py`'nin başına `import numpy as np` ekle (mevcut importların arasına, alfabetik olarak `re`'den önce üçüncü-parti blok olarak). Regex sabitlerini `UNIT_IN_LABEL`'ın altına ekle:

```python
XPM_LEGEND = re.compile(r'/\*\s*legend:\s*"(.*)"\s*\*/')
XPM_AXIS = re.compile(r'/\*\s*([xy])-axis:\s*(.*?)\s*\*/')
XPM_HEADER = re.compile(r'^"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"')
XPM_VALUE = re.compile(r'/\*\s*"(.*?)"\s*\*/')
```

`parse_xvg`'nin hemen altına:

```python
def parse_xpm(path):
    """gmx .xpm matrisini (meta, values, x_ps, y_ps) olarak dondurur.

    UC TUZAK, ucu de GROMACS 2025.4 uzerinde olculerek saptandi (spec 2.4):

    1. Uzun eksenler BIRDEN FAZLA `/* x-axis: */` yorumuna bolunur. Yalnizca
       ilkini okuyan bir ayristirici ekseni sessizce keser.
    2. Piksel satirlari y ekseninin TERSI sirada yazilir. Cevrilmezse matris
       yatayda aynalanir -- self-matris disinda gozle fark edilmez.
    3. Karakter alani CPP genisliginde SABITTIR ve bosluk da gecerli bir
       karakterdir. split() ile ayristirmak CPP>1'de anahtari bozar.
    """
    colors = {}
    axes = {"x": [], "y": []}
    pixel_rows = []
    meta = {}
    width = height = cpp = None
    colors_left = 0

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("/*"):
            m = XPM_LEGEND.search(line)
            if m:
                meta["legend"] = m.group(1)
                continue
            m = XPM_AXIS.search(line)
            if m:
                # TUZAK 1: extend, assign DEGIL.
                axes[m.group(1)].extend(float(v) for v in m.group(2).split())
            continue
        if not line.startswith('"'):
            continue
        if width is None:
            m = XPM_HEADER.match(line)
            if m:
                width, height, ncolors, cpp = (int(g) for g in m.groups())
                colors_left = ncolors
            continue
        if colors_left > 0:
            # TUZAK 3: konumsal dilim.
            key = line[1:1 + cpp]
            vm = XPM_VALUE.search(line)
            if vm is None:
                raise ValueError(f"{path.name}: renk satiri cozulemedi: {line}")
            colors[key] = float(vm.group(1))
            colors_left -= 1
            continue
        pixel_rows.append(line[1:1 + width * cpp])

    if width is None:
        raise ValueError(f"{path.name}: .xpm basligi bulunamadi")
    if len(pixel_rows) != height:
        raise ValueError(
            f"{path.name}: basligi {height} satir diyor, {len(pixel_rows)} bulundu"
        )

    values = np.empty((height, width), dtype=np.float32)
    for r, row in enumerate(pixel_rows):
        for c in range(width):
            key = row[c * cpp:(c + 1) * cpp]
            try:
                values[r, c] = colors[key]
            except KeyError:
                raise ValueError(
                    f"{path.name}: renk tablosunda olmayan karakter {key!r}"
                ) from None
    # TUZAK 2.
    values = values[::-1].copy()

    x_ps = np.asarray(axes["x"], dtype=np.float64)
    y_ps = np.asarray(axes["y"], dtype=np.float64)
    if values.shape != (len(y_ps), len(x_ps)):
        raise ValueError(
            f"{path.name}: matris {values.shape}, eksenler "
            f"({len(y_ps)}, {len(x_ps)}) -- uyusmuyor"
        )

    m = UNIT_IN_LABEL.search(meta.get("legend", ""))
    meta["unit"] = m.group(1) if m else ""
    return meta, values, x_ps, y_ps
```

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_xpm.py -v`
Expected: 7 PASS

- [ ] **Step 5: Commit**

```bash
git add collect_results.py tests/test_xpm.py
git commit -m "feat(collect): gmx .xpm matris ayristiricisi

Uc tuzagi da testle sabitler: cok satirli eksen yorumlari, y ekseninin
tersi sirada yazilan piksel satirlari, CPP genisliginde sabit karakter
alani. Ikincisi cevrilmezse matris sessizce aynalanir."
```

---

## Task 2: `cross_rmsd` eklenti sözleşmesi

**Files:**
- Create: `analysis/cross_rmsd.sh` (mevcut taslağın üstüne yazılır)
- Test: `tests/test_cross_rmsd.py` (yeni)

**Interfaces:**
- Consumes: yok
- Produces: `analysis/cross_rmsd.sh` — `ANALYSIS_NAME="cross_rmsd"`, `ANALYSIS_KIND="matrix"`, `ANALYSIS_NEEDS_INDEX=1`, `ANALYSIS_DEFAULT_BEGIN=10000`, `ANALYSIS_OUTPUTS` = `$REPS`'teki her `<rep>` için `cross_rmsd_<rep>.xpm`. `analysis_run` bu task'ta yalnızca bir iskelet (`return 0`); gerçek gmx çağrısı Task 3'te.

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_cross_rmsd.py` oluştur:

```python
import subprocess

from conftest import run_bash


def _list(mdkit, config):
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(config), "--list"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        if parts and parts[0] == "cross_rmsd":
            return parts
    raise AssertionError(f"--list ciktisinda cross_rmsd yok:\n{r.stdout}")


def test_manifestoda_matrix_kind(mdkit, fake_config):
    parts = _list(mdkit, fake_config)
    assert parts[1] == "matrix"


def test_ciktilar_repsten_turetilir(mdkit, fake_config):
    """fake_config REPS=(rep1 rep2 rep3)."""
    parts = _list(mdkit, fake_config)
    assert parts[3] == ("cross_rmsd_rep1.xpm,cross_rmsd_rep2.xpm,"
                        "cross_rmsd_rep3.xpm")


def test_farkli_reps_farkli_manifesto(mdkit, fake_config, tmp_path):
    """Cikti kumesi config'e bagli; dosyada sabit kodlu degil."""
    alt = tmp_path / "alt_config.sh"
    alt.write_text(
        fake_config.read_text().replace("REPS=(rep1 rep2 rep3)",
                                        "REPS=(repA repB)")
    )
    parts = _list(mdkit, alt)
    assert parts[3] == "cross_rmsd_repA.xpm,cross_rmsd_repB.xpm"


def test_opsiyonel_cikti_yok(mdkit, fake_config):
    """Her replika her matrisi URETIR; kosullu cikti yok, yoksa idempotency
    kontrolu eksik dosya arar ve analiz her kosuda yeniden kosar."""
    parts = _list(mdkit, fake_config)
    assert parts[4].strip() == ""


def test_sozlesme_dogrulamasindan_gecer(mdkit, lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'mdkit_clear_plugin && source "{mdkit}/analysis/cross_rmsd.sh" && '
        f'mdkit_validate_plugin "{mdkit}/analysis/cross_rmsd.sh"'
    )
    assert r.returncode == 0, r.stderr


def test_yardimci_degisken_sizmiyor(mdkit, lib, fake_config):
    """Eklenti ANA KABUGA source edilir; dongu degiskeni orada kalmamali."""
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        f'source "{mdkit}/analysis/cross_rmsd.sh" && '
        f'echo "SIZDI=${{_cross_rmsd_rep:-yok}}"'
    )
    assert "SIZDI=yok" in r.stdout, r.stdout
```

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_cross_rmsd.py -v`
Expected: FAIL — mevcut taslak `ANALYSIS_KIND="timeseries"` ve tek `.xvg` çıktısı ilan ediyor

- [ ] **Step 3: Eklenti sözleşmesini yaz**

`analysis/cross_rmsd.sh` dosyasının TAMAMINI aşağıdakiyle değiştir:

```bash
#!/usr/bin/env bash
# Capraz-RMSD eklentisi. run_analysis.sh tarafindan source edilir.
#
# rep_i baglaminda kosar ve $REPS'teki HER replikaya (kendisi dahil) karsi bir
# frame x frame RMSD matrisi uretir: RECEPTOR_BB uzerine fit, LIGAND_BB uzerinde
# RMSD.
#
# Cevapladigi soru: uc replika ayni konformasyonel bolgeleri mi geziyor?
# Mevcut `rmsd` analizi bunu cevaplayamaz -- her replikayi KENDI referansina
# gore olctugu icin iki replika ayni RMSD degerine FARKLI konformasyonlarda
# ulasabilir.
#
# Kardes replikalara erisim: mdkit_run_isolated analysis_run'i ana kabugun ALT
# KABUGUNDA kosturur (analysis/lib.sh:222), yani $REPS/$TRAJ_NAME/$REF_NAME
# gorunurdur ve kardesler dirname "$rep_dir" ile bulunur. run_analysis.sh
# degistirilmedi.
ANALYSIS_NAME="cross_rmsd"
ANALYSIS_DESC="Replikalar arasi 2D RMSD matrisleri (peptid backbone, oluk-uzeri fit)"
ANALYSIS_KIND="matrix"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000

# Ciktilar $REPS'ten turetilir. Self-matris (rep1'in dizininde
# cross_rmsd_rep1.xpm) iki isi birden gorur: tek replika icindeki metastabil
# durumlari gosterir VE cikti kumesini her replika icin TEKDUZE yapar.
# Tekduzelik olmadan ANALYSIS_OUTPUTS replikaya gore degisirdi (rep1 ->
# rep2,rep3 ama son replika -> bos) ve mdkit_mandatory_outputs tabanli
# idempotency calismazdi: hicbir zorunlu cikti uretmeyen replika her kosuda
# yeniden kosardi.
ANALYSIS_OUTPUTS=()
for _cross_rmsd_rep in ${REPS[@]+"${REPS[@]}"}; do
    ANALYSIS_OUTPUTS+=("cross_rmsd_${_cross_rmsd_rep}.xpm")
done
unset _cross_rmsd_rep

# Ciktilarin hepsi KOSULSUZ uretilir.
ANALYSIS_OPTIONAL_OUTPUTS=()

analysis_run() {
    local rep_dir="$1" out_dir="$2"
    return 0
}
```

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_cross_rmsd.py -v`
Expected: 6 PASS

- [ ] **Step 5: Diğer testlerin bozulmadığını doğrula**

Run: `cd tests && python -m pytest -q`
Expected: hiçbir yeni FAIL (`matrix` kind'lı bir eklentinin varlığı `collect_results.py`'de uyarı üretebilir; bu Task 4'te ele alınır, FAIL üretmemeli)

- [ ] **Step 6: Commit**

```bash
git add analysis/cross_rmsd.sh tests/test_cross_rmsd.py
git commit -m "feat(cross_rmsd): eklenti sozlesmesi, ciktilar REPS'ten turetilir

Self-matris dahil her replikaya karsi bir matris ilan edilir. Tekduze
cikti kumesi olmadan idempotency calismaz: hicbir zorunlu cikti
uretmeyen replika her kosuda yeniden kosardi."
```

---

## Task 3: `cross_rmsd` gmx çağrıları

**Files:**
- Modify: `analysis/cross_rmsd.sh` (`analysis_run` gövdesi + iki yardımcı)
- Test: `tests/test_cross_rmsd.py` (yeni test'ler eklenir)

**Interfaces:**
- Consumes: Task 2'nin sözleşme değişkenleri
- Produces: `<rep_dir>/analysis/cross_rmsd_<peer>.xpm` — `$REPS`'teki her `<peer>` için. Task 4 bu dosyaları tüketir.

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_cross_rmsd.py` sonuna ekle:

```python
CONTRACT = """
source "{lib}"
mdkit_load_config "{config}"
GMX="{gmx}"
REF_PDB="{rep}/check_ref.pdb"
XTC="{rep}/traj_compact_center_dry.xtc"
NDX="{rep}/analysis/index.ndx"
B_PS=0
FORCE=0
mdkit_clear_plugin
source "{mdkit}/analysis/cross_rmsd.sh"
analysis_run "{rep}" "{out}"
"""


def test_uyusmayan_es_sistemi_reddedilir(mdkit, lib, fake_config,
                                         fake_dataset, tmp_path):
    """Capraz rms TEK index'i IKI trajektoriye birden uygular; esin sistemi
    farkliysa matris sessizce cop olur. gmx yolu bilerek gecersiz: koruma
    gmx cagrilmadan ONCE donmeli."""
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    (cx / "rep1" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    (cx / "rep2" / "check_ref.pdb").write_text("ATOM      1  N\n" * 4)
    (cx / "rep3" / "check_ref.pdb").write_text("ATOM      1  N\n" * 10)
    out = tmp_path / "out"
    out.mkdir()
    r = run_bash(CONTRACT.format(
        lib=lib, config=fake_config, gmx="/nonexistent/gmx",
        rep=cx / "rep1", mdkit=mdkit, out=out,
    ))
    assert r.returncode != 0
    assert "uyusmuyor" in (r.stdout + r.stderr).lower()


def test_eksik_es_trajektorisi_reddedilir(mdkit, lib, fake_config,
                                          fake_dataset, tmp_path):
    cx = fake_dataset / "top1_BBB_A0201_pandora"
    (cx / "rep2" / "traj_compact_center_dry.xtc").unlink()
    out = tmp_path / "out2"
    out.mkdir()
    r = run_bash(CONTRACT.format(
        lib=lib, config=fake_config, gmx="/nonexistent/gmx",
        rep=cx / "rep1", mdkit=mdkit, out=out,
    ))
    assert r.returncode != 0
    assert "trajektori" in (r.stdout + r.stderr).lower()


def test_skip_ve_tu_kullanilmiyor(mdkit):
    """spec 2.2: -skip .xpm eksen zamanlarini bozar, -tu -b/-e'yi cevirir."""
    src = (mdkit / "analysis" / "cross_rmsd.sh").read_text()
    assert " -skip" not in src
    assert " -tu" not in src
    assert " -dt " in src


@needs_gmx
def test_gercek_capraz_matris(mdkit, real_config):
    """real_config tek replikali (REPS=(rep1)); self-matris uretilmeli.

    small_rep 0-200 ps'lik 21 frame'dir, bu yuzden -b 0 ve CROSS_RMSD_DT=50
    verilir: eklentinin varsayilani (ANALYSIS_DEFAULT_BEGIN=10000) bu kisa
    trajektoride bos cikti uretirdi."""
    root = subprocess.run(
        ["bash", "-c",
         f'source "{mdkit}/analysis/lib.sh" && '
         f'mdkit_load_config "{real_config}" && printf "%s" "$DATA_ROOT"'],
        capture_output=True, text=True, check=True,
    ).stdout
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-a", "cross_rmsd", "-b", "0", "--all", root],
        capture_output=True, text=True,
        env={**os.environ, "CROSS_RMSD_DT": "50"},
    )
    assert r.returncode == 0, r.stderr
    xpm = next(Path(root).glob("*/rep1/analysis/cross_rmsd_rep1.xpm"), None)
    assert xpm is not None, r.stdout

    sys.path.insert(0, str(mdkit))
    import collect_results
    meta, values, x, y = collect_results.parse_xpm(xpm)
    assert values.shape[0] == values.shape[1] == len(x) == len(y)
    assert meta["unit"] == "nm"
    # Self-matriste kosegen TANIM GEREGI sifir. Ayristiricinin satir cevirmesi
    # ile gmx'in GERCEK yazma sirasinin uyustugunun kaniti -- fixture'lar bunu
    # gosteremez, cunku fixture'i da biz yaziyoruz.
    assert np.allclose(np.diag(values), 0.0, atol=1e-3)
```

`tests/test_cross_rmsd.py` dosyasının başındaki import bloğunu şu hâle getir:

```python
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from conftest import needs_gmx, run_bash
```

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_cross_rmsd.py -v`
Expected: `test_uyusmayan_es_sistemi_reddedilir`, `test_eksik_es_trajektorisi_reddedilir`, `test_skip_ve_tu_kullanilmiyor` FAIL (iskelet `analysis_run` her zaman 0 dönüyor, dosyada `-dt` yok)

- [ ] **Step 3: `analysis_run`'ı uygula**

`analysis/cross_rmsd.sh` içindeki iskelet `analysis_run`'ı şununla değiştir:

```bash
_cross_rmsd_atom_count() {
    # $1 = PDB. Mutlak dogruluk degil, iki replika ARASINDA esitlik onemli.
    local n
    n="$(grep -c -E '^(ATOM|HETATM)' "$1" 2>/dev/null)" || n=0
    printf '%s' "${n:-0}"
}

_cross_rmsd_pairs() {
    # $1 = rep_dir, $2 = out_dir, $3 = gecici .xvg yolu
    local rep_dir="$1" out_dir="$2" tmp_xvg="$3"
    local cx_dir peer peer_dir peer_xtc peer_ref out own_atoms peer_atoms dt
    cx_dir="$(dirname "$rep_dir")"
    dt="${CROSS_RMSD_DT:-200}"
    own_atoms="$(_cross_rmsd_atom_count "$REF_PDB")"

    for peer in ${REPS[@]+"${REPS[@]}"}; do
        peer_dir="$cx_dir/$peer"
        peer_xtc="$peer_dir/$TRAJ_NAME"
        peer_ref="$peer_dir/$REF_NAME"
        out="$out_dir/cross_rmsd_${peer}.xpm"

        if [[ ! -s "$peer_xtc" ]]; then
            echo "es replikanin trajektorisi yok: $peer_xtc" >&2
            return 1
        fi

        # Capraz rms TEK bir index dosyasini IKI trajektoriye birden uygular.
        # Esin sistemi farkliysa (atom sayisi/sirasi) matris SESSIZCE cop olur.
        peer_atoms="$(_cross_rmsd_atom_count "$peer_ref")"
        if [[ "$peer_atoms" != "$own_atoms" ]]; then
            echo "es replikanin sistemi uyusmuyor ($peer: $peer_atoms atom, kendi: $own_atoms)" >&2
            return 1
        fi

        # -skip KULLANILMAZ: .xpm eksen zamanlarini bozar -- ilk zaman dogru,
        #   gerisi 0 yazilir (GROMACS 2025.4'te olculdu). Seyreltme -dt ile.
        # -tu KULLANILMAZ: -b/-e degerlerini de cevirir.
        # -b HER IKI trajektoriye de uygulanir (olculdu), es icin ayrica
        #   kirpma gerekmez.
        # -o acikca verilir: gmx rms -o'suz cagrildiginda rmsd.xvg'yi CALISMA
        #   DIZININE yazar. out_dir'e yazmak da olmaz -- manifestoda ilan
        #   edilmeyen bir .xvg her toplamada uyari uretirdi.
        if [[ "$peer_xtc" == "$XTC" ]]; then
            # Self-matris: -f2 verilmez. Olculen ve kosegeni sifir dogrulanan
            # bicim budur.
            "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
                -m "$out" -o "$tmp_xvg" -b "$B_PS" -dt "$dt" \
                <<< $'RECEPTOR_BB\nLIGAND_BB' \
                || { echo "gmx rms basarisiz: cross_rmsd_${peer}" >&2; return 1; }
        else
            "$GMX" rms -s "$REF_PDB" -f "$XTC" -f2 "$peer_xtc" -n "$NDX" \
                -m "$out" -o "$tmp_xvg" -b "$B_PS" -dt "$dt" \
                <<< $'RECEPTOR_BB\nLIGAND_BB' \
                || { echo "gmx rms basarisiz: cross_rmsd_${peer}" >&2; return 1; }
        fi

        [[ -s "$out" ]] || {
            echo "matris uretilmedi: cross_rmsd_${peer}" >&2
            return 1
        }
    done
    return 0
}

analysis_run() {
    local rep_dir="$1" out_dir="$2" tmp_dir rc
    tmp_dir="$(mktemp -d)" || { echo "gecici dizin acilamadi" >&2; return 1; }
    _cross_rmsd_pairs "$rep_dir" "$out_dir" "$tmp_dir/rms.xvg"
    rc=$?
    rm -rf "$tmp_dir"
    return "$rc"
}
```

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_cross_rmsd.py -v`
Expected: hızlı test'ler PASS; `test_gercek_capraz_matris` gmx+veri varsa PASS, yoksa SKIP

- [ ] **Step 5: Portability test'inin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_portability.py -v`
Expected: 3 PASS

- [ ] **Step 6: Commit**

```bash
git add analysis/cross_rmsd.sh tests/test_cross_rmsd.py
git commit -m "feat(cross_rmsd): capraz ve self RMSD matrislerini uret

Seyreltme -dt ile: -skip, .xpm eksen zamanlarini bozuyor (ilk zaman
dogru, gerisi 0). Es trajektori icin on-kirpma yok -- -b her iki
trajektoriye de uygulaniyor (olculdu). Tek index iki trajektoriye
birden uygulandigi icin atom sayilari once karsilastiriliyor."
```

---

## Task 4: Matris toplama ve `.npz` yazımı

**Files:**
- Modify: `collect_results.py` (`collect` imzası + `.xpm` dalı, `main`)
- Test: `tests/test_collect.py`

**Interfaces:**
- Consumes: `collect_results.parse_xpm` (Task 1), `cross_rmsd_<peer>.xpm` dosyaları (Task 3)
- Produces:
  - `collect_results.peer_replica(output_name, reps) -> str | None`
  - `collect_results.collect(data_root, complex_glob, reps, manifest, matrix_dir=None) -> (timeseries, profile, matrices)` — `matrices` bir özet dict listesi (Task 5'te CSV'ye yazılır), alanları: `complex, replica_i, replica_j, analysis, output, n_x, n_y, unit` ve `values` **içermez**.
  - `<matrix_dir>/<kompleks>_<rep_i>_<çıktı_kökü>.npz` — anahtarlar: `values`, `x_ps`, `y_ps`, `unit`, `complex`, `replica_i`, `replica_j`, `analysis`, `output`.

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_collect.py` sonuna ekle (dosyanın başındaki importlara `import numpy as np` ekle):

```python
MATRIX_XPM = textwrap.dedent("""\
    /* XPM */
    /* legend:  "RMSD (nm)" */
    static char *gromacs_xpm[] = {
    "3 3   3 1",
    "A  c #FFFFFF " /* "0" */,
    "B  c #808080 " /* "0.5" */,
    "C  c #000000 " /* "1" */,
    /* x-axis:  0 100 200 */
    /* y-axis:  0 100 200 */
    "CBA",
    "BAB",
    "ABC"
    """)


def test_peer_replica_bilinen_adlarla_eslesir():
    reps = ["rep1", "rep2", "rep3"]
    assert collect_results.peer_replica("cross_rmsd_rep2.xpm", reps) == "rep2"
    assert collect_results.peer_replica("dssp.xpm", reps) is None


def test_matris_npz_olarak_yazilir(fake_dataset, fake_config, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep2.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep2.xpm": ("cross_rmsd", "matrix")}
    mdir = tmp_path / "matrices"
    ts, pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=mdir,
    )
    assert ts == [] and pr == []
    # 2 kompleks x 3 replika
    assert len(mx) == 6
    f = mdir / "last1_rep1_cross_rmsd_rep2.npz"
    assert f.exists(), sorted(p.name for p in mdir.glob("*"))
    d = np.load(f, allow_pickle=False)
    assert d["values"].shape == (3, 3)
    assert np.allclose(np.diag(d["values"]), 0.0)
    assert np.allclose(d["x_ps"], [0.0, 100.0, 200.0])
    assert str(d["unit"]) == "nm"
    assert str(d["replica_i"]) == "rep1"
    assert str(d["replica_j"]) == "rep2"


def test_matris_kind_artik_uyari_uretmiyor(fake_dataset, capsys, tmp_path):
    """Once .xpm hic gorulmuyordu (*.xvg glob'u); simdi toplaniyor."""
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m",
    )
    err = capsys.readouterr().err
    assert "taninmayan ANALYSIS_KIND" not in err


def test_manifestoda_olmayan_xpm_uyarir(fake_dataset, capsys, tmp_path):
    write_outputs(fake_dataset, {"baskabir.xpm": MATRIX_XPM})
    collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], {},
        matrix_dir=tmp_path / "m",
    )
    err = capsys.readouterr().err
    assert "baskabir.xpm" in err


def test_bozuk_xpm_kosuyu_durdurmaz(fake_dataset, capsys, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": "bozuk icerik\n"})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m",
    )
    assert mx == []
    assert "cross_rmsd_rep1.xpm" in capsys.readouterr().err


def test_matrix_dir_yoksa_matris_atlanir(fake_dataset, tmp_path):
    """collect(), matrix_dir verilmediginde .xpm'leri hic islemez."""
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest)
    assert mx == []
```

`tests/test_collect.py` içinde `collect_results.collect(...)` çağıran **mevcut** test'ler artık üç değer döndüğü için güncellenmeli: her `ts, pr = collect_results.collect(...)` satırını `ts, pr, _mx = collect_results.collect(...)` yap.

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_collect.py -v`
Expected: yeni test'ler FAIL (`peer_replica` yok, `collect` iki değer dönüyor)

- [ ] **Step 3: `collect`'i genişlet**

`collect_results.py`'de `TS_FIELDS`/`PR_FIELDS`'ın altına ekle:

```python
MX_FIELDS = ["complex", "replica_i", "replica_j", "analysis", "output",
             "n_x", "n_y", "min", "mean", "max", "unit"]
```

`collect`'in üstüne ekle:

```python
def peer_replica(output_name, reps):
    """'cross_rmsd_rep2.xpm' -> 'rep2'; eslesme yoksa None.

    Dosya adi kalibini burada yeniden tanimlamiyoruz: config'den gelen
    BILINEN replika adlariyla eslestiriyoruz. Boylece `matrix` destegi
    replikalar-arasi analizlere kilitlenmez -- eslesmeyen bir matris
    (or. residue x zaman) da toplanir, yalnizca replica_j'si bos olur."""
    stem = Path(output_name).stem
    for rep in reps:
        if stem.endswith("_" + rep):
            return rep
    return None
```

`collect`'i şu şekilde değiştir — imza, dönüş değeri ve `.xpm` dalı:

```python
def collect(data_root, complex_glob, reps, manifest, matrix_dir=None):
    timeseries, profile, matrices = [], [], []
    warned = set()
    warned_unlisted = set()
    for cx in sorted(data_root.glob(complex_glob)):
        if not cx.is_dir():
            continue
        cname = cx.name.split("_")[0]
        for rep in reps:
            adir = cx / rep / "analysis"
            if not adir.is_dir():
                continue

            if matrix_dir is not None:
                for xpm in sorted(adir.glob("*.xpm")):
                    entry = manifest.get(xpm.name)
                    if entry is None:
                        if xpm.name not in warned_unlisted:
                            warned_unlisted.add(xpm.name)
                            print(
                                f"mdkit: manifestoda olmayan .xpm atlandi: "
                                f"{xpm.name} (ilk gorulen: {adir})",
                                file=sys.stderr,
                            )
                        continue
                    analysis, kind = entry
                    if kind != "matrix":
                        key = (xpm.name, kind)
                        if key not in warned:
                            warned.add(key)
                            print(
                                f"mdkit: .xpm ciktisi {kind!r} kind'i ile ilan "
                                f"edilmis ({xpm.name}) -- toplanmadan atlandi",
                                file=sys.stderr,
                            )
                        continue
                    # Tek bir bozuk matris 105 dizinlik toplamayi dusurmemeli
                    # (bash tarafindaki hata yalitimiyla ayni gerekce).
                    try:
                        rec = collect_matrix(
                            xpm, cname, rep, reps, analysis, matrix_dir)
                    except Exception as exc:
                        print(f"mdkit: {xpm.name} okunamadi ({adir}): {exc}",
                              file=sys.stderr)
                        continue
                    matrices.append(rec)

            for xvg in sorted(adir.glob("*.xvg")):
                ...  # MEVCUT GOVDE AYNEN KALIR
    return timeseries, profile, matrices
```

`collect`'in altına ekle:

```python
def collect_matrix(xpm, cname, rep, reps, analysis, matrix_dir):
    """Bir .xpm'i .npz olarak yazar ve ozet kaydini dondurur.

    Matrisler uzun-format CSV'ye GIRMEZ: 451x451'lik 315 matris ~64 milyon
    satir ederdi. 1D veri icin dogru olan bicim 2B icin degil."""
    meta, values, x_ps, y_ps = parse_xpm(xpm)
    peer = peer_replica(xpm.name, reps)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    out = matrix_dir / f"{cname}_{rep}_{xpm.stem}.npz"
    np.savez_compressed(
        out,
        values=values, x_ps=x_ps, y_ps=y_ps,
        unit=meta["unit"], complex=cname, replica_i=rep,
        replica_j=peer or "", analysis=analysis, output=xpm.name,
    )
    return {
        "complex": cname, "replica_i": rep, "replica_j": peer or "",
        "analysis": analysis, "output": xpm.name,
        "n_x": int(values.shape[1]), "n_y": int(values.shape[0]),
        "unit": meta["unit"],
    }
```

`main`'i güncelle:

```python
    matrix_dir = out_dir / "matrices"
    timeseries, profile, matrices = collect(
        data_root, complex_glob, reps, manifest, matrix_dir=matrix_dir)
    write_csv(out_dir / "timeseries_long.csv", TS_FIELDS, timeseries)
    write_csv(out_dir / "profile_long.csv", PR_FIELDS, profile)

    print(f"timeseries: {len(timeseries)} satir -> {out_dir / 'timeseries_long.csv'}")
    print(f"profile   : {len(profile)} satir -> {out_dir / 'profile_long.csv'}")
    print(f"matris    : {len(matrices)} dosya -> {matrix_dir}")
```

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_collect.py test_xpm.py -v`
Expected: tümü PASS

- [ ] **Step 5: Commit**

```bash
git add collect_results.py tests/test_collect.py
git commit -m "feat(collect): matrix kind'li .xpm ciktilarini topla

Once .xpm hic gorulmuyordu: glob yalnizca *.xvg idi, yani belgelenen
'taninmayan kind' uyarisi bile tetiklenmiyordu. Matrisler .npz olarak
yazilir; uzun-format CSV 451x451'lik 315 matriste ~64M satir ederdi."
```

---

## Task 5: `matrix_summary.csv`

**Files:**
- Modify: `collect_results.py` (`collect_matrix` istatistikleri + `main` yazımı)
- Test: `tests/test_collect.py`

**Interfaces:**
- Consumes: Task 4'ün `collect_matrix` fonksiyonu
- Produces: `<out_dir>/matrix_summary.csv`, kolonlar `MX_FIELDS` sırasında: `complex, replica_i, replica_j, analysis, output, n_x, n_y, min, mean, max, unit`.

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_collect.py` sonuna ekle:

```python
def test_ozet_istatistikleri(fake_dataset, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep2.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep2.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m")
    rec = [r for r in mx if r["complex"] == "last1"
           and r["replica_i"] == "rep1"][0]
    # Capraz cift (rep1 x rep2): TUM hucreler sayilir.
    # Matris: [[0,.5,1],[.5,0,.5],[1,.5,0]] -> min 0, max 1, mean 4/9
    assert rec["min"] == pytest.approx(0.0)
    assert rec["max"] == pytest.approx(1.0)
    assert rec["mean"] == pytest.approx(4.0 / 9.0, abs=1e-6)


def test_self_matriste_kosegen_haric(fake_dataset, tmp_path):
    """Self-matriste kosegen tanim geregi sifirdir ve min'i anlamsiz kilar."""
    write_outputs(fake_dataset, {"cross_rmsd_rep1.xpm": MATRIX_XPM})
    manifest = {"cross_rmsd_rep1.xpm": ("cross_rmsd", "matrix")}
    _ts, _pr, mx = collect_results.collect(
        fake_dataset, "*_pandora", ["rep1", "rep2", "rep3"], manifest,
        matrix_dir=tmp_path / "m")
    rec = [r for r in mx if r["complex"] == "last1"
           and r["replica_i"] == "rep1"][0]
    assert rec["replica_i"] == rec["replica_j"] == "rep1"
    # Kosegen (uc adet 0) haric: dort 0.5 ve iki 1.0 -> min 0.5
    assert rec["min"] == pytest.approx(0.5)
    assert rec["mean"] == pytest.approx((4 * 0.5 + 2 * 1.0) / 6.0, abs=1e-6)


def test_summary_csv_yazilir(fake_dataset, fake_config, tmp_path):
    write_outputs(fake_dataset, {"cross_rmsd_rep2.xpm": MATRIX_XPM})
    out = tmp_path / "results"
    r = subprocess.run(
        [sys.executable, "collect_results.py", "-c", str(fake_config),
         "-o", str(out)],
        capture_output=True, text=True,
        cwd=str(__import__("pathlib").Path(collect_results.__file__).parent),
    )
    assert r.returncode == 0, r.stderr
    with (out / "matrix_summary.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert rows, r.stdout
    assert set(rows[0]) == set(collect_results.MX_FIELDS)
```

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_collect.py -k "ozet or self_matris or summary" -v`
Expected: FAIL — `KeyError: 'min'` ve `matrix_summary.csv` yok

- [ ] **Step 3: İstatistikleri ve CSV yazımını uygula**

`collect_matrix` içinde `return` bloğunu şununla değiştir:

```python
    # Self-matriste kosegen TANIM GEREGI sifirdir ve min'i anlamsiz kilar;
    # tezde kullanilacak sayi capraz ciftin minimumudur.
    if peer == rep and values.shape[0] == values.shape[1]:
        sample = values[~np.eye(values.shape[0], dtype=bool)]
    else:
        sample = values.ravel()

    return {
        "complex": cname, "replica_i": rep, "replica_j": peer or "",
        "analysis": analysis, "output": xpm.name,
        "n_x": int(values.shape[1]), "n_y": int(values.shape[0]),
        "min": float(sample.min()), "mean": float(sample.mean()),
        "max": float(sample.max()), "unit": meta["unit"],
    }
```

`main` içinde `profile_long.csv` yazımının altına ekle:

```python
    write_csv(out_dir / "matrix_summary.csv", MX_FIELDS, matrices)
```

ve matris `print` satırını güncelle:

```python
    print(f"matris    : {len(matrices)} dosya -> {matrix_dir} "
          f"(+ {out_dir / 'matrix_summary.csv'})")
```

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_collect.py -v`
Expected: tümü PASS

- [ ] **Step 5: Commit**

```bash
git add collect_results.py tests/test_collect.py
git commit -m "feat(collect): matrix_summary.csv -- matris basina min/mean/max

Self-matriste istatistikler kosegen HARIC hesaplanir: kosegen tanim
geregi sifirdir ve min'i anlamsiz kilardi."
```

---

## Task 6: `--matrix` çizim modu

**Files:**
- Modify: `plot_results.py` (`load` toleransı, `plot_matrix`, `main`)
- Test: `tests/test_plot.py`

**Interfaces:**
- Consumes: Task 4'ün `.npz` dosyaları (`values`, `x_ps`, `y_ps`, `unit`, `complex`, `replica_i`, `replica_j`, `analysis`, `output`)
- Produces: `plot_results.plot_matrix(results_dir, out_dir)` → `<out_dir>/matrix/<kompleks>_<analiz>.png`

- [ ] **Step 1: Failing test'leri yaz**

`tests/test_plot.py` sonuna ekle (importlara `import numpy as np` ekle):

```python
def _write_matrix_npz(mdir, cx, rep_i, rep_j, values, unit="nm",
                      analysis="cross_rmsd"):
    mdir.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=np.float32)
    n_y, n_x = values.shape
    np.savez_compressed(
        mdir / f"{cx}_{rep_i}_{analysis}_{rep_j}.npz",
        values=values,
        x_ps=np.arange(n_x, dtype=np.float64) * 100.0,
        y_ps=np.arange(n_y, dtype=np.float64) * 100.0,
        unit=unit, complex=cx, replica_i=rep_i, replica_j=rep_j,
        analysis=analysis, output=f"{analysis}_{rep_j}.xpm",
    )


def test_matrix_modu_figur_uretir(tmp_path):
    mdir = tmp_path / "matrices"
    for i in ["rep1", "rep2"]:
        for j in ["rep1", "rep2"]:
            _write_matrix_npz(mdir, "last1", i, j, [[0.0, 0.5], [0.5, 0.0]])
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert (out / "matrix" / "last1_cross_rmsd.png").exists()


def test_matrix_dizini_yoksa_sessiz(tmp_path):
    """cross_rmsd hic kosulmamis olabilir; mod hata vermemeli."""
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert not (out / "matrix").exists()


def test_bilinmeyen_birim_cevrilmez(tmp_path, capsys):
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0]], unit="nm^2")
    plot_results.plot_matrix(tmp_path, tmp_path / "plots")
    assert "nm^2" in capsys.readouterr().err


def test_tek_kompleksin_hatasi_digerlerini_dusurmez(tmp_path, capsys):
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "iyi", "rep1", "rep1", [[0.0, 0.1], [0.1, 0.0]])
    (mdir / "bozuk_rep1_cross_rmsd_rep1.npz").write_bytes(b"npz degil")
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert (out / "matrix" / "iyi_cross_rmsd.png").exists()
    assert "bozuk" in capsys.readouterr().err


def test_matrix_only_csv_olmadan_calisir(tmp_path):
    """--matrix tek basina verildiginde CSV'ler olmasa da cikmamali."""
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0, 0.2], [0.2, 0.0]])
    r = subprocess.run(
        [sys.executable, "plot_results.py", "--results-dir", str(tmp_path),
         "--matrix"],
        capture_output=True, text=True,
        cwd=str(__import__("pathlib").Path(plot_results.__file__).parent),
    )
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "plots" / "matrix" / "last1_cross_rmsd.png").exists()
```

- [ ] **Step 2: Test'lerin başarısız olduğunu doğrula**

Run: `cd tests && python -m pytest test_plot.py -k matrix -v`
Expected: FAIL — `AttributeError: module 'plot_results' has no attribute 'plot_matrix'`

- [ ] **Step 3: `plot_matrix`'i uygula**

`plot_results.py` importlarına ekle:

```python
import numpy as np  # noqa: E402
```

Sabitlerin altına ekle:

```python
# Isi haritasi icin algisal olarak duzgun, TEK YONLU bir skala. Kategorik
# REP_COLORS paleti burada kullanilmaz: o palet ust uste binen replika
# CIZGILERINI ayirmak icin secilmisti, sirali bir buyuklugu kodlamak icin degil.
MATRIX_CMAP = "viridis"
```

`load`'u toleranslı hâle getir:

```python
def load(results_dir, required=True):
    ts_path = results_dir / "timeseries_long.csv"
    pr_path = results_dir / "profile_long.csv"
    if not ts_path.exists() and not pr_path.exists():
        # --matrix tek basina istendiginde CSV'ler olmayabilir: matrisler
        # ayri bir urun ve cross_rmsd tek basina kosulmus olabilir.
        if required:
            sys.exit(f"sonuc CSV'leri bulunamadi: {results_dir}")
        return pd.DataFrame(), pd.DataFrame()
    ts = pd.read_csv(ts_path) if ts_path.exists() else pd.DataFrame()
    pr = pd.read_csv(pr_path) if pr_path.exists() else pd.DataFrame()
    return ts, pr
```

`plot_compare`'in altına ekle:

```python
def load_matrices(results_dir):
    """results/matrices/*.npz -> {(kompleks, analiz): [kayit]}.

    Bozuk tek bir dosya digerlerini dusurmez; adiyla stderr'e yazilir."""
    mdir = results_dir / "matrices"
    groups = {}
    if not mdir.is_dir():
        return groups
    for f in sorted(mdir.glob("*.npz")):
        try:
            with np.load(f, allow_pickle=False) as d:
                rec = {
                    "values": d["values"],
                    "x_ps": d["x_ps"], "y_ps": d["y_ps"],
                    "unit": str(d["unit"]),
                    "complex": str(d["complex"]),
                    "replica_i": str(d["replica_i"]),
                    "replica_j": str(d["replica_j"]),
                    "analysis": str(d["analysis"]),
                }
        except Exception as exc:
            print(f"matris okunamadi ({f.name}): {exc}", file=sys.stderr)
            continue
        groups.setdefault((rec["complex"], rec["analysis"]), []).append(rec)
    return groups


def plot_matrix(results_dir, out_dir):
    """Kompleks basina N x N isi haritasi izgarasi, ORTAK renk skalasiyla.

    Ortak skala sart: panel basina ayri skala, farkli replika ciftlerini
    gorsel olarak karsilastirilamaz kilardi -- bu figurun tek amaci o
    karsilastirma."""
    groups = load_matrices(results_dir)
    if not groups:
        return
    target = out_dir / "matrix"
    target.mkdir(parents=True, exist_ok=True)

    for (cx, analysis), recs in sorted(groups.items()):
        fig = None
        try:
            rows = sorted({r["replica_j"] for r in recs})
            cols = sorted({r["replica_i"] for r in recs})
            by_cell = {(r["replica_i"], r["replica_j"]): r for r in recs}

            unit = recs[0]["unit"]
            conv, ulabel = scale_and_label(unit)
            if unit != "nm":
                _warn_unknown_unit(f"{cx} ({analysis})", unit)
            vmin = min(float(r["values"].min()) for r in recs) * conv
            vmax = max(float(r["values"].max()) for r in recs) * conv

            fig, axes = plt.subplots(
                len(rows), len(cols), squeeze=False,
                figsize=(2.6 * len(cols) + 1.6, 2.6 * len(rows)),
                sharex=True, sharey=True,
            )
            im = None
            for ri, rep_j in enumerate(rows):
                for ci, rep_i in enumerate(cols):
                    ax = axes[ri][ci]
                    rec = by_cell.get((rep_i, rep_j))
                    if rec is None:
                        ax.set_axis_off()
                        continue
                    x, y = rec["x_ps"], rec["y_ps"]
                    im = ax.imshow(
                        rec["values"] * conv, origin="lower", aspect="auto",
                        vmin=vmin, vmax=vmax, cmap=MATRIX_CMAP,
                        extent=[x[0] * PS_TO_NS, x[-1] * PS_TO_NS,
                                y[0] * PS_TO_NS, y[-1] * PS_TO_NS],
                    )
                    if ri == len(rows) - 1:
                        ax.set_xlabel(f"{rep_i} (ns)")
                    if ci == 0:
                        ax.set_ylabel(f"{rep_j} (ns)")
            if im is not None:
                fig.colorbar(im, ax=axes, label=f"RMSD ({ulabel})",
                             fraction=0.046, pad=0.02)
            fig.suptitle(f"{cx} — {analysis}")
            fig.savefig(target / f"{cx}_{analysis}.png", dpi=150,
                        bbox_inches="tight")
        except Exception as exc:
            print(f"matris cizimi basarisiz ({cx}, {analysis}): {exc}",
                  file=sys.stderr)
        finally:
            if fig is not None:
                plt.close(fig)
```

`main`'i güncelle:

```python
    ap.add_argument("--matrix", action="store_true")
    args = ap.parse_args()

    run_all = not (args.per_complex or args.mean_sd or args.compare
                   or args.matrix)
    needs_csv = args.per_complex or args.mean_sd or args.compare or run_all
    ts, pr = load(args.results_dir, required=needs_csv)
    groups = read_complex_groups(args.config)
    out_dir = args.results_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.per_complex or run_all:
        plot_per_complex(ts, pr, out_dir, groups)
    if args.mean_sd or run_all:
        plot_mean_sd(ts, pr, out_dir, groups)
    if args.compare or run_all:
        plot_compare(ts, out_dir, args.compare_output, groups)
    if args.matrix or run_all:
        plot_matrix(args.results_dir, out_dir)
```

`--matrix` satırı diğer mod bayraklarının hemen altına, `--compare-output`'un üstüne konur.

- [ ] **Step 4: Test'lerin geçtiğini doğrula**

Run: `cd tests && python -m pytest test_plot.py -v`
Expected: tümü PASS

- [ ] **Step 5: Tüm hızlı test paketini koştur**

Run: `cd tests && python -m pytest -q`
Expected: hiç FAIL yok

- [ ] **Step 6: Commit**

```bash
git add plot_results.py tests/test_plot.py
git commit -m "feat(plot): --matrix modu, N x N isi haritasi izgarasi

Ortak renk skalasi sart: panel basina ayri skala replika ciftlerini
karsilastirilamaz kilardi. load() artik --matrix tek basina
verildiginde CSV yoklugunda cikmiyor."
```

---

## Task 7: README

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: Task 1–6'nın tamamı
- Produces: yok (belge)

- [ ] **Step 1: "matrix henuz desteklenmiyor" bölümünü değiştir**

`README.md`'deki `### \`matrix\` henuz desteklenmiyor` başlıklı bölümün (satır ~324–332) tamamını şununla değiştir:

```markdown
### `matrix` — 2B ciktilar

`ANALYSIS_KIND="matrix"` ilan eden bir analiz `.xpm` uretir; `collect_results.py`
bunlari `results/matrices/<kompleks>_<replika>_<cikti>.npz` dosyalarina ve
`results/matrix_summary.csv` ozetine cevirir, `plot_results.py --matrix` de
kompleks basina bir isi haritasi izgarasi cizer. Ilk ornegi `cross_rmsd`.

Matrisler uzun-format CSV'ye **girmez**: 451x451'lik 315 matris ~64 milyon satir
ederdi. `timeseries_long.csv`/`profile_long.csv` icin dogru olan bicim 2B veri
icin degil.

> **`-skip` TUZAGI.** `gmx`'in `-skip`/`-skip2` secenekleri `.xpm`'in eksen
> zaman degerlerini bozar: ilk zaman dogru yazilir, geri kalani `0` olur
> (GROMACS 2025.4'te olculdu). Seyreltme **her zaman `-dt`** (ps) ile yapilir.
> `-dt` frame'leri okuma aninda eler, matris zaten seyreltilmis veriden kurulur
> ve eksenler dogru cikar.

#### Bilinen sinirlar

**Degerler 80 seviyeye yuvarlanmistir.** `.xpm` surekli degerleri `-nlevels`
(varsayilan 80) renk seviyesine indirger. Isi haritasi icin fazlasiyla yeterli,
ama `matrix_summary.csv`'deki `min`/`mean`/`max` sayilari da bu cozunurluktedir.
Tam degerler gerekirse `gmx rms -bin` ham binary dump uretir; arac bunu su an
kullanmaz.

**`-r/--reps` ile kisitlanmis kosu eksik matris uretir.** Cikti manifestosu
`config.sh`'teki `REPS`'ten turetilir, `-r` ise onu yalnizca kosu icin ezer.
Sessiz bir yanlislik degil: eksik dosya toplanmaz, sonraki tam kosu eksikligi
gorup yeniden uretir.

**DSSP gibi farkli semantikli 2B ciktilar** (residue x zaman) bu izgara
ciziminden faydalanmaz; onlar icin ayri bir cizim bicimi gerekir. Toplama
katmani ise genel: `collect_results.py` replika adiyla eslesmeyen bir matris
ciktisini da toplar, yalnizca `replica_j` kolonu bos kalir.
```

- [ ] **Step 2: "Diskte ne nereye yazilir" ağacını güncelle**

`results/` bloğunu şu hâle getir:

```
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
        └── compare_<cikti>.png
```

- [ ] **Step 3: Boru hattı şemasına matris kolunu ekle**

`### Boru hatti` bloğundaki son iki satırı şu hâle getir:

```
                                            results/{timeseries,profile}_long.csv
                                            results/matrices/*.npz + matrix_summary.csv
                                                              │
                                                              ▼
                                                     plot_results.py
                                                              │
                                                              ▼
                                                   results/plots/*.png
```

- [ ] **Step 4: "Uc mod, uc soru" tablosuna dorduncu satiri ekle**

Tablonun başlığını `### Dort mod, dort soru` yap ve şu satırı ekle:

```markdown
| `--matrix` | `plots/matrix/<kompleks>_<analiz>.png` | *Replikalar ayni konformasyonlari mi geziyor?* Kompleks basina N x N isi haritasi izgarasi, ortak renk skalasiyla. Kosegende self-matrisler (tek replika icindeki metastabil durumlar), kosegen disinda capraz ciftler. Koyu bir capraz panel, iki replikanin AYNI bolgeyi ziyaret ettigini soyler. |
```

Bölümün altındaki "Mod verilmezse **ucu birden** kosar" cümlesini
"Mod verilmezse **dordu birden** kosar" yap.

- [ ] **Step 5: Kullanim bolumune ornek ekle**

`## Kullanim` bloğunun sonuna ekle:

```bash
# Replikalar arasi 2D RMSD matrisleri (seyreltme varsayilani 200 ps)
./run_analysis.sh --all -a cross_rmsd /veri/kok
CROSS_RMSD_DT=500 ./run_analysis.sh --all -a cross_rmsd /veri/kok
```

- [ ] **Step 6: Tasarim belgeleri bolumune spec'i ekle**

`## Tasarim belgeleri` listesine üçüncü madde olarak ekle:

```markdown
- **`docs/superpowers/specs/2026-09-21-cross-rmsd-matrix-design.md`** — capraz-RMSD
  matrisleri ve `matrix` katman destegi. `.xpm` biciminin uc tuzagi (cok satirli
  eksen yorumlari, ters sirada yazilan piksel satirlari, sabit genislikli
  karakter alani) ve `-skip`/`-b` davranisi komutla dogrulanmis olarak burada.
```

- [ ] **Step 7: Portability ve tüm test paketini doğrula**

Run: `cd tests && python -m pytest -q`
Expected: hiç FAIL yok

- [ ] **Step 8: Commit**

```bash
git add README.md
git commit -m "docs(mdkit): matris boru hattini belgele, -skip tuzagini yaz

'matrix henuz desteklenmiyor' bolumu kalkti: kind artik uctan uca
tuketiliyor."
```

---

## Uygulama sonrası doğrulama

- [ ] `cd tests && python -m pytest -v` — hızlı testler
- [ ] `cd tests && python -m pytest -v -m slow` — gmx gerektirenler
- [ ] Gerçek veride tek kompleks üzerinde uçtan uca:

```bash
./run_analysis.sh -a cross_rmsd /veri/kok/<bir_kompleks>
python collect_results.py
python plot_results.py --results-dir /veri/kok/results --matrix
```

Beklenen: `results/matrices/` altında 9 `.npz`, `matrix_summary.csv`'de 9 satır,
`results/plots/matrix/<kompleks>_cross_rmsd.png` içinde 3×3 ızgara ve köşegende
sol-alttan sağ-üste uzanan koyu olmayan (sıfır) bir çizgi.
