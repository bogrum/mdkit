# mdkit — Post-MD Analiz Çerçevesi Implementasyon Planı

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rmsd_rmsf.sh`'ı, 35 kompleks × 3 replika üzerinde tekrarlanabilir koşan ve yeni analizlerin dosya ekleyerek takıldığı, çıkarılabilir bir araç setine (`mdkit/`) dönüştürmek.

**Architecture:** Bash katmanı gmx çağrılarını ve orkestrasyonu yapar; her analiz `analysis/<ad>.sh` içinde bir eklenti sözleşmesi uygular (`ANALYSIS_NAME/DESC/KIND/NEEDS_INDEX/DEFAULT_BEGIN/OUTPUTS` + `analysis_run`). Projeye özel tüm varsayımlar tek bir `config.sh`'te toplanır; koddan mutlak yol geçmez. Python katmanı yalnızca `.xvg` → birleşik CSV toplama ve çizim yapar.

**Tech Stack:** Bash 5, GROMACS 2025.4-cuda, Python (anaconda base: numpy 2.4.6, pandas 2.3.3, matplotlib 3.10.6), pytest 8.4.2.

**Spec:** `docs/superpowers/specs/2026-09-21-postmd-analysis-framework-design.md`

## Global Constraints

- **Kök dizin:** tüm yollar `mdsimulations/postmd_analysis/mdkit/` altındadır (repo kökü: `/home/emre/workspace/TUSEB-Bitirme/scratch`).
- **Mutlak yol yasağı:** `config.sh` dışında hiçbir dosyada `/mnt/`, `/usr/local/`, `/home/` geçmez. Task 11'de grep ile test edilir. (spec §3.1, §8.7)
- **Zaman birimi:** gmx çağrılarında `-tu` **kullanılmaz**; araç içinde tüm zamanlar **ps**, tüm mesafeler **nm**. Birim dönüşümü yalnızca çizim katmanında. (spec §2.5)
- **Grup seçimi:** gmx grup sorgularına numara değil **kanonik isim** yazılır: `RECEPTOR`, `AUX`, `LIGAND`, `RECEPTOR_BB`, `LIGAND_BB`. (spec §4.5)
- **Test yorumlayıcısı:** `/home/emre/anaconda3/bin/python` (numpy/pandas/matplotlib/pytest burada). PATH'teki `python3` bir venv'dir ve bilimsel yığını **yoktur** — testlerde kullanılmaz. Testler `mdkit/tests/` içinden `/home/emre/anaconda3/bin/python -m pytest` ile koşulur.
- **`set -e` yasağı:** ana döngüde `set -e` kullanılmaz; hata yalıtımı alt kabuk + manuel durum kontrolüyle yapılır. Tek bir replikanın hatası batch'i düşürmemeli. (spec §11)
- **Git:** Her task bir commit ile biter, **ancak commit komutları kullanıcının açık onayı olmadan çalıştırılmaz** (proje kuralı). Commit mesajlarına atıf/`Co-Authored-By` satırı eklenmez.
- **Gerçek veri:** `/mnt/data/scratch-simulations-TUSEB` — 35 `*_pandora` dizini, her birinde `rep1/rep2/rep3`. `check_ref.pdb` 30 replikada eksiktir; `traj_compact_center_dry.xtc` ve `md_0_10.tpr` 105/105 mevcuttur.

---

## Dosya Yapısı

| Dosya | Sorumluluk |
|---|---|
| `mdkit/config.sh` | Projeye özel **tüm** varsayımlar: yollar, dosya adları, replika listesi, zincir eşlemesi, yorumlayıcılar. Tek istisna dosya. |
| `mdkit/analysis/lib.sh` | Ortak altyapı: config yükleme/doğrulama, gmx keşfi, hedef/replika keşfi, ön koşul kontrolü, `check_ref` self-heal, index kurulumu, idempotency, log. Analiz mantığı içermez. |
| `mdkit/analysis/rmsd.sh` | RMSD eklentisi. Üç zaman serisi üretir. |
| `mdkit/analysis/rmsf.sh` | RMSF eklentisi. İki (veya `--groove-fit` ile üç) residue profili üretir. |
| `mdkit/run_analysis.sh` | Üst script: CLI ayrıştırma, `--list`, `--dry-run`, self-heal onayı, yürütme döngüsü, loglama. Analiz mantığı içermez. |
| `mdkit/collect_results.py` | `.xvg` → `results/{timeseries,profile}_long.csv`. Analiz manifestosunu `run_analysis.sh --list` çıktısından okur (tek kaynak). |
| `mdkit/plot_results.py` | Üç çizim modu: `--per-complex`, `--mean-sd`, `--compare`. |
| `mdkit/README.md` | Kurulum, config, örnek koşu, yeni analiz ekleme. |
| `mdkit/tests/conftest.py` | Ortak fixture'lar: sahte veri ağacı, sahte config, gerçek veriden 21 frame'lik küçük replika. |
| `mdkit/tests/test_*.py` | Task başına test dosyası. |

---

### Task 1: İskelet, `config.sh` ve config yükleme

**Files:**
- Create: `mdsimulations/postmd_analysis/mdkit/config.sh`
- Create: `mdsimulations/postmd_analysis/mdkit/analysis/lib.sh`
- Create: `mdsimulations/postmd_analysis/mdkit/tests/conftest.py`
- Test: `mdsimulations/postmd_analysis/mdkit/tests/test_config.py`

**Interfaces:**
- Consumes: hiçbir şey (ilk task).
- Produces: `mdkit_load_config [config_yolu]` — config'i source eder, zorunlu değişkenleri doğrular; eksikse stderr'e değişken adını yazıp 1 döner. Global `MDKIT_DIR` (mdkit kök dizini). Config değişkenleri: `DATA_ROOT`, `GMX`, `COMPLEX_GLOB`, `REPS` (dizi), `TRAJ_NAME`, `TPR_NAME`, `REF_NAME`, `RESULTS_DIR`, `PYTHON`, `CHAIN_RECEPTOR`, `CHAIN_AUX`, `CHAIN_LIGAND`.

- [ ] **Step 1: `conftest.py` ile ortak fixture'ları yaz**

`mdkit/tests/conftest.py`:

```python
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

MDKIT = Path(__file__).resolve().parents[1]
REAL_ROOT = Path("/mnt/data/scratch-simulations-TUSEB")
REAL_REP = REAL_ROOT / "last10_IMGQQPAPQV_A0201_pandora" / "rep1"
GMX_BIN = Path("/usr/local/gromacs-2025.4-cuda/bin/gmx")

needs_gmx = pytest.mark.skipif(
    not (GMX_BIN.exists() and REAL_REP.exists()),
    reason="gmx veya gercek veri yok",
)


def run_bash(snippet, **kw):
    """Bash parcacigini calistir, CompletedProcess dondur."""
    return subprocess.run(["bash", "-c", snippet], capture_output=True, text=True, **kw)


@pytest.fixture
def mdkit():
    return MDKIT


@pytest.fixture
def lib(mdkit):
    return mdkit / "analysis" / "lib.sh"


@pytest.fixture
def fake_dataset(tmp_path):
    """gmx gerektirmeyen sahte veri agaci: 2 kompleks x 3 replika."""
    root = tmp_path / "data"
    for cx in ["last1_AAA_A0201_pandora", "top1_BBB_A0201_pandora"]:
        for rep in ["rep1", "rep2", "rep3"]:
            d = root / cx / rep
            d.mkdir(parents=True)
            (d / "traj_compact_center_dry.xtc").write_text("x")
            (d / "md_0_10.tpr").write_text("x")
            (d / "check_ref.pdb").write_text("x")
    (root / "not_a_complex").mkdir()
    return root


@pytest.fixture
def fake_config(tmp_path, fake_dataset):
    """fake_dataset'i isaret eden gecerli bir config.sh."""
    cfg = tmp_path / "config.sh"
    cfg.write_text(
        textwrap.dedent(
            f"""\
            DATA_ROOT="{fake_dataset}"
            GMX=""
            COMPLEX_GLOB="*_pandora"
            REPS=(rep1 rep2 rep3)
            TRAJ_NAME="traj_compact_center_dry.xtc"
            TPR_NAME="md_0_10.tpr"
            REF_NAME="check_ref.pdb"
            RESULTS_DIR="{tmp_path}/results"
            PYTHON="{sys.executable}"
            CHAIN_RECEPTOR="A"
            CHAIN_AUX="B"
            CHAIN_LIGAND="C"
            """
        )
    )
    return cfg


@pytest.fixture(scope="session")
def small_rep(tmp_path_factory):
    """Gercek veriden 21 frame'lik kucuk replika kopyasi (oturumda bir kez)."""
    if not (GMX_BIN.exists() and REAL_REP.exists()):
        pytest.skip("gmx veya gercek veri yok")
    d = tmp_path_factory.mktemp("small_rep")
    subprocess.run(
        [
            str(GMX_BIN), "trjconv",
            "-s", str(REAL_REP / "md_0_10.tpr"),
            "-f", str(REAL_REP / "traj_compact_center_dry.xtc"),
            "-o", str(d / "traj_compact_center_dry.xtc"),
            "-b", "0", "-e", "200",
        ],
        input="Protein\n", text=True, capture_output=True, check=True,
    )
    shutil.copy(REAL_REP / "check_ref.pdb", d / "check_ref.pdb")
    os.symlink(REAL_REP / "md_0_10.tpr", d / "md_0_10.tpr")
    return d


@pytest.fixture
def real_config(tmp_path, small_rep):
    """small_rep'i tek replikali tek kompleks gibi gosteren config + veri agaci.

    DIKKAT: icerik GERCEKTEN last10_IMGQQPAPQV_A0201_pandora/rep1 verisidir,
    yalnizca dizin adi 'test1_PEPTIDE_A0201_pandora' olarak maskelenmistir
    (testler kesif mantigini gercek isimden bagimsiz dogrulasin diye).
    Bu yuzden testlerdeki "10 residue", "275 residue", "151 atom" gibi
    beklenen degerler last10'un gercek sayilaridir.
    """
    root = tmp_path / "data"
    rep = root / "test1_PEPTIDE_A0201_pandora" / "rep1"
    rep.mkdir(parents=True)
    for name in ["traj_compact_center_dry.xtc", "check_ref.pdb", "md_0_10.tpr"]:
        os.symlink(small_rep / name, rep / name)
    cfg = tmp_path / "config.sh"
    cfg.write_text(
        textwrap.dedent(
            f"""\
            DATA_ROOT="{root}"
            GMX="{GMX_BIN}"
            COMPLEX_GLOB="*_pandora"
            REPS=(rep1)
            TRAJ_NAME="traj_compact_center_dry.xtc"
            TPR_NAME="md_0_10.tpr"
            REF_NAME="check_ref.pdb"
            RESULTS_DIR="{tmp_path}/results"
            PYTHON="{sys.executable}"
            CHAIN_RECEPTOR="A"
            CHAIN_AUX="B"
            CHAIN_LIGAND="C"
            """
        )
    )
    return cfg
```

- [ ] **Step 2: Başarısız testi yaz**

`mdkit/tests/test_config.py`:

```python
from conftest import run_bash


def test_gecerli_config_yuklenir(lib, fake_config, fake_dataset):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && echo "$DATA_ROOT"'
    )
    assert r.returncode == 0, r.stderr
    assert str(fake_dataset) in r.stdout


def test_eksik_degisken_adiyla_raporlanir(lib, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text('DATA_ROOT="/x"\nREPS=(rep1)\n')
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "TRAJ_NAME" in r.stderr
    assert "REF_NAME" in r.stderr


def test_eksik_python_hata_verir(lib, tmp_path):
    """PYTHON zorunlu: PATH fallback'i yok, cunku PATH'teki python3 bilimsel
    yigina sahip olmayabilir. Erken yakalanmali."""
    bad = tmp_path / "bad.sh"
    bad.write_text(
        'DATA_ROOT=/x\nCOMPLEX_GLOB="*"\nTRAJ_NAME=a\nTPR_NAME=b\nREF_NAME=c\n'
        'RESULTS_DIR=/y\nCHAIN_RECEPTOR=A\nCHAIN_AUX=B\nCHAIN_LIGAND=C\nREPS=(rep1)\n'
    )
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "PYTHON" in r.stderr


def test_bos_reps_hata_verir(lib, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text(
        'DATA_ROOT=/x\nCOMPLEX_GLOB="*"\nTRAJ_NAME=a\nTPR_NAME=b\nREF_NAME=c\n'
        'RESULTS_DIR=/y\nPYTHON=/usr/bin/python3\n'
        'CHAIN_RECEPTOR=A\nCHAIN_AUX=B\nCHAIN_LIGAND=C\nREPS=()\n'
    )
    r = run_bash(f'source "{lib}" && mdkit_load_config "{bad}"')
    assert r.returncode != 0
    assert "REPS" in r.stderr


def test_olmayan_config_anlasilir_hata(lib, tmp_path):
    r = run_bash(f'source "{lib}" && mdkit_load_config "{tmp_path}/yok.sh"')
    assert r.returncode != 0
    assert "bulunamadi" in r.stderr


def test_varsayilan_config_script_yaninda_aranir(lib, mdkit):
    r = run_bash(f'source "{lib}" && mdkit_load_config && echo "$TRAJ_NAME"')
    assert r.returncode == 0, r.stderr
    assert "traj_compact_center_dry.xtc" in r.stdout
```

- [ ] **Step 3: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_config.py -v`
Expected: FAIL — `lib.sh` yok, `mdkit_load_config: command not found`.

- [ ] **Step 4: `config.sh`'i yaz**

`mdkit/config.sh`:

```bash
#!/usr/bin/env bash
# mdkit — projeye ozel TUM varsayimlar burada. Baska bir veri setiyle calistirmak
# icin bu dosyanin bir kopyasini duzenleyip --config ile verin; kod degismez.

# Veri kokü ve cikti
DATA_ROOT="/mnt/data/scratch-simulations-TUSEB"
RESULTS_DIR="$DATA_ROOT/results"

# Yorumlayicilar. GMX bos birakilirsa PATH'ten bulunur.
GMX="/usr/local/gromacs-2025.4-cuda/bin/gmx"
PYTHON="/home/emre/anaconda3/bin/python"

# Dizin ve dosya adlandirma
COMPLEX_GLOB="*_pandora"
REPS=(rep1 rep2 rep3)
TRAJ_NAME="traj_compact_center_dry.xtc"
TPR_NAME="md_0_10.tpr"
REF_NAME="check_ref.pdb"

# Zincir eslemesi. Analiz scriptleri bu harfleri gormez; index kurulumu
# bunlari kanonik adlara cevirir: RECEPTOR / AUX / LIGAND / RECEPTOR_BB / LIGAND_BB
CHAIN_RECEPTOR="A"   # MHC agir zincir
CHAIN_AUX="B"        # beta-2 mikroglobulin
CHAIN_LIGAND="C"     # peptid
```

- [ ] **Step 5: `lib.sh`'in config bölümünü yaz**

`mdkit/analysis/lib.sh`:

```bash
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
```

- [ ] **Step 6: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_config.py -v`
Expected: 6 passed.

- [ ] **Step 7: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add config boundary and config loader"
```

---

### Task 2: Hedef/replika keşfi ve ön koşul kontrolü

**Files:**
- Modify: `mdkit/analysis/lib.sh` (fonksiyon ekle)
- Test: `mdkit/tests/test_discovery.py`

**Interfaces:**
- Consumes: Task 1'den `mdkit_load_config`, config değişkenleri.
- Produces:
  - `mdkit_find_complexes <hedef> <all:0|1>` — satır başına bir kompleks dizini yolu yazar. `all=1` ise `<hedef>` altındaki `$COMPLEX_GLOB` eşleşen dizinler (sıralı), değilse `<hedef>`'in kendisi.
  - `mdkit_rep_status <replika_dizini>` — tek kelime yazar: `OK` | `NO_TRAJ` | `NO_TPR` | `NO_REF`.
  - `mdkit_complex_name <kompleks_dizini>` — kısa ad (`last10_IMG..._pandora` → `last10`).

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_discovery.py`:

```python
from conftest import run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && '


def test_all_modu_kompleksleri_bulur(lib, fake_config, fake_dataset):
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_find_complexes "{fake_dataset}" 1'
    )
    assert r.returncode == 0, r.stderr
    lines = r.stdout.split()
    assert len(lines) == 2
    assert lines[0].endswith("last1_AAA_A0201_pandora")
    assert lines[1].endswith("top1_BBB_A0201_pandora")
    assert not any("not_a_complex" in ln for ln in lines)


def test_tek_kompleks_modu_hedefi_aynen_dondurur(lib, fake_config, fake_dataset):
    target = fake_dataset / "last1_AAA_A0201_pandora"
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_find_complexes "{target}" 0'
    )
    assert r.stdout.strip() == str(target)


def test_rep_status_tam_replikada_ok(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "OK"


def test_rep_status_eksik_ref_bildirir(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "NO_REF"


def test_rep_status_eksik_traj_bildirir(lib, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep2"
    (rep / "traj_compact_center_dry.xtc").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_rep_status "{rep}"')
    assert r.stdout.strip() == "NO_TRAJ"


def test_complex_name_kisaltir(lib, fake_config, fake_dataset):
    cx = fake_dataset / "last1_AAA_A0201_pandora"
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + f'mdkit_complex_name "{cx}"')
    assert r.stdout.strip() == "last1"
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_discovery.py -v`
Expected: FAIL — `mdkit_find_complexes: command not found`.

- [ ] **Step 3: `lib.sh`'e keşif fonksiyonlarını ekle**

`lib.sh` sonuna:

```bash
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
    # /yol/last10_IMGQQPAPQV_A0201_pandora -> last10
    local base
    base="$(basename "${1%/}")"
    printf '%s\n' "${base%%_*}"
}
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_discovery.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add complex/replica discovery and precondition checks"
```

---

### Task 3: gmx katmanı — `check_ref` self-heal ve index kurulumu

**Files:**
- Modify: `mdkit/analysis/lib.sh`
- Test: `mdkit/tests/test_gmx_layer.py`

**Interfaces:**
- Consumes: Task 1-2'nin tümü.
- Produces:
  - `mdkit_resolve_gmx` — `$GMX` çalıştırılabilir değilse PATH'ten bulur, yoksa 1 döner.
  - `mdkit_make_ref <replika_dizini>` — `$TPR_NAME` + `$TRAJ_NAME`'den `-dump 0` ile `$REF_NAME` üretir. Başarıda 0.
  - `mdkit_build_index <replika_dizini> <cikti_dizini>` — `<cikti_dizini>/index.ndx` üretir; grupları kanonik adlara çevirir; beşinin de var olduğunu doğrular.

**Neden `md_0_10.tpr`:** `dry_reference.tpr` zincir ID'lerini taşımaz, üretilen PDB'de chain kolonu boş çıkar ve `chain A` seçimi boş grup döner (spec §2.2). `$TPR_NAME` config'de `md_0_10.tpr`'dir.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_gmx_layer.py`:

```python
import re

import pytest

from conftest import needs_gmx, run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && mdkit_resolve_gmx && '


def test_gmx_bulunamazsa_anlasilir_hata(lib, fake_config):
    r = run_bash(
        f'source "{lib}" && mdkit_load_config "{fake_config}" && '
        'PATH=/nonexistent mdkit_resolve_gmx'
    )
    assert r.returncode != 0
    assert "gmx bulunamadi" in r.stderr


@needs_gmx
def test_make_ref_zincirleri_koruyarak_uretir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_bash(PRE.format(lib=lib, cfg=real_config) + f'mdkit_make_ref "{rep}"')
    assert r.returncode == 0, r.stderr
    text = (rep / "check_ref.pdb").read_text()
    chains = {ln[21] for ln in text.splitlines() if ln.startswith("ATOM")}
    assert chains == {"A", "B", "C"}


@needs_gmx
def test_build_index_kanonik_gruplari_uretir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    r = run_bash(
        PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"'
    )
    assert r.returncode == 0, r.stderr
    names = re.findall(r"^\[ (.*) \]$", (out / "index.ndx").read_text(), re.M)
    for g in ["RECEPTOR", "AUX", "LIGAND", "RECEPTOR_BB", "LIGAND_BB"]:
        assert g in names
    assert not any(n.startswith("ch") for n in names)


@needs_gmx
def test_build_index_grup_boyutlari_dogru(lib, real_config, tmp_path):
    """LIGAND_BB, peptid residue sayisinin 3 kati atom icermeli (N, CA, C)."""
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    run_bash(PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"')

    counts, cur = {}, None
    for line in (out / "index.ndx").read_text().splitlines():
        m = re.match(r"^\[ (.*) \]$", line)
        if m:
            cur = m.group(1)
            counts[cur] = 0
        elif cur:
            counts[cur] += len(line.split())

    assert counts["LIGAND"] == 151        # last10 peptidi, 10 residue
    assert counts["LIGAND_BB"] == 30      # 10 residue x 3 backbone atomu
    assert counts["RECEPTOR_BB"] == 825   # 275 residue x 3
    assert counts["LIGAND_BB"] % 3 == 0


@needs_gmx
def test_build_index_bos_zincirde_hata_verir(lib, real_config, tmp_path):
    """Var olmayan bir zincir harfi verilirse index kurulmamali (bos grup
    sessizce kabul edilirse sonraki analizler anlamsiz sayi uretir)."""
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    cfg = tmp_path / "config_badchain.sh"
    cfg.write_text(
        real_config.read_text().replace('CHAIN_LIGAND="C"', 'CHAIN_LIGAND="Z"')
    )
    r = run_bash(PRE.format(lib=lib, cfg=cfg) + f'mdkit_build_index "{rep}" "{out}"')
    assert r.returncode != 0
    assert not (out / "index.ndx").exists()


@needs_gmx
def test_build_index_bozuk_referansta_hata_verir(lib, real_config, tmp_path):
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    out = rep / "analysis"
    out.mkdir()
    (rep / "check_ref.pdb").unlink()
    (rep / "check_ref.pdb").write_text("bu bir pdb degil\n")
    r = run_bash(
        PRE.format(lib=lib, cfg=real_config) + f'mdkit_build_index "{rep}" "{out}"'
    )
    assert r.returncode != 0
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_gmx_layer.py -v`
Expected: FAIL — `mdkit_resolve_gmx: command not found`.

- [ ] **Step 3: `lib.sh`'e gmx katmanını ekle**

`lib.sh` sonuna:

```bash
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
    local rd="$1"
    "$GMX" trjconv -s "$rd/$TPR_NAME" -f "$rd/$TRAJ_NAME" \
        -o "$rd/$REF_NAME" -dump 0 >/dev/null 2>&1 <<< "Protein"
    if [[ ! -s "$rd/$REF_NAME" ]]; then
        echo "mdkit: $REF_NAME uretilemedi: $rd" >&2
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

    "$GMX" make_ndx -f "$ref" -o "$ndx" >/dev/null 2>&1 <<EOF
chain $CHAIN_RECEPTOR
chain $CHAIN_AUX
chain $CHAIN_LIGAND
"ch$CHAIN_RECEPTOR" & "Backbone"
"ch$CHAIN_LIGAND" & "Backbone"
q
EOF

    if [[ ! -s "$ndx" ]]; then
        echo "mdkit: make_ndx basarisiz: $rep_dir" >&2
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
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_gmx_layer.py -v`
Expected: 6 passed (gmx/veri yoksa 5 skipped).

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add gmx layer with check_ref self-heal and canonical index"
```

---

### Task 4: Log, idempotency ve hata yalıtımı

**Files:**
- Modify: `mdkit/analysis/lib.sh`
- Test: `mdkit/tests/test_logging.py`

**Interfaces:**
- Consumes: Task 1-3.
- Produces:
  - `mdkit_log_init` — `$RESULTS_DIR`'i oluşturur, `$MDKIT_LOG` (= `$RESULTS_DIR/run_log.csv`) değişkenini kurar, dosya yoksa başlık satırını yazar: `timestamp,complex,replica,analysis,status,seconds,error`.
  - `mdkit_log_row <complex> <replica> <analysis> <status> <seconds> [error]` — bir satır ekler. `status` ∈ {`OK`, `SKIP_DONE`, `SKIP_MISSING`, `HATA`}.
  - `mdkit_outputs_present <cikti_dizini> <dosya...>` — hepsi varsa ve boyutu > 0 ise 0.
  - `mdkit_run_isolated <analiz_script> <rep_dir> <out_dir>` — `analysis_run`'ı alt kabukta çalıştırır; hata durumunda `MDKIT_LAST_ERROR`'a son 3 satırı koyar ve 1 döner.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_logging.py`:

```python
import csv

from conftest import run_bash

PRE = 'source "{lib}" && mdkit_load_config "{cfg}" && '


def test_log_init_baslik_yazar(lib, fake_config, tmp_path):
    r = run_bash(PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && echo "$MDKIT_LOG"')
    assert r.returncode == 0, r.stderr
    log = tmp_path / "results" / "run_log.csv"
    assert log.exists()
    assert log.read_text().splitlines()[0] == (
        "timestamp,complex,replica,analysis,status,seconds,error"
    )


def test_log_row_satir_ekler(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd OK 12 && '
          'mdkit_log_row last1 rep2 rmsd HATA 3 "gmx patladi"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert len(rows) == 2
    assert rows[0]["status"] == "OK"
    assert rows[1]["error"] == "gmx patladi"


def test_log_row_virgul_ve_yenisatiri_bozmaz(lib, fake_config, tmp_path):
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 '
          '"$(printf \'a, b\\nc\')"'
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert len(rows) == 1
    assert "a, b" in rows[0]["error"]
    assert "c" in rows[0]["error"]


def test_log_row_cift_tirnagi_standart_kacirir(lib, fake_config, tmp_path):
    """CSV standardi tirnagi ikileyerek kacirir; tek tirnaga cevirmek
    hata mesajini degistirir ve yaniltici olur."""
    run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + 'mdkit_log_init && mdkit_log_row last1 rep1 rmsd HATA 1 '
          "'gmx: \"can not find group\"'"
    )
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert rows[0]["error"] == 'gmx: "can not find group"'


def test_log_init_mevcut_dosyayi_ezmez(lib, fake_config, tmp_path):
    snippet = PRE.format(lib=lib, cfg=fake_config) + 'mdkit_log_init && mdkit_log_row a b c OK 1'
    run_bash(snippet)
    run_bash(snippet)
    lines = (tmp_path / "results" / "run_log.csv").read_text().splitlines()
    assert lines.count("timestamp,complex,replica,analysis,status,seconds,error") == 1
    assert len(lines) == 3


def test_outputs_present_hepsi_varsa_dogru(lib, fake_config, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.xvg").write_text("veri")
    (out / "b.xvg").write_text("veri")
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_outputs_present "{out}" a.xvg b.xvg'
    )
    assert r.returncode == 0


def test_outputs_present_bos_dosyayi_yok_sayar(lib, fake_config, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "a.xvg").write_text("veri")
    (out / "b.xvg").write_text("")
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config) + f'mdkit_outputs_present "{out}" a.xvg b.xvg'
    )
    assert r.returncode != 0


def test_run_isolated_hatayi_yakalar_ve_devam_eder(lib, fake_config, tmp_path):
    bad = tmp_path / "bad.sh"
    bad.write_text('analysis_run() { echo "birinci satir"; echo "patladim" >&2; return 1; }\n')
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_run_isolated "{bad}" /tmp /tmp; echo "rc=$?"; echo "err=$MDKIT_LAST_ERROR"'
    )
    assert "rc=1" in r.stdout
    assert "patladim" in r.stdout


def test_run_isolated_basarida_sifir_doner(lib, fake_config, tmp_path):
    good = tmp_path / "good.sh"
    good.write_text('analysis_run() { echo calisti; return 0; }\n')
    r = run_bash(
        PRE.format(lib=lib, cfg=fake_config)
        + f'mdkit_run_isolated "{good}" /tmp /tmp; echo "rc=$?"'
    )
    assert "rc=0" in r.stdout
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_logging.py -v`
Expected: FAIL — `mdkit_log_init: command not found`.

- [ ] **Step 3: `lib.sh`'e log/idempotency/yalıtım katmanını ekle**

`lib.sh` sonuna:

```bash
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
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_logging.py -v`
Expected: 9 passed.

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add run log, idempotency check and error isolation"
```

---

### Task 5: `run_analysis.sh` — CLI, `--help`, `--list`, `--dry-run`

**Files:**
- Create: `mdkit/run_analysis.sh`
- Modify: `mdkit/analysis/lib.sh` (analiz keşfi + metadata okuma)
- Test: `mdkit/tests/test_cli.py`

**Interfaces:**
- Consumes: Task 1-4'ün tümü.
- Produces:
  - `mdkit_analysis_scripts` — `$MDKIT_DIR/analysis/*.sh` yollarını yazar, `lib.sh` hariç.
  - `mdkit_analysis_meta <script>` — `name<TAB>kind<TAB>desc<TAB>output1,output2` satırı yazar. `collect_results.py` bu formatı okur.
  - `run_analysis.sh --list` — her analiz için bir metadata satırı (TSV). Makine tarafından ayrıştırılabilir; **analiz manifestosunun tek kaynağı budur.**
  - `run_analysis.sh --dry-run <hedef>` — `DRY-RUN <complex>/<rep>/<analiz> -> <ciktilar> (b=<ps>)` satırları.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_cli.py`:

```python
import subprocess

import pytest


def run_cli(mdkit, *args, **kw):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True, **kw,
    )


def test_help_sifirla_doner(mdkit):
    r = run_cli(mdkit, "--help")
    assert r.returncode == 0
    assert "kullanim:" in r.stdout


def test_bilinmeyen_secenek_hata(mdkit):
    r = run_cli(mdkit, "--boyle-bir-sey-yok")
    assert r.returncode == 2
    assert "bilinmeyen secenek" in r.stderr


def test_hedefsiz_cagri_hata(mdkit, fake_config):
    r = run_cli(mdkit, "--config", str(fake_config))
    assert r.returncode == 2
    assert "hedef" in r.stderr.lower()


def test_list_tsv_uretir(mdkit, fake_config):
    r = run_cli(mdkit, "--config", str(fake_config), "--list")
    assert r.returncode == 0, r.stderr
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert len(lines) >= 2
    names = set()
    for ln in lines:
        parts = ln.split("\t")
        assert len(parts) == 4, f"4 alan bekleniyordu: {ln!r}"
        name, kind, desc, outputs = parts
        assert kind in {"timeseries", "profile", "matrix"}
        assert desc.strip()
        assert outputs.strip()
        names.add(name)
    assert {"rmsd", "rmsf"} <= names


def test_bilinmeyen_analiz_adi_hata(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "-a", "yokboylebiranaliz",
        str(fake_dataset / "last1_AAA_A0201_pandora"),
    )
    assert r.returncode == 2
    assert "bilinmeyen analiz" in r.stderr


def test_dry_run_hicbir_sey_yazmaz(mdkit, fake_config, fake_dataset, tmp_path):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "--all", str(fake_dataset)
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.count("DRY-RUN") == 2 * 3 * 2   # 2 kompleks x 3 replika x 2 analiz
    assert "last1/rep1/rmsd" in r.stdout
    assert not (tmp_path / "results" / "run_log.csv").exists()
    assert not (fake_dataset / "last1_AAA_A0201_pandora" / "rep1" / "analysis").exists()


def test_dry_run_tek_kompleks(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == 3 * 2
    assert "last1" not in r.stdout


def test_reps_secenegi_daraltir(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-r", "rep1",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert r.stdout.count("DRY-RUN") == 2
    assert "rep2" not in r.stdout


def test_begin_secenegi_dry_run_ciktisinda_gorunur(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-b", "5000", "-a", "rmsf",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=5000)" in r.stdout


def test_begin_verilmezse_analiz_varsayilani_kullanilir(mdkit, fake_config, fake_dataset):
    r = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-a", "rmsf",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=10000)" in r.stdout
    r2 = run_cli(
        mdkit, "--config", str(fake_config), "--dry-run", "-a", "rmsd",
        str(fake_dataset / "top1_BBB_A0201_pandora"),
    )
    assert "(b=0)" in r2.stdout
```

**Not:** Bu testler Task 7-8'de yazılacak `rmsd.sh` ve `rmsf.sh`'a bağlıdır. Task 5'i bitirmek için `analysis/rmsd.sh` ve `analysis/rmsf.sh`'ın **yalnızca metadata bloklarını** (Step 4) oluşturun; `analysis_run` gövdeleri Task 7-8'de doldurulur.

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_cli.py -v`
Expected: FAIL — `run_analysis.sh` yok.

- [ ] **Step 3: `lib.sh`'e analiz keşfi ve metadata okuma ekle**

`lib.sh` sonuna:

```bash
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
```

- [ ] **Step 4: Analiz metadata iskeletlerini oluştur**

`mdkit/analysis/rmsd.sh`:

```bash
#!/usr/bin/env bash
# RMSD eklentisi. run_analysis.sh tarafindan source edilir.
ANALYSIS_NAME="rmsd"
ANALYSIS_DESC="Peptid ve MHC RMSD zaman serileri (oluk-uzeri, ic, global)"
ANALYSIS_KIND="timeseries"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=0
ANALYSIS_OUTPUTS=(rmsd_pep_on_mhc.xvg rmsd_pep_internal.xvg rmsd_complex_bb.xvg)

analysis_run() {
    echo "rmsd: Task 7'de doldurulacak" >&2
    return 1
}
```

`mdkit/analysis/rmsf.sh`:

```bash
#!/usr/bin/env bash
# RMSF eklentisi. run_analysis.sh tarafindan source edilir.
ANALYSIS_NAME="rmsf"
ANALYSIS_DESC="Residue bazli RMSF profilleri (peptid ic, MHC, opsiyonel oluk-cercevesi)"
ANALYSIS_KIND="profile"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
ANALYSIS_OUTPUTS=(rmsf_pep_self.xvg rmsf_mhc.xvg)

analysis_run() {
    echo "rmsf: Task 8'de doldurulacak" >&2
    return 1
}
```

- [ ] **Step 5: `run_analysis.sh`'ın CLI ve keşif bölümünü yaz**

`mdkit/run_analysis.sh`:

```bash
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
```

- [ ] **Step 6: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_cli.py -v`
Expected: 10 passed.

- [ ] **Step 7: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add runner CLI with --list manifest and --dry-run"
```

---

### Task 6: `run_analysis.sh` — yürütme döngüsü, self-heal onayı, loglama

**Files:**
- Modify: `mdkit/run_analysis.sh` (son iki `echo` satırını gerçek döngüyle değiştir)
- Modify: `mdkit/config.sh` (opsiyonel `POSTMD_SCRIPT` ekle)
- Test: `mdkit/tests/test_runner.py`

**Interfaces:**
- Consumes: Task 1-5'in tümü.
- Produces: Analiz script'lerine geçirilen ortam değişkenleri — `GMX`, `REF_PDB`, `XTC`, `NDX`, `B_PS`, `GROOVE_FIT`. Task 7-8 bunlara güvenir.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_runner.py`:

```python
import csv
import subprocess
import textwrap

from conftest import needs_gmx


def run_cli(mdkit, *args, **kw):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True, **kw,
    )


def _stub_analysis(mdkit, name, body):
    """Gecici bir analiz eklentisi yazar; test sonunda silinir."""
    p = mdkit / "analysis" / f"{name}.sh"
    p.write_text(body)
    return p


def test_eksik_onkosul_skip_missing_loglanir(mdkit, fake_config, fake_dataset, tmp_path):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "traj_compact_center_dry.xtc").unlink()
    run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "rmsd",
            str(fake_dataset / "last1_AAA_A0201_pandora"))
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    skipped = [r for r in rows if r["status"] == "SKIP_MISSING"]
    assert len(skipped) == 1
    assert skipped[0]["replica"] == "rep1"
    assert "NO_TRAJ" in skipped[0]["error"]


def test_bir_replikanin_hatasi_digerlerini_durdurmaz(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() {
            if [[ "$1" == *rep2* ]]; then echo "kasitli hata" >&2; return 1; fi
            echo "veri" > "$2/zz.xvg"
        }
        """))
    try:
        r = run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "zzstub",
                    str(fake_dataset / "last1_AAA_A0201_pandora"))
        assert r.returncode == 0, r.stderr
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        by_rep = {r["replica"]: r["status"] for r in rows}
        assert by_rep == {"rep1": "OK", "rep2": "HATA", "rep3": "OK"}
        assert "kasitli hata" in [r["error"] for r in rows if r["status"] == "HATA"][0]
    finally:
        stub.unlink()


def test_idempotency_ikinci_kosuda_skip_done(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() { echo "veri" > "$2/zz.xvg"; }
        """))
    try:
        args = ["-c", str(fake_config), "-y", "-a", "zzstub",
                str(fake_dataset / "last1_AAA_A0201_pandora")]
        run_cli(mdkit, *args)
        run_cli(mdkit, *args)
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert [r["status"] for r in rows[:3]] == ["OK", "OK", "OK"]
        assert [r["status"] for r in rows[3:]] == ["SKIP_DONE"] * 3
    finally:
        stub.unlink()


def test_force_idempotencyi_ezer(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzstub", textwrap.dedent("""\
        ANALYSIS_NAME="zzstub"
        ANALYSIS_DESC="test icin"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(zz.xvg)
        analysis_run() { echo "veri" > "$2/zz.xvg"; }
        """))
    try:
        args = ["-c", str(fake_config), "-y", "-a", "zzstub",
                str(fake_dataset / "last1_AAA_A0201_pandora")]
        run_cli(mdkit, *args)
        run_cli(mdkit, *args, "--force")
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert [r["status"] for r in rows[3:]] == ["OK"] * 3
    finally:
        stub.unlink()


def test_analiz_scripti_ortam_degiskenlerini_gorur(mdkit, fake_config, fake_dataset, tmp_path):
    stub = _stub_analysis(mdkit, "zzenv", textwrap.dedent("""\
        ANALYSIS_NAME="zzenv"
        ANALYSIS_DESC="ortam testi"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(env.xvg)
        analysis_run() {
            printf 'REF=%s\\nXTC=%s\\nB=%s\\nGF=%s\\n' \\
                "$REF_PDB" "$XTC" "$B_PS" "$GROOVE_FIT" > "$2/env.xvg"
        }
        """))
    try:
        run_cli(mdkit, "-c", str(fake_config), "-y", "-a", "zzenv", "-r", "rep1",
                "-b", "777", "--groove-fit",
                str(fake_dataset / "last1_AAA_A0201_pandora"))
        out = (fake_dataset / "last1_AAA_A0201_pandora" / "rep1" / "analysis" / "env.xvg").read_text()
        assert "check_ref.pdb" in out
        assert "traj_compact_center_dry.xtc" in out
        assert "B=777" in out
        assert "GF=1" in out
    finally:
        stub.unlink()


def test_eksik_ref_onay_reddedilirse_uretilmez(mdkit, fake_config, fake_dataset):
    rep = fake_dataset / "last1_AAA_A0201_pandora" / "rep1"
    (rep / "check_ref.pdb").unlink()
    r = run_cli(mdkit, "-c", str(fake_config), "-a", "rmsd",
                str(fake_dataset / "last1_AAA_A0201_pandora"), input="h\n")
    assert "1 replikada" in r.stdout
    assert not (rep / "check_ref.pdb").exists()


@needs_gmx
def test_eksik_ref_onaylanirsa_uretilir_ve_analiz_devam_eder(mdkit, real_config, tmp_path):
    """DIKKAT: bu test bir STUB eklenti kullanir, `-a rmsd` DEGIL. Task 6'nin
    dogrulamasi gereken sey self-heal akisi ve dongunun kaldigi yerden devam
    etmesidir; rmsd.sh'in govdesi Task 7'ye kadar iskeledir ve `-a rmsd`
    kullanmak bu testi Task 7'ye bagimli kilar."""
    stub = _stub_analysis(mdkit, "zzheal", textwrap.dedent("""\
        ANALYSIS_NAME="zzheal"
        ANALYSIS_DESC="self-heal testi"
        ANALYSIS_KIND="timeseries"
        ANALYSIS_NEEDS_INDEX=0
        ANALYSIS_DEFAULT_BEGIN=0
        ANALYSIS_OUTPUTS=(heal.xvg)
        analysis_run() { echo "veri" > "$2/heal.xvg"; }
        """))
    try:
        rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
        (rep / "check_ref.pdb").unlink()
        r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "zzheal",
                    str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
        assert (rep / "check_ref.pdb").exists(), r.stderr
        rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
        assert any(row["status"] == "OK" for row in rows), rows
    finally:
        stub.unlink()
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_runner.py -v`
Expected: FAIL — `run_log.csv` oluşmuyor; döngü henüz yok.

- [ ] **Step 3: `config.sh`'e opsiyonel `POSTMD_SCRIPT` ekle**

`config.sh` sonuna:

```bash
# Opsiyonel: --full-postmd icin tam post-MD scripti. Bos birakilirsa
# --full-postmd hata verir ve yalnizca hizli yol (-dump 0) kullanilabilir.
POSTMD_SCRIPT="/home/emre/workspace/TUSEB-Bitirme/scratch/mdsimulations/system_prep/post_md_script.sh"
```

- [ ] **Step 4: `run_analysis.sh`'ın son iki `echo` satırını yürütme döngüsüyle değiştir**

`run_analysis.sh` içindeki

```bash
echo "mdkit: ${#complexes[@]} kompleks, ${#REPS[@]} replika, ${#selected[@]} analiz"
echo "mdkit: yurutme dongusu Task 6'da eklenecek"
```

satırlarını şununla değiştir:

```bash
mdkit_resolve_gmx || exit 1
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
                if [[ -z "${POSTMD_SCRIPT:-}" || ! -f "${POSTMD_SCRIPT:-}" ]]; then
                    echo "--full-postmd icin config'de POSTMD_SCRIPT tanimli olmali" >&2
                    exit 2
                fi
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
            mdkit_log_row "$cname" "$rep" "-" SKIP_MISSING 0 "$status"
            continue
        fi

        out_dir="$rep_dir/analysis"
        mkdir -p "$out_dir"

        for script in "${selected[@]}"; do
            # shellcheck source=/dev/null
            source "$script"
            begin_ps="${BEGIN_OPT:-$ANALYSIS_DEFAULT_BEGIN}"

            if [[ "$FORCE" != 1 ]] && \
               mdkit_outputs_present "$out_dir" "${ANALYSIS_OUTPUTS[@]}"; then
                echo "$cname/$rep/$ANALYSIS_NAME: mevcut -- atlaniyor"
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" SKIP_DONE 0
                continue
            fi

            if [[ "${ANALYSIS_NEEDS_INDEX:-0}" == 1 && ! -s "$out_dir/index.ndx" ]]; then
                if ! mdkit_build_index "$rep_dir" "$out_dir"; then
                    mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" HATA 0 "index kurulamadi"
                    continue
                fi
            fi

            # Analiz scriptlerinin gordugu sozlesme degiskenleri
            REF_PDB="$rep_dir/$REF_NAME"
            XTC="$rep_dir/$TRAJ_NAME"
            NDX="$out_dir/index.ndx"
            B_PS="$begin_ps"

            t0=$SECONDS
            if mdkit_run_isolated "$script" "$rep_dir" "$out_dir"; then
                echo "$cname/$rep/$ANALYSIS_NAME: tamam ($((SECONDS - t0))s)"
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" OK "$((SECONDS - t0))"
            else
                echo "$cname/$rep/$ANALYSIS_NAME: HATA -- $MDKIT_LAST_ERROR" >&2
                mdkit_log_row "$cname" "$rep" "$ANALYSIS_NAME" HATA \
                    "$((SECONDS - t0))" "$MDKIT_LAST_ERROR"
            fi
        done
    done
done

echo "mdkit: bitti. log -> $MDKIT_LOG"
```

- [ ] **Step 5: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_runner.py -v`
Expected: 7 passed (gmx/veri yoksa 1 skipped).

- [ ] **Step 6: Önceki testlerin hâlâ geçtiğini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest -v`
Expected: tüm testler passed/skipped, hiç failed yok.

- [ ] **Step 7: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add execution loop with self-heal prompt and per-replica logging"
```

---

### Task 7: `analysis/rmsd.sh`

**Files:**
- Modify: `mdkit/analysis/rmsd.sh` (Task 5'teki iskeletin `analysis_run`'ını doldur)
- Test: `mdkit/tests/test_rmsd.py`

**Interfaces:**
- Consumes: `GMX`, `REF_PDB`, `XTC`, `NDX`, `B_PS` (Task 6'dan); `index.ndx` kanonik grupları (Task 3'ten).
- Produces: `rmsd_pep_on_mhc.xvg`, `rmsd_pep_internal.xvg`, `rmsd_complex_bb.xvg` — ps/nm cinsinden, `@ subtitle` alanında fit/hesap grubu yazılı.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_rmsd.py`:

```python
import subprocess

from conftest import needs_gmx

OUTPUTS = ["rmsd_pep_on_mhc.xvg", "rmsd_pep_internal.xvg", "rmsd_complex_bb.xvg"]


def run_cli(mdkit, *args):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True,
    )


def data_rows(path):
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


def header(path, key):
    for ln in path.read_text().splitlines():
        if ln.startswith(f"@") and key in ln:
            return ln
    return ""


@needs_gmx
def test_uc_cikti_uretilir(mdkit, real_config, tmp_path):
    r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
                str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    assert r.returncode == 0, r.stderr
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    for name in OUTPUTS:
        assert (out / name).exists(), f"{name} uretilmedi. stderr={r.stderr[-800:]}"


@needs_gmx
def test_frame_sayisi_dogru(mdkit, real_config, tmp_path):
    """small_rep 0-200 ps, dt=10 ps -> 21 frame."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    for name in OUTPUTS:
        assert len(data_rows(out / name)) == 21, name


@needs_gmx
def test_zaman_ps_cinsinden(mdkit, real_config, tmp_path):
    """-tu kullanilmadigi icin x ekseni ps olmali: son satir 200 civari."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    f = out / "rmsd_pep_on_mhc.xvg"
    last_x = float(data_rows(f)[-1].split()[0])
    assert abs(last_x - 200.0) < 1e-6
    assert "(ps)" in header(f, "xaxis")


@needs_gmx
def test_subtitle_gruplari_belgeler(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert "LIGAND after lsq fit to RECEPTOR_BB" in header(
        out / "rmsd_pep_on_mhc.xvg", "subtitle")
    assert "LIGAND after lsq fit to LIGAND_BB" in header(
        out / "rmsd_pep_internal.xvg", "subtitle")


@needs_gmx
def test_ilk_frame_sifira_yakin(mdkit, real_config, tmp_path):
    """Referans frame 0 oldugu icin ilk RMSD ~0 olmali."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    first = float(data_rows(out / "rmsd_pep_on_mhc.xvg")[0].split()[1])
    assert first < 0.01


@needs_gmx
def test_begin_secenegi_frameleri_kirpar(mdkit, real_config, tmp_path):
    """-b 100 (ps) -> 100..200 ps arasi 11 frame."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd", "-b", "100",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    rows = data_rows(out / "rmsd_pep_on_mhc.xvg")
    assert len(rows) == 11
    assert abs(float(rows[0].split()[0]) - 100.0) < 1e-6


@needs_gmx
def test_ic_ve_toplam_rmsd_farkli_fit_grubundan_gelir(mdkit, real_config, tmp_path):
    """rmsd_pep_on_mhc (fit RECEPTOR_BB) ve rmsd_pep_internal (fit LIGAND_BB)
    farkli fit gruplari kullanir: ayni sayida frame uretmeliler ama egrileri
    birebir AYNI OLMAMALI. Aynilarsa iki cagriya ayni fit grubu verilmis
    demektir.

    NOT: burada "ic RMSD toplamdan kucuktur" gibi bir FIZIK beklentisi
    dogrulanmaz. gmx rms rijit donusumu FIT grubu uzerinde minimize eder ama
    RMSD'yi HESAP grubu uzerinde raporlar; iki kume farkli oldugu icin
    esitsizlik matematiksel olarak garanti degildir ve baska bir veri setinde
    cikti tamamen dogruyken test kirmiziya donebilir."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsd",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    total = [float(ln.split()[1]) for ln in data_rows(out / "rmsd_pep_on_mhc.xvg")]
    internal = [float(ln.split()[1]) for ln in data_rows(out / "rmsd_pep_internal.xvg")]
    assert len(total) == len(internal)
    assert total != internal
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_rmsd.py -v`
Expected: FAIL — `analysis_run` hâlâ iskelet, "Task 7'de doldurulacak" hatası.

- [ ] **Step 3: `rmsd.sh`'ın `analysis_run`'ını yaz**

`mdkit/analysis/rmsd.sh` içindeki iskelet `analysis_run`'ı şununla değiştir:

```bash
analysis_run() {
    local rep_dir="$1" out_dir="$2"

    # -tu KULLANILMAZ: -tu ns, -b/-e degerlerini de ns'e cevirir (spec 2.5).
    # Araç icinde zaman daima ps, mesafe daima nm.

    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_pep_on_mhc.xvg" -b "$B_PS" <<< $'RECEPTOR_BB\nLIGAND' \
        || { echo "gmx rms basarisiz: rmsd_pep_on_mhc"; return 1; }

    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_pep_internal.xvg" -b "$B_PS" <<< $'LIGAND_BB\nLIGAND' \
        || { echo "gmx rms basarisiz: rmsd_pep_internal"; return 1; }

    "$GMX" rms -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsd_complex_bb.xvg" -b "$B_PS" <<< $'Backbone\nBackbone' \
        || { echo "gmx rms basarisiz: rmsd_complex_bb"; return 1; }

    return 0
}
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_rmsd.py -v`
Expected: 7 passed.

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): implement rmsd analysis with three fit/calc variants"
```

---

### Task 8: `analysis/rmsf.sh`

**Files:**
- Modify: `mdkit/analysis/rmsf.sh`
- Test: `mdkit/tests/test_rmsf.py`

**Interfaces:**
- Consumes: Task 6'nın sözleşme değişkenleri + `GROOVE_FIT`.
- Produces: `rmsf_pep_self.xvg`, `rmsf_mhc.xvg`, ve `GROOVE_FIT=1` ise `rmsf_pep_groovefit.xvg` + yeniden kullanılabilir `traj_fit_mhc.xtc`.

**Kritik ayrım (spec §2.4):** `gmx rmsf` varsayılan olarak fit yapar ve fit'i *seçilen grubun üzerinde* yapar. `LIGAND` seçmek peptidin **iç** esnekliğini verir. Oluk çerçevesindeki esneklik için trajektori önce `RECEPTOR_BB`'ye fitlenip `-nofit` kullanılmalıdır.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_rmsf.py`:

```python
import subprocess

from conftest import needs_gmx


def run_cli(mdkit, *args):
    return subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), *args],
        capture_output=True, text=True,
    )


def data_rows(path):
    return [ln for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


@needs_gmx
def test_iki_profil_uretilir(mdkit, real_config, tmp_path):
    r = run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
                str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    assert r.returncode == 0, r.stderr
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert (out / "rmsf_pep_self.xvg").exists(), r.stderr[-800:]
    assert (out / "rmsf_mhc.xvg").exists()


@needs_gmx
def test_peptid_profili_residue_sayisi_kadar(mdkit, real_config, tmp_path):
    """last10 peptidi IMGQQPAPQV -> 10 residue."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert len(data_rows(out / "rmsf_pep_self.xvg")) == 10


@needs_gmx
def test_mhc_profili_275_residue(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert len(data_rows(out / "rmsf_mhc.xvg")) == 275


@needs_gmx
def test_groove_fit_ucuncu_cikti_ve_fitli_traj_uretir(mdkit, real_config, tmp_path):
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    rep = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1"
    assert (rep / "analysis" / "rmsf_pep_groovefit.xvg").exists()
    assert (rep / "traj_fit_mhc.xtc").exists()
    assert len(data_rows(rep / "analysis" / "rmsf_pep_groovefit.xvg")) == 10


@needs_gmx
def test_groove_ve_self_farkli_hesaplamadan_gelir(mdkit, real_config, tmp_path):
    """Oluk cercevesi profili (RECEPTOR_BB'ye fitli traj + -nofit) ile kendi
    uzerine fitlenmis profil AYNI uzunlukta olmali ama birebir AYNI OLMAMALI.
    Aynilarsa groove-fit yolu self-fit hesabina cokmus demektir.

    GOZLEM (dogrulama degil): bu fixture'da self toplam 0.9583, groove toplam
    1.1121 (~%16 marj) -- yani "oluk cercevesi daha esnek" beklentisi toplamda
    tutuyor. Ama residue bazinda TERS DONUYOR (residue 1: self 0.087 > groove
    0.080) ve gmx rmsf'in fit'i referansa-sapmayi minimize ederken raporlanan
    buyukluk ortalama-etrafindaki-varyanstir; bu yuzden esitsizlik bir kod
    degismezi degil, istatistiksel bir beklentidir ve kapi olarak kullanilmaz."""
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit",
            str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"))
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    self_v = [float(ln.split()[1]) for ln in data_rows(out / "rmsf_pep_self.xvg")]
    groove = [float(ln.split()[1]) for ln in data_rows(out / "rmsf_pep_groovefit.xvg")]
    assert len(self_v) == len(groove)
    assert self_v != groove


@needs_gmx
def test_groove_fit_idempotency_ucuncu_ciktiyi_da_sayar(mdkit, real_config, tmp_path):
    """--groove-fit ile ANALYSIS_OUTPUTS uc elemanli olmali; aksi halde
    ucuncu cikti hic uretilmeden SKIP_DONE verilirdi."""
    import csv
    target = str(tmp_path / "data" / "test1_PEPTIDE_A0201_pandora")
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0", target)
    run_cli(mdkit, "-c", str(real_config), "-y", "-a", "rmsf", "-b", "0",
            "--groove-fit", target)
    out = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora" / "rep1" / "analysis"
    assert (out / "rmsf_pep_groovefit.xvg").exists()
    rows = list(csv.DictReader((tmp_path / "results" / "run_log.csv").open()))
    assert rows[-1]["status"] == "OK"
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_rmsf.py -v`
Expected: FAIL — iskelet `analysis_run` 1 dönüyor.

- [ ] **Step 3: `rmsf.sh`'ı yaz**

`mdkit/analysis/rmsf.sh` dosyasının tamamını şununla değiştir:

```bash
#!/usr/bin/env bash
# RMSF eklentisi. run_analysis.sh tarafindan source edilir.
#
# DIKKAT (spec 2.4): gmx rmsf varsayilan olarak fit yapar ve fit'i SECILEN
# GRUBUN uzerinde yapar. Bu yuzden:
#   LIGAND + -fit            -> peptidin IC esnekligi (kayma/sallanma cikarilmis)
#   RECEPTOR_BB'ye fitli traj + LIGAND + -nofit -> OLUK CERCEVESINDEKI esneklik
ANALYSIS_NAME="rmsf"
ANALYSIS_DESC="Residue bazli RMSF profilleri (peptid ic, MHC, opsiyonel oluk-cercevesi)"
ANALYSIS_KIND="profile"
ANALYSIS_NEEDS_INDEX=1
ANALYSIS_DEFAULT_BEGIN=10000
ANALYSIS_OUTPUTS=(rmsf_pep_self.xvg rmsf_mhc.xvg)

# --groove-fit ile ucuncu cikti eklenir; idempotency kontrolu onu da saymali.
if [[ "${GROOVE_FIT:-0}" == "1" ]]; then
    ANALYSIS_OUTPUTS+=(rmsf_pep_groovefit.xvg)
fi

analysis_run() {
    local rep_dir="$1" out_dir="$2"

    "$GMX" rmsf -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsf_pep_self.xvg" -res -b "$B_PS" <<< 'LIGAND' \
        || { echo "gmx rmsf basarisiz: rmsf_pep_self"; return 1; }

    "$GMX" rmsf -s "$REF_PDB" -f "$XTC" -n "$NDX" \
        -o "$out_dir/rmsf_mhc.xvg" -res -b "$B_PS" <<< 'RECEPTOR' \
        || { echo "gmx rmsf basarisiz: rmsf_mhc"; return 1; }

    if [[ "${GROOVE_FIT:-0}" == "1" ]]; then
        local fitted="$rep_dir/traj_fit_mhc.xtc"
        # Fitli trajektori saklanir: ileride PCA/DCCM ayni dosyayi isteyecek.
        if [[ ! -s "$fitted" ]]; then
            "$GMX" trjconv -s "$REF_PDB" -f "$XTC" -n "$NDX" \
                -fit rot+trans -o "$fitted" <<< $'RECEPTOR_BB\nSystem' \
                || { echo "gmx trjconv -fit basarisiz"; return 1; }
        fi
        "$GMX" rmsf -s "$REF_PDB" -f "$fitted" -n "$NDX" \
            -o "$out_dir/rmsf_pep_groovefit.xvg" -res -nofit -b "$B_PS" <<< 'LIGAND' \
            || { echo "gmx rmsf basarisiz: rmsf_pep_groovefit"; return 1; }
    fi

    return 0
}
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_rmsf.py -v`
Expected: 6 passed.

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): implement rmsf analysis with optional groove-frame fitting"
```

---

### Task 9: `collect_results.py` — `.xvg` → birleşik uzun-format CSV

**Files:**
- Create: `mdkit/collect_results.py`
- Test: `mdkit/tests/test_collect.py`

**Interfaces:**
- Consumes: `run_analysis.sh --list` TSV manifestosu (Task 5), `analysis/lib.sh`'in `mdkit_load_config`'i (Task 1), `rep*/analysis/*.xvg` (Task 7-8).
- Produces:
  - `parse_xvg(path) -> (meta: dict, rows: list[list[float]])` — `meta` anahtarları: `legends` (dict[int, str]), `title`, `subtitle`, `xaxis`, `yaxis` (hepsi opsiyonel).
  - `read_manifest(config) -> dict[str, tuple[str, str]]` — çıktı dosya adı → (analiz adı, KIND).
  - `collect(data_root, complex_glob, reps, manifest) -> (timeseries: list[dict], profile: list[dict])`.
  - `results/timeseries_long.csv` kolonları: `complex,replica,analysis,output,series,time_ps,value,unit`
  - `results/profile_long.csv` kolonları: `complex,replica,analysis,output,series,residue,value,unit`

**Birim kararı:** Değerler `.xvg`'deki ham hâliyle (ps, nm) yazılır, birim `unit` kolonunda belirtilir. Körlemesine nm→Å çevirmek ileride `gmx sasa` (nm²) eklendiğinde sessizce yanlış olurdu. Dönüşüm çizim katmanında.

- [ ] **Step 1: Başarısız testi yaz**

`mdkit/tests/test_collect.py`:

```python
import csv
import subprocess
import sys
import textwrap

import pytest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
import collect_results  # noqa: E402

TIMESERIES_XVG = textwrap.dedent("""\
    # gmx rms tarafindan uretildi
    @    title "RMSD"
    @    xaxis  label "Time (ps)"
    @    yaxis  label "RMSD (nm)"
    @ subtitle "LIGAND after lsq fit to RECEPTOR_BB"
    @TYPE xy
       0.0000000    0.0000005
      10.0000000    0.1457492
      20.0000000    0.1278536
    """)

PROFILE_XVG = textwrap.dedent("""\
    @    title "RMS fluctuation"
    @    xaxis  label "Residue"
    @    yaxis  label "(nm)"
    @TYPE xy
       1   0.1823
       2   0.1490
    """)

MULTI_SERIES_XVG = textwrap.dedent("""\
    @    xaxis  label "Time (ps)"
    @    yaxis  label "RMSD (nm)"
    @ s0 legend "birinci"
    @ s1 legend "ikinci"
       0.0   0.10   0.20
      10.0   0.11   0.21
    """)


def write_outputs(fake_dataset, mapping):
    """mapping: {dosya_adi: icerik} -> her kompleks/replikaya yazar."""
    for cx in fake_dataset.iterdir():
        if not cx.name.endswith("_pandora"):
            continue
        for rep in ["rep1", "rep2", "rep3"]:
            adir = cx / rep / "analysis"
            adir.mkdir(parents=True, exist_ok=True)
            for name, content in mapping.items():
                (adir / name).write_text(content)


def test_parse_xvg_yorum_ve_basliklari_atlar(tmp_path):
    f = tmp_path / "a.xvg"
    f.write_text(TIMESERIES_XVG)
    meta, rows = collect_results.parse_xvg(f)
    assert len(rows) == 3
    assert rows[1] == [10.0, 0.1457492]
    assert meta["yaxis"] == "RMSD (nm)"
    assert meta["subtitle"] == "LIGAND after lsq fit to RECEPTOR_BB"


def test_parse_xvg_legend_okur(tmp_path):
    f = tmp_path / "b.xvg"
    f.write_text(MULTI_SERIES_XVG)
    meta, rows = collect_results.parse_xvg(f)
    assert meta["legends"] == {0: "birinci", 1: "ikinci"}
    assert rows[0] == [0.0, 0.10, 0.20]


def test_manifest_runner_listesinden_okunur(mdkit, fake_config, monkeypatch):
    monkeypatch.setattr(collect_results, "MDKIT", mdkit)
    manifest = collect_results.read_manifest(fake_config)
    assert manifest["rmsd_pep_on_mhc.xvg"] == ("rmsd", "timeseries")
    assert manifest["rmsf_pep_self.xvg"] == ("rmsf", "profile")


def run_collect(mdkit, config):
    return subprocess.run(
        [sys.executable, str(mdkit / "collect_results.py"), "-c", str(config)],
        capture_output=True, text=True,
    )


def test_timeseries_ve_profile_ayri_dosyalara_gider(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "rmsf_pep_self.xvg": PROFILE_XVG,
    })
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr

    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    pr = list(csv.DictReader((tmp_path / "results" / "profile_long.csv").open()))
    assert len(ts) == 2 * 3 * 3   # 2 kompleks x 3 replika x 3 satir
    assert len(pr) == 2 * 3 * 2
    assert ts[0]["analysis"] == "rmsd"
    assert ts[0]["time_ps"] == "0.0"
    assert ts[0]["unit"] == "nm"
    assert pr[0]["residue"] == "1"
    assert pr[0]["analysis"] == "rmsf"


def test_ham_deger_korunur_donusturulmez(mdkit, fake_config, fake_dataset, tmp_path):
    """nm->A donusumu burada YAPILMAZ; ileride sasa (nm^2) eklendiginde
    korlemesine carpma sessizce yanlis olurdu."""
    write_outputs(fake_dataset, {"rmsd_pep_on_mhc.xvg": TIMESERIES_XVG})
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    assert float(ts[1]["value"]) == pytest.approx(0.1457492)


def test_manifestte_olmayan_xvg_yok_sayilir(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {
        "rmsd_pep_on_mhc.xvg": TIMESERIES_XVG,
        "elle_yazilmis_baska.xvg": TIMESERIES_XVG,
    })
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    assert {row["output"] for row in ts} == {"rmsd_pep_on_mhc.xvg"}


def test_cok_serili_xvg_her_seriyi_ayri_satira_yazar(mdkit, fake_config, fake_dataset, tmp_path):
    write_outputs(fake_dataset, {"rmsd_pep_on_mhc.xvg": MULTI_SERIES_XVG})
    run_collect(mdkit, fake_config)
    ts = list(csv.DictReader((tmp_path / "results" / "timeseries_long.csv").open()))
    series = {row["series"] for row in ts}
    assert series == {"birinci", "ikinci"}
    assert len(ts) == 2 * 3 * 2 * 2   # kompleks x replika x satir x seri


def test_bos_veri_seti_sadece_baslik_yazar(mdkit, fake_config, tmp_path):
    r = run_collect(mdkit, fake_config)
    assert r.returncode == 0, r.stderr
    lines = (tmp_path / "results" / "timeseries_long.csv").read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("complex,replica,analysis,output,series,time_ps")
```

- [ ] **Step 2: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_collect.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'collect_results'`.

- [ ] **Step 3: `collect_results.py`'yi yaz**

`mdkit/collect_results.py`:

```python
#!/usr/bin/env python
"""mdkit: rep*/analysis/*.xvg -> results/{timeseries,profile}_long.csv

Analiz manifestosu `run_analysis.sh --list` ciktisindan okunur; boylece hangi
ciktinin hangi analize ve hangi KIND'a ait oldugu tek kaynakta (analysis/*.sh)
kalir ve burada kopyalanmaz.

Degerler .xvg'deki ham haliyle yazilir (zaman ps, mesafe nm) ve birim `unit`
kolonunda belirtilir. Donusum cizim katmaninda yapilir: korlemesine nm->A
carpmak ileride gmx sasa (nm^2) eklendiginde sessizce yanlis olurdu.
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

MDKIT = Path(__file__).resolve().parent

XVG_LEGEND = re.compile(r'@\s+s(\d+)\s+legend\s+"(.*)"')
XVG_LABEL = re.compile(r'@\s+(title|subtitle)\s+"(.*)"')
XVG_AXIS = re.compile(r'@\s+(xaxis|yaxis)\s+label\s+"(.*)"')
UNIT_IN_LABEL = re.compile(r"\(([^)]*)\)")

TS_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "time_ps", "value", "unit"]
PR_FIELDS = ["complex", "replica", "analysis", "output", "series",
             "residue", "value", "unit"]


def parse_xvg(path):
    """xvg dosyasini (meta, satirlar) olarak dondurur."""
    meta = {"legends": {}}
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            m = XVG_LEGEND.match(line)
            if m:
                meta["legends"][int(m.group(1))] = m.group(2)
                continue
            m = XVG_LABEL.match(line)
            if m:
                meta[m.group(1)] = m.group(2)
                continue
            m = XVG_AXIS.match(line)
            if m:
                meta[m.group(1)] = m.group(2)
            continue
        try:
            rows.append([float(p) for p in line.split()])
        except ValueError:
            continue
    return meta, rows


def _bash(snippet):
    r = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr.strip() or "bash cagrisi basarisiz")
    return r.stdout


def read_config(config):
    """Config'i bash tarafindaki dogrulayiciyla okur (tek kaynak)."""
    arg = f'"{config}"' if config else ""
    out = _bash(
        f'source "{MDKIT}/analysis/lib.sh" && mdkit_load_config {arg} && '
        'printf "%s\\n%s\\n%s\\n%s\\n" '
        '"$DATA_ROOT" "$RESULTS_DIR" "$COMPLEX_GLOB" "${REPS[*]}"'
    )
    data_root, results_dir, complex_glob, reps = out.splitlines()[:4]
    return Path(data_root), Path(results_dir), complex_glob, reps.split()


def read_manifest(config):
    """cikti_dosya_adi -> (analiz_adi, KIND)"""
    cmd = ["bash", str(MDKIT / "run_analysis.sh"), "--list"]
    if config:
        cmd += ["--config", str(config)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr.strip() or "analiz manifestosu okunamadi")
    manifest = {}
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        name, kind, _desc, outputs = line.split("\t")
        for out in outputs.split(","):
            manifest[out.strip()] = (name, kind)
    return manifest


def collect(data_root, complex_glob, reps, manifest):
    timeseries, profile = [], []
    for cx in sorted(data_root.glob(complex_glob)):
        if not cx.is_dir():
            continue
        cname = cx.name.split("_")[0]
        for rep in reps:
            adir = cx / rep / "analysis"
            if not adir.is_dir():
                continue
            for xvg in sorted(adir.glob("*.xvg")):
                entry = manifest.get(xvg.name)
                if entry is None:
                    continue
                analysis, kind = entry
                meta, rows = parse_xvg(xvg)
                m = UNIT_IN_LABEL.search(meta.get("yaxis", ""))
                unit = m.group(1) if m else ""
                for row in rows:
                    x, ys = row[0], row[1:]
                    for i, y in enumerate(ys):
                        series = (meta["legends"].get(i)
                                  or meta.get("subtitle")
                                  or f"y{i}")
                        rec = {
                            "complex": cname, "replica": rep,
                            "analysis": analysis, "output": xvg.name,
                            "series": series, "value": y, "unit": unit,
                        }
                        if kind == "timeseries":
                            rec["time_ps"] = x
                            timeseries.append(rec)
                        elif kind == "profile":
                            rec["residue"] = int(x)
                            profile.append(rec)
    return timeseries, profile


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description="mdkit .xvg toplayici")
    ap.add_argument("-c", "--config", type=Path, default=None,
                    help="config.sh (varsayilan: mdkit/config.sh)")
    ap.add_argument("-o", "--out-dir", type=Path, default=None,
                    help="cikti dizini (varsayilan: config'deki RESULTS_DIR)")
    args = ap.parse_args()

    data_root, results_dir, complex_glob, reps = read_config(args.config)
    out_dir = args.out_dir or results_dir
    manifest = read_manifest(args.config)

    timeseries, profile = collect(data_root, complex_glob, reps, manifest)
    write_csv(out_dir / "timeseries_long.csv", TS_FIELDS, timeseries)
    write_csv(out_dir / "profile_long.csv", PR_FIELDS, profile)

    print(f"timeseries: {len(timeseries)} satir -> {out_dir / 'timeseries_long.csv'}")
    print(f"profile   : {len(profile)} satir -> {out_dir / 'profile_long.csv'}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_collect.py -v`
Expected: 8 passed.

- [ ] **Step 5: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add xvg collector producing long-format result tables"
```

---

### Task 10: `plot_results.py` — üç çizim modu

**Files:**
- Create: `mdkit/plot_results.py`
- Test: `mdkit/tests/test_plot.py`

**Interfaces:**
- Consumes: `results/timeseries_long.csv`, `results/profile_long.csv` (Task 9).
- Produces:
  - `compare_order(per_rep_df) -> list[str]` — kompleksleri medyan değere göre artan sırada verir (saf fonksiyon, render'dan bağımsız test edilir).
  - `results/plots/per_complex/<kompleks>_<cikti>.png`
  - `results/plots/mean_sd/<kompleks>_<cikti>.png`
  - `results/plots/compare_<cikti>.png`

- [ ] **Step 1: `dataviz` skill'ini yükle**

Grafik kodunun ilk satırını yazmadan önce `dataviz` skill'ini yükle ve palet/eksen/legend kurallarını oradan al. Aşağıdaki kodda `REP_COLORS` ve `GROUP_COLORS` somut değerlerle verilmiştir; skill'in paleti farklıysa **yalnızca bu iki sabiti** değiştir, kodun yapısı aynı kalır. 35 kategorili `--compare` panelinde okunabilirliği belirleyen şey bu seçimdir.

- [ ] **Step 2: Başarısız testi yaz**

`mdkit/tests/test_plot.py`:

```python
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_results  # noqa: E402

pd = pytest.importorskip("pandas")

TS_HEADER = "complex,replica,analysis,output,series,time_ps,value,unit\n"
PR_HEADER = "complex,replica,analysis,output,series,residue,value,unit\n"


@pytest.fixture
def results_dir(tmp_path):
    d = tmp_path / "results"
    d.mkdir()
    ts_rows, pr_rows = [], []
    for cx, base in [("last1", 0.20), ("top1", 0.10), ("last2", 0.30)]:
        for i, rep in enumerate(["rep1", "rep2", "rep3"]):
            for t in range(0, 50, 10):
                ts_rows.append(
                    f"{cx},{rep},rmsd,rmsd_pep_on_mhc.xvg,LIGAND,{t}.0,"
                    f"{base + 0.01 * i + 0.001 * t},nm"
                )
            for res in range(1, 6):
                pr_rows.append(
                    f"{cx},{rep},rmsf,rmsf_pep_self.xvg,LIGAND,{res},"
                    f"{base + 0.02 * res},nm"
                )
    (d / "timeseries_long.csv").write_text(TS_HEADER + "\n".join(ts_rows) + "\n")
    (d / "profile_long.csv").write_text(PR_HEADER + "\n".join(pr_rows) + "\n")
    return d


def run_plot(mdkit, results_dir, *args):
    return subprocess.run(
        [sys.executable, str(mdkit / "plot_results.py"),
         "--results-dir", str(results_dir), *args],
        capture_output=True, text=True,
    )


def test_compare_order_medyana_gore_siralar():
    df = pd.DataFrame({
        "complex": ["a", "a", "b", "b", "c", "c"],
        "value": [3.0, 3.2, 1.0, 1.1, 2.0, 2.2],
    })
    assert plot_results.compare_order(df) == ["b", "c", "a"]


def test_per_complex_kompleks_basina_figure_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--per-complex")
    assert r.returncode == 0, r.stderr
    out = results_dir / "plots" / "per_complex"
    names = sorted(p.name for p in out.glob("*.png"))
    assert "last1_rmsd_pep_on_mhc.png" in names
    assert "last1_rmsf_pep_self.png" in names
    assert len(names) == 3 * 2          # 3 kompleks x 2 cikti
    assert all((out / n).stat().st_size > 1000 for n in names)


def test_mean_sd_figure_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--mean-sd")
    assert r.returncode == 0, r.stderr
    out = results_dir / "plots" / "mean_sd"
    assert len(list(out.glob("*.png"))) == 3 * 2


def test_compare_tek_panel_uretir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir, "--compare")
    assert r.returncode == 0, r.stderr
    p = results_dir / "plots" / "compare_rmsd_pep_on_mhc.png"
    assert p.exists() and p.stat().st_size > 1000


def test_mod_verilmezse_hepsi_calisir(mdkit, results_dir):
    r = run_plot(mdkit, results_dir)
    assert r.returncode == 0, r.stderr
    assert (results_dir / "plots" / "per_complex").is_dir()
    assert (results_dir / "plots" / "mean_sd").is_dir()
    assert (results_dir / "plots" / "compare_rmsd_pep_on_mhc.png").exists()


def test_eksik_csv_anlasilir_hata(mdkit, tmp_path):
    r = run_plot(mdkit, tmp_path / "yok")
    assert r.returncode != 0
    assert "bulunamadi" in (r.stderr + r.stdout)
```

- [ ] **Step 3: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_plot.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'plot_results'`.

- [ ] **Step 4: `plot_results.py`'yi yaz**

`mdkit/plot_results.py`:

```python
#!/usr/bin/env python
"""mdkit cizim katmani.

Cizim ANALIZ ADINA gore degil, veri SEKLINE gore dallanir: timeseries ve
profile. Yeni bir timeseries analizi (gyrate, sasa) eklendiginde bu dosyaya
dokunmak gerekmez.

Birim donusumu burada yapilir: .xvg'ler ps ve nm cinsindendir (spec 2.5).
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

NM_TO_ANGSTROM = 10.0
PS_TO_NS = 1e-3

# dataviz skill'inin paleti. Sadece bu iki sabit degistirilir.
REP_COLORS = {"rep1": "#4C72B0", "rep2": "#DD8452", "rep3": "#55A868"}
GROUP_COLORS = {"top": "#4C72B0", "last": "#C44E52"}

DEFAULT_COMPARE_OUTPUT = "rmsd_pep_on_mhc.xvg"


def load(results_dir):
    ts_path = results_dir / "timeseries_long.csv"
    pr_path = results_dir / "profile_long.csv"
    if not ts_path.exists() and not pr_path.exists():
        sys.exit(f"sonuc CSV'leri bulunamadi: {results_dir}")
    ts = pd.read_csv(ts_path) if ts_path.exists() else pd.DataFrame()
    pr = pd.read_csv(pr_path) if pr_path.exists() else pd.DataFrame()
    return ts, pr


def _group_color(complex_name):
    return GROUP_COLORS["top" if complex_name.startswith("top") else "last"]


def _panels(ts, pr):
    """(df, x kolonu, x etiketi, x carpani) uclusu uretir."""
    if not ts.empty:
        yield ts, "time_ps", "Zaman (ns)", PS_TO_NS
    if not pr.empty:
        yield pr, "residue", "Residue", 1.0


def plot_per_complex(ts, pr, out_dir):
    """Kompleks basina bir figure; rep1/rep2/rep3 ust uste. Yakinsama denetimi."""
    target = out_dir / "per_complex"
    target.mkdir(parents=True, exist_ok=True)
    for df, xcol, xlabel, xconv in _panels(ts, pr):
        for (cx, output), g in df.groupby(["complex", "output"], sort=True):
            fig, ax = plt.subplots(figsize=(7, 4))
            for rep, gr in g.groupby("replica", sort=True):
                gr = gr.sort_values(xcol)
                ax.plot(gr[xcol] * xconv, gr["value"] * NM_TO_ANGSTROM,
                        label=rep, color=REP_COLORS.get(rep), linewidth=1.0)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Å")
            ax.set_title(f"{cx} — {Path(output).stem}")
            ax.legend(frameon=False)
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            fig.savefig(target / f"{cx}_{Path(output).stem}.png", dpi=150)
            plt.close(fig)


def plot_mean_sd(ts, pr, out_dir):
    """Replika ortalamasi + ±SD seridi. Yayina/teze giden temiz figure."""
    target = out_dir / "mean_sd"
    target.mkdir(parents=True, exist_ok=True)
    for df, xcol, xlabel, xconv in _panels(ts, pr):
        for (cx, output), g in df.groupby(["complex", "output"], sort=True):
            stats = (g.groupby(xcol)["value"]
                       .agg(["mean", "std"])
                       .sort_index()
                       .fillna(0.0))
            x = stats.index.to_numpy() * xconv
            mean = stats["mean"].to_numpy() * NM_TO_ANGSTROM
            sd = stats["std"].to_numpy() * NM_TO_ANGSTROM
            color = _group_color(cx)

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.fill_between(x, mean - sd, mean + sd, color=color, alpha=0.25,
                            linewidth=0)
            ax.plot(x, mean, color=color, linewidth=1.4)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Å")
            ax.set_title(f"{cx} — {Path(output).stem} (n={g['replica'].nunique()} replika)")
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            fig.savefig(target / f"{cx}_{Path(output).stem}.png", dpi=150)
            plt.close(fig)


def compare_order(per_rep):
    """Kompleksleri medyan degere gore artan sirada dondurur (saf fonksiyon)."""
    return (per_rep.groupby("complex")["value"]
                   .median()
                   .sort_values()
                   .index
                   .tolist())


def plot_compare(ts, out_dir, output=DEFAULT_COMPARE_OUTPUT):
    """Tum kompleksler tek panelde, medyana gore sirali boxplot."""
    if ts.empty:
        return
    df = ts[ts["output"] == output]
    if df.empty:
        return
    per_rep = (df.groupby(["complex", "replica"])["value"].mean()
                 .mul(NM_TO_ANGSTROM)
                 .reset_index())
    order = compare_order(per_rep)
    data = [per_rep.loc[per_rep["complex"] == c, "value"].to_numpy() for c in order]

    fig, ax = plt.subplots(figsize=(max(8.0, len(order) * 0.35), 4.5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.6)
    for patch, cx in zip(bp["boxes"], order):
        patch.set_facecolor(_group_color(cx))
        patch.set_alpha(0.75)
        patch.set_edgecolor("#333333")
    for median in bp["medians"]:
        median.set_color("#222222")

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylabel(f"Ortalama {Path(output).stem} (Å)")
    ax.set_title("Kompleksler arasi karsilastirma (replika basina ortalama)")
    ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Line2D([], [], color=c, linewidth=6, alpha=0.75)
               for c in (GROUP_COLORS["top"], GROUP_COLORS["last"])]
    ax.legend(handles, ["top*", "last*"], frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / f"compare_{Path(output).stem}.png", dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="mdkit cizim katmani")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--per-complex", action="store_true")
    ap.add_argument("--mean-sd", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--compare-output", default=DEFAULT_COMPARE_OUTPUT)
    args = ap.parse_args()

    ts, pr = load(args.results_dir)
    out_dir = args.results_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_all = not (args.per_complex or args.mean_sd or args.compare)
    if args.per_complex or run_all:
        plot_per_complex(ts, pr, out_dir)
    if args.mean_sd or run_all:
        plot_mean_sd(ts, pr, out_dir)
    if args.compare or run_all:
        plot_compare(ts, out_dir, args.compare_output)

    print(f"grafikler -> {out_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Testleri koş, geçtiklerini doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_plot.py -v`
Expected: 6 passed.

- [ ] **Step 6: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "feat(mdkit): add plotting layer with per-complex, mean-sd and compare modes"
```

---

### Task 11: README, taşınabilirlik testi ve bilinen sonuçla çapraz doğrulama

**Files:**
- Create: `mdkit/README.md`
- Modify: `mdkit/tests/conftest.py` (gerçek veri yollarını ortam değişkeninden okunabilir yap; `slow` işaretini tanıt)
- Test: `mdkit/tests/test_portability.py`, `mdkit/tests/test_cross_validation.py`

**Depo kısıtı:** `.gitignore` ters-whitelist'tir — yalnızca `*.sh`, `*.py`, `*.md`, `*.csv` (ve birkaç format daha) izlenir. Bu plan boyunca **başka uzantıda dosya oluşturulmaz**.

**Interfaces:**
- Consumes: tüm önceki task'lar.
- Produces: `MDKIT_TEST_DATA_ROOT` ve `MDKIT_TEST_GMX` ortam değişkenleri — laboratuvardaki başka bir makinede testlerin koşabilmesi için.

- [ ] **Step 1: `conftest.py`'deki gerçek veri yollarını ortam değişkenine bağla**

`conftest.py`'deki

```python
REAL_ROOT = Path("/mnt/data/scratch-simulations-TUSEB")
REAL_REP = REAL_ROOT / "last10_IMGQQPAPQV_A0201_pandora" / "rep1"
GMX_BIN = Path("/usr/local/gromacs-2025.4-cuda/bin/gmx")
```

satırlarını şununla değiştir:

```python
# Testler bu makinenin gercek verisine bakar; baska bir makinede kosmak icin
# asagidaki ortam degiskenleri ayarlanir. (mdkit'in KENDISI hicbir mutlak yol
# icermez -- bkz. test_portability.py; bu istisna yalnizca test verisi icindir.)
REAL_ROOT = Path(os.environ.get("MDKIT_TEST_DATA_ROOT",
                                "/mnt/data/scratch-simulations-TUSEB"))
REAL_REP_NAME = os.environ.get("MDKIT_TEST_REP",
                               "last10_IMGQQPAPQV_A0201_pandora/rep1")
REAL_REP = REAL_ROOT / REAL_REP_NAME
GMX_BIN = Path(os.environ.get("MDKIT_TEST_GMX",
                              "/usr/local/gromacs-2025.4-cuda/bin/gmx"))
```

- [ ] **Step 2: Taşınabilirlik testini yaz**

`mdkit/tests/test_portability.py`:

```python
import subprocess


def test_calisma_kodunda_mutlak_yol_yok(mdkit):
    """config.sh disinda hicbir calisma dosyasinda makineye ozel yol olmamali.
    tests/ haric tutulur: testler bu makinenin gercek verisine bakar ve
    yollari MDKIT_TEST_* ortam degiskenleriyle ezilebilir (bkz. conftest.py)."""
    r = subprocess.run(
        ["grep", "-rEn", "/mnt/|/usr/local/|/home/", str(mdkit),
         "--include=*.sh", "--include=*.py",
         "--exclude=config.sh", "--exclude-dir=tests",
         "--exclude-dir=__pycache__"],
        capture_output=True, text=True,
    )
    assert r.stdout == "", f"mutlak yol bulundu:\n{r.stdout}"


def test_proje_adi_kodda_gecmiyor(mdkit):
    """'pandora', 'TUSEB', 'A0201' gibi projeye ozel isimler config'de kalmali."""
    r = subprocess.run(
        ["grep", "-rEni", "pandora|TUSEB|A0201", str(mdkit),
         "--include=*.sh", "--include=*.py",
         "--exclude=config.sh", "--exclude-dir=tests",
         "--exclude-dir=__pycache__"],
        capture_output=True, text=True,
    )
    assert r.stdout == "", f"projeye ozel isim bulundu:\n{r.stdout}"


def test_kopya_configle_calisir(mdkit, fake_config, fake_dataset):
    """Aracin config disinda hicbir seye bagli olmadiginin uctan uca kaniti."""
    r = subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(fake_config),
         "--dry-run", "--all", str(fake_dataset)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "DRY-RUN" in r.stdout
```

- [ ] **Step 3: Testleri koş, başarısız olduklarını doğrula**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_portability.py -v`
Expected: `test_calisma_kodunda_mutlak_yol_yok` FAIL olabilir — `config.sh` dışında sızmış yol varsa. Varsa **kodu düzelt**, testi gevşetme.

- [ ] **Step 4: Bilinen sonuçla çapraz doğrulama testini yaz**

`mdkit/tests/test_cross_validation.py`:

```python
import csv
import subprocess

import pytest

from conftest import GMX_BIN, REAL_REP, REAL_ROOT, needs_gmx

KNOWN_CSV = REAL_ROOT / "rmsf_per_position.csv"

pytestmark = pytest.mark.slow


def _parse_xvg_values(path):
    return [float(ln.split()[1]) for ln in path.read_text().splitlines()
            if ln.strip() and not ln.startswith(("#", "@"))]


@needs_gmx
@pytest.mark.skipif(not KNOWN_CSV.exists(), reason="rmsf_per_position.csv yok")
def test_gmx_cagri_katmani_bilinen_ciktiyi_yeniden_uretir(tmp_path):
    """rmsf_analysis.py'nin urettigi rmsf_per_position.csv, 'Protein' grubuna
    fit edilmis RMSF'dir (-b 10000). Ayni komutu tekrarlayip son 10 residue'yu
    karsilastirmak, gmx cagri katmanini bilinen-dogru bir ciktiya karsi
    dogrular. TAM trajektori gerektirir, birkac dakika surer."""
    out = tmp_path / "rmsf_protein.xvg"
    subprocess.run(
        [str(GMX_BIN), "rmsf",
         "-s", str(REAL_REP / "md_0_10.tpr"),
         "-f", str(REAL_REP / "traj_compact_center_dry.xtc"),
         "-o", str(out), "-res", "-b", "10000"],
        input="1\n", text=True, capture_output=True, check=True,
    )
    computed = [v * 10.0 for v in _parse_xvg_values(out)][-10:]

    known = [
        float(row["rmsf_angstrom"])
        for row in csv.DictReader(KNOWN_CSV.open())
        if row["complex"] == "last10" and row["replica"] == "rep1"
    ]
    assert len(known) == 10
    for got, expected in zip(computed, known):
        assert abs(got - expected) < 0.002, (computed, known)


@needs_gmx
def test_mdkit_rmsf_farkli_bir_buyukluk_olcer(mdkit, real_config, tmp_path):
    """mdkit'in rmsf_pep_self'i LIGAND grubuna fit edilir; rmsf_analysis.py
    ise Protein grubuna fit eder. Bunlar TANIM GEREGI farkli sayilardir ve
    birbirinin yerine kullanilamaz -- bu testi gecmesi o ayrimi kilitler."""
    target = tmp_path / "data" / "test1_PEPTIDE_A0201_pandora"
    subprocess.run(
        ["bash", str(mdkit / "run_analysis.sh"), "-c", str(real_config),
         "-y", "-a", "rmsf", "-b", "0", str(target)],
        capture_output=True, text=True, check=False,
    )
    out = target / "rep1" / "analysis"
    pep_self = _parse_xvg_values(out / "rmsf_pep_self.xvg")

    out_protein = tmp_path / "rmsf_protein_short.xvg"
    subprocess.run(
        [str(GMX_BIN), "rmsf",
         "-s", str(target / "rep1" / "md_0_10.tpr"),
         "-f", str(target / "rep1" / "traj_compact_center_dry.xtc"),
         "-o", str(out_protein), "-res", "-b", "0"],
        input="1\n", text=True, capture_output=True, check=True,
    )
    protein_fit_peptide = _parse_xvg_values(out_protein)[-10:]

    assert len(pep_self) == 10
    assert pep_self != pytest.approx(protein_fit_peptide, abs=1e-6)
```

- [ ] **Step 5: `slow` işaretini `conftest.py` içinde tanıt**

Bu depo `.gitignore`'u ters-whitelist'tir (`*` + `!*.sh`, `!*.py`, `!*.md`); bir
`pytest.ini` dosyası **git tarafından takip edilmez**. İşareti bu yüzden
`conftest.py` içinde tanıtıyoruz — hem izlenir hem de ek dosya gerekmez.

`mdkit/tests/conftest.py` sonuna:

```python
def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: tam trajektori gerektiren, dakikalar suren testler"
    )
```

- [ ] **Step 6: Testleri koş**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest test_cross_validation.py -v -m slow`
Expected: 2 passed (gerçek veri yoksa skipped).

- [ ] **Step 7: `README.md`'yi yaz**

`mdkit/README.md`:

```markdown
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
```

- [ ] **Step 8: Tüm test paketini koş**

Run: `cd mdsimulations/postmd_analysis/mdkit/tests && /home/emre/anaconda3/bin/python -m pytest -v`
Expected: tüm testler passed veya skipped; hiç failed yok.

- [ ] **Step 9: Uçtan uca gerçek koşu (doğrulama planı §8.1-8.2, §8.4)**

```bash
cd mdsimulations/postmd_analysis/mdkit
./run_analysis.sh --all --dry-run /mnt/data/scratch-simulations-TUSEB | tail -20
./run_analysis.sh -y -a rmsd,rmsf /mnt/data/scratch-simulations-TUSEB/last10_IMGQQPAPQV_A0201_pandora
./run_analysis.sh -y -a rmsd,rmsf /mnt/data/scratch-simulations-TUSEB/last10_IMGQQPAPQV_A0201_pandora
```

Doğrula:
- İlk `--dry-run` **180** `DRY-RUN` satırı üretir. (Naif çarpım 35 × 3 × 2 = 210 verir, ama `check_ref.pdb`'si eksik 30 replika analiz başına bir satır yerine tek bir `NO_REF` satırı yazar: 75 sağlam replika × 2 analiz = 150, artı 30 `NO_REF` = 180.)
- İlk gerçek koşuda `rmsd_*.xvg` satır sayısı 10001 (0–100 ns, dt 10 ps), `rmsf_pep_self.xvg` satır sayısı 10.
- İkinci koşuda tüm satırlar `SKIP_DONE`, süre ~0.

- [ ] **Step 10: Commit** *(kullanıcı onayı gerekli)*

```bash
git add mdsimulations/postmd_analysis/mdkit/
git commit -m "docs(mdkit): add README, portability test and cross-validation against known output"
```

---

## Plan Öz-İncelemesi

**1. Spec kapsamı** — spec'teki her bölüm bir task'a bağlandı:

| Spec | Task |
|---|---|
| §3.1 dosya yapısı, mutlak yol yasağı | 1, 11 |
| §3.2 çıktı yapısı | 4, 9 |
| §3.3 `config.sh` sınırı | 1, 6 (POSTMD_SCRIPT) |
| §3.4 eklenti sözleşmesi | 5 (metadata), 7-8 (uygulama) |
| §4.0-4.1 config + gmx keşfi | 1, 3 |
| §4.2-4.3 hedef keşfi, ön koşul | 2 |
| §4.4 `check_ref` self-heal | 3 (fonksiyon), 6 (onay akışı) |
| §4.5 index kurulumu | 3 |
| §4.6 idempotency | 4, 6 |
| §4.7 hata yalıtımı | 4, 6 |
| §4.8 log + durum sözlüğü | 4 |
| §5 CLI | 5, 6 |
| §6.1 rmsd | 7 |
| §6.2 rmsf + groove-fit | 8 |
| §6.3 `rmsf_analysis.py` ile ilişki | 11 (çapraz doğrulama, ayrım kilitlenir) |
| §7.1 collect | 9 |
| §7.2 plot (3 mod) | 10 |
| §8.1-8.6 doğrulama | 2-8 testleri + 11 Step 9 |
| §8.7 taşınabilirlik | 11 |
| §10 çıkarma planı | — (kod değil; `mdkit/` sınırı Task 1'de kuruldu) |

**2. Placeholder taraması** — Task 5 Step 4'teki iskelet `analysis_run` gövdeleri bilinçli geçicidir ve Task 7-8'de değiştirilecekleri hem adımda hem kodda yazılıdır; bunun dışında "TBD"/"sonra doldur" yok. Her kod adımı tam içerik veriyor.

**3. Tip/isim tutarlılığı** — fonksiyon adları task'lar arası kontrol edildi: `mdkit_load_config`, `mdkit_find_complexes`, `mdkit_rep_status`, `mdkit_complex_name`, `mdkit_resolve_gmx`, `mdkit_make_ref`, `mdkit_build_index`, `mdkit_group_size`, `mdkit_log_init`, `mdkit_log_row`, `mdkit_outputs_present`, `mdkit_run_isolated`, `mdkit_analysis_scripts`, `mdkit_analysis_meta`. Çıktı dosya adları Task 5 iskeletinde ilan edilenle Task 7-8'de üretilenle aynı (`rmsd_pep_on_mhc.xvg`, `rmsd_pep_internal.xvg`, `rmsd_complex_bb.xvg`, `rmsf_pep_self.xvg`, `rmsf_mhc.xvg`, `rmsf_pep_groovefit.xvg`) ve Task 9-10 testleri de bu adları kullanıyor. CSV kolon adları Task 9'da tanımlanıp Task 10 testlerinde birebir kullanılıyor.
