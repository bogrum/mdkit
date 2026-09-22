import subprocess
import sys
from pathlib import Path

import numpy as np
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


# --- Fix 1: unit-aware scaling -------------------------------------------

def test_scale_and_label_nm_angstroma_cevirir():
    assert plot_results.scale_and_label("nm") == (10.0, "Å")


def test_scale_and_label_bilinmeyen_birimi_cevirmeden_birakir():
    assert plot_results.scale_and_label("nm^2") == (1.0, "nm^2")


def test_scale_and_label_bos_birim_birimsiz_doner():
    assert plot_results.scale_and_label("") == (1.0, "birimsiz")


@pytest.fixture
def results_dir_odd_unit(tmp_path):
    """timeseries'i nm disinda bir birimle (nm^2, sasa benzeri) tasiyan set."""
    d = tmp_path / "results"
    d.mkdir()
    rows = []
    for rep in ["rep1", "rep2"]:
        for t in range(0, 30, 10):
            rows.append(
                f"only1,{rep},sasa,sasa_total.xvg,LIGAND,{t}.0,"
                f"{1.5 + 0.01 * t},nm^2"
            )
    (d / "timeseries_long.csv").write_text(TS_HEADER + "\n".join(rows) + "\n")
    (d / "profile_long.csv").write_text(PR_HEADER)
    return d


def test_bilinmeyen_birim_cevrilmeden_cizilir_ve_stderrde_ciktiyi_adlandirir(
    mdkit, results_dir_odd_unit
):
    r = run_plot(mdkit, results_dir_odd_unit, "--per-complex")
    assert r.returncode == 0, r.stderr
    p = results_dir_odd_unit / "plots" / "per_complex" / "only1_sasa_total.png"
    assert p.exists() and p.stat().st_size > 1000
    assert "sasa_total.xvg" in r.stderr


# --- Fix 2: error isolation in the render loops --------------------------

def _rows(cx, output, unit="nm", n=3, base=0.1):
    return [
        {
            "complex": cx, "replica": "rep1", "analysis": "rmsd",
            "output": output, "series": "LIGAND",
            "time_ps": float(t * 10), "value": base + 0.001 * t, "unit": unit,
        }
        for t in range(n)
    ]


def test_plot_per_complex_bir_kompleksin_hatasi_digerlerini_durdurmaz(tmp_path):
    output = "rmsd_pep_on_mhc.xvg"
    rows = (_rows("good1", output, base=0.10)
            + _rows("bad/name", output, base=0.15)
            + _rows("good2", output, base=0.20))
    ts = pd.DataFrame(rows)
    pr = pd.DataFrame(columns=[
        "complex", "replica", "analysis", "output", "series",
        "residue", "value", "unit",
    ])
    out_dir = tmp_path / "plots"

    plot_results.plot_per_complex(ts, pr, out_dir)  # not raise

    target = out_dir / "per_complex"
    assert (target / f"good1_{Path(output).stem}.png").exists()
    assert (target / f"good2_{Path(output).stem}.png").exists()


def test_plot_mean_sd_tek_replika_std_nan_ile_bile_calisir(tmp_path):
    output = "rmsd_pep_on_mhc.xvg"
    ts = pd.DataFrame(_rows("solo", output, base=0.15))
    pr = pd.DataFrame(columns=[
        "complex", "replica", "analysis", "output", "series",
        "residue", "value", "unit",
    ])
    out_dir = tmp_path / "plots"

    plot_results.plot_mean_sd(ts, pr, out_dir)  # not raise

    p = out_dir / "mean_sd" / f"solo_{Path(output).stem}.png"
    assert p.exists() and p.stat().st_size > 1000


# --- Fix 3: cok serili ciktilar serileri karistirmamali -------------------

MULTI_SERIES_OUTPUT = "gyrate.xvg"


def _multi_series_rows(cx="only1", reps=("rep1", "rep2", "rep3")):
    """gyrate benzeri iki seri: Rg ~2.0, RgX ~1.0.

    Seriler arasi fark (1.0) replika ici yayilimdan (<=0.02) cok daha
    buyuk secildi; havuzlanmis bir SD bunu sakli tutamaz."""
    rows = []
    for i, rep in enumerate(reps):
        for t in range(0, 40, 10):
            rows.append({
                "complex": cx, "replica": rep, "analysis": "gyrate",
                "output": MULTI_SERIES_OUTPUT, "series": "Rg",
                "time_ps": float(t), "value": 2.0 + 0.01 * i, "unit": "nm",
            })
            rows.append({
                "complex": cx, "replica": rep, "analysis": "gyrate",
                "output": MULTI_SERIES_OUTPUT, "series": "RgX",
                "time_ps": float(t), "value": 1.0 + 0.01 * i, "unit": "nm",
            })
    return rows


def test_line_specs_her_replika_seri_cifti_icin_bir_cizgi_uretir():
    """Seri yok sayilirsa Rg ve RgX ayni replikanin TEK bir cizgisinde
    birlestirilir ve zikzak bir egri cizilir."""
    g = pd.DataFrame(_multi_series_rows())
    specs = plot_results.line_specs(g, "time_ps")

    assert len(specs) == 3 * 2
    assert {(s["replica"], s["series"]) for s in specs} == {
        (r, s) for r in ("rep1", "rep2", "rep3") for s in ("Rg", "RgX")
    }
    # Renk replikayi, cizgi tipi seriyi kodlar.
    for s in specs:
        assert s["color"] == plot_results.REP_COLORS[s["replica"]]
    by_series = {s["series"]: s["linestyle"] for s in specs}
    assert by_series["Rg"] != by_series["RgX"]
    # Hicbir cizgi iki seriyi karistirmaz: her cizgi tek bir sabit degerdedir.
    for s in specs:
        assert len(set(s["y"].round(6))) == 1, s


def test_line_specs_tek_seride_gorunumu_degistirmez():
    """Tek serili ciktida davranis eskisiyle birebir ayni kalmali:
    duz cizgi, etiket yalnizca replika adi."""
    rows = [r for r in _multi_series_rows() if r["series"] == "Rg"]
    specs = plot_results.line_specs(pd.DataFrame(rows), "time_ps")
    assert len(specs) == 3
    assert [s["label"] for s in specs] == ["rep1", "rep2", "rep3"]
    assert {s["linestyle"] for s in specs} == {"-"}


def test_series_stats_sdyi_seri_icinde_hesaplar():
    """Havuzlanmis SD, FARKLI FIZIKSEL BUYUKLUKLER arasindaki yayilimi
    replika degiskenligi gibi gosteren bir band uretir. Burada seri farki
    1.0 nm, replika ici yayilim <=0.02 nm'dir; havuzlanmis SD ~0.5 olurdu."""
    g = pd.DataFrame(_multi_series_rows())
    bands = dict(plot_results.series_stats(g, "time_ps"))

    assert sorted(bands) == ["Rg", "RgX"]
    assert bands["Rg"]["mean"].max() == pytest.approx(2.01, abs=1e-9)
    assert bands["RgX"]["mean"].max() == pytest.approx(1.01, abs=1e-9)
    for name, stats in bands.items():
        assert stats["std"].max() < 0.02, (name, stats)

    pooled = g.groupby("time_ps")["value"].std().max()
    assert pooled > 0.4, pooled     # eski davranisin urettigi anlamsiz band


def test_cok_serili_veri_ucdan_uca_figure_uretir(mdkit, tmp_path):
    d = tmp_path / "results"
    d.mkdir()
    rows = ["{complex},{replica},{analysis},{output},{series},{time_ps},"
            "{value},{unit}".format(**r) for r in _multi_series_rows()]
    (d / "timeseries_long.csv").write_text(TS_HEADER + "\n".join(rows) + "\n")
    (d / "profile_long.csv").write_text(PR_HEADER)

    r = run_plot(mdkit, d)
    assert r.returncode == 0, r.stderr
    assert r.stderr == "", r.stderr
    for sub in ("per_complex", "mean_sd"):
        p = d / "plots" / sub / "only1_gyrate.png"
        assert p.exists() and p.stat().st_size > 1000, sub


# --- Fix 7: top/last ayrimi projeye ozel, config.sh'den okunur -----------

def test_parse_complex_groups_onek_etiket_cozer():
    assert plot_results.parse_complex_groups("top:top* last:last*") == [
        ("top", "top*"), ("last", "last*"),
    ]
    assert plot_results.parse_complex_groups("") == []


def test_parse_complex_groups_bicimsiz_ogeyi_atlar(capsys):
    assert plot_results.parse_complex_groups("top:top* bozuk") == [("top", "top*")]
    assert "bozuk" in capsys.readouterr().err


def test_group_color_yapilandirilmamisken_tek_renk():
    """Gruplama yoksa HER kompleks ayni rengi alir. Eskiden 'top' ile
    baslamayan her kompleks 'last' kirmizisini aliyordu -- baska bir veri
    setinde tum panel kirmizi, efsane 'top*/last*' yaziyordu."""
    assert (plot_results.group_color("top1", ())
            == plot_results.group_color("last1", ())
            == plot_results.group_color("herhangi1", ()))


def test_group_color_yapilandirilmis_gruplari_ayirir():
    groups = [("top", "top*"), ("last", "last*")]
    top = plot_results.group_color("top1", groups)
    last = plot_results.group_color("last1", groups)
    assert top != last
    assert {top, last} == set(plot_results.GROUP_PALETTE[:2])
    # Hicbir gruba uymayan kompleks bir KATEGORI rengi almaz.
    assert plot_results.group_color("xyz1", groups) not in {top, last}


def test_read_complex_groups_configden_okunur(mdkit, fake_config):
    """fake_config COMPLEX_GROUPS tanimlamaz -> gruplama yok.
    Projenin kendi config.sh'i tanimlar -> iki grup."""
    assert plot_results.read_complex_groups(fake_config) == []
    assert plot_results.read_complex_groups(mdkit / "config.sh") == [
        ("top", "top*"), ("last", "last*"),
    ]


def _capture_axes(monkeypatch):
    seen = []
    real = plot_results.plt.subplots

    def fake(*a, **k):
        fig, ax = real(*a, **k)
        seen.append(ax)
        return fig, ax

    monkeypatch.setattr(plot_results.plt, "subplots", fake)
    return seen


def _compare_frame():
    rows = []
    for cx, base in [("last1", 0.20), ("top1", 0.10)]:
        for rep in ["rep1", "rep2"]:
            rows.append({"complex": cx, "replica": rep, "analysis": "rmsd",
                         "output": "rmsd_pep_on_mhc.xvg", "series": "LIGAND",
                         "time_ps": 0.0, "value": base, "unit": "nm"})
    return pd.DataFrame(rows)


def test_compare_gruplama_yoksa_efsane_cizilmez(tmp_path, monkeypatch):
    seen = _capture_axes(monkeypatch)
    plot_results.plot_compare(_compare_frame(), tmp_path, groups=())
    assert (tmp_path / "compare_rmsd_pep_on_mhc.png").exists()
    assert seen[-1].get_legend() is None


def test_compare_gruplar_verilince_efsane_config_etiketlerini_kullanir(
    tmp_path, monkeypatch
):
    seen = _capture_axes(monkeypatch)
    groups = [("top", "ONDE"), ("last", "ARKADA")]
    plot_results.plot_compare(_compare_frame(), tmp_path, groups=groups)
    legend = seen[-1].get_legend()
    assert legend is not None
    assert [t.get_text() for t in legend.get_texts()] == ["ONDE", "ARKADA"]


def _write_matrix_npz(mdir, cx, rep_i, rep_j, values, unit="nm",
                      analysis="cross_rmsd", title=""):
    mdir.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=np.float32)
    n_y, n_x = values.shape
    np.savez_compressed(
        mdir / f"{cx}_{rep_i}_{analysis}_{rep_j}.npz",
        values=values,
        x_ps=np.arange(n_x, dtype=np.float64) * 100.0,
        y_ps=np.arange(n_y, dtype=np.float64) * 100.0,
        unit=unit, complex=cx, replica_i=rep_i, replica_j=rep_j,
        analysis=analysis, output=f"{analysis}_{rep_j}.xpm", title=title,
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


def test_equilibrium_egrisi_tum_eslere_gore_ortalama(tmp_path):
    """rep_i'nin her frame'i icin, TUM eslerin TUM frame'lerine ortalama
    uzaklik.

    Matris duzeni [y=es frame, x=kendi frame] oldugu icin esler axis=0'da
    yigilir. Bu, birlestirilmis 3Nx3N matriste satir ortalamasi almanin
    replika basina karsiligidir.
    """
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0, 1.0], [1.0, 0.0]])
    _write_matrix_npz(mdir, "last1", "rep1", "rep2", [[2.0, 3.0], [4.0, 5.0]])
    recs = plot_results.load_matrices(tmp_path)[("last1", "cross_rmsd")]
    curves = plot_results.equilibrium_curves(recs)
    assert set(curves) == {"rep1"}
    x_ps, curve = curves["rep1"]
    assert np.allclose(x_ps, [0.0, 100.0])
    # kolon 0: 0,1,2,4 -> 1.75 ; kolon 1: 1,0,3,5 -> 2.25
    assert np.allclose(curve, [1.75, 2.25])


def test_equilibrium_replika_basina_ayri_egri(tmp_path):
    mdir = tmp_path / "matrices"
    for i in ["rep1", "rep2"]:
        for j in ["rep1", "rep2"]:
            _write_matrix_npz(mdir, "last1", i, j, [[0.0, 0.5], [0.5, 0.0]])
    recs = plot_results.load_matrices(tmp_path)[("last1", "cross_rmsd")]
    assert set(plot_results.equilibrium_curves(recs)) == {"rep1", "rep2"}


def test_equilibrium_uyusmayan_frame_sayisi_atlanir(tmp_path, capsys):
    """Ayni rep_i'nin matrisleri ayni kendi-frame sayisina sahip olmali;
    degilse yigmak sessizce yanlis olurdu."""
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0, 1.0], [1.0, 0.0]])
    _write_matrix_npz(mdir, "last1", "rep1", "rep2", [[1.0, 2.0, 3.0]])
    recs = plot_results.load_matrices(tmp_path)[("last1", "cross_rmsd")]
    assert plot_results.equilibrium_curves(recs) == {}
    assert "rep1" in capsys.readouterr().err


def test_equilibrium_figuru_uretilir(tmp_path):
    mdir = tmp_path / "matrices"
    for i in ["rep1", "rep2"]:
        for j in ["rep1", "rep2"]:
            _write_matrix_npz(mdir, "last1", i, j, [[0.0, 0.5], [0.5, 0.0]])
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert (out / "matrix" / "last1_cross_rmsd.png").exists()
    assert (out / "matrix" / "last1_cross_rmsd_equilibrium.png").exists()


def test_equilibrium_hatasi_izgarayi_dusurmez(tmp_path, capsys):
    """Izgara ve equilibrium ayri ayri yalitilir: biri patlarsa digeri
    yine de diske yazilmis olmali."""
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0, 1.0], [1.0, 0.0]])
    _write_matrix_npz(mdir, "last1", "rep1", "rep2", [[1.0, 2.0, 3.0]])
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert (out / "matrix" / "last1_cross_rmsd.png").exists()


def test_rolling_mean_ortalanir(tmp_path):
    """Pencere sonuna baglanirsa egri yarim pencere saga kayar ve
    'ne zaman dengelendi' sorusu sistematik olarak GEC cevaplanir."""
    # Tepe cevresi bilerek asimetrik: duz bir tepe tum pencereleri
    # esitler ve argmax testi anlamsizlasir.
    y = np.array([0.0, 1.0, 9.0, 1.0, 0.0])
    start, vals = plot_results.rolling_mean(y, 3)
    assert start == 1                       # pencere ortasi, sonu (2) degil
    assert np.allclose(vals, [10 / 3, 11 / 3, 10 / 3])
    # Yumusatilmis tepe, girdideki tepeyle AYNI indekse dusmeli.
    assert start + int(np.argmax(vals)) == int(np.argmax(y)) == 2


def test_rolling_mean_kisa_seri_bos_doner(tmp_path):
    start, vals = plot_results.rolling_mean(np.array([1.0, 2.0]), 3)
    assert start == 0 and vals.size == 0


def test_window_label_sureyi_de_yazar(tmp_path):
    """Yalnizca frame sayisi yetmez: ayni 5 frame, -dt 200 ile 1 ns,
    -dt 1000 ile 5 ns eder. Okuyucu pencerenin ne kadar zamana denk
    geldigini bilmeli."""
    x_ps = np.array([0.0, 200.0, 400.0])
    assert plot_results.window_label(5, x_ps) == "pencere=5 frame (1 ns)"


def test_window_label_tek_frame_sureyi_atlar(tmp_path):
    """dt cikarilamiyorsa uydurulmaz."""
    assert plot_results.window_label(3, np.array([0.0])) == "pencere=3 frame"


def test_baslik_npzden_tasinir(tmp_path):
    mdir = tmp_path / "matrices"
    _write_matrix_npz(mdir, "last1", "rep1", "rep1", [[0.0, 0.5], [0.5, 0.0]],
                      title="LIGAND_BB RMSD matrix")
    recs = plot_results.load_matrices(tmp_path)[("last1", "cross_rmsd")]
    assert recs[0]["title"] == "LIGAND_BB RMSD matrix"


def test_baslik_olmayan_eski_npz_cokmez(tmp_path):
    """Bu alan sonradan eklendi; onceden yazilmis .npz'ler onu icermez ve
    yeniden cizim bu yuzden patlamamali."""
    mdir = tmp_path / "matrices"
    mdir.mkdir(parents=True)
    np.savez_compressed(
        mdir / "last1_rep1_cross_rmsd_rep1.npz",
        values=np.array([[0.0, 0.5], [0.5, 0.0]], dtype=np.float32),
        x_ps=np.array([0.0, 100.0]), y_ps=np.array([0.0, 100.0]),
        unit="nm", complex="last1", replica_i="rep1", replica_j="rep1",
        analysis="cross_rmsd",
    )
    recs = plot_results.load_matrices(tmp_path)[("last1", "cross_rmsd")]
    assert recs[0]["title"] == ""
    out = tmp_path / "plots"
    plot_results.plot_matrix(tmp_path, out)
    assert (out / "matrix" / "last1_cross_rmsd.png").exists()
