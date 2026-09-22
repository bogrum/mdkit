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
