"""Fast tests of the stage-22 adjoint variants (widegap.adjoint_variants) on the small 4H-SiC
DMOSFET mesh of the stage-17 tests: the E variant re-assembled through the variant machinery equals
the existing exact gradient bit for bit (R_on,ch, soft-min BV, E_ox), and each variant differs
from E only in its own term (checked against independent closed forms).  < 60 s."""

import numpy as np
import pytest
from scipy.special import expit

from widegap.adjoint_channel import ChannelRon
from widegap.adjoint_sic_oxide import OxideAdjoint
from widegap.adjoint_variants import (
    VariantComposite,
    VariantStats,
    fields_gradients,
    objective_gradients,
    ron_gradients,
)
from widegap.channel_mobility import NoguchiParams
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)

DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=12, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
TH = np.array([0.1, -0.05, 0.03, 0.02, 0.1])


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


@pytest.fixture(scope="module")
def cr(sc):
    return ChannelRon(OxideAdjoint(sc, DESIGN, 0), NoguchiParams.from_registry(sc.reg))


@pytest.fixture(scope="module")
def ev_off(sc):
    return OxideAdjoint(sc, DESIGN, 0, cuts_um=(0.0, 0.8), v_blk_v=100.0)


@pytest.fixture(scope="module")
def ron(cr):
    stats = VariantStats()
    return ron_gradients(cr, TH, stats=stats), stats


def _frac(ev, X):
    mod = ev.mod
    psi, pn, _ = mod.split(X)
    psi_c = mod._psi_c(psi)
    return sum(lv.fraction * expit(-(psi_c - pn + lv.log_gamma)) for lv in mod.ion.donors)


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def test_ron_E_bitwise_and_lagged(cr, ron):
    r, stats = ron
    out = r["out"]
    assert np.array_equal(r["grad"]["E"], out["grad"])
    assert np.array_equal(r["grad"]["A"], out["grad_lagged"])  # A is the existing lagged gradient
    assert stats.n_bitwise_checks >= 3
    assert r["out"]["picard"]["converged"]
    assert r["state_sha"]


def test_ron_variants_differ_only_in_their_term(cr, ron):
    r, _ = ron
    ev, out = cr.ev, r["out"]
    X = out["state"]["X"]
    fac = -out["R"] / out["I_scaled"]
    frac = _frac(ev, X)
    vol = ev.mod.vol
    w, wl = out["w"], out["w_lagged"]
    # the test is only meaningful if incomplete ionization and the mobility term both matter
    assert (1.0 - frac)[w[0::3] != 0].max() > 1e-3
    g = r["g_nd"]
    ion_incomplete = fac * w[0::3] * (-vol * frac)
    # B: Poisson (incomplete-ionization) part only
    assert np.allclose(g["B"], ion_incomplete, rtol=1e-12, atol=0)
    # C: E with the ionized fraction replaced by 1, mobility part unchanged
    assert np.allclose(g["C"] - g["E"], fac * w[0::3] * (-vol) * (1.0 - frac), rtol=1e-9,
                       atol=1e-9 * np.abs(g["E"]).max())  # fmt: skip
    mob = g["E"] - g["B"]  # the bulk-mobility part of E
    assert np.abs(mob).max() > 0
    assert np.allclose(
        g["C"] - fac * w[0::3] * (-vol), mob, rtol=1e-9, atol=1e-9 * np.abs(mob).max()
    )
    # A: same dF/dtheta assembly, multiplier without F_mu mu_X
    assert not np.array_equal(w, wl) and _rel(r["grad"]["A"], r["grad"]["E"]) > 1e-4
    # D: ionization-complete Poisson part of the LAGGED multiplier, no mobility part at all
    assert np.allclose(g["D"], fac * wl[0::3] * (-vol), rtol=1e-12, atol=0)


def test_offstate_E_bitwise_and_C_term(ev_off):
    stats = VariantStats()
    o = objective_gradients(ev_off, TH, stats=stats)
    f = fields_gradients(ev_off, TH, stats=stats)
    assert np.array_equal(o["grad"]["E"], o["res"]["grad"])
    assert np.array_equal(f["E_ox"]["grad"]["E"], f["E_ox"]["E_grad_existing"])
    assert stats.n_bitwise_checks >= 4
    # B: zero-current electron rows are Dirichlet -> mobility part exactly zero -> B == E; D == C
    assert np.array_equal(o["grad"]["B"], o["grad"]["E"])
    assert np.array_equal(f["E_ox"]["grad"]["B"], f["E_ox"]["grad"]["E"])
    assert np.array_equal(o["grad"]["D"], o["grad"]["C"])
    assert np.array_equal(f["E_ox"]["grad"]["D"], f["E_ox"]["grad"]["C"])
    # C - E = -lam-weighted (1 - frac) vol term on the Poisson rows
    mod = ev_off.mod
    X = f["X"]
    frac = _frac(ev_off, X)
    lam = f["E_ox"]["lam"]
    d = f["E_ox"]["g_nd"]
    expect = (-lam)[0::3] * (-mod.vol) * (1.0 - frac)
    assert np.abs(expect).max() > 0
    assert np.allclose(d["C"] - d["E"], expect, rtol=1e-9, atol=1e-9 * np.abs(d["E"]).max())
    assert _rel(f["E_ox"]["grad"]["C"], f["E_ox"]["grad"]["E"]) > 1e-9
    # BV: per-cut C - E, then the soft-min combination
    for k, lam_k in enumerate(o["lams"]):
        st = o["res"]["outs"][k]["state"]
        fr = _frac(ev_off, st.X)
        ex = (-lam_k)[0::3] * (-mod.vol) * (1.0 - fr) * (-1.0 / o["di_dv"][k] * o["sc_v"])
        diff = o["g_nd_cuts"]["C"][k] - o["g_nd_cuts"]["E"][k]
        assert np.allclose(diff, ex, rtol=1e-9, atol=1e-9 * np.abs(o["g_nd_cuts"]["E"][k]).max())


def test_composite_E_equals_existing_and_D_function_values_exact(cr, ev_off):
    comp = VariantComposite(ev_off, cr, "D")
    b = comp.objective(TH)
    b0 = ev_off.objective(TH, grad=False)
    assert b["BV"] == pytest.approx(b0["BV"], rel=1e-12)
    r = comp.ron(TH)
    assert r["R"] == pytest.approx(cr.ron(TH, grad=False)["R"], rel=1e-12)
    assert comp.stats.n_bitwise_checks > 0


def test_kkt_residual_toy():
    from widegap.adjoint_variants import kkt_residual

    # max f = x + y  s.t.  c1 = 1 - x >= 0, c2 = 1 - y >= 0, c3 = 5 - x >= 0 (inactive)
    g = np.array([1.0, 1.0])
    cons = [("c", np.array([0.0, 0.0, 4.0]), np.array([[-1.0, 0.0], [0.0, -1.0], [-1.0, 0.0]]))]
    k = kkt_residual(g, cons)
    assert k["n_active"] == 2 and k["residual"] < 1e-12
    assert np.allclose(k["multipliers"], [1.0, 1.0])
    # not stationary: only c1 active, gradient has a y-component
    cons1 = [("c", np.array([0.0]), np.array([[-1.0, 0.0]]))]
    k1 = kkt_residual(g, cons1)
    assert k1["residual"] == pytest.approx(1.0) and k1["ratio"] == pytest.approx(1 / np.sqrt(2))
    # wrong-sign multiplier is not allowed: constraint pointing the other way leaves f's gradient
    cons2 = [("c", np.array([0.0]), np.array([[1.0, 1.0]]))]
    assert kkt_residual(g, cons2)["residual"] == pytest.approx(np.linalg.norm(g))
    assert kkt_residual(g, [])["ratio"] == 1.0
