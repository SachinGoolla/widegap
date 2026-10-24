"""Stage 26a: CSL layer, retrograde body, gate overhang / contact start, series terms (E1-E3)."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from widegap.adjoint_sic_corner import RoundedSiCDMOS, smooth_step
from widegap.offstate import CorrectedCornerAdjoint
from widegap.series import contact_resistance_sp, source_contact_length_cm, substrate_resistance_sp
from widegap.sic import DMOSDesign, SiC4H, SiCDMOS, TracedRegistry

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")
REPO = Path(__file__).resolve().parents[1]

BASE = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
EXT = DMOSDesign(
    x_jfet=0.5, ch_len=0.3, src_w=2.0, bc_w=0.1, pad_w=0.05, d_src=0.25, d_body=0.6,
    w_drift=2.0, t_drain=0.3, t_ox=0.052, n_drift=1e16, n_jfet=1e17, n_body=5e17, n_src=1e20,
    n_drain=1e19, t_csl=0.3, n_csl=1e17, n_body_b=1e18, l_ov=0.5, src_gap=0.8,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


# ---- E1: defaults are bit-identical to the pre-stage-26 code -----------------------------------
# digests produced with the code of commit eeaa03e (before the extensions) on BASE, level 0
GOLD = {
    "SiCDMOS": (504.04582402362985, -514.0880548599459, 2801725.1449261364, 227639.11282910584,
                32120.79888970869, 2512.1008743461525, 152812.20181244102),
    "RoundedSiCDMOS": (504.04582402362985, -514.0880548599459, 2932619.9423425696,
                       267996.039000813, 30975.478392979057, 4400.891211904948, 804460.4348277581),
}  # fmt: skip


@pytest.mark.parametrize("cls,kw", [(SiCDMOS, {}), (RoundedSiCDMOS, dict(r_um=0.1))])
def test_E1_defaults_bit_identical_small_mesh(sc, cls, kw):
    b = cls(sc, BASE, 0, **kw)
    m = b.dev.mesh
    i = lambda a: float((a * np.arange(len(a))).sum())  # noqa: E731
    eq = b.equilibrium()
    got = (float(m.x.sum()), float(m.y.sum()), i(b.nd), i(b.na), i(b.mu_edge), float(eq.X.sum()),
           i(eq.X))  # fmt: skip
    assert got == pytest.approx(GOLD[cls.__name__], rel=1e-12, abs=0)
    assert b.prm.t_csl == 0.0 and b.prm.a_body_b == 0.0 and b.prm.x_split == ()
    assert not b.dev.m_csl.any()


def test_E1_stage21_initial_on_and_off_state():
    """Stage-21 initial state: BV_obj (h) and channel R_on (h/2) at theta = 0 reproduce
    results/stage21/results.json to <= 1e-12 relative."""
    ref = json.loads((REPO / "results/stage21/results.json").read_text())["initial"]
    spec = importlib.util.spec_from_file_location("stage21_run", REPO / "benchmarks/stage21/run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage21_run"] = mod
    spec.loader.exec_module(mod)
    s17 = mod.S17
    r = s17.make_cr(1).ron(mod.TH0, grad=False)["R"]
    bv = s17.make_ev_off().objective(mod.TH0, grad=False)["BV"]
    assert r == pytest.approx(ref["R_channel_h2"], rel=1e-12)
    assert bv == pytest.approx(ref["bv_obj"], rel=1e-12)


# ---- E2: analytic CSL and retrograde doping ----------------------------------------------------
def _cell_xy(dev):
    cc = dev.mesh.cell_centres
    return cc[:, 0], cc[:, 1]


def test_E2_csl_layer_and_retrograde_cells_match_analytic(sc):
    b = RoundedSiCDMOS(sc, EXT, 0, r_um=0.1)
    dev, p = b.dev, b.prm
    xc, yc = _cell_xy(dev)
    semi = dev.reg.semi_cell
    # CSL cells: exactly n_csl between -(d_body + t_csl) and -d_body
    csl = semi & (yc >= -(p.d_body + p.t_csl)) & (yc < -p.d_body)
    assert csl.any() and np.array_equal(csl, dev.m_csl)
    don = dev.cell_donors(np.zeros(5))
    assert np.max(np.abs((don - dev.src_cell)[csl] / p.n_csl - 1.0)) <= 1e-12  # + n+ tail
    # drift cells below the CSL carry the drift profile with eta = -1 at the CSL bottom
    drift = dev.m_drift
    assert yc[drift].max() < -(p.d_body + p.t_csl)
    assert (
        np.max(np.abs((don - dev.src_cell)[drift] / dev.drift_level(np.zeros(5))[drift] - 1.0))
        <= 1e-12
    )
    assert dev.drift_profile(np.zeros(5), [dev.y_dtop])[0] == pytest.approx(p.n0_drift, rel=1e-14)
    # mesh nodes sit exactly at the layer interfaces
    ys = dev.mesh.y
    for y0 in (-p.d_body, -(p.d_body + p.t_csl)):
        assert np.min(np.abs(ys - y0)) <= 1e-12 * abs(y0)
    # retrograde acceptors: log-linear x erfc rounding, cell by cell
    t = np.clip(-yc / p.d_body, 0.0, 1.0)
    a_y = p.a_body * (p.a_body_b / p.a_body) ** t
    r = b.r_s
    acc = np.where(semi, a_y * smooth_step(xc - p.x_jfet, r) * smooth_step(yc + p.d_body, r), 0.0)
    assert np.max(np.abs(dev.cell_acceptors() - acc)) <= 1e-12 * p.a_body_b
    # surface and bottom values of the profile (in cm^-3)
    n0 = sc.n0
    assert float(dev.body_acceptor(0.0)) * n0 == pytest.approx(5e17, rel=1e-12)
    assert float(dev.body_acceptor(-p.d_body)) * n0 == pytest.approx(1e18, rel=1e-12)
    # nodal doping = semiconductor control-volume average of the cell values (linear map)
    na_nodes = dev.M @ dev.cell_acceptors()
    assert np.max(np.abs(na_nodes - b.na)) <= 1e-12 * p.a_body_b
    nd_nodes = dev.M @ don
    assert np.max(np.abs(nd_nodes - b.nd)) <= 1e-12 * p.b_src


def test_E2_box_device_nodal_values_and_uniform_limit(sc):
    d = EXT.__class__(**{**EXT.__dict__, "n_body_b": 5e17})  # retrograde with N_b = N_s
    d0 = EXT.__class__(**{**EXT.__dict__, "n_body_b": None})
    a = SiCDMOS(sc, d, 0)
    b = SiCDMOS(sc, d0, 0)
    assert np.max(np.abs(a.na - b.na)) <= 1e-12 * a.prm.a_body
    c = SiCDMOS(sc, EXT, 0)  # box device: retrograde body, CSL
    dev, p = c.dev, c.prm
    xs = dev.mesh.node_xy
    deep = (xs[:, 0] > p.x_jfet + 0.0) & (xs[:, 1] < -0.4 * p.d_body) & (xs[:, 1] > -0.6 * p.d_body)
    deep &= xs[:, 0] < p.x_src0
    y = xs[deep, 1]
    # interior nodes: control volume average of a(y) over a mesh column differs from a(y_node)
    # only at second order in the cell size
    ana = p.a_body * (p.a_body_b / p.a_body) ** (-y / p.d_body)
    assert np.max(np.abs(c.na[deep] / ana - 1.0)) < 0.05
    csl_nodes = (xs[:, 1] < -p.d_body - 0.3 * p.t_csl) & (xs[:, 1] > -p.d_body - 0.7 * p.t_csl)
    assert csl_nodes.any()
    assert np.max(np.abs(c.nd[csl_nodes] / p.n_csl - 1.0)) <= 1e-12


# ---- gate overhang, contact start, V_GS ---------------------------------------------------------
def test_gate_overhang_and_contact_start_nodes(sc):
    b = SiCDMOS(sc, EXT, 0)
    p, dev = b.prm, b.dev
    xy = dev.mesh.node_xy
    ls = lambda v: float(sc.length_s(v * 1e-4))  # noqa: E731
    assert p.gate_x1 == pytest.approx(p.x_src0 + ls(0.5), rel=1e-12)
    gate_x = xy[dev.nodes["gate"], 0]
    assert gate_x.max() == pytest.approx(p.gate_x1, rel=1e-12)  # a mesh node sits at the gate end
    src_x = xy[dev.nodes["source"], 0]
    assert src_x.min() == pytest.approx(p.x_src0 + ls(0.8), rel=1e-12)
    assert source_contact_length_cm(b) * 1e4 == pytest.approx(EXT.src_w - 0.8, rel=1e-9)


def test_V_GS_20_is_settable_on_the_ron_path(sc):
    ev = CorrectedCornerAdjoint(sc, BASE, 0, r_um=0.1, cuts_um=(), tau_v=5.0)
    assert ev.V_GS == 15.0
    ev.V_GS = 20.0
    s, bc = ev.on_state()
    assert s.report.converged
    i20 = abs(ev.mod.contact_current(s, ev.dev.nodes["drain"]))
    ev2 = CorrectedCornerAdjoint(sc, BASE, 0, r_um=0.1, cuts_um=(), tau_v=5.0)
    s2, _ = ev2.on_state()
    i15 = abs(ev2.mod.contact_current(s2, ev2.dev.nodes["drain"]))
    assert i20 > i15


# ---- E3: series terms --------------------------------------------------------------------------
def test_E3_series_terms_match_their_formulas(sc):
    rho, x_max, l_c = 1e-3, 5e-4, 2.96e-4
    assert contact_resistance_sp(rho, x_max, l_c) == rho * x_max / l_c
    assert contact_resistance_sp(0.0, x_max, l_c) == 0.0
    assert substrate_resistance_sp(0.02, 350e-4) == 0.02 * 350e-4
    with pytest.raises(ValueError):
        contact_resistance_sp(rho, x_max, 0.0)
    # R = V x_max / I: a series resistor rho/L_c (ohm cm per cm depth) adds exactly rho x_max/L_c
    v, i = 0.1, 3.7e-3  # V, A per cm depth
    r_cell = v * x_max / i
    r_lin_c = rho / l_c  # ohm cm
    r_tot = (v / i + r_lin_c) * x_max  # specific resistance of the half cell with the contact
    assert r_tot - r_cell == pytest.approx(contact_resistance_sp(rho, x_max, l_c), rel=1e-14)


def test_E3_contact_length_comes_from_the_model(sc):
    # default (1 Debye-length gap, no node at the gap): the contact starts at the first mesh node
    b = SiCDMOS(sc, BASE, 0)
    gap_um = float(sc.length_cm(1.0)) * 1e4
    x0 = b.dev.mesh.node_xy[b.dev.nodes["source"], 0].min()
    exp = float(sc.length_cm(b.prm.x_src1 - x0)) * 1e4
    assert source_contact_length_cm(b) * 1e4 == pytest.approx(exp, rel=1e-12)
    assert source_contact_length_cm(b) * 1e4 <= BASE.src_w - gap_um + 1e-9
