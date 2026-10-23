"""Fast tests of the stage-15b field-line machinery (widegap.fieldline).  < 60 s."""

import numpy as np
import pytest

from widegap.fieldline import (
    EXIT,
    LOW_FIELD,
    STAGNATION,
    STEP_LIMIT,
    FieldLineSet,
    NodalField,
    build_lines,
    contour_seeds,
    junction_contour,
    line_integral,
    trace,
)
from widegap.fv.avalanche import Chynoweth, uniform_field_integral
from widegap.offstate import UniformPN2D
from widegap.sic import SiC4H, TracedRegistry, tangent_predictor

pytestmark = pytest.mark.filterwarnings(
    "ignore::RuntimeWarning", "ignore::scipy.sparse.linalg.MatrixRankWarning"
)


def test_uniform_field_line_matches_closed_form():
    """Linear psi on a nonuniform tensor grid: straight lines, chord length exact, I_n closed form."""
    x = np.array([0.0, 0.3, 0.45, 1.0, 1.7, 2.0])
    y = np.array([-1.5, -1.1, -0.4, -0.2, 0.0])
    ex, ey = 0.35, -0.8
    psi = -(ex * x[:, None] + ey * y[None, :])  # E = -grad psi = (ex, ey)
    fld = NodalField(x, y, psi)
    ch = Chynoweth(a_n=6.0, b_n=2.0, a_p=2.0, b_p=3.0)
    seed = np.array([[1.1, -0.7]])
    lines, codes = build_lines(fld, seed, ds=0.004, max_steps=5000)
    pts = lines[0]
    assert (codes == EXIT).all()  # both ends leave the rectangle
    d = np.array([ex, ey]) / np.hypot(ex, ey)
    ts = [
        (b - s) / dd
        for s, dd, lo, hi in zip(seed[0], d, fld.lo, fld.hi, strict=True)
        for b in ((hi,) if dd > 0 else (lo,))
        if dd != 0
    ]
    t_fw = min(ts)  # distance to the boundary along +E
    ts = [
        (b - s) / (-dd)
        for s, dd, lo, hi in zip(seed[0], d, fld.lo, fld.hi, strict=True)
        for b in ((hi,) if -dd > 0 else (lo,))
        if dd != 0
    ]
    t_bw = min(ts)
    # p-side end first: the trace along +E (toward lower psi) is the start of the line
    np.testing.assert_allclose(pts[0], seed[0] + t_fw * d, atol=1e-12)
    np.testing.assert_allclose(pts[-1], seed[0] - t_bw * d, atol=1e-12)
    i_n, length = line_integral(fld, pts, ch)
    e = np.hypot(ex, ey)
    exact = uniform_field_integral(ch.alpha_n(np.array([e]))[0], ch.alpha_p(np.array([e]))[0],
                                   t_fw + t_bw)  # fmt: skip
    assert length == pytest.approx(t_fw + t_bw, rel=1e-12)
    assert i_n == pytest.approx(exact, rel=1e-5)  # trapezoid error O(ds^2)


def test_low_field_stop_and_step_limit():
    x = np.linspace(0.0, 1.0, 6)
    y = np.linspace(-1.0, 0.0, 6)
    # quadratic psi: |E| = 2|y + 1| vanishes on the bottom edge; lines along -E end in the low field
    psi = (y[None, :] + 1.0) ** 2 + 0.0 * x[:, None]
    fld = NodalField(x, y, psi)
    seeds = np.array([[0.5, -0.5]])
    _, _, cf = trace(fld, seeds, ds=2e-4, sign=+1.0, max_steps=5000)
    assert cf[0] == LOW_FIELD  # runs into the zero of E at y = -1 (|E| < 1e-3 max|E| within ds)
    _, _, cb = trace(fld, seeds, ds=2e-4, sign=-1.0, max_steps=5000)
    assert cb[0] == EXIT  # climbs to the surface y = 0
    _, _, c_lim = trace(fld, seeds, ds=2e-4, sign=+1.0, max_steps=3)
    assert c_lim[0] == STEP_LIMIT


class _Rotation:
    """Analytic circular field E = (-y, x): lines are circles, |E| = r."""

    lo = np.array([-2.0, -2.0])
    hi = np.array([2.0, 2.0])
    emax = 1.0

    def vec(self, p):
        p = np.atleast_2d(p)
        return np.column_stack([-p[:, 1], p[:, 0]])


def test_rk4_fourth_order():
    errs = []
    for n in (10, 20, 40):
        ds = 1.0 / n  # one radian of arc along the unit circle
        nodes, _, code = trace(_Rotation(), np.array([[1.0, 0.0]]), ds, +1.0, n)
        assert code[0] == STEP_LIMIT
        errs.append(float(np.hypot(*(nodes[-1, 0] - np.array([np.cos(1.0), np.sin(1.0)])))))
    assert errs[0] / errs[1] == pytest.approx(16.0, rel=0.2)
    assert errs[1] / errs[2] == pytest.approx(16.0, rel=0.2)


def test_junction_contour_seeds():
    x = np.linspace(0.0, 2.0, 9)
    y = np.linspace(-1.0, 0.0, 6)
    net = np.broadcast_to(y[None, :] + 0.37, (len(x), len(y))).copy()  # zero at y = -0.37
    chains = junction_contour(x, y, net)
    assert len(chains) == 1
    seeds, total = contour_seeds(chains, 8)
    assert total == pytest.approx(2.0)
    np.testing.assert_allclose(seeds[:, 1], -0.37, atol=1e-12)
    np.testing.assert_allclose(seeds[:, 0], (np.arange(8) + 0.5) * 0.25, atol=1e-12)


@pytest.fixture(scope="module")
def sc():
    return SiC4H(TracedRegistry.from_json())


def test_uniform_pn_lines_are_vertical(sc):
    """Laterally uniform p+n (12c-C5 structure): all field lines are vertical and the line integral
    reproduces the column ionization integral."""
    ch = sc.niwa()
    pn = UniformPN2D(sc, 1e16, level=0)
    mod = pn.model
    bc_of = pn.offstate_bc
    s = mod.solve_bc(mod.initial_guess(bc_of(0.0)), bc_of(0.0), max_step=20.0)
    assert s.report.converged
    x, v = s.X, 0.0
    h = pn.v_est / 3.0
    while v < 0.25 * pn.v_est:  # tangent-predictor continuation, step halving on a failed solve
        s = mod.solve_bc(tangent_predictor(mod, x, bc_of, v, h), bc_of(v + h), max_step=20.0)
        if not s.report.converged:
            h *= 0.5
            assert h > 1e-6 * pn.v_est
            continue
        x, v = s.X, v + h
    um = float(sc.l_d) * 1e4
    fls = FieldLineSet.build(pn.mesh, pn.nd, pn.na, ch, um, cut_x_um=[0.0], n_s=12, ds_um=0.005)
    out = fls.evaluate(x, keep_lines=True)
    for p in out["lines"]:
        assert np.ptp(p[:, 0]) < 1e-9
    assert (out["codes"] != STEP_LIMIT).all()
    i_col = pn.integral(x, ch)
    assert out["imax"] == pytest.approx(i_col, rel=0.05)
    assert out["I"].max() - out["I"].min() < 1e-6 * out["I"].max()


class _Sink:
    """Point sink E = -(r - r0) (|E| = distance to r0): lines converge on r0."""

    lo = np.array([-2.0, -2.0])
    hi = np.array([2.0, 2.0])
    emax = 1e-3  # low-field threshold 1e-6: the orbit of a fixed step ds = 0.01 stays above it

    def vec(self, p):
        return -(np.atleast_2d(p) - np.array([0.3, -0.2]))


def test_sink_stops_by_stagnation_before_step_limit():
    seeds = np.array([[0.8, 0.1]])
    nodes, n_nodes, code = trace(_Sink(), seeds, 0.01, +1.0, 5000)
    assert code[0] == STAGNATION
    assert n_nodes[0] < 5000
