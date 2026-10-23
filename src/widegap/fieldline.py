r"""Stage 15b: field-line ionization integrals (post-hoc re-evaluation of the avalanche BV).

Spec: docs/benchmarks/stage15_caveat_closure.md, section 15b (pre-registered; AUTHORITATIVE).

* **Nodal field.**  ``E = -grad psi`` is reconstructed at the nodes of the SiC part (``y <= 0``) of
  the tensor mesh by second-order finite differences on the nonuniform grid (``numpy.gradient``,
  ``edge_order=2``; one-sided at the boundary), and interpolated bilinearly inside each cell.
* **Field lines.**  ``dr/ds = +-E/|E|`` with fixed-step RK4 (step ``ds``), traced from every seed in
  both directions.  A trace stops (i) on leaving the SiC rectangle (the last step is shortened to
  end exactly on the boundary) or (ii) at the first node where ``|E| < e_frac * max|E|``
  (``e_frac = 1e-3``; that node is kept).  A line still running after ``max_steps`` is reported
  (``step_limit``) and is a failure of criterion F5.
* **Step.**  The spec leaves ``ds`` free.  A line that ends in a field-reversal sink of the nodal
  field (the p-body next to the n+ source junction, where dE/dy ~ q N_A / eps) keeps a
  ``|E| ~ (dE/dy) ds`` above the low-field threshold and orbits the sink until ``max_steps``
  (F5 violation) unless the stagnation rule (net displacement over the last 200 steps < 10 ds,
  spec Amendment 2) stops the trace first; the integral up to the stop is kept.  For the 13d initial design this happened at ``ds`` = 0.01, 0.005 and 0.0025 um
  (marched states of 180-200 V) and not at 0.00125 or 0.000625 um (40 contour seeds + cut
  crossings, all marched states); the default is 0.00125 um and F2 halves it.
* **Orientation.**  Electrons drift toward higher psi, so every line is ordered by increasing psi:
  index 0 is the p-side (hole-injecting) end, as for the vertical cuts of
  :mod:`widegap.adjoint_sic_corner` (anode at the first node, depth increasing).
* **Integral.**  ``I_n = int alpha_n exp(-int_0^s (alpha_n - alpha_p) ds') ds`` by the trapezoidal
  rule on the RK4 nodes (chord lengths between nodes), with the coefficient set of the evaluator
  (``ev.ch``, the registry Niwa set) evaluated on the interpolated ``|E|``.
* **Seeds.**  ``N_s`` points uniformly spaced in arc length along the metallurgical junction
  (``N_D - N_A = 0``, piecewise-linear zero contour of the nodal net doping on the SiC cells;
  seeds at arc positions ``(k + 1/2) L / N_s``), plus the junction crossings of the vertical cut
  lines (bilinear interpolation of the net doping between mesh columns, as ``DenseCuts``).
* **Field-line BV.**  Smallest ``V`` with ``max_lines I_n = 1``: the stage-13b ``dense_bv`` recipe
  (20 V zero-current march with tangent predictor, then regula falsi on ``ln I`` clipped to
  [2 %, 98 %] of the bracket, stopped when the bracket is <= ``tol_v`` = 0.02 V; BV from the final
  bracket by log-linear interpolation).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from widegap.adjoint_sic import EvalError
from widegap.adjoint_sic_corner import march_states, state_at
from widegap.sic import tangent_predictor

E_FRAC = 1e-3  # low-field stopping threshold (fraction of the global max |E|)
V_STEP_V = 20.0  # march step [V], as dense_bv
TOL_V = 0.02  # bracket tolerance [V], as dense_bv (tol_v)
CLIP = 0.02  # secant point kept within [CLIP, 1 - CLIP] of the bracket, as dense_bv
MAX_REFINE = 400  # as dense_bv

EXIT, LOW_FIELD, STEP_LIMIT, STAGNATION = 0, 1, 2, 3
STOP_NAMES = {
    EXIT: "left_sic",
    LOW_FIELD: "low_field",
    STAGNATION: "stagnation",
    STEP_LIMIT: "step_limit",
}
STAG_STEPS = 200  # stagnation: net displacement over the last STAG_STEPS steps < STAG_FACTOR * ds
STAG_FACTOR = 10.0


# ------------------------------------------------------------------------------ the field
class NodalField:
    """Bilinear-interpolated nodal field ``E = -grad psi`` on the SiC rectangle.

    ``x`` (nx+1) and ``y`` (ny+1, ascending, last = 0) are the SiC grid; ``psi`` is (nx+1, ny+1).
    """

    def __init__(self, x, y, psi):
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        psi = np.asarray(psi, float)
        self.ex = -np.gradient(psi, self.x, axis=0, edge_order=2)
        self.ey = -np.gradient(psi, self.y, axis=1, edge_order=2)
        self.lo = np.array([self.x[0], self.y[0]])
        self.hi = np.array([self.x[-1], self.y[-1]])
        self.emax = float(np.max(np.hypot(self.ex, self.ey)))

    @classmethod
    def from_mesh(cls, mesh, psi_nodes):
        """SiC part (``y <= 0``) of the 2D tensor mesh; ``psi_nodes`` in the mesh node order."""
        psi = np.asarray(psi_nodes, float).reshape(mesh.nx + 1, mesh.ny + 1)
        j = mesh.y <= 1e-12 * max(1.0, float(np.max(np.abs(mesh.y))))
        return cls(mesh.x, mesh.y[j], psi[:, j])

    def vec(self, p):
        """(n, 2) field at points ``p`` (clipped into the rectangle)."""
        p = np.clip(np.atleast_2d(p), self.lo, self.hi)
        x, y = self.x, self.y
        i = np.clip(np.searchsorted(x, p[:, 0], side="right") - 1, 0, len(x) - 2)
        j = np.clip(np.searchsorted(y, p[:, 1], side="right") - 1, 0, len(y) - 2)
        fx = (p[:, 0] - x[i]) / (x[i + 1] - x[i])
        fy = (p[:, 1] - y[j]) / (y[j + 1] - y[j])

        def bil(f):
            return (
                (1 - fx) * (1 - fy) * f[i, j] + fx * (1 - fy) * f[i + 1, j]
                + (1 - fx) * fy * f[i, j + 1] + fx * fy * f[i + 1, j + 1]
            )  # fmt: skip

        return np.column_stack([bil(self.ex), bil(self.ey)])

    def mag(self, p):
        return np.hypot(*self.vec(p).T)


# ------------------------------------------------------------------------------ tracing
def trace(field, seeds, ds, sign, max_steps, e_frac=E_FRAC):
    """RK4 trace of ``dr/ds = sign E/|E|`` from every seed (vectorized over the seeds).

    ``field`` needs ``vec(points)``, ``lo``, ``hi``, ``emax``.  Returns
    ``(nodes, n_nodes, code)``: ``nodes`` (K, m, 2) padded with NaN, line ``i`` being
    ``nodes[:n_nodes[i], i]`` (the seed first), ``code`` the stop rule of each line."""
    pos = np.array(seeds, float).reshape(-1, 2)
    m = len(pos)
    thr = e_frac * field.emax
    active = np.ones(m, bool)
    code = np.full(m, STEP_LIMIT)
    n_nodes = np.ones(m, int)
    hist = [pos.copy()]
    ring = np.empty((STAG_STEPS + 1, m, 2))  # positions of the last STAG_STEPS + 1 steps
    ring[0] = pos

    def unit(p):
        e = field.vec(p)
        mag = np.hypot(e[:, 0], e[:, 1])
        return sign * e / np.where(mag > 0, mag, 1.0)[:, None]

    for k_step in range(1, max_steps + 1):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        p = pos[idx]
        low = np.hypot(*field.vec(p).T) < thr
        code[idx[low]] = LOW_FIELD
        active[idx[low]] = False
        idx, p = idx[~low], p[~low]
        if idx.size == 0:
            break
        k1 = unit(p)
        k2 = unit(p + 0.5 * ds * k1)
        k3 = unit(p + 0.5 * ds * k2)
        k4 = unit(p + ds * k3)
        d = ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        new = p + d
        tol = 1e-12 * (field.hi - field.lo)  # round-off sliding along a boundary is not an exit
        v_hi, v_lo = new > field.hi + tol, new < field.lo - tol
        out = np.any(v_hi | v_lo, axis=1)
        if out.any():
            dd = d[out]
            with np.errstate(divide="ignore", invalid="ignore"):
                f_hi = np.where(v_hi[out], (field.hi - p[out]) / dd, np.inf)
                f_lo = np.where(v_lo[out], (field.lo - p[out]) / dd, np.inf)
            frac = np.clip(np.minimum(f_hi, f_lo).min(axis=1), 0.0, 1.0)
            new[out] = p[out] + frac[:, None] * dd
            code[idx[out]] = EXIT
            active[idx[out]] = False
        new = np.clip(new, field.lo, field.hi)
        pos[idx] = new
        ring[k_step % (STAG_STEPS + 1), idx] = new
        if k_step >= STAG_STEPS:
            old = ring[(k_step - STAG_STEPS) % (STAG_STEPS + 1), idx]
            stag = np.hypot(*(new - old).T) < STAG_FACTOR * ds
            stag &= ~out
            code[idx[stag]] = STAGNATION
            active[idx[stag]] = False
        rec = np.full((m, 2), np.nan)
        rec[idx] = new
        hist.append(rec)
        n_nodes[idx] += 1
    return np.stack(hist), n_nodes, code


def build_lines(field, seeds, ds, max_steps, e_frac=E_FRAC):
    """Both-direction field lines, each ordered from the p-side end (lowest psi): the trace along
    ``+E`` (toward lower psi) reversed, the seed, then the trace along ``-E``.

    Returns ``(lines, codes)``: a list of (n_i, 2) arrays and an (m, 2) array with the stop rule
    of the (p-side, n-side) end."""
    fw, nf, cf = trace(field, seeds, ds, +1.0, max_steps, e_frac)
    bw, nb, cb = trace(field, seeds, ds, -1.0, max_steps, e_frac)
    lines = [
        np.vstack([fw[: nf[i], i][::-1], bw[1 : nb[i], i]])
        for i in range(len(np.atleast_2d(seeds)))
    ]
    return lines, np.column_stack([cf, cb])


def line_integral(field, pts, ch):
    """(I_n, length) of one line by the trapezoidal rule (chord lengths between the nodes)."""
    seg = np.hypot(*np.diff(pts, axis=0).T)
    if len(seg) == 0:
        return 0.0, 0.0
    e = field.mag(pts)
    an, ap = ch.alpha_n(e), ch.alpha_p(e)
    d = an - ap
    big_a = np.concatenate(([0.0], np.cumsum(0.5 * (d[:-1] + d[1:]) * seg)))
    f = an * np.exp(-big_a)
    return float(np.sum(0.5 * (f[:-1] + f[1:]) * seg)), float(seg.sum())


# ------------------------------------------------------------------------------ seeds
def junction_contour(x, y, net):
    """Zero contour of the nodal ``net`` (nx+1, ny+1) with linear interpolation on the grid edges.

    Returns a list of polylines (each (n, 2)), oriented to start at the smaller (x, y) end and
    sorted by their start point (deterministic)."""
    nx1, ny1 = net.shape
    pos = net > 0
    pts: dict = {}
    adj: dict = {}

    def key_pt(kind, i, j):
        k = (kind, i, j)
        if k not in pts:
            if kind == "H":
                a, b = net[i, j], net[i + 1, j]
                t = a / (a - b)
                pts[k] = (x[i] + t * (x[i + 1] - x[i]), y[j])
            else:
                a, b = net[i, j], net[i, j + 1]
                t = a / (a - b)
                pts[k] = (x[i], y[j] + t * (y[j + 1] - y[j]))
        return k

    def link(k1, k2):
        adj.setdefault(k1, []).append(k2)
        adj.setdefault(k2, []).append(k1)

    for i in range(nx1 - 1):
        for j in range(ny1 - 1):
            c = (pos[i, j], pos[i + 1, j], pos[i + 1, j + 1], pos[i, j + 1])  # BL BR TR TL
            if len({*c}) == 1:
                continue
            edges = []  # in order around the cell: bottom, right, top, left
            if c[0] != c[1]:
                edges.append(("b", key_pt("H", i, j)))
            if c[1] != c[2]:
                edges.append(("r", key_pt("V", i + 1, j)))
            if c[2] != c[3]:
                edges.append(("t", key_pt("H", i, j + 1)))
            if c[3] != c[0]:
                edges.append(("l", key_pt("V", i, j)))
            if len(edges) == 2:
                link(edges[0][1], edges[1][1])
            else:  # saddle: the centre value (mean of the corners) decides the connection
                kk = dict(edges)
                centre = 0.25 * (net[i, j] + net[i + 1, j] + net[i + 1, j + 1] + net[i, j + 1])
                if (centre > 0) == c[0]:  # corners BL, TR connected: cut off BR and TL
                    link(kk["b"], kk["r"])
                    link(kk["t"], kk["l"])
                else:
                    link(kk["b"], kk["l"])
                    link(kk["r"], kk["t"])
    seen: set = set()
    chains = []

    def walk(start):
        chain, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nxt = [n for n in adj[cur] if n != prev and n not in seen]
            if not nxt:
                return chain
            prev, cur = cur, nxt[0]
            seen.add(cur)
            chain.append(cur)

    for k in sorted(adj):
        if len(adj[k]) == 1 and k not in seen:
            chains.append(walk(k))
    for k in sorted(adj):
        if k not in seen:
            chains.append(walk(k))
    out = []
    for ch in chains:
        p = np.array([pts[k] for k in ch])
        if (p[0, 0], p[0, 1]) > (p[-1, 0], p[-1, 1]):
            p = p[::-1]
        out.append(p)
    out.sort(key=lambda p: (round(p[0, 0], 12), round(p[0, 1], 12)))
    return out


def contour_seeds(chains, n_s):
    """``n_s`` points at arc positions ``(k + 1/2) L / n_s`` along the concatenated polylines;
    returns (points, total length)."""
    seg = [np.hypot(*np.diff(c, axis=0).T) for c in chains]
    lens = np.array([s.sum() for s in seg])
    total = float(lens.sum())
    if total <= 0:
        raise ValueError("no junction contour in the SiC region")
    s_k = (np.arange(n_s) + 0.5) * total / n_s
    off = np.concatenate(([0.0], np.cumsum(lens)))
    out = np.empty((n_s, 2))
    comp = np.clip(np.searchsorted(off, s_k, side="right") - 1, 0, len(chains) - 1)
    for c in range(len(chains)):
        sel = comp == c
        if not sel.any():
            continue
        cum = np.concatenate(([0.0], np.cumsum(seg[c])))
        t = np.clip(s_k[sel] - off[c], 0.0, cum[-1])
        out[sel, 0] = np.interp(t, cum, chains[c][:, 0])
        out[sel, 1] = np.interp(t, cum, chains[c][:, 1])
    return out, total


def cut_crossings(x, y, net, x_cuts):
    """Junction crossings of vertical lines at ``x_cuts`` (bilinear net doping between columns)."""
    out = []
    for xc in np.asarray(x_cuts, float):
        il = int(np.clip(np.searchsorted(x, xc, side="right") - 1, 0, len(x) - 2))
        w = float((xc - x[il]) / (x[il + 1] - x[il]))  # as line_operator of DenseCuts
        col = (1.0 - w) * net[il] + w * net[il + 1]
        for j in range(len(y) - 1):
            if (col[j] > 0) != (col[j + 1] > 0):
                t = col[j] / (col[j] - col[j + 1])
                out.append((float(xc), float(y[j] + t * (y[j + 1] - y[j]))))
    return np.array(out).reshape(-1, 2)


# ------------------------------------------------------------------------------ evaluator
@dataclass
class FieldLineSet:
    """Seeds + tracing parameters for one design on one mesh; evaluates ``I_n`` of all lines."""

    mesh: object
    ch: object
    seeds: np.ndarray
    kind: np.ndarray  # 0 = contour seed, 1 = cut crossing
    ds: float  # step [scaled length]
    max_steps: int
    um_per_unit: float
    e_frac: float = E_FRAC
    info: dict = field(default_factory=dict)
    stats: dict = field(
        default_factory=lambda: {"n_eval": 0, "n_lines_traced": 0, "step_limit": 0,
                                 "stops": {n: 0 for n in STOP_NAMES.values()}}
    )  # fmt: skip

    @classmethod
    def build(cls, mesh, nd, na, ch, um_per_unit, cut_x_um=(), n_s=400, ds_um=0.00125,
              max_steps=None, e_frac=E_FRAC):  # fmt: skip
        j = mesh.y <= 1e-12 * max(1.0, float(np.max(np.abs(mesh.y))))
        net = (np.asarray(nd) - np.asarray(na)).reshape(mesh.nx + 1, mesh.ny + 1)[:, j]
        x, y = mesh.x, mesh.y[j]
        chains = junction_contour(x, y, net)
        s_c, total = contour_seeds(chains, n_s)
        s_x = cut_crossings(x, y, net, np.asarray(cut_x_um, float) / um_per_unit)
        seeds = np.vstack([s_c, s_x])
        kind = np.r_[np.zeros(len(s_c), int), np.ones(len(s_x), int)]
        ds = ds_um / um_per_unit
        if max_steps is None:
            max_steps = int(np.ceil(4.0 * ((x[-1] - x[0]) + (y[-1] - y[0])) / ds))
        info = {"n_contour_seeds": int(len(s_c)), "n_cut_seeds": int(len(s_x)),
                "contour_length_um": total * um_per_unit, "n_contour_chains": len(chains),
                "ds_um": ds_um, "max_steps": int(max_steps), "e_frac": e_frac}  # fmt: skip
        return cls(mesh, ch, seeds, kind, ds, int(max_steps), um_per_unit, e_frac, info)

    def lines(self, psi_nodes):
        fld = NodalField.from_mesh(self.mesh, psi_nodes)
        lines, codes = build_lines(fld, self.seeds, self.ds, self.max_steps, self.e_frac)
        return fld, lines, codes

    def evaluate(self, X, keep_lines=False):
        """``I_n`` of every line at the full state vector ``X`` (psi = X[0::3])."""
        fld, lines, codes = self.lines(np.asarray(X)[0::3])
        res = [line_integral(fld, p, self.ch) for p in lines]
        i_n = np.array([r[0] for r in res])
        st = self.stats
        st["n_eval"] += 1
        st["n_lines_traced"] += 2 * len(lines)
        for c, name in STOP_NAMES.items():
            st["stops"][name] += int(np.sum(codes == c))
        st["step_limit"] = st["stops"]["step_limit"]
        out = {"I": i_n, "length": np.array([r[1] for r in res]), "codes": codes,
               "imax": float(i_n.max()), "argmax": int(np.argmax(i_n))}  # fmt: skip
        if keep_lines:
            out.update({"field": fld, "lines": lines})
        return out

    def imax(self, X):
        return self.evaluate(X)["imax"]

    def line_record(self, ev_out, k):
        """Seed position and end points [um] of line ``k`` (limiting-line report)."""
        u = self.um_per_unit
        pts = ev_out["lines"][k] * u
        return {"index": int(k), "kind": "cut_crossing" if self.kind[k] else "contour",
                "seed_um": (self.seeds[k] * u).tolist(), "p_end_um": pts[0].tolist(),
                "n_end_um": pts[-1].tolist(), "length_um": float(ev_out["length"][k]) * u,
                "I_n": float(ev_out["I"][k]),
                "stop_p_n": [STOP_NAMES[int(c)] for c in ev_out["codes"][k]]}  # fmt: skip


# ------------------------------------------------------------------------------ BV search
def dmos_states(ev, theta, fls: FieldLineSet, v_step_v=V_STEP_V):
    """Zero-current march (``march_states``) until ``max I_n >= 1``; (states, state_at(v_s))."""
    states = march_states(ev, theta, lambda X: np.array([fls.imax(X)]), v_step_v=v_step_v,
                          extend=1.0)  # fmt: skip
    return states, lambda v: state_at(ev, states, v)


def diode_states(pn, fls: FieldLineSet, v_step_s, v_max_s, max_step=20.0):
    """The same march for :class:`widegap.offstate.UniformPN2D` (corrected off state)."""
    mod = pn.model
    bc_of = pn.offstate_bc
    s = mod.solve_bc(mod.initial_guess(bc_of(0.0)), bc_of(0.0), max_step=max_step)
    if not s.report.converged:
        raise EvalError("zero-bias solve did not converge")
    v, x, h = 0.0, s.X.copy(), float(v_step_s)
    states = []
    while True:
        t = tangent_predictor(mod, x, bc_of, v, 1.0) - x
        states.append((v, x.copy(), t))
        if fls.imax(x) >= 1.0:
            break
        if v > v_max_s:
            raise EvalError("I = 1 not reached below the voltage ceiling")
        while True:
            s = mod.solve_bc(x + h * t, bc_of(v + h), max_step=max_step)
            if s.report.converged:
                break
            h *= 0.5
            if h < 1e-9:
                raise EvalError(f"march failed at V = {v}")
        v, x = v + h, s.X.copy()

    def at(vs):
        k = max(i for i, st in enumerate(states) if st[0] <= vs + 1e-12)
        v0, x0, t0 = states[k]
        if abs(vs - v0) < 1e-12:
            return x0
        sol = mod.solve_bc(x0 + (vs - v0) * t0, bc_of(vs), max_step=max_step)
        if not sol.report.converged:
            raise EvalError(f"refinement solve failed at V = {vs}")
        return sol.X

    return states, at


def refine_bv(states, at, fls: FieldLineSet, volt_s, volt_phys, tol_v=TOL_V):
    """Regula falsi on ``ln max I`` (clipped to [CLIP, 1 - CLIP] of the bracket) until the bracket
    is <= ``tol_v`` [V]; BV by log-linear interpolation of the final bracket.  The root is bracketed
    by the last two marched states."""
    if len(states) < 2:
        raise EvalError("I >= 1 already at zero bias")
    lni = lambda X: float(np.log(max(fls.imax(X), 1e-300)))  # noqa: E731
    (v0, x0, _), (v1, x1, _) = states[-2], states[-1]
    l0, l1 = lni(x0), lni(x1)
    if not (l0 < 0.0 <= l1):
        raise EvalError("root not bracketed by the marched states")
    tol_s = float(volt_s(tol_v))
    n_solves, converged = 0, False
    for _ in range(MAX_REFINE):
        if v1 - v0 <= tol_s:
            converged = True
            break
        vn = v0 + (0.0 - l0) * (v1 - v0) / (l1 - l0)
        vn = float(np.clip(vn, v0 + CLIP * (v1 - v0), v1 - CLIP * (v1 - v0)))
        ln = lni(at(vn))
        n_solves += 1
        if ln < 0.0:
            v0, l0 = vn, ln
        else:
            v1, l1 = vn, ln
    v_root = v0 + (0.0 - l0) * (v1 - v0) / (l1 - l0)
    return {"bv_s": float(v_root), "bv": float(volt_phys(v_root)), "converged": converged,
            "n_refine_solves": n_solves, "n_march_states": len(states),
            "bracket_V": [float(volt_phys(v0)), float(volt_phys(v1))]}  # fmt: skip


__all__ = [
    "CLIP",
    "E_FRAC",
    "EXIT",
    "LOW_FIELD",
    "STEP_LIMIT",
    "STOP_NAMES",
    "TOL_V",
    "V_STEP_V",
    "FieldLineSet",
    "NodalField",
    "build_lines",
    "contour_seeds",
    "cut_crossings",
    "diode_states",
    "dmos_states",
    "junction_contour",
    "line_integral",
    "refine_bv",
    "trace",
]
