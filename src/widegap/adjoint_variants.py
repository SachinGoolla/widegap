r"""Stage 22: adjoint gradient variants that drop individual terms of the exact adjoint.

Spec: docs/benchmarks/stage22_adjoint_simplifications.md (pre-registered; AUTHORITATIVE).

Nothing in the exact adjoint (:mod:`widegap.adjoint_channel`, :mod:`widegap.adjoint_sic`,
:mod:`widegap.adjoint_sic_oxide`) is modified.  Each variant is built from the SAME converged
state as the exact gradient E by masking one input of the unchanged
:meth:`widegap.adjoint_sic.SiCAdjoint.dF_dnd`:

* **A** (R_on only): drop ``F_mu mu_X`` (lagged channel mobility): the adjoint system is solved
  with ``jcore`` instead of ``jcore + F_mu mu_X`` (the multiplier ``w_lagged`` of the existing
  ``detail`` output); everything else as E.
* **B** (R_on only): ``d mu_bulk / d N := 0`` in ``dF/dtheta`` (``ev.dmu_dnd`` replaced by zeros;
  E already zeroes it on the channel edges).
* **C** (all figures of merit): complete-ionization doping derivative in ``dF/dtheta``,
  ``d N_D^+ / d N_D := 1`` (the Poisson-row ionized fraction ``sum_l f_l sigma(.)`` is replaced
  by 1: ``mod.ion.donors`` is emptied for the call, which is the ``frac = ones`` branch of
  ``dF_dnd``); ``d N_A^- / d N_A := 1`` has no effect because ``N_A`` does not depend on theta.
  The mobility part of ``dF/dtheta`` (it depends on ``N_D + N_A`` as a total, not on the ionized
  fraction) is untouched.
* **D**: A + B + C.

Off-state figures of merit (soft-min BV, E_ox) have zero-current electron rows (Dirichlet), so
the adjoint multiplier vanishes there and the mobility part of ``dF/dtheta`` is exactly zero:
A and B do not apply (``B`` is nevertheless evaluated to record that its difference is 0), and
D equals C.

Because the variants call the unchanged ``dF_dnd`` on the unchanged state, "E re-assembled"
equals the gradient returned by the existing code bit for bit; this is checked at every
assembly (``VariantStats.n_bitwise_checks``) and raises on a mismatch.
"""

from __future__ import annotations

from contextlib import contextmanager

import numpy as np

from widegap.adjoint2d import solve_transposed

VARIANTS = ("E", "A", "B", "C", "D")
APPLIES = {"A": ("R",), "B": ("R",), "C": ("R", "BV", "E_ox"), "D": ("R", "BV", "E_ox")}


class VariantCheckError(AssertionError):
    """E re-assembled through the variant machinery differs from the existing exact gradient."""


@contextmanager
def complete_ionization(ev):
    """``d N_D^+/d N_D := 1`` in ``ev.dF_dnd`` (donor ionization levels hidden for the call)."""
    saved = ev.mod.ion.donors
    ev.mod.ion.donors = []
    try:
        yield
    finally:
        ev.mod.ion.donors = saved


@contextmanager
def bulk_mobility_derivative(ev, dmu):
    """``ev.dmu_dnd`` replaced by ``dmu`` for the call."""
    saved = ev.dmu_dnd
    ev.dmu_dnd = dmu
    try:
        yield
    finally:
        ev.dmu_dnd = saved


def dF_dnd_variant(ev, w, X, dmu=None, complete_ion=False):
    """``ev.dF_dnd(w, X)`` with the bulk-mobility derivative ``dmu`` (None = as set) and/or the
    complete-ionization donor derivative.  With the defaults this IS ``ev.dF_dnd(w, X)``."""
    if dmu is None and not complete_ion:
        return ev.dF_dnd(w, X)
    if complete_ion and dmu is None:
        with complete_ionization(ev):
            return ev.dF_dnd(w, X)
    with bulk_mobility_derivative(ev, dmu):
        if complete_ion:
            with complete_ionization(ev):
                return ev.dF_dnd(w, X)
        return ev.dF_dnd(w, X)


class VariantStats:
    """Counters of the E-re-assembly checks (Q2 evidence)."""

    def __init__(self):
        self.n_bitwise_checks = 0
        self.n_assemblies = 0

    def as_dict(self):
        return {"n_bitwise_checks": self.n_bitwise_checks, "n_assemblies": self.n_assemblies}


def _same(a, b, what):
    if not np.array_equal(a, b):
        raise VariantCheckError(f"{what}: re-assembled E differs from the existing gradient")


# ------------------------------------------------------------------------------ R_on,ch
def ron_gradients(cr, theta, variants=VARIANTS, stats: VariantStats | None = None, check=True):
    """R_on,ch and the gradient of every requested variant at ONE converged channel state.

    Returns dict with ``R``, ``I_scaled``, ``grad`` {variant: dR/dtheta}, ``g_nd`` {variant: ...},
    ``state_sha`` (hash of the converged X), ``guard``, ``out`` (the raw ``ChannelRon.ron``
    output)."""
    stats = stats or VariantStats()
    ev, cmap = cr.ev, cr.cmap
    out = cr.ron(theta, grad=True, detail=True)  # E and A (lagged) from the existing code
    fac = -out["R"] / out["I_scaled"]
    X = out["state"]["X"]
    saved_mu = ev.mod.mu_edge
    ev.mod.mu_edge = np.asarray(out["picard"]["mu_edge"], dtype=float).copy()  # converged state
    try:
        masked = ev.dmu_dnd.copy()
        masked[cmap.edges] = 0.0  # as in ChannelRon.ron
        zero = np.zeros_like(ev.dmu_dnd)
        w, w_l = out["w"], out["w_lagged"]

        def gnd(wv, dmu, ion):
            return fac * dF_dnd_variant(ev, wv, X, dmu=dmu, complete_ion=ion)

        stats.n_assemblies += 1
        spec = {"E": (w, masked, False), "A": (w_l, masked, False), "B": (w, zero, False),
                "C": (w, masked, True), "D": (w_l, zero, True)}  # fmt: skip
        g_nd = {v: gnd(*spec[v]) for v in variants}
        if check:
            _same(gnd(*spec["E"]), out["g_nd"], "R_on E")
            _same(gnd(*spec["A"]), out["g_nd_lagged"], "R_on A (lagged)")
            stats.n_bitwise_checks += 2
    finally:
        ev.mod.mu_edge = saved_mu
    # mu_edge restored to the bulk mobility exactly as ChannelRon.ron leaves it
    grads = {v: ev.dnd_dtheta.T @ g_nd[v] for v in variants}
    if check and "E" in grads:
        _same(grads["E"], out["grad"], "R_on E grad")
        stats.n_bitwise_checks += 1
    return {
        "R": out["R"],
        "I_scaled": out["I_scaled"],
        "grad": grads,
        "g_nd": g_nd,
        "state_sha": _sha(X),
        "guard": out["guard"],
        "out": out,
    }


def _sha(x):
    import hashlib

    return hashlib.sha256(np.ascontiguousarray(x, dtype=float).tobytes()).hexdigest()[:16]


# ------------------------------------------------------------------------------ soft-min BV
def objective_gradients(ev, theta, variants=("E", "B", "C", "D"), stats: VariantStats | None = None,
                        check=True):  # fmt: skip
    """Soft-min BV objective and the gradient of every requested variant, assembled from the
    converged per-cut states of one ``ev.objective(theta)`` call (no re-solve)."""
    stats = stats or VariantStats()
    res = ev.objective(theta, grad=True)
    sc_v = float(ev.sc.vt)
    mod = ev.mod
    zero = np.zeros_like(ev.dmu_dnd)
    gk = {v: [] for v in variants}
    lams, dvs = [], []
    stats.n_assemblies += 1
    for o in res["outs"]:
        st = o["state"]
        lam = solve_transposed(mod.jacobian_bc(st.X, st.rows), st.g_psi)
        di_dv = float(lam[st.rows] @ st.dvals)
        lam[st.rows] = 0.0
        lams.append(lam.copy())
        dvs.append(di_dv)
        for v in variants:
            g = dF_dnd_variant(ev, -lam, st.X, dmu=zero if v in ("B", "D") else None,
                               complete_ion=v in ("C", "D"))  # fmt: skip
            gk[v].append(-g / di_dv * sc_v)
        if check:
            e = dF_dnd_variant(ev, -lam, st.X)
            _same(-e / di_dv * sc_v, o["g_nd"], "BV cut E")
            stats.n_bitwise_checks += 1
    w = res["w"]
    g_nd = {v: sum(wk * g for wk, g in zip(w, gk[v], strict=True)) for v in variants}
    # same summation order as CornerAdjoint.objective: soft-min weights times per-cut gradients
    grads = {v: sum(wk * (ev.dnd_dtheta.T @ g) for wk, g in zip(w, gk[v], strict=True))
             for v in variants}  # fmt: skip
    if check and "E" in grads:
        _same(grads["E"], res["grad"], "BV_obj E")
        stats.n_bitwise_checks += 1
    return {
        "BV": res["BV"],
        "bv_k": res["bv_k"],
        "w": w,
        "grad": grads,
        "g_nd": g_nd,
        "res": res,
        "lams": lams,
        "di_dv": dvs,
        "sc_v": sc_v,
        "g_nd_cuts": gk,
    }


# ------------------------------------------------------------------------------ E_ox (fields)
def fields_gradients(ev, theta, which=("E_ox",), variants=("E", "B", "C", "D"),
                     stats: VariantStats | None = None, check=True):  # fmt: skip
    """Smooth-max field measures at V_blk and the gradient of every requested variant, from the
    converged state of one ``ev.fields(theta)`` call."""
    stats = stats or VariantStats()
    o = ev.fields(theta, grad=True, which=which)
    X = o["X"]
    mod = ev.mod
    psi = X[0::3]
    rows, _ = ev._off_bc(float(ev.sc.volt_s(o["v_volt"])))
    jac = mod.jacobian_bc(X, rows)
    zero = np.zeros_like(ev.dmu_dnd)
    out = {"X": X}
    stats.n_assemblies += 1
    for name in which:
        _, gpsi = ev.measures[name].value_grad(psi)
        g = np.zeros(len(X))
        g[0::3] = gpsi
        lam = solve_transposed(jac, g)
        lam[rows] = 0.0
        g_nd = {}
        for v in variants:
            g_nd[v] = dF_dnd_variant(ev, -lam, X, dmu=zero if v in ("B", "D") else None,
                                     complete_ion=v in ("C", "D"))  # fmt: skip
        if check:
            _same(dF_dnd_variant(ev, -lam, X), o[name]["g_nd"], f"{name} E")
            stats.n_bitwise_checks += 1
        out[name] = {"S": o[name]["S"], "hard_max": o[name]["hard_max"], "g_nd": g_nd,
                     "grad": {v: ev.dnd_dtheta.T @ g_nd[v] for v in variants},
                     "E_grad_existing": o[name]["grad"], "lam": lam}  # fmt: skip
    return out


# ------------------------------------------------------------------------------ KKT residual
def kkt_residual(g_obj, constraints, tol=1e-6):
    r"""Stationarity residual of ``max f s.t. c_i >= 0`` at a design, from gradients ``g_obj`` of
    ``f`` and the constraint set ``constraints = [(name, values (m,), jac (m, n)), ...]``.

    Active constraints are those with ``c_i <= tol``.  With ``A`` the matrix of their gradients
    (as columns) the first-order condition for a maximum is ``g_obj + A lam = 0, lam >= 0``;
    the residual is ``min_{lam >= 0} || g_obj + A lam ||`` (non-negative least squares).
    Returns residual, residual / ||g_obj||, the number of active constraints and multipliers."""
    from scipy.optimize import nnls

    g = np.asarray(g_obj, dtype=float)
    cols, names = [], []
    for name, vals, jac in constraints:
        vals = np.atleast_1d(np.asarray(vals, dtype=float))
        jac = np.atleast_2d(np.asarray(jac, dtype=float))
        for i in np.flatnonzero(vals <= tol):
            cols.append(jac[i])
            names.append((name, int(i)))
    if not cols:
        res = float(np.linalg.norm(g))
        lam = np.zeros(0)
    else:
        a = np.stack(cols, axis=1)
        lam, _ = nnls(a, -g, maxiter=50 * max(a.shape))
        res = float(np.linalg.norm(g + a @ lam))
    return {"residual": res, "ratio": res / float(np.linalg.norm(g)), "n_active": len(cols),
            "active": names, "multipliers": lam.tolist(), "tol": float(tol),
            "grad_norm": float(np.linalg.norm(g))}  # fmt: skip


# ------------------------------------------------------------------------------ SLSQP adapter
class VariantComposite:
    """Drop-in for ``Eval17`` (BV_obj / E_ox on the off-state evaluator, R_on,ch on the channel
    model) whose gradients are those of ``variant``; function values are the exact ones."""

    def __init__(self, ev_off, cr, variant="D", check=True):
        self.ev_off, self.cr, self.sc = ev_off, cr, ev_off.sc
        self.variant, self.check = variant, check
        self.stats = VariantStats()

    def objective(self, th, grad=True):
        if not grad:
            return self.ev_off.objective(th, grad=False)
        o = objective_gradients(self.ev_off, th, (self.variant,), self.stats, self.check)
        res = dict(o["res"])
        res["grad"] = o["grad"][self.variant]
        return res

    def fields(self, th, grad=True, which=("E_ox",), v_volt=None):
        if not grad:
            return self.ev_off.fields(th, grad=False, which=which, v_volt=v_volt)
        o = fields_gradients(self.ev_off, th, which, (self.variant,), self.stats, self.check)
        return {k: {"S": o[k]["S"], "hard_max": o[k]["hard_max"], "grad": o[k]["grad"][self.variant]}
                for k in which}  # fmt: skip

    def ron(self, th, grad=True):
        if not grad:
            return {"R": self.cr.ron(th, grad=False)["R"]}
        o = ron_gradients(self.cr, th, (self.variant,), self.stats, self.check)
        return {"R": o["R"], "grad": o["grad"][self.variant]}


__all__ = [
    "APPLIES",
    "VARIANTS",
    "VariantCheckError",
    "VariantComposite",
    "VariantStats",
    "bulk_mobility_derivative",
    "complete_ionization",
    "dF_dnd_variant",
    "fields_gradients",
    "kkt_residual",
    "objective_gradients",
    "ron_gradients",
]
