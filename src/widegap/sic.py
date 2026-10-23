r"""Stage 12: sourced-parameter 4H-SiC device physics (physical -> scaled mapping).

Spec: docs/benchmarks/stage12_sic_devices.md (pre-registered; AUTHORITATIVE).

Every numerical material input is read through :class:`TracedRegistry`, a thin wrapper
around :class:`widegap.materials.registry.Registry` that (i) enforces the stage-7 provenance
rules (an unusable entry raises) and (ii) records each read.  :meth:`SiC4H.traceability`
returns the machine-generated table (name, value, unit, source DOI or owner approval) of
criterion S1.  The only numbers that do not come from the registry are

* documented physical constants (``K_B``, ``Q``, ``H_PLANCK``, ``M0``, ``EPS0``; SI 2019 /
  CODATA 2018), and
* the two spec-fixed choices ``T = 300 K`` and ``N_0 = 1e16 cm^-3`` (spec section 1).

Both groups appear in the table with source ``"constant"`` / ``"spec"``.

Unit conversions (explicit; each is covered by a test).  With ``V_T = k_B T/q`` [V], the
Debye length ``L_D = sqrt(eps_0 eps_r V_T / (q N_0))`` [cm] (``eps_r = eps_r_parallel_c``) and
the reference mobility ``mu_0`` [cm^2/(V s)]:

    length      x_s     = x[cm] / L_D
    potential   psi     = V / V_T
    field       E_s     = E[V/cm] * L_D / V_T            (unit V_T / L_D)
    alpha       alpha_s = alpha[1/cm] * L_D              (unit 1 / L_D)
    mobility    mu_s    = mu[cm^2/(V s)] / mu_0
    velocity    v_s     = v[cm/s] / (mu_0 V_T / L_D)     (so mu_s E_s / v_s = mu E / v)
    density     n_s     = n[cm^-3] / N_0
    sheet       Q_s     = Q[cm^-2] / (N_0 L_D)           (per elementary charge)
    energy      eps     = eps[eV] / V_T                  (E in units of k_B T)

Flat band of the n+ poly-Si gate (derivation, code ``V_FB`` convention of the stage-8.1 MOS
solver ``psi_gate = V_G - V_FB``; psi is measured so that ``n = n_i exp(psi - phi_n)``,
i.e. ``E_i(x) = E_i0 - q psi(x)``):

1. Vacuum level continuity: ``E_vac(x) = E_c0 + chi_SiC - q psi(x)`` with
   ``E_c0 = E_i0 + (E_c - E_i)_0`` (psi = 0 level), valid in the semiconductor and (no
   dipole assumed) in the oxide up to the gate.
2. Gate: ``E_F,gate = E_F,semi - q V_G`` (source/body at 0 V, ``E_F,semi = 0``) and
   ``E_vac,gate = E_F,gate + Phi_gate``.  Hence the gate potential in the psi scale is
   ``psi_gate = [E_c0 + chi_SiC - Phi_gate + q V_G]/q``.
3. Body bulk (psi_b, E_F = 0): ``E_c0 = q psi_b + (E_c - E_F)_bulk``.
4. n+ poly-Si: ``E_F,gate = E_c(Si)`` (degenerate), so ``Phi_gate - chi_SiC =
   (chi_Si - chi_ox) - (chi_SiC - chi_ox) = Phi_C,Si - Phi_C,SiC`` (conduction-band
   barriers to SiO2; Afanas'ev et al. 1996: 3.15 and 2.7 eV).
5. Therefore ``psi_gate = V_G + psi_b - V_FB`` with

       V_FB = (Phi_C,Si - Phi_C,SiC)/q - (E_c - E_F)_bulk/q            (spec formula),
       (E_c - E_F)_bulk = k_B T ln(N_c / n_b),   n_b = n_i exp(psi_b),

   where ``psi_b`` is the incomplete-ionization neutral root (stage 7/9a) and ``n_i``,
   ``N_c`` the stage-7 statistics.  The solver convention is
   ``psi_gate = V_G - V_FB,code`` with ``V_FB,code = V_FB - psi_b V_T``; combining with
   ``(E_c - E_F)_bulk = k_B T (ln(N_c/n_i) - psi_b)`` gives the *doping-independent*
   ``V_FB,code = [(Phi_C,Si - Phi_C,SiC) - k_B T ln(N_c/n_i)]/q`` (scaled by 1/V_T),
   because ``(E_c-E_i)_0 = k_B T ln(N_c/n_i)`` fixes the psi = 0 reference.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.special import expit

from widegap.devices.dmosfet import DMOSFET, DMOSParams
from widegap.fv.avalanche import ionization_integrals
from widegap.fv.drift_diffusion import DDParams
from widegap.fv.equilibrium import control_volumes
from widegap.fv.reverse_bias import RBReport
from widegap.fv2d.mos import Dirichlet
from widegap.fv2d.physics import (
    InterfaceCharge,
    Ionization,
    Mobility,
    PhysicsMOSDD,
    ct_factor,
)
from widegap.materials.ionization import (
    K_B,
    Level,
    Q,
    ionized_acceptors,
    ionized_donors,
    log_gamma_acceptor,
    log_gamma_donor,
    neutral_potential,
)
from widegap.materials.registry import Registry

H_PLANCK = 6.62607015e-34  # J s, exact (SI 2019)
M0 = 9.1093837015e-31  # kg, CODATA 2018
EPS0 = 8.8541878128e-12  # F/m, CODATA 2018 (vacuum permittivity)
T_K = 300.0  # spec section 1
N0_CM3 = 1.0e16  # spec section 1
EPS_MACH = np.finfo(float).eps
REGISTRY_PATH = Path(__file__).resolve().parent / "materials" / "data" / "sic4h.json"


# ---------------------------------------------------------------- traceability
class TracedRegistry:
    """Registry wrapper that records every parameter read (criterion S1)."""

    def __init__(self, reg: Registry):
        self.reg = reg
        self.used: dict[str, object] = {}

    @classmethod
    def from_json(cls, path=REGISTRY_PATH):
        return cls(Registry.from_json(path))

    def value(self, name: str) -> float:
        p = self.reg.get(name)  # raises ProvenanceError / KeyError: no fallback
        self.used[name] = p
        return p.value

    def usable(self, name: str) -> bool:
        return self.reg.usable(name)

    def rows(self) -> list[dict]:
        out = []
        for name in sorted(self.used):
            p = self.used[name]
            rec = p.record
            if p.param_class == "owner-approved-convention":
                src = f"owner approval: {rec['approved_by']} on {rec['approved_on']}"
            else:
                src = f"doi:{p.doi}"
            out.append(
                {
                    "name": name,
                    "value": p.value,
                    "unit": p.unit,
                    "source": src,
                    "class": p.param_class,
                    "locator": p.locator,
                    "orientation": p.orientation,
                    "registry_usable": True,
                }
            )
        return out


CONSTANT_ROWS = [
    ("K_B", K_B, "J/K", "constant: SI 2019 (exact)"),
    ("Q", Q, "C", "constant: SI 2019 (exact)"),
    ("H_PLANCK", H_PLANCK, "J s", "constant: SI 2019 (exact)"),
    ("M0", M0, "kg", "constant: CODATA 2018"),
    ("EPS0", EPS0, "F/m", "constant: CODATA 2018"),
    ("T", T_K, "K", "spec: stage12 section 1"),
    ("N0", N0_CM3, "cm^-3", "spec: stage12 section 1"),
]


# ---------------------------------------------------------------- statistics (stage 7)
def _nc_parabolic(m_rel, t_k):
    """2 (2 pi m k T / h^2)^{3/2} in cm^-3 for one parabolic valley (spin included)."""
    return 2.0 * (2.0 * np.pi * m_rel * M0 * K_B * t_k / H_PLANCK**2) ** 1.5 * 1e-6


def band_gap(reg, t_k: float) -> float:
    """E_g(T) = E_GX(T) + E_ex, simplified Passler model (stage 7)."""
    e0, eps, th, pp = (reg.value(n) for n in ("EGX_0", "passler_eps", "passler_theta", "passler_p"))
    egx = e0 - 0.5 * eps * th * ((1.0 + (2.0 * t_k / th) ** pp) ** (1.0 / pp) - 1.0)
    return float(egx + reg.value("E_exciton"))


def derived_statistics(reg, t_k: float):
    """(n_i, N_c, N_v) in cm^-3, stage-7 construction (Persson & Lindefelt DOS masses,
    parabolic bands; second CB minimum Boltzmann-weighted, three valence bands)."""
    kt = K_B * t_k / Q
    nc = reg.value("nu_c_band1") * _nc_parabolic(reg.value("m_c_band1"), t_k)
    nc += (
        reg.value("nu_c_band2")
        * _nc_parabolic(reg.value("m_c_band2"), t_k)
        * np.exp(-reg.value("Ec_band2_offset") / kt)
    )
    nv = _nc_parabolic(reg.value("m_v_band1"), t_k)
    nv += _nc_parabolic(reg.value("m_v_band2"), t_k) * np.exp(reg.value("Ev_band2") / kt)
    nv += _nc_parabolic(reg.value("m_v_band3"), t_k) * np.exp(reg.value("Ev_band3") / kt)
    eg = band_gap(reg, t_k)
    ni = np.sqrt(nc * nv) * np.exp(-eg / (2 * kt))
    return float(ni), float(nc), float(nv)


# ---------------------------------------------------------------- impact ionization
@dataclass(frozen=True)
class OktoCrowell:
    r"""Niwa (2015) Okuto-Crowell coefficients in scaled units (Table II, Eqs. 8-9):

        alpha(E) = a (1 + c (T-T0)) E^gamma exp{-[b (1 + d (T-T0)) / E]^delta},

    E in V/cm, alpha in 1/cm.  The registry carries ``a`` in cm^-1 and no ``gamma`` entry; the
    dimensionally consistent reading ``gamma = 0`` is used (flagged in R3).  Scaled
    evaluation: E_s -> E = E_s V_T/L_D, alpha_s = alpha L_D.  Interface-compatible with
    :class:`widegap.fv.avalanche.Chynoweth` (``alpha_n/p``, ``dalpha_n/p``)."""

    an: float
    bn: float
    dn: float
    ap: float
    bp: float
    dp: float
    e_unit: float  # V/cm per scaled field unit  = V_T / L_D[cm]
    l_d: float  # cm per scaled length unit

    def _alpha(self, e, a, b, d):
        ae = np.abs(np.asarray(e, dtype=float)) * self.e_unit  # V/cm
        out = np.zeros_like(ae)
        m = ae > 0.0
        with np.errstate(over="ignore", under="ignore"):
            out[m] = a * np.exp(-((b / ae[m]) ** d))
        return out * self.l_d

    def _dalpha(self, e, a, b, d):
        """d alpha_s / d|E_s| = alpha_s * delta (b/E)^delta / |E_s|."""
        ae = np.abs(np.asarray(e, dtype=float))
        al = self._alpha(ae, a, b, d)
        out = np.zeros_like(ae)
        m = (al > 0.0) & (ae > 0.0)
        with np.errstate(over="ignore", under="ignore"):
            out[m] = al[m] * d * (b / (ae[m] * self.e_unit)) ** d / ae[m]
        return out

    def alpha_n(self, e):
        return self._alpha(e, self.an, self.bn, self.dn)

    def alpha_p(self, e):
        return self._alpha(e, self.ap, self.bp, self.dp)

    def dalpha_n(self, e):
        return self._dalpha(e, self.an, self.bn, self.dn)

    def dalpha_p(self, e):
        return self._dalpha(e, self.ap, self.bp, self.dp)


@dataclass(frozen=True)
class ChynowethPhys:
    """Hatakeyama (2004) alpha = a exp(-b/E) with a [1/cm], b [V/cm]; scaled evaluation.
    Same interface as :class:`OktoCrowell` (and ``avalanche.Chynoweth``)."""

    an: float
    bn: float
    ap: float
    bp: float
    e_unit: float
    l_d: float

    def _alpha(self, e, a, b):
        ae = np.abs(np.asarray(e, dtype=float)) * self.e_unit
        out = np.zeros_like(ae)
        m = ae > 0.0
        with np.errstate(over="ignore", under="ignore"):
            out[m] = a * np.exp(-b / ae[m])
        return out * self.l_d

    def _dalpha(self, e, a, b):
        ae = np.abs(np.asarray(e, dtype=float))
        al = self._alpha(ae, a, b)
        out = np.zeros_like(ae)
        m = (al > 0.0) & (ae > 0.0)
        out[m] = al[m] * (b / self.e_unit) / ae[m] ** 2
        return out

    def alpha_n(self, e):
        return self._alpha(e, self.an, self.bn)

    def alpha_p(self, e):
        return self._alpha(e, self.ap, self.bp)

    def dalpha_n(self, e):
        return self._dalpha(e, self.an, self.bn)

    def dalpha_p(self, e):
        return self._dalpha(e, self.ap, self.bp)


# ---------------------------------------------------------------- the mapping
@dataclass
class SiC4H:
    """Physical -> scaled mapping for 4H-SiC built only from the registry (spec section 1)."""

    reg: TracedRegistry
    t_k: float = T_K
    n0: float = N0_CM3
    vsat_set: str = "khan_cooper"  # saturation-velocity source for the DD mobility (stage 15f)
    vsat_beta: float | None = None  # explicit Caughey-Thomas exponent (required for "ishikawa")
    vt: float = field(init=False)  # V
    eps_r: float = field(init=False)
    l_d: float = field(init=False)  # cm
    n_i: float = field(init=False)  # cm^-3
    n_c: float = field(init=False)
    n_v: float = field(init=False)
    e_g: float = field(init=False)  # eV
    delta: float = field(init=False)
    mu0: float = field(init=False)  # cm^2/(V s), reference mobility

    def __post_init__(self):
        r = self.reg
        self.vt = K_B * self.t_k / Q
        self.eps_r = r.value("eps_r_parallel_c")  # vertical (E || c) devices
        eps_cm = EPS0 * self.eps_r * 1e-2  # F/cm
        self.l_d = float(np.sqrt(eps_cm * self.vt / (Q * self.n0)))
        self.n_i, self.n_c, self.n_v = derived_statistics(r, self.t_k)
        self.e_g = band_gap(r, self.t_k)
        self.delta = self.n_i / self.n0
        # reference mobility: the largest registered low-field electron mobility (perp c)
        self.mu0 = r.value("mu_n_max_perp_c")

    # ----- unit conversions ------------------------------------------------
    @property
    def nu_c(self):
        return self.n_c / self.n0

    @property
    def nu_v(self):
        return self.n_v / self.n0

    @property
    def e_unit(self):
        """V/cm per scaled field unit (V_T / L_D)."""
        return self.vt / self.l_d

    def length_s(self, x_cm):
        return np.asarray(x_cm, dtype=float) / self.l_d

    def length_cm(self, x_s):
        return np.asarray(x_s, dtype=float) * self.l_d

    def volt_s(self, v):
        return np.asarray(v, dtype=float) / self.vt

    def volt_phys(self, psi):
        return np.asarray(psi, dtype=float) * self.vt

    def field_s(self, e_v_per_cm):
        return np.asarray(e_v_per_cm, dtype=float) * self.l_d / self.vt

    def field_phys(self, e_s):
        return np.asarray(e_s, dtype=float) * self.e_unit

    def alpha_s(self, alpha_per_cm):
        return np.asarray(alpha_per_cm, dtype=float) * self.l_d

    def mu_s(self, mu_cm2):
        return np.asarray(mu_cm2, dtype=float) / self.mu0

    def saturation(self):
        """(v_sat [cm/s], Caughey-Thomas exponent) of the selected set ``self.vsat_set``.

        * ``"khan_cooper"`` (default, unchanged since stage 12): Khan & Cooper, perp c, one doping
          (``vsat_n_perp_c``, ``ct_beta_n_perp_c``).
        * ``"ishikawa"`` (stage 15f, Amendment 1): isotropic ``v_sat = 1.7e7`` cm/s (Ishikawa thesis
          Table 3.1; ``vsat_n_perp_c_ishikawa``, ``vsat_n_par_c`` must agree) with the exponent
          given explicitly as ``vsat_beta``.  The thesis fit exponents at N_D <= 1e16 are < 1, which
          the stage-9a Caughey-Thomas implementation does not support (beta >= 1), so there is no
          implicit default: ``vsat_beta`` is a MODEL CHOICE, not a registry value.  A missing or
          unverified v_sat entry raises (no fallback).

        The set applies to electrons and, as before, to holes (``vsat_p = vsat_n``)."""
        r = self.reg
        if self.vsat_set == "khan_cooper":
            return r.value("vsat_n_perp_c"), r.value("ct_beta_n_perp_c")
        if self.vsat_set == "ishikawa":
            v_perp, v_par = r.value("vsat_n_perp_c_ishikawa"), r.value("vsat_n_par_c")
            if v_perp != v_par:
                raise ValueError(f"Ishikawa v_sat is not isotropic: perp {v_perp}, par {v_par}")
            if self.vsat_beta is None:
                raise ValueError(
                    "vsat_set='ishikawa' needs an explicit vsat_beta (a model choice >= 1; the "
                    "thesis alpha < 1 is outside the model, see stage-15 Amendment 1)"
                )
            return v_perp, float(self.vsat_beta)
        raise ValueError(f"unknown vsat_set {self.vsat_set!r} (use 'khan_cooper' or 'ishikawa')")

    def vsat_s(self, v_cm_s):
        return np.asarray(v_cm_s, dtype=float) / (self.mu0 * self.vt / self.l_d)

    def dens_s(self, n_cm3):
        return np.asarray(n_cm3, dtype=float) / self.n0

    def sheet_s(self, q_cm2):
        return np.asarray(q_cm2, dtype=float) / (self.n0 * self.l_d)

    def energy_s(self, e_ev):
        return np.asarray(e_ev, dtype=float) / self.vt

    def current_unit_a_per_cm2(self, dim=2):
        """q N0 mu0 V_T / L_D * L_D^(dim-1)... returned per unit (scaled) contact length:
        J0 = q N0 mu0 V_T/L_D [A/cm^2]; a 2D scaled current multiplies J0 L_D [A/cm]."""
        return Q * self.n0 * self.mu0 * self.vt / self.l_d

    # ----- dopant levels (stage 7) ----------------------------------------------
    def donor_levels(self):
        """Götz two-site N donor (h / k), registered site fraction and degeneracy."""
        r = self.reg
        g = r.value("g_donor")
        fh = r.value("N_site_fraction_h")
        return (
            Level(fh, log_gamma_donor(g, r.value("dE_N_donor_h"), self.t_k, self.delta, self.nu_c)),
            Level(
                1 - fh,
                log_gamma_donor(g, r.value("dE_N_donor_k"), self.t_k, self.delta, self.nu_c),
            ),
        )

    def acceptor_energy_ev(self, na_cm3):
        """Koizumi (2009) Eq. (7): dE_A = dE_A0 - alpha N_A^(1/3)."""
        return self.reg.value("dE_Al_acceptor") - self.reg.value("dE_Al_acceptor_alpha") * (
            float(na_cm3) ** (1.0 / 3.0)
        )

    def acceptor_levels(self, na_cm3):
        de = self.acceptor_energy_ev(na_cm3)
        return (
            Level(
                1.0,
                log_gamma_acceptor(
                    self.reg.value("g_acceptor"), de, self.t_k, self.delta, self.nu_v
                ),
            ),
        )

    # ----- mobility (Ishikawa Caughey-Thomas, Hall; flagged) ----------------------
    def mobility_cm2(self, n_total_cm3, carrier="n", orient="perp"):
        """mu = (mu_max-mu_min)(T/300)^d1 / [1 + (N/N_ref)^a (T/300)^d2] + mu_min (T/300)^d3."""
        r = self.reg
        c, o = carrier, ("par" if orient == "par" else "perp")
        sfx = f"{o}_c"
        mx = r.value(f"mu_{c}_max_{sfx}")
        mn = r.value(f"mu_{c}_min_{sfx}")
        nref = r.value(f"N_{c}_ref_{sfx}")
        al = r.value(f"alpha_{c}_{sfx}")
        t = self.t_k / 300.0
        if c == "n":
            d1, d2, d3 = (r.value(f"delta{k}_n_{sfx}") for k in (1, 2, 3))
        else:  # hole fits: thesis deltas are in the registry validity text only; at T = 300 K
            d1 = d2 = d3 = 0.0  # all (T/300)^d factors equal 1 (t == 1 below), so d_i is moot
            if abs(t - 1.0) > 1e-12:
                raise ValueError("hole mobility temperature exponents are not registry values")
        n = np.asarray(n_total_cm3, dtype=float)
        return (mx - mn) * t**d1 / (1.0 + (n / nref) ** al * t**d2) + mn * t**d3

    # ----- coefficient sets ---------------------------------------------------------
    def niwa(self) -> OktoCrowell:
        r = self.reg
        dt = self.t_k - r.value("oc_T0")
        fa_n = 1.0 + r.value("oc_c_n") * dt
        fb_n = 1.0 + r.value("oc_d_n") * dt
        fa_p = 1.0 + r.value("oc_c_p") * dt
        fb_p = 1.0 + r.value("oc_d_p") * dt
        return OktoCrowell(
            r.value("oc_a_n") * fa_n,
            r.value("oc_b_n") * fb_n,
            r.value("oc_delta_n"),
            r.value("oc_a_p") * fa_p,
            r.value("oc_b_p") * fb_p,
            r.value("oc_delta_p"),
            self.e_unit,
            self.l_d,
        )

    def hatakeyama(self) -> ChynowethPhys:
        r = self.reg
        return ChynowethPhys(
            r.value("ii_a_n_par_c"),
            r.value("ii_b_n_par_c"),
            r.value("ii_a_p_par_c"),
            r.value("ii_b_p_par_c"),
            self.e_unit,
            self.l_d,
        )

    # ----- flat band ----------------------------------------------------------------
    def vfb(self, na_body_cm3):
        """(V_FB [V] spec formula, V_FB,code [V_T], psi_b) for an Al-doped body (module doc)."""
        r = self.reg
        dphi = r.value("PhiC_Si_SiO2") - r.value("PhiC_4HSiC_SiO2")  # eV
        lev = self.acceptor_levels(na_body_cm3)
        psi_b = neutral_potential(0.0, na_body_cm3 / self.n0, self.delta, (), lev)
        n_b = self.n_i * np.exp(psi_b)
        ec_ef = self.vt * np.log(self.n_c / n_b)  # (E_c - E_F)_bulk / q  [V]
        v_fb = dphi - ec_ef
        v_code = (dphi - self.vt * np.log(self.n_c / self.n_i)) / self.vt
        return float(v_fb), float(v_code), float(psi_b)

    def traceability(self) -> list[dict]:
        rows = self.reg.rows()
        for name, val, unit, src in CONSTANT_ROWS:
            rows.append(
                {
                    "name": name,
                    "value": val,
                    "unit": unit,
                    "source": src,
                    "class": "constant" if src.startswith("constant") else "spec",
                    "locator": None,
                    "orientation": None,
                    "registry_usable": None,
                }
            )
        rows.append(
            {
                "name": "gamma_OktoCrowell",
                "value": 0.0,
                "unit": "1",
                "source": "structural: no registry entry; alpha 'a' is in 1/cm so E^gamma "
                "must be absent (flagged R3)",
                "class": "structural",
                "locator": None,
                "orientation": None,
                "registry_usable": None,
            }
        )
        return rows


# ---------------------------------------------------------------- 1D reverse-bias diode
class IonizedReverseDiode:
    r"""Zero-current reverse-bias electrostatics with incomplete ionization (1D).

    Flat quasi-Fermi levels (stage 4): ``phi_n = 0`` (cathode), ``phi_p = -V_R`` (anode),

        -psi'' = delta (e^{phi_p - psi} - e^{psi - phi_n})
                 + N_D^+(psi - phi_n) - N_A^-(psi - phi_p),

    with the stage-7 ionized fractions evaluated at ``psi - phi`` (so at ``V_R = 0`` this is
    exactly :class:`widegap.fv.ionized_equilibrium.IonizedJunction`).  The contact potentials
    are the incomplete-ionization neutral roots.  Interface-compatible with
    :class:`widegap.fv.reverse_bias.Diode1D` (``x``, ``solve``, ``element_field``) so
    :func:`widegap.fv.avalanche.avalanche_bv` applies unchanged.

    Newton: box-method FV, analytic Jacobian, step limiter ``max_step``; convergence
    requires the residual to be below its round-off tolerance and the last step below
    ``1e-8 + 1e3 eps max|psi|`` (round-off floor of psi itself, |psi| up to 1e5 V_T at kV).
    The initial guess is the abrupt-depletion profile at the requested bias (valid for a
    p+ / n-drift structure), so large bias jumps need no continuation.
    """

    def __init__(self, x, nd, na, delta, donors, acceptors):
        self.x = np.asarray(x, dtype=float)
        self.nd = np.asarray(nd, dtype=float)
        self.na = np.asarray(na, dtype=float)
        self.delta = float(delta)
        self.donors, self.acceptors = tuple(donors), tuple(acceptors)
        self.vol = control_volumes(self.x)
        self.g = 1.0 / np.diff(self.x)
        self.psi_anode_eq = float(
            neutral_potential(self.nd[0], self.na[0], delta, self.donors, self.acceptors)
        )
        self.psi_cathode_eq = float(
            neutral_potential(self.nd[-1], self.na[-1], delta, self.donors, self.acceptors)
        )
        self._guess_cache: dict = {}

    # charge as a function of psi at flat quasi-Fermi levels (phi_p, phi_n)
    def _rho(self, psi, phi_n, phi_p):
        return (
            self.delta * (np.exp(phi_p - psi) - np.exp(psi - phi_n))
            + ionized_donors(psi - phi_n, self.nd, self.donors)
            - ionized_acceptors(psi - phi_p, self.na, self.acceptors)
        )

    def _drho(self, psi, phi_n, phi_p):
        return (
            -self.delta * (np.exp(phi_p - psi) + np.exp(psi - phi_n))
            + self._d_ion(psi - phi_n, self.nd, self.donors, +1)
            + self._d_ion(psi - phi_p, self.na, self.acceptors, -1)
        )

    @staticmethod
    def _d_ion(z, dens, levels, sign):
        """d/dz of N_D^+(z) (sign=+1) or of -N_A^-(z) (sign=-1) wrt z = psi - phi."""
        out = np.zeros_like(z)
        for lv in levels:
            if sign > 0:  # N_D^+ = f nd / (1 + exp(z + l))
                s = expit(z + lv.log_gamma)
                out = out - lv.fraction * dens * s * (1.0 - s)
            else:  # -N_A^- = -f na / (1 + exp(-z + l))
                t = expit(lv.log_gamma - z)
                out = out - lv.fraction * dens * t * (1.0 - t)
        return out

    def residual(self, psi, v_r):
        flux = -np.diff(psi) * self.g
        f = np.zeros_like(psi)
        rho = self._rho(psi, 0.0, -v_r)
        f[1:-1] = (flux[1:] - flux[:-1]) - self.vol[1:-1] * rho[1:-1]
        return f

    def jacobian_banded(self, psi, v_r):
        dr = self._drho(psi, 0.0, -v_r)
        ab = np.zeros((3, len(psi)))
        ab[1, 1:-1] = self.g[:-1] + self.g[1:] - self.vol[1:-1] * dr[1:-1]
        ab[0, 2:] = -self.g[1:]
        ab[2, :-2] = -self.g[:-1]
        ab[1, 0] = ab[1, -1] = 1.0
        return ab

    def tolerance(self, psi, v_r):
        ab = self.jacobian_banded(psi, v_r)
        a = np.abs(psi)
        s = np.abs(ab[1]) * a
        s[:-1] += np.abs(ab[0, 1:]) * a[1:]
        s[1:] += np.abs(ab[2, :-1]) * a[:-1]
        mag = self.delta * (np.exp(-v_r - psi) + np.exp(psi)) + np.abs(self.nd) + np.abs(self.na)
        return 1e2 * EPS_MACH * float(np.max((s + self.vol * mag)[1:-1]))

    def depletion_guess(self, v_r):
        """Abrupt-depletion guess for a p+ (x < x_j) / n-drift structure, junction at the
        first node with nd > na.  psi = psi_p - V on the p side, psi_n on the n side, parabola
        in between from the drift charge ``N_D`` (complete ionization estimate)."""
        x = self.x
        j = int(np.argmax(self.nd > self.na))
        xj = x[j]
        psi_p = self.psi_anode_eq - v_r
        psi_n = self.psi_cathode_eq
        nd_drift = float(self.nd[min(j + 1, len(x) - 1)])
        na_p = float(self.na[0])
        drop = psi_n - psi_p  # total potential drop to be taken up
        # one-sided width, small p-side width by charge balance
        w_n = np.sqrt(2.0 * drop / (nd_drift * (1.0 + nd_drift / na_p)))
        w_p = w_n * nd_drift / na_p
        psi = np.where(x < xj, psi_p, psi_n)
        xn = (x - xj) / w_n
        inside_n = (x >= xj) & (xn < 1.0)
        psi = np.where(inside_n, psi_n - 0.5 * nd_drift * (w_n - (x - xj)) ** 2, psi)
        inside_p = (x < xj) & (x > xj - w_p)
        psi = np.where(inside_p, psi_p + 0.5 * na_p * (x - (xj - w_p)) ** 2, psi)
        return psi

    def solve(self, v_r, psi0=None, max_iter=200, max_step=2.0):
        if psi0 is not None and abs((self.psi_anode_eq - psi0[0]) - v_r) <= 5.0:
            psi = np.array(psi0, dtype=float)  # warm start only for a small bias change
        else:
            psi = self.depletion_guess(v_r)
        psi[0] = self.psi_anode_eq - v_r
        psi[-1] = self.psi_cathode_eq
        rep = RBReport()
        last_step = np.inf
        for _ in range(max_iter + 1):
            f = self.residual(psi, v_r)
            rn = float(np.max(np.abs(f)))
            tol = self.tolerance(psi, v_r)
            rep.residual_history.append(rn)
            rep.tolerance = tol
            if rn <= tol and last_step <= 1e-8 + 1e3 * EPS_MACH * float(np.max(np.abs(psi))):
                rep.converged = True
                break
            if rep.iterations == max_iter:
                break
            dpsi = solve_banded((1, 1), self.jacobian_banded(psi, v_r), -f)
            step = float(np.max(np.abs(dpsi)))
            if step > max_step:
                dpsi *= max_step / step
                rep.limited_steps += 1
            last_step = float(np.max(np.abs(dpsi)))
            psi = psi + dpsi
            rep.iterations += 1
        return psi, rep

    def element_field(self, psi):
        return -np.diff(psi) / np.diff(self.x)

    def max_field(self, psi):
        return float(np.max(np.abs(self.element_field(psi))))


# ---------------------------------------------------------------- 12a design helpers
def graded_mesh(lo, hi, h0, growth, hmax):
    """Nodes on [lo, hi]: first cell h0, geometric growth ratio ``growth`` up to ``hmax``."""
    x = [lo]
    h = h0
    while x[-1] < hi - 1e-12 * (hi - lo):
        x.append(min(x[-1] + h, hi))
        h = min(h * growth, hmax)
    return np.array(x)


def refine_bisect(x, level):
    x = np.asarray(x, dtype=float)
    for _ in range(level):
        x = np.sort(np.r_[x, 0.5 * (x[:-1] + x[1:])])
    return x


def triangular_breakdown(nd_s, ch, n_pts=3001):
    """Design estimate (not a result): abrupt one-sided junction, triangular field
    E(x) = E_m (1 - x/W_d), W_d = E_m / N_D (scaled; dE/dx = N_D/N_0).  Returns
    (E_m, W_d, V_bd) in scaled units with I_n(E_m) = 1 for coefficient set ``ch``."""

    def integral(e_m):
        w = e_m / nd_s
        x = np.linspace(0.0, w, n_pts)
        xm = 0.5 * (x[1:] + x[:-1])
        with np.errstate(all="ignore"):
            val = ionization_integrals(x, e_m * (1.0 - xm / w), ch)[0]
        return val if np.isfinite(val) else np.inf

    lo, hi = 10.0, 1500.0  # scaled field bracket: ~7e4 .. 1e7 V/cm
    if integral(hi) < 1.0:
        raise ValueError("no breakdown field below the search ceiling")
    e_m = float(brentq(lambda e: integral(e) - 1.0, lo, hi, xtol=1e-12, rtol=1e-12))
    w_d = e_m / nd_s
    return e_m, w_d, 0.5 * e_m * w_d


def build_pn_diode(sc: SiC4H, nd_cm, w_d_s, na_cm=1e19, level=0, rule=None):
    """p+/n 1D diode (anode p+ for x < 0, junction at x = 0, n drift of width ``w_d_s`` scaled).

    Mesh rule (fixed up front, all lengths in L_D; ``rule`` overrides for tests):
    junction spacing  h0 = w_p / 40  with w_p = W_dep N_D / N_A the p+ depletion width at BV;
    drift spacing capped at W_dep/120 (``w_d_s`` here is the non-punch-through *width*; the
    cap uses its depletion fraction 1/rule['width_factor']); geometric growth 1.08;
    p+ contact region 2 L_D thick (cells grow up to 1 L_D).  ``level`` bisections.
    """
    r = {"w_dep": w_d_s / 1.25, "growth": 1.08, "h0_div": 40.0, "hmax_div": 120.0, "lp": 2.0}
    if rule:
        r.update(rule)
    w_dep = r["w_dep"]
    w_p = w_dep * nd_cm / na_cm
    h0 = w_p / r["h0_div"]
    xr = graded_mesh(0.0, w_d_s, h0, r["growth"], w_dep / r["hmax_div"])
    xl = -graded_mesh(0.0, r["lp"], h0, r["growth"], 1.0)[::-1]
    x = refine_bisect(np.r_[xl[:-1], xr], level)
    nd = np.where(x > 0, sc.dens_s(nd_cm), 0.0)
    na = np.where(x < 0, sc.dens_s(na_cm), 0.0)
    j = int(np.argmin(np.abs(x)))  # junction node: half donor, half acceptor (CV average)
    nd[j], na[j] = 0.5 * sc.dens_s(nd_cm), 0.5 * sc.dens_s(na_cm)
    diode = IonizedReverseDiode(x, nd, na, sc.delta, sc.donor_levels(), sc.acceptor_levels(na_cm))
    return diode, {"h0": h0, "w_p": w_p, "w_dep": w_dep, "n_nodes": len(x)}


# ---------------------------------------------------------------- 2D vertical DMOSFET
class AnisoMobilityDD(PhysicsMOSDD):
    """Stage-9a physics with a fixed per-edge low-field electron mobility factor.

    ``mu_edge[e] = mu_low,e / mu_n`` (``prm.mu_n`` is the reference mobility, scaled 1).  On a
    tensor mesh each edge is axis-aligned, so an edge along x (E perp c, basal plane) can take
    the perp-c mobility and an edge along y (E || c) the par-c mobility exactly.  The
    Caughey-Thomas factor uses the edge's own low-field mobility,
    ``m_e = [1 + (mu_e F / v_sat)^beta]^(-1/beta)``, and the flux is ``mu_edge * m_e * j_0``;
    ``mu_edge`` is a fixed constant, so the analytic Jacobian of the stage-9a class stays exact
    (both ``m`` and ``dm/dphi`` are multiplied by it)."""

    def __init__(self, *args, mu_edge, **kw):
        super().__init__(*args, **kw)
        self.mu_edge = np.asarray(mu_edge, dtype=float)

    def _mob_factors(self, X):
        _, pn, pp = self.split(X)
        a, b, h = self.a, self.b, self.mesh.h
        mob, prm = self.mob, self.prm
        mn, kn = ct_factor(pn[b] - pn[a], h, prm.mu_n * self.mu_edge, mob.vsat_n, mob.beta_n)
        mp, kp = ct_factor(pp[b] - pp[a], h, prm.mu_p, mob.vsat_p, mob.beta_p)
        return self.mu_edge * mn, mp, self.mu_edge * kn, kp


@dataclass
class DMOSDesign:
    """Physical dimensions [um] and dopings [cm^-3] of the half-cell (see run.py)."""

    x_jfet: float
    ch_len: float
    src_w: float
    bc_w: float
    pad_w: float
    d_src: float
    d_body: float
    w_drift: float
    t_drain: float
    t_ox: float
    n_drift: float
    n_jfet: float
    n_body: float
    n_src: float
    n_drain: float
    mesh_kw: dict = field(default_factory=dict)


class SiCDMOS:
    """Vertical 4H-SiC DMOSFET half-cell in registry-sourced scaled units (stage 12b)."""

    def __init__(
        self, sc: SiC4H, design: DMOSDesign, level=0, dit_scale=1.0, qf_cm2=0.0, dit_width_ev=0.2
    ):
        self.sc = sc
        self.design = d = design
        um = 1e-4
        ls = lambda v: float(sc.length_s(v * um))  # noqa: E731
        ds = lambda v: float(sc.dens_s(v))  # noqa: E731
        x_jfet = ls(d.x_jfet)
        x_src0 = x_jfet + ls(d.ch_len)
        x_src1 = x_src0 + ls(d.src_w)
        x_bc = x_src1 + ls(0.5 * d.bc_w)
        prm = DMOSParams(
            x_jfet=x_jfet, x_src0=x_src0, x_src1=x_src1, x_bc=x_bc, x_max=x_bc + ls(d.pad_w),
            gate_x1=x_src0, x_ch_mid=x_jfet + 0.5 * ls(d.ch_len),
            d_src=ls(d.d_src), d_body=ls(d.d_body), w_drift=ls(d.w_drift),
            t_drain=ls(d.t_drain), t_ox=ls(d.t_ox),
            a_body=ds(d.n_body), b_src=ds(d.n_src), b_drain=ds(d.n_drain),
            n0_drift=ds(d.n_drift), n_jfet=ds(d.n_jfet),
            eps_ratio=sc.reg.value("eps_r_SiO2_thermal") / sc.eps_r,
            delta=sc.delta, mu_n=1.0, mu_p=1.0, **d.mesh_kw,
        )  # fmt: skip
        self.prm = prm
        self.dev = dev = DMOSFET(prm, level)
        self.level = level
        q = dev.cell_doping(np.zeros(5))
        nd_c, na_c = np.maximum(q, 0.0), np.maximum(-q, 0.0)
        self.nd, self.na = dev.M @ nd_c, dev.M @ na_c
        self.donors = sc.donor_levels()
        self.acceptors = sc.acceptor_levels(d.n_body)
        mesh = dev.mesh
        # per-edge low-field electron mobility (Ishikawa, Hall; perp c along x, par c along y)
        # N = total ionizable impurity density (N_D + N_A): the registry fit is defined for the
        # total donor density of n-type material and is undefined for minority electrons in the
        # p-body channel (flagged in R3).
        nd_cm = 0.5 * (self.nd[mesh.ea] + self.nd[mesh.eb] + self.na[mesh.ea] + self.na[mesh.eb])
        nd_cm = nd_cm * sc.n0
        mu_perp = sc.mobility_cm2(nd_cm, "n", "perp")
        mu_par = sc.mobility_cm2(nd_cm, "n", "par")
        self.mu_edge = sc.mu_s(np.where(mesh.edge_dir == 0, mu_perp, mu_par))
        mu_p = float(sc.mu_s(sc.mobility_cm2(d.n_body, "p", "perp")))
        vsat_cm, beta = sc.saturation()
        vs = float(sc.vsat_s(vsat_cm))
        self.mob = Mobility(vsat_n=vs, vsat_p=vs, beta_n=beta, beta_p=beta)
        self.dp = DDParams(delta=sc.delta, mu_n=1.0, mu_p=mu_p, recombination=False)
        # interface: D_it uniform on [E_c - width, E_c] (above midgap, E_i = intrinsic level)
        ec_ei = sc.vt * np.log(sc.n_c / sc.n_i)  # eV
        self.e_c_minus_e_i_ev = float(ec_ei)
        dit_cm = sc.reg.value("Dit_Ec_minus_0p2eV_NO_0001") * dit_scale
        self.itf = InterfaceCharge(
            q_f=float(sc.sheet_s(qf_cm2)),
            dit=float(sc.sheet_s(dit_cm * sc.vt)),  # cm^-2 eV^-1 * eV per V_T
            e1=float(sc.energy_s(ec_ei - dit_width_ev)),
            e2=float(sc.energy_s(ec_ei)),
        )
        self.v_fb_phys, self.vfb_code, self.psi_body = sc.vfb(d.n_body)
        self.model = AnisoMobilityDD(
            dev.reg, self.nd - self.na, self.dp, float(np.log(sc.delta)),
            eps_ratio=prm.eps_ratio,
            ionization=Ionization(self.nd, self.na, self.donors, self.acceptors),
            mobility=self.mob, interface=self.itf, mu_edge=self.mu_edge,
        )  # fmt: skip
        self._blocking = dev.blocking_edge_mask()
        self._count_solves()

    def bc(self, vg, vs, vd) -> Dirichlet:
        n, mod = self.dev.nodes, self.model
        bc = Dirichlet()
        mod.ohmic(bc, n["source_body"], vs)
        mod.ohmic(bc, n["drain"], vd)
        bc.gate(n["gate"], vg, self.vfb_code)
        return bc

    # ----- off state (zero current, stage-4 electrostatics in 2D) -------------------------
    def offstate_bc(self, vd, vg=0.0) -> Dirichlet:
        """Zero-current reverse bias: all quasi-Fermi levels prescribed (flat, as stage 4).

        phi_n = V_D on the n-type nodes (drift, JFET, drain) and 0 on p-type nodes and the n+
        source box (those are connected to the source through the p-body / channel);
        phi_p = 0 on every semiconductor node (the contacts keep their ohmic values).  Minority
        carrier densities are < delta-level wherever the prescribed level differs from the
        true one, so the electrostatics equal those of the full drift-diffusion state at zero
        current.  The solver then only updates psi (phi rows are identity)."""
        dev, p = self.dev, self.prm
        bc = self.bc(vg, 0.0, vd)
        xy = dev.mesh.node_xy
        semi = dev.reg.semi_node
        in_src = (
            (xy[:, 0] >= p.x_src0 - 1e-9)
            & (xy[:, 0] <= p.x_src1 + 1e-9)
            & (xy[:, 1] >= -p.d_src - 1e-9)
        )
        n_type = self.nd > self.na
        contact = np.zeros(len(xy), dtype=bool)
        for k in ("source_body", "drain"):
            contact[dev.nodes[k]] = True
        free = np.flatnonzero(semi & ~contact)
        phi_n = np.where(n_type[free] & ~in_src[free], vd, 0.0)
        bc.phi(free, phi_n, 0.0)
        return bc

    # ----- solves, bookkeeping ---------------------------------------------------------
    def _count_solves(self):
        self.stats = {"attempts": 0, "failed_attempts": 0}
        orig = self.model.solve_bc

        def counted(x0, bc, **kw):
            sol = orig(x0, bc, **kw)
            self.stats["attempts"] += 1
            self.stats["failed_attempts"] += int(not sol.report.converged)
            return sol

        self.model.solve_bc = counted

    def equilibrium(self):
        bc0 = self.bc(0.0, 0.0, 0.0)
        mod = self.model
        return mod.solve_bc(mod.initial_guess(bc0), mod.equilibrium_bc(bc0))

    def currents(self, sol):
        """Terminal currents (scaled) of the source, drain and body contact node sets."""
        n, mod = self.dev.nodes, self.model
        return {
            "I_source": mod.contact_current(sol, n["source"]),
            "I_drain": mod.contact_current(sol, n["drain"]),
            "I_body": mod.contact_current(sol, n["body"]),
        }

    def current_a_per_cm(self, i_scaled):
        """Scaled 2D contact current -> A per cm of device depth: I_s * q N_0 mu_0 V_T
        (J_0 = q N_0 mu_0 V_T / L_D [A/cm^2] times one scaled length L_D)."""
        sc = self.sc
        return float(i_scaled) * Q * sc.n0 * sc.mu0 * sc.vt

    # ----- BV by the ionization integral along a vertical cut -----------------------------
    def peak_field(self, psi, x_min=None):
        """Largest |E| (scaled) over the blocking semiconductor edges (optionally only edges
        with midpoint x >= ``x_min``); returns (E, x, y, dir) of the edge midpoint."""
        mod, mesh = self.model, self.dev.mesh
        e = np.abs(psi[mod.b] - psi[mod.a]) / mesh.h
        keep = self._blocking
        if x_min is not None:
            keep = keep & (0.5 * (mesh.node_xy[mod.a, 0] + mesh.node_xy[mod.b, 0]) >= x_min)
        e = np.where(keep, e, 0.0)
        k = int(np.argmax(e))
        mid = 0.5 * (mesh.node_xy[mod.a[k]] + mesh.node_xy[mod.b[k]])
        return float(e[k]), float(mid[0]), float(mid[1]), int(mesh.edge_dir[k])

    def column_nodes(self, x_target):
        """Node column (x snapped to the mesh) and its semiconductor nodes, top (y = 0) first."""
        mesh = self.dev.mesh
        x_col = float(mesh.x[int(np.argmin(np.abs(mesh.x - x_target)))])
        xy = mesh.node_xy
        sel = np.flatnonzero((xy[:, 0] == x_col) & (xy[:, 1] <= 0.0))
        return x_col, sel[np.argsort(-xy[sel, 1])]

    def column_integrals(self, psi, nodes, ch):
        """(I_n, I_p, A_W) along the column, anode (p-body, y = 0) at index 0 as in stage 9b."""
        depth = -self.dev.mesh.node_xy[nodes, 1]
        e = -np.diff(psi[nodes]) / np.diff(depth)
        return ionization_integrals(depth, e, ch)

    def offstate_step(self, x, v, v_new, max_step=20.0):
        pred = tangent_predictor(self.model, x, self.offstate_bc, v, v_new - v)
        return self.model.solve_bc(pred, self.offstate_bc(v_new), max_step=max_step)

    def offstate_bv(self, ch, dv_s, v_max_s, x_eq, form="n", max_pass=4, x_min=None):
        """BV: root of I_form(V_D) = 1 on a FIXED vertical node column (continuity of I(V)).

        1. March V_D up (tangent predictor + Newton, step halving) until the integral along the
           column through the current peak-field point reaches 1; store all states.
        2. Fix the column at the peak of the bracketing state; Brent on the stored bracket
           with the same predictor/corrector.  If the peak at the root is on another column,
           repeat with that column (``max_pass``); report the final consistency.
        ``x_min`` restricts the peak-field search to edges with x >= x_min (planar cut).
        Returns a dict with the BV, the column, the peak-field point and diagnostics.
        """
        k = 0 if form == "n" else 1
        mod = self.model
        x0 = mod.solve_bc(x_eq, self.offstate_bc(0.0), max_step=20.0)
        if not x0.report.converged:
            raise RuntimeError("off-state zero-bias solve did not converge")
        states = [(0.0, x0.X.copy())]
        v, h = 0.0, dv_s
        col = None
        while True:
            v_try = min(v + h, v_max_s)
            s = self.offstate_step(states[-1][1], v, v_try)
            if not s.report.converged:
                h *= 0.5
                if h < 1e-3:
                    raise RuntimeError(f"off-state march failed at V_D = {v:.1f} V_T")
                continue
            v = v_try
            states.append((v, s.X.copy()))
            psi = s.X[0::3]
            _, xp, _, _ = self.peak_field(psi, x_min)
            col, nodes = self.column_nodes(xp)
            if self.column_integrals(psi, nodes, ch)[k] >= 1.0:
                break
            if v >= v_max_s:
                raise ValueError("I = 1 not reached below the voltage ceiling")
        history = []
        for _ in range(max_pass):
            _, nodes = self.column_nodes(col)
            vals = [self.column_integrals(x[0::3], nodes, ch)[k] for _, x in states]
            kk = next((i for i, val in enumerate(vals) if val >= 1.0), None)
            if kk is None or kk == 0:
                # the new column has not broken down yet within the stored states: extend
                raise RuntimeError("stored states do not bracket I = 1 on the chosen column")
            v_lo, x_lo = states[kk - 1]
            v_hi = states[kk][0]
            cache = {}

            def g(v_, x_lo=x_lo, v_lo=v_lo, nodes=nodes, cache=cache):
                sol = self.offstate_step(x_lo, v_lo, v_)
                if not sol.report.converged:
                    raise RuntimeError(f"Brent evaluation failed at V_D = {v_:.2f} V_T")
                cache[v_] = sol.X.copy()
                return self.column_integrals(sol.X[0::3], nodes, ch)[k] - 1.0

            bv = float(brentq(g, v_lo, v_hi, xtol=1e-9, rtol=4 * EPS_MACH))
            sol = self.offstate_step(x_lo, v_lo, bv)
            psi_bv = sol.X[0::3]
            pk = self.peak_field(psi_bv, x_min)
            col_new, _ = self.column_nodes(pk[1])
            i_n, i_p, a_w = self.column_integrals(psi_bv, nodes, ch)
            history.append({"x_col": col, "bv": bv, "peak_x": pk[1], "peak_col": col_new})
            # consistent: the peak is on this column, or this column's own peak field equals the
            # global peak to 1e-3 (a planar junction has a plateau of near-equal fields)
            e_col = float(
                np.max(np.abs(np.diff(psi_bv[nodes]) / np.diff(-self.dev.mesh.node_xy[nodes, 1])))
            )
            if col_new == col or e_col >= (1.0 - 1e-3) * pk[0]:
                return {
                    "bv": bv, "x_col": col, "form": form, "I_n": i_n, "I_p": i_p, "A_W": a_w,
                    "peak_E": pk[0], "peak_x": pk[1], "peak_y": pk[2], "peak_dir": pk[3],
                    "consistent": True, "passes": history, "X": sol.X.copy(),
                    "n_states": len(states),
                }  # fmt: skip
            col = col_new
        return {
            "bv": bv, "x_col": col, "form": form, "I_n": i_n, "I_p": i_p, "A_W": a_w,
            "peak_E": pk[0], "peak_x": pk[1], "peak_y": pk[2], "peak_dir": pk[3],
            "consistent": False, "passes": history, "X": sol.X.copy(), "n_states": len(states),
        }  # fmt: skip


def tangent_predictor(model, x, bc_of, v, dv):
    """First-order continuation predictor dX/dV = J^{-1} d(bc)/dV (Euler step of size ``dv``)."""
    return tangent_predictor_bc(model, x, bc_of(v), bc_of(v + 1.0), dv)


def tangent_predictor_bc(model, x, bc_a, bc_b, frac):
    """Euler predictor X + frac * J(X)^{-1} (g_b - g_a) for Dirichlet data moving linearly from
    ``bc_a`` to ``bc_b`` (frac = fraction of that move; row-equilibrated sparse solve)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve

    rows, va = model._all_bc(bc_a)
    rows_b, vb = model._all_bc(bc_b)
    if not np.array_equal(rows, rows_b):
        raise ValueError("Dirichlet row sets of the two bias points differ")
    rhs = np.zeros(len(x))
    rhs[rows] = vb - va
    jac = model.jacobian_bc(x, rows)
    rs = np.asarray(abs(jac).max(axis=1).todense()).ravel()
    return x + frac * spsolve(sp.csc_matrix(sp.diags(1.0 / rs) @ jac), rhs / rs)


__all__ = [
    "tangent_predictor",
    "tangent_predictor_bc",
    "AnisoMobilityDD",
    "DMOSDesign",
    "SiCDMOS",
    "IonizedReverseDiode",
    "OktoCrowell",
    "ChynowethPhys",
    "SiC4H",
    "TracedRegistry",
    "band_gap",
    "derived_statistics",
]
