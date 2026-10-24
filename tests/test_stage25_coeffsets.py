"""Fast tests for stage 25: coefficient-set switch (``ev.ch``) and ionization-set sanity."""

import numpy as np

from widegap.adjoint_sic import SiCAdjoint
from widegap.fv.avalanche import ionization_integrals
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip


def test_default_is_niwa_and_sets_differ():
    sc = SiC4H(TracedRegistry.from_json())
    ev = SiCAdjoint(sc, DESIGN, 0)
    assert ev.ch == sc.niwa()  # default unchanged (frozen dataclass equality)
    e = np.array([2.0, 2.5, 3.0]) / sc.e_unit * 1e6
    an_n, an_h = sc.niwa().alpha_n(e), sc.hatakeyama().alpha_n(e)
    assert np.all(an_n > 0) and np.all(an_h > 0) and not np.allclose(an_n, an_h)


def test_ionization_integral_responds_to_set():
    sc = SiC4H(TracedRegistry.from_json())
    x = np.linspace(0.0, 1.0, 40) * float(sc.length_s(1e-4))
    e = np.full(39, 2.5e6 / sc.e_unit)
    i_n = ionization_integrals(x, e, sc.niwa())[0]
    i_h = ionization_integrals(x, e, sc.hatakeyama())[0]
    assert i_n > 0 and i_h > 0 and abs(i_n / i_h - 1.0) > 1e-3
