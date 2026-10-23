"""Stage 21: the 13d machinery on the deeper (d_body = 1.4 um) body_design geometry (small mesh)."""

import numpy as np
import pytest

from widegap.adjoint_sic_oxide import WINDOWS_UM, OxideAdjoint
from widegap.bodyredesign import body_design, column_index_matrix, source_columns
from widegap.sic import DMOSDesign, SiC4H, TracedRegistry

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

BASE = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4,
    d_src=0.3, d_body=0.8, w_drift=3.0, t_drain=0.3, t_ox=0.05,
    n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(nx=(2, 1, 2, 2, 2, 1), ny_drain=1, ny_drift=5, ny_body=4, ny_surf=5, ny_ox=1,
                 q_surf=0.6),
)  # fmt: skip
SHIFT = 0.6


def test_deeper_device_keeps_x_and_shifts_y_and_windows_map():
    sc = SiC4H(TracedRegistry.from_json())
    deep = body_design(BASE, 1.4)
    assert deep.d_body == 1.4 and deep.mesh_kw["ny_body"] > BASE.mesh_kw["ny_body"] - 1
    shifted = {
        "W1": ((0.5, 1.1), (-1.1 - SHIFT, -0.5 - SHIFT)),
        "W2": ((1.3, 2.0), (-1.2 - SHIFT, -0.6 - SHIFT)),
        "W3": WINDOWS_UM["W3"],
    }
    e0 = OxideAdjoint(sc, BASE, 0)
    e1 = OxideAdjoint(sc, deep, 0, windows_um=shifted)
    assert np.array_equal(e0.mesh.x, e1.mesh.x)  # x unchanged
    dy = float(sc.length_s(SHIFT * 1e-4))
    assert e1.mesh.y.min() == pytest.approx(e0.mesh.y.min() - dy, rel=1e-9)
    info = e1.window_info()
    assert all(info[k]["n_edges"] > 0 for k in ("W1", "W2", "W3", "E_ox"))
    assert info["E_ox"]["n_edges"] == e0.window_info()["E_ox"]["n_edges"]
    # the stage-20 barrier columns exist on the deep device
    ev = OxideAdjoint(sc, deep, 0)
    ev.set_theta(np.zeros(5))
    cols = source_columns(ev)
    idx = column_index_matrix(ev, cols)
    assert len(cols) >= 1 and idx.shape[0] == len(cols)
