"""Stage 18: accumulation-surface mobility option of the channel map (edges, default path, registry)."""

from __future__ import annotations

import numpy as np
import pytest

from widegap.channel_mobility import ACC_NAME, ChannelMap, acc_mobility_from_registry
from widegap.materials.registry import Registry
from widegap.sic import DMOSDesign, SiC4H, SiCDMOS, TracedRegistry

DESIGN = DMOSDesign(
    x_jfet=0.8, ch_len=0.8, src_w=1.0, bc_w=0.8, pad_w=0.4, d_src=0.3, d_body=0.8, w_drift=3.0,
    t_drain=0.3, t_ox=0.05, n_drift=5e16, n_jfet=1e17, n_body=3e17, n_src=1e19, n_drain=1e19,
    mesh_kw=dict(ny_surf=10, q_surf=0.6),
)  # fmt: skip


@pytest.fixture(scope="module")
def stack():
    reg = TracedRegistry.from_json()
    sc = SiC4H(reg)
    dev = SiCDMOS(sc, DESIGN, 0)
    mu_acc = acc_mobility_from_registry(reg)
    return reg, sc, dev, mu_acc


def _maps(sc, dev, mu_acc):
    kw = dict(d_ch_um=20.0e-3, gate_x_max=dev.prm.gate_x1)
    base = ChannelMap(sc, dev.model, dev.nd, dev.na, **kw)
    acc = ChannelMap(sc, dev.model, dev.nd, dev.na, mu_acc_cm2=mu_acc, **kw)
    return base, acc


def test_registry_read_not_hardcoded(stack):
    reg, _, _, mu_acc = stack
    assert mu_acc == 17.5 == reg.reg.get(ACC_NAME).value
    assert ACC_NAME in reg.used  # read through the (traced) registry
    # an unverified copy of the entry is refused: no fallback value
    ent = dict(next(e for e in reg.reg._raw.values() if e["name"] == ACC_NAME))
    ent["verified_by"] = None
    with pytest.raises(Exception, match="verif|Provenance|verified"):
        acc_mobility_from_registry(TracedRegistry(Registry([ent])))


def test_default_is_bit_identical(stack):
    _, sc, dev, mu_acc = stack
    base, acc = _maps(sc, dev, mu_acc)
    assert base.mu_acc_cm2 is None and len(base.acc_edges) == 0 and not base.col_acc.any()
    # stage-16 fields untouched by the option
    assert np.array_equal(base.edges, acc.edges) and np.array_equal(base.cols, acc.cols)
    assert np.array_equal(base.mu_bulk, acc.mu_bulk)
    mu_col = np.full(len(base.cols), 80.0)
    m0, m1 = base.edge_mobility(mu_col), acc.edge_mobility(mu_col)
    sel = np.ones(len(m0), bool)
    sel[acc.acc_edges] = False
    assert np.array_equal(m0[sel], m1[sel])  # everything outside the accumulation edges
    # default: bulk everywhere except the stage-16 channel edges (input untouched)
    keep = np.setdiff1d(np.arange(len(m0)), base.edges)
    assert np.array_equal(m0[keep], base.mu_bulk[keep])
    with pytest.raises(ValueError):
        ChannelMap(sc, dev.model, dev.nd, dev.na, 20.0e-3, mu_acc_cm2=mu_acc)  # no gate limit


def test_edge_selection_jfet_surface_only(stack):
    _, sc, dev, mu_acc = stack
    base, acc = _maps(sc, dev, mu_acc)
    mesh = dev.dev.mesh
    p = dev.prm
    ix = np.arange(mesh.n_nodes) // (mesh.ny + 1)
    cols = np.flatnonzero(acc.col_acc)
    assert len(cols) > 2 and len(acc.acc_edges) > 0
    surf = acc.surf_node[cols]
    assert np.all(dev.nd[surf] > dev.na[surf])  # n-type surface
    assert np.all(mesh.x[cols] <= p.gate_x1 + 1e-9)  # under the gate
    assert np.all(mesh.x[cols] <= p.x_jfet + 1e-9)  # JFET surface (before the p-body channel)
    assert cols[0] == 0 and np.array_equal(cols, np.arange(len(cols)))  # contiguous from x = 0
    assert not (acc.col_acc & base.col_p).any()  # never a p-type (inversion) column
    assert not acc.col_acc[mesh.x >= p.x_src0 - 1e-9].any()  # n+ source side excluded
    e = acc.acc_edges
    a, b = dev.model.a[e], dev.model.b[e]
    assert np.all(mesh.edge_dir[e] == 0)  # lateral
    assert np.all(acc.in_ch[a] & acc.in_ch[b])  # both nodes within d_ch
    assert np.all(acc.col_acc[ix[a]] & acc.col_acc[ix[b]])
    assert not np.intersect1d(e, acc.edges).size  # disjoint from the stage-16 channel edges
    # applying the option sets the constant (scaled) mobility on exactly those edges
    mu_col = np.full(len(acc.cols), 80.0)
    mu, m0 = acc.edge_mobility(mu_col), base.edge_mobility(mu_col)
    assert np.all(mu[e] == sc.mu_s(17.5))
    other = np.setdiff1d(np.arange(len(mu)), e)
    assert np.array_equal(mu[other], m0[other])  # nothing else changes
