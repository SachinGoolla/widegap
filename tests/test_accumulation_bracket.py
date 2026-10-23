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


# ------------------------------------------------------------------ stage 18b (excess mode)
def _excess(sc, dev, mu_acc):
    return ChannelMap(sc, dev.model, dev.nd, dev.na, 20.0e-3, gate_x_max=dev.prm.gate_x1,
                      mu_acc_cm2=mu_acc, acc_mode="excess")  # fmt: skip


def test_18b_defaults_and_all_mode_bit_identical(stack):
    _, sc, dev, mu_acc = stack
    base, acc = _maps(sc, dev, mu_acc)
    assert base.acc_mode == acc.acc_mode == "all" and base.nfb_node is None
    mu_col = np.full(len(base.cols), 80.0)
    # the new optional n argument is ignored outside excess mode
    n = np.full(dev.model.mesh.n_nodes, 1.0)
    assert np.array_equal(acc.edge_mobility(mu_col), acc.edge_mobility(mu_col, n))
    assert np.array_equal(base.edge_mobility(mu_col), base.edge_mobility(mu_col, n))
    with pytest.raises(ValueError):
        ChannelMap(sc, dev.model, dev.nd, dev.na, 20.0e-3, gate_x_max=dev.prm.gate_x1,
                   mu_acc_cm2=17.5, acc_mode="bogus")  # fmt: skip
    with pytest.raises(ValueError):
        _excess(sc, dev, mu_acc).edge_mobility(mu_col)  # needs n


def test_18b_excess_limits(stack):
    _, sc, dev, mu_acc = stack
    ex = _excess(sc, dev, mu_acc)
    e = ex.acc_edges
    assert len(e) > 0 and np.all(ex.nfb_edge > 0)
    mu_col = np.full(len(ex.cols), 80.0)
    nodes = np.zeros(dev.model.mesh.n_nodes)
    a, b = dev.model.a[e], dev.model.b[e]
    # depletion / flat band (n <= n_fb): bulk mobility
    nodes[:] = 0.5 * np.median(ex.nfb_node[np.isfinite(ex.nfb_node)])
    nodes[a], nodes[b] = ex.nfb_edge, ex.nfb_edge
    mu = ex.edge_mobility(mu_col, nodes)
    assert np.array_equal(mu[e], ex.mu_bulk[e])
    nodes[a], nodes[b] = 0.3 * ex.nfb_edge, 0.3 * ex.nfb_edge
    assert np.array_equal(ex.edge_mobility(mu_col, nodes)[e], ex.mu_bulk[e])
    # strong accumulation: tends to mu_acc; always between mu_acc and mu_bulk
    nodes[a], nodes[b] = 1e6 * ex.nfb_edge, 1e6 * ex.nfb_edge
    mu = ex.edge_mobility(mu_col, nodes)[e]
    assert np.allclose(mu, sc.mu_s(17.5), rtol=1e-4)
    nodes[a], nodes[b] = 3.0 * ex.nfb_edge, 3.0 * ex.nfb_edge
    mu = ex.edge_mobility(mu_col, nodes)[e]
    assert np.allclose(mu, ex.mu_bulk[e] * (1 / 3) + sc.mu_s(17.5) * (2 / 3), rtol=1e-12)
    # mu_acc := mu_bulk reduces to bulk for any density
    bk = ChannelMap(sc, dev.model, dev.nd, dev.na, 20.0e-3, gate_x_max=dev.prm.gate_x1,
                    mu_acc_cm2="bulk", acc_mode="excess")  # fmt: skip
    nodes[bk.model.a[e]] = np.linspace(0.5, 50, len(e)) * bk.nfb_edge
    nodes[bk.model.b[e]] = nodes[bk.model.a[e]]
    assert np.allclose(bk.edge_mobility(mu_col, nodes)[e], bk.mu_bulk[e], rtol=1e-14)
    # other edges untouched
    other = np.setdiff1d(np.arange(len(bk.mu_bulk)), e)
    assert np.array_equal(bk.edge_mobility(mu_col, nodes)[other],
                          _maps(sc, dev, mu_acc)[0].edge_mobility(mu_col)[other])  # fmt: skip


def test_18b_flat_band_density_is_ionized_donors(stack):
    _, sc, dev, mu_acc = stack
    ex = _excess(sc, dev, mu_acc)
    nfb = ex.nfb_node[np.isfinite(ex.nfb_node)] * sc.n0
    # JFET N_D = 1e17: incomplete ionization -> n_fb a bit below N_D, same order
    assert np.all(nfb < 1.0e17 * 1.001) and np.all(nfb > 0.5e17)
