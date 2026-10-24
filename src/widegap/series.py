r"""Stage 26: lumped series terms added to R_on,sp AFTER the solve (never inside the PDE).

The half-cell R_on,sp is ``R = V_DS x_max / I_D`` with ``I_D`` the drain current per cm of device
depth and ``x_max`` the half-pitch (``adjoint_sic.SiCAdjoint.ron``).  A resistance per unit depth
(ohm cm) times ``x_max`` is therefore a specific resistance (ohm cm^2) of the half cell:

* source contact (specific contact resistivity ``rho_c`` [ohm cm^2] over the modelled source
  contact length ``L_c`` [cm], per cm depth):  ``Delta R_c = (rho_c / L_c) x_max``  -- the
  half-pitch over contact-length rule.  No current crowding / transfer-length correction (the
  n+ sheet resistance is already in the mesh; ``L_c`` is geometric).
* substrate:  ``R_sub = rho_sub t_sub`` [ohm cm^2] (a vertical 1D resistor over the whole cell
  area, so no geometric factor).

``source_contact_length_cm`` is the length of the top-node set the model grounds as source
(``DMOSFET.nodes["source"]``: from the first node at ``x >= x_src0 + src_contact_gap`` to
``x_src1``).
"""

from __future__ import annotations


def contact_resistance_sp(rho_c_ohm_cm2: float, x_max_cm: float, l_c_cm: float) -> float:
    """Delta R_c = rho_c x_max / L_c  [ohm cm^2]."""
    if l_c_cm <= 0.0:
        raise ValueError("source contact length must be positive")
    return float(rho_c_ohm_cm2) * float(x_max_cm) / float(l_c_cm)


def substrate_resistance_sp(rho_sub_ohm_cm: float, t_sub_cm: float) -> float:
    """R_sub = rho_sub t_sub  [ohm cm^2]."""
    return float(rho_sub_ohm_cm) * float(t_sub_cm)


def source_contact_length_cm(dmos) -> float:
    """Modelled source-contact length of a :class:`widegap.sic.SiCDMOS` (cm)."""
    xy = dmos.dev.mesh.node_xy
    nodes = dmos.dev.nodes["source"]
    x0 = float(xy[nodes, 0].min())
    return float(dmos.sc.length_cm(dmos.prm.x_src1 - x0))
