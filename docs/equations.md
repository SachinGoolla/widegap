# Governing equations, scaling and sign conventions

This document grows stage by stage. Each section names the stage that
introduced it and the code that implements it.

## 1. Electrostatics (stages 1–2)

Poisson's equation for the electrostatic potential $\phi$ (SI units):

$$
-\nabla\cdot(\varepsilon\nabla\phi)=\rho,\qquad \rho=q\,(p-n+N_D^+-N_A^-).
$$

The electric field is $\mathbf E=-\nabla\phi$. Positive $\rho$ curves $\phi$
downward.

### 1.1 Equilibrium carrier statistics (stage 2)

Boltzmann statistics with the quasi-Fermi level as the potential reference
($\phi_F=0$):

$$
n=n_i\,e^{\phi/V_T},\qquad p=n_i\,e^{-\phi/V_T},\qquad V_T=k_BT/q.
$$

Stage 2 assumes complete ionization ($N_D^+=N_D$, $N_A^-=N_A$). This is
**not** adequate for Al acceptors in 4H-SiC near room temperature. That
physics is deferred to stage 7.

### 1.2 Scaling (stage 2, `widegap.fv.equilibrium`)

Pick a reference doping $N_0$. Then:

| Quantity | Unit |
|---|---|
| potential $\psi=\phi/V_T$ | $V_T$ |
| length | $L_D=\sqrt{\varepsilon V_T/(qN_0)}$ (extrinsic Debye length) |
| densities | $N_0$ |
| field | $V_T/L_D$ |

With $\delta=n_i/N_0$ and $C=(N_D-N_A)/N_0$, Poisson's equation becomes

$$
-\psi''=\delta\,(e^{-\psi}-e^{\psi})+C .
$$

The charge-neutral (ohmic contact) potential is
$\psi_{eq}=\operatorname{asinh}(C/(2\delta))$. The built-in potential of an
abrupt junction is $V_{bi}=\psi_{eq}(C_n)-\psi_{eq}(C_p)$.

### 1.3 Depletion approximation for the peak field (stage 2, criterion S5)

For an abrupt junction with $N_A=aN_0$ on the p side and $N_D=bN_0$ on the
n side, the depletion approximation with the $2V_T$ majority-carrier-tail
correction (S. M. Sze, *Physics of Semiconductor Devices*, 2nd ed., Wiley,
1981, ch. 2) gives the total depletion width

$$
W^2=\frac{2\varepsilon}{q}\frac{N_A+N_D}{N_AN_D}\,(V_{bi}-2V_T).
$$

Charge balance gives $x_p=W\,N_D/(N_A+N_D)$. The peak field is
$E_{max}=qN_Ax_p/\varepsilon$. In the scaled units of §1.2 this becomes

$$
E_{max}=\sqrt{\frac{2ab}{a+b}\,(V_{bi}-2)} .
$$

This is an approximation, used here only as a 10% plausibility check (S5).
It is not an exact reference.
