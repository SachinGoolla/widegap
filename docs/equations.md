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

## 2. Reverse-bias electrostatics and the breakdown index (stage 4)

### 2.1 Zero-current approximation

At reverse bias $V_R$ with negligible current, the quasi-Fermi potentials
are flat: $\phi_n=0$ (cathode) and $\phi_p=-V_R$ (anode). Then
$n=\delta e^{\psi-\phi_n}$, $p=\delta e^{\phi_p-\psi}$, and
$-\psi''=\delta(e^{\phi_p-\psi}-e^{\psi-\phi_n})+C$ (`widegap.fv.reverse_bias`).

### 2.2 Exact first integral

Multiply by $\psi'$ and integrate. Since $E=-\psi'$ vanishes in a neutral
bulk,
$$\tfrac12E(x_0)^2=\int_{\psi(x_0)}^{\psi_{bulk}}\rho(\psi)\,d\psi .$$
For a one-sided p⁺n junction ($N_A=aN_0$, $N_D=bN_0$) this gives, on the two
sides,
$\tfrac12E(0)^2=b\Delta\psi_n-b+p(0)=a\Delta\psi_p-a+p(0)$ (minority
terms dropped), with $p(0)=a\,e^{-\Delta\psi_p}$. Hence
$$\frac{d(E^2)}{dV}=2\,\frac{ab-b\,p(0)}{a+b}\quad(\text{for }b=1),$$
which reduces to the depletion-approximation value $2ab/(a+b)$ only when
$p(0)\ll a$. In stage 4, $p(0)\approx19$ for $a=100$, so this
correction is about 20%.

### 2.3 Breakdown index (Zhu & Andrei 2019, eq. 1)

$B=\max_e|E_e|/E_{crit}$ over element fields $E_e=-(\psi_{e+1}-\psi_e)/h_e$.
The field-based BV solves $B(V_R)=1$. Because the maximum sits on the
element at the junction, where $|E'|=|\rho(0)|$, the discrete maximum
converges at $O(h)$.

## 3. Flat band of an n⁺ poly-Si gate on 4H-SiC (stage 12, `widegap.sic`)

**Convention.** The stage-8.1 MOS solver imposes $\psi_{gate}=V_G-V_{FB}^{code}$ on
the gate edge of the oxide, with $\psi$ the electrostatic potential in $V_T$
normalised so that $n=n_i e^{\psi-\phi_n}$, i.e. the intrinsic level is
$E_i(\mathbf x)=E_{i0}-q\,\psi(\mathbf x)\,V_T$ (energies in eV, potentials in V
below) and the source/body Fermi level is $E_F=0$.

**Derivation.**

1. Vacuum-level continuity (no interface dipole): in the semiconductor and
   throughout the oxide
   $E_{vac}(\mathbf x)=E_{c0}+\chi_{SiC}-q\psi(\mathbf x)$, where
   $E_{c0}=E_{i0}+(E_c-E_i)_0$ is the conduction-band edge at $\psi=0$ and
   $\chi_{SiC}$ the SiC electron affinity.
2. Gate: $E_{F,gate}=E_F-qV_G=-qV_G$ and $E_{vac,gate}=E_{F,gate}+\Phi_{gate}$, hence
   $$\psi_{gate}=V_G+\frac{E_{c0}+\chi_{SiC}-\Phi_{gate}}{q}.$$
3. Neutral p-body: $E_F=0$ at $\psi=\psi_b$ gives
   $E_{c0}=q\psi_b+(E_c-E_F)_{bulk}$.
4. n⁺ poly-Si: $E_{F,gate}=E_c(\mathrm{Si})$, so
   $\Phi_{gate}-\chi_{SiC}=(\chi_{Si}-\chi_{ox})-(\chi_{SiC}-\chi_{ox})
   =\Phi_{C,Si}-\Phi_{C,SiC}$, the conduction-band barriers of Si and SiC to SiO₂
   (Afanas'ev et al., J. Appl. Phys. 79, 3108 (1996): 3.15 eV and 2.7 eV).
5. Therefore
   $$\psi_{gate}=V_G+\psi_b-V_{FB},\qquad
   V_{FB}=\frac{\Phi_{C,Si}-\Phi_{C,SiC}}{q}-\frac{(E_c-E_F)_{bulk}}{q},$$
   the usual flat-band voltage, with $(E_c-E_F)_{bulk}=k_BT\ln(N_c/n_b)$,
   $n_b=n_i e^{\psi_b}$ (stage-7 statistics, incomplete-ionization neutral root $\psi_b$).
6. In the solver convention $V_{FB}^{code}=V_{FB}-\psi_bV_T$ (in $V_T$ units
   $V_{FB}/V_T-\psi_b$). Since $(E_c-E_F)_{bulk}=k_BT\,(\ln(N_c/n_i)-\psi_b)$,
   $$V_{FB}^{code}=\frac{(\Phi_{C,Si}-\Phi_{C,SiC})-k_BT\ln(N_c/n_i)}{k_BT},$$
   independent of the body doping: $\psi=0$ already fixes $(E_c-E_i)_0$.

For the stage-12 body ($N_A=3\cdot10^{17}\,\mathrm{cm^{-3}}$ Al, $T=300$ K) this gives
$V_{FB}\approx-2.6$ V ($V_{FB}^{code}\approx-45.3\,V_T$); the numbers are recorded in
`results/stage12/results.json` (`12b.flat_band`). Assumptions flagged in the stage-12
results: $E_{F,gate}=E_c(\mathrm{Si})$ (degenerate n⁺ poly-Si), Si-face/dry-oxide offsets,
$\Phi_{C,SiC}=\Phi_V-E_g$ derived by the offset source, no interface dipole.
