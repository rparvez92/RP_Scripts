# Hand-flux Born SIDIS virtual-photon cross-section convention

Version: `RP_VIRTUAL_PHOTON_XSEC_HAND_V1`

This document defines the convention implemented by
`RP_sidis_reduced_kinematics.py` and `RP_extract_sidis_reduced_xsec.py`. It is
based on the actual `calc_semi_xsec.f` call to `peepiX` in
`simc_gfortran/semi_physics.f`. Neither Fortran source is modified by this
analysis step.

All energies and momenta below are in GeV. `Q2` denotes positive (Q^2).

## 1. Cross section returned by SIMC

`peepiX` first constructs

\[
\frac{d^5\sigma}
 {dE'\,d\Omega_{e'}\,dz\,dP_T^2\,d\phi_{pq}}
=
(\hbar c)^2 10^4
\frac{4\alpha^2E'^2}{Q^4}
\left[
\frac{F_2^A}{\nu}\cos^2\frac{\theta_e}{2}
+\frac{2F_1^A}{M_{\rm struck}}\sin^2\frac{\theta_e}{2}
\right]
\mathcal H_\pi,
\]

where

\[
\mathcal H_\pi=
R_h(A,z)\,
\frac{d\sigma_{\rm frag}}{dz}\,
\frac{b(z)e^{-b(z)P_T^2}}{2\pi}.
\]

It then multiplies by

\[
J_{\rm SIMC}=
\left|\frac{\partial(z,P_T^2)}
{\partial(E_\pi,\cos\theta_{pq})}\right|
=\frac{2P_\pi^2\cos\theta_{pq}}{\nu}
\]

and returns

\[
\frac{d^6\sigma_{\rm lab}}
 {dE'\,d\Omega_{e'}\,dE_\pi\,d\Omega_\pi}
\]

in `microbarn/GeV^2/sr^2` after the output conversion performed by
`calc_semi_xsec.f`. The final SIMC differential is therefore in (E_\pi),
not (P_\pi).

The current SIMC SIDIS model has no explicit \(\cos\phi_{pq}\) or
\(\cos2\phi_{pq}\) terms. Its central model cross section is consequently
constant across the requested phi centers.

## 2. Required pion-variable Jacobian

The desired Hand-flux decomposition uses

\[
\frac{d^6\sigma}
 {dE'\,d\Omega_{e'}\,dP_\pi\,d\Omega_\pi}.
\]

Because

\[
E_\pi=\sqrt{P_\pi^2+m_\pi^2},
\qquad
\frac{dE_\pi}{dP_\pi}=\frac{P_\pi}{E_\pi},
\]

the required transformation is

\[
\frac{d^6\sigma}
 {dE'\,d\Omega_{e'}\,dP_\pi\,d\Omega_\pi}
=
\frac{P_\pi}{E_\pi}
\frac{d^6\sigma}
 {dE'\,d\Omega_{e'}\,dE_\pi\,d\Omega_\pi}.
\]

The dimensionless conversion Jacobian is therefore

\[
J_{E_\pi\rightarrow P_\pi}=P_\pi/E_\pi.
\]

SIMC's internal (J_{\rm SIMC}) must not be inverted for this convention.
Doing so would produce a differential in (z,P_T^2,\phi_{pq}), not the
required differential in (P_\pi,\Omega_\pi).

## 3. Central kinematics and target mass

The central stationary-nucleon quantities are

\[
\nu=\frac{Q^2}{2M_px},\qquad E'=E-\nu,
\]

\[
W^2=M_p^2+2M_p\nu-Q^2,
\]

\[
E_\pi=z\nu,\qquad P_\pi=\sqrt{E_\pi^2-m_\pi^2}.
\]

Both `calc_semi_xsec.f` and full SIMC set `targ%Mtar_struck=Mp` for the
semi-inclusive model. Consequently this central convention uses

\[
M=M_p
\]

for both LH2 and LD2. LD2 enters the model through its proton-plus-neutron
PDF and inclusive-structure-function content and its nuclear factors; the
deuteron mass is not substituted into the central SIDIS electron kinematics.

If LD2 Fermi motion is enabled, full SIMC changes event-level (x) through
(P\cdot q) and applies an additional event-level flux factor. That
correction remains part of the radiated SIMC yield and is not absorbed into
the central Hand flux reported by this table.

## 4. Virtual-photon polarization

\[
\varepsilon=
\left[
1+2\left(1+\frac{\nu^2}{Q^2}\right)
\tan^2\frac{\theta_e}{2}
\right]^{-1}.
\]

The implementation verifies this against

\[
\varepsilon=
\frac{1-y-\gamma^2y^2/4}
 {1-y+y^2/2+\gamma^2y^2/4},
\quad y=\frac{\nu}{E},\quad \gamma=\frac{2M_px}{Q}.
\]

## 5. Hand flux

The adopted flux is exactly

\[
\Gamma_{\rm Hand}=
\frac{\alpha}{2\pi^2}
\frac{E'}{E}
\frac{1}{Q^2}
\frac{1}{1-\varepsilon}
\frac{W^2-M_p^2}{2M_p}.
\]

The Hand equivalent photon energy is

\[
K_{\rm Hand}=\frac{W^2-M_p^2}{2M_p}=\nu(1-x).
\]

The central flux has units `GeV^-1`.

## 6. Extracted Born virtual-photon cross section

The laboratory data cross section is obtained using the existing ratio
method,

\[
\sigma_{\rm lab,Born}^{\rm data}=
\frac{Y_{\rm data}}{Y_{\rm MC}^{\rm rad}}
\sigma_{\rm lab,Born}^{\rm SIMC,central}.
\]

The virtual-photon cross section is

\[
\boxed{
\frac{d^3\sigma_{\gamma^*,\rm Born}}
 {dP_\pi\,d\Omega_\pi}
=
\frac{P_\pi/E_\pi}{\Gamma_{\rm Hand}}
\frac{d^6\sigma_{\rm lab,Born}}
 {dE'\,d\Omega_{e'}\,dE_\pi\,d\Omega_\pi}
}.
\]

Its units are `microbarn/GeV/sr`. With central kinematic uncertainties
neglected, the statistical uncertainty is multiplied by the same factor.

The reconstruction identity is

\[
\sigma_{\rm lab,Born}^{(dE_\pi)}=
\frac{E_\pi}{P_\pi}\Gamma_{\rm Hand}
\sigma_{\gamma^*,\rm Born}.
\]

## 7. Born and radiative treatment

`calc_semi_xsec.f` calls `peepiX` directly and does not invoke SIMC's
radiation generator. Its central value is therefore the Born model cross
section at the supplied vertex kinematics.

When radiation is enabled in full SIMC, `generate_rad` changes event
kinematics and multiplies `main%gen_weight` by the radiative weight. The MC
yield is thus a radiative convolution of the Born `peepiX` model. In the
ratio method above, the radiated MC yield supplies the radiative correction,
while the central multiplier remains Born. No second radiative correction is
applied in this extraction stage.

## 8. Azimuthal decomposition

The extracted observable follows

\[
\begin{aligned}
\frac{d^3\sigma_{\gamma^*}}
 {dP_\pi\,d\Omega_\pi}
={}&\frac{d^3\sigma_T}{dP_\pi\,d\Omega_\pi}
+\varepsilon\frac{d^3\sigma_L}{dP_\pi\,d\Omega_\pi}\\
&+\sqrt{2\varepsilon(1+\varepsilon)}
\frac{d^3\sigma_{LT}}{dP_\pi\,d\Omega_\pi}\cos\phi_{pq}\\
&+\varepsilon
\frac{d^3\sigma_{TT}}{dP_\pi\,d\Omega_\pi}\cos2\phi_{pq}.
\end{aligned}
\]

The flux and Jacobian are evaluated from the central kinematics of each
epsilon-family member and analysis (P_T) bin. No event-averaged flux is
claimed by this convention.
