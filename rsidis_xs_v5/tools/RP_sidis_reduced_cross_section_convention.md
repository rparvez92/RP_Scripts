# SIMC-native reduced SIDIS cross-section convention

Version: `RP_REDUCED_XSEC_SIMC_LAB_V1`

This note records the convention implemented by
`RP_sidis_reduced_kinematics.py`. It is derived from the central, non-Fermi
calculation in `simc_gfortran/semi_physics.f`; it does not modify SIMC or
`calc_semi_xsec`.

All energies and momenta below are expressed in GeV. The symbol `Q2` in the
Python implementation denotes the positive quantity \(Q^2\) in GeV\(^2\).

## Central SIMC laboratory cross section

SIMC calculates

\[
\frac{d^6\sigma_{\rm lab}}
 {dE'\,d\Omega_{e'}\,dE_h\,d\Omega_h}
=
(\hbar c)^2 10^4
\frac{4\alpha^2E'^2}{Q^4}
\left[
W_2\cos^2\frac{\theta_e}{2}
+2W_1\sin^2\frac{\theta_e}{2}
\right]
\mathcal M_h
J_h,
\]

with

\[
J_h=
\left|\frac{\partial(z,P_T^2)}
{\partial(E_h,\cos\theta_{pq})}\right|
=\frac{2p_h^2\cos\theta_{pq}}{\nu}
\]

for the forward central kinematics used by the analysis. The factor \(10^4\)
converts fm\(^2\) to microbarns after \(\hbar c\) is expressed in GeV fm.

The central kinematics are

\[
\nu=\frac{Q^2}{2Mx},\quad E'=E-\nu,\quad
Q^2=4EE'\sin^2\frac{\theta_e}{2},
\]

\[
E_h=z\nu,\quad p_h=\sqrt{E_h^2-m_h^2},\quad
P_T=p_h\sin\theta_{pq}.
\]

The constants in the Python module reproduce `simc_gfortran/constants.inc`:

- \(\alpha=1/137.0359895\),
- \(\hbar c=0.197327053\) GeV fm,
- \(M=0.93827231\) GeV,
- \(m_{\pi^\pm}=0.13957018\) GeV.

## Transformation to the transverse/longitudinal basis

SIMC uses

\[
W_1=\frac{F_1}{M},\qquad W_2=\frac{F_2}{\nu}.
\]

Define

\[
\gamma=\frac{2Mx}{Q},\qquad
F_T=2xF_1,\qquad
F_L=(1+\gamma^2)F_2-2xF_1,
\]

and

\[
\varepsilon=
\left[1+2\left(1+\frac{\nu^2}{Q^2}\right)
\tan^2\frac{\theta_e}{2}\right]^{-1}.
\]

Direct algebra gives

\[
W_2\cos^2\frac{\theta_e}{2}
+2W_1\sin^2\frac{\theta_e}{2}
=
\frac{\cos^2(\theta_e/2)}
{\nu\varepsilon(1+\gamma^2)}
\left(F_T+\varepsilon F_L\right).
\]

The same \(\varepsilon\) is obtained from

\[
y=\frac{\nu}{E},\qquad
\varepsilon=
\frac{1-y-\gamma^2y^2/4}
{1-y+y^2/2+\gamma^2y^2/4}.
\]

## Reduced normalization

The convention therefore defines

\[
K_{\rm base}=
(\hbar c)^2 10^4
\frac{4\alpha^2E'^2}{Q^4}J_h,
\]

\[
A_{\rm lab}=
\frac{\cos^2(\theta_e/2)}
{\nu\varepsilon(1+\gamma^2)},
\]

and

\[
K_{\rm red}=K_{\rm base}A_{\rm lab}.
\]

The reduced cross section is

\[
\sigma_{\rm red}
=\frac{\sigma_{\rm lab}}{K_{\rm red}}
=\left(F_T+\varepsilon F_L\right)\mathcal M_h.
\]

When uncertainty in the kinematic normalization is neglected, its statistical
uncertainty is

\[
\delta\sigma_{\rm red}
=\frac{\delta\sigma_{\rm lab}}{K_{\rm red}}.
\]

This version applies only to the central stationary-nucleon convention. The
event-dependent Fermi-motion flux correction used by some deuterium SIMC runs
is intentionally outside its scope.
