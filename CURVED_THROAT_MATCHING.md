# Curved throat modes and finite-mouth matching

This note supplies a fixed-background mode problem for replacing the constant cylindrical leads. It assumes an existing throat; it neither constructs that throat nor fixes its quantum state. Units are c = hbar = 1, the time dependence is exp(-iEt), and E refers to one common static clock.

## Published starting point

[Maldacena, Milekhin and Popov, Appendix C, (C.10)-(C.18)](https://arxiv.org/pdf/1807.04726) give the spherical metric, vierbein and rescaling Psi = exp(-sigma/2) psi eta / R. Their angular zero modes obey a free two-dimensional equation. A single four-dimensional Weyl field has one radial chirality per angular zero mode, rather than two independent counterpropagating zero channels.

[Freivogel, Fumagalli and Tomasevic, Appendix C, HTML equations (271)-(295)](https://arxiv.org/html/2606.12528v2#A3) retain nonzero angular eigenvalues in the reduced action. The flat kinetic coefficient and the angular term in (293)-(294) imply the optical potential below. Their later Weyl-rescaled gravitational variables must not be substituted for the physical lapse without converting frames. Their interacting collective charged boson is a different object from this single-fermion mode.

For flux magnitude 2 pi q, q > 0 integer, the self-adjoint unit-sphere Dirac spectrum is lambda_0 = 0 with multiplicity q and lambda_n = sqrt(n(n+q)), n >= 1, with multiplicity q+2n for each angular sign. This convention follows [Borokhov, Kapustin and Wu, section 4.1 and appendix](https://arxiv.org/pdf/hep-th/0206054); their section 4 flux integer is our q. The angular differential operator written without -i is anti-Hermitian: its eigenvalue convention must be converted before using real lambda_n.

## Our radial Hamiltonian and current normalization

Take a spherical, static, purely magnetic background, with A_t = A_x = 0 and an x-independent angular bundle:

\[
 ds^2=N(x)^2(-dt^2+dx^2)+r(x)^2d\Omega^2,\qquad
 \Psi={\chi(t,x,\Omega)\over r\sqrt N}.
\]

The rescaling cancels the radial spin-connection term N'/(2N)+r'/r. It also gives the canonical norm integral int dx dOmega chi-dagger chi, since the physical spatial measure is N r^2 dx dOmega. For each nonzero angular pair, constant phase choices put the equation in the form

\[
 i\partial_t f_n=H_n f_n,\quad H_n=-i\sigma_3\partial_x+V_n\sigma_1,
 \qquad V_n(x)={N(x)\lambda_n\over r(x)}.
\]

Consequently the stationary system and an independent squared-operator check are

\[
 f_n'=(iE\sigma_3+V_n\sigma_2)f_n,\qquad
 H_n^2=-\partial_x^2+V_n^2+V_n'\sigma_2.
\]

For real E the transfer matrix satisfies T-dagger sigma_3 T = sigma_3. The conserved radial current per unit common coordinate time is f-dagger sigma_3 f. The proper-area-integrated local current is this quantity divided by N. Thus a radius-normalized physical spinor amplitude acquires a factor sqrt(N_cut) when converted into canonical f; a local flux amplitude must not silently be treated as a common-clock flux amplitude.

For lambda_0 = 0, a right-moving channel satisfies f_+(x_+) = exp(i E L_opt) f_+(x_-), L_opt = x_+-x_-. The opposite chirality travels the reverse way when that field is present. No curvature barrier is generated for this isolated, minimally coupled zero channel. This statement does not include scattering at an unmatched mouth, interaction with a dynamical gauge field, or nonspherical channel mixing.

## One concrete analytic benchmark

Use the leading global AdS2 x S2 profile, with rho = tan(sigma), t = ell tau and x = ell sigma:

\[
 r=r_e,\qquad N={r_e\over\ell\cos(x/\ell)},\qquad
 V_n={\lambda_n\over\ell}\sec(x/\ell),\quad |x|\leq\ell\sigma_c<\pi\ell/2.
\]

Its partner potentials are ell^(-2)[lambda_n^2 sec^2(sigma) +/- lambda_n sec(sigma) tan(sigma)]. For example, q=1 and rho_c=10 specify a finite profile for the first nonzero pair lambda_1=sqrt(2); ell sets the common-clock scale. Integrate its 2 x 2 system, checking real-energy current conservation, refinement and the exact zero-mode phase. A scattering or Riccati formulation avoids overflow in a long evanescent interval. This is a proposed calculation, not a saved numerical result.

At symmetric finite cuts,

\[
 L_{\rm opt}=2\ell\arctan\rho_c,\qquad
 N_c={r_e\over\ell}\sqrt{1+\rho_c^2}.
\]

If the prescribed profile is calibrated to N_c=1 in the external clock, ell>r_e is required and sigma_c=arccos(r_e/ell), so L_opt=2 ell arccos(r_e/ell). This is an algebraic calibration, not a junction solution. In the MMP near-throat approximation the controlled overlap has 1 << rho_c << ell/r_e; N_c=1 instead reaches rho_c approximately ell/r_e, near the opening. Continuing the leading product metric that far does not establish the physical collar. These scales follow [MMP (5.18)-(5.20)](https://arxiv.org/pdf/1807.04726).

Reflecting conditions at the formal AdS endpoints would define a different spectral problem from finite transparent cuts. The limiting length pi ell is an optical light-crossing length, not the throat's proper spatial length.

## Sewing the modes to the exterior

At each cut one must identify the physical radius and induced clock, the pulled-back magnetic connection up to a gauge transformation, and the oriented tangent/normal spin frame. Without a singular interface interaction the physical spinor is related by its spin/gauge transition matrix U. With equal cut radii, canonical amplitudes satisfy chi_throat = sqrt(N_throat/N_exterior) U chi_exterior; a genuine matching in one common time has equal induced lapse. Local energy is E/N. Reversing a port's outward normal also reverses which characteristic is called incoming.

U is obtained from the angular overlap after the frame and bundle identification; it is not generally diagonal. The actual point-pair field at a finite exterior cut differs from an isolated monopole. Therefore the present cylinder basis cannot be assumed to give the exact curved junction even when its flux norm is correct. A nonspherical collar would also produce radial derivative couplings between angular modes.

For a single transparently identified zero channel the conditional propagation factor is

\[
 z(E)=e^{i\varphi_{\rm end}}e^{iE L_{\rm opt}},\qquad
 z(i\kappa)=e^{i\varphi_{\rm end}}e^{-\kappa L_{\rm opt}},\quad\kappa>0.
\]

The endpoint phase is fixed only after specifying spin/gauge-frame matching, any longitudinal Wilson line and the global spin structure. An endpoint rephasing changes the matrix representation, while the complete sewn-loop holonomy is invariant. Testing lengths 1 and 4 and phases 0 and pi is consequently a useful conditional sewing check, not a determination of the physical phase or curved mouth geometry.

An imaginary-axis exterior map and this propagation factor supply ingredients of a closed resolvent. They do not determine the renormalized stress: that requires the full angular/frequency content, a specified global quantum state and spin structure, controlled ultraviolet subtraction, and consistent variation with respect to geometry and gauge field. Selecting N and r does not solve those field equations. No stable support medium, formation process, or backward-time path follows from this benchmark.
