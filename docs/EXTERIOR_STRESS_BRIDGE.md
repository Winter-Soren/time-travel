# From the certified exterior probe to a stress calculation

Increment120. The [119 certificate](../results/dihole_probe_certificate_results.json)
bounds one compact resolvent matrix element in the exact classical balanced
dihole, for $q=512000$, $\mu=398908.5$ and scaled $\kappa=1$.
It supplies neither a stress tensor nor negative energy. This note specifies a
quantum state, the next metric response that can be bounded, and the remaining
steps. The background is Emparan's uncut black-hole exterior, with Melvin
asymptotics ([equations8--16](https://arxiv.org/pdf/hep-th/9906160)); no horizon
has been replaced by a throat.

## A definite exterior state

Let $H_D$ be the natural self-adjoint, full massless charged Dirac Hamiltonian,
including both chiralities and every azimuthal sector. Use the exact balanced
parameters, integer-$q$ magnetic bundle and time coordinate of the
[domain proof](DIHOLE_DIRAC_DOMAIN.md). The charge-included connection is
$A=q\,a(\rho,z)d\varphi$ in the outer patch. Fix the gauge-invariant quasifree
ground state by its occupation covariance

$$
C_-=\mathbf1_{(-\infty,0)}(H_D),\qquad C_+=I-C_-.
$$

**Zero-energy modes are empty.** This is a convention, not a claim that the
kernel vanishes. The canonical two-point kernels are the distribution kernels of

$$
K_\pm(t,t')=e^{-i(t-t')H_D}C_\pm.
$$

Here $K_-$ represents $\langle\Psi^\dagger(t',y)\Psi(t,x)\rangle$ with its
indices ordered as an operator kernel; $K_+$ is the opposite CAR ordering.
Recover physical spinor/cospinor kernels using the inverse canonical Hilbert map
and Dirac adjoint. A single fixed-$\mu$ kernel is not the four-dimensional state.

The Hadamard hypotheses can be checked on this **open exterior**. Its complete
optical spatial metric makes the ultrastatic conformal spacetime globally
hyperbolic: bounded time intervals bound optical travel distance, and complete
Riemannian closed balls are compact. The physical conformal metric has the same
causal structure. Time translations are complete, everywhere timelike, and
preserve the smooth magnetic connection. The charged Dirac equation, doubled
with its conjugate bundle, has a normally hyperbolic square; its two-point
functions are bisolutions and satisfy CAR locality. The spectral state is a
ground state. Thus Sahlmann--Verch's
[passivity theorem5.1](https://arxiv.org/pdf/math-ph/0002021) gives the microlocal
spectrum condition, and their [Dirac theorem5.8 and remark(ii)](https://arxiv.org/pdf/math-ph/0008029)
give Hadamard form. Conjugation exchanges the magnetic bundle with its dual.
The charged ultrastatic construction is explicit in
[Zahn, section3, pp13--14 and footnote12](https://arxiv.org/pdf/1210.4031): a
time-independent spatial gauge connection and spectral projections produce the
Hadamard covariances. Our massless conformal map transfers this construction
to the physical static metric without changing frequency signs.
No spectral gap or global positive lower bound on the lapse is assumed.
This establishes Hadamard behavior at exterior points, **not regular extension
across either horizon**. A bounded occupation change supported in $\ker H_D$
has a kernel annihilated by $H_D$ in both variables. The sum of their squares
is elliptic on the spatial product, so this kernel is locally smooth.

There is an exact resolvent bridge, with explicit infrared and ultraviolet
cutoffs. If $P_0=\mathbf1_{\{0\}}(H_D)$, spectral calculus gives

$$
C_-=\mathop{\mathrm{s\!-!lim}}_{\alpha\downarrow0,\,\Lambda\uparrow\infty}
\left\{\frac{I-P_0}{2}-\frac1{2\pi}\int_\alpha^\Lambda
\left[(H_D-i\omega)^{-1}+(H_D+i\omega)^{-1}\right]d\omega\right\}.
$$

For a nonzero spectral value, the integral reduces to an arctangent and tends
to half its sign. This is a strong operator limit, not permission to interchange
the limit with a trace, coincidence limit or derivative.

## The next bounded metric response

Write $h_{\rm sp}$ for the spatial metric and $h$ for the fixed compact source.
Hold $h_{\rm sp}$, the magnetic connection, canonical Hilbert map and baseline
$E_q=v_*/L$ fixed. Vary only

$$
g_\epsilon=-N^2(1+\epsilon f)^2dt^2+h_{\rm sp},\qquad
f\in C_c^\infty(X,\mathbb R),\quad \|f\|_\infty\le1,\quad |\epsilon|<1.
$$

This is an off-shell metric variation. Substituting a new lapse into the old
Weyl-coordinate metric formulas would also change the spatial metric and would
be a different variation.

Use the baseline affine coordinates $(u,v)$ and the three bumps

$$
f_c=\zeta(u)\zeta(v)\exp\!\left[-\frac{(u-c)^2+(v-c)^2}{8}\right],
\qquad c=-2,0,2.
$$

Define $s(t)=e^{-1/t}$ for $t>0$ and $s(t)=0$ otherwise,
$S(t)=s(t)/(s(t)+s(1-t))$. Set $\zeta(x)=1$ for $|x|\le9$,
$\zeta(x)=1-S((|x|-9)/2)$ for $9<|x|<11$, and zero beyond.
This is smooth, including at the joins. The support has $\rho>.484$, so extension
by zero defines a smooth axisymmetric function on the physical exterior.
On the entire trial square $[-9,9]^2$, the bump is exactly the analytic Gaussian;
its one-complex-coordinate ellipse modulus is below $e^{.75^2/8}<1.073$.
At $\kappa=1$, the overlap integrands remain within the 119 common
$3\times10^{25}$ analytic modulus budget.
The separate choice $f=1$ is only a uniform clock-rescaling control.

For the massless Hamiltonian, the first-order Leibniz rule gives **exactly**

$$
H_\epsilon=(1+\epsilon f)^{1/2}H(1+\epsilon f)^{1/2}
=H+\epsilon V,\qquad V=\tfrac12(fH+Hf).
$$

Here $H$ may be the scaled probe restriction. The compact commutator $[H,f]$ is
bounded; multiplication by $(1+\epsilon f)^{\pm1/2}$ preserves $D(H)$.
Consequently this is a common-domain self-adjoint analytic family. The same
fixed spatial Hilbert norm is used throughout.

Let $R=(H-i\kappa)^{-1}$, $\psi=Rh$, $\chi=R^*h$, and
$G=\langle h,Rh\rangle$, with inner products antilinear in the first slot.
Resolvent differentiation gives $\delta G=-\langle\chi,V\psi\rangle$.
Using $H\psi=h+i\kappa\psi$ and $H\chi=h-i\kappa\chi$ yields

$$
\boxed{\delta G=-\tfrac12\bigl(\langle\chi,fh\rangle+
\langle h,f\psi\rangle\bigr)-i\kappa\langle\chi,f\psi\rangle.}
$$

Only bounded multiplication overlaps remain. For graph-domain approximants
$\widetilde\psi,\widetilde\chi$ with certified errors $e_\psi,e_\chi$, the overlap
formula's error is at most

$$
\frac{\|fh\|}{2}(e_\psi+e_\chi)
+\kappa\left(e_\chi\|f\widetilde\psi\|
+e_\psi\|f\widetilde\chi\|+e_\psi e_\chi\right).
$$

Residual norms divided by $\kappa$ supply these solution-error bounds.
Quadrature and arithmetic errors in the overlaps must also be enclosed.
The global error bounds cover the exact solutions outside the compact trial
support; those contributions are not set to zero.
This observable is a **packet lapse susceptibility**, independent of choosing
the ground state. Its sign is not an energy-density or null-energy test.

## What a stress result still requires

A renormalized stress calculation must reconstruct the full matrix two-point
function, including both chiralities, all azimuthal modes and the complete
frequency range. In point-split notation its structure is

$$
\langle T_{ab}\rangle_{\rm ren}
=\lim_{y\to x}\mathcal D_{ab}\bigl(W_\omega-H_{\rm Had}\bigr)(x,y)
+C_{ab}(x).
$$

Here $\mathcal D_{ab}$ includes the Dirac stress derivatives and spin/gauge
parallel transport; $H_{\rm Had}$ is the local charged-field parametrix.
The renormalization length, local gravitational/gauge counterterms and Ward
identities must be fixed consistently. Matter stress need not be separately
conserved in a nonzero external electromagnetic field. See
[Dappiaggi--Hack--Pinamonti](https://arxiv.org/pdf/0904.0612) for Dirac
point-splitting and [Zahn, sections3--4](https://arxiv.org/pdf/1210.4031) for
background gauge fields and renormalization freedom.

The [calculation](../scripts/dihole_lapse_susceptibility.py) now encloses the three
[specified susceptibilities](../results/dihole_lapse_susceptibility_results.json).
Stress then requires controlled differentiated kernels and subtracted mode
sums; backreaction additionally requires Einstein--Maxwell consistency.
Neither the 119 scalar disk nor this metric derivative proves a supported
throat, stability or travel to yesterday.
