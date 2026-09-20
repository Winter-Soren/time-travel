# Dirac subtraction-scale matching and the remaining contact definition

Increment136. The logarithmic conversion between the convention in
[133](SEWN_SOURCE_STATE.md) and the heat/zeta convention in
[135](PHYSICAL_MASS_REGULATOR.md) can be fixed. The complete charged
stress conversion cannot yet be assigned numerical finite constants:
133 has not specified an explicit corrected point-split stress/current
pair. This is a missing local definition, not an exterior PDE problem.

## The logarithmic scale is determined

Keep133's convention

$$
G_{H,\ell}(x,x')=\frac1{8\pi^2}
\left[\frac{U(x,x')}{\sigma}+V(x,x')\log\frac{\sigma}{\ell^2}\right],
\qquad \sigma=\tfrac12(\hbox{geodesic distance})^2,
$$

understood here on the Euclidean section for the auxiliary Laplace-type
spinor Green kernel. Applying the Dirac operator gives the first-order
parametrix. This step alone does not choose its stress contact correction.

A flat massive kernel fixes the otherwise easily missed constants. For
$r^2=2\sigma$ and positive test mass $m$,

$$
G_m(r)=\frac{mK_1(mr)}{4\pi^2r}
=\frac1{8\pi^2\sigma}
+\frac{m^2}{16\pi^2}
\left[\log\frac{m^2\sigma}{2}+2\gamma_E-1\right]+o(1).
$$

Subtracting $G_{H,\ell}$ leaves
$m^2[\log(m^2\ell^2/2)+2\gamma_E-1]/(16\pi^2)$.
Independently, the diagonal zeta kernel is

$$
\operatorname{FP}_{z=0}\left[
\mu_R^{2z}\frac{(m^2)^{1-z}}{16\pi^2z(z-1)}\right]
=\frac{m^2}{16\pi^2}
\left[\log\frac{m^2}{\mu_R^2}-1\right].
$$

Their equality requires

$$
\boxed{\mu_H=\frac{\sqrt2e^{-\gamma_E}}{\ell}}
\simeq\frac{0.794023616383283}{\ell}.
$$

This fixes the universal logarithmic normalization, independently of the
chosen test mass. The same relation is found for neutral Dirac stress on
AdS in [Ambrus--Winstanley, equations41 and52--56](https://arxiv.org/pdf/1505.04962),
after their specified stress correction. That example is not a theorem
equating arbitrary charged stress prescriptions.

Changing the logarithm to $\log(2\sigma/\ell^2)$ would instead give
$\mu_H=2e^{-\gamma_E}/\ell$. Zahn uses $\Gamma=-2\sigma$ and
$\log(-\Gamma/\Lambda^2)$, so matching133 requires
$\Lambda=\sqrt2\ell$; see
[Zahn, section3.1 and equation41](https://arxiv.org/pdf/1210.4031).
Setting either scale equal to $1/\ell$ without conversion changes a finite
local term.

## A computable action change, not the completed matching

In135's convention the parity-even Euclidean fermion action is
$\Gamma_\zeta(\mu_R)=\tfrac12\zeta'_{\mathscr P/\mu_R^2}(0)$.
For a fixed zero-mode convention and compact metric variation,

$$
\Gamma_\zeta(\mu_2)-\Gamma_\zeta(\mu_1)
=\frac{\log(\mu_2^2/\mu_1^2)}{32\pi^2}
\int\sqrt{g_E}\,a_4,
$$

where the equality is for the local scale dependence; no finite global
determinant on the noncompact manifold is assumed. Possible zero-mode
projector terms require a separate infrared convention and are not local
UV counterterms.

Thus, if the starting scale is $\mu_0=1/D_0$ and $\ell_*=D_0$, the
determined change is

$$
\Delta\Gamma_{\rm scale}
=b\int\sqrt{g_E}\,a_4,\qquad
b=\frac{\log2-2\gamma_E}{32\pi^2}
=-0.001460557999898805.
$$

The four-dimensional coefficient in135 can be rewritten as

$$
a_4=\frac{R^2}{30}-\frac{R_{ab}R^{ab}}{10}
-\frac7{360}\mathcal E_4-\frac{\Delta R}{30}
+\frac23F_{ab}F^{ab},
$$

with $\mathcal E_4=R_{abcd}R^{abcd}-4R_{ab}R^{ab}+R^2$ and
charge-included $F$. Euler and total-derivative terms have zero compact
integrated variation here. The corresponding Euclidean coefficients of
$R^2$, $\mathrm{Ric}^2$, $F^2$ are therefore
$b/30$, $-b/10$, $2b/3$ respectively.

Reusing only the order128 action derivatives in
[133's saved local responses](../results/dihole_local_source_results.json)
gives

$$
\frac{\delta\int\sqrt{g_E}a_4}{\mathcal T}
=\frac{7.706927314195081}{D_0},\qquad
\Delta E'_{\rm scale}
=\frac{-0.011256414343386234}{D_0}.
$$

Here $E=\Gamma_E/\mathcal T$. The saved Lorentzian action has
$E=-S_L/\mathcal T$, so its action coefficients have the opposite sign.
These numbers are algebraic reuse of existing floating quadrature, not
new enclosures or a vacuum-source result. They give only the scale part
of a convention change. Compensating gravitational/gauge coupling shifts
must accompany it; it cannot be added as physical support.

## What133 still leaves unspecified

The singular coefficients and scale do not themselves prescribe the
finite differentiated coincidence limit. One needs the particular
first-order spinor symmetrization, the off-shell Hilbert stress
bidifferential operator, and an explicit local correction pair
$(K_{ab},K^a)$ for stress and current.

Zahn's section3.1 supplies a current-preserving symmetrization, and
Proposition4.1 proves current conservation. Section4.2 explicitly
specializes its stress argument to vanishing gauge curvature. Neither
that argument nor the neutral Dirac correction in
[Dappiaggi--Hack--Pinamonti, equation34 and theorem4.1](https://arxiv.org/pdf/0904.0612)
is an implemented proof of the required charged identity

$$
\nabla^aT_{ab}=F_{ba}j^a.
$$

The exact unresolved freedom can be exhibited without any spectrum.
For the Lorentzian action $I_F=c_F\int\sqrt{-g}\,F^2$, define variations
by $\delta I_F=\int\sqrt{-g}(\tfrac12\Delta T^{ab}\delta g_{ab}
+\Delta j^a\delta A_a)$. Then

$$
\Delta T_{ab}=c_F(g_{ab}F^2-4F_{ac}F_b{}^c),\qquad
\Delta j^b=-4c_F\nabla_aF^{ab}.
$$

This pair obeys the joint Ward identity for every $c_F$. Variations of
$R^2$ and $\mathrm{Ric}^2$ supply the other two homogeneous freedoms.
They generally affect the paired null response. Pure metric contact
terms $g_{ab}C$, by contrast, have zero null contraction and zero paired
volume variation.

133's instruction to set the *remaining additions* to zero fixes these
freedoms only relative to a specified Ward-corrected reference. That
reference is not yet given as a formula. Its prohibition on freely added
smooth **bisolutions** does not settle this: permitted local parametrix
corrections need not be bisolutions. See Zahn, remark3.11 and equation35.
This does not invalidate133's intended normalization; it prevents assigning
its conversion constants from the present text alone.

## The finite next calculation

Fix the named parametrix and off-shell stress/current operators, compute
their local joint Ward defect from the charged Hadamard jets, and give
one explicit correcting pair. Compare that pair with the variational
heat/zeta prescription at $\mu_H$. Express the remaining difference as
the action variations of $R^2$, $\mathrm{Ric}^2$ and $F^2$; then133's
zero-addition rule has a definite reference and fixes the three constants.
This requires local coefficient algebra, not an exterior solve.

The following alternative fixes an operational reference without waiting
for that comparison. It does not modify133 or assign its charged contact
constants.

## Explicit alternative reference

For subsequent calculations, choose the parity-even, covariant physical
four-dimensional heat/zeta action at the fixed scale
$\mu_H=\sqrt2e^{-\gamma_E}/D_0$. This is a declared mathematical model,
not a proved identification with133 or a measurement of couplings. Keep
the same static spectral state: negative energies occupied and exact
zero modes empty. Keep $\mu_H$, proper regulator masses and the
renormalized cosmological and Einstein--Maxwell couplings fixed during
the compact variation. Add no further finite $R^2$, $R_{ab}R^{ab}$ or
charge-included $F^2$ action terms in this chosen scheme.

Let $\mathscr P_\epsilon$ be the squared physical Euclidean Dirac
operator, $m_j^2=jM^2$, $c=(1,-3,3,-1)$, and formally write

$$
\Gamma_{\rm PV}^{(M)}[g_\epsilon,A]
=-\frac12\sum_{j=0}^3c_j\log\det(\mathscr P_\epsilon+jM^2).
$$

Only relative compact variations, per unit Euclidean Killing time
$\mathcal T$, are required; this notation does not assert a finite
global determinant. Define the renormalized paired response by

$$
\boxed{\displaystyle
\mathcal Q_{\rm HK}[f]=\lim_{M\to\infty}
\left.\partial_\epsilon\right|_0\!\left\{
\frac{\Gamma_{\rm PV}^{(M)}[g_\epsilon,A]}{\mathcal T}
-\frac1{32\pi^2\mathcal T}\int\sqrt{g_{\epsilon,E}}\,
\bigl[\beta_0M^4a_0+\beta_1M^2a_2+\lambda_Ma_4\bigr]
\right\},}
$$
$$
\beta_0=\frac{9\log3-12\log2}{2},\qquad
\beta_1=6\log2-3\log3,\qquad
\lambda_M=\log\frac{M^2}{\mu_H^2}+\log\frac38
=\log(M^2D_0^2)+2\gamma_E+\log\frac3{16}.
$$

The $a_{2j}$ are the **physical four-dimensional** coefficients in135,
evaluated on $g_{\epsilon,E},A$; all three local terms are subtracted.
The paired family is $N_\epsilon=e^{\epsilon f}N$,
$h_{ss,\epsilon}=e^{-2\epsilon f}h_{ss}$, with the other spatial
coefficients and magnetic connection fixed. Its volume density
$N\sqrt h$ is unchanged, so the $a_0=4$ variation vanishes exactly.
Once the indicated limits and differentiation exist, $\Gamma_E/\mathcal T=E$
gives

$$
\mathcal Q_{\rm HK}[f]
=\int fN\sqrt h\,
\bigl(\langle T_{\hat0\hat0}\rangle_{\rm HK}
+\langle T_{\hat s\hat s}\rangle_{\rm HK}\bigr).
$$

The full vectorlike Dirac action, proper scalar PV masses and local
counterterms are gauge and diffeomorphism invariant. Compact background
transformations conjugate the regulated operators; differentiating the
invariant relative action therefore yields
$\nabla_a j^a=0$ and $\nabla^aT_{ab}=F_{ba}j^a$.
There is no vector or diffeomorphism anomaly here; a trace anomaly can
remain. See [Bastianelli--Martelli, section2.2](https://arxiv.org/pdf/1610.02304)
for the covariant Dirac PV mass and
[Bastianelli--Broccoli, section4.2.1](https://arxiv.org/pdf/1808.03489)
for the vector Ward identity. This continuum action argument neither
implements a charged point-split correction nor proves Ward identities
for a truncated angular or absorbing-boundary calculation.

This fixes the local ultraviolet normalization. Existence and
differentiability of the noncompact, massless relative action remain
separate infrared requirements. A squared determinant alone does not
settle a derivative at a zero crossing; no thermal prescription is
silently substituted for the declared zero-mode convention. At fixed
$M$, the angular, frequency, spatial and physical-exterior limits must
be controlled before $M\to\infty$, or justified by uniform estimates.
The scale-only number above is **not** added to the finite-$M$135 pilot:
the selected scale already occurs once in $\lambda_M$. Scheme changes
require compensating coupling changes, never a fit to the desired
supporting stress. No physical renormalized source has yet been computed
by making this definition.

## Executable local subtraction and fixed-mass cancellation

The [local counterterm calculator](../scripts/dihole_pv_counterterms.py)
reuses133's saved invariant variations; its [results](../results/dihole_pv_counterterms_results.json)
include an independent high-precision auxiliary Green-kernel scale check.
In units $D_0=1$, the local response to subtract is

$$
\mathcal Q_{\rm loc}(M)=0.00031154963588733757\,M^2
+0.0244023436787409\,[\log M^2+2\gamma_E+\log(3/16)].
$$

This is the large-mass local subtraction in the declared scheme, not an
approximation of the whole source at $MD_0=1$. It has not been applied to
the finite angular/absorbing pilot. Its role is to define the subtraction
once the physical fixed-$M$ limits are available.

A useful exact identity for the *same* positive physical Euclidean operator
$\mathscr P$ is

$$
\sum_{j=0}^3c_j(\mathscr P+jM^2)^{-1}
=\frac{6M^6}{\mathscr P(\mathscr P+M^2)
(\mathscr P+2M^2)(\mathscr P+3M^2)}.
$$

With a legitimate relative trace and the canonical Hilbert-space variation
included, cyclicity gives

$$
\delta\Gamma_{\rm PV}
=-3M^6\operatorname{Tr}\left[
\delta\mathscr P\,
\{\mathscr P(\mathscr P+M^2)(\mathscr P+2M^2)
(\mathscr P+3M^2)\}^{-1}\right].
$$

No commutation of $\delta\mathscr P$ and $\mathscr P$ is assumed. The
metric variation has differential order2 and the resolvent product order-8:
the localized ultraviolet operator has order-6, below the trace-class
threshold-4 in four dimensions. Equivalently, the localized variation heat
trace behaves as $O(s^{-3})$, while $(1-e^{-sM^2})^3=O(s^3)$, making the
short-proper-time variation integrand integrable.

This ultraviolet statement does not cure the $1/\lambda$ behavior at zero
spectral parameter. A common positive infrared shift makes the displayed
resolvent identity bounded; removing it needs the actual massless and
zero-mode analysis. Nor does the identity give an angular-tail error bar:
local elliptic and angular spectral-density bounds are still required.
The spherical angular projection does not commute with the actual curved
operator. Convergence of projections of a trace-class response does not,
by itself, bound a compressed-operator determinant calculation.

The four first-order absorbing-pilot determinants have not been shown to
arise from a common physical $\mathscr P+jM^2$ boundary/domain realization.
Replacing them by this rational formula without that derivation would
change the problem. Squaring the Hamiltonian also retains $m[H_0,N]$;
it is not the substitution $H_0^2+m^2I$.
