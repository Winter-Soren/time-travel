# Conforming trials across both physical mouth cuts

Increment138. This gives global trial sections for the common physical
Euclidean form in [137](PHYSICAL_EUCLIDEAN_OPERATOR.md), using the exact
geometry, spin lifts and magnetic bundle in [133](SEWN_SOURCE_STATE.md).
The exterior is the actual exterior of the two round cuts. Trial fields
vanish outside a compact support; no absorbing outer boundary is imposed.
This construction supplies admissible fields, not a source or remainder
certificate.

## Coordinates, frames and phase

Let $R_c=.05D$, $k=D/2$, $s_L=-.1D$ and $s_R=.1D$. The cut maps are

$$
x_R(\theta)=(R_c\sin\theta,k+R_c\cos\theta),\qquad
x_L(\theta)=(R_c\sin\theta,-k-R_c\cos\theta).
$$

Thus the lower global polar angle is $\pi-\theta$. With the rotating
exterior spin frame and the bridge's increasing-$s$ frame, use

$$
U_R=e^{-i\theta\sigma_y/2},\qquad U_L=e^{+i\theta\sigma_y/2},\qquad
\mathbb U_j=\operatorname{diag}(U_j,U_j).
$$

For the exterior flux amplitude $v_e=\sqrt P\,\Psi_\mu$, the corresponding
bridge flux trace, whose angular measure is $d\theta$, is

$$
\tau_j(\theta)=R_c\sqrt{\sin\theta}\,
\mathbb U_j^\dagger v_e(x_j(\theta)).
$$

The azimuthal factor is $e^{i\mu\varphi}/\sqrt{2\pi}$ with the fixed
global $\mu\in\mathbb Z+1/2$. Both Weyls have the same charge and gauge
patch. The lower orientation is handled by $x_L$ and $U_L$; do not also
reverse $q$ or $\mu$, or insert an additional modal parity into this
pointwise lift.

Fix an oriented handle phase by a charge-included flat connection
$A_s=\Theta'(s)$, with $\Theta(s_L)=0$ and $\Theta(s_R)=\phi$. In the gauge
where $A_s$ is retained, the traces are $\tau_L,\tau_R$. To reuse the
radial operator with $A_s=0$, set $v_b^0=e^{-i\Theta}v_b$; its traces are

$$
v_b^0(s_L,\theta)=\tau_L(\theta),\qquad
v_b^0(s_R,\theta)=e^{-i\phi}\tau_R(\theta).
$$

Use one convention throughout. The same constant phase multiplies both
Weyls. It is a fixed bundle transition, not charge conjugation.

## An explicit lift with analytic derivatives

Choose any four-component exterior test amplitude $v_e(\rho,z)$ that is
piecewise smooth and globally $C^1$ on the closure of the actual exterior,
with compact support away from $\rho=0$. Smooth test amplitudes are included.
For a directly
implementable family, use tensor Legendre polynomials times

$$
B(x,y)=[(1-x^2)(1-y^2)]^2\quad (|x|,|y|<1),
\qquad B=0\text{ otherwise},
$$

where $x,y$ affinely map a finite rectangle in $(\rho,z)$ to $[-1,1]^2$
and its lower radial edge is positive. Multiply by each constant spinor
component and restrict to the exterior of the two balls. The double
zeros make the zero extension $C^1$ and $H^1$; only first derivatives
are needed in the form. A rectangle crossing both mouths can describe
one field connecting them through the exterior.

Set $t=(s-s_L)/(s_R-s_L)$. The simple linear lift is

$$
\boxed{v_b^0(s,\theta)=(1-t)\tau_L(\theta)
+t e^{-i\phi}\tau_R(\theta).}
$$

It has the required two traces and explicit radial derivative. Smooth
cutoff blends could replace the linear weights, but are unnecessary
for form conformity. Add independent bridge test fields with zero
traces when enriching the space.

For implementation without $\sqrt{\sin\theta}$, write
$b=v_b^0/\sqrt{\sin\theta}$ and
$\widetilde\tau_j=R_c\mathbb U_j^\dagger v_e(x_j)$.
The same linear lift applies to $b$, now with measure $ds\,d(\cos\theta)$,
integrating the angular variable from $-1$ to1. Its Weyl angular operator is

$$
\widetilde{\mathcal K}
=-i\sigma_1(\partial_\theta+\tfrac12\cot\theta)
+\sigma_2\frac{NP(\mu-qa)}{\sin\theta}.
$$

The exact trace derivatives are, with
$\mathbb S_y=\operatorname{diag}(\sigma_y,\sigma_y)$,

$$
\widetilde\tau_R'
=R_c\mathbb U_R^\dagger
\left[R_c(\cos\theta\,\partial_\rho-\sin\theta\,\partial_z)
+\tfrac i2\mathbb S_y\right]v_e(x_R),
$$
$$
\widetilde\tau_L'
=R_c\mathbb U_L^\dagger
\left[R_c(\cos\theta\,\partial_\rho+\sin\theta\,\partial_z)
-\tfrac i2\mathbb S_y\right]v_e(x_L).
$$

For the $d\theta$ version, additionally differentiate
$\sqrt{\sin\theta}$, giving $\tau_j'/\sqrt{\sin\theta}
=\widetilde\tau_j'+\tfrac12\cot\theta\,\widetilde\tau_j$.
The exterior support condition makes both traces identically zero near
the poles. If its radial support starts at $\rho_{\min}<R_c$, the support
endpoints include $\cos\theta=\pm\sqrt{1-(\rho_{\min}/R_c)^2}$;
piecewise quadrature must respect these and any other support crossings.
The coordinate cotangent factors then introduce no axis singularity.

## Why the joined field is in the form domain

Recover the physical bridge spinor by dividing $v_b$ by
$R(s)\sqrt{P(s,\theta)\sin\theta}$. At either cut this agrees with the
rotated exterior physical spinor, including the declared phase. The
sections are piecewise $H^1$ on smooth regions and their physical traces
agree exactly. Integrating a weak first derivative over both regions
therefore cancels its interface term: the global section is $H^1$.
Normal derivatives need not agree. There is no delta-function term in
the first-order Dirac action.

This is also visible in the current. Write
$\Sigma=\operatorname{diag}(\sigma_3,-\sigma_3)$. The exterior outward
normal points into each removed ball. Its characteristic matrices pull
back to

$$
\mathbb U_R^\dagger\alpha_{n_e,R}\mathbb U_R=-\Sigma,
\qquad
\mathbb U_L^\dagger\alpha_{n_e,L}\mathbb U_L=+\Sigma.
$$

The exterior boundary measure in flux variables is
$\rho\,dl=R_c^2\sin\theta\,d\theta$, exactly canceled by the two trace
normalizations. The bridge Hamiltonian Green form is

$$
\mathcal B_b(v,w)=-i\int d\theta\,
[v^\dagger\Sigma w]_{s_L}^{s_R}.
$$

The exterior supplies its negative at both cuts. Unit phases cancel
between two fields in the same phase domain. The Euclidean Green
identity cancels for the same reason, since $\beta$ commutes with the
sewing matrices. Thus the joined trial belongs to the closed first-order
domain and to the norm-square form domain of $P_\kappa$ in137.
It need not belong to the second-order **operator** domain: applying
$P_\kappa$ pointwise and ignoring derivative jumps would be invalid.

Assemble the form and Euclidean norm by adding the two physical regions:

$$
q_\kappa=q_{\kappa,e}+q_{\kappa,b},\qquad
n=n_e+n_b,
$$

with exterior weights $\rho P$ and bridge weights $P$ in their respective
flux coordinates. The massive identity $q_{\kappa,m}=q_\kappa+m^2n$
then uses global, canceled boundary terms. The paired metric variation
has $f=0$ at the cuts, so137's identification $J_\epsilon=e^{\epsilon f/2}$
preserves both traces. Its $f'$ contribution must still be included.

## Enrichment and the limit still required

A useful sequence expands the exterior support, refines its bump/polynomial
space, includes every azimuthal sector, and adds a dense family of
zero-trace bridge fields. Exterior smooth-core approximation plus the lift
then leaves a bridge remainder with zero traces, which can be approximated
in $H^1_0$. The elementary lift is used for smooth traces here; no bounded
map from arbitrary $H^{1/2}$ traces to $H^1$ is claimed for this formula.

Support away from the axis is not an intrinsic form-density obstruction.
For a smooth compact section, a logarithmic cutoff between transverse
radii $\delta^2$ and $\delta$ has derivative energy
$\int_{\delta^2}^{\delta}|\chi'|^2\rho\,d\rho=1/|\log\delta|$.
In regular spin/gauge charts this tends to zero, so shrinking the excluded
axis neighborhood recovers such core sections in the form norm.

These observations specify how an exhaustive trial sequence can be built.
A fixed rectangle, axis exclusion, finite set of modes and linear bridge
lift are not that sequence. Finite matrix positivity, Green identities or
refinement alone do not bound its missing exterior, angular, frequency or
infrared contributions. A physical source sign requires those remainders
and136's local subtraction; it cannot be inferred from admissible trial
fields alone.

## Executed global pilot

The [implementation](../scripts/dihole_sewn_form.py) reuses137's physical
form and derivative on these global sections. Six scalar exterior bumps
(two radial by three axial polynomials), each with four spinor components,
give24 global trials per azimuthal sector. The exact round-domain
[quadrature helper](../scripts/dihole_sewn_quadrature.py) uses18,432 exterior
nodes at order24; the bridge uses6,912 nodes. Three sectors
$\mu=-1/2,1/2,3/2$ and phases $0,\pi$ are evaluated at $MD=1$.

The form and its derivative are exactly quadratic in real Fourier frequency.
Assembling their three matrix coefficients makes the33-frequency sweep
inexpensive. Cholesky whitening uses a fixed Euclidean norm. The response
is evaluated by the common-$P$ rational function from136, avoiding
subtraction of four nearly equal matrix inverses. The full four-component
trace is even in frequency; its $-1/2$ fibre response integrates with
$1/\pi$ on positive frequencies, with no further Weyl doubling.

All nine connected checks pass in71.35 seconds, including the six geometric
quadrature checks. The global Dirac skew error is below9.38e-13 in the
24-dimensional spaces and1.25e-12 in the enriched spaces. Proper-mass form
agreement is below1.09e-13. The metric derivative agrees with finite
differences to1.05e-10 relative error. Higher quadrature changes the form by
1.23e-12 and its derivative by8.19e-10; none is an outward error enclosure.

The original derivative check subtracted the large unchanged exterior form
at $\epsilon=\pm10^{-5}$, losing precision in the small compact change and
failing the2e-8 relative threshold (up to3.61e-7 observed). The exact
[original source](../provenance/dihole_sewn_form_pre_fd_local.py) and
[failed diagnostic](../results/dihole_sewn_form_pre_fd_local_results.json)
are preserved. Because $f=f'=0$ in the exterior, differencing only the
bridge tests the identical full derivative and removes this cancellation.
The tolerance and field equations were unchanged. Every saved physical
frequency response and basis-refinement value is byte-value equal to the
first run's JSON numbers.

The important physical-computation finding is **large basis sensitivity**.
On $\mu=1/2$, expanding the scalar polynomials to three by four gives48
nested global trials. The original24-dimensional form and norm are recovered
exactly as matrix subblocks under the same quadrature. The finite-window
responses over $\kappa D\in[.01,10]$ are:

| Phase | 24 global trials | 48 global trials | Ratio |
| --- | ---: | ---: | ---: |
| $0$ | $1.9059213387810255\times10^{-6}$ | $7.946254442718241\times10^{-6}$ | $4.16925$ |
| $\pi$ | $1.8319366359477733\times10^{-6}$ | $8.102834444067308\times10^{-6}$ | $4.42310$ |

Units are $\hbar c/D$. The enriched17-to33-frequency-grid sensitivity is
1.0043% and0.9775%, respectively. These positive values belong to finite
common-form matrices, not a controlled approximation of the renormalized
source. They cannot be combined numerically with136's absorbing pilot or
used to subtract136's large-$M$ local terms at a fixed basis.

See [results](../results/dihole_sewn_form_results.json) and
[figure](../figures/dihole_sewn_form_progress.png). The next enrichment must
include independent zero-trace bridge fields as well as these trace lifts;
otherwise the bridge field is still artificially determined by the exterior
polynomial traces. Use137's known compact bridge family, assembled together
with the lift on compatible split angular quadrature. Then assess the full
physical exterior/axis and spectral remainders. No physical source sign,
new law, stability result or backward passage has been established.
