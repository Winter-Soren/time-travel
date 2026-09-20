# A localized-field hypothesis and its first rejection tests

This is a separate speculative matter branch, started now instead of waiting
for every refinement of the Dirac calculation. It tests whether a localized
negative gradient coefficient can avoid an immediate linear instability.
It does **not** supply a physical material, a new law, or a covariant source
for the existing bridge. No novelty claim is made for these mathematical tools.

The motivation is the published ghost-condensate effective-field framework:
a higher spatial derivative can provide a positive fourth-power dispersion
term. Gravity and nonlinear evolution remain essential additional tests.
See [Arkani-Hamed et al., equations 2.17–18 and 7.8](https://arxiv.org/pdf/hep-th/0312099),
[Creminelli et al., equation 41 and its cutoff qualification](https://arxiv.org/pdf/hep-th/0606090),
and the [nonlinear caustic discussion, section 3.5](https://arxiv.org/pdf/hep-ph/0507120).
These papers do not construct the localized profile assumed below.

## Defined hypothesis

Use natural units and a real perturbation $\pi$ with prescribed coefficients:

$$
S_2=\frac12\int dt\,d^3x\,
\left[K\dot\pi^2-p_g(\mathbf x)|\nabla\pi|^2-a(\nabla^2\pi)^2\right],
\qquad K>0,\quad a>0.
$$

$p_g$ is a gradient coefficient, **not** the background pressure. First take
a planar comparison profile $p_g=-\sigma<0$ for $|x|<L$ and
$p_g=p_{\rm out}>0$ outside. The domain is all of space; there are no reflecting
walls. An infinite homogeneous negative region has
$\omega^2=(-\sigma k^2+a k^4)/K$ and is unstable at sufficiently small $k$.
Localization changes that conclusion under the following explicit condition.

## Exact slab criterion, including transverse modes

Define the auxiliary one-dimensional operator

$$T=-a\partial_x^2+p_g(x).$$

Its first zero-energy even solution is proportional to $\cos(kx)$ inside
and a decaying exponential outside. Continuity of the function and its first
derivative gives $k\tan(kL)=\kappa$, where
$k=\sqrt{\sigma/a}$ and $\kappa=\sqrt{p_{\rm out}/a}$. The ground level
therefore crosses zero at

$$
\boxed{L_{\rm crit}=\sqrt{\frac a\sigma}
\arctan\sqrt{\frac{p_{\rm out}}\sigma}}.
$$

$T\ge0$ precisely when $L\le L_{\rm crit}$. For transverse momentum $q$,
integration by parts gives the full spatial quadratic form

$$
E_q=\langle\pi',T\pi'\rangle
+q^2\left[a\|\pi'\|^2+\langle\pi,T\pi\rangle\right]
+a q^4\|\pi\|^2.
$$

Thus the same condition controls every transverse momentum. Conversely, a
negative direction of $T$ can be approximated by a derivative of a compact
function: subtract a distant broad bump to enforce zero integral, with
arbitrarily small added form cost. This produces a negative $q=0$ direction.
At equality the statement is nonnegativity, not a positive spectral gap or
nonlinear stability. The auxiliary eigenvalues are not physical growth rates.

A smooth compact three-dimensional profile satisfying
$p_{\rm blob}(\mathbf x)\ge p_{\rm slab}(x)$ inherits nonnegativity, since
the difference in spatial energy is
$\int(p_{\rm blob}-p_{\rm slab})|\nabla\pi|^2\ge0$.
This avoids treating the discontinuous comparison profile as a physical wall.
If the step is used directly, its form domain is $H^2$; the physical
fourth-order interface conditions involve $\pi,\pi',a\pi''$ and
$p_g\pi'-a\pi'''$, not the auxiliary matching conditions applied to $\pi$.

## A compact blob without a positive exterior coefficient

A positive clock-field coefficient outside would itself need a justified
background and stress. A separate conservative comparison avoids assuming
that exterior for free. Suppose $p_g\ge-\sigma$ inside a ball of radius $R$
and $p_g\ge0$ outside. The componentwise three-dimensional Hardy inequality
and the integrated Hessian identity imply

$$
\int_{B_R}|\nabla\pi|^2
\le R^2\int\frac{|\nabla\pi|^2}{r^2}
\le4R^2\int|\nabla\nabla\pi|^2
=4R^2\int|\Delta\pi|^2.
$$

For smooth compact tests, Hardy's constant follows by integrating
$\nabla\cdot(\mathbf x/r^2)=1/r^2$ against each squared derivative and
applying Cauchy–Schwarz; the result extends by form closure. Hence

$$E[\pi]\ge(a-4\sigma R^2)\|\Delta\pi\|^2.$$

The sufficient condition $R\le\tfrac12\sqrt{a/\sigma}$ works even when
$p_g=0$ outside. It is **not necessary**; exceeding this conservative radius
does not prove instability. This is genuinely a three-dimensional comparison,
so the vanishing positive-exterior slab window must not be used to reject it.

## Checked examples and conditional resource scale

[The executable calculation](../scripts/localized_field_hypothesis.py) passes
13 checks, including independent ground-level matching, transverse-form
integrals, the three-dimensional Hessian/Hardy control, and the gravitational
rejection test below. Passing these checks does not mean the medium is viable.
For $a=\alpha/\Lambda^2$, $\alpha=K=p_{\rm out}=1$ and $\sigma=10^{-4}$,
$L_{\rm crit}\Lambda=156.079666$. An illustrative smooth slab-comparable blob
can have outer half-width $100/\Lambda$, core $80/\Lambda$, and transition
width $20/\Lambda$. The separate zero-exterior sufficient radius is
$50/\Lambda$; an example radius $40/\Lambda$ retains a bending margin of 0.36.
These are coefficient hierarchies, not measured properties or verified EFT
strong-coupling scales. See [results](../results/localized_field_hypothesis_results.json)
and [visual](../figures/localized_field_hypothesis_progress.png).

Only if a complete background realizes the normalization
$\mathcal L=\Lambda^4\mathcal P(X)$,
$X=-(\partial\phi)^2/(2\Lambda^4)>0$, does its two-derivative contribution obey

$$
p_g=\mathcal P_X,\quad K=\mathcal P_X+2X\mathcal P_{XX},\quad
\rho+P_{\rm phys}=2\Lambda^4X\mathcal P_X.
$$

Higher-derivative terms may add background stress. Within the exact slab
window, the conditional magnitude obeys

$$
|\rho+P_{\rm phys}|<\frac{\pi^2}{2}\,
\frac{\alpha X\Lambda^2}{L^2}.
$$

If a specified geometry required $B\bar M_{\rm Pl}^2/L^2$ from this contribution,
it would require $\Lambda/\bar M_{\rm Pl}>\sqrt{2B/(\pi^2\alpha X)}$.
For $B=\alpha=X=1$ this is about 0.45016. This is a conditional scale tension,
not a gravitational no-go theorem; $B$ has not been matched to our bridge.
Do not apply the exact slab bound as a necessary condition for all compact
blobs, or turn the sufficient Hardy criterion into a necessary cutoff bound.

## A concrete completion fails its gravitational test

We tested the minimal covariant candidate

$$
S=\int\sqrt{-g}\,d^4x\left[
\frac{\bar M_{\rm Pl}^2}{2}R+
\frac{\Lambda^4}{4}(X-1)^2-\alpha\Lambda^2X\Theta^2\right],
\quad u_\mu=-\frac{\partial_\mu\phi}{\sqrt{2\Lambda^4X}},
\quad\Theta=\nabla_\mu u^\mu.
$$

For the zero-shift static clock $\phi=\sqrt2\Lambda^2t$,
$X=N^{-2}$ and $\Theta=0$. The last term then has zero background stress;
the remaining coefficients are $p_g=(X-1)/2$, $K=(3X-1)/2$.
Healthy time kinetics and negative NEC contribution require
$1/3<X<1$, or $1<N<\sqrt3$. This condition does not solve Einstein's
equations for such a lapse profile. The exterior $X=1$ has zero background
stress and the decoupled spatial coefficient zero.

That exterior is nevertheless unstable once gravity is included. Write
$P=\bar M_{\rm Pl}^2$, $c=\alpha\Lambda^2$, scalar spatial metric
$\gamma_{ij}=e^{2\zeta}\delta_{ij}$, lapse $1+n$ and shift $\partial_i\beta$.
The actual quadratic action has the constraints

$$n=\frac{P}{\Lambda^4}\Delta\zeta,\qquad
\Delta\beta=(3+P/c)\dot\zeta.$$

Symbolic elimination in the executable reproduces

$$
\mathcal L_{\rm red}=(3P+P^2/c)\dot\zeta^2
+P|\nabla\zeta|^2-\frac{P^2}{\Lambda^4}(\Delta\zeta)^2,
$$

and therefore

$$
\boxed{\omega^2=
\frac{\alpha/\Lambda^2}{1+3\alpha\Lambda^2/\bar M_{\rm Pl}^2}
k^2\left(k^2-\frac{\Lambda^4}{\bar M_{\rm Pl}^2}\right)}.
$$

The band $0<k<\Lambda^2/\bar M_{\rm Pl}$ is unstable. This agrees in the
weak-gravity limit with the published Jeans mechanism, with paper parameters
$M^4=2\Lambda^4$, $\bar M^2=2\alpha\Lambda^2$; the retained denominator above
comes from our explicit action and constraint calculation. At the illustrative
$\Lambda/\bar M_{\rm Pl}=0.01$, $\alpha=1$, the maximal growth rate is about
$4.99925\times10^{-5}\Lambda$. It can be slow; this is not an immediate-collapse
claim or a completed finite-operation analysis.

**Decision: reject this particular action with an infinite static ghost-vacuum
exterior as a globally linearly stable medium.** Making the negative region
smaller does not remove an instability already present far outside it.
The result is not an exclusion of all NEC-violating theories.

## Why tuning the same static constraint terms cannot fix this

This is our algebraic test of a defined class, not a universal no-go theorem.
For a constant background, let $v=\dot\zeta$, $d=\Delta\zeta$ and let $y$
collect nondynamical scalar constraints. Consider

$$
\mathcal L=A v^2+2v u^Ty+y^TM(k)y+2d w^Ty+D|\nabla\zeta|^2.
$$

Assume real, time-independent, even spatial coefficients,
$M(k)=M_0+O(k^2)$ with invertible $M_0$, and finite couplings. Eliminating
$y$ gives kinetic coefficient $K(k)=A-u^TM^{-1}u$ and fourth-derivative
coefficient $W(k)=w^TM^{-1}w$. The mixed term
$-2u^TM^{-1}w\,\dot\zeta\Delta\zeta$ is a time boundary term. Consequently,

$$
\omega^2=\frac{-Dk^2+W(k)k^4}{K(k)}
=-\frac{D}{K_0}k^2+O(k^4).
$$

With $D>0$ and healthy $K_0>0$, regular finite changes to those constraints
cannot remove the longest unstable wavelengths. If
$K(k)=K_0+K_2k^2+\cdots$, the quartic dispersion coefficient is
$W_0/K_0+DK_2/K_0^2$; its sign still cannot overcome the leading term near zero.
Negative $K_0$ introduces a ghost. Singular $M_0$, vanishing $K_0$, extra
propagating fields, a different leading gradient term, finite boundaries,
and time-dependent or curved backgrounds lie outside this test.

## Changing the background: finite growth, then another obstruction

Expansion changes an actual assumption of the static test. To isolate what
redshifting could accomplish, consider the explicitly hypothetical equation

$$
\ddot z+3H\dot z+A k(t)^2[k(t)^2-k_J^2]z=0,\quad
k(t)=k_J e^{-Ht},\quad H,A>0,
$$

with $z(0)=1$, $\dot z(0)=0$. The constant coefficients here have **not** been
derived for an expanding completion of the action. Set
$\tau=Ht$, $r=\sqrt A k_J^2/(2H)$. Then
$z''+3z'=b(\tau)z$ with
$b=4r^2(e^{-2\tau}-e^{-4\tau})\ge0$, and the exact integral equation is

$$
z(\tau)=1+\frac13\int_0^\tau[1-e^{-3(\tau-s)}]b(s)z(s)\,ds.
$$

Positivity implies $z\ge1$, $z'\ge0$. Since $\int_0^\infty b=r^2$,
Gronwall bounds the monotone limit:

$$\boxed{1+r^2/3\le z(\infty)\le e^{r^2/3}}.$$

Thus this toy has finite amplification, not an absence of amplification or
positive instantaneous energy. For example, $H\ge\sqrt A k_J^2/2$ bounds
this initial mode's amplification by $e^{1/3}<1.396$. An independent DOP853
calculation gives $z(\infty)\simeq1.371885020551$ at $r=1$; this numerical
value is not an interval certificate. The analytic bound is the result used
here. It does not assert a uniform physical-mode bound for a completed medium.

The same bare action fails a prerequisite even before applying that toy.
On a homogeneous clock and flat Friedmann background, let
$H=\dot a/(aN)$, $X=\dot\phi^2/(2\Lambda^4N^2)$. The homogeneous Lagrangian is

$$
L=a^3N\left[-3\bar M_{\rm Pl}^2H^2+
\Lambda^4\mathcal P(X)-9\alpha\Lambda^2XH^2\right].
$$

For constant $H>0$, constant $X>0$, the lapse equation and conserved
shift current respectively give

$$
(3\bar M_{\rm Pl}^2+27\alpha\Lambda^2X)H^2
=\Lambda^4(2X\mathcal P_X-\mathcal P),\qquad
\mathcal P_X=9\alpha H^2/\Lambda^2.
$$

Combining them yields

$$
(3\bar M_{\rm Pl}^2+9\alpha\Lambda^2X)H^2=-\Lambda^4\mathcal P.
$$

The independent scale-factor equation gives the same condition. With the
bare $\mathcal P=(X-1)^2/4\ge0$ and $\alpha>0$, it forbids this expanding
constant-clock background. These homogeneous variations were checked
independently with symbolic algebra; they are additional analytic checks,
separate from the 13 checks in the saved executable result.

**Decision:** the damping toy cannot rescue this unchanged action. A positive
constant potential, represented by $\mathcal P\mapsto\mathcal P-v$, is one
explicit possible change, with an additional energy density $v\Lambda^4$.
For this changed action, a constant solution would require

$$
X=1+18\alpha H^2/\Lambda^2,\qquad
v\Lambda^4=(3\bar M_{\rm Pl}^2+9\alpha\Lambda^2)H^2+243\alpha^2H^4.
$$

Its coupled perturbations must be derived anew; neither that energy source
nor a localized boundary can be assumed available. This exterior condition
does not supply the sought localized NEC-violating region.

## Remaining physical requirement

The missing ingredient is a background that generates $p_g(\mathbf x)$.
Even flat stationary pure $\mathcal P(X)$ with
$\phi=ct+\psi(x)$ and vanishing asymptotic spatial current gives
$\mathcal P_X\psi'=0$. It cannot freely supply the negative-inside,
positive-outside profile. Higher derivatives, gravity or additional fields
must be included in their actual background equations and stress.

The regular static constraint class above does not remove the exterior
instability. A background or field-content change must be stated, not silently
assumed. Background equations and boundary stress must be checked before
another coefficient sweep. A generic
$(\Box\phi)^2$ term cannot be interpreted as an exact all-frequency theory
while silently discarding higher time derivatives. Nonlinear stability,
causality, a usable cutoff, material generation and the full Einstein source
remain open. Neither this branch nor the Dirac calculation establishes a
stable medium or travel from 13 September 2027 to 12 September 2027.
