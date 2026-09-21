# An expanding exterior with its energy source included

This follows the [rejected static completion](LOCALIZED_FIELD_HYPOTHESIS.md).
The change is explicit: add a positive constant vacuum-energy density and
derive both the background and its coupled perturbations again. Expansion
limiting ghost-condensate growth is an established idea; see
[Arkani-Hamed et al., section 8](https://arxiv.org/pdf/hep-th/0312099).
The exact coefficients below are derived for our displayed truncated action.
No new law, novel mechanism, physical material or time machine is claimed.

## Action and background

In natural units, take

$$
S=\int d^4x\sqrt{-g}\left[
\frac{P}{2}R+\frac{\Lambda^4}{4}(X-1)^2-\rho_\Lambda
-\alpha\Lambda^2X\Theta^2\right],\quad
P=\bar M_{\rm Pl}^2,\quad \alpha>0,
$$

where $X=-(\partial\phi)^2/(2\Lambda^4)$ and
$\Theta=\nabla_\mu[-\partial^\mu\phi/\sqrt{2\Lambda^4X}]$.
For the spatially flat background $a(t)=e^{Ht}$,
$\phi=\sqrt{2\Lambda^4X_0}\,t$, the actual background equations give

$$
X_0=1+18\alpha H^2/\Lambda^2,\qquad
\rho_\Lambda=(3P+9\alpha\Lambda^2)H^2+243\alpha^2H^4.
$$

The added energy is a required source, not a freely available laboratory
resource. Removing it returns the earlier obstruction to constant expansion.

## Scalar constraints and the term that changes in an expanding metric

Use unitary gauge $N=1+n$, $\gamma_{ij}=a^2e^{2\zeta}\delta_{ij}$ and scalar
shift $N^i$. Define

$$c=\alpha\Lambda^2X_0,\quad V=\Lambda^4X_0^2,\quad
Z=cV+H^2(P+3c)(P+12c),\quad d=\Delta\zeta/a^2.$$

The lapse and shift constraints, for a nonzero scalar Fourier mode, are

$$
n=\frac{P[H(P+9c)\dot\zeta+c\,d]}{Z},\qquad
\partial_iN^i=\frac{P+3c}{c}\dot\zeta-\frac{H(P+9c)}{c}n.
$$

After eliminating them, the action contains
$-2J\dot\zeta\Delta\zeta/a^2$, with
$J=P^2H(P+9c)/Z$. This is not just an endpoint term: with the expanding
integration measure,

$$\int dt\,d^3x\,a\dot\zeta\Delta\zeta
=\frac H2\int dt\,d^3x\,a|\nabla\zeta|^2$$

up to endpoints. The reduced action is therefore

$$
S_2=\int dt\,d^3x\,a^3\left[
K\dot\zeta^2+G\frac{|\nabla\zeta|^2}{a^2}
-W\frac{(\Delta\zeta)^2}{a^4}\right],
$$

$$
K=\frac{P[(P+3c)V+108c^2H^2]}{Z},\quad
G=\frac{Pc[V+6H^2(P+6c)]}{Z},\quad W=\frac{P^2c}{Z}.
$$

All three coefficients are positive. The positive sign of $G$ in the action
is a negative contribution to squared frequency, so it does not establish
positive spatial energy. The transverse traceless tensor action remains
$P\int a^3[\dot h_{ij}^2-a^{-2}(\partial_kh_{ij})^2]/8$.
For $\gamma_{ij}=a^2(e^h)_{ij}$, $\operatorname{tr}h=0$, its determinant is
unchanged and $\Theta=3H$ at unperturbed lapse and shift; the added operator
provides no tensor kinetic or gradient term. This is a quadratic statement.

The executable derives the scalar quadratic action from an exact geometric
ADM ansatz before using the background equations or eliminating constraints.
An independent three-dimensional derivation agrees. Limits $H\to0$,
$c\to0$, and $P\to\infty$ supply distinct checks; the last agrees with an
independent fixed-metric scalar expansion.

## What expansion controls

Each mode obeys

$$
\ddot\zeta_k+3H\dot\zeta_k+
A k_{\rm phys}^2(k_{\rm phys}^2-k_J^2)\zeta_k=0,
$$

$$
A=\frac{Pc}{(P+3c)V+108c^2H^2},\qquad
k_J^2=\frac{V+6H^2(P+6c)}{P},\qquad k_{\rm phys}=k/a(t).
$$

The coefficients are now derived from the expanding action. They are not
the static coefficients with damping added. The instantaneous negative band
persists. Let a mode enter it at $t=0$, with $\zeta(0)=1$ and
$\dot\zeta(0)=0$, and define $r=\sqrt A k_J^2/(2H)$. The positive Volterra
kernel from the preceding report proves

$$1+r^2/3\le\zeta(\infty)\le e^{r^2/3}.$$

For arbitrary finite band-entry data, the same absolute-value argument gives
$|\zeta(t)|\le(|\zeta(0)|+|\dot\zeta(0)|/(3H))e^{r^2/3}$.
These are linear mode bounds from band entry onwards. They do not certify
nonlinear evolution, arbitrary ultraviolet initial data, or a finite apparatus.

At the illustrative $\Lambda/\bar M_{\rm Pl}=0.01$, $\alpha=1$, minimizing
$r$ numerically within the stated scan interval gives

| Quantity | Value |
| --- | ---: |
| $H/\Lambda$ | 0.00408125816 |
| $X_0$ | 1.00029982003 |
| $k_J/\Lambda$ | 0.01414425551 |
| $r$ | 0.02450224259 |
| Analytic fractional-amplification bound, evaluated numerically | 0.00020013999 (0.020014%) |
| Integrated fractional amplification | 0.00020013312 |
| Required $\rho_\Lambda/\Lambda^4$ | 0.49985002270 |

This is a numerical optimum within the tested parameter family, not a
certified global optimization or a measured material scale. The cutoff ratios
are illustrative; small frequencies relative to $\Lambda$ alone do not prove
the true strong-coupling cutoff. At fixed cutoff, increasing expansion without
limit does not improve the bound, because the dispersion changes too.

All 15 executable checks pass. Two independent ODE resolutions differ by at
most $2.7\times10^{-15}$ in the reported late-mode estimates. These numerical
values are not interval enclosures; the Volterra inequality is analytic.
See [source](../scripts/expanding_field_completion.py),
[saved results](../results/expanding_field_completion_results.json), and
[figure](../figures/expanding_field_completion_progress.png).

## First total-null-stress calculation

For a scalar disturbance depending on one spatial coordinate, define
$q=\dot\zeta-Hn$, $s=\partial_xN^x$, $d=\partial_x^2\zeta/a^2$.
In the physical orthonormal frame of the foliation, with future normal $u$
and longitudinal unit vector $e$, the linear Ricci projections are

$$
\delta R_{uu}+\delta R_{ee}=-2(\dot q+Hs+d),\qquad
\delta R_{ue}=-2\partial_xq/a.
$$

These were independently checked by expanding the full four-dimensional
Christoffels and Ricci tensor, including the frame perturbation; all three
individual projection residuals vanished exactly. These are additional
analytic checks, separate from the saved executable's 15 checks.
Using $k_\pm=u\pm e$ with $-u\cdot k_\pm=1$, the total on-shell source is
$T_\pm=P\delta R_{\mu\nu}k_\pm^\mu k_\pm^\nu$ to first order. The constraints
give $n=(J\dot\zeta+Wd)/P$, $s=(K\dot\zeta-Jd)/P$, so

$$
T_\pm=-2\left[G\ddot\zeta+HK\dot\zeta-HW\dot d+Gd
\ \pm\frac2a(G\partial_x\dot\zeta-HW\partial_xd)\right].
$$

At band entry, let $\zeta=\epsilon\cos(k_{\rm com}x)$,
$\dot\zeta=0$ and $k_{\rm com}/a=k_J$. The mode equation gives
$\ddot\zeta=0$, but the physical Laplacian still evolves:
$\dot d=-2Hd$. Hence

$$
T_\pm=2\epsilon k_J^2\left[
(G+2H^2W)\cos(k_{\rm com}x)
\pm2HWk_J\sin(k_{\rm com}x)\right].
$$

Both directions are negative at a cosine trough for $\epsilon>0$. At the
illustrative optimum in the table,

$$T_\pm/\Lambda^4\simeq\epsilon[0.00559592136923\cos(k_{\rm com}x)
\pm0.00276837291271\sin(k_{\rm com}x)].$$

For $\epsilon=0.001$, the trough has $T_\pm/\Lambda^4\simeq-5.596\times10^{-6}$.
The associated lapse perturbation is about $1.199\times10^{-6}$. This is
total linear null stress, not an inference from $\mathcal P_X$ alone:
the local first-order clock value remains $X\simeq1.00029742>1$.
The disturbance is sinusoidal, with compensating positive regions and no
finite transverse extent. It establishes a linear sign-changing source
response, not a localized supporting material or a nonlinear bridge solution.

## What this does not supply

The total homogeneous source obeys $\rho+p=-2P\dot H=0$ on this background.
Its two-derivative part alone gives
$\Lambda^4X_0(X_0-1)>0$; the higher-derivative background stress must be
included. A homogeneous expanding exterior is therefore not the required
localized negative-null-energy source. The curvature disturbance used above
is a field mode, not an abiotic object passing through a medium.

The next physical test is to construct on-shell localized disturbances and
calculate their total null stress, flux and evolution from the same action.
A prescribed $P_X$ profile is insufficient. Matching to the proposed bridge,
backreaction, a realizable energy source, nonlinear stability and chronology
change remain unresolved. Travel from 13 September 2027 to 12 September 2027
is still unproved.
