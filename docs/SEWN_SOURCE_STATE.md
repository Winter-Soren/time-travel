# Quantum state and source on the smooth sewn bridge

Increment133. This extends the [uncut domain proof](DIHOLE_DIRAC_DOMAIN.md)
and [stress prescription](EXTERIOR_STRESS_BRIDGE.md) to the exact smooth
[129 bridge](../scripts/dihole_stretched_collar.py). The construction below
defines a complete static test-field background, a self-adjoint charged Dirac
Hamiltonian and a Hadamard ground state. It does **not** solve its supporting
Einstein--Maxwell equations or identify the absorbing numerical boundary with
the physical exterior.

## Exact geometry and smooth sewing

Use the linked exact balanced parameters of118, not independently rounded
parameters: $D=1$, $k=D/2$, $Q=D/5$,
$M=(\sqrt{k^2+4Qk}-k)/2$, $a_g=\sqrt{k^2-M^2}$ and
$B=M/[(k+a_g)(M+k)]$. Remove the two coordinate balls of radius
$R_c=.05D$ around $(0,0,\pm k)$ and attach $S^2\times[-s_c,s_c]$.
Put $R_0=.025D$, $\delta=2(R_c-R_0)$ and a central coordinate length
$L_c=.1D$, so $s_c=L_c/2+\delta$. On each collar,

$$
u=|s|-L_c/2,\qquad
R(s)=R_0+\delta\int_0^{u/\delta}S(t)\,dt,
\quad S(t)=\frac{e^{-1/t}}{e^{-1/t}+e^{-1/(1-t)}}.
$$

Extend $S$ by0 and1 at its flat endpoints; set $R=R_0$ centrally.
Use the actual upper-mouth fields evaluated at
$(\rho,z)=(R\sin\theta,k+R\cos\theta)$:

$$
g=-N^2dt^2+P^2(ds^2+R^2d\theta^2)
  +\frac{R^2\sin^2\theta}{N^2}d\varphi^2.
$$

At the right cut this equals the exterior germ to every derivative:
$R=R_c$, $R'=1$, and all higher derivatives vanish. At the left cut use
$(\rho,z)=(R\sin\theta,-k-R\cos\theta)$ and the corresponding reversed
normal. Reflection symmetry of the exact dihole gives the same $N,P,a$.
The flat jets also make both cylinder joins smooth. This is a prescribed
metric, **not** the degenerate coordinate pullback with radial factor
$P^2R'^2$.

At both poles $NP=1$, while the coefficients are smooth in $\theta^2$ or
$(\pi-\theta)^2$. Consequently the angular circumference/radius ratio
tends to $2\pi$ and the metric extends smoothly over the axes. They are
included points, not boundaries. The bridge is compact, $R\ge R_0>0$, and
$N,P$ are smooth and positive there. The **whole sewn spatial manifold**,
including the exterior, has topology $\mathbb R^3\#(S^2\times S^1)$;
the compact bridge alone is $S^2\times[-s_c,s_c]$.

## Magnetic bundle, spin structure and loop phase

For integer $q$, extend the charge-included connection by
$A^o=q\,a(R,\theta)d\varphi$. Its polar values are $a=0,1$; the second
regular patch is

$$
A^c=q(a-1)d\varphi,\qquad \Psi^c=e^{-iq\varphi}\Psi^o.
$$

The vanishing coefficients are quadratic at their respective poles. The
transition is single-valued and the cross-sectional flux is $2\pi q$.
Opposite outward mouth fluxes therefore join without changing charge or
introducing a magnetic sheet. All connection jets match at the cuts.
The generally nonzero Maxwell current inside the prescribed collars is
smooth; the free Dirac construction does not require Maxwell's source-free
equation.

Fix the exterior spin structure and choose the bridge-frame transitions,
away from the polar coordinate singularities, as

$$
U_R=e^{-i\theta\sigma_y/2},\qquad U_L=e^{+i\theta\sigma_y/2}.
$$

These proper frame rotations, completed by the regular polar spin charts,
specify one spin structure. The unique spin structure on each overlap
$S^2$ permits either lift; reversing one endpoint's overall sign selects
the other handle spin structure. Keep the displayed choice fixed.

An additional phase $\phi$ is implemented by a flat connection
$\phi\eta(s)ds$, with smooth compact support inside the neck and
$\int\eta(s)ds=1$. It adds no curvature. Gauging it away within the neck
leaves a constant endpoint transition and the opposite phases on the two
transmission directions used in131. Thus $\phi=0,\pi$ are explicit
global alternatives, not arbitrary changes of modal basis. The fixed
lower-frame parity $(-1)^{n+1}$ is not an additional adjustable holonomy.

Axial rotations preserve this connection and the spin sewing. The118 lift
still gives global $\mu\in\mathbb Z+1/2$, with local exponents $\mu$ and
$\mu-q$ in the two gauge patches. Its Fourier projections reduce the
full operator. One retained $\mu$ is not the complete field.

## Optical completeness and Hamiltonian

Set $\widehat h=N^{-2}h$. Outside a compact set this is the original
optical metric, whose meridional factor obeys $P/N\ge1$. An escaping path
therefore has infinite optical length. The removed horizon ends have been
replaced by a smooth compact positive-metric region with no boundary;
there are no other escaping ends. Hence $(\Sigma,\widehat h)$ is complete.
The ultrastatic spacetime $-dt^2+\widehat h$ is globally hyperbolic, and
its conformal metric $g$ has the same causal structure. This supplies a
global time function, not a chronology-changing spacetime.

For each fixed spin/phase choice, the smooth unitary magnetic connection
twists a Dirac bundle on this complete manifold. Its symmetric elliptic
operator $D_A$ is essentially self-adjoint on $C_c^\infty$ by
[Bandara--Saratchandran, section3.1 and theorem3.5](https://arxiv.org/pdf/1606.06190).
Equivalently, complete-metric cutoffs with $|d\chi_R|\le C/R$ eliminate
deficiency vectors using $[D_A,\chi_R]=-i\sigma(d\chi_R)$. Its domain is
the maximal distributional graph domain $\{\psi\in L^2:D_A\psi\in L^2\}$;
no cut or axis boundary condition remains.

The massless conformal map and canonical amplitudes are

$$
\widehat\Psi=N^{3/2}\Psi,\qquad
\|\widehat\Psi\|_{L^2(\widehat h)}=\|\Psi\|_{L^2(h)},
$$
$$
u_e=P\sqrt{\rho/N}\,\Psi_\mu,\qquad
u_b=PR\sqrt{\sin\theta/N}\,\Psi_\mu,\qquad
w=\sqrt{N/P}\,u_b.
$$

Here the Fourier factor is $e^{i\mu\varphi}/\sqrt{2\pi}$; $u_e,u_b$
have measures $d\rho dz,dsd\theta$, respectively. At either cut,
$w=\sqrt{R_cN/P}\,U_{R,L}^\dagger u_e$, up to the declared flat phase.
The exterior inner-normal flux matrices become $-\sigma_3$ on the right
and $+\sigma_3$ on the left, so seam currents cancel. The full massless
Dirac Hamiltonian is $H_D=D_A\oplus(-D_A)$ with the same charge and sewing.
These are exact continuum identifications, not claims about polygonal
facets or truncated modal Gram matrices.

## A definite state

Choose the gauge-invariant quasifree ground state with occupation
$C_-=\mathbf1_{(-\infty,0)}(H_D)$ and $C_+=I-C_-$: **all exact zero modes
are empty**. Its canonical covariances are
$e^{-i(t-t')H_D}C_\pm$. This defines a state without a spectral gap.

The smooth static connection has no time component. The charged
ultrastatic spectral construction therefore applies directly; see
[Zahn, section3 and footnote12](https://arxiv.org/pdf/1210.4031).
Alternatively, complete timelike translations, the ground-state spectrum,
the field equation and CAR locality meet
[Sahlmann--Verch theorem5.1](https://arxiv.org/pdf/math-ph/0002021);
their [Dirac theorem5.8](https://arxiv.org/pdf/math-ph/0008029)
identifies the resulting microlocal condition with Hadamard form.
For the charged field one includes the conjugate bundle in the CAR
description. Smooth positive conformal transfer preserves that local
property. Changing the bounded occupation on $\ker H_D$ adds a locally
smooth kernel: $H_{D,x}^2+H_{D,y}^2$ is elliptic and annihilates it.
Thus the stated zero-mode convention is admissible. These conclusions
hold everywhere on the sewn manifold; there is no horizon to extend across.

## The source that132 has not yet computed

Choose the following baseline mathematical model, not measured physical
couplings. Use the locally covariant charged-Dirac Hadamard prescription,
with subtraction length $\ell_*=D_0$ fixed at the baseline and held fixed
throughout the $\epsilon$ variation. The logarithmic parametrix convention
is $\log(\sigma_+/\ell_*^2)$, with its standard charged Hadamard
coefficients and no freely added smooth bisolution. Here $\sigma_+$ is
Synge's half-squared geodesic distance with the causal boundary prescription.

Normalize the current and stress jointly with the required local
Ward-identity correction:
$\nabla_a\langle j^a\rangle=0$ and
$\nabla^a\langle T_{ab}\rangle=F_{ba}\langle j^a\rangle$.
Matter stress is not separately conserved in this magnetic background.
The explicit conserved-stress calculation in
[Zahn, section4.2](https://arxiv.org/pdf/1210.4031) assumes vanishing gauge
curvature; it is not an implemented nonzero-$F$ correction for our model.
The joint charged correction must be included in the actual point-split
calculation.

After that correction, set the remaining finite action additions
$\int\sqrt{-g}(c_1R^2+c_2R_{ab}R^{ab}+c_FF_{ab}F^{ab})$ to
$c_1=c_2=c_F=0$ at this scale, using the **charge-included** $F=dA$.
Keep the renormalized Einstein--Maxwell couplings fixed consistently.
This legitimate normalization choice does not remove the anomaly or the
fixed Ward terms. Scheme or scale changes require compensating changes
of the gravitational and Maxwell couplings; none may be fitted to obtain
support. These conditions define the baseline convention, while the
charged coincidence-limit correction and finite kernel remain to be
computed. In this prescription $\mathcal R$, define

$$
\mathcal Q_\phi^{\mathcal R}[f]
=\int_\Sigma fN\sqrt h\,
 (\langle T_{\hat0\hat0}\rangle_\phi^{\mathcal R}
 +\langle T_{\hat s\hat s}\rangle_\phi^{\mathcal R})\,d^3x.
$$

The stress is obtained from the differentiated full two-point kernel minus
the charged Hadamard parametrix, with the local terms required by
$\mathcal R$ and its Ward identities. This compact integral is the primary
observable; a finite total vacuum energy on the noncompact space is not
assumed. When the ground-state energy variation exists, it equals
$\partial_\epsilon E_{\rm ren}$ for
$N_\epsilon=e^{\epsilon f}N$, $h_{ss,\epsilon}=e^{-2\epsilon f}h_{ss}$,
with angular metric and $A$ fixed. Differentiability through zero-energy
crossings and interchange with infrared limits require separate control.

The actual130 mask gives the required non-Maxwell source
$-.171015352199415D/(8\pi G)$ in the Einstein--Maxwell convention.
Maxwell is already subtracted; finite higher-curvature terms and any other
current-carrying matter must be included consistently before identifying
this number with the Dirac contribution alone.

The132 determinant computes only a finite approximation to an interaction
part, $-(2/\pi)\int\Re\log\det(I-BR_{\rm ext})d\kappa$.
Its open one-pass factors are not a specified self-adjoint reference state.
A justified decomposition would still require that reference's metric
derivative, local subtraction terms, both chiralities, every $\mu$, both
frequency tails and the physical exterior/angular limits. Flat-holonomy
contrasts cancel common local terms but not the need for an absolute
source; the required geometric target is identical in both phases.

The absorbing DG boundary is not the domain proved here. The concrete next
source calculation must recover the compactly smeared differentiated
kernel of this sewn operator, subtract locally in a declared scheme, and
control its omitted modes and frequencies. State existence is established;
its required negative source, Maxwell current, backreaction and stability
remain unproved.
