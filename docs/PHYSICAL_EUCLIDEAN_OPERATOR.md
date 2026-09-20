# A common physical Euclidean Dirac operator

Increment137. This constructs the common positive operator needed by the
fixed-mass PV identity in [136](DIRAC_SCHEME_MATCHING.md). The background is
the exact smooth sewn manifold, magnetic bundle and chosen spin/flat-phase
sewing of [133](SEWN_SOURCE_STATE.md). The ultraviolet reference is now136's
explicit physical heat/zeta action; equality to133's point-split prescription
is not assumed. No absorbing numerical boundary is introduced below.

## Physical geometry, Hilbert spaces and domains

Write the static metric as $g=-N^2dt^2+h$, with $N>0$, and its Euclidean
section as $g_E=N^2d\tau^2+h$ on $\mathbb R_\tau\times\Sigma$. The connection
is spatial and unitary, with the same integer charge, polar patches and
handle phase as133. The four-component field is vectorlike.

The **physical** spatial metric $h$ is complete as well as the optical metric.
Indeed, the exact positive balanced exterior has
$N=\Lambda\ge1-2Mr/\Sigma_g$, where
$\Sigma_g=r^2-a_g^2\cos^2\theta$. This lower bound tends uniformly to1
at infinity; the lapse itself need not do so. The remaining terms in the numerator of
$\Lambda$ are nonnegative; see [Emparan, equations8--10](https://arxiv.org/pdf/hep-th/9906160).
The bridge and remaining bounded exterior are smooth with positive lapse,
so $N\ge c_N>0$ globally. Hence $h=N^2\widehat h\ge c_N^2\widehat h$,
and133's optical completeness implies physical completeness. A finite-length
Euclidean curve has a finite-length spatial projection, which converges in
$h$; near its limit $N$ has a positive lower bound, so its time coordinate
also cannot escape. Thus $g_E$ is complete. There is no horizon or Euclidean
bolt in the sewn geometry. The compact paired variations preserve this result.

Use Hamiltonian and Euclidean spatial norms

$$
\mathcal H=L^2(\Sigma,d\mathrm{vol}_h;S_4\otimes L_q),\qquad
\mathcal H_E=L^2(\Sigma,N\,d\mathrm{vol}_h;S_4\otimes L_q).
$$

The full Euclidean norm additionally integrates $d\tau$. Fix Hermitian
Hamiltonian matrices $\alpha^a,\beta$, with $\beta^2=I$ and
$\{\beta,\alpha^a\}=0$. In the project's chiral representation,
$H_0=\operatorname{diag}(h_W,-h_W)$ and
$\beta=\left(\begin{smallmatrix}0&I\\I&0\end{smallmatrix}\right)$.
Here $h_W$ includes the physical lapse connection:

$$
H_0=-iN\alpha^i\left(\nabla_i^{h,A}
+\tfrac12\partial_i\log N\right),\qquad \{\beta,H_0\}=0.
$$

Choose Euclidean Clifford matrices $\gamma^0=\beta$,
$\gamma^i=-i\beta\alpha^i$. The skew-Hermitian physical Dirac expression is

$$
\boxed{\mathcal D_E
=\beta N^{-1}\partial_\tau
+\gamma^i\left(\nabla_i^{h,A}+\tfrac12\partial_i\log N\right)
=N^{-1}\beta(\partial_\tau+H_0).}
$$

The convention differs by a factor of $i$ from a self-adjoint Euclidean
Dirac operator. The local Clifford-connection square is the usual
Laplace-type operator; compare [Vassilevich, equations3.26--3.28](https://arxiv.org/pdf/hep-th/0306138).
On the complete smooth twisted Dirac bundle, $-i\mathcal D_E$ and its square
are essentially self-adjoint on compact smooth sections by
[Bandara--Saratchandran, section3.1 and theorem3.5](https://arxiv.org/pdf/1606.06190).
Their hypotheses concern the complete metric and compatible smooth unitary
connection, not a bounded global gauge potential. Define

$$
\mathscr P=\overline{\mathcal D_E}^{\,\dagger}
\overline{\mathcal D_E}\ge0.
$$

All proper scalar masses use this same realization:
$(\mathcal D_E+m)^\dagger(\mathcal D_E+m)=\mathscr P+m^2$.
The first-order Euclidean domains are identical for finite $m$, since $mI$
is bounded. This does not assert identical maximal domains for the physical
Hamiltonians $H_m=H_0+mN\beta$, whose multiplication coefficient $N$ is unbounded.

## Fourier fibre and a form that can be assembled

Use $e^{i\kappa\tau}$, so $\partial_\tau=i\kappa$ and
$A_\kappa=i\kappa+H_0$. Existing $H_m-i\kappa$ conventions correspond to
$\kappa\mapsto-\kappa$. The Euclidean fibre and its positive form are

$$
\mathcal D_\kappa=N^{-1}\beta A_\kappa,\qquad
q_\kappa[u,v]=\langle A_\kappa u,N^{-1}A_\kappa v\rangle_{\mathcal H},
\qquad n[u,v]=\langle u,Nv\rangle_{\mathcal H}.
$$

Weighted adjunction gives $\mathcal D_\kappa^{\dagger_E}=-\mathcal D_\kappa$.
The unitary map $U_Nu=N^{1/2}u$ from $\mathcal H_E$ to $\mathcal H$ yields

$$
\widehat{\mathcal D}_\kappa
=D_h+i\kappa\beta/N,
\qquad D_h=\gamma^i\nabla_i^{h,A}.
$$

The spatial Dirac operator $-iD_h$ is self-adjoint on complete $h$, and
$\kappa\beta/N$ is bounded. Thus every real-$\kappa$ fibre has a definite
closure with a common spatial first-order domain; $q_\kappa$ is its closed
norm-square form. Axial projections reduce it with the global labels and
gauge transitions of133. The associated $P_\kappa$ is the Fourier fibre
of $\mathscr P$. This conclusion does not give a positive spectral gap:
even $\kappa\ne0$ supplies no uniform lower bound through $\kappa^2/N^2$.

For real $m$,

$$
\mathcal D_\kappa+m=N^{-1}\beta(i\kappa+H_m),\qquad
\boxed{q_{\kappa,m}=q_\kappa+m^2 n.}
$$

The cross term vanishes because
$(H_0-i\kappa)\beta+\beta(H_0+i\kappa)=0$ on the core. Closure gives the
same form identity. For conforming trial sections $\phi_a$, assemble

$$
G_{ab}=n[\phi_a,\phi_b],\qquad Q_{ab}=q_\kappa[\phi_a,\phi_b],\qquad
P_G=G^{-1/2}QG^{-1/2}.
$$

Every regulator is the **same** matrix pencil $Q+jM^2G$. Integrate the
actual differential actions in $Q$, including their components outside
the trial space. Squaring a projected first-order matrix generally omits
those components and gives a different form. This distinction is testable
with a basis that is not invariant under the differential operator.

## Bridge flux variables and the exact paired variation

For one Weyl field, the133 bridge flux amplitude is
$v_b=R\sqrt{P\sin\theta}\,\Psi_\mu$, with measure $ds\,d\theta$ before
weights. Set $W=P/N$ and

$$
T_W=-i\sigma_3\partial_s+\mathcal K/R,\qquad
\mathcal K=-i\sigma_1\partial_\theta
+\sigma_2\frac{NP(\mu-qa)}{\sin\theta}.
$$

Use the regular polar spin/gauge domain, not arbitrary endpoint conditions
on this coordinate expression. For the full field let
$T_4=\operatorname{diag}(T_W,-T_W)=F+\mathcal A$,
$F=-i\Sigma\partial_s$, $\Sigma=\operatorname{diag}(\sigma_3,-\sigma_3)$,
and $\mathcal A=\operatorname{diag}(\mathcal K/R,-\mathcal K/R)$.
Then

$$
H_0=W^{-1}T_4,\quad
\mathcal D_\kappa=P^{-1}\beta(i\kappa W+T_4),\quad
n_b[v]=\int P|v|^2,\quad
q_{\kappa,b}[v]=\int\left|P^{-1/2}(i\kappa W+T_4)v\right|^2.
$$

For $N_\epsilon=e^{\epsilon f(s)}N$ and
$h_{ss,\epsilon}=e^{-2\epsilon f(s)}h_{ss}$, keep the other spatial
coefficients and $A$ fixed. The physical volume $N_\epsilon\sqrt{h_\epsilon}$
is unchanged. In the varying flux coordinates, however,

$$
W_\epsilon=e^{-2\epsilon f}W,\qquad
T_{4,\epsilon}=F+e^{-\epsilon f}\mathcal A,\qquad
P_\epsilon=N_\epsilon W_\epsilon=e^{-\epsilon f}P.
$$

Here $P_\epsilon$ denotes the Euclidean flux density; it is not a replacement
of $P$ in all angular metric components. The unitary fixed-space identification
is $v_\epsilon=J_\epsilon v$, $J_\epsilon=e^{\epsilon f/2}$, from baseline
$L^2(P\,dsd\theta)$ to $L^2(P_\epsilon\,dsd\theta)$. Therefore

$$
q_{\kappa,\epsilon}[J_\epsilon v]=\|L_\epsilon v\|_{L^2(dsd\theta)}^2,
$$
$$
\boxed{L_\epsilon=P^{-1/2}\left[
e^{-\epsilon f}i\kappa W
+e^{\epsilon f}\left(F-\tfrac{i\epsilon}{2}f'\Sigma\right)
+\mathcal A\right],\qquad
\dot L=P^{-1/2}\left[-i\kappa fW+fF-\tfrac i2f'\Sigma\right].}
$$

The $f'$ term comes from differentiating $J_\epsilon v$ and cannot be
dropped. The angular operator cancels from $\dot L$, but remains in $L$:
$\dot q[u,v]=\langle\dot Lu,Lv\rangle+\langle Lu,\dot Lv\rangle$.
The norm matrix $G$ and proper mass term $m^2G$ are constant under this
identification. Compact smooth variations give equivalent form graph norms
by local elliptic estimates and unchanged coefficients outside the support.

For the exterior amplitude $v_e=\sqrt P\,\Psi_\mu$, the corresponding
Hamiltonian and Euclidean weights are $\rho W$ and $\rho P$. Use

$$
T_{e,W}=-i\sigma_x(\partial_\rho+1/(2\rho))
-i\sigma_z\partial_z+\sigma_y NP(\mu-qa)/\rho
$$

and the same formulas with measure $\rho\,d\rho dz$. At each mouth the
physical trace relation is
$v_b=R_c\sqrt{\sin\theta}\,U_{R,L}^\dagger v_e$, with133's fixed spin lift
and flat phase. These are internal sewing conditions, not absorbing ends.

## Relation to first-order determinants and what remains

On the smooth core there is the exact congruence

$$
\widehat{\mathcal D}_\kappa+m
=\beta N^{-1/2}(H_m+i\kappa)N^{-1/2}.
$$

It is not similarity to $H_m+i\kappa$. In a common finite representation
that respects this factorization exactly,

$$
\log\det(P_\kappa+m^2)
=2\Re\log\det(H_m+i\kappa)-2\log\det N.
$$

The mass-independent lapse factor cancels in the common PV sum because
$\sum c_j=0$. This algebra does not authorize separate multiplicative
factorizations of infinite-dimensional zeta determinants. The continuum
comparison requires the common domains, Hilbert identification, and a
jointly regulated relative trace permitting the cyclic operations.
In particular $N^{1/2}$ is unbounded at infinity. Boundary terms and the
local normalization must remain those of the chosen physical action.

Under those conditions the parity-even full-Dirac normalization is
$-\tfrac12\int d\kappa/(2\pi)\,\log\det(P_\kappa+m^2)$, equivalently
$-1/\pi\int_0^\infty\Re\log\det(H_m-i\kappa)\,d\kappa$ after the common
normalization. There is no further Weyl doubling. The136 local subtraction
is still required when the proper regulator mass tends to infinity.

A form-supported calculation on a finite bridge with zero traces converges,
at fixed cuts, to an artificial Dirichlet form problem. It does not replace
the natural full sewn exterior. Genuine conforming global trials must match
the bundle/spin traces and exhaust the physical domain, or come with a
controlled exterior resolvent/form reduction. Existing absorbing DG matrices
have neither been squared into this realization nor proved equivalent to it.

The immediate implementable object is $Q$, $G$ and $\dot Q$ from the displayed
forms, with all masses sharing $Q+jM^2G$. This settles operator and measure
normalization. Spatial/exterior, angular and frequency remainders remain to
be bounded. The infinite-volume relative trace and infrared limit, including
derivatives at zero crossings, are separate: retain133's empty exact-zero
prescription and do not silently replace it by a finite-temperature limit.
No renormalized supporting stress or physical throat solution follows from
this construction alone.

## Remainders tied to the physical source test

These are our conditional estimates for the common operator above, not
evaluated error bars. They identify what the next global calculation must
measure. All traces are per unit Euclidean Killing time and the compact
metric variation is represented on the fixed Euclidean Hilbert space.
Write its symmetric form as $\delta P=L_f^\dagger J_fL_f$, where $L_f$ is a
first-order row operator with compact spatial support and
$\|J_f\|\le C_f$. Define a positive and a signed local spectral measure:

$$
\mu_f(B)=\operatorname{Tr}_{/\mathcal T}
 [L_fE_P(B)L_f^\dagger],\qquad
\nu_f(B)=\operatorname{Tr}_{/\mathcal T}
 [J_fL_fE_P(B)L_f^\dagger],\qquad |\nu_f|\le C_f\mu_f.
$$

Whenever the indicated relative trace exists, the exact136 resolvent
identity gives

$$
Q_M=-\int_0^\infty w_M(\lambda)\,d\nu_f(\lambda),\qquad
w_M(\lambda)=\frac{3M^6}
 {\lambda(\lambda+M^2)(\lambda+2M^2)(\lambda+3M^2)}.
$$

For a **uniform proved bound** $\mu_f([0,\lambda])\le A\lambda^3$ for
every $\lambda\ge\lambda_*$, integration by parts, dropping the nonpositive
lower-end term, and $-w_M'(\lambda)\le12M^6\lambda^{-5}$ imply

$$
|Q_M^{\lambda\ge\lambda_*}|\le\frac{12C_fAM^6}{\lambda_*}.
$$

The cutoff $\lambda_*$ has units of mass squared. A local Weyl asymptotic
alone supplies neither this uniform inequality nor its constant $A$.
The bound is on the exact spectral projection of $P$, not on eigenvalues
of a separately compressed determinant.

In the infrared $w_M(\lambda)\sim1/(2\lambda)$. Absolute convergence of
the signed scalar response at zero requires and, away from zero atoms,
is equivalent to

$$
\int_{(0,\eta)}\lambda^{-1}\,d|\nu_f|(\lambda)<\infty.
$$

A stronger sufficient condition uses $C_f\mu_f$. For example, if
$\mu_f([0,\lambda])\le B\lambda^{1+\epsilon}$ at small $\lambda$, with
$\epsilon>0$, then $w_M\le1/(2\lambda)$ and integration by parts give

$$
|Q_M^{0<\lambda<\eta}|
\le\frac{C_fB(1+\epsilon)}{2\epsilon}\eta^\epsilon.
$$

No such threshold estimate has been proved here. Exact zero atoms and
differentiation through zero crossings retain the separate state convention;
a positive finite-trial eigenvalue does not settle them.

The paired variation offers a more efficient possible UV route. Since
$\delta a_0=0$, the signed inserted heat trace is formally

$$
h_f(s)=\operatorname{Tr}_{/\mathcal T}(\delta P e^{-sP})
=-\frac1s\delta\operatorname{Tr}_{/\mathcal T}e^{-sP}
\sim-\frac{a_2'}{16\pi^2s^2}-\frac{a_4'}{16\pi^2s}+O(1),
$$

where136 already supplies $a_2'=0.11400999281577749D$ and
$a_4'=7.706927314195081/D$. If a remainder estimate establishes
$|h_f(s)|\le H_2s^{-2}$ throughout $0<s\le s_0$, then

$$
\left| -\tfrac12\int_0^{s_0}(1-e^{-sM^2})^3h_f(s)\,ds\right|
\le\frac{H_2M^6s_0^2}{4}.
$$

This follows directly from $1-e^{-sM^2}\le sM^2$. The asymptotic coefficients
do not prove $H_2$, and this signed cancellation is not an estimate of the
total-variation spectral measure. The heat/resolvent framework is given in
[Vassilevich, sections1--2](https://arxiv.org/pdf/hep-th/0306138); the displayed
conditional bounds are derived here for this response.

For one Dirac field define the final, dimensionless physical response by
$Q_{\rm HK}=(\hbar c/D)C$. At the prescribed dimensionless geometry, the
existing masked non-Maxwell source target is

$$
Q_{\rm required}=-\frac{0.171015352199415\,c^4D}{8\pi G},\qquad
\boxed{C_{\rm required}=-0.006804484661784583\,(D/\ell_P)^2},\quad
\ell_P^2=\hbar G/c^3.
$$

Only after global fixed-$M$ control,136's local subtraction and regulator
removal can an interval for $C$ be compared to that target. This is a
necessary smeared Einstein-equation test, not a sufficient throat solution:
other stress components, induced current, backreaction and stability remain.
No finite compact-trial response in137, or absorbing total in136, is an
estimate of $C$ with a known remainder. The next calculation therefore needs
globally sewn trial traces or a controlled physical exterior reduction,
not more refinement of the isolated bridge Dirichlet problem.

The focused [implementation](../scripts/dihole_euclidean_form.py) exposes
`trial_data(bridge, radial_modes=3, radial_order=64)` and
`assemble(data, kappa, epsilon=0., mass=0.)`, returning `q`, `norm`,
`derivative` and `direct_massive` among its controls. Its
[saved results](../results/dihole_euclidean_form_results.json) use actual
bridge coefficients and compact $H^1$ sine trials in three azimuthal sectors.
The assembly evaluates multiplication by $P^{-1}$ pointwise, retains the
$f'$ term, and compares the direct massive form with $Q+m^2G$. These are
finite form/variation checks. The trials vanish at the mouth cuts, and
neither this calculation nor its refinement contains the sewn exterior.
