# Natural Dirac domain on the uncut balanced dihole exterior

Increment118. This note supplies the operator-domain identification needed by
increment117. It concerns a massless charged test field on the **uncut classical
black-hole exterior**, with the regular symmetry axis included. It establishes a
unique self-adjoint Hamiltonian and its fixed-azimuthal restrictions. It does not
establish an operator on a replacement throat, a quantum state, or a stress tensor.

## 1. Geometry and the complete optical manifold

Use the positive balanced branch of Emparan's *Black diholes*, equations8--16,
with azimuthal period $2\pi$ ([primary paper](https://arxiv.org/pdf/hep-th/9906160)).
These are the fields in [balanced_dihole.py](../scripts/balanced_dihole.py).
Write $a_g$ for the geometric separation parameter, reserving $a(\rho,z)$ for the
unit-flux gauge coefficient. The parameters obey

$$
0<M<k,\quad a_g^2=k^2-M^2,\quad r_+=M+k,\quad
Q=Mr_+/k,\quad B=(k-a_g)/(Mr_+)>0.
$$

The exact target calibration used here is

$$
D=1,\quad Q=\tfrac15,\quad k=D/2,\quad
M=\frac{\sqrt{k^2+4Qk}-k}{2},\quad
a_g=\sqrt{k^2-M^2},\quad B=\frac{M}{(k+a_g)(M+k)}.
$$

Stored floating-point parameters approximate these linked exact constants.
Treating every stored value as an independently exact parameter could introduce
a small conical imbalance, outside the smooth-axis proof below. The domain
conclusion concerns the exact balanced target. Any certified curved-response
calculation must enclose evaluation error in these constants as well as in the
coefficient functions and integrals.

The spatial manifold is $X=\mathbb R^3\setminus\{p_+,p_-\}$, where
$p_\pm=(0,0,\pm k)$. The removed points represent horizon ends, not regular
points to be filled in. In cylindrical coordinates away from the axis,

$$
g=-N^2dt^2+P^2(d\rho^2+dz^2)+\frac{\rho^2}{N^2}d\varphi^2.
$$

Set

$$
d_\pm=\sqrt{\rho^2+(z\mp k)^2},\quad
r=M+\frac{d_++d_-}{2},\quad c=\frac{d_--d_+}{2k},
$$
$$
\Sigma=r^2-a_g^2c^2,\qquad
\mathfrak h=(r-M)^2-k^2c^2=d_+d_-,\qquad
\frac PN=\frac{\Sigma^2}{\mathfrak h^2}.
$$

All coefficients are smooth and positive away from the horizon ends after the
regular axis is included. Indeed, in the published expression for $N$, the
numerator is a sum of nonnegative terms and vanishes only at those ends. Also

$$
\Sigma-\mathfrak h=2Mr-M^2(1-c^2)>0.
$$

The optical spatial metric is

$$
\widehat h=N^{-2}h
=\alpha^2(d\rho^2+dz^2)+\beta^2d\varphi^2,
\qquad \alpha=P/N,\quad \beta=\rho/N^2.
$$

Its completeness follows from three explicit checks.

* **Spatial infinity:** $\alpha\ge1$. A path escaping to infinity in $X$ has
  unbounded $\sqrt{\rho^2+z^2}$ and therefore infinite optical length, regardless
  of the shrinking azimuthal circles in the Melvin region.
* **Either horizon:** if $\epsilon=d_+\to0$, then $d_-\to2k$,
  $\Sigma\to r_+^2-a_g^2=2Qk$, and
  $\alpha\epsilon^2\to Q^2$, uniformly in the direction of approach. Thus, in
  a sufficiently small punctured neighborhood,
  $ds_{\widehat h}\ge C\epsilon^{-2}|d\epsilon|$ with $C>0$; the end has infinite
  optical distance. The lower end has the same estimate.
* **Axis:** the axis is part of $X$, not a boundary. On the outer segments,
  $NP=1$. On the middle segment, balance gives $a_g+BMr_+=k$ and hence
  $N=\mathfrak h/\Sigma$, $P=\Sigma/\mathfrak h$, again $NP=1$.
  Away from $p_\pm$, the coefficient functions are smooth in $\rho^2$.
  Consequently $\beta/\rho=\alpha$ at the axis and their difference is
  $O(\rho^2)$: the optical metric extends smoothly in Cartesian coordinates,
  without a conical defect.

Every escape from compact subsets either approaches a horizon end or spatial
infinity; the remaining bounded regions, including the axis, have a smooth
positive metric. These checks prove completeness of $(X,\widehat h)$.

## 2. Spin structure, magnetic bundle, and conformal Hilbert map

$X$ inherits a spin structure from $\mathbb R^3$. It is unique because
$H^1(X;\mathbb Z_2)=0$. Use its rank-two spatial spin bundle $S$ for one Weyl
chirality. The other chirality has the opposite Hamiltonian sign below; their
direct sum is the massless four-component Dirac field.

The project uses the charge-included connection
$A^o=q\,a(\rho,z)d\varphi$, where
$a=-\widehat A_\varphi/(2Q)$ and $\widehat A_\varphi$ is Emparan's potential.
Its axis limits are $a=0$ on the outer segments and $a=1$ on the middle segment.
For an **integer** $q$, introduce the middle-axis patch

$$
A^c=A^o-q\,d\varphi,\qquad
\Psi^c=e^{-iq\varphi}\Psi^o,
\qquad \nabla^A=\nabla^S-iA.
$$

$A^o=O(\rho^2)d\varphi$ near the outer axis and
$A^c=O(\rho^2)d\varphi$ near the middle axis, so both are smooth connections in
their respective Cartesian neighborhoods. Their overlap is off the axis, where
$e^{-iq\varphi}$ is single-valued precisely for integer $q$. This constructs a
Hermitian line bundle $L_q$ with unitary connection and Chern numbers $+q,-q$
on outward spheres around the upper and lower ends. Use the twisted Dirac bundle
$S\otimes L_q$. A globally single-valued potential is neither required nor used.

Write $g=N^2\widehat g$, with $\widehat g=-dt^2+\widehat h$. Under the standard
identification of conformally related spin frames, the four-dimensional massless
Dirac identity is

$$
\slashed D_g^{,A}\Psi
=N^{-5/2}\slashed D_{\widehat g}^{,A}(N^{3/2}\Psi).
$$

The spinor weight follows from the spin-connection transformation; see
Godazgar--Nicolai, equations9 and13--15
([primary paper](https://arxiv.org/pdf/1801.01728)). The magnetic one-form has
conformal weight zero, so the same identity holds after minimal coupling.
There is no electric potential here. Thus
$\widehat\Psi=N^{3/2}\Psi$ obeys an ultrastatic equation
$i\partial_t\widehat\Psi=D_A\widehat\Psi$ for the chosen chirality.

This transformation is unitary for the conserved one-particle norms:

$$
d\mathrm{vol}_{h}=\frac{P^2\rho}{N}\,d\rho\,dz\,d\varphi,
\qquad
d\mathrm{vol}_{\widehat h}=\frac{P^2\rho}{N^4}\,d\rho\,dz\,d\varphi,
\qquad
\|N^{3/2}\Psi\|_{L^2(\widehat h)}=\|\Psi\|_{L^2(h)}.
$$

No global upper or lower bound on $N$ is needed for this equality of weighted
spaces. Since $N$ is positive and smooth on $X$, the map also sends compact
smooth cores bijectively to one another.

## 3. Self-adjointness and the natural domain

Apply Bandara--Saratchandran, section3.1 and Theorem3.5
([primary paper](https://arxiv.org/pdf/1606.06190)). Its Dirac-bundle hypotheses
are met: the base has the complete smooth Riemannian metric just proved; the
bundle metric is Hermitian; Clifford multiplication is compatible with it; and
the spin connection tensored with the smooth unitary magnetic connection is a
compatible Clifford connection. The resulting Dirac operator is elliptic and
symmetric on $C_c^\infty(X,S\otimes L_q)$. The theorem therefore gives essential
self-adjointness on that core. It needs no global bound on the gauge potential.

Equivalently, the cutoff proof has no hidden horizon boundary condition. On a
complete optical manifold choose compact Lipschitz cutoffs tending to1 with
$|d\chi_R|_{\widehat h}\le C/R$. For a deficiency vector
$D_A^*\psi=\pm i\psi$, local ellipticity and compact approximation permit
$\chi_R\psi$ in the minimal graph domain. Since
$[D_A,\chi_R]=-i\,\sigma(d\chi_R)$,

$$
\|\chi_R\psi\|
\le\|(D_A\mp i)\chi_R\psi\|
\le(C/R)\|\psi\|.
$$

Both deficiency spaces vanish. The connection drops out of this commutator.
The unique closed domain is

$$
\mathcal D_A=\{\psi\in L^2(X,S\otimes L_q;\widehat h):
D_A\psi\in L^2\text{ in the bundle-distributional sense}\}.
$$

This is a graph-domain statement, not an assertion that derivatives and each
unbounded coefficient term separately belong to global $L^2$.

## 4. Exact reduced Hamiltonian and the global meaning of $\mu$

Choose the rotating cylindrical orthonormal frame, assigning
$(\sigma_x,\sigma_y,\sigma_z)$ to $(\rho,\varphi,z)$. Its spin lift changes sign
after a $2\pi$ rotation. Consequently the local spinor components are
antiperiodic in $\varphi$, although the geometric section is single-valued.
In this frame the optical operator is

$$
D_A=-\frac{i\sigma_x}{\alpha}
 \left(\partial_\rho+\tfrac12\partial_\rho\log(\alpha\beta)\right)
-\frac{i\sigma_z}{\alpha}
 \left(\partial_z+\tfrac12\partial_z\log(\alpha\beta)\right)
-\frac{i\sigma_y}{\beta}(\partial_\varphi-iA_\varphi).
$$

The half-logarithmic terms follow from the Levi-Civita spin connection of
$\alpha^2(d\rho^2+dz^2)+\beta^2d\varphi^2$. Conjugation by the volume half-density
$\alpha\sqrt\beta$ gives

$$
(\alpha\sqrt\beta)D_A(\alpha\sqrt\beta)^{-1}
=-i\sigma_x(v\partial_\rho+v_\rho/2)
 -i\sigma_z(v\partial_z+v_z/2)
 +b\sigma_y(-i\partial_\varphi-A_\varphi),
\qquad v=N/P,\quad b=N^2/\rho.
$$

For example, the radial connection changes by
$\frac12\partial_\rho\log(\alpha\beta)
-\partial_\rho\log(\alpha\sqrt\beta)
=-\frac12\partial_\rho\log\alpha$; this is exactly $v_\rho/2$ after multiplying
by $v$. There is no additional $1/(2\rho)$ term in the canonical operator.

Fix the global axial lift by taking, in the outer gauge and rotating frame,

$$
(U_\theta\psi)^o(\varphi)=\psi^o(\varphi-\theta),\qquad
(U_\theta\psi)^c(\varphi)=e^{-iq\theta}\psi^c(\varphi-\theta).
$$

These agree on the overlap with the transition above and extend smoothly across
the regular axis. They form a strongly continuous unitary action with
$U_{2\pi}=-I$, $U_{4\pi}=I$. Its generator is

$$
J^o=-i\partial_\varphi,\qquad J^c=-i\partial_\varphi+q.
$$

Thus one global eigenvalue $\mu\in\mathbb Z+\tfrac12$ corresponds to exponent
$e^{i\mu\varphi}$ in the outer patch and $e^{i(\mu-q)\varphi}$ in the central
patch. The gauge-invariant angular combination agrees:
$\mu-qa=(\mu-q)-q(a-1)$.

The lifted rotations preserve the connection and optical metric, preserve the
compact smooth core, and commute there with $D_A$. Essential self-adjointness
then makes them commute with its closure and resolvent. In particular,

$$
\Pi_\mu=\frac1{4\pi}\int_0^{4\pi}e^{i\mu\theta}U_\theta\,d\theta
$$

is a reducing projection. Its restriction has self-adjoint domain
$\Pi_\mu\mathcal D_A$ and core
$\Pi_\mu C_c^\infty(X,S\otimes L_q)$. Averaging over this compact rotation group
preserves compact support; applying the projection to graph approximants proves
the core assertion.

Use $e^{i\mu\varphi}/\sqrt{2\pi}$ for the outer-gauge Fourier factor. For its
physical amplitude $\Psi_\mu(\rho,z)$ the canonical two-dimensional spinor is

$$
\boxed{u=\alpha\sqrt\beta\,\widehat\Psi_\mu
       =P\sqrt{\rho/N}\,\Psi_\mu},\qquad
\|u\|^2=\int_0^\infty\!\int_{\mathbb R}|u|^2\,dz\,d\rho.
$$

The restriction is therefore unitarily represented by exactly

$$
\boxed{H_\mu=-i\sigma_x(v\partial_\rho+v_\rho/2)
 -i\sigma_z(v\partial_z+v_z/2)+\sigma_y b(\mu-qa)},
$$

the expression in [dihole_null_scaling.py](../scripts/dihole_null_scaling.py).
Its domain is the image of the reducing global graph domain above. It is not
defined by imposing an arbitrary boundary condition at $\rho=0$. As an axis
check on smooth core sections, in a regular gauge with local exponent $\nu$,
Cartesian spin components have integer azimuthal powers $\nu-1/2,\nu+1/2$;
their radial factors have the corresponding regular Taylor behavior. The
central value is $\nu=\mu-q$, not $\mu$. No pointwise Taylor condition is asserted
for every vector in the closed graph domain.

## 5. Consequence for increment117 and limitations

The chosen gauge, rotating spin frame, chirality, and axial lift now specify a
global realization matching increment117's local Pauli matrices and half-integer
labels. Compact $C^2$ trial fields supported away from the axis and horizon ends
reconstruct to compact sections in its graph domain by local smoothing. Hence
the exact imaginary-frequency residual and dual-response inequalities apply to
this **particular uncut classical exterior**, without an unspecified
self-adjoint-extension assumption. For the saved $q=512000$, $\mu=398908.5$
sector, the central-patch exponent is $-113091.5$, still half-integer.

This completes the mathematical domain identification; it does not certify the
quadrature of the curved residuals. Those numerical estimates become certified
response disks only after their integrals and arithmetic are enclosed. Physical
frequency is $iE_q\kappa$, and physical Green responses divide by $E_q$ (and by
the source norm squared if the source is normalized).

No unresolved spin/gauge alignment remains for the explicit conventions in this
note. A different gauge or axial lift requires transforming both coefficients
and local labels. A throat replacement or finite-cut boundary is a different
operator problem: the optical completeness and domain proved here do not extend
to it automatically. No real-axis limit, vacuum choice, renormalized stress,
backreaction, stable wormhole, or chronology change is established.
