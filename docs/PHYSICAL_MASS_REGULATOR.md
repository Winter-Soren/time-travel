# Physical Dirac masses and Pauli--Villars subtraction

Increment135. This analytic note specifies a regulator candidate on the
[complete smooth sewn background](SEWN_SOURCE_STATE.md). It does not report
numerical135 results, a completed Ward-identity calculation, or supporting
stress. The local subtraction below is derived in a declared heat/zeta
convention; its finite matching to133's Hadamard convention remains missing.

## Physical operator, domain and exterior ports

For the positive-Weyl Hamiltonian $H_0$, a constant **proper scalar mass** $m$
gives the full four-component operator

$$
H_m=\begin{pmatrix}H_0&mN\\mN&-H_0\end{pmatrix},\qquad
H_m^2=\begin{pmatrix}
H_0^2+m^2N^2&m[H_0,N]\\
-m[H_0,N]&H_0^2+m^2N^2
\end{pmatrix}.
$$

Both Weyls have the same charge, magnetic connection, spin structure and
loop phase. Mass couples them; it does not conjugate the charge. Replacing
the massless optical square by $H_0^2+m^2$ omits both lapse factors and the
commutator.

All masses use the same canonical Hilbert space and compact smooth core.
Complete-optical cutoffs prove essential self-adjointness separately for
each $m$: the smooth Hermitian mass potential cancels in the deficiency
Green identity, leaving the same vanishing cutoff commutator as in133.
The maximal domains are $\{u\in L^2:H_mu\in L^2\}$ distributionally.
The unbounded asymptotic lapse prevents assuming these graph domains are
identical for different masses. No artificial boundary condition may be
changed between regulators.

Multiplying **both** equations of $(H_m-i\kappa)u=0$ by $i$ gives diagonal
terms $iH_0+\kappa,-iH_0+\kappa$ and two equal off-diagonal terms $+imN$.
Dividing by $v=N/P$ gives frequency weight $P/N$ and mass coefficient $imP$.
Thus, in the exterior DG weak measure $\rho\,d\rho dz$, define

$$
W_{ab}=\int\rho(P/N)\,\overline\phi_a\phi_b\,d\rho dz,
\qquad P_{ab}=\int\rho P\,\overline\phi_a\phi_b\,d\rho dz.
$$

If $A_0$ contains the conservative positive-Weyl operator plus self-adjoint
upwind loss, the full matrix is

$$
\mathcal A_m=\begin{pmatrix}
A_0+\kappa W&imP\\imP&A_0^\dagger+\kappa W
\end{pmatrix}.
$$

The adjoint reverses the conservative terms and retains the loss. Since
$P=P^\dagger$, the mass blocks contribute nothing to the Hermitian part
or absorption ledger. The same regular-axis active degrees of freedom apply.

For $R=C(A_0+\kappa W)^{-1}G$, use input lift
$\operatorname{diag}(G,C^\dagger)$ and receiver
$\operatorname{diag}(C,G^\dagger)$. At $m=0$ their return is exactly
$\operatorname{diag}(R,R^\dagger)$ for real $\kappa$.
Lower incoming channels are the original outgoing channels, and vice versa.
This dual ordering is essential; it is not charge conjugation. Any separate
lower spinor rephasing must also transform the mass blocks.

The [massive bridge helper](../scripts/dihole_massive_bridge.py) uses the
same physical mass through $mP$ in flux coordinates. Under the compact
paired variation $N_\epsilon=e^{\epsilon f}N$,
$h_{ss,\epsilon}=e^{-2\epsilon f}h_{ss}$, the proper mass stays fixed:
$W_\epsilon=e^{-2\epsilon f}W$ and the entire transverse coupling,
including $mP$, scales by $e^{-\epsilon f}$.

## Regulator moments and the flat calibration

Take $m_j^2=jM^2$, $j=0,1,2,3$, with $c=(1,-3,3,-1)$. Then

$$
\sum c_j=\sum c_jm_j^2=\sum c_jm_j^4=0,\qquad
\sum c_j e^{-sm_j^2}=(1-e^{-sM^2})^3=O(s^3).
$$

For the positive **physical four-dimensional Euclidean** charged Dirac
square $\mathscr P$, the parity-even action is

$$
\Gamma_{\rm PV}=-\frac12\sum_jc_j\log\det(\mathscr P+m_j^2)
=\frac12\int_0^\infty\frac{ds}{s}
(1-e^{-sM^2})^3\operatorname{Tr}e^{-s\mathscr P}.
$$

The $s^3$ factor regularizes the four-dimensional ultraviolet integral at
fixed $M$. It does not remove infrared or infinite-volume problems: use
relative actions or compact variations. This factorization does not apply
to $H_0^2+m^2$ on the optical slice. Covariant Dirac-mass PV preserves the
vector gauge symmetry; regulator masses still generate the trace anomaly.
See [Bastianelli--Martelli, section2.2](https://arxiv.org/pdf/1610.02304)
and [Bastianelli--Broccoli, section4.2.1](https://arxiv.org/pdf/1808.03489).
The negative coefficients are subtraction weights, not exotic matter or
a physical state.

For134's constant flat segment, let $w=P/N$ and
$t_j^2=t_0^2+P^2jM^2$. Each power-subtracted bulk derivative has integral

$$
\int_0^\infty fL\left[2\kappa w-
\frac{2\kappa^2w^2+t_j^2}{\sqrt{\kappa^2w^2+t_j^2}}\right]d\kappa
=-\frac{fLt_j^2}{2w}.
$$

Its signed PV sum vanishes exactly because $\sum c_jt_j^2=0$.
This removes the flat bulk regulator artifact at every finite $M$.
Finite-boundary/Casimir contributions need not vanish. Passing this
calibration alone does not establish curved-space covariance.

## Full-Dirac determinant normalization

The symmetric vacuum convention is $E=-\operatorname{Tr}|H_m|/2$.
For an already full-Dirac, canonically normalized determinant, it gives

$$
E'_{\rm PV}=-\frac1\pi\int_0^\infty
\sum_jc_j\,\partial_\epsilon\Re\log
D_{{\rm full},m_j}(i\kappa)\,d\kappa.
$$

Include the direct/reference normalization and the sewn interaction
determinant. At $m=0$, adjoint sewing in the common physical Hilbert space
makes the real full-Dirac log twice the one-Weyl log, recovering132's
$-2/\pi$. This is not an exact identity of the finite polygon pilot:
its cut transforms $F_i,F_o$ need not be isometries. Although the raw
bridge blocks are $B,B^\dagger$, their projected blocks are
$B_R=F_iBF_o^{-1}$ and $B_L=F_oB^\dagger F_i^{-1}$, generally
$B_L\ne B_R^\dagger$. Test the raw adjoint identity and independently
projected two-Weyl sewing; do not attribute this projection defect to ODE
tolerance. The physical full-Dirac prefactor remains $-1/\pi$, with no
further chirality, antiparticle or azimuthal-sector doubling. A raw finite
generalized pencil still needs its canonical measure normalization; PV
weights do not authorize omitting terms inconsistently between masses.

## Local large-mass subtraction and remaining limits

Use sphere-positive curvature, $\Delta=\nabla^a\nabla_a$, and
charge-included $F=dA$, all for the physical Euclidean metric. Define
$\operatorname{tr}e^{-s\mathscr P}(x,x)\sim
(4\pi s)^{-2}(a_0+sa_2+s^2a_4+\cdots)$. The rank-four coefficients are

$$
a_0=4,\qquad a_2=-R/3,\qquad
a_4=\frac{5R^2-8R_{ab}R^{ab}-7R_{abcd}R^{abcd}-12\Delta R}{360}
+\frac23F_{ab}F^{ab}.
$$

These follow from the Dirac connection and general Laplace-type formulas
in [Vassilevich, equations3.27--3.28 and4.26--4.28](https://arxiv.org/pdf/hep-th/0306138).
The three-dimensional optical curvature identity used in121 is inapplicable.

In a **declared zeta/heat finite-part convention** at scale $\mu_R$, the
local terms to subtract as $M\to\infty$ are

$$
\Gamma_{\rm loc}(M)=\frac1{32\pi^2}\int\sqrt g\,
\left[\beta_0M^4a_0+\beta_1M^2a_2+
\left(\log\frac{M^2}{\mu_R^2}+\log\frac38\right)a_4\right],
$$
$$
\beta_0=-\tfrac12\sum_{j=1}^3c_jj^2\log j
=\tfrac12(9\log3-12\log2),\qquad
\beta_1=\sum_{j=1}^3c_jj\log j=6\log2-3\log3.
$$

Here $\mu_R$ is a renormalization scale, not the azimuthal label.
The remaining mass-ratio sum is $\sum_{j=1}^3c_j\log j=\log(8/3)$.
The displayed formula specifies the local asymptotic subtraction, not
existence of a global determinant or a uniform remainder on this noncompact
geometry.

Counterterms therefore include cosmological, Einstein, curvature-squared
and Maxwell terms. Keep Euler and total-derivative terms in local
bookkeeping; their compact integrated bulk variations vanish without
boundaries. PV moment cancellation does not remove the anomaly.

Subtracting $\Gamma_{\rm loc}$ is **not proven to match133's prescription**.
Its finite local conversion must be fixed at the declared Hadamard length
$\ell_*=D_0$, including joint charged Ward normalization and the specified
zero finite $R^2$, $R_{ab}R^{ab}$ and $F^2$ additions. Einstein--Maxwell
couplings remain fixed consistently. No coefficient may be fitted for
support. The relevant finite-renormalization freedom is discussed in
[Zahn, sections4.2--4.3](https://arxiv.org/pdf/1210.4031).

At each fixed $M$, combine the regulator contributions consistently and
control frequency tails, spatial resolution, angular towers, every
azimuthal sector and the physical exterior limit before increasing $M$.
A fixed angular truncation or absorbing DG exterior is not a covariant
limit. The needed controls are massless full-Dirac recovery, the exact flat
PV cancellation, physical mass/commutator agreement, common sewing and
phase, and ultimately the physical heat coefficients and joint Ward
identities. These checks do not yet supply the finite renormalized null
source required by the bridge geometry.
