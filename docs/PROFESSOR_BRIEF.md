# Technical assessment: Dirac routing near a magnetic null ring

**Proposed brief for first methodological feedback.** This note assesses
scattering near a magnetic null as possible matching data for fermions in a
specified curved exterior. It makes no claim of novelty, wormhole support or
chronology change.

The background is the externally balanced, magnetic Einstein–Maxwell dihole
of [Emparan, *Black diholes*, equations 8–16](https://arxiv.org/pdf/hep-th/9906160).
It contains oppositely charged extremal horizons and has Melvin asymptotics;
it is not a traversable throat. Our dimensionless calibration is $Q/D=0.2$,
where $D$ is the Weyl coordinate separation, not the proper horizon separation.
The motivation comes from [Maldacena, Milekhin and Popov](https://arxiv.org/abs/1807.04726):
charged massless fermions can provide negative Casimir energy in their
four-dimensional traversable-wormhole construction, which does not violate
ambient causality. That construction does not validate this different, static
exterior.

An earlier project energy ansatz assumed fixed exterior optical lengths,
closed antiperiodic Dirac channels, transparent open channels in their
infinite-line vacuum, and a positive classical opening cost. It has a finite
minimum only when the closed-channel fraction $f>1/4$. This is a conditional
statement about that ansatz, not a wormhole-support theorem. Classical flux
connectivity cannot simply be substituted for $f$. That missing quantum
connection motivates the present calculation. Tunnelling between Dirac
zero-mass lines is established subject matter; see
[Tudorovskiy and Katsnelson](https://arxiv.org/abs/1206.2244).

[Anglin and Schulz](https://arxiv.org/abs/1612.05126) analytically solve a
right-angle channel splitter for a checkerboard step mass. Its mass profile
differs from the smooth saddle below; the splitting mechanism itself is
established physics.

In natural units, for one chirality, quantized flux $q>0$ and half-integer
azimuthal label $\mu$, write the metric as
$ds^2=-N^2dt^2+P^2(d\rho^2+dz^2)+\rho^2N^{-2}d\phi^2$ and set
$c=N/P$, $b=N^2/\rho$, $a=A_\phi/q$. In the flat-norm spinor representation,
the physical Hamiltonian is

$$
H=-i\sigma_x(c\partial_\rho+\tfrac12c_\rho)
  -i\sigma_z(c\partial_z+\tfrac12c_z)
  +\sigma_y b(\mu-qa).
$$

At the equatorial null ring $(\rho_*,0)$, both first derivatives of $a$ vanish.
The source-free Maxwell equation
$\partial_\rho(ba_\rho)+\partial_z(ba_z)=0$ therefore implies
$a_{zz}=-a_{\rho\rho}$ there. With $a_2=a_{\rho\rho}/2>0$,

$$
a=a_*+a_2[(\rho-\rho_*)^2-z^2]+O((|\rho-\rho_*|+|z|)^3).
$$

This Hessian relation is exact; using its quadratic Taylor polynomial is a
local approximation. Balancing kinetic and quadratic terms gives

$$
L_q=\left(\frac{c_*}{qb_*a_2}\right)^{1/3},\qquad
E_q=\frac{c_*}{L_q},\qquad
\eta=\frac{b_*(\mu-qa_*)}{E_q}.
$$

Let $s=(\rho-\rho_*)/L_q$, $t=z/L_q$, $u=(s+t)/\sqrt2$ and
$v=(s-t)/\sqrt2$. The leading dimensionless operator is

$$
H_0=-iA_U\partial_u-iA_V\partial_v+\sigma_y(\eta-2uv),\qquad
A_U=\frac{\sigma_x+\sigma_z}{\sqrt2},\quad
A_V=\frac{\sigma_x-\sigma_z}{\sqrt2}.
$$

The actual operator has coefficient corrections of order $q^{-1/3}$ on
bounded rescaled regions, including the kinetic connection terms. This is
not a bound on the global scattering error. Both derivatives occur at leading
order, so independent frozen transverse slices cannot resolve this crossing.
For the chosen geometry, $L_q\simeq0.940916Dq^{-1/3}$; at $q=1$ the local
approximation is not controlled over the exterior. Changing $q$ at fixed
couplings also requires consistent physical rescaling.

The implemented pilot solves $H_0\psi=0$ on $[-R,R]^2$ using piecewise-linear
discontinuous Galerkin elements. Each face is attached to a **frozen**,
translation-invariant linear-mass lead. Its transverse channels are exact
Hermite functions: each $v$ face has one incoming chiral Gaussian, each $u$
face one outgoing Gaussian, and massive pairs have
$\lambda_j=2\sqrt{Rj}$. In the fixed channel phases, the decaying massive
relation at zero energy is $a_{\rm in}=-i a_{\rm out}$. Zero modes are
longitudinal scattering states. Whole-line normalization and finite-face
Gram matrices are kept distinct.

The output is a complex $2\times2$ matrix from $(v_+,v_-)$ inputs to
$(u_+,u_-)$ outputs. Raw squared amplitudes retain numerical losses; they
are not rescaled to sum to one. The driver separately reports interior jump
loss, incoming-trace mismatch, unresolved outgoing components, current-balance
error and the conditioning of boundary feedback. Mesh, retained-mode and
box-size sensitivities must accompany any routing conclusion.

For $\eta=1$, $R=2.5$, mesh spacing $h=0.0625$ and three massive pairs per
face, the [saved finite calculation](../results/null_saddle_scattering_results.json)
gives, in the stated input/output order and fixed spin phases,

$$
S\simeq i\begin{pmatrix}
-0.988650955 & 0.149611923\\
 0.149611923 & 0.988650955
\end{pmatrix}.
$$

For input $v_+$, the raw output powers are $0.977430711$ and $0.022383728$,
totalling $0.999814439$. DG loss is $0.000167392$ and unresolved outgoing
norm is $0.0000181694$; these are numerical deficits, not physical absorption.
At $\eta=0$, $h=0.125$, each output is $0.499353826$, without normalization
to $0.5$. The reflected $\eta=-1$ case exchanges the preferred output.

The maximum change in a power is $0.00113608$ on refining $h=0.125$ to
$0.0625$, $0.000143830$ on increasing $R$ to $3$ at the fine spacing, and
$0.0000157071$ on adding a fourth massive pair. These are measured
sensitivities, **not certified error bars**. Nine implementation checks pass.
Separate AI-agent checks verified the Hessian reduction, lead phases and
Gram conversion; there has been no external expert review. The exact
reflection relation holds to $2.2\times10^{-15}$ and current/power
accounting to $7.6\times10^{-15}$. A manufactured plane-wave control also
improves under mesh refinement. The result demonstrates routing in this
finite frozen-lead model; it does not determine the infinite-domain matrix.

The remaining scientific gaps are substantial. The frozen leads need a
controlled connection to the true bilinear tails and then to the curved
dihole exterior; finite-channel and discretization error estimates are still
needed. Complex phases and exterior reflections matter for that matching.
A physical vacuum choice, renormalized stress tensor, self-consistent
gravitational backreaction, throat attachment, formation and stability remain
unproved. No chronology-changing trajectory follows, including travel from
13 September 2027 to 12 September 2027.

**Focused question for the professor:** Is treating this locally derived Dirac
saddle as a two-input, two-output scattering node, followed by controlled
matching to the dihole exterior, a mathematically sound and scientifically
worthwhile next step?

From the repository root, using the existing Python environment:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
.\.venv\Scripts\python.exe -B -W error -X utf8 scripts/null_saddle_leads.py
.\.venv\Scripts\python.exe -B -W error -X utf8 scripts/null_saddle_scattering.py
```

The first command checks the analytic lead implementation; the second runs
the finite-box driver. The derivation and prior assumptions are documented
in [the research record](RESEARCH.md), with the
[local scaling calculation](../scripts/dihole_null_scaling.py),
[lead implementation](../scripts/null_saddle_leads.py), and
[scattering driver](../scripts/null_saddle_scattering.py).
The saved results include source, input and field hashes. The outstanding
readiness gaps are a controlled lead-to-exterior matching procedure and
continuum error estimates, before any physical return or support claim.

A [subsequent finite-arm calculation](../results/null_saddle_matched_results.json)
propagates the changing Hermite widths/centers and the source/detector maps.
At the same fine mesh with four pairs it gives powers 0.977383541 and
0.022443584. Ten coupled checks pass; a separate short-transfer audit also
agrees. Tested remote-cap sensitivity is small; mesh error is larger. The
four full-line arm domains overlap at corners, so this improvement still
does not certify the complete plane or curved exterior. See the
[updated derivation](RESEARCH.md) for assumptions and quantitative comparisons.
