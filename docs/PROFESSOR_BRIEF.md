# Technical assessment: a bounded Dirac response in a balanced dihole exterior

**Proposed brief for first methodological feedback.** We have an a posteriori
numerical enclosure for one compact quantum-field response in a specified curved
exterior. This is an AI-assisted calculation with independent agent checks,
not external peer review. No novelty, supporting material, stable wormhole or
chronology change is claimed.

## Physical question and background

The motivation is whether a controlled exterior fermion response can eventually
supply useful input to a quantum stress or throat-matching calculation.
[Maldacena, Milekhin and Popov](https://arxiv.org/abs/1807.04726) construct a
four-dimensional traversable wormhole supported by charged-fermion Casimir
energy. Their result does not validate the different geometry used here.

Our background is the externally balanced magnetic Einstein–Maxwell solution
in [Emparan, *Black diholes*, equations8–16](https://arxiv.org/pdf/hep-th/9906160).
It has two oppositely charged extremal horizons and Melvin asymptotics. It is a
known classical black-hole exterior, not a traversable throat. Set the Weyl
coordinate separation D=1 and charge radius Q=1/5. The geometric parameters
are linked exact radicals; stored floating values are not treated as independent
exact parameters. Changing quantized flux at fixed Q/D requires consistent
coupling calibration or physical rescaling.

For one massless chirality and half-integer angular label mu, the canonical
flat-norm Hamiltonian is

$$
H_\mu=-i\sigma_x(v\partial_\rho+v_\rho/2)
      -i\sigma_z(v\partial_z+v_z/2)
      +\sigma_y\frac{N^2}{\rho}(\mu-q a),\qquad v=N/P.
$$

Here a is the charge-included gauge coefficient divided by q. The
[natural-domain proof](DIHOLE_DIRAC_DOMAIN.md) includes the regular axis and
magnetic gauge patches. The optical metric is complete: infinity and the
horizon ends have infinite optical distance. The conformal spinor map is
unitary, and a [complete Dirac-bundle theorem](https://arxiv.org/pdf/1606.06190)
gives essential self-adjointness. The lifted axial symmetry supplies the
reducing fixed-mu sector. This conclusion concerns the uncut exact exterior;
it does not extend automatically to a replacement throat.

## Concrete result

A local saddle reduction near the magnetic null supplied an initial trial.
Two corrections in the actual curved coefficients then reduced its residual.
The present result encloses that trial's full continuum residual directly;
it does not infer accuracy solely from basis or grid convergence.

Use q=512000, mu=398908.5 and the compact unnormalized Gaussian source h, with
cutoff5/7 in each rescaled coordinate. The trial has128 Hermite modes per
coordinate and cutoff7/9, wholly away from the axis and horizons. The coordinate
scale and center are specified dyadics; the energy scale is the exact ratio
Eq=v*/L. Define Hq=T^-1 Hmu T/Eq using the unitary coordinate map.
At scaled imaginary frequency i, the response is

$$
G=\langle h,(H_q-i)^{-1}h\rangle.
$$

For A=chi_t a, the exact reflection construction B=-sigma_z A(v,u) supplies
the adjoint trial. The reflected residual has the same norm. The dual identity
gives an exact disk bound with center iC and radius at most

$$
2.785503521209716\times10^{-6}.
$$

The imaginary-axis projection is

$$
\boxed{\operatorname{Im}G\in
[0.41900184413867064,\;0.4190074151728139]}.
$$

The source norm squared is enclosed in
[0.9999999999999585,0.9999999999999833], rather than assumed to equal1.
Physical responses for the unitary-transformed source divide by Eq.
This is an imaginary-frequency resolvent response, not a transmission probability
or a stress component.

All numerical components are covered: exact-parameter evaluation, distinct
interval-validated Gauss nodes and weights, Hermite functions and derivatives,
outward interval contractions, coefficient cancellation, pairwise sums, and
an analytic integration remainder. The certificate assumes IEEE754 elementary
arithmetic with gradual underflow and the documented mpmath interval behavior;
it is not a formal verification of those libraries or hardware.

The full driver takes about23 seconds. Its controls include exact-rational
arithmetic cases, Gauss/Hermite identities, and an independent100-digit
reconstruction of N/P derivatives and the original gauge potential. The
[results](../results/dihole_probe_certificate_results.json) record all source
hashes, parameter definitions, enclosures and limitations. Earlier local
scattering pilots and their unresolved matching errors remain in the
[research record](RESEARCH.md).

## Question for methodological feedback

**Which quantum-state prescription and additional response kernels would make
these certified exterior probes useful for a renormalized fermion stress
calculation or a controlled throat-matching problem?** In particular, we seek
advice on the frequency/angular completion and ultraviolet subtraction before
attempting any self-consistent support claim.

One compact source, one angular sector and one imaginary frequency cannot
determine the vacuum stress. The state, full response kernels, renormalization,
backreaction, throat attachment, formation and stability remain unresolved.
No trajectory to an earlier date follows.

Reproduce from the repository root using the existing environment:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
.\.venv\Scripts\python.exe -B -W error -X utf8 scripts/dihole_probe_certificate.py
```

This runs the helper controls, complete interval integration and progress
figure. No email has been sent; this document is a reviewable technical brief.
