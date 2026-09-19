# New attempt: test the onset of time loops, not just the throat radius

**Verdict:** An exact classical time-loop background passes the vacuum field
equations, but fails a bounded test-field response diagnostic. The tested repairs
do not fix that diagnostic. This is **not** a completed time-travel medium, and
it is **not** a theorem excluding every time-machine proposal.

## What changed from the previous attempt

The previous membrane calculation tested a wormhole's uniform radius. Its local
bending calculation then produced a negative kinetic coefficient. That result
was a probe warning, not a full coupled-gravity ghost theorem. Retuning its
radius would not settle the missing question.

Here the supporting configuration is instead a **vacuum gravitational field**
in the compact-core family introduced by Amos Ori [1]. The background contains
no exotic membrane in its core and has a chronology-respecting region followed
by an extension with closed timelike curves. This is closer to the requested
"initially no time travel; time travel develops later" scenario than an eternal
rotating universe. It does not by itself supply a factory, an Earth exterior,
or a selectable one-day journey.

The original geometry is published work, not our invention. The calculation
below independently checks it and constructs an exact scalar-field diagnostic
and a trial modification. Novelty of the diagnostic has not been established.

## 1. Exact vacuum background

Use dimensionless coordinates `(T,x,y,z)`, with `z` periodic with period `2*pi`,
and metric signature `(-,+,+,+)`:

\[
 ds^2=dx^2+dy^2-2\,dT\,dz+
 \left[\frac a2(x^2-y^2)-T\right]dz^2.
\]

An overall length scale can multiply this metric. None of the numerical times
below has been identified with seconds, a calendar day, or a year.

`verify_geometry.py` calculates the connection and all Ricci components directly.
It finds, exactly,

\[
 \det g=-1,\qquad R_{ab}=0.
\]

This verifies the vacuum Einstein equation inside the core. It does not solve
the surrounding material region or its evolution. The original paper explicitly
leaves the surrounding material's specific field description and the stability
question unresolved [1].

Define the chronology-adapted coordinate used in [1],

\[
 t=T-\frac a2(x^2-y^2)+e(x^2+y^2).
\]

Direct calculation gives

\[
 (\nabla t)^2=t+[(2e-a)^2-e]x^2+[(2e+a)^2-e]y^2.
\]

For `e=1/16` and `0<a<1/8`, both square-bracket coefficients are negative.
Thus `t` is a time function throughout the `t<0` region: future-directed causal
curves cannot close there. This parameter choice directly satisfies the needed
inequalities; it need not use every auxiliary sufficient restriction in [1].

On a circular path at fixed `t,x,y`,

\[
 ds^2=[e(x^2+y^2)-t]dz^2.
\]

At the central circle, it becomes null at `t=0` and timelike for `t>0`.
The determinant never vanishes. Thus a smooth extension of this idealized
background contains closed timelike curves. Ordinary Cauchy evolution only
establishes the chronology-respecting domain up to its boundary, not unique
predictions for everything beyond the Cauchy horizon.

**Important:** This is a statement about a specified mathematical spacetime,
not a demonstration of preparing its initial data in a laboratory.

## 2. A smooth field disturbance that can be arbitrarily small initially

Test two real canonical scalar fields, written as one complex combination
`phi=phi1+i*phi2`. Both have positive kinetic energy. A single real component
also suffices for an unbounded-energy subsequence; the pair makes the energy
formula nonoscillatory.

Let

\[
 S=T+b x^2+d y^2,
\]

where

\[
 4b^2-b-a/2=0,\qquad 4d^2-d+a/2=0.
\]

Choose the small roots

\[
 b=\frac{1-\sqrt{1+8a}}8,\qquad
 d=\frac{1-\sqrt{1-8a}}8.
\]

The four-dimensional metric then gives the exact identities

\[
 (\nabla S)^2=S,\quad
 \Box S=\kappa=1+2b+2d,\quad
 \nabla S\cdot\nabla z=-1,\quad
 (\nabla z)^2=\Box z=0.
\]

Consequently the complete massless scalar equation closes on the family

\[
 \phi=e^{ikz}u(S),\qquad
 S u''+(\kappa-2ik)u'=0,
\]

where integer nonzero `k` ensures periodicity. The exact nonconstant solution is

\[
 \phi=\epsilon e^{ikz}(-S)^{s+2ik},\qquad
 s=-2(b+d)=\frac{\sqrt{1+8a}+\sqrt{1-8a}-2}{4}.
\]

Concavity of the square root implies `s<0` for `0<a<1/8`; the flat `a=0`
limit has `s=0`. The symbolic script verifies the full four-dimensional wave
operator on this field, not merely an assumed ordinary differential equation.

The field is smooth wherever `S<0`. On `t=-1` in a bounded core, `S` is bounded
away from zero and the field and its derivatives are finite. Multiplying it by
arbitrarily small `epsilon` makes the initial disturbance arbitrarily small.
The exact solution is not claimed to have finite total energy on the untruncated
infinite background. It is a compact-core test; exterior completion and the
nonlinear gravitational constraints remain separate tasks.

A numerical integral on `t=-1`, `x^2+y^2<=1`, and `0<=z<2*pi`, for `a=.02`,
`k=1`, gives initial scalar energy `19.50949564166929 * epsilon^2` in normalized
model units. Two quadrature resolutions agree to displayed precision. The
induced volume and the normal-observer energy, not coordinate energy alone,
are used in this check.

## 3. The divergence is measured by a regular observer

On the central axis `x=y=0`, take a timelike observer with

\[
 u^a=\frac{(1,0,0,1)}{\sqrt{2+T}},\qquad -1\leq T<0.
\]

This vector has norm `-1`. It follows the path `z=T+constant` and reaches the
central null circle in finite proper time. Its acceleration remains finite:

\[
 a^a a_a=\frac{(T+3)^2}{4(T+2)^3}
 \longrightarrow \frac9{32}.
\]

Thus the effect below is not an observer made infinitely accelerated by an
unfortunate coordinate choice.

For the two canonical scalar components, the energy measured by this observer is

\[
 \rho=\epsilon^2(-T)^{2s}
 \left[
 \frac{k^2(T+2)}{T^2}
 +\frac{s^2(T^2+2T+2)}{2T^2(T+2)}
 \right].
\]

In particular,

\[
 \rho\sim\epsilon^2\left(2k^2+\frac{s^2}{2}\right)(-T)^{2s-2}
 \quad\hbox{as }T\to0^-.
\]

The exponent is at most `-2`. For any nonzero amplitude and `k`, the measured
energy diverges in this fixed-background solution. The field need not start
with infinite energy or a singularity on the initial compact region.

**What this proves:** there is an explicit small-initial-data disturbance that
does not remain uniformly bounded as this background approaches the central
null circle. Therefore this background fails the requested bounded test-field
stability diagnostic.

**What this does not prove:** the final outcome after the disturbance's gravity
is included. The stress enters the Einstein equation at second order in
`epsilon`. The test-field approximation must cease to be valid as that stress
becomes important. A nonlinear coupled solution might form a singularity,
change the horizon, or behave differently. We have not simulated that outcome.
Related instability results in the flat Misner limit are known [2]; their
nonlinear conclusions are not silently transferred to this curved core.

## 4. First repair: change the curvature strength

The exponent formula covers the entire interval `0<a<1/8`, not only sampled
points. Since `s<=0`, retuning `a` inside this family does not remove the
massless instability witness. The vacuum equations and causal change can
continue to hold while this stability diagnostic fails.

## 5. Second repair: give the disturbance a positive mass

The massive scalar equation remains exactly separable in the same family:

\[
 S u''+(\kappa-2ik)u'-m^2u=0.
\]

Set `s_log=-log(-S)` and `p=s+2ik`. The equation integrated numerically is

\[
 u_{s_{log}s_{log}}+p u_{s_{log}}+m^2e^{-s_{log}}u=0.
\]

The mass term becomes subleading near the horizon. An independent exact basis is

\[
 u_1={}_0F_1(;1-p;m^2S),\qquad
 u_2=(-S)^p\,{}_0F_1(;1+p;m^2S).
\]

The second solution retains the problematic near-horizon exponent for every
finite mass. There is a fine-tuned regular branch, but stability against general
small disturbances cannot be obtained by discarding the other allowed branch.

### Executed numerical tests

- `a`: 0, .001, .01, .02, .04, .06, .08, .10, .12.
- `k`: 1, 2, 4, 8.
- `m`: 0, .1, 1, 10.
- 144 mode evolutions, of which 128 use nonzero curvature; 36 are massless and
  108 massive.
- Start at `S=-1`; stop at `S=-10^-12`. Nothing is integrated through the horizon.
- Initial data: `u=1`, `du/ds_log=-p`. These are diagnostic seeds, not randomly
  sampled universe histories.
- DOP853 tolerances: relative `5e-11`, absolute `5e-13`.
- All modes compared against a separately evaluated 60-digit hypergeometric
  solution at four sample times.
- Selected massive curved-background run also checked with implicit Radau.

Results from the executed run:

| Check | Result |
|---|---|
| Energy amplification range | `1.18893e23` to `5.02358e26` |
| Largest scaled discrepancy against exact basis | `1.424e-9` or less |
| Largest asymptotic energy-slope discrepancy | `1.50e-10` or less |
| DOP853/Radau scaled state discrepancy | `3.70e-11` or less |
| Regular-branch control | Bounded; energy falls from 1 to about 0.5 |

For `a=.02`, `k=1`, `m=0`, the amplification is about `2.18646e24`.
These are ratios in the linear model. They are **not** energy predictions for
a real apparatus: extrapolating normalized unit-amplitude test fields to huge
energies ignores the point where backreaction invalidates the approximation.

No computation here is called a full spacetime simulation. These are ODE
evolutions of an exact invariant field sector on a fixed four-dimensional metric.

## 6. Third repair: a variable, positive kinetic coupling

Try to suppress the growth using

\[
 S_{\phi}=-\frac12\int\sqrt{-g}\,Z(S)(\nabla\phi)^2\,d^4x,
 \qquad Z(S)=(-S)^{-q}>0.
\]

This is deliberately a generous diagnostic: `Z` is treated as a prescribed
scalar background. Its own field equation, source and gravity have **not** been
supplied. Therefore it is not yet a complete covariant material model. It
already fails the desired response test even under that favorable simplification.

The mode equation becomes

\[
 S u''+(\kappa-q-2ik)u'+\frac{iqk}{S}u=0.
\]

Power-law exponents satisfy

\[
 \lambda^2-(p+q)\lambda+iqk=0.
\]

For a nonconstant mode, the leading measured-energy exponent is

\[
 E_j=2\operatorname{Re}\lambda_j-q-2.
\]

The sum of both exponents is fixed:

\[
 E_1+E_2=2s-4<0.
\]

For `q!=0`, neither root vanishes, so at least one physical, nonconstant branch
has divergent energy for **every real q**. At `q=0` there is a constant regular
branch but the original divergent branch remains. Repeated roots do not remove
the negative power and may add logarithms.

This is an analytic obstruction within the proposed power-law repair, not merely
a failed parameter search. A supplementary algebraic check of 2,916 combinations
(`9 a * 4 k * 81 q`) found a nonconstant divergent branch in all of them.
These are algebraic checks, **not** 2,916 time-domain simulations.

## 7. Other simple changes do not alter this linear diagnostic

Adding a regular positive quartic scalar potential around `phi=0` has no linear
term in its field equation. It therefore cannot eliminate this linear
instability witness merely by making the potential stiffer. Its nonlinear
outcome has not been calculated.

A conventional nonminimal term `xi R phi^2` also does not change the linear
scalar equation on this vacuum background, because `R=0`. This observation is
restricted to that term and background; it is not a statement about all modified
gravity theories or higher-curvature couplings.

## 8. What survives and what fails

**Established here:**

1. The published compact-core metric satisfies the exact vacuum equations.
2. Its negative-time region admits a time function, while an extension contains
   closed timelike circles.
3. A specified, initially smooth test-field family has unbounded energy near the
   central chronology-onset circle as measured by a regular timelike observer.
4. Curvature retuning, a finite scalar mass, and the tested positive power-law
   kinetic coupling do not provide robust bounded response in this family.

**Not established:**

- Fully nonlinear gravitational evolution of these disturbances.
- A stable final time machine or a theorem excluding every final CTC geometry.
- Formation in a laboratory from accessible matter.
- A specified, healthy material source for the exterior or for variable `Z`.
- A path to the requested Earth date, inclusion of a human, or safe passage.

The intended success condition has **not** been met. The value of this attempt
is a concrete failure mechanism at time-loop formation, rather than another
assumed time-shift formula. A replacement must address that mechanism within a
complete coupled field theory, not delete the unstable mode by fiat.

## 9. Reproduce the work

From this directory, with Python, SymPy, NumPy, SciPy, and mpmath available:

```bash
python verify_geometry.py
python check_initial_energy.py
OPENBLAS_NUM_THREADS=1 python run_modes.py
python test_kinetic_repair.py
python verify_outputs.py
```

The numerical results record exact package versions and tolerances. The JSON
outputs contain all cases, not only favorable examples. The scripts have no
network dependency and do not use any external credentials.

## Primary sources

[1] A. Ori, *A Class of Time-Machine Solutions with a Compact Vacuum Core*,
Physical Review Letters 95, 021101 (2005).
https://doi.org/10.1103/PhysRevLett.95.021101
Preprint: https://arxiv.org/abs/gr-qc/0503077
The core metric, its causal structure, and unresolved source/stability issues
come from this paper. The numerical outputs here are our own computations.

[2] P. Denaro and G. Dotti, *Strong cosmic censorship and Misner spacetime*,
Physical Review D 92, 024017 (2015).
https://arxiv.org/abs/1505.01092
Relevant comparison for the flat limit; not used as a substitute for a nonlinear
analysis of the curved case.

[3] S. Mironov, V. Rubakov and V. Volkova, *In hot pursuit of a stable wormhole
in beyond Horndeski theory*, Physical Review D 107, 104061 (2023).
https://arxiv.org/abs/2212.05969
Screened as an alternative source theory. Its stated high-momentum stability is
not full stability or a time-machine construction.

[4] S. Mironov and V. Volkova, *Complete stability for spherically symmetric
backgrounds in beyond Horndeski theory*, version 2 (17 January 2025).
https://arxiv.org/abs/2404.06297v2
Gives stronger sufficient criteria and explicitly notes that previous examples
do not meet the full set. Failing a sufficient criterion alone would not prove
instability; we did not claim otherwise or simulate those backgrounds here.

[5] J. Dietz, A. Dirmeier and M. Scherfner, *Geometric Analysis of Ori-type
Spacetimes* (2012). https://arxiv.org/abs/1201.1929
Further mathematical context for the geometry, not a stability certificate.
