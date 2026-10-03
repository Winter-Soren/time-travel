# Direction assessment and proposed redirect (external review, 3 October 2026)

**Status of this note.** One reviewer, one session, skim-level reading of
`TASK_STATE.md` (checkpoint193), `README.md`, `docs/PROFESSOR_BRIEF.md`,
`docs/ATOMIC_TRANSPORT_TARGET.md`, `docs/THROAT_SOURCE_BUDGET.md`,
`docs/TILTED_CLOCK_ADMISSIBILITY.md` and the saved 193 result. No simulation was
rerun and the interval proofs were **not** audited. Statements about the
literature are from memory (arXiv was unreachable); verify them before relying
on them. This is a recommendation, not a result.

## What the repository has actually established

- A very large body of careful numerical analysis of **one auxiliary number**:
  an error ceiling for the renormalized Dirac-field response in the balanced
  magnetic dihole exterior, at flux $q\approx4.6\times10^{24}$ (1 nm target).
- Ceiling went from $4\times10^{14}$ to $7.8\times10^{-5}$ (checkpoints 191-193)
  and is still ~4000x the (uncertified) numerical reference value.
- The only locally computed supporting piece has the **wrong sign** and is ~10^7
  times too small versus the required $-5.83$ (increment 162 / `README.md`).
- No supporting source, no traversable throat, no stability, no atom transit,
  no clock gain, no new law. The repository says this itself; that honesty is
  its strongest feature.

## Why I recommend changing course

1. **The goal is not on the path.** Even a perfect, stable, traversable throat is
   a *spatial* shortcut. Backward travel needs a second, separate ingredient: a
   time offset between the mouths (twin-paradox style) plus surviving the
   resulting chronology horizon. The repo already records that its static
   geometry has a global time function and cannot give backward passage. So
   finishing the stress certificate cannot, by itself, reach 12 Sept 2027.
2. **Hardest, least decisive step first.** The programme spends its effort on
   rigorous error bounds (checkpoints ~175-193 are almost all bound-tightening)
   before any cheap test of whether the *time-shift* step is compatible with
   the vacuum stress. If chronology protection kills that step, the throat
   certificate is moot. Do kill-tests in order of cost, cheapest first.
3. **Progress metric is misleading.** "1.76 billion times smaller ceiling" is
   measurable progress in a proof technique, not physics progress. The dominant
   remaining uncertainty is the reference quadrature and the omitted
   curved-field/global contribution, not the constants being tightened.
4. **The 1 nm target forces extreme parameters.** Order-of-magnitude check
   (my arithmetic, using the repo's $q$ and $D$ and reading $q$ as $q$ flux
   quanta $h/e$; treat as +/- a factor of a few): flux $\approx1.9\times10^{10}$ Wb
   through a ~1 nm cross-section gives $B\sim10^{28}$ T, about $10^{18}$ times
   the Schwinger field, energy density $\sim10^{62}$ J/m^3, and an exterior
   mass scale $Dc^2/G\sim3\times10^{18}$ kg (asteroid mass in a nm).
   Separately, the repo's own tidal diagnostic is ~280x hydrogen's binding
   energy, and its `q=1` alternative needs $\alpha_R\sim3\times10^{-52}$. The
   payload size is an engineering choice; it should not be the thing that
   dictates where the fundamental physics is tested.
5. **Nothing here is a new law.** All work is established QFT/GR applied to a
   hypothetical background. If new physics is the aim, it has to be stated as a
   concrete postulate or Lagrangian and then tested; none is proposed yet.
6. **The branches that were tried already match known lessons.** The tilted-clock
   scalar branch is excluded by the stable-causality theorem (correctly found
   in `TILTED_CLOCK_ADMISSIBILITY.md`); the Ori-type branch failed a
   test-field diagnostic. Both say: a field that provides a global time
   function cannot make closed causal curves, so any route must go through
   geometry/topology change at a chronology horizon, where semiclassical
   gravity is exactly the open question.

## Is the dihole + Maldacena-Milekhin-Popov idea a sound base?

As a *starting point*, yes: it is, to my recollection, the one 4D traversable
wormhole with a controlled-looking quantum source (charged massless fermions
in a strong magnetic flux). Two caveats to check against the paper: it needs
matter (many light charged fermions) absent from the Standard Model, and, as I
recall, the wormhole is long, so it is not a shortcut relative to the exterior
path, which matters directly for any time-machine construction. As a
*route to a time machine* it is unproven and may be the wrong object.

## Proposed redirect for the automation

**Stop** extending the dihole error-certificate chain (checkpoint 194+ as
currently planned: retaining signed first-order cancellations). Archive it in
`TASK_STATE.md` as "paused: inconclusive, not on critical path". Do not delete
any saved evidence (see `AGENTS.md`).

**Start: a chronology-protection stress test on solvable models**, with gates
that can falsify the whole programme cheaply.

- **Gate 1 (time shift vs. vacuum stress).** In 1+1D and then 2+1D/3+1D
  conformal/free-field toy spacetimes with a wormhole or Misner-type
  identification and a tunable mouth time offset $\Delta T$, compute
  $\langle T_{\mu\nu}\rangle_{\rm ren}$ by image sums/Hadamard methods as the
  offset approaches the CTC onset. Controls: published Misner-space and
  Kim-Thorne-type results; the known divergence at the chronology horizon
  (Kay-Radzikowski-Wald). Output: does stress stay sub-Planckian up to and
  past onset for any parameter window, and with what scaling in
  $\Delta T/$(mouth size)?
- **Gate 2 (quantum-inequality budget).** For any window found, check it
  against Ford-Roman-type quantum inequalities and ANEC-integral requirements
  for the throat; report required vs. allowed negative energy as a scale-free
  ratio.
- **Gate 3 (only if 1 and 2 pass).** Return to a physical source (dihole or
  another) at the scale where the approximation is controlled, not at 1 nm.
- **Gate 4 (only then).** Payload survival and clock comparison.

Rules: every gate gets a pre-stated pass/fail criterion written *before*
running; a negative result is a deliverable; each checkpoint must state what
decision it moves, not only what bound it tightened; keep new status short and
do not grow `TASK_STATE.md` further without pruning.

**Separate, cheap, genuinely achievable side track: future travel.** Forward
time travel is real physics (time dilation; optical clocks resolve ~1 cm of
gravitational redshift). The repo never computed a clock difference. A
two-page calculation of attainable forward offsets for a given trajectory
would give a verified result and an honest benchmark for the "medium" claim.

**Optional lab-scale track.** Quantum simulations of closed-timelike-curve
channels (Deutsch and postselected CTCs, photonic experiments) and of
traversable-wormhole dynamics are real and reproducible, but they simulate
the *rules*; they do not create a medium that carries a person or atom back.
The existing `FICTION_PHYSICS.md` null-benchmark framing is a good start.

## Honest expectation

The most likely outcome of the redirected programme is a documented
chronology-protection obstruction or an inconclusive window. That is still
more informative than another tightened ceiling. A sub-Planckian, stable,
CTC-forming window, if one appeared, would be a notable result and would then
deserve independent expert review before any claim.
