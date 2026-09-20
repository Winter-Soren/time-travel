# Fiction-inspired time travel: mechanisms and two diagnostic calculations

**15 September 2026 — bounded research note.** These stories suggest useful
questions about transport, information and consistency. They supply no physical
medium. **A stable tiny-probe route from 13 September 2027 to 12 September 2027
remains unproved.** The two calculations below are small, reproducible checks,
not construction proposals.

Context: the authoritative top of `TASK_STATE.md` and latest available opening
section of `docs/RESEARCH.md` were read. The current source work remains conditional:
no complete, stable throat-plus-exterior solution has been established. This
note adds information/consistency diagnostics without changing that work.

## Six story mechanisms and their physical requirements

The transport apparatus and the rule for changing history are separate ideas;
the table includes both. Mathematical analogues below are **our comparisons**,
not explanations supplied or validated by the stories. Physics references
P1–P6 are linked immediately after the table.

| Story mechanism and verified reference | What moves or changes in the story | Nearest mathematical analogue; what a model would need | Missing physical assumptions |
|---|---|---|---|
| **Marvel comics: Doom's time platform.** *Fantastic Four* #5 (1962): Doom sends Reed, Ben and Johnny to the pirate era for Blackbeard's treasure. [Marvel retrospective](https://www.marvel.com/articles/comics/kirby-100-meet-doctor-doom). | **Matter:** whole travelers and carried objects participate in the earlier period. | An event-to-event transport map; a traversable wormhole with a suitable mouth time offset is one conditional GR analogue (P1). The story does not identify its platform as that geometry. | A derived metric and source, preparation, accessible endpoints, causal orientation, payload survival and stability. Naming an arrival date specifies none of them. |
| **Marvel comics: changing a counterpart's history.** *Marvel Two-in-One* #50 (1979): Ben cures an earlier counterpart but remains the Thing. Reed explains that intervention creates another reality. [Original ending page, visually checked](https://tombrevoort.com/wp-content/uploads/2024/08/m2in1-v1-050-p033.jpg), [editor Tom Brevoort's account](https://tombrevoort.com/2024/08/04/bhoc-marvel-two-in-one-50/). | **Matter plus an intervention:** Ben and the cure travel; the counterpart's future differs without curing the visitor. | Histories with a branch label, or quantum relative states as a limited analogy (P2). A transport law would have to connect histories, not merely name them. | Physically accessible alternate histories and an interaction that selects one. Attribute this rule specifically to **#50**, not every Marvel continuity. It does not establish alteration of the departure history. |
| **Steins;Gate: D-Mail.** In the 2011 anime, episode 6, “Butterfly Effect's Divergence,” develops past-directed messages. [Licensed episode listing](https://www.crunchyroll.com/watch/GRWE7W5QR/butterfly-effects-divergence); the publisher separately confirms the modified-microwave e-mail premise for *ELITE*. [Spike Chunsoft](https://www.spike-chunsoft.com/games/steinsgate-elite/). | **Information:** a message influences earlier decisions; a body need not be transported. | A communication channel whose earlier receiver statistics depend on a later freely chosen message. Ordinary causal quantum operations provide a null benchmark (P3); postselection provides a potentially misleading imitation (P6). | An actual backward channel with nonzero usable capacity. Entanglement, a microwave, or calling a field a black hole does not supply it. See calculation 1. |
| **Steins;Gate: memory time leap.** Episode 11, “Dogma in Event Horizon,” introduces transport of memories. [Licensed synopsis](https://tv.apple.com/ca/episode/dogma-in-event-horizon/umc.cmc.75fk12vx6vqqi7odbbvfyj05o?showId=umc.cmc.7bzye9l3t2ww0rpuv91symto7), [Crunchyroll description](https://www.crunchyroll.com/watch/GYNQ4NG0Y/dogma-in-event-horizon). | **Encoded memories**, rather than the traveler's whole body. This note makes no claim about transfer of an entire unknown quantum brain state. | An encode → channel → decode/write process. Quantum teleportation transfers a quantum state using entanglement **and classical communication** (P4), but does not provide the backward channel. | A valid memory encoding and restoration process, sufficient capacity/fidelity, and backward delivery. Classical records can be copied; universal perfect copying of unknown quantum states is forbidden in ordinary linear quantum mechanics (P5). |
| **Rick and Morty: interdimensional relocation.** S1E6, “Rick Potion #9”: Rick and Morty replace dead counterparts in another reality and leave the damaged reality behind. [Official clip](https://www.youtube.com/watch?v=FIbIwA0koBg), [episode metadata](https://tv.apple.com/us/episode/rick-potion-9/umc.cmc.20nmj9mj6vbl1v6kkvci2zd88?showId=umc.cmc.12dp30hnvyq5fbm9716puu8zc). | **Matter and retained experience:** living travelers occupy another history. | A spatial wormhole is an analogy for traversing an aperture (P1); a branch label is an analogy for the destination (P2). Neither supplies the other. | A physical connection to an accessible alternate universe, destination selection and a specified time relation. A green portal and a different history do not establish a negative time offset. |
| **Rick and Morty: time fracture and recombination.** S2E1, “A Rickle in Time”: uncertainty fractures time; the resolution depicts multiple versions and collar-assisted recombination. [Official fracture description](https://www.adultswim.com/videos/rick-and-morty/time-fractures), [official resolution clip](https://www.youtube.com/watch?v=9aUowoykENE). | **Macroscopic world/character states:** multiple depicted versions resolve to one scene, not just a transferred message. | Coherent superposition and interference provide a narrow quantum analogy (P2). A model would need explicit dynamics for splitting, preserving coherence and recombining the relevant degrees of freedom. | Control over the characters and their environment, plus a physical coupling and consistent records. Human indecision is not an established macroscopic branch-control mechanism. Panel counts are not measured quantum states. |

### Physics references and distinctions

- **P1 — Conditional spacetime geometry:** [Morris, Thorne and Yurtsever,
  *Phys. Rev. Lett.* 61, 1446 (1988)](https://doi.org/10.1103/PhysRevLett.61.1446).
  A spatial wormhole is not automatically a time machine. A suitable time
  offset and an exterior return path must make a closed causal route possible.
  A closed **timelike** curve is a future-directed massive-particle worldline
  returning to its starting event. The paper's argument is conditional on
  creating and maintaining a traversable wormhole.
- **P2 — Relative states:** [Everett, *Rev. Mod. Phys.* 29, 454
  (1957)](https://doi.org/10.1103/RevModPhys.29.454). A decomposition such as
  `|Psi> = sum_i c_i |record_i>|environment_i>` is not a portal prescription.
  Our inference: writing a branch index or a superposition supplies no
  interaction that lets an observer choose and enter another universe.
- **P3 — Causal operations:** [Beckman, Gottesman, Nielsen and Preskill,
  *Phys. Rev. A* 64, 052309 (2001), Sections I–II](https://preskill.caltech.edu/pubs/preskill-2001-causal.pdf).
  Their distinction between trace-preserving operations and selected outcomes
  motivates calculation 1. Its temporal application and example are derived here.
- **P4 — Teleportation:** [Bennett et al., *Phys. Rev. Lett.* 70, 1895
  (1993)](https://doi.org/10.1103/PhysRevLett.70.1895). The receiver requires the
  sender's classical measurement result. Entanglement alone is not a usable
  earlier-delivery service; transferring a state is not transporting its original matter.
- **P5 — No universal quantum copier:** [Wootters and Zurek, *Nature* 299,
  802–803 (1982)](https://www.nature.com/articles/299802a0). This restriction
  concerns arbitrary unknown quantum states, not ordinary copying of known
  classical records. It is not a general proof against time travel.
- **P6 — Different hypothetical loop rules:** [Deutsch, *Phys. Rev. D* 44,
  3197 (1991)](https://doi.org/10.1103/PhysRevD.44.3197) and
  [Lloyd et al., *Phys. Rev. Lett.* 106, 040403 (2011), author manuscript,
  Eq. 1 and pp. 1–2](https://arxiv.org/pdf/1005.2219).
  The latter's postselected experiment simulates a proposed consistency rule.
  Selecting successful trials later is not the physical power to force their
  outcome or send an operational message earlier.

Neither a named fictional substance, a Möbius shape nor a branching diagram
is evidence for a medium. Likewise, the
[quantum switch](https://arxiv.org/abs/0912.0195) controls the order of operations;
that fact supplies no day-long backward communication channel.

**Canon scope and evidence limits.** The Marvel rows refer to comics, not film
adaptations; no MCU Quantum Realm or Pym-particle premise is imported.
The Steins;Gate rows use the 2011 anime references; *ELITE* is separately
identified corroboration. Exact byte budgets, leap limits and a universal
worldline rule are omitted because they were not established from the accessed
primary material. Crunchyroll's indexed official descriptions were readable,
but direct opens returned an unsupported-browser page; Apple TV and the game
publisher provided accessible corroboration. Marvel's retrospective was
available through indexed official text but returned HTTP 403 on direct open.
The #50 ending itself was read visually. The Rick and Morty checks used
official clip captions and selected frames, not full-episode viewing; no
precise collar operating law is inferred.

## Candidate calculation 1 — Can the earlier record carry a later choice?

**Priority: recommended as a lightweight acceptance test for any information-only
proposal.** This checks a necessary operational property, not how to build it.

**Assumptions:** ordinary forward-ordered quantum mechanics, the Born rule,
an earlier recorded outcome `a`, and a later independently randomized choice
`b`. The preparation and early measurement do not depend on `b`. Include all
outcomes, losses and failures. Let a complete early instrument `{M_a}` produce
the unnormalized state `sigma_a = M_a(rho)`, with `sum_a Tr(sigma_a)=1`.
A complete later operation `E_b` is
completely positive and trace-preserving (CPTP). Then

\[
p(a\mid\mathrm{do}(b))=\operatorname{Tr} E_b(\sigma_a)
=\operatorname{Tr}\sigma_a,\qquad
D_{\rm early}=\frac12\sum_a|p(a\mid\mathrm{do}(0))-p(a\mid\mathrm{do}(1))|=0.
\]

Here `do(b)` means an independently set choice, not a value chosen after
examining `a`. Arbitrary forward propagation can be included in `E_b`.
This is our direct trace-preservation derivation using P3's framework; it
does not cover every proposed modification of quantum mechanics.

A later selected outcome `s` instead gives

\[
p(a\mid b,s)=\frac{\operatorname{Tr}E_{s\mid b}(\sigma_a)}
{\sum_{a'}\operatorname{Tr}E_{s\mid b}(\sigma_{a'})}.
\]

**Computable negative control:** choose independent fair bits `a,b`, then retain
only trials with `s=1` iff `a=b`. Exactly `p(s=1|b)=1/2`,
`D_selected=1`, but `D_early=0`. The apparent perfect prediction exists only
after the future filter is known. Even constant acceptance probability does
not remove this selection bias.

**Output and consistency test:** report `D_early`, `D_selected`, acceptance
rates and the full trial counts separately. For a numerical quantum proposal,
check completeness `sum_s K_(s|b)^dagger K_(s|b)=I` and positivity before
interpreting a signal. The exact negative control must recover `(0,1,1/2)`.
For a real experiment, preregister the statistic, freeze the early records
before choosing `b`, and test the unfiltered data against independently
randomized labels. In a finite experiment or Monte Carlo simulation, empirical
TV is usually positive under the null; positivity alone is neither an error
nor evidence of a backward signal. An excess beyond expected sampling and
numerical error demands an audit of timing, leakage, exclusions and independence.
A strictly positive TV between exact distributions would contradict the stated
assumptions. A null result rejects that channel claim
to the test's sensitivity, not all possible theories of retrocausality.

**13 → 12 September 2027:** put the early receiver record on 12 September
and choose the message on 13 September, at the same local clock time on both
days. The required lag is `-86,400 s`. Substituting these dates in the equations
creates no channel. Even a verified information channel would leave transport
of the tiny abiotic probe unproved.

## Candidate calculation 2 — Which loop consistency rule, and which stability?

**Priority: secondary conceptual check.** Useful if a proposed probe controller
eventually couples to a backward path. It does not advance the source solution
by itself. P6 supplies alternative hypothetical consistency rules; the minimal
unmonitored loop and noise calculation here are our explicit examples.

For a loop qubit with no interacting exterior system, let `X` be the Pauli
bit-flip matrix. A Deutsch-style boundary condition is

\[
\tau=X\tau X,\qquad
\tau=\frac{I+r_xX}{2},\quad -1\leq r_x\leq1.
\]

The consistency equation alone admits infinitely many fixed states, including
pure `X` eigenstates, in this **unmonitored** example. Deutsch's additional
maximum-entropy selection postulate would select `I/2` from this family;
that is an extra assumption, as distinguished by
[Lloyd et al., p. 3](https://arxiv.org/pdf/1005.2219#page=3). By contrast, a
definite classical bit cannot satisfy `a=1-a`. Simply iterating the quantum
map from `|0><0|` alternates between two states: existence is not convergence.

**Our added noise assumption:** for `0<epsilon<=1`, define the ordinary
depolarizing extension

\[
F_\epsilon(\tau)=(1-\epsilon)X\tau X+\epsilon I/2.
\]

Its Bloch-vector linear part has eigenvalues
`(1-epsilon, -(1-epsilon), -(1-epsilon))`; the unique fixed state is `I/2`
and the contraction factor is `1-epsilon`. This is stability of a specified
iteration. The noise is an assumed operation, not a derived stabilizer of
a spacetime geometry or a method to preserve a probe's stored information.

**Compare the noiseless postselection rule:** prepare a Bell pair, act on one
half with `U_theta=exp(-i theta X)`, then project back onto that same Bell state.
There is no exterior input coupling in this example. The accepted amplitude is

\[
k=\tfrac12\operatorname{Tr}U_\theta=\cos\theta,\qquad
p_{\rm accept}=|k|^2=\cos^2\theta.
\]

At `theta=pi/2`, conjugation is the same bit flip used above, but acceptance
is zero: the normalized postselected output is undefined. Do not replace a
zero denominator by a fallback state. A Deutsch fixed state and a nonzero
postselection probability are different requirements. This Bell-overlap
calculation is a forward experiment unless a backward link is independently
supplied. Adding an uncoupled exterior register gives `K=cos(theta) I_ext`:
whenever accepted, its normalized output equals its input. This particular
example carries no message even within the conditional model.

**Outputs and consistency tests:** fixed-point residual, multiplicity,
iteration eigenvalues and acceptance probability. Test pure and mixed initial
states, recover nonuniqueness at zero noise, and preserve the exact zero-success
case. The check below gives a maximum noisy iteration residual `3.56e-14`
after 300 steps at `epsilon=0.1`; the bit-flip Bell acceptance is analytically
zero. Neither result measures gravitational or material stability.

**13 → 12 September 2027:** assigning the two ends of the consistency condition
to those dates would *assume* the missing backward connection. The equations
contain no duration, stress-energy, apparatus or delivery mechanism and cannot
establish that date route. They expose incompatible model assumptions before
we spend resources on a larger conditional simulation.

## One runnable check for the two examples

Run the following Python block with the project's Python, or any ordinary
Python 3 interpreter. It uses only the standard library and writes no files.

```python
from fractions import Fraction as F
from math import cos, pi, isclose

# Independent early record a and late choice b; filter only after both exist.
joint = {(a, b): F(1, 4) for a in (0, 1) for b in (0, 1)}
raw = [[joint[a, b] / F(1, 2) for a in (0, 1)] for b in (0, 1)]
keep = [sum(joint[a, b] for a in (0, 1) if a == b) / F(1, 2)
        for b in (0, 1)]
selected = [[(joint[a, b] if a == b else 0) / (keep[b] * F(1, 2))
             for a in (0, 1)] for b in (0, 1)]
tv = lambda p, q: sum(abs(x - y) for x, y in zip(p, q)) / 2
assert tv(*raw) == 0 and tv(*selected) == 1
assert keep == [F(1, 2), F(1, 2)]

# Algebraic Bloch-vector iteration, with no physical backward connection.
def step(r, eps):
    return tuple((1 - eps) * sign * v for sign, v in zip((1, -1, -1), r))

assert step((1, 0, 0), 0) == (1, 0, 0)
assert step((-1, 0, 0), 0) == (-1, 0, 0)
assert step((0, 0, 1), 0) != (0, 0, 1)
assert step(step((0, 0, 1), 0), 0) == (0, 0, 1)
residual = 0.0
for seed in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0.2, 0.3, 0.4)):
    r = seed
    for _ in range(300):
        r = step(r, 0.1)
    assert max(map(abs, r)) < 1e-12
    residual = max(residual, max(abs(a - b) for a, b in zip(step(r, 0.1), r)))
assert cos(pi / 2)**2 < 1e-30  # Analytically zero; floating-point pi is approximate.
assert isclose(cos(pi / 3)**2, 0.25, abs_tol=1e-15)
print("PASS: raw TV=0; selected TV=1; acceptance=1/2")
print(f"PASS: maximum noisy loop residual={residual:.3e}")
```

**Research recommendation:** retain calculation 1 as an early rejection test
for information-only proposals; use calculation 2 to make any chosen loop
rule explicit. Neither replaces the quantum-source and geometry work. Ordinary
future time dilation remains the permitted fallback; arriving at an earlier
date is a separate, still unsolved objective.
