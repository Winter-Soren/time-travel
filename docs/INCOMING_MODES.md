# Incoming spinor modes: a normalized collar prescription

20 September 2026. Bounded follow-up to `docs/RESEARCH.md` increment87. Only this
report, `scripts/incoming_modes.py`, and `results/incoming_modes_results.json` are owned here.
No historical propagation, closure, or field simulation was rerun.

**Result:** an isolated monopole collar supplies exact positive-energy,
unit-flux reference modes, including the spherical `k=0` mode. A closed-form
positive-energy pulse supplies unit *time-integrated* injection. These are
usable incident lead data. They become a physical two-mouth incoming
scattering state only after matching and an outgoing exterior solution.

## 1. Three different objects

| Object | What is fixed | What it does not establish |
|---|---|---|
| Local test packet | Initial spinor and a chosen evolution operator/domain | Positive energy, throat preparation, or incident asymptotics |
| Positive-energy projection | `f_+=P_+ f / ||P_+ f||`, `P_+=1_(0,infinity)(H)` for a specified self-adjoint `H` | Incident lead, unit incoming flux, or scattering boundary conditions |
| Incoming scattering state | Energy, incident channel with unit flux, and no other incident channels, plus outgoing response | Physical preparation of the assumed throat/background |

The increment87 real lower-component Gaussian belongs to the first row.
The main task independently reports `<H>~1.24e-16` and
`<H^2>=17.8916077633` for its 256-grid initial state; those numbers were not
recomputed here. Nonzero spectral width and zero mean exclude a state
supported only at strictly positive energy.

There is an analytic reason to expect symmetric spectral weight: in this
reduced representation `H_mu` is purely imaginary, so `K H_mu K=-H_mu` for
complex conjugation `K`, provided the operator domain is invariant under K.
A real vector then has equal positive/negative spectral weights, apart
from any zero-energy weight. This is a statement about that reduced
operator/domain, not an assumed identity for every numerical discretization.

For the *frozen*, z-independent guide `W(rho,z)=W(rho,0)`, a normalized
transverse zero mode `u(rho)(0,1)` obeys `H=i partial_z`. Thus the carrier
for energy `E0>0` is `exp(-i E0 z)`, not `exp(+i E0 z)`. Multiplying by
`exp[-z^2/(4 sigma^2)]` leaves negative-energy weight

\[
 P_- = \tfrac12\operatorname{erfc}(\sqrt2\,\sigma E_0).
\]

At `sigma=.12`, it is `.5`, `.11506967`, and `.00015910859` for
`E0=0,5,15`. These are exact frozen-guide results; this Fourier projection
is not the spectral projector of the z-dependent dipole operator. A
periodic-box projector would define a different preparation again.

## 2. Exact monopole mode in the project's conventions

Use `hbar=c=1`, electric charge one, `H=+sigma.(p-A)`, `A_t=0`, and
the standard Cartesian curl. About the upper mouth at `z=d/2`, set

\[
 r^2=\rho^2+(z-d/2)^2,\quad \rho=r\sin\theta,
 \quad \tau=\sin^2(\theta/2),\quad A^N_\phi=q\tau,
 \quad \mu=k+\tfrac12,\quad k=0,\ldots,q-1.
\]

Here `q>0` and N is the north regular patch. The monopole field is outward.
Write `B_k=B(k+1,q-k)` and

\[
 \xi_+(\theta)=\binom{\cos(\theta/2)}{\sin(\theta/2)},\qquad
 h_k(\theta)=\sqrt{\frac{\sin\theta}{2B_k}}
       \tau^{k/2}(1-\tau)^{(q-1-k)/2}.
\]

The following is our conversion and normalization in the increment87
cylindrical representation:

\[
 \boxed{f_{E,k}(r,\theta)=r^{-1/2}h_k(\theta)\xi_+(\theta)
                 e^{iEr},\qquad E>0.}
\]

The full spinor is
`Psi=exp(-i sigma_z phi/2) exp(i mu phi) f/sqrt(2 pi rho)`;
time dependence is `exp(-iEt)`. Its sphere spinor has norm one:

\[
 \chi_k=\frac{e^{-i\sigma_z\phi/2}e^{i\mu\phi}}{\sqrt{4\pi B_k}}
 \tau^{k/2}(1-\tau)^{(q-1-k)/2}\xi_+,
 \quad \int d\Omega\,\chi_k^\dagger\chi_k=1,
 \quad \Psi_{E,k}=r^{-1}\chi_k e^{iEr-iEt}.
\]

Direct substitution uses `partial_theta log(h_k)=(mu-q tau)/sin(theta)`
and `sigma_r xi_+=xi_+`; radial and spin-connection terms cancel, giving
`H_mu f_E=E f_E`. The angular probability is precisely
`tau^k(1-tau)^(q-1-k) d tau / B_k`. The spherical `k=0` mode is therefore
normalizable even though its midplane transverse counterpart is marginal.

The underlying angular zero modes and their chiral radial propagation are
published in [Maldacena, Milekhin and Popov, Eq. 2.10 and Appendix C,
C.15–C.18](https://arxiv.org/html/1807.04726v3#A3). The gauge/sign conversion,
flux normalization and direct residual tests here are our calculations;
do not copy the paper's printed magnetic-field sign into this Hamiltonian.

With normal pointing *away from the mouth*, the probability flux is

\[
 \mathcal J_r=\int_0^\pi r\,d\theta\,
 f^\dagger(\sigma_x\sin\theta+\sigma_z\cos\theta)f=+1.
\]

There is no missing `rho` weight: it was absorbed in the reduced spinor.
The factor `r dtheta` is the arc measure of the boundary in the reduced
half-plane. A stationary mode has unit flux, not finite total norm.
At an excised inner sphere the exterior-domain normal is `-e_r`, so this
same mode has domain-outward flux `-1`: it enters the exterior.

## 3. Opposite mouth, patches, and chirality

In a north patch of a negative-charge mouth, use

\[
 A^N_\phi=-q\tau,\quad\mu_N=-(k'+\tfrac12),\quad
 \xi_-=\binom{\sin(\theta/2)}{-\cos(\theta/2)},\quad
 f_{E,k'}=r^{-1/2}h_{k'}\xi_-e^{-iEr}.
\]

It has `H_mu f=E f` and radial flux `-1`. For the given Weyl block the
zero-mode direction is fixed by the charge; one cannot reverse it just
by reversing the positive-energy carrier. The other Weyl block `H=-sigma.Pi`
reverses the physical propagation direction at fixed monopole charge.

Near the lower mouth the project's global potential is
`A_global=q-q tau`. Thus `A_N=A_global-q`, `mu_N=mu_global-q`, and
the same global `mu=k+1/2` corresponds to `k'=q-k-1`. This reverses the
local angular index; setting `k'=k` would match the wrong channel.
Here theta at both mouths is measured from the common positive z axis,
and the azimuthal spin rotation uses the same fixed Cartesian frame.
The lower-mouth angular density is therefore
`tau^(q-1-k)*(1-tau)^k/B(q-k,k+1)`, the reflection of the upper density
under `tau_lower=1-tau_upper`. Returning to the global lower patch multiplies
the spinor by `exp(+i q phi)`. At the lower excision surface `n_ext=-e_r`,
so the mode's radial flux `-1` becomes exterior-outward flux `+1`, flowing
into the lower throat. At the upper surface the signs are reversed.

For either charge sign `s`, the south patch is `A_S=A_N-s q`,
`mu_S=mu_N-s q`. The spinor transition is `Psi_S=exp(-i s q phi) Psi_N`.
The expression `mu-A_phi` and all currents agree on overlap. The
half-angle spin rotation and half-integer mu together make the Cartesian
components single-valued. Use patchwise limits at the axes, not values of
a singular `dphi` coordinate. Relative constant spinor phases are a basis
choice that must be retained when reporting complex scattering amplitudes.

## 4. A strictly positive-energy pulse

Let the chosen reference lead have optical coordinate x increasing away
from its mouth and propagation sign s. Choose any `a>0` and

\[
 a_E=2a^{3/2}E e^{-aE}\,\mathbf1_{E>0},\qquad
 \int_0^\infty |a_E|^2dE=1,
\]
\[
 g(u)=\frac1{\sqrt{2\pi}}\int_0^\infty a_E e^{-iEu}dE
     =\frac{2a^{3/2}}{\sqrt{2\pi}(a+iu)^2}.
\]

Replace `exp(i s E x-iEt)` by `g(t-s x)` in the reference mode. Then
`J_r(t)=s |g(t-s x)|^2`, `integral dt |J_r|=1`, and the packet has norm
one on a complete optical lead. Its norm on a finite collar is only the
fraction present there at that time. Also

\[
 \langle E\rangle=\frac{3}{2a},\quad
 \operatorname{var}E=\frac{3}{4a^2},\quad
 P(E>E_c)=e^{-2aE_c}[1+2aE_c+2(aE_c)^2].
\]

This pulse is a test preparation, with a specified spectrum; `a` must be
chosen in the actual clock units. It has noncompact temporal tails.
Abruptly switching it on, truncating its spatial tail, or applying a time
window destroys exact positive-frequency support and requires a quantified
preparation error. The exact isolated monopole zero mode exists at all E;
the pulse's high-energy tail must still be controlled if an approximate
low-energy collar treatment is used.

For a strictly bounded energy test instead, take
`a_E=sqrt(30/Delta^5)*(E-Emin)*(Emax-E)` on `0<Emin<E<Emax`, zero elsewhere,
where `Delta=Emax-Emin`. Its squared integral is one. The gamma-shaped
pulse above is the one evaluated numerically in this artifact.

## 5. Curved collar and genuinely needed matching data

For a static spherical collar

\[
 ds^2=-N(l)^2dt^2+dl^2+R(l)^2d\Omega^2,
 \qquad dx=dl/N,
\]

the zero-mode spinor in a consistent orthonormal frame is
`Psi=chi_k g(t-s x)/(R sqrt(N))`. Its probability norm is
`integral dx |g|^2`; radial flux per unit Killing time t is
`N R^2 integral dOmega j_hat_r=s |g|^2`.
This follows by inserting the conformal rescaling of [MMP, C.15–C.17 and
Eq. 2.12](https://arxiv.org/html/1807.04726v3#A3). It needs no arbitrary
transmission factor in an exactly spherical, massless, purely magnetic
collar. The collar contributes the optical propagation phase; an assumed
origin extension through a point monopole is not a throat boundary condition.

Some of the spherical data already exist. Read-only inspection of the
finest `cutoff=100` record in `results/dirac_return_closure_results.json` gives
`A=6165936382.795819`, `N_c=6104275059.254434`,
`x_c=1.5708225328305907`, and isotropic radius `98.99996752543319` at
areal radius `R=100` in re units. The exterior-normalized lapse is
`N_c/A=.9899996821709947`. If `t_infinity=A*t_c`, also use
`x_infinity=A*x_c`, `E_infinity=E_c/A`; the phase `E*x` is invariant.
Do not feed the large central-clock lapse into an infinity-normalized
boundary formula. These are conditional saved spherical-model data, not a
completed nonspherical matching solution.

The remaining inputs divide naturally:

* **To define a reference injection now:** select its mouth, global k,
  Weyl block, spectrum and phase/time origin; specify matching sphere,
  outward normal, gauge patch, tetrad rotation and common clock. The
  saved collar can supply N, R and x under its stated approximations.
* **To call it the physical finite-mouth incident channel:** provide or
  solve the nonspherical magnetic/metric collar corrections and their
  overlap with these modes, including excited/evanescent channels where
  required. Nonzero `A_t`, mass or inter-Weyl coupling would change the
  radial problem and must be specified, not silently set by normalization.
* **To extract an exterior scattering matrix:** prescribe incident
  amplitudes only and solve for outgoing response at both mouths and
  infinity. A flux-preserving basis with fixed complex phases is needed;
  arbitrary simultaneous Dirichlet data on both spinor components are
  not an incident-mode boundary condition.
* **To close a throat-plus-exterior quantum state:** additionally supply
  the throat transfer matrix/phase, spin structure and Wilson line and
  state occupations. They are not needed for a conditional single-pass
  exterior S-matrix, but cannot be inferred from its return probability.

## 6. The reference mode's finite-dipole defect is explicit

On an upper sphere of radius b, write `Z=d+b cos(theta)>0` and
`r_other=sqrt(rho^2+Z^2)`. The exact point-pair potential differs from the
isolated monopole by `delta A=q/2*(Z/r_other-1)`. Hence

\[
 (H_{\rm pair}-E)f_E=\sigma_y\delta W f_E,
 \qquad \delta W=\frac{q\rho}{2r_{\rm other}(r_{\rm other}+Z)}.
\]

This rationalized formula is well behaved near the axis. The angular RMS
defect is `sqrt(<delta W^2>_Beta(k+1,q-k))`. At fixed q,k and `b/d -> 0`,

\[
 \epsilon_H\sim\frac{qb}{2d^2}
 \sqrt{\frac{(k+1)(q-k)}{(q+1)(q+2)}}.
\]

Our small test at `q=31,k=4,d=1` gives RMS defects
`.1080115650,.0547078712,.0275313397` for `b=.02,.01,.005`, approaching
the derived leading expression. These have units of energy. They are not
reflection coefficients or certified scattering-error bounds. In particular
`epsilon_H/E` can be large at low E despite a geometrically small collar.
This q=31 demonstration is not a matching certification at the project's
q=100000 or at its saved physical-model collar.

The isolated spherical excited angular gap is `sqrt(q+1)/b`, from
`lambda_n^2=n(n+q)`; [Borokhov, Kapustin and Wu, Section 4 and Appendix](https://arxiv.org/pdf/hep-th/0206054).
Comparing to this local gap diagnoses one source of mixing only. It does
not establish a global gap for the dipole exterior or rescue a planar
zero-mode projection. Here the defect/gap ratios decrease as b squared;
no propagation or accumulated-phase error has been bounded.

## 7. Reproduction and next dependency

```powershell
.\.venv\Scripts\python.exe -B -W error -X utf8 scripts/incoming_modes.py
```

The original reference-mode run exited 0 in 1.94 seconds. Eight controls pass: independent
45-digit differentiation of both charge-sign modes, implementation/formula
agreement, gauge covariance,
direct angular flux quadrature (including q=1 and q=100000), positive
spectrum normalization/moments, integrated pulse flux, independent Fourier
quadrature and the small-collar partner-field asymptotic.

* Maximum relative stationary Dirac residual: `3.38083e-41`.
* Maximum implementation/formula relative discrepancy: `3.84e-11`.
* Maximum tested unit-flux error: `7.77e-11` (large-q double precision).
* Pulse integrated flux: `.9999999999999989`;
  maximum Fourier-transform discrepancy: `1.71e-13`.
* Original reference-only source SHA256: `354d7e818c374250fcc77cf5c0771eb6c8ac78ef7637a2e11f6396cdf9c9ca0b`.

An initial results-write failed on a NumPy boolean; explicit conversion to
Python bool fixed serialization. No mathematical check was relaxed.

**Exact next dependency:** adopt a definite matched collar surface and its
gauge/frame/clock map, then apply this incident mode to the open-boundary
solver while resolving the explicit partner-field defect. Check complex
amplitudes under displacement of that surface with the predicted optical
phase removed. The companion `docs/OPEN_BOUNDARY_SCATTERING.md` task owns the
outgoing-boundary and flux-accounting calculation. It received the formulas.

This completes the incoming *reference-mode* calculation. It supplies no
full scattering solution, escape rate, vacuum stress, stable material,
new physical law, or passage from 13 September 2027 to 12 September 2027.

## 8. Actual saved cutoff100 collar: one bounded magnetic matching term

The requested follow-up evaluates the defect at the existing finest saved
collar, rather than extrapolating the q=31 demonstration. The original
`collar_mode`, `positive_pulse`, and `checks` functions are unchanged;
`saved_collar_checks` reads the existing result and adds six controls.
No closure or historical propagation is rerun.

The selected record is `order=128, refinement=2, cutoff=100` in
`results/dirac_return_closure_results.json`, SHA256
`6baf84864641725f0e2378a368e99eeef5beb0dd08b053a84ea21b674d2fce03`.
It supplies

```text
q = 100000
R = 100                         areal radius, in re
b = 98.99996752543319            saved isotropic radius, in re
d = 4188251800.784477            prescribed separation, in re
b/d = 2.3637539535442946e-8
A = 6165936382.795819            common central-to-exterior clock
N_infinity = N_c/A = .9899996821709947
```

The distance remains the saved model's prescribed coordinate parameter;
it is not a measured separation or a newly solved geometric distance.
Here b, not the areal R, is inserted into the *flat point-pair* potential.
All following energies are in `re^-1`, with `hbar=c=1`. Physical SI energies
would multiply by `hbar*c/re_physical`; no physical re is specified here.

### Exact angular bound and stable expectation

Set `e=b/d`, `c=cos(theta)=1-2 tau`, and

\[
 D=\sqrt{1+2ec+e^2},\qquad K=\frac2{D(D+1+ec)},\qquad
 C=\frac{qb}{2d^2}.
\]

Then `delta W=C sqrt(tau(1-tau))*K`. For every angle and `0<b<d`,
`1-e <= D <= 1+e` and `1-e <= 1+ec <= 1+e`. Consequently

\[
 \frac1{(1+e)^2}\le K\le\frac1{(1-e)^2},\qquad
 \boxed{\sup_\theta |\delta W|\le\frac{qb}{4(d-b)^2}
             =1.4109431438961075\times10^{-13}.}
\]

This is an angular-uniform bound on the multiplication operator at this
sphere, valid for any angular spinor. It is not a bound throughout the
exterior volume. Its factor of length is correct: A_phi and q are
dimensionless, while `delta W=-delta A_phi/rho` has inverse-length units.

Using the exact Beta identity

\[
 \langle\tau(1-\tau)K^2\rangle_{\mathrm{Beta}(k+1,q-k)}
 =M_k\langle K^2\rangle_{\mathrm{Beta}(k+2,q-k+1)},\quad
 M_k=\frac{(k+1)(q-k)}{(q+1)(q+2)},
\]

gives `epsilon_k=C sqrt(M_k)*sqrt(<K^2>_tilted)` without integrating
numbers of order `1e-30`. In particular the exact channel RMS is enclosed by

\[
 \frac{C\sqrt{M_k}}{(1+e)^2}\le\epsilon_k
 \le\frac{C\sqrt{M_k}}{(1-e)^2}.
\]

Thus the leading Beta moment has a rigorous relative enclosure of about
`+/-4.72751e-8`, independent of the numerical integration.
For an independent check, expand

\[
 K^2=1-4ec+e^2(23c^2/2-3/2)+O(e^3).
\]

Under the tilted Beta distribution, with `S=q+3`,
`<c>=(q-2k-1)/S` and
`<c^2>=<c>^2+4(k+2)(q-k+1)/(S^2(S+1))`.
These known moments reproduce the directly integrated `K^2-1` within
`3.31e-22` in the tested channels. The truncated series is an independent
numerical cross-check; the inequalities above supply the rigorous bound.

### Numerical size relative to the gap and saved energies

The saved 128 NS loop samples supply `E_nu=pi/(A L_nu)`, not `2pi/(A L_nu)`:

```text
minimum sampled E = 4.757400615445382e-14
maximum sampled E = 1.3334509095484067e-10
weighted RMS E     = 1.0319386260789103e-10
```

The RMS agrees with `sqrt(w)/(2A)` from the saved state. These are comparison
scales across the stored loop quadrature. They do not assign a single
energy to a finite-width spherical k mode. The smallest stored value is
not a lower bound on the continuous incoming spectrum or the omitted
far-loop endpoint. In particular it cannot justify a uniform low-energy
relative-error bound.

The following table uses the flat residual and flat local spherical gap
`Delta_flat=sqrt(q+1)/b=3.1942368776078194` consistently:

| k | RMS defect epsilon_k | epsilon_k / Delta_flat | epsilon_k / sampled Emin | epsilon_k / sampled ERMS |
|---:|---:|---:|---:|---:|
| 0 | 8.9234532720e-16 | 2.79361e-16 | .0187570 | 8.64727e-6 |
| 16 | 3.6789396881e-15 | 1.15174e-15 | .0773309 | 3.56508e-5 |
| 50000 | 1.4109360225e-13 | 4.41713e-14 | 2.96577 | .00136727 |
| 99999 | 8.9234541157e-16 | 2.79361e-16 | .0187570 | 8.64727e-6 |

The middle channel was chosen as `k=q/2`. Its mirror `k=q/2-1` has the
same leading Beta moment but the opposite odd angular correction. The
polar channels differ slightly because their distance to the partner
differs. Their relative corrections to the leading formula are approximately
`-4.72732e-8` and `+4.72732e-8`; the middle-channel correction is `4.72318e-13`.

This magnetic term is tiny compared with the *local* excited angular gap,
yet need not be tiny relative to every low incident energy. Neither ratio
is a scattering error. Also `xi_+^dagger sigma_y xi_+=0`: the norm of this
residual is not a first-order diagonal energy shift of the zero mode.
No global gap, transmitted amplitude, or accumulated phase follows.

### Optional spherical tetrad conversion, kept separate

If one embeds the same angular potential correction into the assumed
spherical metric at this collar, the azimuthal orthonormal length is
`R sin(theta)`, rather than `b sin(theta)`. Therefore this *magnetic term
alone* would convert as

\[
 \epsilon_{\rm local}=(b/R)\epsilon_{\rm flat},\qquad
 \epsilon_\infty=N_\infty(b/R)\epsilon_{\rm flat},\quad
 N_\infty b/R=.9800993638511766.
\]

The spherical local gap is `sqrt(q+1)/R=3.162293471517152`; its
infinity-normalized value is `3.1306695317333917`. The defect/gap ratio
is unchanged under this consistent conversion. All defect/energy ratios
in the preceding table multiply by `.9800993638511766` if this conditional
embedding is adopted. The conditional uniform infinity-energy bound is
`1.382864477762754e-13`.

This conversion is not a solution of curved nonspherical Maxwell equations
or a calculation of metric, spin-connection, radial-derivative, polarization,
or throat-transfer corrections. Those are separate missing matching terms.
Using the saved isotropic coordinate provides a concrete model evaluation;
it does not certify the map between two complete geometries.

### Verification and handoff

The full small command in section 7 now passes the original eight controls
plus six saved-collar controls, with warnings treated as errors (exit 0,
1.63 seconds). The six check tilted-Beta normalization, the independent
moment expansion, exact per-channel/uniform bounds, RMS energy consistency,
the tetrad/gap conversion, and finite positive scales.

An initial adaptive integration missed `1e-9` probability in a very long
tail interval of the polar Beta distributions. The normalization and
independent-moment checks caught it. Geometrically subdividing those tail
intervals resolved the failure; tolerances were unchanged, and no tail
probability was discarded or renormalized. Final normalization error is
at most `4.67e-15`. Quadrature error estimates are diagnostics, not interval
certificates. The positive pulse still has the algebraic tails described
in section 4; its preparation was not truncated or otherwise changed.

Current source SHA256:
`0e0266ba4e25fac86b96c1a94edc14d68347ec048194f6f12ed3e4972ed1d455`.
Input provenance, individual quadrature diagnostics, all channel bounds,
and both sets of energy ratios are saved in the `saved_collar` JSON field.

**Next dependency:** use this selected R=100/b=98.9999675 collar and common
clock in the outgoing-boundary calculation, with the gauge/frame map from
sections 2–5. The partner magnetic residual now has a concrete bound at the
matching surface. The actual nonspherical collar solution and the exterior
response to this residual remain uncomputed. It is not an S-matrix error
bound, physical injection certification, or evidence of past-directed travel.
