# Open-boundary spinful exterior: flux contract and calibrated local boundary

20 September 2026. Bounded companion to increment87, read together with
`docs/INCOMING_MODES.md` when available. This report derives the boundary accounting
for the actual reduced Weyl operator and supplies a small executable calibration.
It does **not** compute the two-mouth scattering matrix, a physical escape
probability, a vacuum stress, or a stable time-travel medium.

## 1. Operator, domain and exact flux identity

Use the conventions already established in `scripts/spinful_return.py` and increment87:

\[
 H_\mu=-i\sigma_x\partial_\rho-i\sigma_z\partial_z+\sigma_y W,
 \qquad W=(\mu-A_\phi)/\rho,\qquad \mu\in\mathbb Z+\tfrac12.
\]

The reduced spinor has measure \(d\rho\,dz\), because
\(\Psi=e^{-i\sigma_z\phi/2}e^{i\mu\phi}f/\sqrt{2\pi\rho}\).
For an oriented curve bounding the meridional exterior domain, define
\(\Gamma_n=n_\rho\sigma_x+n_z\sigma_z\). Integration over the revolved
physical surface cancels the factor \(2\pi\rho\), so the reduced surface
flux is exactly \(\int f^\dagger\Gamma_n f\,ds\), not an extra
\(\rho\)-weighted integral. These are flat-exterior conventions; a curved
collar needs its own tetrad, conserved inner product and clock conversion.

Direct integration by parts, with real W, gives

\[
 \langle f,Hg\rangle-\langle Hf,g\rangle
 =-i\int_{\partial D}f^\dagger\Gamma_n g\,ds .
\]

For a scalar complex absorbing potential (CAP), \(H_\eta=H-i\eta I\),
\(\eta\ge0\), the local and integrated identities are

\[
 \partial_t(f^\dagger f)+\partial_\rho(f^\dagger\sigma_xf)
 +\partial_z(f^\dagger\sigma_zf)=-2\eta f^\dagger f,
\]
\[
 N_D(T)-N_D(0)+\int_0^T\!dt\int_{\partial D}j_n\,ds
 +2\int_0^T\!dt\int_D\eta |f|^2\,d\rho dz=0.
\]

The W term cancels from this identity but remains in evolution, including the
spin terms \(\sigma_z W_\rho-\sigma_x W_z\) in \(H^2\).
Nothing here substitutes a scalar midplane transmission barrier.

Let \(P_\pm=(I\pm\Gamma_n)/2\). Then
\(j_n=\|P_+f\|^2-\|P_-f\|^2\). Only incoming characteristics can be
prescribed independently in a first-order evolution. The homogeneous local
condition \(P_-f=0\) makes the boundary dissipative, but **does not in general
make it transparent**. In particular, even a pure oblique outgoing plane wave
has a nonzero \(P_-\) component. Characteristic decomposition is an exact
accounting identity, not a separation into physical scattering channels.

## 2. Incoming, reflected, returned, exterior and residual probability

Excise both mouth interiors, leaving ports \(\Sigma_A,\Sigma_B\), and put a
large observation surface \(\Sigma_R\) outside the interaction region.
All normals point outward **from the exterior computational domain**. At
either excised mouth this is opposite to the radial-from-mouth normal.
Thus the companion's upper +q monopole mode with \(\sigma_r=+1\) enters the
exterior (negative exterior-normal current), whereas the lower -q mode with
\(\sigma_r=-1\) leaves the exterior into the receiving mouth.

For an injection at A, distinguish:

| Quantity | Meaning and required measurement |
|---|---|
| Incoming \(I_A\) | Prepared, unit-flux incoming channel packet from A; its energy distribution and phase must be supplied. |
| Reflected \(P_{AA}\) | Flux back into A in its physical outgoing lead modes, including all open channels, not just the injected profile. |
| Returned \(P_{BA}\) | Flux from the exterior into B in receiving lead modes. It is a mouth-to-mouth exterior return, not automatically a completed throat loop. |
| Outer signed flux \(F_R\) | Time-integrated **net** outward current through \(\Sigma_R\), measured before the absorber. At finite radius/time it can include recrossing and reactive/evanescent contributions. |
| Residual \(N_D\) | Norm remaining inside the selected physical observation domain at the stopping time. Include slowly departing packets and any bound or trapped component. |
| Numerical loss \(A_\eta\) | Integral of \(2\eta|f|^2\) in a numerical layer. This is not physical escape. |

In asymptotic leads with resolved, flux-orthogonal incoming/outgoing channels,
and no other incoming data, the physical observation-domain ledger is

\[
 N_D(T)+P_{AA}(T)+P_{BA}(T)+F_R(T)=N_D(0)+I_A(T).
\]

For general finite cuts, use the signed surface-current identity above instead
of forcing each surface contribution into a positive event probability. Counting
every outward crossing as a distinct particle double-counts recrossings.
If a throat feeds amplitude back into the exterior, that is additional incoming
flux; the one-pass ledger cannot silently count it as a new source-free trial.

If the region outside \(\Sigma_R\) is a buffer plus CAP and has no port,
its separate ledger is

\[
 F_R(T)=N_{\rm layer}(T)-N_{\rm layer}(0)+A_\eta(T)+F_{\rm far}(T).
\]

This identity explains why absorption and physical escape differ: they differ
even numerically by stored norm and outer-boundary flux, and the CAP itself can
send artificial reflected amplitude back across \(\Sigma_R\). Do not add both
\(F_R\) and \(A_\eta\) to the same global ledger.

Physical \(P_{\rm escape}\) requires the outgoing scattering flux at infinity
(or a justified asymptotic limit of \(F_R\)), with absence of incoming infinity
data, radius/time convergence, correct long-range propagation and absorber
reflection control across the occupied spectrum. For a normalized packet and
complete asymptotic decomposition, this gives
\(1=P_{AA}+P_{BA}+P_{\rm escape}+P_{\rm bound}\).
Finite-time residual is not proof of a bound state. The nonnormalizable planar
k=0 threshold mode in increment87 supplies none of these probabilities.

## 3. Preserve complex amplitudes at fixed reference surfaces

At a stationary lead, choose modes with a specified gauge, tetrad and phase,
normalized by the indefinite surface form

\[
 \int_\Sigma\chi_a^{s\dagger}\Gamma_n\chi_b^{s'}ds
 =s\,\delta_{ss'}\delta_{ab},\quad s=+1\ ({\rm outgoing}),-1\ ({\rm incoming}).
\]

In a complete propagating-mode expansion, amplitudes are
\(a_b^s=s\int\chi_b^{s\dagger}\Gamma_n f\,ds\). Evanescent fields have
zero individual asymptotic flux but can affect matching; do not discard them
at a finite interaction-region cut or assume this diagonal formula describes
an arbitrary finite-surface basis.

Store the complex matrices \(r_{AA}(E),t_{BA}(E),t_{\infty A}(E)\), as well
as the reference surfaces, gauge transitions, spin frames and common Killing
energy convention. With complete lossless scattering leads, their incoming-A
columns satisfy \(r^\dagger r+t_{BA}^\dagger t_{BA}
+t_{\infty A}^\dagger t_{\infty A}=I\). A CAP destroys exact physical
unitarity and therefore cannot by itself establish that identity.

If \(a(E)\) is normalized in the same flux/energy convention, a stationary
packet return probability is \(\int\|t_{BA}(E)a(E)\|^2dE\). The waveform
and time delay require the **complex** product before Fourier synthesis. A
shift of a reference plane changes propagation phase while leaving probability
unchanged. Compare complex amplitudes only after removing the known phase for
that shift, and transform modes and fields together under gauge/frame changes.

A later vacuum calculation additionally needs the throat scattering/matching
map, the chosen state and spin/Wilson phases, both energy sectors or a justified
analytic continuation, all relevant continuum and bound spectral information,
and renormalization. A single positive-energy packet or scalar \(|t|^2\) is
insufficient. The present implementation deliberately covers positive energies
only; it does not yet supply a vacuum spectral determinant.

## 4. Exact local outgoing map and why a characteristic exit reflects

Freeze W in a planar exterior patch at a +z face and Fourier transform the
tangential rho coordinate with momentum k. This gives an exact constant-
coefficient **spinor** problem

\[
 H_k=-i\sigma_z\partial_x+k\sigma_x+W\sigma_y,
 \quad E>0,\quad p=\sqrt{E^2-k^2-W^2}.
\]

For propagating modes use p>0. The outgoing solution satisfies

\[
 f_2=\beta f_1,\qquad \beta={k+iW\over E+p}.
\]

Below threshold choose Im p>0, giving a decaying evanescent wave rather than
an artificially propagating channel. Threshold unit-flux normalization is
singular and is rejected by the code. A time-domain exact condition uses
the causal frequency-dependent map; its frequency/tangential dependence
means it is generally nonlocal in time and along the face. Freezing W is
only a boundary-symbol calibration for the varying dipole problem, not its
exact exterior map. The \(\mu/\rho\) centrifugal term persists even where
the magnetic field has decayed, and cannot simply be set to zero globally.

For a flat product lead whose W varies along the tangent s but is independent
of normal x, the corresponding exact outgoing map is **operator-valued**.
Define \(Q=-i\partial_s+iW(s)\) on a specified tangential domain, with its
adjoint \(Q^\dagger=-i\partial_s-iW(s)\) and the compatible boundary
conditions. Then

\[
 H=\begin{pmatrix}-i\partial_x&Q^\dagger\\Q&i\partial_x\end{pmatrix},
 \qquad f_2=Q\,[E+\sqrt{E^2-Q^\dagger Q}]^{-1}f_1.
\]

The spectral square root has nonnegative real part on propagating channels
and positive imaginary part on evanescent channels, selected by the outgoing
limit for E>0. The operator order shown matters. In particular
\(Q^\dagger Q=-\partial_s^2+W^2+W'\): replacing it pointwise by
\(k^2+W(s)^2\) loses the derivative term and can miss transverse zero modes.
Both Q and its adjoint sector are needed to characterize the full lead.
Curved boundaries and varying frames introduce geometry/frame corrections;
normal-dependent W requires the actual exterior evolution/outgoing map, not
this product-lead formula applied pointwise. This formula is a derivation and
integration requirement, not an implemented general exterior boundary solver.

The implemented flux-normalized columns are

\[
 u_+={1\over\sqrt{2p(E+p)}}\binom{E+p}{k+iW},\qquad
 u_-={1\over\sqrt{2p(E+p)}}\binom{k-iW}{E+p}.
\]

They obey \(u_s^\dagger\sigma_z u_{s'}=s\delta_{ss'}\). Setting
\(f_2=0\), as the local characteristic condition would, yields

\[
 r=-{k+iW\over E+p},\qquad R=|r|^2={E-p\over E+p}
 ={1-v_n\over1+v_n},\quad v_n=p/E.
\]

Thus \(v_n=0.8,0.2,0.05\) gives spurious reflection
\(1/9,2/3,19/21\). A stable dissipative boundary is not automatically a
good scattering boundary.

## 5. Executable absorber benchmark and independent checks

`scripts/open_boundary_spinor.py` uses the installed NumPy/SciPy stack and imports
only `operator_identity` from `scripts/spinful_return.py`. No periodic evolution,
historical suite or spatial eigensolver is rerun. The standalone command is

```powershell
.\.venv\Scripts\python.exe -B -W error scripts/open_boundary_spinor.py
```

The finite layer has \(\eta(x)=\eta_0(x/L)^3\), 0<x<L, and the deliberately
reflecting termination \(f_2(L)=e^{i\vartheta}f_1(L)\), for which j(L)=0.
Backward integration of the original two-component stationary ODE gives

\[
 f'=[i(E+i\eta)\sigma_z+k\sigma_y-W\sigma_x]f.
\]

At x=0 the code decomposes f into the exact u+/u- basis and divides by its
incoming amplitude. It stores complex r=[real,imag], R and reflection phase.
Adaptive quadrature of \(2\eta f^\dagger f\) independently checks
\(1-R-A_\eta-j(L)=0\). It never defines absorption as 1-R merely to make
the identity pass. This is a stationary **flux ratio**, not a normalized
finite-time packet with a residual probability.

For E=2, k=0.6, W=0.8, eta0=1 and wall phase zero:

| Width L | Artificial reflection R | Complex r at x=0 |
|---:|---:|---:|
| 4 | 0.00664687501085 | 0.0377669866 + 0.0722532334 i |
| 8 | 0.0000721277519853 | -0.00638911832 + 0.00559525863 i |
| 16 | 8.78709716607e-9 | 0.0000368549256 - 0.0000861905542 i |

At the same E, eta0 and L=16, with k=1.17575507654 and W=1.56767343538
(v_n=0.2), **R=0.0243320430643**. This failed-quality grazing example is
retained explicitly. The non-grazing calibration does not validate the CAP
over the physical incoming-mode spectrum.

All 13 checks passed in the saved run:

- Existing symbolic square of the full arbitrary-W operator; random complex
  boundary-projector flux identity error 2.67e-15.
- Spinor eigenvector/unit-flux identities and analytic characteristic-reflection
  formula; exact outgoing-map reflection below 5.2e-31 in the sampled modes.
- Exact evanescent example E=2, k=1, W=2: p=i and beta=(4+3i)/5.
  Symbolic operator residual and normal current are exactly zero, and
  exp(ipx)=exp(-x). The numerical outgoing ratio matches this exact beta.
  `flux_modes` rejects both the threshold E=2,k=0,W=2 and this evanescent
  example, which cannot be normalized as propagating unit-flux channels.
- For k=W=0 the exact wall-backed result is
  \(r=e^{i\vartheta+2iEL-2\int\eta dx}\). Complex-amplitude error 1.63e-12.
- Constant-CAP ODE compared with an independent 2x2 matrix exponential:
  complex-amplitude error 3.99e-13. This also tests the nonzero k and W signs.
- Largest independent loss-balance residual 1.45e-11; tightening ODE tolerance
  changes the L=16 complex r by 1.59e-15.
- Free propagation has t=exp(ipL) and r=0: amplitude/phase error 4.85e-16.

These are numerical checks, not interval-certified global error bounds. The
initial command exposed NumPy-boolean JSON serialization, fixed by explicit
conversion to native bool; numerical equations and thresholds were unchanged.
The code is small backwards shooting for bounded layer strengths; for very
thick/strong layers its unnormalized growing solution can overflow or become
ill-conditioned. A production boundary sweep would use stabilized impedance
or scattering propagation and still perform spectral/angular calibration.

## 6. Regular axis, gauge patches and missing strong-field matching

For the point-dipole gauge in increment87, Aphi approaches q on the axis
between mouths and 0 outside that segment. On each smooth axis segment use
the regular patch \(A_\phi^{reg}=A_\phi-A_{axis}\),
\(\mu_{reg}=\mu-A_{axis}\), so W is invariant. Integer patch changes
transform Aphi and mu together. Away from excised foci, smooth Cartesian
spinors require

\[
 f_+=O(\rho^{|\mu_{reg}-1/2|+1/2}),\quad
 f_-=O(\rho^{|\mu_{reg}+1/2|+1/2}),
\]

with the regular series/parity structure and no logarithmic branch. Their
axis boundary form vanishes. The coordinate axis is not an absorbing exit.
An arbitrary wall or CAP at a small positive rho is not this regularity
condition and would require its own axis-limit convergence. The local planar
benchmark has no axis and makes no claim to implement these global patches.

The isolated-monopole collar modes supplied by the companion define useful
conditional incoming/receiving data. They do not determine the continuation
through a finite, curved, strong-field mouth. Required physical inputs remain:

1. Finite matching surfaces and the actual exterior/collar metric and magnetic
   field, with a controlled overlap between approximations.
2. Matching of tetrads, gauge patches, redshifted energy, and flux-normalized
   spinor modes across those surfaces, including any excited/evanescent modes.
3. A throat transfer/scattering map with its complex phase, state/spin/Wilson
   choice and possible reflection; point-source regularity does not fix it.
4. Asymptotic infinity data and a boundary calibration covering the actual
   packet spectrum, especially near thresholds and grazing/long-tail modes.

Specify incoming mode amplitudes through a consistent incoming boundary map;
do not clamp the full spinor at a mouth while attempting to measure its freely
generated outgoing part. Replacing the throat by a perfectly absorbing port
would be a selected one-pass model, not a derived physical throat matching.

## 7. Primary-source context and integration contract

The flux formulas, outgoing symbol, characteristic-reflection calculation and
benchmark above are derivations for this project's operator. Their scope is
separate from the following published results:

- [Maldacena, Milekhin and Popov, *Traversable wormholes in four dimensions*,
  Section 3 and Appendices C/E](https://arxiv.org/html/1807.04726v3) supplies
  the magnetic-fermion setting. Increment87 already records the curl/sign
  convention adjustment. Its field-line estimate is not a substitute for the
  full spinful boundary problem derived here.
- [Hammer, Poetz and Arnold, *A dispersion and norm preserving finite
  difference scheme with transparent boundary conditions for the Dirac
  equation in (1+1)D*](https://arxiv.org/abs/1302.5587) derives discrete
  transparent conditions for constant-coefficient exteriors by a Z transform.
  Their time convolution and consistency with the interior discretization
  support using an actual outgoing map rather than a local zero-component
  rule. Their discrete scheme has not been implemented in this component.
- [Pinaud, *Absorbing layers for the Dirac equation*](https://www.math.colostate.edu/~pinaud/menu_data/publis/paper_Dirac.pdf)
  constructs Dirac PMLs and compares them with absorbing potentials and
  transport conditions. It distinguishes continuum matching from discretized
  reflection/stability, including fermion-doubling concerns. The present CAP
  is not a PML, and importing a PML into a Fourier split method would need
  separate stability and dispersion checks.

Integration can reuse `boundary_flux` for arbitrary meridional normals and
use `flux_modes`/`modal_amplitudes` only on compatible frozen planar reference
leads. For the actual exterior, combine the incoming task's physical modes
with a compatible flux-preserving discretization, regular-axis patches and
an exterior boundary map or calibrated layer. Record separate port currents,
pre-layer outer flux, layer norm/loss and interior residual, and extract
complex modal amplitudes before absorption. Compare moved boundaries after
phase alignment. No claimed physical return/escape result precedes those steps.

Only this report, `scripts/open_boundary_spinor.py` and its JSON are owned by this task.
`scripts/spinful_return.py`, shared plots, `docs/RESEARCH.md`, `TASK_STATE.md`, and incoming
task files were not edited. The specific proposed passage from 13 September
2027 to 12 September 2027 remains unproved.
