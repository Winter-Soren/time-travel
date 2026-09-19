# Quantum-source completion: keep the fermions on their returning magnetic loops

Research increment, 15 September 2026. Inputs: the authoritative top of
`TASK_STATE.md`, the finite-boundary, Polyakov-source and magnetic-return-path
sections of `RESEARCH.md`, and the saved `quantum_taper` definitions. Those
files and all simulation/results files were read only. No Python jobs were run.

## Recommendation

Use the **adiabatic lowest-Landau-level fermion action on complete magnetic
loops in an axisymmetric throat-plus-exterior**, varying the metric and gauge
field together. This supplies a definite leading angular stress and a
magnetization current. Match its omitted modes and local terms to the original
four-dimensional theory before attempting an exact semiclassical junction.

This is a defensible next source calculation, **not an established global
solution**. A radius profile alone does not change the free zero-mode action
into an angular-pressure source. The positive pressure demanded by the saved
fixed-background taper is not a quantum constitutive law.

Conventions below: signature `(-+++)`, `hbar=c=1`, one massless Dirac field of
electric charge one, magnetic flux magnitude `2*pi*q`, `q>0`; gauge action
`-F^2/(4*g^2)`. The symbol `g` in that expression is the gauge coupling.
`R` is sphere radius, `R_2` is two-dimensional scalar curvature.

## 1. What the primary equations establish

MMP uses Einstein-Maxwell theory with charged massless fermions. There are
`q` angular zero modes. Their normalization is `chi=exp(-sigma/2)*psi*eta/R`,
and their remaining equation is free in conformal coordinates, including
when `R` varies. The exterior field lines satisfy

```text
A_phi = (q/2)(cos(theta_1)-cos(theta_2)),
cos(theta_1)-cos(theta_2) = nu,
L_nu = pi*ell + d*f(nu)                         [leading matched regions].
```

MMP includes both circle Casimir energy and the conformal anomaly; it uses
overlap expansions, not exact electrovacuum joined across a finite quantum
cut. Massive-mode local stresses are omitted at leading order. Field-line
localization fails for the most extended trajectories, around `nu*q~1`.
The antiperiodic spin structure and gauge Wilson line are part of the state
specification. These are published approximations, not a full numerical
four-dimensional stress tensor. [MMP, Eqs. 2.1, 3.13-3.17, 5.24-5.27,
5.44; Appendices C, E, F](https://arxiv.org/pdf/1807.04726)

The newer scattering calculation independently displays the varying-radius
mode action and a gauge-coupled collective mode. It studies probes of an
assumed background; it does not close our stress-tensor matching problem.
[Freivogel, Fumagalli and Tomasevic, v2, 14 September 2026,
Eqs. 4.49-4.54 and 4.78-4.82](https://arxiv.org/pdf/2606.12528)

## 2. Radius dependence: which terms can actually give angular pressure?

For a spherical reduction define the quantum action variation by

```text
delta Gamma = (1/2) integral sqrt(-h) t^ab delta h_ab
              + integral sqrt(-h) J_R delta R,
t_ab = 4*pi*R^2*T_ab,       p_perp = J_R/(8*pi*R).
```

Diffeomorphism invariance then gives, in a purely magnetic spherical state,

```text
D_a t^a_b = 8*pi*R*p_perp*partial_b R.                    (1)
```

**Our inference from the reduction:** for the canonically normalized free
zero modes, `Gamma_0=q*W_D[h;state]` has no explicit `R`. Thus
`p_perp=0` even for nonconstant `R`. The `R^-2` dilution of four-dimensional
radial stress is already included in (1); differentiating that factor alone
does not generate a new pressure. This statement applies to the truncated
action and its chosen regulator, not to the complete renormalized Dirac field.

There are two real ways to change this conclusion within spherical symmetry:

**Higher modes.** In terms of the Hermitian angular Dirac operator,

```text
lambda_0 = 0,                   degeneracy q,
lambda_(n,+/-) = +/-sqrt(n*(n+q)),
degeneracy of each sign = q+2*n, n=1,2,...,
m_n(R) = abs(lambda_n)/R.                                 (2)
```

The spectrum is the monopole-spinor spectrum, with the flux integer renamed
`q` here. [Borokhov, Kapustin and Wu, Section 4 and Appendix,
pp. 9 and 20-21](https://arxiv.org/pdf/hep-th/0206054)
Its use in a radius-dependent determinant is our proposed extension:

```text
p_perp,high = -[1/(8*pi*R^2*sqrt(-h))]
                 sum_modes m_n * delta Gamma_n/delta m_n
               + p_perp,matching.                         (3)
```

The sum counts both nonzero angular eigenvalue signs. Equation (3) is a
functional derivative, including any derivative dependence on `m_n(x)`.
One must determine `Gamma_n` and its state; the equation does not predict a
pressure sign. A large throat gap `m_1*R=sqrt(q+1)` suppresses nonlocal
transmission by these modes, but does not erase their local vacuum
polarization or gauge-coupling renormalization. A derivative expansion fails
where `abs(grad m)/m^2` ceases to be small.

**Matched local terms and interactions.** For illustration, a term

```text
Gamma_C = C * integral sqrt(-h) log(R/R_*) R_2
```

would give, by direct variation,

```text
t_ab,C = 2*C*(D_a D_b log R - h_ab*Box log R),
p_perp,C = C*R_2/(8*pi*R^2),
t^a_a,C = -2*C*Box log R.                                 (4)
```

The commutator of two covariant derivatives verifies (1) for (4).
Thus adding angular stress also changes base stress and anomaly. `C` is
shown to expose that linkage; it is **not a parameter to fit to the taper**.
Dilaton-coupled spinors can produce local curvature-dilaton terms, with a
total-derivative anomaly coefficient requiring a renormalization condition.
Their precise measure and coupling must be matched to this charged reduction.
[van Nieuwenhuizen, Nojiri and Odintsov, Eqs. 6-11](https://arxiv.org/pdf/hep-th/9901119)

Likewise, summing independently renormalized two-dimensional determinants
need not reproduce the four-dimensional subtraction. Sutton explicitly
demonstrates this for scalar spherical reduction; its correction formula is
not a fermion result to import unchanged. [Sutton, Section 1](https://arxiv.org/pdf/hep-th/0003290)

For the spherical U(1) interaction, dimensional reduction and bosonization
give

```text
e_2^2(R) = g^2/(4*pi*R^2),
m_S^2(R) = q*e_2^2/pi = q*g^2/(4*pi^2*R^2).               (5)
```

This is one collective bosonic mode; the other `q-1` combinations remain
massless in this approximation. A radius-dependent determinant of that
mode can supply angular response. It does not gap all `q` channels or
justify extinguishing their vacuum energy. Equation (5) uses our Maxwell
normalization; the collective-mode mechanism is published in
[Freivogel et al., Eqs. 4.73-4.82](https://arxiv.org/pdf/2606.12528).

## 3. Derived leading source in the nonspherical returning field

This section is our static variational extension of the magnetic-loop
approximation. It assumes a smooth positive lapse, negligible mode mixing,
and adiabatic confinement along each complete field line. Let

```text
ds^2 = -N^2 dt^2 + gamma_ij dx^i dx^j,
b = B/abs(B),          n_B = abs(B)/(2*pi),
dN_modes = abs(dPhi)/(2*pi),
dx_nu = dl/N,          L_nu = closed_integral dl/N.
```

`dN_modes` counts flux channels; `N` is the lapse. Each channel is counted
once per complete loop. In axisymmetry the azimuthally distributed quantum
modes are represented by this flux measure, not by an extra factor of `2*pi`.

The candidate low-energy functional is

```text
Gamma_LLL = integral dN_modes W_D[h_nu; global state],
h_nu = -N^2 dt^2 + dl^2 = N^2(-dt^2+dx_nu^2).             (6)
```

For the antiperiodic ground state, its finite Casimir contribution is

```text
E_C = -(pi/6) integral dN_modes / L_nu,
rho_C = -pi*n_B/(6*N^2*L_nu^2),
T_C^(mu nu) = rho_C*(u^mu*u^nu + b^mu*b^nu),
u = N^-1*partial_t.                                       (7)
```

Metric variation is at fixed magnetic flux form. In particular,
`delta L=integral [b^i*b^j*delta gamma_ij/(2*N)
-delta N/N^2] dl`; varying `E_C` gives (7). Integration back gives
`integral N*rho_C*dV=E_C`, checking normalization and redshift factors.

The geometric two-dimensional anomaly must also be retained. With primes
denoting `d/dx_nu`, `sigma=log N`, `Omega_nu=2*pi/L_nu`, the complete
leading per-channel expressions in this regulator are

```text
epsilon_2 = (2*sigma'' - sigma'^2 - Omega_nu^2)/(24*pi*N^2),
p_2       = (             -sigma'^2 - Omega_nu^2)/(24*pi*N^2),
rho=n_B*epsilon_2,       P=n_B*p_2,
T_LLL^(mu nu)=rho*u^mu*u^nu+P*b^mu*b^nu,
-rho+P = n_B*R_2[h_nu]/(24*pi).                           (8)
```

This reproduces the existing Polyakov source when all channels are radial
and share a common circumference. It is a leading flux-sheet description,
not a computation of the full four-dimensional anomaly.

### Angular stress and the current must be derived together

In a local orthonormal mouth frame, (8) gives

```text
T_rtheta = P*b_r*b_theta,  T_thetatheta = P*b_theta^2,
T_phiphi = P*b_phi^2.                                     (9)
```

For poloidal axisymmetric field lines `b_phi=0`. No spherical isotropic
`p_perp` captures this tensor. Moreover `p_2<0` in (8), so bending alone
does not produce the positive isotropic pressure requested by the taper.
For a weak companion field, `b_theta=-epsilon*sin(theta)+O(epsilon^2)`:
mixed stress begins at first order; angular diagonal stress begins at
second order. A linear dipole cannot repair a zeroth-order pressure mismatch.

Put `a=D log N`, `kappa=(b.D)b`, and `a_perp=a-b*(a.b)`. Flux conservation
`D.(n_B*b)=0` and the two-dimensional static Ward equation give

```text
f = spatial_projection(div T_LLL) = P*kappa + rho*a_perp,
j_perp = B cross f / abs(B)^2,
j_perp cross B = f.                                       (10)
```

Thus a prescribed bent field with `j=0` is generally not a conserved combined
source. Equation (10) is independently obtainable from (6), not merely a
force fitted afterward: under a transverse displacement `delta X` of a
loop, its Casimir length changes by

```text
delta L = -closed_integral (kappa+a_perp).delta X * dl/N.
```

In a poloidal axisymmetric patch choose `A=psi*dphi` and
`B=grad psi cross e_phi/varpi`, with `varpi=sqrt(gamma_phiphi)`.
Cross products and gradients use the spatial metric. At fixed flux label,
`delta X=-delta psi*grad psi/abs(grad psi)^2`. Varying
`Gamma_Cas=-integral dt E_C` then gives (10), with a purely azimuthal current;
its divergence vanishes by axisymmetry. Variation of the full induced-metric
functional gives the separate `P` and `rho` weights in (10).
This uses variations supported away from channel endpoints and fixed
loop holonomy. Nonminimal spin/gradient corrections and mode mixing are
outside (6). Maxwell's equation must subsequently be solved with this
current; the original prescribed field need not solve it.

## 4. What counts as a global state and a completed source

The circle state is admissible within the ideal two-dimensional closed-channel
model if `N>0`, every retained `L_nu` is finite, and each phase comes from a
common four-dimensional spin structure and gauge connection. Adiabatic spin
transport must be included in that phase. Independent adjustable phases or
spatial functions `Omega^2(x)` are not this construction.

With channel mixing, use the actual throat and exterior propagation operators,
including their spin and gauge phases:

```text
det[1-U_ext(omega)*U_th(omega)] = 0.                       (11)
```

The diagonal lossless limit gives the antiperiodic frequencies
`omega_k=2*pi*(k+1/2)/L_nu`. Leaky channels require the continuum scattering
problem and its vacuum at infinity; a resonance is not an exact periodic
normal mode. A full four-dimensional Hadamard state, with the correct local
short-distance singularity, cannot be certified from a finite zero-mode
truncation. The full self-adjoint Dirac problem and its ultraviolet matching
remain necessary.

The required final identities are

```text
delta Gamma_q = (1/2) integral sqrt(-g4) T_q^(mu nu) delta g_munu
                + integral sqrt(-g4) j_q^mu delta A_mu,
div j_q=0,                 div T_q = F . j_q,
G_mu_nu = 8*pi*G*(T_EM_mu_nu + T_q_mu_nu),
delta(S_Maxwell+Gamma_q)/delta A_mu=0.                     (12)
```

**Implementation sign guard:** with `+j^mu*delta A_mu` in (12) and
`S_Maxwell=-integral sqrt(-g4)*F^2/(4*g^2)` at constant renormalized `g`,
integration by parts gives

```text
nabla_mu F^(mu nu) = -g^2*j^nu,
equivalently nabla_mu F^(nu mu) = g^2*j^nu.
```

This agrees with `nabla_mu T_q^mu_nu=F_nu_mu*j^mu`. Preserve the index
order and current convention when implementing Maxwell's equation.

The four-dimensional trace must be the anomaly of the specified theory,
with common renormalization conditions for gravitational and Maxwell local
terms. A heat-kernel subtraction fixes that test; setting it equal to the
LLL trace in (8) is insufficient. [Vassilevich, Eqs. 3.27-3.28,
4.28, and 7.8](https://arxiv.org/pdf/hep-th/0306138)

## 5. Checks performed in this increment

* **Varying spherical radius:** the normalization cancels its derivatives
  from the free zero-mode kinetic action; (1) reduces to the existing Ward
  identity with `p_perp=0`. A radius-dilaton term instead obeys all three
  linked variations in (4).
* **Throat normalization:** for `N=R_e/(ell*cos(chi))`, `x=ell*chi` and
  `L=pi*ell`, (8) gives per-channel traceless base stress
  `(t_tt+t_xx)/2=-1/(8*pi*ell^2)`. At `L=2*pi*ell` it vanishes.
  A direct six-point JavaScript check at `ell=3`, `chi=0,.3,1.1` gives
  dimensionless discrepancies below `2.1e-17`.
* **No-flux spectrum:** (2) has no zero modes at `q=0`; its nonzero
  eigenvalues reduce to `+/-n` with degeneracy `2*n` per sign. Reversing
  flux changes zero-mode chirality, not the counting or mass magnitudes.
* **Nonspherical Ward check:** differentiated the Cartesian stress in (7)
  directly using a five-point stencil. Test field: opposite unit monopoles
  at `(0,0,+/-1)`, `N=exp[.03*(x^2+y^2)+.02*z]`,
  `L=5+psi^2`, `psi=(z-1)/r_+-(z+1)/r_-`. This deliberately simple
  loop-constant `L` tests the identity, not physical optical matching.
  Points `(.7,0,.2)`, `(1.3,.4,.6)`, `(2.2,-.3,1.5)` and steps
  `1e-3`, `5e-4` give maximum normalized Ward residuals `5.86e-12`
  and `5.44e-13`; the force/current cross-product residual is below
  `1.58e-12`. Normalization is `max(abs(rho),max_i abs(f_i))` in these
  dimensionless test units. No geometry or Maxwell solution is claimed.
* **Returning-field tail:** in the formal dipole/adiabatic approximation,
  `B=O(q*d/r^3)` and a loop reaching `r` has `L=Omega(r)`, so
  `rho_C=O(q*d/r^5)` and its far energy integral converges. This removes
  the spherical constant-state `1/R^2` tail at that level only. The
  approximation fails near `r~q*d`; it is not a finite-q ADM-mass proof.

**Independent review, reported by the coordinator:** conditional PASS for
the radius/local-term Ward identity, Casimir/redshift normalization,
spherical Polyakov limit, full transverse force, axisymmetric gauge
variation, angular spectrum and collective mode. This validates those
relations within the stated approximation; it does not establish a
completed four-dimensional quantum source or a joined geometry.

## 6. One concrete next calculation

**Evaluate (6)-(10), including the geometric anomaly, on the existing
weak-gravity magnetic-return paths over a compact, adiabatically valid flux
range.** Reuse the saved path geometry and optical lengths. Compute the
angular/mixed stress and the first variation with respect to `A_phi`, then
independently differentiate the resulting tensor to test
`div T_q-F.j_q=0`. Report where the field-line adiabatic condition fails,
and retain the flux-range boundary terms rather than interpreting their
omission as a physical source taper.

Acceptance: recover the spherical Polyakov limit and the existing
Casimir-energy length derivative, and obtain the same transverse force
from gauge variation and direct tensor divergence. Compute the current
from the action's gauge variation independently: defining it from the
tensor divergence alone would not constitute this check. This modest source
calculation can precede a new coupled geometry solve; evaluating it on a
fixed background establishes neither backreaction nor a junction.

**Unresolved:** common four-dimensional renormalization; higher-mode and
spin-gradient response near the mouth; finite-q leakage and global
holonomies; self-consistent Maxwell and metric corrections; one smooth
physical attachment surface; dynamical stability and formation. Controlled
small corrections cannot be assumed to supply the taper's pressure,
enhanced by roughly `1/alpha` in the nearly constant-radius test window.
No full state/geometry completion, material recipe or route from
13 September 2027 to 12 September 2027 has been established.

## 7. Next increment: linear Maxwell response to the verified loop current

Added 15 September 2026; Sections 1-6 above retain the reviewed proposal.
The new `magnetic_loop_source` block now tests its compact Casimir-plus-anomaly
shape variation against tensor divergence. Saved maximum errors are
`9.94e-12` for the Ward identity, `1.45e-8` for action variation and
`1.16e-8` for anomaly variation. The coordinator reports independent reviews
passed within this approximation. I inspected the saved definitions and
current/geometry implementation; I did not rerun its computation.

The profiles use `q=100000`, flux labels `nu=.5,1,1.5`, and the band
`[.5,1.5]`, containing half the flux. Their full reference circumferences
are `L_nu/d=41.21073658045383+f(nu)`. These are flat exterior profiles;
the separate lapse tests hold an unvaried loop complement fixed. Neither
those tests nor the amplified `mu=.01` metric define a completed source.

### 7.1 Operator, normalization and the first response

Restrict to a static, axisymmetric, purely poloidal magnetic sector with no
electric field or frame dragging:

```text
ds^2=-N^2 dt^2+h_AB dy^A dy^B+varpi^2 dphi^2,
A=psi(y)*dphi,       F_Aphi=partial_A psi,
j_hat_phi=varpi*j^phi,
K^AB=N*sqrt(det h)*h^AB/varpi,
J=N*varpi*sqrt(det h)*j^phi=N*sqrt(det h)*j_hat_phi.
```

Direct substitution in the sign convention of (12) gives

```text
L_g psi := partial_A(K^AB*partial_B psi) = -g^2*J.         (13)
```

Here `g` is the constant renormalized gauge coupling; `L_g` denotes the
metric-dependent operator. It is elliptic wherever `N,varpi>0` and `h`
is positive definite. Equation (13) is our specialization of Maxwell's
equations, not the fluid or force-free closure of a plasma model. Its sign
and flat limit agree with [Ioka and Sasaki, Eqs. 2.4-2.5, 5.27,
6.19-6.20](https://arxiv.org/pdf/gr-qc/0302106), replacing their `4*pi`
by our `g^2` and converting coordinate to orthonormal current.

For the saved conformally flat spatial metric, write its conformal factor
as `a` here to distinguish it from the gauge potential:

```text
gamma_ij=a^2*delta_ij,  varpi=a*rho,
partial_rho[N/(a*rho)*psi_rho]+partial_z[N/(a*rho)*psi_z]
    = -g^2*N*a^2*j_hat_phi,

Delta_* psi + grad log(N/a).grad psi = -g^2*a^3*rho*j_hat_phi,
Delta_* = partial_rho^2-rho^-1*partial_rho+partial_z^2.    (14)
```

All gradients in the second line of (14) are in the flat meridional
coordinates. In flat space the orientation is
`B_rho=-psi_z/rho`, `B_z=psi_rho/rho`, so
`(curl B)_phi=-Delta_*psi/rho=g^2*j_hat_phi`.

**Use the saved orientation:** with mouths at `z=+/-d/2`,
`nu=cos(theta_lower)-cos(theta_upper)`, set `psi_0=q*nu/2`.
Writing displacement vectors `R_+=(rho,0,z-d/2)` and
`R_-=(rho,0,z+d/2)`, this gives the implemented field
`B_0=(q/2)*(R_+/abs(R_+)^3-R_-/abs(R_-)^3)` in physical coordinates.
Do not reverse `psi_0` while retaining the saved current sign.

The flat first response on a specified domain is

```text
Delta_* delta psi = -g^2*rho*j_hat_phi,0.                 (15)
```

In the saved dimensionless coordinates `hat y=y/d`, the stored current
is `j_hat_phi,0*d^3`, with **no factor of q**. If the unknown is
`delta u=2*delta psi/q`, (15) becomes

```text
Delta_hat_* delta u=-(2*g^2/q)*hat rho*(j_hat_phi,0*d^3).
```

This normalization exposes the formal `g^2/q` relative field response;
it is not a bound near a boundary, a weak-field zero, or an ill-conditioned
coupled mode. No value of `g` is selected here.
As an algebraic sign check, `psi=rho^4` gives `Delta_*psi=8*rho^2`
and `(curl B)_phi=-8*rho`, hence `j_hat_phi=-8*rho/g^2`.

For a metric perturbation, or a background not solving homogeneous Maxwell,
retain the residual and metric contribution:

```text
L_0 delta psi = -[L_0 psi_0+g^2*J_0]
               -partial_A(delta K^AB*partial_B psi_0)
               -g^2*delta J.                            (16)
```

The saved old flat point field is not a source-free solution on the
nonconstant-lapse test metric. At first weak-field order `N/a=1-2U`, its
additional geometric forcing in the `Delta_*` form is
`+2*grad U.grad psi_0`. The existing `weak_magnetic_response` calculation
already studies this separate geometric correction. Reuse its consistent
order if needed; do not count that forcing twice or identify it with the
new quantum current.

### 7.2 Physical boundary data are still missing

At fixed source, (15) determines a particular response plus a homogeneous
solution. Total magnetic charges and three sampled channels do not determine
that homogeneous solution. A band calculation must evaluate the known
source functional throughout its computational patch, not interpolate
three profiles as though they specify the entire exterior.

| Location | Data or matching required | Merely diagnostic without that data |
|---|---|---|
| Two mouth surfaces | Their actual embeddings, induced metric/lapse, angular flux distribution, and the throat's relation between flux and normal field derivatives; the two mouths and their quantum transport are linked. | `delta psi=0` on fixed cuts clamps the angular flux pattern. Setting a normal derivative to zero assumes a particular tangential magnetic response. Neither follows from fixed total charge. |
| Infinity | Fixed total flux, absence or specification of an applied field, and the asymptotic quantum state/current. For a genuinely source-free far region, select decaying magnetic multipoles. | `delta psi=0` at a finite outer radius; or a dipole-only Robin condition while discarding higher multipoles or currents beyond that radius. |
| Retained flux-band edges | Matching to omitted channels: continuity of the field and its appropriate normal derivative unless an actual surface current is supplied. Edge locations move with the perturbed flux. | Zero extension of the current, reflective edges, or clamping both band edges. These change the problem; no physical completion follows. |
| Symmetry axis, if included | Regular field, patchwise potential constants and fixed monopole-flux jumps. In a regular perturbation gauge, `delta psi=O(rho^2)`. | Allowing new jumps of the axis constants changes magnetic charge or introduces an unmodelled string. |

These conditions follow by integrating `dF=0` and the Maxwell equation
across a physical interface. With one normal pointing from throat to
exterior, no magnetic sheet charge, and tangential surface electric current
`k^a`, they read

```text
[pullback F]=0,       n_mu [F^(mu a)] = -g^2*k^a.         (17)
```

No surface current means continuity of the normal field-tensor contraction,
not that either side vanishes. In three-dimensional language it gives
continuous normal `B` and `n cross [B]=g^2*k`.
Specifying both Dirichlet and Neumann data independently generally
overdetermines the exterior elliptic equation; a throat response map must
supply their relation on the same embedded surfaces. MMP's
[overlapping-region construction, Section 4](https://arxiv.org/pdf/1807.04726)
does not provide these exact finite-cut boundary data.

For a source-free flat region beyond a sphere, a decaying mode has
`delta psi_l proportional to r^-l*(1-mu^2)*P_l'(mu)`, `l>=1`.
Its exact modal outer condition is
`partial_r delta psi_l+(l/r)*delta psi_l=0`; the `l=1` condition alone
is not the whole exterior map. Omitted quantum currents prevent assuming
that source-free region without justification.

The retained edges are labels `psi_a=q/4`, `psi_b=3*q/4`, not material walls.
For their normal displacement `eta`,
`delta psi+eta.grad psi_0=0` at fixed label. Keeping those curves fixed
and setting `delta psi=0` instead imposes a clamped diagnostic boundary.
Elliptic magnetic response can cross these surfaces even though the
idealized fermions stay on their own flux channels.

### 7.3 Vary the full loops before claiming self-consistency

A source held fixed in (15) defines a first response. It can be the first
term in a formal expansion with a single constant `epsilon` multiplying
the quantum action; that bookkeeping parameter is not a spatial taper.
It is not an exact solution of `J=J[psi,N,h;state]` at `epsilon=1`.

Carry each flux label, its spin/gauge holonomy and mode measure with its
loop. For the transverse displacement `eta` and Eulerian metric variations,
the optical-length variation is

```text
delta L_nu = integral dl/N *
                [0.5*b^i*b^j*delta gamma_ij-delta N/N
                 -(kappa+a_perp).eta]
             + [b.eta/N]_segment ends.                   (18)
```

Here `a_perp` retains its acceleration definition in (10).
Sum (18) over the complete loop, including the throat and every mouth
segment. Artificial cut contributions must cancel between neighboring
segments; genuine moving-interface response must be retained. The saved
compact bump killed these terms only for its particular test.
The anomaly action has its own endpoint work, recorded in
`magnetic_loop_source.endpoint_terms`; optical-length endpoint work alone
does not exhaust the boundary variation.
Consequently `delta Omega_nu^2=-2*Omega_nu^2*delta L_nu/L_nu`.
The anomaly term, local magnetic density, tangent and curvature also change.
Moving to a new spatial point and calling the old helper with its old
`f(nu(point))` would generally change channel labels while retaining an
unperturbed optical geometry; that is not this variation.

At fixed metric and with a specified common state/boundary completion,
define the static response kernel by independent gauge variation:

```text
Pi(y,y') = delta J(y)/delta psi(y')
          = (1/(2*pi))*delta^2(Gamma_q/time)
                            /[delta psi(y)*delta psi(y')],
[L_0+g^2*Pi] delta psi = -R_0,
R_0 = L_0 psi_0+g^2*J[psi_0].                            (19)
```

Integration over `y'` in (19) is understood. For a geometry correction add
the `delta K` term and metric derivative of `J` from (16), together with
the Einstein equations. The loop lengths make `Pi` nonlocal. It is the
Hessian of the stated static action under the same boundary conditions;
its positivity, invertibility and dynamical stability have not been shown.
Keeping just `J_0` is the first-response approximation, not evaluation of
this kernel.

### 7.4 One acceptance test for the next coupled calculation

Use a declared finite patch within the retained band, with a stated
diagnostic boundary condition or supplied interface response, plus a
specified variational rule for how the loops join their uncomputed
complement. Dirichlet data for the magnetic field alone do not supply
that optical/state rule. Once those data exist, compute a
Newton direction `v` from (19), then independently reevaluate the gauge
derivative of the full loop action on `psi_0+h*v`, moving the loops and
updating their complete optical lengths and anomaly. With the same metric,
state and boundary rule, test

```text
e = (L_0+g^2*Pi)*v+R_0,
R[psi_0+h*v]-(1-h)*R_0 = h*e+O(h^2).                    (20)
```

Monitor the linear-solve residual `e` and the current-evaluation and
discretization errors separately. When `h*e` and those errors are negligible,
expect quadratic or faster decrease: halving `h` gives approximately a
factor of four while a nonzero quadratic term is resolved, or faster
decrease if that term vanishes. Exact quartering is not required.
Normalize against a fixed initial source/residual scale. Obtain the
reevaluated current from action variation, not from Maxwell's equation or
the stress-divergence expression being checked. Equation (20) tests the
source sign, linear solve and nonlocal feedback together. Passing it is
conditional verification of that patch problem; physical mouth/omitted-band
data, a global state, backreaction and stability remain unresolved.

**Independent review, reported by the coordinator:** conditional PASS for
Section 7's operator, weights/signs, dimensionless current normalization,
interface conditions and Newton/Hessian factor. The reviewer also confirmed
the Ioka-Sasaki primary citation. This reviews the formulation within its
stated scope; it does not establish a solved Maxwell response or physical
boundary/state completion.

No new Maxwell solve or heavy simulation was performed for this addendum.

## 8. What the patch response establishes, and one physical matching step

Addendum, 16 September 2026. I inspected `magnetic_patch_response_checks`
and its saved `magnetic_patch_response` result; I did not rerun them.
The coordinator reports independent checks of the source and continuum
barrier. The interpretation below uses those results; the proposed sourced
mouth calculation has not yet been performed or independently reviewed.

### 8.1 A converged response to the chosen boundary data

The frozen-current problem is
`Delta_* w=-rho*jhat_phi`, with `delta_u=(2*g^2/q)*w` and
`2*g^2/q=7.8956835209e-8`. Coordinates are in units of separation `d`.
Both rectangles impose artificial zero Dirichlet data:

| Rectangle | Grids | Finest centre `delta B_z*d^2/q` |
|---|---|---:|
| `[0.7,1.1] x [-0.2,0.2]` | 33, 65, 129 | `-2.5798041e-14` |
| `[0.6,1.2] x [-0.3,0.3]` | 97, 193 | `-5.8896104e-14` |

The centre is the same point `(rho,z)=(0.9,0)`. Grid restriction changes
decrease by a factor `3.9944`; finest linear residuals are below `7.1e-12`.
At matching grid spacings, the boundary-induced difference changes only
`7.69e-11` under refinement, against a difference of about `1.62e-6` in `w`.
Thus the approximately **2.28-fold centre-field change** is resolved
dependence on the selected boundary problem. It does not estimate the
physical answer by increasing the box.

The absolute perturbation is tiny under these conditions. The maximum
principle bounds the flux correction `w`, not its field derivatives.
Moreover, `w_h=epsilon*rho^2` solves the homogeneous equation exactly and
adds a uniform axial field through changed boundary data. These results
establish a local first response, with no stability conclusion or physical
mouth/state completion. The previously computed weak-metric forcing is a
separate contribution.

### 8.2 One next calculation: the sourced RN mouth transfer relation

**Derive the first axisymmetric coupled dipole transfer relation across a
Reissner-Nordstrom (RN) mouth, retaining unknown boundary amplitudes and
uncomputed source moments explicitly.** Reuse `coupled_exterior` and
`coupled_boundary`. This adds the missing quantum forcing to their existing
homogeneous operator and reconstruction, without assigning arbitrary mouth
values.

The target, in a declared RN normalization, is

```text
(F H_1')'-U_1 H_1 = S_1[T_q0,j_q0],
F=1-2M/r+Q^2/r^2,       Y=(H_1,F H_1')^T,
Y' = [[0,1/F],[U_1,0]] Y + (0,S_1)^T,

Y(r_in) = T(r_in,r_out) Y(r_out)
          + integral_(r_out)^(r_in) T(r_in,s) (0,S_1(s))^T ds,
T(a,b)=Phi(a) Phi(b)^(-1).                              (21)
```

Here `Q` is geometrized RN charge, distinct from the flux integer `q`;
`U_1` is the saved coupled polar-dipole potential, related by duality to
the magnetic sector. Its published **homogeneous** equation is
[Burko, Eq. 87](https://arxiv.org/pdf/gr-qc/9507012). Deriving `S_1` from
the sourced Einstein-Maxwell constraints, including the duality convention,
is the new work. Substituting the flat current into a test-Maxwell RN
equation would omit the coupled metric response and quantum stress.

For homogeneous radial solutions `h_a,h_b`, the columns of `Phi` are
`(h_a,F*h_a')^T` and `(h_b,F*h_b')^T`. The saved exact extremal basis can
be reused only where that approximation is controlled: its `M=Q` background
differs from the slightly `M<Q` optical-matching branch. Use one consistent
background and common clock throughout (21). Its second component is a
coordinate radial momentum; conversion to physical normal data still uses
the embedded surface and clock factors.

Derive the source-dependent terms in the metric/field reconstruction as
part of this calculation. The existing homogeneous reconstruction does not
automatically satisfy matter-sourced equations. Obtain both `T_q0` and
`j_q0` from the same loop action/state and preserve its Ward identity.
Full-sphere dipole projection requires the omitted angular channels:
retain their unknown moments and any band-edge terms explicitly. Extending
the retained-band source by zero would invent a different source.

### 8.3 What a matched region supplies

There are useful, distinct overlaps. Around the upper positive mouth, the
known two-charge exterior supplies, for `r_e << r << d`,

```text
B_0 = q/(2*r^2)*e_r - q/(2*d^2)*e_z + O(q*r/d^3).
```

Here `r` is local distance at leading weak-field order. This fixes the
leading companion incident field and its direction. Reconstructing that
field and its consistent metric perturbation in the master variables can
fix the corresponding incident combination; it does not identify the
saved growing amplitude with a measured magnetic field by assertion.
Decay at infinity about one isolated mouth would discard this companion
term. The quantum field generated in the uncomputed exterior also remains
part of the outer matching data.

At the throat end, the leading overlap has
`1 << rho_th << ell/r_e`,
`r-r_e approximately r_e^2*rho_th/ell` and `t=ell*tau` in MMP's clock
convention. RN connects these two overlaps; the flat region does not
extend through the throat. These leading relations follow from
[MMP, Section 4, Section 5.3 and Appendix B](https://arxiv.org/pdf/1807.04726).
They supply asymptotic geometry and normalization, but do not determine
exact finite-surface stresses or a complete quantum state. Black-hole
horizon regularity cannot replace the traversable-throat response.

### 8.4 Minimum data still needed for a physical number

- **Common background:** compatible mass, flux/radius and clock relations,
  mouth embeddings and their displacement, with companion/inertial terms
  kept to the chosen quasi-static order.
- **Physical response at both ends:** a linked two-mouth relation between
  metric/flux values and normal derivatives, plus the completed outer
  incident/source response. Continuity relates these data; it does not
  set both to zero.
- **Source/state moments:** complete unperturbed optical lengths, spin/gauge
  holonomies and state choice, omitted-channel contributions and the local
  renormalization terms needed for `T_q0,j_q0`. Their variations, including
  `delta L_nu`, are additionally required for nonlocal feedback, beyond
  this first sourced response.

Equation (21) can be derived with these missing quantities as inputs.
Validate it by reconstructing the sourced four-dimensional field equations
and Ward identity; splitting its source integral at an intermediate
matching surface must compose to the same transfer for unchanged physical
data. That checks matching without moving zero-value walls. This dipole
increment leaves the known spherical throat-pressure mismatch and the
other angular sectors unresolved.
