# Regularized receiver: source energy and feasibility

**15 September 2026 — bounded source assessment; independently reviewed within the stated source assumptions.**

**Current follow-up:** Section 7's total initial-energy bound and exact self/mutual-energy formulas have passed independent analytic review. The earlier reviewed radiation derivation is preserved. The coordinator has confirmed completed spectral refinements: the finite-band radiation integral is numerically converged, not interval-certified; its independent source-code and saved-result audit has passed.

## Recommendation

Keep the selected 3.5 ns receiver as a **prescribed Maxwell source**; the completed spectral refinements reported in Section 7 do not establish a physical receiver. For this exact isolated source in vacuum, the earlier analytic derivation below establishes

\[
\boxed{W_{\rm rad}\ge0.3676\ {\rm GJ}.}
\]

This is a conservative, finite-frequency-band lower bound, not an estimate of total radiation. It uses an inequality for the exact spatial Fourier integral; it does **not** replace the short pulse with a long-wavelength multipole formula. The previously established **6.3552 GJ** exterior electrostatic-energy bound remains a separate stored-energy quantity. Adding these two numbers does not establish driver energy.

The coordinator confirms final trajectory-team verification of the selected 3.5 ns discharge: the downstream proton crossing has 75.522401 keV kinetic energy, a 6.814349 m transverse offset and two temporary axial reversals. All trajectory-team checks, including higher-order production and independent verification, passed. Those results are attributed to that team; this report's source audit is independent of the proton trajectory and does not interpret a plane crossing as capture. No material realization or backward calendar-time route follows.

## 1. Source actually assessed

Read-only inputs: the latest receiver entries of `TASK_STATE.md`, the receiver sections of `RESEARCH.md`, and these functions in [future_time.py](<future_time.py>): `_electrostatic_torus_pair_kernel`, `_receiver_discharge_waveform`, `_switched_torus_receiver_kernel`, `_switched_receiver_definition`, `_regularized_receiver_definition`, `_regularized_receiver_source_checks`, and `_receiver_duration_source_checks`.

The source file was changing under its owner's trajectory work. Its observed SHA256 after the source-audit reads was `50696c8bdb9e045a12ccf842cbc1198c7672b890da230197fdcdbc00fbae338b`; the equations and constants below define this assessment independently of later numerical refinements. The coordinator subsequently supplied the final trajectory-team artifacts, whose SHA256 hashes were confirmed by read-only checks here: `future_time.py` — `eee735ed1eef15af31182e25e9ef5ad6324aa4b863a8a43f36ea54a9962dae07`; `future_time_results.json` — `688185af47ae5368369281852d99f06b86c6a970e5a2ef18e7cdad0be4e7b1cd`. These final hashes identify the trajectory team's verified artifacts; they do not imply a trajectory rerun or a new source audit here. No simulation or result file was edited or rerun.

Use receiver-centered cylindrical coordinates. Let

\[
a_i=10\ {\rm m},\quad a_o=20\ {\rm m},\quad b=0.1\ {\rm m},\quad
Q=18.77026815794324\ {\rm C},\quad T=3.5\ {\rm ns}.
\]

The inner and outer surfaces carry opposite charges with normalized measure

\[
d\nu={2\cos^2\alpha\over4\pi^2}\,d\alpha\,d\phi,\qquad
\mathbf x_a=(a+b\cos\alpha)\mathbf e_r(\phi)+b\sin\alpha\,\mathbf e_z.
\]

Thus \(\rho=q(t)\int[\delta^3(\mathbf x-\mathbf x_i)-\delta^3(\mathbf x-\mathbf x_o)]d\nu\), and the actual connecting current is

\[
\mathbf J(\mathbf x,t)=\dot q(t)\mathbf j(\mathbf x),\qquad
\mathbf j=-\int d\nu\,\mathbf e_r\int_{a_i+b\cos\alpha}^{a_o+b\cos\alpha}
\delta^3(\mathbf x-r\mathbf e_r-b\sin\alpha\mathbf e_z)\,dr.
\]

The endpoint terms in its divergence cancel \(\partial_t\rho\) exactly. The nonzero volume current in the otherwise net-neutral transfer region requires mobile carriers if physically realized; continuity alone does not supply them.

For \(s=(t-t_0)/T\in[0,1]\), the implemented waveform is

\[
q(t)=Q\left[1-\int_0^s g(v)\,dv\right],\qquad
g(s)=630s^4(1-s)^4,\qquad \dot q=-Q g(s)/T.
\]

Outside this interval, \(q=Q\) before and \(q=0\) after. The onset translates the radiation pulse and does not affect its integrated energy.

| Quantity | Value for the specified source | Meaning |
|---|---:|---|
| Peak transfer current | 13.19784480 GA | \(630Q/(256T)\) |
| Supremum current density | 1.323981442 GA/m² | Bounded volume current, not a conductor rating |
| \(\int dt\,d^3x\,J^2\) | \(1.078246315\times10^{11}\) A² s/m | Conditional Ohmic diagnostic |
| Initial field energy outside 20.2 m | At least 6.355200925 GJ | Positive electrostatic multipole bound |
| Receiver radiation | At least 0.3676 GJ | New analytic bound below |

The first three values follow independently by substituting \(T=3.5\) ns into the existing exact formulas, including \(\int j^2dV=0.590461304917929\ {\rm m}^{-1}\) and \(\int\dot q^2dt=(4410/2431)Q^2/T\). For a **specified constant scalar resistivity** \(\eta\), \(W_{\rm Ohm}=\eta\int J^2dVdt\). This is conditional arithmetic, not a heating prediction: no constitutive law, temperature, conductor geometry or current response has been solved.

## 2. Exact radiation expression

For time dependence \(e^{-i\omega t}\), the exact outgoing Maxwell far field is the transverse spatial Fourier transform of the current. This is [Wang et al., *Fourier method for identifying electromagnetic sources with multi-frequency far-field data*, Eq. (2.4)](https://arxiv.org/pdf/1801.03263). It retains the factor \(e^{-ik\mathbf n\cdot\mathbf x}\), with \(k=\omega/c\); the far-observer limit is not an electrically small-source approximation.

Define \(\widehat{\dot q}(\omega)=\int e^{i\omega t}\dot q(t)dt\), \(P_{\mathbf n}=I-\mathbf n\mathbf n^T\), and \(\widetilde{\mathbf j}(k\mathbf n)=\int e^{-ik\mathbf n\cdot\mathbf x}\mathbf j(\mathbf x)d^3x\). Parseval's identity and the vacuum Poynting flux give our pulse-energy specialization:

\[
W_{\rm rad}={\mu_0\over16\pi^3c}\int_0^\infty d\omega\,
\omega^2|\widehat{\dot q}(\omega)|^2
\int d\Omega\,|P_{\mathbf n}\widetilde{\mathbf j}(k\mathbf n)|^2.\tag{1}
\]

The coefficient uses positive frequencies, the stated Fourier convention, and real physical fields; no additional sinusoidal time-average factor belongs in (1).

Axisymmetry reduces the spatial integral exactly. Choose \(\mathbf n=(\sin\beta,0,\cos\beta)\). The azimuthal integral gives \(\widetilde{\mathbf j}=iH(k,\beta)\mathbf e_x\), where

\[
H(k,\beta)=\int_0^{2\pi}{2\cos^2\alpha\over2\pi}\,
\cos(kb\sin\alpha\cos\beta)
\int_{r_i}^{r_o} J_1(kr\sin\beta)\,dr\,d\alpha,\qquad
r_a=a+b\cos\alpha.\tag{2}
\]

Here \(J_0,J_1\) are ordinary Bessel functions, unrelated to current density. For \(k\sin\beta\ne0\), the radial integral equals

\[
{J_0(kr_i\sin\beta)-J_0(kr_o\sin\beta)\over k\sin\beta}.
\]

Consequently \(|P_{\mathbf n}\widetilde{\mathbf j}|^2=\cos^2\beta\,H^2\), and

\[
W_{\rm rad}={\mu_0\over8\pi^2c}\int_0^\infty d\omega\,
\omega^2|\widehat{\dot q}|^2
\int_0^\pi \sin\beta\cos^2\beta\,H^2\,d\beta.\tag{3}
\]

All inter-torus and intra-torus interference is inside the squared amplitude. Radiation vanishes along the axis and in the equatorial plane for this exact source symmetry. Equations (1)–(3) include charge radiation through the transverse projection and continuity; adding a separate charge-radiation term would double count it.

## 3. Positive radiation bound without a short-source approximation

Set \(R_s=a_o+b=20.1\) m and restrict the **positive** integral in (1) to

\[
0\le\omega\le\Omega={2c\over R_s},\qquad
{\Omega\over2\pi}=4.747607124\ {\rm MHz}.
\]

The actual pulse has \(cT=1.049273603\) m. We do not assume its spectrum is restricted to this band. Discarding higher frequencies only lowers the energy.

For \(0\le x\le2\), the alternating Bessel series, or the integral identity obtained from \(\sin y\ge y-y^3/6\) for \(y\ge0\), gives the rigorous inequalities

\[
J_1(x)\ge{x\over2}-{x^3\over16}\ge{x\over4}.
\]

Every argument \(kr\sin\beta\) in (2) lies in that interval. Also
\(\cos(kb\sin\alpha\cos\beta)\ge\cos(2b/R_s)>0\). The weighted mean of \(\cos\alpha\) is zero, so with \(D_a=a_o^2-a_i^2=300\ {\rm m}^2\),

\[
H(k,\beta)\ge \cos(2b/R_s)\,{k\sin\beta\over4}
\left\langle{r_o^2-r_i^2\over2}\right\rangle
=\cos(2b/R_s)\,{kD_a\sin\beta\over8}.\tag{4}
\]

The time weight \(g\) is nonnegative and integrates to one. After removing the phase of the pulse midpoint, its transform has real part at least \(\cos(\omega T/2)\). Since \(\Omega T/2=0.05220266682<\pi/2\),

\[
|\widehat{\dot q}(\omega)|\ge Q\cos(\Omega T/2),\quad 0\le\omega\le\Omega.\tag{5}
\]

Use \(\int d\Omega\sin^2\beta\cos^2\beta=8\pi/15\), \(\int_0^\Omega\omega^4d\omega=\Omega^5/5\), and the exact initial quadrupole coefficient \(M_2=QD_a/2=2815.540223691486\ {\rm C\,m^2}\). Equations (1), (4), and (5) imply

\[
W_{\rm rad}\ge
{\mu_0M_2^2\Omega^5\over2400\pi^2c^3}
\cos^2(2b/R_s)\cos^2(\Omega T/2)
\simeq3.676261579\times10^8\ {\rm J}.\tag{6}
\]

This is our derivation for the stated source, not a published apparatus result. The quadrupole coefficient in (6) is a convenient exact geometric moment inside an inequality for the full source transform. No omitted multipole or high-frequency contribution is assigned a sign by an approximation: the omitted **spectral energy** is nonnegative. The numerical value is ordinary floating-point evaluation of the analytic bound; **0.3676 GJ** is conservatively rounded down.

## 4. Stored energy, driver work, and proton work

The initial electrostatic multipole expansion remains legitimate outside a sphere enclosing the static source. It gives \(U_{0,r>20.2}\ge3k_eM_2^2/(10\,20.2^5)=6.355200925\) GJ. This does not impose a long-wavelength approximation on the subsequent radiation.

For the receiver alone, initially static and finally discharged, Poynting's theorem gives

\[
W_{\rm impressed\to EM}\equiv-\int dt\,d^3x\,\mathbf J\cdot\mathbf E_{\rm receiver}
=W_{\rm rad}-U_0.\tag{7}
\]

Here \(U_0\) is the **total** initial receiver field energy, not its exterior lower bound. The sign in (7) permits field-energy recovery. Neither lower bound determines the sign or magnitude of (7). A real driver's electrical input must additionally account for its carrier/material energy changes and losses. The distinction between an imposed current and a solved conducting boundary, and the global nature of the energy balance, are explicit in [Jackson, *How an antenna launches its input power into radiation*, Eq. (1) and Sections I–III](https://arxiv.org/pdf/physics/0506053).

The gamma10 proton begins with about **1.35295 nJ** of kinetic energy. Even the radiation lower bound is about \(2.72\times10^{17}\) times that amount. This compares scales; it is not a per-proton operating-cost or efficiency prediction. Proton work is separately \(\int e\mathbf E\cdot\mathbf v\,dt\). Equation (7) concerns the prescribed receiver self-field and omits proton backreaction and other source cross terms. A driver, shields, return paths or other time-dependent currents may change radiation by interference and require a new combined-source calculation; this isolated-source bound is not universal across receiver designs.

## 5. What can be said about material and mechanical demands

The prescribed surface density is

\[
\sigma_a(\alpha)=\pm{Q\cos^2\alpha\over2\pi^2b(a+b\cos\alpha)}.
\]

At the inner surface maximum, \(|\sigma|=0.9605180743\) C/m². Vacuum Maxwell boundary conditions give \(E_n^+-E_n^-=\sigma/\epsilon_0\), so at least one side has

\[
|E_n|\ge {|\sigma|\over2\epsilon_0}=5.424089109\times10^{10}\ {\rm V/m},
\qquad
u_E\ge{\sigma^2\over8\epsilon_0}=1.302483906\times10^{10}\ {\rm J/m^3}.
\]

These are local field and energy-density bounds. **They do not give a lower bound on the sheet's net mechanical traction**, which depends on both one-sided fields through the Maxwell-stress jump. For example, opposite fields with zero average can have a nonzero jump while their force on the sheet cancels.

If one instead assumes an ideal conductor with zero interior field and the *same local surface density*, its outward electrostatic pressure would be \(\sigma^2/(2\epsilon_0)=52.0994\) GPa at that point. This is a **conditional conductor comparison**, not the load of the current model: the prescribed charge distribution was never solved for equipotential torus surfaces. A physical assessment needs that boundary solution and the actual supporting geometry before applying a material strength or charge-emission model.

The transfer distance is 10 m, with a vacuum light-crossing time of 33.3564 ns; the outer diameter takes 134.0938 ns. A 3.5 ns current change everywhere cannot be attributed to a switch that first communicates its decision across the apparatus during that interval. Distributed prearranged drive could be causal, but no such feed, local energy supply or synchronization response is specified. Nor does the formula prove individual carriers traverse 10 m in 3.5 ns: pre-existing carriers can move locally. The current construction proves conservation and retarded fields, not a carrier or circuit solution.

## 6. Earlier calculation plan and checks — spectrum results now in Section 7

1. **Integrate (3)** for the unchanged 3.5 ns source. Use the analytic radial Bessel integral, one poloidal average, and streamed frequency/observation-angle quadrature. The temporal factor is the single bounded integral \(-Qe^{i\omega t_0}\int_0^1g(s)e^{i\omega Ts}ds\). This avoids a three-dimensional field grid and proton integration. Use a series or the integral form where the Bessel difference suffers cancellation.
2. **Certify the omitted frequency tail**, separately from quadrature convergence. One sufficient bound follows from four integrations by parts: with \(C_4=\int_0^1|g^{(4)}(s)|ds\), \(|\widehat{\dot q}|\le QC_4/(\omega T)^4\). Since \(\int|\mathbf j|dV=a_o-a_i=L\),
   \[
   W_{\omega>\omega_{\max}}\le {\mu_0L^2Q^2C_4^2\over20\pi^2cT^8\omega_{\max}^5}.
   \]
   This is a sufficient, possibly loose tail bound. The first four endpoint terms vanish because \(g,g',g'',g'''=0\) at both ends. No spectral-tail calculation was run here.
3. **Then specify a conductor/drive boundary problem** and compute its required work and stress. Equations (2)–(7) are benchmarks for a proposed realization; matching the existing proton path alone does not satisfy them. Keep full receiver radiation, initial stored energy, material losses and proton work as distinct ledger entries.

Checks performed here used scalar JavaScript arithmetic and small source quadratures only; no extra Python process or large simulation was started. At nine points, \(k/(\Omega/c)=0.001,0.2,1\) and \(\beta=0.25,0.7,1.2\), direct azimuthal integration of the original radial-line source agrees with its Bessel reduction within \(9.8\times10^{-15}\) m. Both used 48 midpoint nodes per integrated variable. All samples satisfy (4). These checks test the reduction and normalization; the inequalities supply the bound between samples. Temporal normalization follows exactly from \(630B(5,5)=1\).

The following dependency-free JavaScript reproduces the reported arithmetic and leaves a small runnable consistency check. It can be pasted into Node.js; it does not import the simulation or write files.

```javascript
const Q=18.77026815794324, c=299792458, eps=8.8541878188e-12;
const T=3.5e-9, b=.1, ai=10, ao=20, Rs=ao+b;
const mu=1/(eps*c*c), omega=2*c/Rs, M2=Q*(ao*ao-ai*ai)/2;
const radiation=mu*M2*M2*omega**5/(2400*Math.PI**2*c**3)
  *Math.cos(2*b/Rs)**2*Math.cos(omega*T/2)**2;
const sigma=Q/(2*Math.PI**2*b*(ai-b));
const peak=630*Q/(256*T);
const j2=.59046130491792895*(4410/2431)*Q*Q/T;
const exterior=3*M2*M2/(40*Math.PI*eps*20.2**5);
function check(ok) { if (!ok) throw Error('Source accounting check failed'); }
check(omega*T/2<Math.PI/2 && 2*b/Rs<Math.PI/2);
check(radiation>3.676e8 && radiation<3.677e8);
check(Math.abs(peak/1.3197844798553839e10-1)<1e-12);
check(Math.abs(j2/1.0782463148506638e11-1)<1e-12);
check(Math.abs(exterior/6.355200924575164e9-1)<1e-12);
console.log({radiationLowerJ: radiation, exteriorStoredLowerJ: exterior,
  peakCurrentA: peak, currentSquaredIntegral: j2,
  minimumOneSidedNormalFieldVm: sigma/(2*eps),
  conditionalConductorPressurePa: sigma*sigma/(2*eps)});
```

**Independent review:** a separate derivation confirmed the Fourier/Parseval coefficient, angular and temporal inequalities, and the numerical lower bound. Nine independent Cartesian source-transform checks agree with the Bessel reduction within `8.9e-15 m`. This verifies the bound for the specified isolated source, not total radiation or driver feasibility.

**Earlier completion boundary:** the analytic radiation bound, exact full-spectrum integrand, energy ledger and conditional mechanical assessment were reviewed before the spectral computation. Section 7 now records the completed spectral refinements. A physical driver's work, conductor equilibrium and physical capture remain uncomputed. The coordinator has integrated these findings into the model review.

## 7. Follow-up: total initial electrostatic energy

**Independently reviewed analytic result:** for the same prescribed weighted surfaces, the **total** initial electrostatic energy satisfies

\[
\boxed{6.3552\ {\rm GJ}\le U_0\le1.0663\ {\rm TJ}.}\tag{8}
\]

The upper bound includes both self-energies and their negative mutual interaction. It does not identify the exterior floor with total energy. It needs no grid, conductor assumption or thin-filament approximation.

The coordinator confirms a **numerically converged finite-band radiation integral** of \(W_{{\rm rad},X<80}=3.2740472480895234\times10^{12}\) J. Frequency refinement changes the integral by **4.6797 J**; doubling the observer resolution and using 64 poloidal nodes changes it by **0.1167 J**. The analytic omitted-tail upper bound is **21,190.217592327 J** at the reported cutoff \(X=80\). These completed refinements are the radiation team's results; its independent source-code and saved-result audit has passed. The finite-band integral is **not interval-certified**, and the tail bound is not a bound on its quadrature error.

The converged finite-band value exceeds the conservative total-energy upper bound of 1.0663 TJ by more than **2.2077 TJ numerically**. Since omitted spectral energy is nonnegative, this supports **positive** work by the impressed current into the electromagnetic field in equation (7). This is a numerical sign conclusion with an analytic bound on \(U_0\), not an interval-certified lower bound on work or a physical driver's energy budget.

### 7.1 Exact self/mutual energy and its units

Let \(w(\alpha)=2\cos^2\alpha\), \(r_a(\alpha)=a+b\cos\alpha\), and \(z(\alpha)=b\sin\alpha\). Averaging the Coulomb kernel over the relative azimuth yields

\[
\begin{aligned}
L_{ab}^2&=(r_a+r_b')^2+(z-z')^2,\\
d_{ab}^2&=(r_a-r_b')^2+(z-z')^2,\\
k_{ab}'&=d_{ab}/L_{ab},\qquad k_{ab}^2=1-(k_{ab}')^2,\\
G_{ab}(\alpha,\alpha')
&={1\over2\pi}\int_0^{2\pi}
{d\varphi\over\sqrt{d_{ab}^2+4r_a r_b'\sin^2(\varphi/2)}}
={2K(k_{ab})\over\pi L_{ab}}.
\end{aligned}\tag{9}
\]

Primes on \(r_b',z'\) mean evaluation at \(\alpha'\), not differentiation; \(k'\) is the complementary elliptic modulus. \(K\) uses the **modulus** convention in [NIST DLMF, Eqs. 19.2.4 and 19.2.8](https://dlmf.nist.gov/19.2). A numerical library using the parameter convention must receive \(m=k^2\). To avoid rounding \(m\) to one near a self diagonal, compute the complement \(d^2/L^2\) directly. An alternative is \(K(k)=\pi/[2\operatorname{AGM}(1,k')]\), from [DLMF Eq. 19.8.5](https://dlmf.nist.gov/19.8#E5), where AGM is the arithmetic-geometric mean.

Define normalized poloidal averages by \(\langle h\rangle=(2\pi)^{-1}\int_0^{2\pi}h\,d\alpha\). Then

\[
S_{ab}=\langle\langle w(\alpha)w(\alpha')G_{ab}(\alpha,\alpha')\rangle\rangle,
\qquad
U_0={k_e Q^2\over2}\,[S_{ii}+S_{oo}-2S_{io}].\tag{10}
\]

The kernels \(G,S\) have units m⁻¹; \(k_eQ^2\) has units J m, so (10) is in joules. The minus sign follows from the opposite torus charges. These are Coulomb energies of the actual prescribed angular charge measures, not capacitances of equipotential conductors.

### 7.2 Removing the integrable self singularity

For one torus, put \(\delta=\alpha-\alpha'\). Its meridional separation is exactly \(d_{aa}=2b|\sin(\delta/2)|\). The [convergent complementary-modulus expansion, DLMF Eq. 19.12.1](https://dlmf.nist.gov/19.12#E1), starts with \(K(k)=\log(4/k')+O(k'^2\log(1/k'))\). Thus the ring-averaged diagonal singularity is logarithmic and integrable; deleting a diagonal sample is not a valid quadrature rule for it.

Use the symmetric subtraction

\[
\begin{aligned}
R_a(\alpha,\alpha')&=G_{aa}(\alpha,\alpha')+
{\log(2|\sin(\delta/2)|)\over\pi\sqrt{r_a r_a'}},\\
R_a(\alpha,\alpha)&={\log(8r_a/b)\over\pi r_a},\qquad
f_a(\alpha)={w(\alpha)\over\sqrt{r_a(\alpha)}}.
\end{aligned}\tag{11}
\]

The coefficient \(2/(\pi L_{aa})\) differs from \(1/(\pi\sqrt{r_a r_a'})\) only by \(O(\delta^2)\). Therefore \(R_a\) extends continuously through the diagonal, with only an \(O(\delta^2\log|\delta|)\) non-smooth remainder. This is an exact subtraction, not an approximation to the energy.

With \(f_{a,n}=\langle f_a e^{-in\alpha}\rangle\), the Fourier series \(-\log(2|\sin(\delta/2)|)=\sum_{n\ge1}\cos(n\delta)/n\) gives

\[
S_{aa}=\langle\langle w w'R_a\rangle\rangle
+{1\over\pi}\sum_{n=1}^{\infty}{|f_{a,n}|^2\over n}.\tag{12}
\]

The Fourier coefficients use the normalized average above. The sum is over positive \(n\) only: adding another factor of two would double count the logarithmic term. The mean logarithm is zero. The smooth positive radius makes \(f_a\) analytic and its Fourier sum rapidly convergent. The mutual term in (10) has no singularity because the surfaces are separated by at least 9.8 m.

**Convergence method for the source owner:** apply periodic quadrature to \(R_a\), supplying its exact diagonal value, and evaluate the logarithmic part by (12) or equivalent logarithmic product quadrature. Refine self and mutual terms separately at \(N,2N,4N\); vary Fourier truncation separately if it is not already negligible. The residual still contains \(\delta^2\log|\delta|\), so do not assume exponential convergence of the whole self integral. For very small \(k'\), evaluate the complementary-modulus series with the logarithm separated to avoid subtractive loss. An interval or explicit remainder analysis is required if a certified numerical interval, beyond bound (8), is desired.

### 7.3 Rigorous upper bound with scalar arithmetic

For two points on one torus of major radius \(a\), write \(A=a-b>0\), and reduce both angle differences to \([-\pi,\pi]\). Their exact distance obeys

\[
\begin{aligned}
D^2&=4b^2\sin^2(\delta_\alpha/2)+4rr'\sin^2(\delta_\phi/2)\\
&\ge {4\over\pi^2}\,[b^2\delta_\alpha^2+A^2\delta_\phi^2].
\end{aligned}\tag{13}
\]

This uses \(rr'\ge A^2\) and \(\sin(|\delta|/2)\ge|\delta|/\pi\). The weighted angular autocorrelation is exactly

\[
\langle w(\alpha)w(\alpha-\delta_\alpha)\rangle
=1+\tfrac12\cos(2\delta_\alpha)\le\tfrac32.
\]

After applying (13), the distance bound no longer depends on the average poloidal angle, so this autocorrelation can be used without replacing the actual surface measure. Integrating the remaining rectangle gives

\[
\begin{aligned}
S_{aa}&\le {3\over4}\int_0^1\!\int_0^1
{du\,dv\over\sqrt{b^2u^2+A^2v^2}}\\
&={3\over4}\left[{\operatorname{asinh}(A/b)\over A}
+{\operatorname{asinh}(b/A)\over b}\right]\equiv B(a).
\end{aligned}\tag{14}
\]

The corner singularity in this rectangle is integrable. Numerically, \(B(10)=0.476384510078266\) m⁻¹ and \(B(20)=0.263308571969715\) m⁻¹. Dropping the negative mutual energy already gives \(U_0\le1.171126806\) TJ. A slightly tighter bound uses the maximum distance between the two surfaces:

\[
|\mathbf x_i-\mathbf x_o'|\le a_i+a_o+2b=30.2\ {\rm m}
\quad\Longrightarrow\quad S_{io}\ge1/(30.2\ {\rm m}).
\]

Hence

\[
U_0\le {k_eQ^2\over2}\left[B(10)+B(20)-{2\over30.2\ {\rm m}}\right]
\simeq1.0662751195\times10^{12}\ {\rm J}<1.0663\ {\rm TJ}.\tag{15}
\]

Equation (15) is an analytic bound for **all** initial field energy. The final 1.0663 TJ is rounded outward. Together with the prior exterior floor it proves (8), independently of the numerical integration below.

### 7.4 Light numerical checks and handoff

Scalar JavaScript evaluation of (9)–(12), using the AGM and 16/32/64 periodic poloidal nodes, produced:

| Poloidal nodes N | Inner self S_ii (m⁻¹) | Outer self S_oo (m⁻¹) | Mutual S_io (m⁻¹) | Total U0 (J) |
|---:|---:|---:|---:|---:|
| 16 | 0.216764971822932 | 0.119411274624214 | 0.0536605427181682 | 362,337,326,220.991 |
| 32 | 0.216764971664213 | 0.119411274604382 | 0.0536605427181682 | 362,337,325,938.299 |
| 64 | 0.216764971644884 | 0.119411274601967 | 0.0536605427181683 | 362,337,325,903.873 |

The finest change is 34.426 J; this is a convergence diagnostic, **not a certified error bar**. Each Fourier sum used positive modes through \(N/2\). Three separate direct azimuthal midpoint integrations with 8,192 scalar nodes agree with (9) within \(4.3\times10^{-15}\) relative: \((a,\alpha;a',\alpha')=(10,0.2;10,2.1),(20,0.8;20,3.7),(10,0.3;20,2.5)\), with radii in metres and angles in radians. All three numerical total energies lie inside analytic interval (8).

These small source checks ran in the existing JavaScript tool runtime without a Python process, large arrays, source edits or result-file writes. The numerical value **about 0.36234 TJ** and its convergence method have been sent to the coordinator for the simulation owner's independent implementation and review. Combining this light check with the numerically converged radiation value suggests about **2.91171 TJ** in (7); this remains a diagnostic subtraction, not a production \(U_0\) result. The sign comparison above uses the analytic \(U_0\) upper bound instead of this light numerical estimate, while retaining the finite-band radiation qualification.

The following standalone JavaScript reproduces the bound without numerical quadrature:

```javascript
{
  const Q=18.77026815794324, eps=8.8541878188e-12, b=.1;
  const ke=1/(4*Math.PI*eps);
  const B=a=>.75*(Math.asinh((a-b)/b)/(a-b)+Math.asinh(b/(a-b))/b);
  const upper=ke*Q*Q/2*(B(10)+B(20)-2/30.2);
  if (!(upper>1.0662e12 && upper<1.0663e12)) throw Error('Energy bound check failed');
  console.log({totalInitialEnergyUpperJ:upper,
    numericalMarginUsingRoundedRadiationJ:3.274047248e12-1.0663e12});
}
```

**Independent Section 7 review:** the coordinator confirms PASS for the saved elliptic-kernel normalization, logarithmic diagonal subtraction, positive-mode Fourier sum without an extra factor of two, total-energy upper bound, units and stated scope. The light \(U_0\approx0.36234\) TJ calculation remains an unfrozen convergence diagnostic. The work conclusion remains conditional on the radiation team's numerical quadrature.

**Updated completion boundary:** the ideal-source total-energy upper bound and integrable exact formulas are independently reviewed; the light convergence check remains unfrozen. The source owner retains all production integration and result-file changes. The radiation team's finite-band integral has completed numerical refinements, with an analytic omitted-tail bound; it is not interval-certified, and the independent source-code and saved-result audit has passed. No conductor, driver efficiency, source preparation or capture realization has been established.
