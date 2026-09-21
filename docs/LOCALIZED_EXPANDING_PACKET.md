# A spherical packet: negative linear stress, no supported throat

The [coupled expanding theory](EXPANDING_FIELD_COMPLETION.md) now has a
genuinely three-dimensional, finite-energy perturbation. Its linear source
has a negative core and positive surrounding regions. Twelve linear checks
pass. A separate four-check calculation **rejects the chosen amplitude as a
quantitatively controlled linear prediction of compactness**: geometric
quadratic terms exceed that small signal. Neither result supplies a medium.

## Packet and evolution

Use the previously derived coefficients at $\Lambda/\bar M_{\rm Pl}=0.01$,
$\alpha=1$, $H/\Lambda=0.00408125816086$. Initially $\dot\zeta=0$ and
$\zeta\simeq-\epsilon\exp[-r^2/(2R^2)]$, with $\epsilon=0.001$,
$R=1/k_J=70.70008026/\Lambda$. Specifically,

$$
\zeta(r,t)=-\epsilon\int_0^{10}dq\,w(q)f(q,t)j_0(qr/R),\qquad
\int_0^{10}w(q)dq=1.
$$

The positive weight is proportional to $q^2e^{-q^2/2}$, unchanged through
$q=8$ and smoothly tapered to zero at10. The taper and all its endpoint
derivatives vanish where appropriate. Thus the spatial packet has rapidly
decaying tails, not a hard wall or compact spatial support. Its spectrum is
finite; the largest initial wavenumber is about $0.14144\Lambda$. This alone
does not establish the true effective-theory cutoff.

With $\tau=Ht$ and $r_g=\sqrt A k_J^2/(2H)$, each mode evolves by

$$f''+3f'+4r_g^2q^2e^{-2\tau}(q^2e^{-2\tau}-1)f=0,
\quad f(0)=1,\quad f'(0)=0.$$

The calculation uses the continuous radial transform and ordinary mode
equations, with no reflecting spatial box. Both spectral quadrature and ODE
resolution are refined. The largest relative change among the monitored
fields is $6.5\times10^{-12}$; this is numerical agreement, not an interval
error enclosure. Evolution is sampled through $Ht=4$.

## Total stress, flux and gravity

Let $Q=\dot\zeta-Hn$, $s=\partial_iN^i$, $S_r=\partial_rN^r$ and
$d=\Delta\zeta/a^2$. Direct four-dimensional Ricci/frame expansion gives

$$
B\equiv\delta R_{uu}+\delta R_{ee}
=-2\dot Q+\dot s-\dot S_r+H(s-3S_r)
-\frac{\zeta_{rr}+\Delta\zeta}{a^2}
+\frac{\Delta n-n_{rr}}{a^2},\qquad
\delta R_{ue}=-2Q_r/a.
$$

The full regular spherical expression reduces on shell to
$B=-2(\dot Q+Hs+d)$. The implementation checks both expressions, the
regular-center limits and the angular Einstein equation. The total null
stresses are $T_\pm=P(B\mp4Q_r/a)$, for null directions normalized by
$-u\cdot k_\pm=1$. The energy density and flux are

$$\delta\rho=P(-2d+6HQ-2Hs),\qquad j=T_{ue}=-2P Q_r/a.$$

These include the full linear source of the action; they are not inferred
from the two-derivative $\mathcal P_X$ contribution alone. Initially,

| Central quantity | Value in units of $\Lambda^4$ |
| --- | ---: |
| Both null stresses | $-1.6799273741\times10^{-5}$ |
| Energy density | $-1.7985965512\times10^{-5}$ |
| Null stress at $Ht=4$ | $-1.3205903481\times10^{-9}$ |

The initial density has a positive surrounding shell; the signed radial
integral cancels to $4.8\times10^{-14}$ of the absolute integral numerically.
For the smooth full-space linear packet, the zero spatial mode of the
Laplacian and velocity vanishes, explaining the compensation analytically.

For gravitational response, use spherical areal radius
$\mathcal R=ae^\zeta r$ and geometric Misner--Sharp compactness
$C=1-g^{ab}\partial_a\mathcal R\partial_b\mathcal R$.
The historical spherical-mass framework is due to
[Misner and Sharp](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B571);
the following expanding-background subtraction and flux checks are derived
here for this perturbation, not taken from their ideal-fluid equations.
Subtract the de Sitter value **at the same areal radius**:

$$\delta C_{\mathcal R}=C-H^2\mathcal R^2
=2H(ar)^2(Q-N^r/r)-2r\zeta_r+O(\epsilon^2).$$

The associated mass integral contains $\delta\rho-Har\,j$, not density
alone. This flux-corrected integral agrees with the boundary compactness to
$5.9\times10^{-13}$ of the absolute integral. The largest sampled linear
$|\delta C_{\mathcal R}|$ is $2.670\times10^{-6}$, subject to the decisive
truncation problem below.

## No throat in this packet ansatz

A spatial neck requires $\partial_r\mathcal R=ae^\zeta(1+r\zeta_r)=0$.
Initially the Gaussian has $r\zeta_r=\epsilon(r/R)^2e^{-r^2/(2R^2)}>0$.
The smooth spectral approximation obeys a stronger all-time bound:

$$|r\zeta_r|\le2\epsilon M<0.0284,\qquad
M=e^{0.025^2/3}[1+5(\log10+1/3)]<14.2.$$

To see this, $q\le1$ is controlled directly by the positive Volterra kernel.
For $1<q<10$, mode energy decreases before band entry, giving
$|f'|\le2r_gq\sqrt{q^2-1}\le5$ and
$|f|\le1+5\log10$ at entry. The general-entry Volterra estimate then gives
$|f|\le M$. Finally $|y j_0'(y)|\le2$ and normalized positive weights give
the derivative bound. It covers the continuum linear packet, not just the
sampled radial grid. The sampled maximum is much smaller, about0.000736.

Thus this linear packet ansatz has no spatial neck. Its physical width also
expands by $e^4=54.60$ over the simulated interval; the nominal width equals
the background cosmological horizon at $Ht\simeq1.243$. A cosmological
marginal sphere is not a wormhole mouth or a stationary apparatus.

## The geometric truncation test fails at the chosen amplitude

The small metric fields do not by themselves validate a tiny gravitational
signal. Evaluate the exact compactness formula on the reconstructed fields
$N=1+n$, $\gamma_{ij}=a^2e^{2\zeta}\delta_{ij}$ and $N^r$:

$$
F=\frac{\dot\zeta-(N^r/r)(1+r\zeta_r)-Hn}{1+n},\qquad
C-H^2\mathcal R^2=\mathcal R^2(2HF+F^2)-2r\zeta_r-(r\zeta_r)^2.
$$

This evaluates geometry exactly for those reconstructed fields. It is **not**
a nonlinear solution of the matter/gravity constraints or evolution.

| Amplitude | Maximum geometric remainder / maximum linear signal |
| --- | ---: |
| 0.001 | 1.2449 |
| 0.0001 | 0.12447 |
| 0.00001 | 0.012447 |
| 0.000001 | 0.0012447 |

At amplitude0.001, the remainder is $3.324\times10^{-6}$, exceeding the
linear maximum. The largest discrepancy already occurs at $Ht=0$, $r/R=1.05$,
inside the background horizon. The remainder decreases quadratically with
amplitude as expected. This does not bound the missing second-order physical
fields. Lower amplitude weakens the stress and gravitational signal too;
it does not establish stronger support or nonlinear stability.

**Decision:** retain the linear localized stress result and no-neck bound,
but reject amplitude0.001 as a quantitatively controlled linear compactness
prediction. The next step is constraint-compatible nonlinear initial data.
No new law, available material, stable bridge or backward passage is proved.

## Nonlinear constraints prepared for the next calculation

For specified $\zeta(r)$ and $\dot\zeta(r)$, put $\beta=N^r$,
$A_0=H+\dot\zeta-\beta\zeta_r$, $E_r=A_0-\beta_r$,
$E_t=A_0-\beta/r$, $E=E_r+2E_t$, $Q_E=E_r^2+2E_t^2-E^2$,
$X=X_0/N^2$ and $c_0=\alpha\Lambda^2X_0$. The full radial action gives

$$
\frac P2\left(R^{(3)}-\frac{Q_E}{N^2}\right)
+\Lambda^4\left[\frac{(X-1)^2}{4}-X(X-1)\right]-\rho_\Lambda
+\frac{3c_0E^2}{N^4}=0,
$$

$$
\partial_r\left(\frac{PE_t}{N}+\frac{c_0E}{N^3}\right)
-\frac P N(\zeta_r+1/r)(E_r-E_t)=0,
$$

where $R^{(3)}=a^{-2}e^{-2\zeta}[-4\Delta\zeta-2\zeta_r^2]$.
The lapse constraint is quadratic in $y=N^{-2}$:

$$
(3c_0E^2-3\Lambda^4X_0^2/4)y^2
+(-PQ_E/2+\Lambda^4X_0/2)y+PR^{(3)}/2+\Lambda^4/4-\rho_\Lambda=0.
$$

The branch must connect continuously to $N=1$ at infinity. Smooth-center
data have $\beta=b_1r+O(r^3)$ and even $N,\zeta$; regularity does not set
$\beta_r(0)=0$. Independent derivation checks the homogeneous background and
recovers the previous linear constraints. A nonlinear solver and evolution
have not yet been executed; these equations define the next physical test.

Reproduction: [packet calculation](../scripts/localized_expanding_packet.py),
[packet results](../results/localized_expanding_packet_results.json),
[saved fields](../data/localized_expanding_packet_state.npz),
[geometry test](../scripts/localized_packet_geometry.py),
[geometry results](../results/localized_packet_geometry_results.json), and
[figure](../figures/localized_expanding_packet_progress.png).
