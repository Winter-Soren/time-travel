# Optical heat subtraction and lapse counterterms

Increment121. These are local ultraviolet subtraction data for one full,
four-component massless charged Dirac field on the exact balanced, uncut
classical dihole exterior. They connect the
[lapse susceptibility](EXTERIOR_STRESS_BRIDGE.md) to required counterterm
variations. They do not determine finite vacuum energy or supporting stress.
The [driver](../scripts/dihole_optical_heat.py) and
[saved results](../results/dihole_optical_heat_results.json) evaluate the actual
optical geometry and check the formulas below.

## Operator, normalization and local coefficients

Use the natural domain established in the
[exterior domain proof](DIHOLE_DIRAC_DOMAIN.md). The optical spatial metric is
$g=N^{-2}h_{\rm sp}$, with

$$
g=\alpha^2(d\rho^2+dz^2)+\beta^2d\varphi^2,
\qquad \alpha=P/N,\qquad \beta=\rho/N^2.
$$

All curvatures, contractions and volume elements in this note are optical.
Choose sphere-positive scalar curvature and $\Delta=\nabla^i\nabla_i$,
which is negative on nonconstant flat Fourier modes. The connection is
$\nabla_A=\nabla_S-iA$, with the **charge already included**:
$A=q\,a(\rho,z)d\varphi$, $F=dA$, $q=512000$. Thus

$$
F^2=F_{ij}F^{ij}
=\frac{2q^2(a_\rho^2+a_z^2)}{\alpha^2\beta^2},\qquad
\Delta R=\frac{\partial_\rho(\beta R_\rho)+
\partial_z(\beta R_z)}{\alpha^2\beta}.
$$

There is no additional charge or antiparticle multiplicity. Both chiralities
are included once; these are not coefficients of a single azimuthal sector.
For $P=H_D^2$, define the diagonal trace per optical volume by

$$
k(s;x)=\operatorname{tr}_4 e^{-sP}(x,x)
\sim(4\pi s)^{-3/2}\left[a_0+s a_2+s^2a_4+\cdots\right].
$$

The result is

$$
\boxed{a_0=4,\qquad a_2=-\frac R3,\qquad
a_4=\frac{R^2-3|\mathrm{Ric}|^2-\Delta R}{30}+\frac23F^2.}
$$

To check signs, write $P=-(\nabla_A^2+\mathcal E)$ with
$\mathcal E=-RI/4-i[\gamma^i,\gamma^j]F_{ij}/4$. Then

$$
\operatorname{tr}\mathcal E=-R,\qquad
\operatorname{tr}\mathcal E^2=R^2/4+2F^2,\qquad
\operatorname{tr}\Omega_{ij}\Omega^{ij}
=-|\mathrm{Riem}|^2/2-4F^2.
$$

Substitution into the general Laplace-type coefficients gives
$(5R^2-8|\mathrm{Ric}|^2-7|\mathrm{Riem}|^2-12\Delta R)/360+2F^2/3$.
The three-dimensional identity
$|\mathrm{Riem}|^2=4|\mathrm{Ric}|^2-R^2$ gives the boxed expression.
See [Vassilevich, equations3.27--3.28 and4.26--4.28](https://arxiv.org/pdf/hep-th/0306138).

The local Laplacian term must remain. With compact smearing $w$,
$\int w\Delta R\,dV=\int R\Delta w\,dV$, generally nonzero.
The exact flat uniform-field check is $4sB\coth(sB)$; its $s^2$ coefficient
is $4B^2/3$, since $F^2=2B^2$. Unit $S^2\times\mathbb R$ gives
$a_2=-2/3$, $a_4=-1/15$.

At the saved dyadic central coordinate $\rho\simeq0.6672354660$, $z=0$,
the driver finds $a_2\simeq0.0583718091$ and
$a_4\simeq-0.0238900203$, in units $D^{-2}$ and $D^{-4}$ respectively.
The tiny nonzero $F^2\simeq4.85\times10^{-24}$ reflects evaluation near the
magnetic null. These numbers are subtraction coefficients, not energy signs.
Precision comparisons, independent curvature contractions and exact benchmarks
check the implementation; they are not interval enclosures or heat-remainder
bounds.

## A specified proper-time convention

The parity-even Euclidean fermion effective action has the sign
$-\tfrac12\log\det(-\partial_\tau^2+P)$; see
[de Berredo-Peixoto, equations3--6](https://arxiv.org/pdf/hep-th/0108223).
After the time trace, its ultraviolet energy prescription is

$$
\mathcal E_{\varepsilon,s_0}
=\frac1{4\sqrt\pi}\int_\varepsilon^{s_0}s^{-3/2}
\operatorname{Tr}e^{-sP}\,ds,
$$

interpreted locally, relatively, or variationally when the unperturbed trace
is infinite. A determinant normalization constant is omitted. This is not
the same regulator as a directly cut off occupied-energy sum.
Here $\varepsilon$ has dimensions length squared. The local divergences are

$$
\mathcal E_{\rm div}(x)=\frac1{32\pi^2}
\left[\frac{a_0}{2\varepsilon^2}+\frac{a_2}{\varepsilon}
+a_4\log\frac{s_0}{\varepsilon}\right].
$$

No finite global heat trace is assumed on this noncompact exterior.
The long-proper-time part and finite remainder require separate control.

## Actual lapse variation

Hold the physical spatial metric, magnetic connection and time normalization
fixed. Let $N_\eta=N(1+\eta f)$ with real $f\in C_c^\infty$,
$\|f\|_\infty\le1$, $|\eta|<1$. Then
$g_\eta=(1+\eta f)^{-2}g$. The conformal spinor identification followed by
the volume-unitary map gives, on the fixed Hilbert space,

$$
H_\eta=\sqrt{1+\eta f}\,H_D\sqrt{1+\eta f}
=H_D+\frac\eta2\{f,H_D\}.
$$

The equality is exact by the first-order Leibniz rule. Write
$\delta=\partial_\eta|_0$. The curvature variations are
$\delta dV=-3f\,dV$, $\delta R=2fR+4\Delta f$,
$\delta R_{ij}=\nabla_i\nabla_jf+g_{ij}\Delta f$,
$\delta F^2=4fF^2$, and

$$
\delta(\Delta R)=4f\Delta R+2R\Delta f
+3\nabla f\cdot\nabla R+4\Delta^2f.
$$

Consequently the **local density** variations are

$$
\frac{\delta(dVa_0)}{dV}=-3fa_0,\qquad
\frac{\delta(dVa_2)}{dV}=-fa_2-\frac43\Delta f,
$$

$$
\frac{\delta(dVa_4)}{dV}
=fa_4-\frac{6\mathrm{Ric}:\mathrm{Hess}f
+3\nabla f\cdot\nabla R+4\Delta^2f}{30}.
$$

Define $A_{2j}[f]=\int f a_{2j}\,dV$. Integrating the compactly supported
variation over the whole manifold, contracted Bianchi gives
$6\int\mathrm{Ric}:\mathrm{Hess}f=-3\int\nabla R\cdot\nabla f$;
the Laplacian integrals vanish. Therefore

$$
\boxed{\delta A_{2j}=(2j-3)A_{2j}[f],\qquad j=0,1,2.}
$$

This denotes a finite coefficient variation even when $A_{2j}$ itself diverges.
An additional spatial window would retain boundary/derivative terms.
Where the relative heat variation and cyclicity are justified,
$P'=H_DfH_D+(Pf+fP)/2$ yields the independent check

$$
\delta\operatorname{Tr}e^{-sP}
=2s\partial_sK_f(s),\qquad K_f(s)=\operatorname{Tr}(f e^{-sP}).
$$

This is not the variation of a separately windowed trace
$\operatorname{Tr}(w e^{-sP})$. It also does not assert
$P_\eta=(1+\eta f)^2P$ as operators.
The general heat-coefficient scaling argument appears in
[Vassilevich, equation4.7](https://arxiv.org/pdf/hep-th/0306138).

At **fixed** proper-time cutoffs, integration by parts gives

$$
\boxed{\delta\mathcal E_{\varepsilon,s_0}
=\mathcal E_{\varepsilon,s_0}[f]
+\frac1{2\sqrt\pi}[s^{-1/2}K_f(s)]_\varepsilon^{s_0},}
$$

where $\mathcal E[f]$ means the same integral with insertion $f$.
The cutoff term cannot be discarded. In particular,

$$
\delta\mathcal E_{\rm div}
=\frac1{32\pi^2}\left[-\frac{3A_0[f]}{2\varepsilon^2}
-\frac{A_2[f]}{\varepsilon}+A_4[f]\log\frac{s_0}{\varepsilon}\right].
$$

Changing a metric-dependent cutoff introduces further terms. The driver checks
the local identities by directly perturbing the actual optical metric.

## Remaining physical calculation

Three asymptotic coefficients do not bound the finite heat kernel or its
differentiated coincidence limit. The next requirement is a controlled finite
nonlocal kernel, including all angular modes and proper times, with these
subtractions. Physical renormalized stress additionally requires a consistent
covariant prescription, finite gravitational/gauge counterterms, and the
conformal anomaly; multiplying an optical answer by a lapse power is insufficient.
See [Zahn, sections3--4](https://arxiv.org/pdf/1210.4031) for the charged-field
renormalization framework and [Fulling](https://arxiv.org/pdf/quant-ph/0302117)
for the distinction between heat asymptotics and finite vacuum energy.
No supporting stress, self-consistent throat or chronology change follows here.

Reproduce the saved calculation from the repository root with
`.\.venv\Scripts\python.exe -B scripts/dihole_optical_heat.py`.
