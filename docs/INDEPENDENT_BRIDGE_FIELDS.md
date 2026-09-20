# Independent interior fields in the sewn physical form

Increment139 removes a restriction in [138](GLOBAL_SEWN_FORM.md): its bridge
field was a linear interpolation of its two exterior traces. The
[new bridge helper](../scripts/dihole_enriched_bridge.py) adds137's regular
monopole profiles times radial sine functions with zero values at both cuts.
Their zero extensions to the exterior are valid global first-order form
trials. They supply independent interior freedom without changing the spin,
magnetic or phase sewing, geometry, state or136 action reference.

The trace lifts and compact functions share radial quadrature and three
angular panels. The polar panels matter: exterior trace functions vanish
there, but the independent interior fields generally do not. Actual
differential actions are evaluated before integration, including mixed
matrix elements between both families. Adding two separately computed
determinants would omit these interactions.

## One matrix assembly gives nested physical trial spaces

Let $L_\kappa=X+\kappa Y$ be137's physical Euclidean residual and
$\delta L_\kappa=\dot X+\kappa\dot Y$ its fixed-Hilbert-space variation.
The [driver](../scripts/dihole_bridge_enrichment.py) integrates these
pointwise products in batches:

$$
Q_\kappa=X^\dagger X+
\kappa(X^\dagger Y+Y^\dagger X)+\kappa^2Y^\dagger Y,
$$
$$
\delta Q_\kappa=
\dot X^\dagger X+X^\dagger\dot X
+\kappa(\dot X^\dagger Y+Y^\dagger\dot X+
\dot Y^\dagger X+X^\dagger\dot Y)
+\kappa^2(\dot Y^\dagger Y+Y^\dagger\dot Y).
$$

The norm matrix $G$ is fixed by the paired metric identification. The
unchanged exterior adds only to the block for the global trace-lift
columns; independent compact columns are zero there. Principal submatrices
select nested trial spaces **before** each space's Cholesky whitening.
Thus one maximal assembly per sector and phase gives several radial,
angular and exterior-polynomial comparisons.

This is a finite restriction of the physical quadratic form. Positivity
or a negative metric response of that matrix does not certify its
uncomputed infinite-dimensional counterpart. The basis still has fixed
compact exterior support and excludes an axis neighborhood.

## Frequency tails with an explicit finite-matrix scope

In the orthonormal finite basis write
$P_\kappa=A\kappa^2+B\kappa+C$ and
$\delta P_\kappa=A_1\kappa^2+B_1\kappa+C_1$. The exact common-mass fibre
response, for dimension $n$, is

$$
F_\kappa=-3M^6\operatorname{Tr}
\frac{\delta P_\kappa}
 {P_\kappa(P_\kappa+M^2)(P_\kappa+2M^2)(P_\kappa+3M^2)}.
$$

The notation uses cyclicity of the trace, not commutation of
$\delta P_\kappa$ with $P_\kappa$. Obtain outward bounds
$a\le\lambda_{\min}(A)$, $b\ge\|B\|$,
$c\ge\max(0,-\lambda_{\min}(C))$ and
$d_2\ge\|A_1\|$, $d_1\ge\|B_1\|$, $d_0\ge\|C_1\|$.
If $\chi=a-b/K-c/K^2>0$, then $P_\kappa\ge\chi\kappa^2I$ for either
sign of $\kappa$ with $|\kappa|\ge K$. It follows that the combined
positive/negative frequency tail, including the per-time measure, obeys

$$
\boxed{\left|\int_{|\kappa|\ge K}\frac{d\kappa}{2\pi}F_\kappa\right|
\le\frac{3M^6n}{\pi\chi^4}
\left(\frac{d_2}{5K^5}+\frac{d_1}{6K^6}+\frac{d_0}{7K^7}\right).}
$$

The [tail helper](../scripts/dihole_finite_frequency_tail.py) reuses the
existing outward binary64 interval primitives. Numerical eigenvectors are
only preconditioners: interval congruence and Gram bounds certify the
needed eigenvalue bounds; Frobenius bounds control the other coefficients.
The helper has independent scalar and dense complex controls.

Its inputs are the **exact binary64 Hermitian canonical matrices** saved
in the [matrix artifact](../data/dihole_bridge_enrichment_state.npz), with
explicit sector/phase keys. The certificate does not enclose the geometric
quadrature, floating whitening relative to the unrounded Galerkin matrix,
the intermediate numerical integral, or physical continuum errors.
Those limitations must not be hidden by calling the entire response
certified. In particular this finite-dimensional frequency tail is not137's
physical spectral or infrared remainder.

The numerical integral now starts at zero. Positive definiteness of each
finite $C$ makes that endpoint regular in this calculation; it is not a
proof of a gap or differentiability in the full noncompact massless theory.
The full two-Weyl finite space is invariant under the chirality swap, so
the middle integral uses $\pi^{-1}\int_0^K F_\kappa\,d\kappa$. The tail
helper already includes its $1/\pi$ and is not multiplied by it again.

## Interpretation

The proper regulator remains $MD=1$. All responses below are finite-matrix
results with that regulator; no covariant local subtraction or regulator
removal has occurred. The required physical source in137 cannot be compared
to them as though these missing limits were controlled. In particular,
136's large-$M$ local terms must not be subtracted at a fixed trial space
and presented as a physical source.

## Executed comparisons and next direction

The three tested sectors are $\mu=-1/2,1/2,3/2$, with maximal
full-component dimensions84,90,84. All use48 exterior-polynomial columns and three radial sines for
three paired angular modes; the central sector also has its upper zero mode.
Seven connected checks pass in42.57 seconds. The bridge helper's six
independent small controls and the tail helper's scalar/dense controls pass.

The polynomial assembly agrees with137's independent single-frequency
implementation to3.59e-16 for $Q$ and6.11e-16 for its derivative. Centered
metric variation agrees to2.33e-10. The old138 global-only frequency responses
are recovered within5.55e-14 relative error. Full coupled higher-quadrature
changes are3.21e-11 for $Q$,3.42e-12 for $G$,8.38e-11 for the derivative;
its sample fibre response changes by5.98e-14. These are numerical checks,
not bounds on the omitted physical fields.

Adding independent interior modes **changes the sign of this finite-model
response**. For the central sector on the old window $\kappa D\in[.01,10]$:

| Exterior columns | Paired angular modes | Radial modes | Phase0 response | Phase $\pi$ response |
| ---: | ---: | ---: | ---: | ---: |
|24|0|0|$+1.90592134\times10^{-6}$|$+1.83193664\times10^{-6}$|
|24|1|1|$-6.82301968\times10^{-4}$|$-6.82237734\times10^{-4}$|
|24|2|2|$-7.16938445\times10^{-4}$|$-7.16927443\times10^{-4}$|
|24|2|3|$-5.60800679\times10^{-4}$|$-5.60766044\times10^{-4}$|
|24|3|2|$-7.16951301\times10^{-4}$|$-7.16940167\times10^{-4}$|
|24|3|3|$-5.60805603\times10^{-4}$|$-5.60770881\times10^{-4}$|
|48|3|3|$-5.60738782\times10^{-4}$|$-5.60856968\times10^{-4}$|

Units are $\hbar c/D$, proper regulator $MD=1$. The radial2-to3 change at
fixed angular2 is21.778% and21.782%. The angular2-to3 change at fixed radial2
is0.00179% and0.00177%. With all three interior radial/angular modes retained,
exterior24-to48 changes the response by0.01192% and0.01535%. These are
particular refinements, not a proof that other angular/exterior errors are
small. They identify radial resolution as the next useful refinement.

The largest central matrices integrate numerically from zero to $\kappa D=80$
to approximately $-0.000614837285578$ and $-0.000614971524909$. Their16-to24
Gauss-order sensitivities are below6.77e-13. The separately certified
finite-matrix $|\kappa D|>80$ tails are below7.452e-8 for each phase. For
the other two tested sectors each corresponding tail is below6.230e-8.
The old window omitted a material low-frequency contribution; the new
finite-matrix integration includes it. This does not establish a physical
infrared limit, and the numerical middle integration is not enclosed.

The canonical arrays stored under `s<sector>_p<phase>_canonical_q` and
`..._canonical_derivative` were independently reloaded. Their outward
coefficient bounds and tail results reproduce the JSON certificates exactly.
Sector order is $(-1/2,1/2,3/2)$ and phase order is $(0,\pi)$. Raw $Q,G,\delta Q$
and weak-Dirac polynomials are also stored, allowing further nested selections
without another geometry assembly. The saved maximal radial order is only3;
adding higher radial modes needs new mixed matrix elements.

See [results](../results/dihole_bridge_enrichment_results.json),
[matrix state](../data/dihole_bridge_enrichment_state.npz), and
[figure](../figures/dihole_bridge_enrichment_progress.png).
Next: increase radial sine order while retaining the same exterior and angular
families, compare the coupled response, and develop its omitted-radial-mode
control. After that, physical exterior/axis and angular completeness, full
infrared behavior, local subtraction and regulator removal still remain.
The present negative finite response is not a new physical law, an obtainable
supporting material, a stability result or a time-travel proof.
