"""A global heat-operator error for collar-supported data, not a source bound.

Convert187's columns to rows, apply Schur, then use the already specified
global nonnegative operator in Duhamel's identity. No exterior is truncated.
"""
import json

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, upper
from dihole_heat_defect import unpack
from dihole_local_transport import packed
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_localized_heat_error'
DEGREES = dict(leading=3, constant=0, varying_field=3, quadratic=2, quartic=4)
POWERS = dict(leading=1, constant=2, varying_field=3, quadratic=4, quartic=6)
CUT_DEGREES = dict(gradient=1, constant=0, quadratic=2)
CUT_POWERS = dict(gradient=1, constant=2, quadratic=4)


def enclosure(iv, previous, record, fraction):
    if fraction.a<=0 or fraction.b>1: raise ValueError('Require0<fraction<=1.')
    eps = unpack(iv, previous['constants']['epsilon'])
    if eps.b>=1: raise ArithmeticError('Metric comparison requires epsilon<1.')
    ratio = (1+eps)/(1-eps)
    R = [ratio**(iv.mpf(2)+iv.mpf(k)/2) for k in range(5)]
    # Mixed-degree187 groups use their larger degree, since R_k increases.
    interior = {k:unpack(iv,v)*fraction**(iv.mpf(POWERS[k])/2)
                for k,v in record['interior_terms'].items()}
    # For fraction<1 use the original geometric beta lower endpoint instead:
    # the packed beta is an UPPER bound and cannot certify exponential decay.
    geom = json.loads(project_path('dihole_local_transport_results.json').read_text())
    heat = json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    s0 = iv.mpf(next(r['cutoff_s0_over_D_squared'] for r in heat['records'] if r['q']==record['q']))
    H = iv.mpf(geom['real_geometry']['inverse_metric_norm_upper'])
    exact_beta = iv.mpf(2)**-64/(32*H*s0)
    decay = iv.exp(-exact_beta*(1/fraction-1))
    cut = {k:unpack(iv,v)*fraction**(iv.mpf(CUT_POWERS[k])/2)*decay
           for k,v in record['cutoff_terms'].items()}
    col = sum(interior.values())+sum(cut.values())
    row_parts = {k:v*R[DEGREES[k]] for k,v in interior.items()}
    row_cut = {k:v*R[CUT_DEGREES[k]] for k,v in cut.items()}
    row = sum(row_parts.values())+sum(row_cut.values())
    eta = iv.sqrt(col*row)
    chart_Q = iv.exp(iv.mpf(1)/100)*iv.sqrt(R[0])
    # Orthogonal hemisphere source projections; defect outputs may overlap.
    # Collar output supports are disjoint, giving a direct sum across collars.
    global_error = iv.sqrt(2)*eta
    triangle = 1+iv.sqrt(2)*chart_Q
    useful = bool(global_error.b < triangle.a)
    selected = global_error if useful else triangle
    return dict(q=record['q'], time_fraction=float(fraction.a),
        row_moment_factors=[packed(iv,v) for v in R],
        interior_row_terms={k:packed(iv,v) for k,v in row_parts.items()},
        cutoff_row_terms={k:packed(iv,v) for k,v in row_cut.items()},
        column_integral_upper=packed(iv,col), row_integral_upper=packed(iv,row),
        chart_integrated_L2_defect_upper=packed(iv,eta),
        chart_Q_operator_upper=packed(iv,chart_Q),
        localized_global_heat_error_upper=packed(iv,global_error),
        trivial_triangle_error_upper=packed(iv,triangle),
        selected_heat_error_upper=packed(iv,selected), improves_triangle=useful,
        differentiated_diagonal_error=None, physical_source=None, atom_transit=None)


def controls():
    """Varying-source Gaussian integrals and an exact exterior-coupled control."""
    mp.dps = 80
    ratios = []
    for eps in map(mp.mpf, ('.0033', '.2')):
        for k in range(5):
            R = ((1+eps)/(1-eps))**(2+mp.mpf(k)/2)
            # g(y)=scale(|y|)*I. Its variation remains inside the comparison.
            scale = lambda r: 1+eps*mp.sin(r)
            mass = 2/mp.gamma(mp.mpf(k)/2+2)*mp.quad(
                lambda r: scale(r)**(2+mp.mpf(k)/2)*r**(k+3)*mp.exp(-scale(r)*r*r),
                [0, 1, 3, 8, mp.inf])
            ratios.append(float(mass/R))
    P = mp.matrix([[3,mp.mpf('.3'),mp.mpf('.4')],
                   [mp.mpf('.3'),4,mp.mpf('.2')],[mp.mpf('.4'),mp.mpf('.2'),5]])
    A = mp.diag([2,3]); J = mp.matrix([[1,0],[0,1],[0,0]])
    B = P*J-J*A
    eigenvalues, U = mp.eigsy(P)
    UB = U.T*B
    errors, leakage, bound_ratios = [], [], []
    for t in map(mp.mpf, ('.01','.1','.3')):
        H = U*mp.diag([mp.exp(-t*v) for v in eigenvalues])*U.T
        Q = J*mp.expm(-t*A)
        integral = mp.matrix(3,2)
        for i in range(3):
            for j in range(2):
                integral[i,j] = mp.fsum(U[i,l]*UB[l,j]*mp.exp(-eigenvalues[l]*t)
                    *mp.expm1((eigenvalues[l]-A[j,j])*t)/(eigenvalues[l]-A[j,j])
                    for l in range(3))
        errors.append(float(mp.norm(H*J-Q+integral)))
        leakage.append(float(mp.norm((H*J)[2,:])))
        # Frobenius norm bounds the operator norm; both semigroups contract.
        bound_ratios.append(float(mp.norm(H*J-Q)/(t*mp.norm(B))))
    # Overlapping outputs show sqrt2 is necessary for orthogonal inputs.
    overlap = mp.matrix([[1,1],[0,0]])
    split_norm = mp.sqrt(max(mp.eigsy(overlap.T*overlap,eigvals_only=True)))
    return dict(varying_source_gaussian_max_ratio=max(ratios),
        global_P_min_eigenvalue=float(min(eigenvalues)),
        exterior_coupled_Duhamel_error=max(errors),
        exterior_leakage_min=min(leakage), Duhamel_bound_max_ratio=max(bound_ratios),
        orthogonal_source_factor_error=float(abs(split_norm-mp.sqrt(2))),
        omitted_exterior_would_fail=min(leakage)>0)


def run():
    iv = MPIntervalContext(); iv.dps = 100
    path = project_path('dihole_frozen_magnetic_results.json')
    previous = json.loads(path.read_text())
    inputs = {**previous['input_sha256'], path.name:sha(path),
        'dihole_frozen_magnetic.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing187 inputs required.')
    records = [enclosure(iv,previous,r,iv.mpf(1)/n) for r in previous['records'] for n in (1,2,4,16)]
    control = controls()
    checks = dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        row_not_below_column=all(unpack(iv,r['row_integral_upper']).a>=unpack(iv,r['column_integral_upper']).b for r in records),
        varying_source_gaussian=control['varying_source_gaussian_max_ratio']<=1,
        nonnegative_global_control=control['global_P_min_eigenvalue']>0,
        global_Duhamel_identity=control['exterior_coupled_Duhamel_error']<1e-65,
        global_Duhamel_bound=control['Duhamel_bound_max_ratio']<1,
        exterior_retained=control['omitted_exterior_would_fail'],
        orthogonal_hemisphere_factor=control['orthogonal_source_factor_error']<1e-65,
        atomic_bounds_useful=all(r['improves_triangle'] for r in records[4:]),
        low_flux_loose_bound_retained=all(not r['improves_triangle'] for r in records[:4]),
        smaller_times_improve=all(unpack(iv,records[i+1]['localized_global_heat_error_upper']).b<
            unpack(iv,records[i]['localized_global_heat_error_upper']).a for i in (0,1,2,4,5,6)),
        no_source_or_transit_claimed=all(r['physical_source'] is None and r['atom_transit'] is None
            and r['differentiated_diagonal_error'] is None for r in records))
    result = dict(passed=all(checks.values()), checks=checks, controls=control, records=records,
        source_sha256=sha(__file__), input_sha256=inputs, bound_encoding=previous['bound_encoding'],
        scope=['Global heat operator H=e^-tP of the already specified nonnegative physical Euclidean operator; no exterior truncation.',
            'Initial-data projection Pi is the two collar cores, each split into measurable orthogonal hemisphere source projections.',
            'Qloc uses original187 local kernels, source projections only, and extension of compact output supports into the global physical Hilbert space.',
            '||H(t)Pi-Qloc(t)||_(L2->L2) is bounded uniformly0<t<=T by Duhamel, integrated column/row Schur bounds and sqrt2 for hemispheres.',
            'This is not pointwise, differentiated-diagonal, trace-norm or supporting-source control, and is not atom propagation.',
            'One hydrogen atom, trial1nm opening. Support, stability, survival and clock gain remain unproved.'],
        sources=['https://web.math.ku.dk/grubb/notes/heat.pdf'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']: raise ArithmeticError('Localized heat checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,atomic=[dict(
        fraction=r['time_fraction'],error_log10=r['localized_global_heat_error_upper']['log10_upper'])
        for r in records[4:]]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = result['records'][4:]
    values = [10**r['localized_global_heat_error_upper']['log10_upper'] for r in rows]
    fig,ax = plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['T = s0','T = s0/2','T = s0/4','T = s0/16'],values,color='#386e94')
    ax.set(ylabel='Global heat-operator error ceiling on collar data',
           title='Atomic parameters: the local estimate controls a global error')
    for i,v in enumerate(values): ax.text(i,v,f'{v:.5f}',ha='center',va='bottom')
    ax.set_ylim(0,max(values)*1.18)
    fig.supxlabel('Auxiliary Euclidean heat parameter, not a physical travel duration.\nExterior propagation retained. Supporting stress and atom passage remain unproved.',fontsize=9)
    for ext in ('png','pdf'): fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__': run()
