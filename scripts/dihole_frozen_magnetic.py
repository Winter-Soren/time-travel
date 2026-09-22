"""Frozen magnetic first defect: local column L1, not a full heat/source error.

Keep the constant magnetic operator exact, factor differences in D0, and
whiten at the source before bounding the actual variable coefficients.
"""
import json

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import unpack
from dihole_local_transport import dprod, packed
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_frozen_magnetic'


def constants(iv, prior, frame, local, magnetic):
    real = local['real_geometry']
    M, H = (iv.mpf(real[k]) for k in ('metric_norm_upper', 'inverse_metric_norm_upper'))
    G1 = 8*iv.mpf(real['first_component_upper'])
    radius = iv.mpf(2)**-32
    eps = H*G1*radius
    if eps.b >= 1: raise ArithmeticError('Frozen metric comparison fails.')
    a = 1/iv.sqrt(1-eps)
    alpha = 2*iv.sqrt(M*H)*unpack(iv, frame['frame_bounds']['internal_inverse_frame_operator_upper'][1])
    fields = prior['jets']['total_order_derivative_upper']
    E, I = (list(map(iv.mpf, frame['frame_bounds'][k])) for k in ('coframe', 'inverse_coframe'))
    V, Vi = (list(map(iv.mpf, fields[k])) for k in ('half_density', 'inverse_half_density'))
    spin = dprod(dprod(I, I), E[1:])
    density = dprod(I, dprod(V[1:], Vi))
    geom = [768*x+16*y for x, y in zip(spin, density)]
    S1, S2 = map(iv.mpf, prior['step_derivative_upper'][1:3])
    return dict(radius=radius, epsilon=eps, a=a, alpha=alpha,
        metric_difference=H*iv.sqrt(H)*G1/(1-eps),
        b0=geom[0], b1=2*iv.sqrt(H)*geom[1],
        curvature_gradient_per_q=16*H*iv.sqrt(H)*iv.mpf(fields['gauge_per_abs_q'][2]),
        magnetic_per_q=iv.mpf(magnetic['geometry_constants']['magnetic_unit_flux_upper_D2']),
        cutoff_gradient=2*iv.sqrt(H)*S1/radius,
        cutoff_hessian=4*H*(S2+S1)/radius**2,
        cutoff_exponent=radius**2/(32*H))


def first_defect(iv, c, q, s0):
    """Positive majorants and their time integrals; no row estimate assumed."""
    a, alpha, b0, b1 = (c[k] for k in ('a', 'alpha', 'b0', 'b1'))
    u = iv.mpf(1)/100
    coth = 1+u*u/3  # beta*coth(beta) <= 1+beta^2/3, 0<=beta<=u.
    ell = coth+u
    L = abs(q)*c['curvature_gradient_per_q']
    b = abs(q)*c['magnetic_per_q']
    if (b*s0).b > u.a: raise ArithmeticError('Magnetic time window exceeded.')
    b2 = a*L/3
    b11, b12 = 2*a*L/3, alpha*L/3
    p1 = 2*b*(a+1)*alpha
    V0 = b0*b0+4*a*b1
    V1 = p1+4*a*b11
    V2 = 2*b0*b2+4*a*b12
    V4 = b2*b2
    rootpi = iv.sqrt(iv.pi)
    cm = rootpi*(c['metric_difference']*(15*(coth*coth+u*u)/4+3*coth)
                 +iv.mpf(3)/2*a*ell*(2*alpha+b0))
    ch = 15*rootpi*a*ell*b2+iv.mpf(3)/2*rootpi*V1
    integrals = dict(leading=2*cm*iv.sqrt(s0), constant=V0*s0,
        varying_field=iv.mpf(2)/3*ch*s0*iv.sqrt(s0),
        quadratic=4*V2*s0**2, quartic=32*V4*s0**3)
    integrals = {k: iv.exp(u)*v for k, v in integrals.items()}
    interior = sum(integrals.values())
    # Free Gaussian tail splitting: e^-r²/(4t) <= e^-d²/(8t)e^-r²/(8t).
    # The widened normalized 4D Gaussian contributes 4*2^(k/2) to moment k.
    Gchi, Hchi = c['cutoff_gradient'], c['cutoff_hessian']
    tail = iv.exp(-c['cutoff_exponent']/s0)
    cut = dict(gradient=2*a*a*Gchi*ell*(iv.mpf(3)/2*rootpi)*4*iv.sqrt(2)*iv.sqrt(s0),
        constant=(4*a*a*Hchi+4*a*alpha*Gchi+2*a*b0*Gchi)*4*s0,
        quadratic=2*a*b2*Gchi*8*8*s0**2/2)
    cut = {k: iv.exp(u)*tail*v for k, v in cut.items()}
    cutoff = sum(cut.values())
    return dict(q=q, s0=upper(s0), magnetic_product_upper=upper(b*s0),
        interior_terms={k: packed(iv, v) for k, v in integrals.items()},
        interior_column_time_integral_upper=packed(iv, interior),
        cutoff_terms={k: packed(iv, v) for k, v in cut.items()},
        cutoff_column_time_integral_upper=packed(iv, cutoff),
        cutoff_exponent_over_s0=packed(iv, c['cutoff_exponent']/s0),
        total_column_time_integral_upper=packed(iv, interior+cutoff),
        total_less_than_one=bool((interior+cutoff).b < 1),
        coefficients={k: packed(iv, v) for k, v in dict(b2=b2, b11=b11, b12=b12,
            p1=p1, V0=V0, V1=V1, V2=V2, V4=V4, Cminus=cm, Chalf=ch).items()},
        row_bound=None, iterated_heat_remainder=None, physical_source=None, atom_transit=None)


def controls():
    """Independent differential identities, with a varying field and frame."""
    import sympy as sp
    x, y, b, h, v = sp.symbols('x y b h v', real=True)
    coords = (x, y)
    gamma = (sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]))
    aa = ((1+h*x)*gamma[0], gamma[1])
    a0 = (-b*y/2, b*x/2)
    delta = (-v*x*y/3, v*x*x/3)
    bg = h*gamma[0]/2
    B = bg-sp.I*sum((aa[i]*delta[i] for i in range(2)), sp.zeros(2))
    F = sp.Matrix([[0, b], [-b, 0]])
    comm = lambda a, c: a*c-c*a
    D0 = lambda f, i: f.diff(coords[i])-sp.I*a0[i]*f
    D = lambda f: sum((aa[i]*D0(f, i) for i in range(2)), sp.zeros(2, 1))+B*f
    pauli = sp.I/4*sum((comm(gamma[i], gamma[j])*F[i, j]
        for i in range(2) for j in range(2)), sp.zeros(2))
    P0 = lambda f: -sum((D0(D0(f, i), i) for i in range(2)), sp.zeros(2, 1))+pauli*f
    first = [-sum((aa[i]*aa[k].diff(coords[i]) for i in range(2)), sp.zeros(2))
             -aa[k]*B-B*aa[k] for k in range(2)]
    zero = sp.I/4*sum(((comm(aa[i], aa[j])-comm(gamma[i], gamma[j]))*F[i, j]
        for i in range(2) for j in range(2)), sp.zeros(2))
    zero -= sum((aa[i]*B.diff(coords[i]) for i in range(2)), sp.zeros(2))+B*B
    f = sp.Matrix([1+x*x+y, x*y+2*y*y])
    principal = ((1+h*x)**2, sp.Integer(1))
    factored = -sum(((principal[i]-1)*D0(D0(f, i), i) for i in range(2)), sp.zeros(2, 1))
    factored += sum((first[i]*D0(f, i) for i in range(2)), sp.zeros(2, 1))+zero*f
    residual = (-D(D(f))-P0(f)-factored).applyfunc(sp.simplify)
    chi = 1+x*y+x*x
    cutoff = sum((-2*principal[i]*chi.diff(coords[i])*D0(f, i)
        -principal[i]*chi.diff(coords[i], 2)*f+first[i]*chi.diff(coords[i])*f
        for i in range(2)), sp.zeros(2, 1))
    cutoff_residual = (-D(D(chi*f))+chi*D(D(f))-cutoff).applyfunc(sp.simplify)
    mp.dps = 90
    errors, derivative_errors = [], []
    def kernel(t, z, field, spin):
        beta = field*t
        Q = field*mp.coth(beta) if field else 1/t
        pref = beta/mp.sinh(beta) if field else 1
        return pref*mp.exp(-(Q*(z[0]**2+z[1]**2)+(z[2]**2+z[3]**2)/t)/4-t*spin)/(4*mp.pi*t)**2
    z = list(map(mp.mpf, ('.11', '-.07', '.05', '.03')))
    for field in (mp.mpf(0), mp.mpf(2)):
        t = mp.mpf('.005')
        Q = mp.diag([field*mp.coth(field*t) if field else 1/t]*2+[1/t]*2)
        F0 = mp.matrix([[0, field, 0, 0], [-field, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        C = Q-1j*F0
        for spin in (-field, field):
            f0 = kernel(t, z, field, spin)
            lap = mp.mpc(0)
            for i in range(4):
                fun = lambda w: kernel(t, [w if k==i else z[k] for k in range(4)], field, spin)
                A = -(F0*mp.matrix(z))[i]/2
                d1 = mp.diff(fun, z[i])-1j*A*f0
                d2 = mp.diff(fun, z[i], 2)-2j*A*mp.diff(fun, z[i])-A*A*f0
                lap += d2
                derivative_errors.extend((abs(d1+(C*mp.matrix(z))[i]*f0/2),
                    abs(d2-((C*mp.matrix(z))[i]**2/4-C[i,i]/2)*f0)))
            errors.append(abs(mp.diff(lambda s: kernel(s, z, field, spin), t)-lap+spin*f0))
    beta = mp.mpf('.01')
    # Check a nondiagonal constant metric in the original coordinates too.
    E = mp.matrix([[2, mp.mpf('.1'), 0, 0], [mp.mpf('.1'), 3, mp.mpf('.2'), 0],
                   [0, mp.mpf('.2'), 4, mp.mpf('.1')], [0, 0, mp.mpf('.1'), 5]])
    gx = (E.T*E)**-1
    Fx = E.T*F0*E
    x0 = list(E**-1*mp.matrix(z))
    raw = lambda t, x: mp.det(E)*kernel(t, list(E*mp.matrix(x)), field, -field)
    kx = raw(t, x0)
    A = -Fx*mp.matrix(x0)/2
    first = [mp.diff(lambda w:raw(t, [w if k==i else x0[k] for k in range(4)]), x0[i]) for i in range(4)]
    lap = mp.mpc(0)
    for i in range(4):
        for j in range(4):
            if i==j:
                second = mp.diff(lambda w:raw(t, [w if k==i else x0[k] for k in range(4)]), x0[i], 2)
            else:
                second = mp.diff(lambda w,v:raw(t, [w if k==i else v if k==j else x0[k] for k in range(4)]),
                                 (x0[i],x0[j]), (1,1))
            lap += gx[i,j]*(second+1j*Fx[j,i]*kx/2-1j*A[j]*first[i]-1j*A[i]*first[j]-A[i]*A[j]*kx)
    metric_error = abs(mp.diff(lambda s:raw(s,x0),t)-lap-field*kx)
    Qx = E.T*Q*E
    gaussian_mass = mp.det(E)*beta/mp.sinh(beta)*mp.exp(beta)/(4*mp.pi*t)**2*mp.pi**2/mp.sqrt(mp.det(Qx/4))
    normalization = abs(gaussian_mass-(1+mp.tanh(beta)))
    return dict(operator_factorization=residual==sp.zeros(2, 1),
        varying_field_preserved=sp.simplify(sp.diff(a0[1]+delta[1], x)-sp.diff(a0[0]+delta[0], y)-b-v*x)==0,
        constant_interior_defect_zero=factored.subs({h:0, v:0})==sp.zeros(2, 1),
        cutoff_commutator=cutoff_residual==sp.zeros(2, 1),
        frozen_PDE_error=float(max(errors)), covariant_derivative_error=float(max(derivative_errors)),
        nondiagonal_metric_PDE_error=float(metric_error),
        frozen_column_normalization_error=float(normalization))


def run():
    iv = MPIntervalContext(); iv.dps = 100
    path = project_path('dihole_radial_cancellation_results.json')
    previous = json.loads(path.read_text())
    inputs = {**previous['input_sha256'], path.name:sha(path),
        'dihole_radial_cancellation.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing186 inputs required.')
    prior, frame, local, magnetic, heat = (json.loads(project_path(n+'_results.json').read_text()) for n in
        ('dihole_regular_jets', 'dihole_positive_frame', 'dihole_local_transport',
         'dihole_magnetic_remainder', 'dihole_regulator_heat_tail'))
    c = constants(iv, prior, frame, local, magnetic)
    records = [first_defect(iv, c, r['q'], iv.mpf(r['cutoff_s0_over_D_squared'])) for r in heat['records']]
    control = controls()
    checks = dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        local_metric_comparison=bool(c['epsilon'].b<1),
        support_within_chart=bool(c['radius'].b<iv.mpf(local['normal_geometry']['coordinate_clearance']).a),
        radial_gauge_support=bool(c['radius']==iv.mpf(2)**-32),
        magnetic_windows=all(r['magnetic_product_upper']<=.01 for r in records),
        operator_factorization=control['operator_factorization'],
        varying_field_preserved=control['varying_field_preserved'],
        constant_interior_defect_zero=control['constant_interior_defect_zero'],
        cutoff_commutator=control['cutoff_commutator'],
        frozen_PDE=control['frozen_PDE_error']<1e-75,
        covariant_derivatives=control['covariant_derivative_error']<1e-75,
        nondiagonal_metric_PDE=control['nondiagonal_metric_PDE_error']<1e-75,
        frozen_normalization=control['frozen_column_normalization_error']<1e-75,
        atomic_first_bound_useful=records[-1]['total_less_than_one'],
        low_flux_failure_retained=not records[0]['total_less_than_one'],
        no_full_remainder_or_transport_claimed=all(r['iterated_heat_remainder'] is None
            and r['physical_source'] is None and r['atom_transit'] is None for r in records))
    result = dict(passed=all(checks.values()), checks=checks, controls=control,
        constants={k:packed(iv,v) for k,v in c.items()}, records=records,
        source_sha256=sha(__file__), input_sha256=inputs, bound_encoding=previous['bound_encoding'],
        scope=['Fixed-source canonical Lebesgue kernel, source-centred radial gauge, frozen whitened coordinates; D=1.',
            'Exact constant magnetic kernel includes its Pauli endomorphism; actual field derivatives, geometry, spin and density remain in the defect.',
            'Uniform local column-L1 first-defect time integral only, including the local cutoff. Not a global Schur or iterated heat remainder.',
            'The earlier B4 certificate is preserved and is not replaced by this dimensionless first-defect estimate.',
            'No source, supported opening, atom passage or time travel established. One hydrogen atom and trial1nm opening remain the target.'],
        sources=['https://arxiv.org/abs/0810.4889'])
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if not result['passed']: raise ArithmeticError('Frozen magnetic checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True, checks=checks, controls=control,
        bounds=[dict(q=r['q'], total_log10=r['total_column_time_integral_upper']['log10_upper'],
            cutoff_log10=r['cutoff_column_time_integral_upper']['log10_upper']) for r in records]), indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r = result['records'][-1]
    values = [10**x['log10_upper'] for x in r['interior_terms'].values()]
    fig, ax = plt.subplots(figsize=(8,4.5), layout='constrained')
    labels = ['Leading geometry', 'Constant', 'Varying field', 'Quadratic', 'Quartic']
    ax.barh(labels, values, color=['#386e94', '#9bafbc', '#b66f35', '#9bafbc', '#9bafbc'])
    ax.set(xlabel='Contribution to first-defect column integral (dimensionless)',
        title='One-atom target: a useful first field-error bound')
    ax.text(.98,.04, 'Total < 0.050\nLocal cutoff < 10^-15700', transform=ax.transAxes,
        ha='right', va='bottom')
    fig.supxlabel('This is one intermediate bound, not a transit probability or a complete error bound.\nSupport, stability, atom survival and clock gain remain unproved.', fontsize=9)
    for ext in ('png','pdf'): fig.savefig(project_path(NAME+'_progress.'+ext), dpi=170)
    plt.close(fig)


if __name__=='__main__': run()
