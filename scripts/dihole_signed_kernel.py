"""Signed positive-time frozen reference on the actual field; quadrature pilot.

This is neither a quadrature enclosure nor a renormalized supporting source.
Spatial-only cutoff: static coefficients permit exact Gaussian time integration.
Historical 187 kernels/results are unchanged; this is a new reference choice.
"""
import itertools
import json
from time import perf_counter

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext
from dihole_heat_defect import unpack
from dihole_local_transport import packed
from dihole_pointwise_dirac import PointwiseDirac, self_check
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_signed_kernel'


def gaussian_rule(n):
    if n < 2:
        raise ValueError('At least two nodes required for Gaussian second moments.')
    x, w = mp.gauss_quadrature(n, 'hermite')
    return [(mp.sqrt(2)*x[j], w[j]/mp.sqrt(mp.pi)) for j in range(n)]


def unit_rule(n):
    x, w = mp.gauss_quadrature(n, 'legendre')
    return [((x[j]+1)/2, w[j]/2) for j in range(n)]


def freeze(point, q, t):
    if q <= 0 or t <= 0:
        raise ValueError('Positive flux and heat time required.')
    E, I = point['E'][1:4,1:4], point['I'][1:4,1:4]
    field = q*I.T*point['F'][1:4,1:4]*I
    b = mp.sqrt(-sum((field*field)[j,j] for j in range(3))/2)
    if not b:
        raise ValueError('This pilot requires a nonzero magnetic field.')
    beta = b*t
    c = beta/mp.tanh(beta)
    perp = -field*field/b**2
    normal_map = mp.eye(3)+(1/mp.sqrt(c)-1)*perp
    precision = mp.eye(3)+(c-1)*perp
    prefactor = mp.det(point['E'])/(64*mp.pi**2*(t*q)**2)*(beta/mp.sinh(beta))**2/c
    return dict(E=E, I=I, F=field, b=b, beta=beta, c=c, perp=perp,
                normal_map=normal_map, precision=precision, prefactor=prefactor,
                t=t, q=q, point=point)


def density(frozen, output, z_scaled, cutoff=1, cutoff_gradient=None, delta_A=None):
    """Time-integrated signed density / q^2 for one spatial Gaussian draw.

    v and w below are sqrt(t) times their physical coefficients. The common
    source radial-gauge phase cancels before this expression is evaluated.
    """
    t = frozen['t']; root_t = mp.sqrt(t)
    gradient = mp.zeros(3,1) if cutoff_gradient is None else cutoff_gradient
    delta_A = mp.zeros(3,1) if delta_A is None else delta_A
    p = (root_t*gradient-cutoff*frozen['E'].T
         *(frozen['precision']-1j*t*frozen['F'])*z_scaled/2
         -1j*cutoff*root_t*delta_A)
    inverse = output['I'][1:4,1:4]
    v = inverse*p+root_t*cutoff*output['divI'][1:4,0]/2
    f = output['f']
    w = mp.matrix([f*inverse[0,0]*p[0]+root_t*cutoff/2
        *(output['f_u']*inverse[0,0]+f*output['d_u_inverse']), 0, 0])
    spatial = mp.re((v.H*w)[0])
    time = -f*cutoff**2*(frozen['point']['N']/output['N'])**2/4
    pauli = mp.re(1j*mp.sinh(2*frozen['beta'])
                  *(v.H*(frozen['F']/frozen['b'])*w)[0])
    return -16*frozen['prefactor']*(mp.cosh(2*frozen['beta'])*(spatial+time)+pauli)


def uniform_density(frozen):
    """Exact uniform-coefficient control, including canonical source volume."""
    b, t, q = (frozen[k] for k in ('b','t','q'))
    perpendicular = sum(frozen['F'][0,j]**2 for j in (1,2))
    x = 2*t*b
    return (mp.det(frozen['point']['E'])*frozen['point']['f']*perpendicular
            /(4*mp.pi**2*(2*t)*b*q**2)*(1/mp.tanh(x)-x/mp.sinh(x)**2))


def shifted_radius(evaluator, u, radial, du, rule):
    """Integrate the short radial displacement, retaining the actual smooth step."""
    target = u+du
    L = mp.mpf(1)/20
    if target <= 0 or target >= L:
        return evaluator.radial(target)
    if not 0 < u < L:
        return evaluator.radial(target)
    radius = radial[0]+du*mp.fsum(w*evaluator.step((u+v*du)/L)[0] for v,w in rule)
    rate, first, _ = evaluator.step(target/L)
    return radius, rate, first/L


def source_expectation(evaluator, u, mu, pole, q, t, normal, line):
    y = mp.matrix([u, mp.sqrt((1-mu)/(1+mu)), 0])
    radial = evaluator.radial(u)
    source = evaluator.point(*y, pole=pole, radial=radial)
    frozen = freeze(source,q,t)
    reference_point = dict(source, divI=mp.zeros(4,1), f_u=mp.mpf(0), d_u_inverse=mp.mpf(0))
    values = []; max_displacement = mp.mpf(0)
    radius = mp.mpf(2)**-32
    for entries in itertools.product(normal, repeat=3):
        standard = mp.matrix([e[0] for e in entries])
        weight = mp.fprod(e[1] for e in entries)
        z = frozen['normal_map']*standard
        d = mp.sqrt(t)*frozen['I']*z
        distance = mp.norm(d)
        max_displacement = max(max_displacement,distance)
        reference = density(frozen,reference_point,z)
        x = y+d
        if distance >= radius or not 0 <= x[0] <= mp.mpf(1)/20:
            values.append(-weight*reference)
            continue
        if distance <= radius/2:
            cutoff, gradient = mp.mpf(1), mp.zeros(3,1)
        else:
            s, first, _ = evaluator.step(2*distance/radius-1)
            cutoff, gradient = 1-s, -2*first/radius*d/distance
        output = evaluator.point(*x,pole=pole,
            radial=shifted_radius(evaluator,u,radial,d[0],line))
        # Exact curvature line-integral formula, evaluated by declared quadrature.
        delta = mp.zeros(3,1)
        for v,w in line:
            location = y+v*d
            middle = evaluator.point(*location,pole=pole,
                radial=shifted_radius(evaluator,u,radial,v*d[0],line))
            derivative = sum((d[k]*middle['dF'][k+1][1:4,1:4] for k in range(3)),mp.zeros(3))
            delta += q*w*(1-v*v)/2*derivative.T*d
        actual = density(frozen,output,z,cutoff,gradient,delta)
        values.append(weight*(actual-reference))
    # The reference is quadratic in the Gaussian coordinates, hence the chosen
    # normal rule integrates it exactly. Subtraction preserves small differences.
    uniform = uniform_density(frozen)
    return uniform+mp.fsum(values), max_displacement, uniform


def integrate(evaluator,q,T,n_u,n_mu,n_normal=3,n_line=4):
    start = perf_counter()
    normal, line = gaussian_rule(n_normal),unit_rule(n_line)
    total = []; corrections = []; maximum = mp.mpf(0)
    radius,L = mp.mpf(2)**-32,mp.mpf(1)/20
    # Explicitly retain both sides of the input buffer; no exterior source drop.
    for left,right in ((-2*radius,0),(0,L),(L,L+2*radius)):
        for v,wu in unit_rule(n_u):
            u = left+(right-left)*v
            for mu,wm in unit_rule(n_mu):
                for pole in (-1,1):
                    value,distance,uniform = source_expectation(evaluator,u,mu,pole,q,T/2,normal,line)
                    weight = 4*mp.pi*wu*wm*(right-left)/(1+mu)**2
                    total.append(weight*value)
                    corrections.append(weight*(value-uniform))
                    maximum = max(maximum,distance)
    return dict(n_u=n_u,n_mu=n_mu,n_normal=n_normal,n_line=n_line,dps=mp.dps,
        signed_reference_per_q2=str(mp.fsum(total)),
        actual_minus_uniform_reference_per_q2=str(mp.fsum(corrections)),
        maximum_sampled_coordinate_displacement=str(maximum),
        seconds=perf_counter()-start,quadrature_error_certified=False)


def uniform_integral(evaluator,q,T,n_u,n_mu):
    """Cheap source control variate; no claim that it is the full reference."""
    values=[]
    for v,wu in unit_rule(n_u):
        u=v/20; radial=evaluator.radial(u)
        for mu,wm in unit_rule(n_mu):
            for pole in (-1,1):
                point=evaluator.point(u,mp.sqrt((1-mu)/(1+mu)),0,pole,radial)
                values.append(mp.pi*wu*wm/(5*(1+mu)**2)*uniform_density(freeze(point,q,T/2)))
    return dict(n_u=n_u,n_mu=n_mu,value_per_q2=str(mp.fsum(values)))


def quadrature_controls(q,T):
    rows=[]
    for precision,ng,nl in ((60,3,4),(60,4,4),(60,3,6),(80,3,4)):
        with mp.workdps(precision):
            evaluator=PointwiseDirac()
            value,_,_=source_expectation(evaluator,mp.mpf('.021'),mp.mpf('.37'),1,
                q,T/2,gaussian_rule(ng),unit_rule(nl))
            rows.append(dict(dps=precision,n_normal=ng,n_line=nl,value=str(value)))
    with mp.workdps(80):
        evaluator=PointwiseDirac(); u=mp.mpf('.021'); du=mp.mpf('1e-13')
        actual=shifted_radius(evaluator,u,evaluator.radial(u),du,unit_rule(4))
        expected=evaluator.radial(u+du)
        radial_error=mp.norm(mp.matrix(actual)-mp.matrix(expected))
    return dict(point_refinements=rows,short_radius_integral_absolute_error=str(radial_error),
        scope='Finite diagnostics at one interior source, not an integral error enclosure.')


def signed_error(iv,row,T):
    """190 two-leg error about the EXACT new reference, not its quadrature."""
    gamma = iv.sqrt(1+80*iv.mpf(2)**-32)
    a,b = (unpack(iv,row['geometry'][k]) for k in ('insertion_graph_a','insertion_graph_b'))
    F0,F1,Q1 = (unpack(iv,row[k]) for k in
        ('value_HS_error_upper','Dirac_HS_error_upper','frozen_Dirac_HS_upper'))
    X0,X1 = (unpack(iv,x) for x in row['selected_heat_HS_upper'])
    error = 2*T*(gamma*F1*(a*X1+b*X0)+gamma**2*Q1*(a*F1+b*F0)+F1*(a*F1+b*F0))
    return packed(iv,error/iv.mpf(row['q'])**2)


def controls():
    errors=[]; normal=gaussian_rule(3)
    for fields in ((2,0,0),(0,0,2),(1,2,3)):
        curvature=mp.zeros(4)
        for (i,j),v in zip(((1,2),(1,3),(2,3)),fields):
            curvature[i,j],curvature[j,i]=v,-v
        point=dict(E=mp.diag([2,3,4,5]),I=mp.diag([mp.mpf(1)/v for v in (2,3,4,5)]),
            F=curvature,N=mp.mpf(2),divI=mp.zeros(4,1),f=mp.mpf('.7'),f_u=0,d_u_inverse=0)
        frozen=freeze(point,mp.mpf(11),mp.mpf('.03'))
        observed=mp.fsum(mp.fprod(e[1] for e in entries)*density(frozen,point,
            frozen['normal_map']*mp.matrix([e[0] for e in entries]))
            for entries in itertools.product(normal,repeat=3))
        errors.append(abs(observed-uniform_density(frozen)))
    # Direct four-dimensional Clifford trace verifies the Pauli sign independently.
    from dihole_uv_bound import gamma_matrices
    gamma=[mp.matrix(g.tolist()) for g in gamma_matrices()]
    field=mp.zeros(4); field[1:4,1:4]=frozen['F']
    V=1j/4*sum(((gamma[a]*gamma[b]-gamma[b]*gamma[a])*field[a,b]
                for a in range(4) for b in range(4)),mp.zeros(4))
    v=mp.matrix([1+2j,3-1j,mp.mpf('.2'),-.5j]); w=mp.matrix([.7j,2,1-1j,3])
    gv=sum((v[a]*gamma[a] for a in range(4)),mp.zeros(4))
    gw=sum((w[a]*gamma[a] for a in range(4)),mp.zeros(4))
    matrix=gv.H*gw*mp.expm(-2*frozen['t']*V)
    direct=sum(matrix[j,j] for j in range(4))
    scalar=4*mp.cosh(2*frozen['beta'])*(v.H*w)[0]+4j*mp.sinh(2*frozen['beta'])/frozen['b']*(v.H*field*w)[0]
    # Independent variable-field gauge control: A_2=B*x_1+g*x_1^2/2+partial_2 phi,
    # A_1=partial_1 phi, phi=h*x_1^2*x_2. Direct line-gauge differentiation must
    # reproduce the curvature formula including its orientation and sign.
    B,g,h=map(mp.mpf,('1.3','.7','.4'))
    y=mp.matrix([mp.mpf('.2'),mp.mpf('-.3')]); d=mp.matrix([mp.mpf('.4'),mp.mpf('-.1')]); x=y+d
    def potential(x): return mp.matrix([2*h*x[0]*x[1],B*x[0]+g*x[0]**2/2+h*x[0]**2])
    def phase(x):
        displacement=x-y
        return mp.fsum(w*((displacement.T*potential(y+v*displacement))[0]) for v,w in unit_rule(3))
    phase_gradient=mp.matrix([mp.diff(lambda z:phase(mp.matrix([z,x[1]])),x[0]),
                              mp.diff(lambda z:phase(mp.matrix([x[0],z])),x[1])])
    frozen_A=mp.matrix([-(B+g*y[0])*d[1]/2,(B+g*y[0])*d[0]/2])
    delta_A=mp.matrix([-g*d[0]*d[1]/3,g*d[0]**2/3])
    gauge_error=mp.norm(potential(x)-phase_gradient-frozen_A-delta_A)
    return dict(uniform_field_max_absolute_error=str(max(errors)),
        Clifford_trace_absolute_error=str(abs(direct-scalar)),
        variable_field_pure_gauge_absolute_error=str(gauge_error),
        uniform_parallel_field_absolute_error=str(errors[1]))


def run():
    mp.dps=60
    point_checks=self_check()
    control=controls()
    base_path=project_path('dihole_heat_source_graph_results.json')
    base=json.loads(base_path.read_text())
    heat_path=project_path('dihole_regulator_heat_tail_results.json')
    heat=json.loads(heat_path.read_text())
    inputs={**base['input_sha256'],base_path.name:sha(base_path),
        'dihole_heat_source_graph.py':base['source_sha256'],
        'dihole_pointwise_dirac.py':sha(project_path('dihole_pointwise_dirac.py')),
        heat_path.name:sha(heat_path)}
    if not base['passed'] or not heat['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing input provenance required.')
    row=max((r for r in base['records'] if r['heat_time_fraction']==1),key=lambda r:r['q'])
    q=row['q']; hr=next(r for r in heat['records'] if r['q']==q)
    T=mp.mpf(hr['cutoff_s0_over_D_squared'])
    evaluator=PointwiseDirac()
    records=[]
    for nu,nm in ((6,4),(10,6),(14,8)):
        record=integrate(evaluator,q,T,nu,nm)
        records.append(record)
        print(json.dumps(record),flush=True)
    uniform_records=[uniform_integral(evaluator,q,T,nu,nm) for nu,nm in ((16,12),(32,24),(64,48))]
    for record in records:
        record['control_variate_reference_per_q2']=str(mp.mpf(uniform_records[-1]['value_per_q2'])
            +mp.mpf(record['actual_minus_uniform_reference_per_q2']))
    refinements=quadrature_controls(q,T)
    iv=MPIntervalContext();iv.dps=100
    error=signed_error(iv,row,iv.mpf(hr['cutoff_s0_over_D_squared']))
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        pointwise_operator=point_checks['passed'],
        uniform_field_control=mp.mpf(control['uniform_field_max_absolute_error'])<mp.mpf('1e-50'),
        Clifford_trace=mp.mpf(control['Clifford_trace_absolute_error'])<mp.mpf('1e-50'),
        variable_field_gauge=mp.mpf(control['variable_field_pure_gauge_absolute_error'])<mp.mpf('1e-50'),
        short_radius_quadrature=mp.mpf(refinements['short_radius_integral_absolute_error'])<mp.mpf('1e-60'),
        all_results_finite=all(mp.isfinite(mp.mpf(r['signed_reference_per_q2'])) for r in records),
        pilot_not_claimed_certified=all(not r['quadrature_error_certified'] for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,
        pointwise_checks={k:(v if isinstance(v,bool) else str(v)) for k,v in point_checks.items()},
        q=q,T=str(T),records=records,uniform_source_refinements=uniform_records,
        quadrature_diagnostics=refinements,error_about_exact_reference_per_q2_upper=error,
        reference_cutoff='chi(|x_sp-y_sp|), radius 2^-32; time displacement uncut.',
        source_domain='Both collars, u in [-2r,.05+2r], both source hemispheres.',
        signed_full_heat_interval=None,renormalized_source=None,atom_transit=None,
        source_sha256=sha(__file__),input_sha256=inputs,
        scope=['Numerical signed reference with actual output coefficients and curvature line integral.',
            'Spatial cutoff replaces the historical 4D reference only in this new calculation.',
            'Exact-reference error radius retains the full exterior; quadrature error is NOT enclosed.',
            'This heat time is a Euclidean calculation parameter, not elapsed travel time.',
            'No supporting source, stability, atom transport, clock gain or new physical law.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:
        raise ArithmeticError('Signed-kernel controls failed; results retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'],checks=checks,error=error)),flush=True)


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3),layout='constrained')
    records=result['records']
    axes[0].plot([r['n_u'] for r in records],
        [float(r['signed_reference_per_q2'])*1e8 for r in records],'o-',color='#247886',label='Direct quadrature')
    estimate=float(records[-1]['control_variate_reference_per_q2'])
    axes[0].axhline(estimate*1e8,color='#73488e',linestyle='--',label='With refined control variate')
    axes[0].set(xlabel='Radial source nodes (angular nodes also increase)',
        ylabel='Signed reference / q squared, in units of 1e-8',title='Positive numerical reference')
    axes[0].legend(fontsize=8)
    axes[1].bar(['Reference magnitude','Full-field error ceiling'],
        [float(mp.log10(abs(estimate))),result['error_about_exact_reference_per_q2_upper']['log10_upper']],
        color=['#247886','#b05740'])
    axes[1].axhline(0,color='#888888',linewidth=.7)
    axes[1].set(ylabel='log10 of normalized magnitude',title='The full response sign remains unknown')
    fig.supxlabel('Positive Euclidean heat time; numerical quadrature is not an enclosure.\nThis is not supporting energy, an atom-transport result, or a time-travel demonstration.',fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
