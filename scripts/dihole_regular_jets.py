"""Certified high jets of the actual collar in regular local frames.

These are coefficient bounds, not a heat-defect bound or a supporting source.
Lengths and coordinate components use D=1. Gauge bounds are per abs(q).
"""
import json
import math
from types import SimpleNamespace

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, calibration, field_function, lower, upper
from dihole_slice_channel import sha
from dihole_stretched_collar import smooth_step
from project_paths import project_path


NAME = 'dihole_regular_jets'
ORDER = 20


def complex_cover(iv, radial=24, angular=32):
    """Cauchy polydiscs over the whole real buffer, including both axes."""
    guards = dict(root_real=math.inf, NP_real=math.inf, W_real=math.inf)

    def root(z):
        # mpmath's interval complex sqrt is unavailable; this is the
        # principal root on Re(z)>0, evaluated with real interval squares.
        if z.real.a <= 0:
            raise ArithmeticError('Complex root branch is not separated.')
        guards['root_real'] = min(guards['root_real'], lower(z.real))
        u = iv.sqrt((iv.sqrt(z.real**2+z.imag**2)+z.real)/2)
        return iv.mpc(u, z.imag/(2*u))

    function = field_function(SimpleNamespace(sqrt=root))
    parameters = calibration(iv)
    rr, rx = iv.mpf(1)/1000, iv.mpf(1)/100
    bounds = dict.fromkeys(('N', 'P', 'N_inverse', 'P_inverse', 'a'), 0.)
    for i in range(radial):
        R = iv.mpc(iv.mpf([iv.mpf(1)/40+iv.mpf(9)*i/(200*radial)-rr,
                          iv.mpf(1)/40+iv.mpf(9)*(i+1)/(200*radial)+rr]), iv.mpf([-rr, rr]))
        for j in range(angular):
            x = iv.mpc(iv.mpf([-1+iv.mpf(2)*j/angular-rx,
                              -1+iv.mpf(2)*(j+1)/angular+rx]), iv.mpf([-rx, rx]))
            f, W, a = function(R, x, *parameters)[:3]
            if f.real.a <= 0 or W.real.a <= 0:
                raise ArithmeticError('Complex metric factors are not separated from zero.')
            guards['NP_real'] = min(guards['NP_real'], lower(f.real))
            guards['W_real'] = min(guards['W_real'], lower(W.real))
            N, P = root(f/W), root(f*W)
            for name, value in zip(bounds, (N, P, 1/N, 1/P, a)):
                bound = upper(abs(value))
                if not math.isfinite(bound):
                    raise ArithmeticError('A complex coefficient enclosure is not finite.')
                bounds[name] = max(bounds[name], bound)
    return dict(radial_tiles=radial, angular_tiles=angular, boxes=radial*angular,
                radius_R='1/1000', radius_x='1/100', bounds=bounds, positive_guards=guards)


def step_bounds(iv, order):
    """Cauchy on endpoint discs and a fixed middle disc; see report proof."""
    C = iv.exp(iv.mpf(5)/3)/(1-iv.exp(-iv.mpf(5)/6))
    result = [1.]
    for n in range(1, order+1):
        endpoint = C*(2*n/iv.exp(1))**n
        result.append(upper(math.factorial(n)*max(iv.mpf(10)**n, endpoint)))
    return result


def chart_bounds(iv, base):
    """Two oriented stereographic disks |y|^2<=2; all polar quotients regular."""
    N, P, Ni, Pi, a = (iv.mpf(base[k]) for k in ('N','P','N_inverse','P_inverse','a'))
    eta, rx = iv.mpf(1)/2000, iv.mpf(1)/100
    displacement = 4*eta+2*eta**2
    denominator = 1-displacement
    x_displacement = 2*displacement/denominator
    if x_displacement.b >= (rx/2).a:
        raise ArithmeticError('Chart polydiscs leave the inner analytic x tube.')
    Y = iv.sqrt(2)+eta
    rhi, rlo = iv.mpf(71)/1000, iv.mpf(24)/1000
    E = 2*rhi/denominator*(Ni+4*(P+Ni)*Y**2/(denominator*rx))
    Ei = (3+displacement)/(2*rlo)*(N+4*(Pi+N)*Y**2/(denominator*rx))
    constants = dict(coframe=max(N,P,E), inverse_coframe=max(Ni,Pi,Ei),
        metric=max(N*N,P*P,2*E*E), inverse_metric=max(Ni*Ni,Pi*Pi,2*Ei*Ei),
        gauge_per_abs_q=4*a*Y/(denominator*rx),
        half_density=2*P*rhi/denominator,
        inverse_half_density=Pi*(3+displacement)/(2*rlo),
        volume_density=4*P*P*rhi*rhi/denominator**2,
        inverse_volume_density=Pi*Pi*(3+displacement)**2/(4*rlo*rlo))
    return dict(radius_y='1/2000', squared_real_disk_radius=2,
                complex_x_displacement_upper=upper(x_displacement),
                complex_denominator_lower=lower(denominator),
                coefficient_modulus_upper={k:upper(v) for k,v in constants.items()})


def composed_bounds(iv, steps, chart):
    """Finite Faà di Bruno majorant; no analyticity of the sewn R(s) assumed."""
    L, rr, eta = iv.mpf(1)/20, iv.mpf(1)/1000, iv.mpf(1)/2000
    radii = [iv.mpf(0), iv.mpf(1)]
    radii += [iv.mpf(steps[n-1])/L**(n-1) for n in range(2, ORDER+1)]
    a = [radii[n]/(math.factorial(n)*rr) for n in range(ORDER+1)]
    b = [iv.mpf(1)]
    for n in range(1, ORDER+1):
        b.append(sum(a[j]*b[n-j] for j in range(1,n+1)))
    factors = [max(math.factorial(n)*math.factorial(k-n)*b[n]/eta**(k-n)
                   for n in range(k+1)) for k in range(ORDER+1)]
    fields = {name:[upper(iv.mpf(M)*factor) for factor in factors]
              for name,M in chart['coefficient_modulus_upper'].items()}
    masks = [upper(iv.mpf(steps[n+1])/(2*L**n)) for n in range(ORDER+1)]
    return dict(radius_derivative_upper=[upper(v) for v in radii[1:]],
                radial_normalized_majorant=[upper(v) for v in b],
                total_order_derivative_upper=fields, mask_derivative_upper=masks,
                mixed_formula='abs(d_s^n d_y1^j d_y2^k F)/(n! j! k!) <= M*b[n]*2000^(j+k), n+j+k<=20; b values are saved upward.',
                time_derivatives='All background coefficients and mask are independent of Euclidean time.')


def operator_bounds(iv, jets, q):
    """Canonical Dirac coefficients, P=-D^2 and the fixed null-probe variation."""
    fields=jets['total_order_derivative_upper']
    E,I,A,V,Vi=(list(map(iv.mpf,fields[k])) for k in
        ('coframe','inverse_coframe','gauge_per_abs_q','half_density','inverse_half_density'))
    F=list(map(iv.mpf,jets['mask_derivative_upper']))
    def product(a,b):
        return [sum(math.comb(n,j)*a[j]*b[n-j] for j in range(n+1))
                for n in range(min(len(a),len(b)))]
    def scale(c,a):return [c*x for x in a]
    def add(*rows):return [sum(values) for values in zip(*rows)]
    derivative=lambda a:a[1:]
    if product([1,2,6],[3,4,10]) != [3,10,44]:
        raise ArithmeticError('Total-derivative product control failed.')
    J=add(scale(abs(q),A),product(derivative(V),Vi))
    a=scale(4,I)
    b=add(scale(768,product(product(I,I),derivative(E))),scale(16,product(I,J)))
    dE,dI=product(F,E),product(F,I)
    da=scale(4,dI)
    db=add(scale(1536,product(product(I,dI),derivative(E))),
           scale(768,product(product(I,I),derivative(dE))),scale(16,product(dI,J)))
    second=product(a,a)
    first=add(scale(4,product(a,derivative(a))),scale(2,product(a,b)))
    zero=add(scale(4,product(a,derivative(b))),product(b,b))
    dsecond=scale(2,product(a,da))
    dfirst=add(scale(4,product(da,derivative(a))),scale(4,product(a,derivative(da))),
               scale(2,product(da,b)),scale(2,product(a,db)))
    dzero=add(scale(4,product(da,derivative(b))),scale(4,product(a,derivative(db))),scale(2,product(b,db)))
    return dict(q=q,derivatives_through=18,
        convention='P=sum_{mu,nu}P2[mu,nu]*partial_mu*partial_nu+sum_mu P1[mu]*partial_mu+P0; second coefficients use ordered pairs, not collected mixed derivatives.',
        variation='e0 -> exp(epsilon*f)*e0, e1 -> exp(-epsilon*f)*e1; angular coframe, gauge one-form and physical volume fixed.',
        coefficient_operator_norm_upper={
        name:[upper(x) for x in row[:19]] for name,row in
        (('P_second',second),('P_first',first),('P_zero',zero),
         ('deltaP_second',dsecond),('deltaP_first',dfirst),('deltaP_zero',dzero))})


def controls(iv, cover, steps, chart):
    mp.dps=85
    function=field_function(mp); parameters=calibration(mp)
    step=lambda t:1/(1+mp.exp(1/t-1/(1-t)))
    derivative_ratios=[]; value_errors=[]
    for t in map(mp.mpf, ('.01','.05','.15','.3','.5','.8','.95')):
        value_errors.append(abs(float(step(t))-float(smooth_step(float(t))[0])))
        for n in (1,2,4,8,12,16,20,21):
            derivative_ratios.append(float(abs(mp.diff(step,t,n))/mp.mpf(steps[n])))
    inverse_errors=[]; metric_errors=[]; polar_errors=[]; modulus_ratios=[]
    for R in map(mp.mpf, ('.025','.047','.07')):
        for pole in (-1,1):
            f,W,a,*_=function(R,pole,*parameters)
            polar_errors.append(float(max(abs(f-1),abs(a-(1-pole)/2))))
            for u,v in (('0','0'),('.2','-.4'),('1','1')):
                y=mp.matrix([mp.mpf(u),mp.mpf(v)]);t=(y.T*y)[0]
                x=pole*(1-t)/(1+t)
                raw=function(R,x,*parameters);f,W,a=raw[:3]
                N,P=mp.sqrt(f/W),mp.sqrt(f*W)
                E=2*R/(1+t)*(mp.eye(2)/N+(P-1/N)/t*y*y.T) if t else 2*R/N*mp.eye(2)
                Ei=(1+t)/(2*R)*(N*mp.eye(2)+(1/P-N)/t*y*y.T) if t else N/(2*R)*mp.eye(2)
                expected=4*R*R/(1+t)**2*(mp.eye(2)/N**2+(P**2-1/N**2)/t*y*y.T) if t else 4*R*R/N**2*mp.eye(2)
                inverse_errors.append(float(mp.norm(E*Ei-mp.eye(2))))
                metric_errors.append(float(mp.norm(E.T*E-expected)))
                gauge=pole*(a-(1-pole)/2)/t*mp.matrix([-y[1],y[0]]) if t else mp.matrix([0,0])
                for name,values in (('coframe',list(E)+[N,P]),('inverse_coframe',list(Ei)+[1/N,1/P]),('gauge_per_abs_q',list(gauge))):
                    modulus_ratios.append(float(max(abs(z) for z in values)/chart['coefficient_modulus_upper'][name]))
    # Independent coefficient convolution for a known rational composition.
    trial=[mp.mpf(0),mp.mpf(2),mp.mpf(3),mp.mpf(5)]
    b=[mp.mpf(1)]
    for n in range(1,9): b.append(sum(trial[j]*b[n-j] for j in range(1,min(n,3)+1)))
    exact=mp.taylor(lambda z:1/(1-2*z-3*z*z-5*z**3),0,8)
    return dict(step_derivative_max_bound_ratio=max(derivative_ratios), step_value_max_error=max(value_errors),
                frame_inverse_max_error=max(inverse_errors), frame_metric_max_error=max(metric_errors),
                polar_identity_max_error=max(polar_errors), chart_modulus_max_sample_ratio=max(modulus_ratios),
                composition_coefficient_error=float(max(abs(x-y) for x,y in zip(b,exact))))


def run():
    iv=MPIntervalContext();iv.dps=65
    prior_path=project_path('dihole_magnetic_remainder_results.json')
    prior=json.loads(prior_path.read_text())
    inputs={**prior['input_sha256'],prior_path.name:sha(prior_path),
            'dihole_magnetic_remainder.py':prior['source_sha256']}
    if not prior['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing inputs required.')
    cover=complex_cover(iv);steps=step_bounds(iv,ORDER+1)
    chart=chart_bounds(iv,cover['bounds']);jets=composed_bounds(iv,steps,chart)
    operators=[operator_bounds(iv,jets,row['q']) for row in prior['records']]
    check=controls(iv,cover,steps,chart)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        positive_complex_branches=all(v>0 for v in cover['positive_guards'].values()),
        whole_sphere_and_buffer=cover['boxes']==768 and chart['squared_real_disk_radius']==2,
        chart_displacement=chart['complex_x_displacement_upper']<.005,
        step_matches_actual=check['step_value_max_error']<1e-14,
        step_derivative_controls=check['step_derivative_max_bound_ratio']<1,
        frame_inverse=check['frame_inverse_max_error']<1e-75,
        frame_metric=check['frame_metric_max_error']<1e-75,
        balanced_poles=check['polar_identity_max_error']<1e-75,
        chart_sample_controls=check['chart_modulus_max_sample_ratio']<1,
        composition_control=check['composition_coefficient_error']<1e-75,
        finite_high_jets=all(len(row)==ORDER+1 and all(math.isfinite(v) and v>0 for v in row)
                             for row in jets['total_order_derivative_upper'].values()),
        finite_operator_jets=all(len(row)==19 and all(math.isfinite(v) and v>0 for v in row)
            for r in operators for row in r['coefficient_operator_norm_upper'].values()))
    result=dict(passed=all(checks.values()),checks=checks,controls=check,complex_cover=cover,
        step_derivative_upper=steps,regular_charts=chart,jets=jets,max_derivative_order=ORDER,
        operator_bounds=operators,
        source_sha256=sha(__file__),input_sha256=inputs,
        missing=dict(normal_coordinate_transport=None,heat_defect_B4=None,point_evaluation_constant=None,
                     curved_heat_remainder=None,target_source=None,stability=None,atom_transit=None),
        scope=['Directed high-derivative coefficient bounds in two overlapping regular local frames, over the actual buffered collars, including the joins and both polar axes.',
               'The collar is smooth but not analytic at its joins. Analytic Cauchy bounds are used for the exterior field functions; finite composition bounds are used for the sewn radial profile.',
               'Gauge coefficients are local, per abs(q); the nonzero-flux connection is not represented by one global potential.',
               'The coefficient bounds do not evaluate the normal-coordinate transport, heat defect, full supporting source, stability, atom passage or clock gain.'],
        sources=['https://dlmf.nist.gov/1.9.E31','https://arxiv.org/abs/1607.05152'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Regular high-jet checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,complex=cover,chart=chart,controls=check),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.5,4.5),layout='constrained')
    for name,label in (('coframe','Regular local frame'),('inverse_coframe','Inverse frame'),('gauge_per_abs_q','Local gauge potential / |q|')):
        ax.semilogy(range(ORDER+1),result['jets']['total_order_derivative_upper'][name],label=label)
    ax.set(xlabel='Total derivative order in (s, y1, y2), D = 1',ylabel='Uniform component-wise derivative ceiling',
           title='Actual collar derivatives bounded through order 20')
    ax.legend();ax.grid(alpha=.2)
    fig.supxlabel('Both polar charts and smooth joins included. Conservative bounds, not measured gradients.\nThe heat defect and supporting source remain uncomputed; atom passage is not established.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
