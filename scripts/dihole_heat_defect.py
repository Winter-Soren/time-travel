"""An explicit Q4 heat-defect ceiling, not a supporting stress or atom transit.

All derivatives use normalized coordinate Frechet norms; all lengths use D=1.
The Gaussian's radial powers are retained before taking four-dimensional L2.
"""
import json
import math

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_local_transport import compose, dprod, packed, pmul
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_heat_defect'


def unpack(iv, value):
    return iv.mpf(int(value['mantissa'])) * iv.mpf(2)**value['binary_exponent']


def derivative(values, order=1):
    return [math.factorial(k+order)//math.factorial(k)*values[k+order]
            for k in range(len(values)-order)]


def physical_coefficients(iv, prior, record, prefix='P'):
    fields = prior['jets']['total_order_derivative_upper']
    v, vi = (list(map(iv.mpf, fields[k])) for k in ('half_density', 'inverse_half_density'))
    saved = record['coefficient_operator_norm_upper']
    p2, p1, p0 = (list(map(iv.mpf, saved[prefix+'_'+k])) for k in ('second', 'first', 'zero'))
    v1, v2 = dprod(v[1:], vi), dprod(v[2:], vi)
    physical = [p2, [a+8*b for a,b in zip(p1, dprod(p2,v1))],
                [a+4*b+16*c for a,b,c in zip(p0,dprod(p1,v1),dprod(p2,v2))]]
    return [[2**n*x/math.factorial(n) for n,x in enumerate(row)] for row in physical]


def apply_operator(coefficients, values, order):
    p2, p1, p0 = coefficients
    return [16*a+4*b+c for a,b,c in zip(pmul(p2, derivative(values,2), order),
            pmul(p1, derivative(values), order), pmul(p0, values, order))]


def graph_norm(coefficients, values):
    p2, p1, p0 = coefficients
    p0 = p0.copy(); p0[0] += 1
    for order in (4,2,0):
        values = apply_operator((p2,p1,p0), values, order)
    return values[0]


def distance_cutoff(iv, inverse, metric_upper, rho, step):
    # d^2 = <initial velocity, g(source) initial velocity>.
    # Keep the term linear in d: replacing it by a constant loses the t->0 bound.
    A = [iv.mpf(0)] + [metric_upper*sum(inverse[j]*inverse[n-j] for j in range(1,n))
                           for n in range(1,9)]
    B = [iv.mpf(0)] + [2*iv.sqrt(metric_upper)*inverse[n] for n in range(1,9)]
    inner, outer = rho/4, rho/2
    squared = [iv.mpf(0)] + [A[n]+B[n]*outer for n in range(1,9)]
    root = [iv.mpf(0)]
    coefficient = iv.mpf(1)
    for k in range(1,9):
        coefficient *= (iv.mpf(1)/2-(k-1))/k
        root.append(abs(coefficient)*inner**(1-2*k))
    distance = compose(root, squared, 8)
    profile = [iv.mpf(0)] + [iv.mpf(step[k])*(4/rho)**k/math.factorial(k) for k in range(1,9)]
    cutoff = compose(profile, distance, 8); cutoff[0] = iv.mpf(1)
    return A, B, cutoff


def gaussian_polynomials(iv, A, B, order=7):
    # Each key (a,b) multiplies d^a t^-b in the normalized Gaussian derivative.
    rows = [{(0,0): iv.mpf(1)}]
    for n in range(1,order+1):
        row = {}
        for k in range(1,n+1):
            for (a,b), value in rows[n-k].items():
                for shift, factor in ((0,A[k]), (1,B[k])):
                    if factor == 0:
                        continue
                    key = (a+shift,b+1)
                    row[key] = row.get(key,iv.mpf(0)) + value*factor*k/(4*n)
        if any(2*b-a > n for a,b in row):
            raise ArithmeticError('Gaussian radial scaling was lost.')
        rows.append(row)
    return rows


def time_supremum(iv, exponent, beta, cutoff):
    """sup_{0<t<=cutoff} t^exponent exp(-beta/t), with beta positive."""
    if beta.a <= 0 or cutoff.a <= 0:
        raise ValueError('Positive Gaussian gap and time cutoff required.')
    if exponent >= 0 or beta.a >= (-exponent*cutoff).b:
        return cutoff**exponent*iv.exp(-beta/cutoff)
    p = -exponent
    # The unrestricted positive-time maximum is also an upper bound if
    # interval uncertainty overlaps the endpoint/critical-point decision.
    return (p/beta)**p*iv.exp(-p)


def assemble(iv, prior, local, record, transport, point, heat_record, magnetic):
    normal, real = local['normal_geometry'], local['real_geometry']
    rho = iv.mpf(normal['normal_metric_radius_lower'])
    theta = iv.mpf(normal['orthonormal_volume_density_upper'])
    inverse = [iv.mpf(0)] + [unpack(iv,x) for x in normal['inverse_normalized_frechet_upper']]
    A,B,chi = distance_cutoff(iv,inverse,iv.mpf(real['metric_norm_upper']),rho,prior['step_derivative_upper'])
    gaussian = gaussian_polynomials(iv,A,B)
    p = physical_coefficients(iv,prior,record)
    uj = [[unpack(iv,x) for x in row] for row in transport['u_normalized_frechet_upper']]
    s0 = iv.mpf(heat_record['cutoff_s0_over_D_squared'])
    # ||d^a exp(-d^2/4t)||_L2 <= pi sqrt(theta*(a+1)!) (2t)^((a+2)/2).
    # Including t^4 (4pi t)^-2 gives exponent 3+a/2-b >=0 through order6.
    moments = []
    powers = set()
    for row in gaussian[:7]:
        value = iv.mpf(0)
        for (a,b), coefficient in row.items():
            twice_power = 6+a-2*b
            if twice_power < 0:
                raise ArithmeticError('Interior defect is not uniformly integrable.')
            powers.add(twice_power)
            value += coefficient*iv.pi*iv.sqrt(theta*math.factorial(a+1))*iv.mpf(2)**(iv.mpf(a+2)/2)*s0**(iv.mpf(twice_power)/2)
        moments.append(value/(4*iv.pi)**2)
    amplitude = pmul(chi,apply_operator(p,uj[4],6),6)
    interior_jets = pmul(amplitude,moments,6)
    interior = graph_norm(p,interior_jets)

    # Annular commutator: 16 p2 chi'' G u +32 p2 chi'(G'u+Gu')+4 p1 chi'Gu.
    beta, outer = rho**2/64, rho/2
    annulus_jets = [iv.mpf(0)]*7
    chi1, chi2 = derivative(chi), derivative(chi,2)
    h0 = [16*a+4*b for a,b in zip(pmul(p[0],chi2,6),pmul(p[1],chi1,6))]
    h1 = [32*x for x in pmul(p[0],chi1,6)]
    time_powers = set()
    for j, u in enumerate(uj):
        f0 = [a+b for a,b in zip(pmul(h0,u,6),pmul(h1,derivative(u),6))]
        f1 = pmul(h1,u,6)
        for n in range(7):
            for k in range(n+1):
                for shift, amplitude_value in ((0,f0[n-k]),(1,f1[n-k]*(k+1))):
                    for (a,b), coefficient in gaussian[k+shift].items():
                        exponent = j-2-b; time_powers.add(exponent)
                        annulus_jets[n] += amplitude_value*coefficient*outer**a*time_supremum(iv,exponent,beta,s0)/(4*iv.pi)**2
    # Physical 4-volume of the normal support, including the time-difference axis.
    volume = theta*iv.pi**2*outer**4/2
    annulus_jets = [iv.sqrt(volume)*x for x in annulus_jets]
    annulus = graph_norm(p,annulus_jets)
    column = interior+annulus
    # Canonical trace: 2 collars, each length1/20, two hemisphere disks of area pi.
    coordinate_source_volume = iv.pi/5
    source_factor = 4*coordinate_source_volume*iv.sqrt(iv.mpf(real['volume_upper']))
    Cev = unpack(iv,point['canonical_point_evaluation_upper'])
    B4 = source_factor*Cev*column
    # Diagonal Gaussian Hessian cancels for the exact volume-preserving probe.
    # These are Q4 coefficients; the O(t^2) defect alone does not identify
    # higher full heat-asymptotic coefficients. Lower coefficients stay exact.
    dp = physical_coefficients(iv,prior,record,prefix='deltaP')
    physical_source_volume = coordinate_source_volume*iv.mpf(real['volume_upper'])
    diagonal = [4*physical_source_volume*apply_operator(dp,u,0)[0] for u in uj[3:5]]
    Cstar = diagonal[0]+s0*diagonal[1]+16*iv.pi**2*B4
    Cgeometry = Cstar+iv.mpf(magnetic['magnetic_C_upper_D5'])
    matching = [dict(nominal_M_squared_s0=row['nominal_M_squared_s0'],
        remainder_envelope_per_q2_upper=packed(iv,Cgeometry*iv.mpf(row['C_geometry_error_coefficient_upper_D4'])/record['q']**2))
        for row in magnetic['matching_rows']]
    return dict(q=record['q'],cutoff_s0_over_D_squared=float(heat_record['cutoff_s0_over_D_squared']),
        gaussian_gap_over_s0_lower=lower(beta/s0),
        normal_support_volume_upper=packed(iv,volume),
        cutoff_normalized_frechet_upper=[packed(iv,x) for x in chi],
        squared_distance_constant_upper=[packed(iv,x) for x in A],
        squared_distance_linear_upper=[packed(iv,x) for x in B],
        interior_time_twice_powers=sorted(powers),annular_time_powers=sorted(time_powers),
        interior_L2_derivative_upper=[packed(iv,x) for x in interior_jets],
        annular_L2_derivative_upper=[packed(iv,x) for x in annulus_jets],
        physical_graph_interior_upper=packed(iv,interior),physical_graph_annulus_upper=packed(iv,annulus),
        physical_graph_column_upper=packed(iv,column),
        canonical_source_factor_upper=packed(iv,source_factor),
        heat_defect_interior_upper=packed(iv,source_factor*Cev*interior),
        heat_defect_annulus_upper=packed(iv,source_factor*Cev*annulus),
        heat_defect_B4_uniform_upper=packed(iv,B4),
        integrated_defect_coefficient_upper=packed(iv,16*iv.pi**2*B4),
        parametrix_b8_absolute_upper=packed(iv,diagonal[0]),
        parametrix_b10_absolute_upper=packed(iv,diagonal[1]),
        exact_low_order_remainder_Cstar_upper=packed(iv,Cstar),
        curved_remainder_Cgeometry_upper=packed(iv,Cgeometry),
        remainder_matching_diagnostics=matching,
        gaussian_polynomial_terms=[len(row) for row in gaussian],
        physical_source=None,atom_transit=None)


def controls(iv):
    mp.dps=80
    A=[iv.mpf(0),iv.mpf(0),iv.mpf(1)]+[iv.mpf(0)]*6
    B=[iv.mpf(0),iv.mpf(2)]+[iv.mpf(0)]*6
    g=gaussian_polynomials(iv,A,B)
    ratios=[]
    for n in range(8):
        for x in (mp.mpf('.2'),mp.mpf('.7')):
            t=mp.mpf('.3')
            actual=abs(mp.diff(lambda z:mp.exp(-z*z/(4*t)),x,n))/math.factorial(n)
            bound=sum(mp.mpf(c._mpi_[0])*x**a/t**b for (a,b),c in g[n].items())*mp.exp(-x*x/(4*t))
            ratios.append(actual/bound)
    beta,cutoff=iv.mpf(1)/50,iv.mpf(1)/2
    supremum_ratios=[]
    for power in range(-9,5):
        b=mp.mpf(1)/50;s=mp.mpf(1)/2
        critical=s if power>=0 else min(s,b/(-power))
        exact=critical**power*mp.exp(-b/critical)
        bound=time_supremum(iv,power,beta,cutoff)
        supremum_ratios.append(exact/mp.mpf(bound._mpi_[1]))
    # Constant-coefficient graph action on exp(k*x), including derivative factors.
    k=iv.mpf(3); f=[k**n/math.factorial(n) for n in range(7)]
    p=[[iv.mpf(c)]+[iv.mpf(0)]*6 for c in (2,3,5)]
    actual=graph_norm(p,f); expected=(1+16*2*k*k+4*3*k+5)**3
    graph_error=abs(actual-expected)
    # Independent radial Gaussian integral in four dimensions.
    t=mp.mpf('.2')
    integral=mp.sqrt(2*mp.pi**2*mp.quad(lambda r:r**3*mp.exp(-r*r/(2*t)),[0,1,mp.inf]))
    gaussian_integral_error=abs(integral-2*mp.pi*t)
    paired_gaussian_trace_error=abs((2*iv.mpf(3)/2)*2-(2*iv.mpf(3)/5)*5)
    return dict(gaussian_derivative_max_ratio=float(max(ratios)),
        time_supremum_max_ratio=float(max(supremum_ratios)),
        constant_graph_error=upper(graph_error),gaussian_integral_error=float(gaussian_integral_error),
        paired_gaussian_trace_error=upper(paired_gaussian_trace_error))


def run():
    iv=MPIntervalContext();iv.dps=90
    path=project_path('dihole_local_transport_results.json');local=json.loads(path.read_text())
    inputs={**local['input_sha256'],path.name:sha(path),'dihole_local_transport.py':local['source_sha256']}
    if not local['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    heat=json.loads(project_path('dihole_regulator_heat_tail_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    records=[assemble(iv,prior,local,record,transport,point,hr,mr) for record,transport,point,hr,mr in
             zip(prior['operator_bounds'],local['heat_transport'],local['point_evaluation'],heat['records'],magnetic['records'])]
    control=controls(iv)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        flux_records_aligned=all(a['q']==b['q']==c['q']==d['q']==e['q'] for a,b,c,d,e in
            zip(prior['operator_bounds'],local['heat_transport'],local['point_evaluation'],heat['records'],magnetic['records'])),
        two_flux_cases=len(records)==2,
        nonnegative_interior_time_powers=all(min(r['interior_time_twice_powers'])>=0 for r in records),
        annular_powers_covered=all(min(r['annular_time_powers'])>=-9 for r in records),
        atomic_endpoint_supremum=records[-1]['gaussian_gap_over_s0_lower']>9,
        finite_packed_defect=all(int(r['heat_defect_B4_uniform_upper']['mantissa'])>0 for r in records),
        gaussian_derivative_control=control['gaussian_derivative_max_ratio']<=1+1e-14,
        time_supremum_control=control['time_supremum_max_ratio']<=1+1e-14,
        graph_operator_control=control['constant_graph_error']<1e-70,
        gaussian_integral_control=control['gaussian_integral_error']<1e-70,
        paired_gaussian_cancellation=control['paired_gaussian_trace_error']<1e-70,
        finite_curved_remainder=all(int(r['curved_remainder_Cgeometry_upper']['mantissa'])>0 for r in records),
        physical_source_not_claimed=all(r['physical_source'] is None and r['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,
        convention='Uniform B4 ceiling for0<t<=s0, D=1. Four-dimensional physical column norm, canonical source trace per unit Killing time; rank4 and both collars included.',
        bound_encoding=local['bound_encoding'],
        missing=dict(signed_parametrix_b8=None,signed_parametrix_b10=None,
                     certified_numerical_lower_coefficients=None,finite_target_source=None,stability=None,atom_transit=None),
        scope=['The full local Q4 defect is bounded, including cutoff, Gaussian derivatives and physical half-density conversion.',
               'This is a conservative upper bound, not the true defect size or a signed supporting stress.',
               'Cstar subtracts exact lower coefficients; old numerical a2/a4/a6 estimates are not substituted into this certificate.',
               'Cgeometry uses the valid but loose triangle bound Cstar+Cmag. Matching diagnostics bound only this remainder, not the complete source.',
               'One hydrogen atom and a trial1nm areal opening remain the target; no atom transport or time shift is established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:
        raise ArithmeticError('Heat-defect checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,
        records=[dict(q=r['q'],interior_log10=r['physical_graph_interior_upper']['log10_upper'],
            annulus_log10=r['physical_graph_annulus_upper']['log10_upper'],
            B4_log10=r['heat_defect_B4_uniform_upper']['log10_upper']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row=result['records'][-1]
    keys=['heat_defect_interior_upper','heat_defect_annulus_upper','heat_defect_B4_uniform_upper']
    labels=['Interior contribution','Cutoff contribution','Full B4 bound']
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
    ax.bar(labels,[row[k]['log10_upper'] for k in keys],color=['#4769b2','#8b5897','#c17b37'])
    ax.set(ylabel='log10 of conservative ceiling, D = 1',title='The full heat-defect bound is finite, but extremely loose')
    fig.supxlabel('Atomic-flux reference. Worst-case error ceilings, not actual errors or energies.\nNo supporting stress, stable opening or atom passage is established.',fontsize=9)
    for ext in ('png','pdf'):
        fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
