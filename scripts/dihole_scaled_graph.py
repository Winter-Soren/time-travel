"""Tune analytic graph/Fourier scales without changing the physical model.

Every selected pair is an exact power of two and its bound is evaluated with
interval arithmetic. The exploratory search is not an optimality certificate.
"""
import json
import math

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_heat_defect import apply_operator, physical_coefficients, unpack
from dihole_local_transport import dirac_sequences, packed
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_scaled_graph'


def point_polynomials(iv, prior, record):
    """W_j(mu), represented as coefficients of mu^(p-6), p=0..6."""
    a,b,_ = dirac_sequences(iv,prior,record['q'])
    h=iv.mpf(2)**-39;L=iv.mpf(20)
    A=L*(b[0]+40*a[0]/h)+20/h
    c=[iv.mpf(1),25*iv.sqrt(5)/216,iv.mpf(4)/27,iv.mpf(1)/8,
       iv.mpf(4)/27,25*iv.sqrt(5)/216,iv.mpf(1)]
    U={}
    for p in range(7):
        row=[iv.mpf(0)]*7;row[p]=c[p];U[0,p]=row
    for ell in range(1,7):
        for p in range(7-ell):
            row=[L*x+A*y for x,y in zip(U[ell-1,p+1],U[ell-1,p])]
            for j in range(1,ell):
                row=[x+L*math.comb(ell-1,j)*(2**(j+1)*a[j]*y+2**j*b[j]*z)
                     for x,y,z in zip(row,U[ell-j,p],U[ell-1-j,p])]
            U[ell,p]=row
    W=[]
    for n in range(7):
        row=[iv.mpf(0)]*7
        for k in range(n+1):
            factor=math.comb(n,k)*2**k*math.factorial(k)*(20/h)**k
            row=[x+factor*y for x,y in zip(row,U[n-k,0])]
        W.append(row)
    saved=record['coefficient_operator_norm_upper']
    insertion=[iv.mpf(saved[k][0])/denominator for k,denominator in
               [('deltaP_second',12),('deltaP_first',30),('deltaP_zero',48)]]
    coefficients={}
    for n,weight in ((0,1),(2,3),(4,3),(6,1)):
        for p,value in enumerate(W[n]):
            if value==0:continue
            for k,delta in enumerate(insertion):
                key=(p-6,4-n-k)
                coefficients[key]=coefficients.get(key,iv.mpf(0))+weight*value*delta
    return coefficients,c,W


def evaluate(iv, coefficients, mu_exponent, fourier_exponent):
    return sum(value*iv.mpf(2)**(p*mu_exponent+k*fourier_exponent)
               for (p,k),value in coefficients.items())


def choose_scales(iv, coefficients):
    """Small finite search; the final interval evaluation supplies validity."""
    logs=[(p,k,upper(iv.log(value,2))) for (p,k),value in coefficients.items()]
    def score(a,b):
        values=[c+p*a+k*b for p,k,c in logs];peak=max(values)
        return peak+math.log2(sum(2**(x-peak) for x in values))
    candidates=[(score(a,b),a,b) for a in range(-32,225,8) for b in range(-32,225,8)]
    _,a,b=min(candidates)
    candidates.extend((score(x,y),x,y) for x in range(a-8,a+9) for y in range(b-8,b+9))
    value,a,b=min(candidates)
    return a,b,dict(coarse_exponent_limits=[-32,224],coarse_step=8,fine_radius=8,
        selected_mu_exponent=a,selected_fourier_exponent=b,approximate_log2_bound=value,
        optimum_claimed=False)


def record_bound(iv, prior, operator, defect, magnetic):
    point,spectral,W=point_polynomials(iv,prior,operator)
    coefficients=physical_coefficients(iv,prior,operator)
    values=[unpack(iv,a)+unpack(iv,b) for a,b in
            zip(defect['interior_L2_derivative_upper'],defect['annular_L2_derivative_upper'])]
    graph=[values[0]]
    for n in (4,2,0):
        values=apply_operator(coefficients,values,n);graph.append(values[0])
    source_factor=unpack(iv,defect['canonical_source_factor_upper'])
    combined={}
    for (p,k),value in point.items():
        for j,g in enumerate(graph):
            key=(p+6-2*j,k)
            combined[key]=combined.get(key,iv.mpf(0))+source_factor*value*math.comb(3,j)*g
    a,b,search=choose_scales(iv,combined)
    B4=evaluate(iv,combined,a,b)
    B4_packed=packed(iv,B4);B4_upper=unpack(iv,B4_packed)
    old=unpack(iv,defect['heat_defect_B4_uniform_upper'])
    unscaled=evaluate(iv,combined,0,0)
    Cev=evaluate(iv,point,a,b)
    graph_selected=sum(math.comb(3,j)*iv.mpf(2)**((6-2*j)*a)*g for j,g in enumerate(graph))
    diagonal=unpack(iv,defect['parametrix_b8_absolute_upper'])+iv.mpf(defect['cutoff_s0_over_D_squared'])*unpack(iv,defect['parametrix_b10_absolute_upper'])
    Cstar=diagonal+16*iv.pi**2*B4
    Cgeometry=Cstar+iv.mpf(magnetic['magnetic_C_upper_D5'])
    matching=[dict(nominal_M_squared_s0=row['nominal_M_squared_s0'],
        remainder_envelope_per_q2_upper=packed(iv,Cgeometry*iv.mpf(row['C_geometry_error_coefficient_upper_D4'])/defect['q']**2))
        for row in magnetic['matching_rows']]
    return dict(q=defect['q'],mu_exact_power_of_two=a,fourier_exact_power_of_two=b,
        selected_Cev_upper=packed(iv,Cev),selected_graph_column_upper=packed(iv,graph_selected),
        graph_power_column_upper=[packed(iv,x) for x in graph],
        spectral_constants_upper=[packed(iv,x) for x in spectral],
        point_polynomial=[dict(mu_power=p,fourier_power=k,coefficient_upper=packed(iv,v)) for (p,k),v in sorted(point.items())],
        defect_polynomial=[dict(mu_power=p,fourier_power=k,coefficient_upper=packed(iv,v)) for (p,k),v in sorted(combined.items())],
        unscaled_refined_B4_upper=packed(iv,unscaled),heat_defect_B4_uniform_upper=B4_packed,
        old_to_new_ceiling_log10_ratio_lower=lower(iv.log(old/B4_upper,10)),
        exact_low_order_remainder_Cstar_upper=packed(iv,Cstar),
        curved_remainder_Cgeometry_upper=packed(iv,Cgeometry),
        remainder_matching_diagnostics=matching,
        source_factor_upper=packed(iv,source_factor),search=search,
        selected_factorization_ratio=upper(source_factor*Cev*graph_selected/B4),
        improvements_verified=bool(B4.b<unscaled.a and unscaled.b<old.a),
        physical_source=None,atom_transit=None)


def controls(iv):
    mp.dps=85
    c=[mp.mpf(1),25*mp.sqrt(5)/216,mp.mpf(4)/27,mp.mpf(1)/8,
       mp.mpf(4)/27,25*mp.sqrt(5)/216,mp.mpf(1)]
    critical_errors=[]
    for p in range(1,6):
        x=mp.mpf(p)/(6-p)
        critical_errors.append(abs(x**(mp.mpf(p)/2)/(1+x)**3-c[p]))
    # Independent radial quadrature at several Fourier scales in dimension4.
    fourier_ratios=[]
    for lam in (mp.mpf('.2'),mp.mpf(1),mp.mpf(7)):
        for k in range(3):
            actual=mp.sqrt(mp.quad(lambda r:r**(3+2*k)/(lam*lam+r*r)**6,[0,lam,mp.inf])/(8*mp.pi**2))
            exact=lam**(k-4)/(4*mp.pi*mp.sqrt(30 if k==1 else 20))
            fourier_ratios.append(abs(actual/exact-1))
    # Check the spectral multiplier inequality on nonnegative eigenvalues.
    graph_ratio=mp.mpf(0)
    for mu in (mp.mpf('.2'),mp.mpf(1),mp.mpf(5)):
        for eigenvalue in (mp.mpf(0),mp.mpf('.1'),mp.mpf(3),mp.mpf(90)):
            for p in range(7):
                numerator=eigenvalue**(mp.mpf(p)/2) if p else mp.mpf(1)
                graph_ratio=max(graph_ratio,numerator/(mu*mu+eigenvalue)**3/(c[p]*mu**(p-6)))
    toy={(-2,0):iv.mpf(1),(2,0):iv.mpf(1),(0,-2):iv.mpf(1),(0,2):iv.mpf(1)}
    a,b,_=choose_scales(iv,toy)
    return dict(spectral_critical_error=float(max(critical_errors)),fourier_scaling_error=float(max(fourier_ratios)),
                spectral_grid_max_ratio=float(graph_ratio),known_search_minimum=(a==0 and b==0))


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_heat_defect_results.json');defect=json.loads(path.read_text())
    inputs={**defect['input_sha256'],path.name:sha(path),'dihole_heat_defect.py':defect['source_sha256']}
    if not defect['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing inputs required.')
    prior=json.loads(project_path('dihole_regular_jets_results.json').read_text())
    magnetic=json.loads(project_path('dihole_magnetic_remainder_results.json').read_text())
    records=[record_bound(iv,prior,operator,row,mr) for operator,row,mr in
             zip(prior['operator_bounds'],defect['records'],magnetic['records'])]
    control=controls(iv)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        aligned_flux=len(records)==2 and all(r['q']==p['q']==m['q'] for r,p,m in
            zip(records,prior['operator_bounds'],magnetic['records'])),
        both_certified_ceilings_improve=all(r['improvements_verified'] for r in records),
        dyadic_scales=all(isinstance(r[k],int) for r in records for k in ('mu_exact_power_of_two','fourier_exact_power_of_two')),
        positive_polynomial=all(int(c['coefficient_upper']['mantissa'])>0 for r in records for c in r['defect_polynomial']),
        factorization_control=all(abs(r['selected_factorization_ratio']-1)<1e-12 for r in records),
        spectral_extremum_control=control['spectral_critical_error']<1e-70,
        fourier_scaling_control=control['fourier_scaling_error']<1e-70,
        spectral_grid_control=control['spectral_grid_max_ratio']<=1+1e-14,
        search_control=control['known_search_minimum'],
        search_not_optimality_claim=all(r['search']['optimum_claimed'] is False for r in records),
        physical_source_not_claimed=all(r['physical_source'] is None and r['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        source_sha256=sha(__file__),input_sha256=inputs,bound_encoding=defect['bound_encoding'],
        convention='D=1. mu is an auxiliary graph-norm scale and lambda is an auxiliary Fourier evaluation scale; neither is the physical regulator, a field parameter or an opening size.',
        scope=['The full181 defect, geometry, physical field, normal neighbourhood, time interval and quantum payload are unchanged.',
               'The positive Laurent polynomial certifies every selected dyadic pair. The finite numerical search does not certify the global minimum.',
               'The improvement is a ratio of conservative ceilings, not a measured physical improvement or an error reduction of known size.',
               'The matching diagnostics bound only the curved remainder. Neither auxiliary scale changes the saved target regulator M.',
               'Cstar subtracts exact lower coefficients; the full supporting source, stability, atomic passage and clock gain remain unproved.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Scaled graph checks failed; evidence retained.')
    plot(result,defect)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,records=[dict(q=r['q'],mu_power=r['mu_exact_power_of_two'],
        fourier_power=r['fourier_exact_power_of_two'],B4_log10=r['heat_defect_B4_uniform_upper']['log10_upper'],
        improvement_log10=r['old_to_new_ceiling_log10_ratio_lower'],Cstar_log10=r['exact_low_order_remainder_Cstar_upper']['log10_upper']) for r in records]),indent=2))


def plot(result,prior):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    row=result['records'][-1];old=prior['records'][-1]
    ax.bar(['Previous ceiling','Adjusted proof scales'],[old['heat_defect_B4_uniform_upper']['log10_upper'],
        row['heat_defect_B4_uniform_upper']['log10_upper']],color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of B4 ceiling, D = 1',title='A tighter error bound for the same proposed opening')
    fig.supxlabel('Atomic-flux reference. Bounds remain too loose to establish physical support.\nProof scales changed; the field, opening and atom did not.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
