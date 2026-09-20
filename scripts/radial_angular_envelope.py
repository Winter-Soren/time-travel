"""Analytic infinite-angular majorant for the DIRECT omitted source trace.

The estimate uses passive environments and the exact symmetric product core.
It does not bound feedback into retained channels, numerical error, or the
frequency integral. Nonnegative compact bumps are bounded by upper rectangles;
the secant-mass barrier is bounded from below by monotone Riccati comparison.
No continuum-error or interval-arithmetic certificate is claimed.
"""
from project_paths import project_path
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from angular_defect_bound import estimate
from complete_angular_source import _decode, port_labels
from throat_null_response import PROBES, bump


def y_lower(lam, k, s):
    omega = math.hypot(lam, k)
    t = math.tanh(omega*s)
    return lam*t/(omega+k*t)


def moments(alpha, first):
    """Exact geometric moments sum_(n>=first) (n+1/2)^(2,3) exp(-alpha*n)."""
    if not math.isfinite(alpha) or alpha <= 0 or first < 1:
        raise ValueError('Positive exponent and first massive level required')
    r, gap, p = math.exp(-alpha), -math.expm1(-alpha), first+.5
    s0, s1 = 1/gap, r/gap**2
    s2, s3 = r*(1+r)/gap**3, r*(1+4*r+r*r)/gap**4
    factor = math.exp(-alpha*first)
    return factor*(p*p*s0+2*p*s1+s2), factor*(p**3*s0+3*p*p*s1+3*p*s2+s3)


def coefficients(x, rho, probe, cutoff=4, bins=32):
    if not (math.isfinite(x) and x>0 and math.isfinite(rho) and rho>0):
        raise ValueError('Positive finite scaled frequency and throat cut required')
    if not isinstance(cutoff,int) or cutoff<0 or not isinstance(bins,int) or bins<1:
        raise ValueError('Integer nonnegative cutoff and positive bin count required')
    if probe not in range(len(PROBES)):
        raise ValueError('Unknown compact probe')
    c = math.atan(rho)
    length, first = 2*c, cutoff+1
    lam0 = math.sqrt(first*(first+1))
    primitive = lambda s: math.asinh(math.tan(s-c))
    y = lambda s: y_lower(lam0,x,s)
    grid = length*np.array([0.,1/8,1/4,1/2,3/4,7/8,1.])
    barrier = sum(y(float(a))*(primitive(float(b))-primitive(float(a))) for a,b in zip(grid[:-1],grid[1:]))
    ccut = math.hypot(1.,rho)
    delta_inverse_coefficient = ccut/(2*x*y(length/2)**2*(-math.expm1(-lam0*ccut*length)))
    spec = PROBES[probe]
    center, width = c*spec['center_fraction'], c*spec['width_fraction']
    edges = np.linspace(center-width,center+width,bins+1)
    rows = []
    for left,right in zip(edges[:-1],edges[1:]):
        a,b = float(left+c),float(right+c)
        closest = min(max(center,float(left)),float(right))
        weight = float(right-left)*bump(closest,c,probe)
        max_sec = 1/math.cos(max(abs(float(left)),abs(float(right))))
        rows.append(dict(a=a,b=b,weight=weight,C=max_sec,
            ar=2*y(a)*(primitive(length)-primitive(b)),
            al=2*y(length-b)*(primitive(a)-primitive(0.)),
            A0=2*x*(1+1/y(a)),A1=max_sec*(1+1/(2*y(a)))))
    return dict(x=x,rho_c=rho,probe=probe,cutoff=cutoff,first=first,bins=bins,
        length=length,ell=.2*ccut,lambda0=lam0,barrier=barrier,
        delta_inverse_coefficient=delta_inverse_coefficient,force_bins=rows)


def envelope(data):
    k,L,A = data['x'],data['length'],data['barrier']
    s2,s3 = moments(A,data['first'])
    reflected = transmitted = shift_null = 0.
    for row in data['force_bins']:
        f,C = row['weight'],row['C']
        for alpha,distance in ((row['al'],row['a']),(row['ar'],L-row['b'])):
            r2,r3 = moments(alpha,data['first'])
            reflected += f*math.exp(-2*k*distance)*(4*k*r2+C*r3)
        transmitted += 2*f*math.exp(-k*L)*(row['A0']*s2+row['A1']*s3)
        shift_null += 4*f*k*math.exp(-k*L)*s2
    factor = 8*data['delta_inverse_coefficient']/math.pi
    pressure = factor*(reflected+transmitted)
    return dict(direct_pressure_bound=pressure,direct_null_bound=pressure+factor*shift_null,
        scaled_pressure_bound=pressure/data['ell'],scaled_null_bound=(pressure+factor*shift_null)/data['ell'],
        reflection_part=factor*reflected,transmission_part=factor*transmitted,
        null_shift_addition=factor*shift_null)


def mode_weighted_bound(data,n):
    """Single paired mode ||dot B||_1/delta; no multiplicity included."""
    if n<data['first']:
        raise ValueError('Mode must lie in the declared omitted range')
    lam = math.sqrt(n*(n+1))
    k,L = data['x'],data['length']
    tau = math.exp(-k*L-lam*data['barrier'])
    value = 0.
    for row in data['force_bins']:
        reflection = math.exp(-2*k*row['a']-lam*row['al'])+math.exp(-2*k*(L-row['b'])-lam*row['ar'])
        value += row['weight']*((4*k+lam*row['C'])*reflection+2*tau*(row['A0']+row['A1']*lam))
    return value*data['delta_inverse_coefficient']*lam


def checks():
    series_error = 0.
    for alpha in (.025,.1,1.,5.,12.):
        for first in (3,5,9):
            n = np.arange(first,20001,dtype=float)
            exact = moments(alpha,first)
            for power,value in zip((2,3),exact):
                numerical = float(np.sum((n+.5)**power*np.exp(-alpha*n)))
                series_error = max(series_error,abs(numerical/value-1))
    comparisons = []
    for rho in (2.,32.):
        for x in (.044152948044253765,.44391506147705256,4.628335060152456):
            for probe in range(2):
                data = coefficients(x,rho,probe)
                bound = envelope(data)
                partial = 4/math.pi*sum((2*n+1)*mode_weighted_bound(data,n) for n in range(5,129))
                finer = envelope(coefficients(x,rho,probe,bins=64))
                c,L,lam = math.atan(rho),data['length'],data['lambda0']
                continuous_barrier = quad(lambda s:y_lower(lam,x,s)/math.cos(s-c),0,L,
                    epsabs=1e-12,epsrel=1e-12)[0]
                comparisons.append(dict(rho_c=rho,x=x,probe=probe,closed_bound=bound['direct_pressure_bound'],
                    finite_exact_lambda_majorant=partial,finer_bin_bound=finer['direct_pressure_bound'],
                    transmission_lower_barrier=data['barrier'],continuous_lower_barrier=continuous_barrier))
    tests=dict(geometric_series=series_error<1e-12,
        analytic_sum_dominates_finite_exact_lambda_sum=all(r['finite_exact_lambda_majorant']<=r['closed_bound']*(1+2e-13) for r in comparisons),
        nested_force_partition_tightens=all(r['finer_bin_bound']<=r['closed_bound']*(1+2e-13) for r in comparisons),
        monotone_mass_barrier=all(0<r['transmission_lower_barrier']<=r['continuous_lower_barrier']+2e-12 for r in comparisons))
    assert all(tests.values()),tests
    return tests,dict(maximum_geometric_series_relative_error=series_error,comparisons=comparisons)


def finite_feedback_diagnosis(prior):
    """Separate known N2->N4 direct trace, feedback and Gram projection terms."""
    answer=[]
    for group in prior['geometries'][1:]:
        for phase_index,phase in enumerate(('0','pi')):
            totals=np.zeros((3,2));rows=[]
            for sector in group['sectors']:
                sample=sector['samples'][1];mat=sample['matrices'];mu=sector['mu']
                labels_in,labels_out=port_labels(mu,4);count=len(labels_in)
                pin=[i for i,(_,n) in enumerate(labels_in) if n<=2]
                pout=[i for i,(_,n) in enumerate(labels_out) if n<=2]
                il=np.r_[pin,np.asarray(pout)+count];ol=np.r_[pout,np.asarray(pin)+count]
                record=estimate(_decode(mat['S_connected']),_decode(mat['S_reference']),
                    _decode(mat['phases'][phase_index]['B']),_decode(mat['phases'][phase_index]['pressure']),il,ol)
                contributions=[]
                for j,p in enumerate(record['probes']):
                    low=(sample['lower_cut_response']['kernels'][phase]['pressure'][j] if abs(sector['m'])<=2 else 0.)
                    contributions.append([p['direct_high'],p['low_feedback'],p['raw_compressed']-low])
                contribution=np.asarray(contributions).T;totals+=contribution
                rows.append(dict(m=sector['m'],direct_high=contribution[0].tolist(),
                    feedback=contribution[1].tolist(),gram=contribution[2].tolist()))
            expected=group['N4_minus_N2'][1]['kernels'][phase]['pressure']
            error=float(np.max(abs(totals.sum(axis=0)-expected)))
            assert error<2e-12
            answer.append(dict(rho_c=group['geometry']['rho_c'],phase=phase,x=group['sums']['4'][1]['x'],
                direct_high=totals[0].tolist(),feedback=totals[1].tolist(),gram=totals[2].tolist(),
                saved_cut_difference=expected,reconstruction_error=error,sectors=rows))
    return answer


def run():
    inputs=('throat_null_response.py','angular_defect_bound.py','complete_angular_source.py',
            'complete_angular_modes.py','azimuthal_modes.py','long_throat_source_results.json')
    hashes={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in inputs}
    tests,controls=checks()
    prior=json.loads(project_path('long_throat_source_results.json').read_text(encoding='utf-8'))
    feedback=finite_feedback_diagnosis(prior)
    tests['finite_feedback_decomposition']=max(r['reconstruction_error'] for r in feedback)<2e-12
    x_values=(1e-6,.044152948044253765,.44391506147705256,4.628335060152456,12.,24.,48.,96.)
    rows=[]
    for rho in (2.,8.,32.):
        for x in x_values:
            for probe in range(2):
                data=coefficients(x,rho,probe)
                rows.append(dict(x=x,rho_c=rho,probe=PROBES[probe]['name'],**envelope(data),
                    transmission_barrier=data['barrier'],delta_inverse_coefficient=data['delta_inverse_coefficient']))
    low=[]
    for group in prior['geometries']:
        geo=group['geometry'];epsilon=1e-6/geo['ell']
        low.append(dict(rho_c=geo['rho_c'],physical_epsilon=epsilon,
            zero_block_integral_bounds=[4*p['I_f']*epsilon/(math.pi*geo['optical_length']) for p in geo['probes']],
            scope='Pressure and each null zero-block trace, including mixing; massive blocks are not covered.'))
    result=dict(passed=all(tests.values()),tests=tests,controls=controls,records=rows,
        finite_feedback_diagnosis=feedback,low_frequency_zero_block=low,
        parameters=dict(cutoff=4,force_bins=32,mass_bins=[0,1/8,1/4,1/2,3/4,7/8,1],q=1,b=.2),
        input_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        meaning='Pointwise full-Dirac connected-minus-reference DIRECT omitted angular trace majorant; phase independent.',
        unresolved=['Retained-mode feedback through omitted channels','Low/high frequency integration',
                    'Numerical discretization errors','Physical supporting geometry'],
        numerical_scope='Analytic inequalities evaluated in floating point; no outward-rounded interval certificate.')
    project_path('radial_angular_envelope_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,middle=[r for r in rows if r['x']==x_values[2]]),indent=2))
    return result


if __name__=='__main__':
    run()
