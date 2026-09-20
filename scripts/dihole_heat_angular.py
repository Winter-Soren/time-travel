"""Complete the angular sum of the frozen122 finite spatial model.

The tail certificate concerns the saved exact dyadic matrix polynomial only.
It does not enclose quadrature, continuum, boundary or finite eigensolve error.
"""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh
from mpmath.ctx_iv import MPIntervalContext

from dihole_finite_heat import physical_setup, form_matrices, point_evaluation, heat_sector
import dihole_interval_arithmetic as ia
from dihole_quadrature_remainder import dyadic, outward
from project_paths import project_path


def gershgorin_lower(matrix):
    lo,hi=matrix
    magnitude=np.maximum(abs(lo),abs(hi)).copy()
    np.fill_diagonal(magnitude,0.)
    row_sum=ia.sum_interval(ia.point(magnitude),axis=1)
    return float(np.min(ia.sub(ia.point(np.diag(lo)),row_sum)[0]))


def polynomial(matrices,j):
    return ia.add(ia.add(ia.point(matrices[0]),ia.mul(ia.point(float(j)),ia.point(matrices[1]))),
                  ia.mul(ia.point(float(j*j)),ia.point(matrices[2])))


def anchor_certificate(matrices,anchor,energy_lower,derivative_lower,quadratic_lower):
    """Outward congruence and Gershgorin proof; LAPACK supplies trial bases only."""
    assert all(np.array_equal(m,m.T) for m in matrices)
    a=gershgorin_lower(ia.point(matrices[2]))
    assert a>=quadratic_lower>0
    certificates=[]
    for sign in (-1,1):
        derivative=ia.add(ia.mul(ia.point(float(sign)),ia.point(matrices[1])),
                          ia.mul(ia.point(float(2*anchor)),ia.point(matrices[2])))
        d=gershgorin_lower(derivative)
        assert d>=derivative_lower>0
        numeric=matrices[0]+sign*anchor*matrices[1]+anchor**2*matrices[2]
        eigenvalues,Q=eigh(numeric,check_finite=False)
        gram=ia.matmul(ia.point(Q.T),ia.point(Q))
        gram_lower=gershgorin_lower(gram)
        shifted=ia.sub(polynomial(matrices,sign*anchor),ia.point(energy_lower*np.eye(len(Q))))
        transformed=ia.matmul(ia.point(Q.T),ia.matmul(shifted,ia.point(Q)))
        shifted_lower=gershgorin_lower(transformed)
        assert gram_lower>0 and shifted_lower>0,(gram_lower,shifted_lower)
        certificates.append(dict(sign=sign,numerical_anchor_eigenvalue=float(eigenvalues[0]),
            gram_gershgorin_lower=gram_lower,shifted_congruence_gershgorin_lower=shifted_lower,
            derivative_gershgorin_lower=d))
        print(json.dumps(dict(stage='anchor_certificate',sign=sign,shifted_lower=shifted_lower)),flush=True)
    return dict(anchor=anchor,energy_lower=energy_lower,derivative_lower=derivative_lower,
        quadratic_lower=quadratic_lower,quadratic_gershgorin_lower=a,signs=certificates)


def tail_bound(certificate,evaluation,times,normalization):
    """Both signs, including the unsummed anchors; exact dyadic input matrices."""
    iv=MPIntervalContext();iv.dps=75
    ell,d,a=(dyadic(iv,certificate[key]) for key in ('energy_lower','derivative_lower','quadratic_lower'))
    weights=ia.sum_interval(ia.square(ia.point(evaluation)),axis=(1,2))[1]
    bounds=[]
    for t,factors in zip(times,normalization):
        tau=dyadic(iv,t)
        series=2*iv.exp(-tau*ell)/(1-iv.exp(-tau*(d+a)))
        bounds.append([outward(series*dyadic(iv,w)*dyadic(iv,f))['upper']
                       for w,f in zip(weights,factors)])
    return np.array(bounds)


def controls():
    # Noncommuting exact-dyadic Gram family; checks congruence and both signs.
    X=np.array([[1.,.25],[0.,2.],[.5,0.]])
    Y=np.array([[.5,0.],[.25,1.],[0.,.5]])
    cross=X.T@Y
    matrices=(X.T@X,cross+cross.T,Y.T@Y)
    cert=anchor_certificate(matrices,16,20.,.5,.05)
    for j in (-31,-16,16,29):
        exact=(X+j*Y).T@(X+j*Y)
        interval=polynomial(matrices,j)
        assert np.all(interval[0]<=exact) and np.all(exact<=interval[1])
    evaluation=np.array([[[1.,.25],[.5,1.]]]);times=[.5,1.]
    bound=tail_bound(cert,evaluation,times,np.ones((2,1)))[:,0]
    direct=np.zeros(2)
    for j in list(range(-80,-15))+list(range(16,81)):
        values,vectors=eigh((X+j*Y).T@(X+j*Y))
        weights=np.sum((evaluation[0]@vectors)**2,axis=0)
        direct+=np.exp(-np.outer(times,values))@weights
    assert np.all(direct<=bound)
    # Interval Gershgorin cannot mistake a negative direction for positivity.
    assert gershgorin_lower(ia.point(np.diag([1.,-1.])))<0
    return dict(passed=True,noncommuting_direct_tail=direct.tolist(),tail_upper=bound.tolist())


def run():
    started=time.perf_counter();checks=controls()
    previous_path=project_path('dihole_finite_heat_results.json')
    assert hashlib.sha256(previous_path.read_bytes()).hexdigest()=='9b345db0f031f660bc46fbcc5329b513c32e93def1a0ee95298d11b366d7bf11'
    previous=json.loads(previous_path.read_text())
    assert hashlib.sha256(project_path('dihole_finite_heat.py').read_bytes()).hexdigest()==previous['source_sha256']
    for name,sha in previous['input_sha256'].items():
        assert hashlib.sha256(project_path(name).read_bytes()).hexdigest()==sha,name
    baseline=next(c for c in previous['cases'] if c['count']==18 and c['radius']==4.)
    coefficient,target=physical_setup()
    assert target==previous['target']
    matrices,quadrature=form_matrices(18,4.,60,coefficient)
    evaluation=point_evaluation([p['uv'] for p in target['points']],18,4.)
    times=np.array(previous['scaled_proper_times'])
    factors=((4*np.pi*times/target['energy']**2)**1.5)[:,None]*np.array(
        [2/(2*np.pi*target['length']**2*p['optical_volume_density']) for p in target['points']])
    old_heat=np.array(baseline['sector_heat'])
    reconstruction_errors=[]
    for j in (-384,0,384):
        got=heat_sector(matrices,j,evaluation,times)[0]
        error=float(np.max(abs(got-old_heat[j+384])))
        assert error<1e-12,error
        reconstruction_errors.append(error)
    cert=anchor_certificate(matrices,2048,160.,.05,.0003)
    tail=tail_bound(cert,evaluation,times,factors)
    offsets=np.array([s*j for j in range(385,2048) for s in (-1,1)])
    heat=[];retained=[];cutoff=40/min(times)
    # ponytail: dense finite spatial matrices; spatial continuum is a separate calculation.
    for index,j in enumerate(offsets):
        matrix=matrices[0]+j*matrices[1]+j*j*matrices[2]
        values,vectors=eigh(matrix,subset_by_value=(-np.inf,cutoff),driver='evr',check_finite=False)
        assert not len(values) or values[0]>-1e-8
        # Unlike the frozen122 helper, a sector with no eigenvalue below cutoff is valid.
        amplitude=np.sum((evaluation@vectors)**2,axis=1)
        heat.append(np.exp(-np.outer(times,values))@amplitude.T)
        retained.append(len(values))
        if index%256==0:
            print(json.dumps(dict(stage='added_angular_sectors',completed=index+1,total=len(offsets),
                seconds=time.perf_counter()-started)),flush=True)
    heat=np.array(heat)
    total_old=np.sum(old_heat,axis=0)*factors
    summaries=[]
    for cut in (384,512,768,1024,1536,2047):
        selected=abs(offsets)<=cut
        added=np.sum(heat[selected],axis=0)*factors
        final_shell=(abs(offsets)>cut-16)&selected
        summaries.append(dict(mode_cut=cut,normalized_heat=(total_old+added).tolist(),
            added_to_384=added.tolist(),final_16_offset_shell=(np.sum(heat[final_shell],axis=0)*factors).tolist()
            if cut>384 else None))
    eig_bound=(2*2047+1)*np.outer(np.exp(-times*cutoff),np.sum(evaluation*evaluation,axis=(1,2)))*factors
    state=project_path('dihole_heat_angular_state.npz')
    np.savez_compressed(state,M0=matrices[0],M1=matrices[1],M2=matrices[2],evaluation=evaluation,
        normalization=factors,offsets=offsets,sector_heat=heat,retained_eigenpairs=np.array(retained))
    names=['dihole_finite_heat.py','dihole_finite_heat_results.json','dihole_interval_arithmetic.py',
           'dihole_quadrature_remainder.py','dihole_heat_angular_state.npz']+list(previous['input_sha256'])
    result=dict(passed=True,controls=checks,target=target,count=18,radius=4.,quadrature=quadrature,
        scaled_heat_parameters=times.tolist(),reconstructed_baseline_heat_errors=reconstruction_errors,
        reused_angular_sectors=769,added_angular_sectors=len(offsets),summed_angular_sectors=4095,
        positive_spectrum_cutoff=cutoff,empty_retained_sectors=int(np.sum(np.array(retained)==0)),
        summaries=summaries,tail_certificate=cert,normalized_infinite_angular_tail_upper=tail.tolist(),
        normalized_finite_matrix_eigen_omission_bound=eig_bound.tolist(),
        scope=['The interval certificate covers |j|>=2048 for the saved exact dyadic matrix polynomial, evaluation and normalization.',
               'LAPACK supplies trial bases to the tail proof; positivity is verified independently by outward interval products.',
               'Finite angular contributions remain ordinary floating eigensolves; the whole displayed sum is not interval certified.',
               'Coefficient quadrature, spatial truncation and exterior boundary errors are not covered by the angular certificate.',
               'No renormalized stress, stable throat, supporting material or backward passage is established.'],
        arithmetic_assumptions='Inherited IEEE754 interval-ufunc and mpmath interval contract; exact saved dyadic inputs only.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in dict.fromkeys(names)},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_heat_angular_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps({key:result[key] for key in ('passed','runtime_seconds','summaries','normalized_infinite_angular_tail_upper')},indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,4.9))
    cuts=[r['mode_cut'] for r in result['summaries']]
    for j,label in enumerate(('At null','c = -1','c = +1')):
        left.plot(cuts,[r['normalized_heat'][0][j] for r in result['summaries']],'o-',label=label)
        shell=np.array([r['final_16_offset_shell'][0][j] for r in result['summaries'][1:]])
        positive=shell>0
        right.semilogy(np.array(cuts[1:])[positive],shell[positive],'o-',label=label)
    left.set(xlabel='Angular offset cutoff',ylabel='Normalized heat diagonal',title='Complete finite-model sum at heat parameter 0.5')
    right.set(xlabel='Angular offset cutoff',ylabel='Last 16 offsets on both sides',title='Measured angular shell contributions')
    for axis in (left,right):axis.grid(alpha=.2);axis.legend()
    fig.suptitle('Angular completion at fixed spatial basis and boundary',fontsize=14)
    right.text(.04,.06,'Zero retained contributions omitted;\neigenvalue omission is accounted separately.',transform=right.transAxes,fontsize=8)
    fig.text(.5,.02,'All offsets beyond 2047: certified contribution < 1.50e-33 for these matrices. Spatial/exterior errors and stress remain unresolved.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.94))
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_heat_angular_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
