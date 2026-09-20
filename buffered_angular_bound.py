"""Explicit q=1 buffered angular/source bounds on saved exterior maps.

The continuum theorem requires true canonical-flux compressions of the
exterior operator. Saved finite-facet maps only illustrate its constants;
their discretization error and the frequency threshold remain unbounded.
No old solve or result is rerun. See RESEARCH increment104 for the proof.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from angular_defect_bound import defects, norm, response
from complete_angular_source import _decode, port_labels
from radial_angular_envelope import coefficients, envelope


def partner_bounds(a, d=1.):
    if not (math.isfinite(a) and math.isfinite(d) and 0<a<d/2):
        raise ValueError('Disjoint positive collars require 0<a<d/2')
    fmax=1/(4*(d-a)**2)
    axis=a*fmax*(1+a*(4*d-a)/(2*(d-a)**2))
    return dict(potential=a*fmax,rotation_axis=axis,rotation_vector=math.sqrt(2)*axis)


def graph_constant(kappa, U, V, a, b=.2):
    """Column-summed full-exterior angular graph bound, q=1 only.

    U=||X||_F, V=||D X||_F. Piecewise-linear nested cutoffs have
    chi2(b)=1 and h1=h2=2/(a-b). No derivative of the outer return is used.
    """
    if not all(math.isfinite(v) for v in (kappa,U,V,a,b)) or not (kappa>0 and U>=0 and V>=0 and 0<b<a<.5):
        raise ValueError('Positive frequency, nonnegative norms, .0<b<a<.5 required')
    p=partner_bounds(a); h=2/(a-b); c0=1-1/math.sqrt(2)
    u0=U/math.sqrt(2*kappa)
    energy=(kappa+p['potential']+h)**2*u0*u0+2*U*V/b
    u1=a*math.sqrt(energy/c0)
    result=math.sqrt(V*V+2*u1*(p['rotation_vector']*u0+h*u1))
    return dict(a=a,b=b,h1=h,h2=h,c0=c0,U=U,V=V,U0=u0,U1=u1,
                C_X=result,**p)


def source_bound(kappa, U, V, rank_columns, N, c_delta, a):
    if not isinstance(rank_columns,int) or rank_columns<0 or not isinstance(N,int) or N<0 or c_delta<=0:
        raise ValueError('Structural column count, angular cutoff and defect bound required')
    g=graph_constant(kappa,U,V,a)
    lam=math.sqrt((N+1)*(N+2))
    return dict(**g,rank_bound=rank_columns,lambda_first_omitted=lam,
                c_delta=c_delta,feedback_bound=math.sqrt(rank_columns*c_delta/lam)*g['C_X']/math.pi)


def rotation_checks():
    """Independent finite physical rotations test the Cartesian commutator."""
    e=np.array([0.,0.,1.]);worst=0.;rows=[]
    def field(x):
        D=np.linalg.norm(x+e)
        return np.cross(e,x)/(2*D*(D+1+x[2]))
    directions=np.random.default_rng(104).normal(size=(32,3))
    for a in (.2,.3,.4,.49):
        observed=0.;axial=0.
        for direction in directions:
            x=a*direction/np.linalg.norm(direction);D=np.linalg.norm(x+e)
            F=1/(2*D*(D+1+x[2]));Fz=-F*(1/D**2+(1/D+1)/(D+1+x[2]))
            for i,u in enumerate(np.eye(3)):
                w=np.cross(u,e)
                exact=F*np.cross(w,x)+Fz*np.dot(w,x)*np.cross(e,x)
                h=1e-5;forward=Rotation.from_rotvec(h*u);back=Rotation.from_rotvec(-h*u)
                finite=(forward.apply(field(back.apply(x)))-back.apply(field(forward.apply(x))))/(2*h)
                worst=max(worst,float(np.linalg.norm(finite-exact)))
                observed=max(observed,float(np.linalg.norm(exact)))
                if i==2:axial=max(axial,float(np.linalg.norm(exact)))
        rows.append(dict(a=a,**partner_bounds(a),sampled_rotation_max=observed,axial_commutator=axial))
    return dict(maximum_rotation_difference=worst,rows=rows)


def resolvent_checks():
    """Nonnormal returns and distinct port permutations test the ordering."""
    rng=np.random.default_rng(104);rows=[]
    for size,low,source in ((5,3,1),(9,5,3),(15,10,6)):
        z=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
        s=.97*z/norm(z)
        vals=np.linspace(.15,.95,size)*np.exp(1j*rng.uniform(-math.pi,math.pi,size))
        p=rng.permutation(size);q=rng.permutation(size)
        b=np.diag(vals)[np.ix_(p,q)];s=s[np.ix_(q,p)]
        il=np.flatnonzero(p<low);ol=np.flatnonzero(q<low)
        ik=np.flatnonzero(p<source);ok=np.flatnonzero(q<source)
        ih=np.flatnonzero(p>=low);oh=np.flatnonzero(q>=low)
        b0=np.zeros_like(b);b0[np.ix_(il,ol)]=b[np.ix_(il,ol)]
        f,f0=response(s,b),response(s,b0)
        w,v=defects(b)
        rhs=f[:,ih]@b[np.ix_(ih,oh)]@s[np.ix_(oh,il)]@np.linalg.inv(np.eye(low)-b[np.ix_(il,ol)]@s[np.ix_(ol,il)])
        ident=float(np.max(abs((f-f0)[:,il]-rhs)))
        off=norm(w[np.ix_(ok,ok)]@f[np.ix_(ok,ih)]@v[np.ix_(ih,ih)])
        rows.append(dict(size=size,identity_error=ident,unknown_weighted_offdiagonal_norm=off))
    return rows


def saved_source_inputs(prior, N, K, phase, probe):
    """Aggregate structural columns/Frobenius norms across both Weyls and m.

    Raw compressions are used because the theorem requires exact canonical
    projection. They are not substituted for the earlier finite-Gram N2 run.
    """
    group=next(g for g in prior['geometries'] if g['geometry']['rho_c']==32.)
    U2=V2=0.;rank=0;known=0.;sectors=[]
    for sector in group['sectors']:
        if abs(sector['m'])>K:continue
        sample=sector['samples'][1];mat=sample['matrices'];mu=sector['mu']
        li,lo=port_labels(mu,4);li,lo=li+lo,lo+li
        ni=np.array([n for _,n in li]);no=np.array([n for _,n in lo])
        il=np.flatnonzero(ni<=N);ol=np.flatnonzero(no<=N)
        ki=np.flatnonzero(ni[il]<=K);ko=np.flatnonzero(no[ol]<=K)
        full_s=_decode(mat['S_connected']);full_b=_decode(mat['phases'][phase]['B'])
        full_d=_decode(mat['phases'][phase]['pressure'])[probe]
        s=full_s[np.ix_(ol,il)];b=full_b[np.ix_(il,ol)];der=full_d[np.ix_(il,ol)]
        if norm(s)>1+2e-12:raise ValueError('Saved map is not passive in the assumed norm')
        w,_=defects(b);wk=w[np.ix_(ko,ko)]
        masked=np.zeros((len(il),len(ko)),complex)
        masked[ki,:]=der[np.ix_(ki,ko)]
        X=np.linalg.solve(np.eye(len(il))-b@s,masked@np.linalg.inv(wk))
        mass=np.sqrt(ni[il]*(ni[il]+1))
        U=float(np.linalg.norm(X));V=float(np.linalg.norm(mass[:,None]*X))
        U2+=U*U;V2+=V*V;rank+=len(ko)
        # Independent finite N4-versus-N K-supported feedback witness.
        b0=np.zeros_like(full_b);b0[np.ix_(il,ol)]=b
        dm=np.zeros_like(full_d)
        iK=np.flatnonzero(ni<=K);oK=np.flatnonzero(no<=K)
        dm[np.ix_(iK,oK)]=full_d[np.ix_(iK,oK)]
        feedback=float(np.trace((response(full_s,full_b)-response(full_s,b0))@dm).real/math.pi)
        known+=feedback
        sectors.append(dict(m=sector['m'],U=U,V=V,columns=len(ko),connected_norm=norm(s),finite_feedback=feedback))
    return dict(U=math.sqrt(U2),V=math.sqrt(V2),rank_columns=rank,finite_N4_minus_N_K_feedback=known,sectors=sectors)


def run():
    paths=('long_throat_source_results.json','angular_defect_bound.py','radial_angular_envelope.py',
           'complete_angular_source.py','complete_angular_modes.py','azimuthal_modes.py')
    hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
    prior=json.loads(Path(paths[0]).read_text(encoding='utf-8'))
    group=next(g for g in prior['geometries'] if g['geometry']['rho_c']==32.)
    x=group['sums']['4'][1]['x'];ell=group['geometry']['ell'];kappa=x/ell
    rotations=rotation_checks();identities=resolvent_checks();records=[]
    for N in (2,4):
        cd=coefficients(x,32.,0,cutoff=N)['delta_inverse_coefficient']
        for K in (0,1,2):
            for ip,phase in enumerate(('0','pi')):
                for probe in range(2):
                    inputs=saved_source_inputs(prior,N,K,ip,probe)
                    bounds=[source_bound(kappa,inputs['U'],inputs['V'],inputs['rank_columns'],N,cd,a) for a in (.3,.4,.49)]
                    best=min(bounds,key=lambda r:r['feedback_bound'])
                    direct=envelope(coefficients(x,32.,probe,cutoff=K))['direct_pressure_bound']
                    records.append(dict(N=N,K=K,phase=phase,probe=probe,**inputs,collars=bounds,
                        best_feedback_bound=best['feedback_bound'],best_a=best['a'],twice_direct_K_bound=2*direct,
                        conditional_full_source_bound=best['feedback_bound']+2*direct,
                        sampled_full_N4_kernel=group['sums']['4'][1]['kernels'][phase]['pressure'][probe]))
    fixed=[]
    cdelta=coefficients(x,32.,0,cutoff=4)['delta_inverse_coefficient']
    for N in (4,8,16,64,256,1024,1000000):
        g=graph_constant(kappa,1.,math.sqrt(6),.4)
        fixed.append(dict(N=N,weighted_mixing_bound=math.sqrt(cdelta/math.sqrt((N+1)*(N+2)))*g['C_X']))
    tests=dict(rotation_derivative=rotations['maximum_rotation_difference']<3e-10,
        commutator_samples_below_analytic_bound=all(r['sampled_rotation_max']<=r['rotation_axis'] and r['axial_commutator']==0 for r in rotations['rows']),
        distinct_port_resolvent_identity=max(r['identity_error'] for r in identities)<2e-12,
        offdiagonal_julia_contraction=max(r['unknown_weighted_offdiagonal_norm'] for r in identities)<=1+2e-12,
        structural_rank_and_complete_sectors=all(r['rank_columns']==2+4*r['K']*(r['K']+2) and len(r['sectors'])==2*r['K']+1 for r in records),
        known_finite_feedback_below_conditional_bound=all(abs(r['finite_N4_minus_N_K_feedback'])<=r['best_feedback_bound'] for r in records),
        fixed_input_tail_bound_decreases=all(a['weighted_mixing_bound']>b['weighted_mixing_bound'] for a,b in zip(fixed[:-1],fixed[1:])),
        inputs_preserved=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest for p,digest in hashes.items()))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,parameters=dict(q=1,b=.2,d=1.,rho_c=32.,x=x,kappa=kappa,ell=ell),
        rotation_controls=rotations,resolvent_controls=identities,records=records,fixed_input_cutoff_bounds=fixed,
        conclusions=dict(source_sign_certified=False),
        scope=['The explicit graph/source theorem is conditional on true continuum canonical-flux compressions and a well-posed passive exterior.',
               'Saved finite-facet matrices illustrate the constants, not their continuum values; raw canonical compression is distinct from prior finite-Gram reprojection.',
               'The fixed-input bound decreases with N, but replacing K by N does not prove convergence. The actual finite dressed X is included in each source bound.',
               'Two direct envelopes bound the full and compressed relative high-K traces. All multiplicities and both Weyls are already included once.',
               'These are pointwise positive-frequency estimates. They do not bound frequency tails, physical junctions, absolute smooth-collar stress or stability.',
               'Analytic inequalities are evaluated in floating point, not outward-rounded interval arithmetic.'],
        input_sha256=hashes,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=['https://arxiv.org/pdf/hep-th/0206054','https://arxiv.org/pdf/1803.09329'])
    if not result['passed']:raise AssertionError(tests)
    Path(__file__).with_name('buffered_angular_bound_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,phase_pi_K2=[r for r in records if r['phase']=='pi' and r['K']==2],fixed_input_cutoff_bounds=fixed),indent=2))
    return result


if __name__=='__main__':
    run()
