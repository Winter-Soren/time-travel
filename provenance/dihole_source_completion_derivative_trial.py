"""Exact full-source completion and a global lapse-weight commutator bound.

The full weighted smoothing traces are NOT evaluated. Finite matrices check
the identities; the geometric bounds concern the actual prescribed metric.
"""
import json
import math

import numpy as np
import sympy as sp
from scipy.integrate import quad

from balanced_dihole import symbolic_fields
from dihole_channel_bounds import MPIntervalContext, calibration, lower, upper
from dihole_slice_channel import sha
from dirac_green_mass_limit import benchmark
from project_paths import project_path


NAME='dihole_source_completion'


def lapse_expression(context):
    R,c,m,a,b=sp.symbols('R c m a b',real=True)
    r=R+m;s=1-c*c;delta=R*R-m*m-a*a;sigma=r*r-a*a*c*c
    w=(r*r-a*a)**2+delta*a*a*s
    N=(delta+a*a*s+2*b*m*r*a*s+b*b*w*s/4)/sigma
    fn=sp.lambdify((R,c,m,a,b),(N,sp.diff(N,R),sp.diff(N,c)),
                   modules=[{'sqrt':context.sqrt}],cse=True)
    return fn


def weight_bound(iv, floor):
    k,m,a,B,_=calibration(iv);cn=iv.mpf(floor);epsilon=iv.mpf(1)/40
    h0=epsilon*(2*k-epsilon);fn=lapse_expression(iv)
    maxima=[0.,0.]
    for i in range(16):
        R=iv.mpf([str(.5+i/32),str(.5+(i+1)/32)])
        for j in range(16):
            c=iv.mpf([str(-1+j/8),str(-1+(j+1)/8)])
            _,Nr,Nc=fn(R,c,m,a,B)
            for n,value in enumerate((Nr,Nc)):
                maxima[n]=max(maxima[n],abs(lower(value)),abs(upper(value)))
    nr,nc=map(iv.mpf,maxima)
    near=iv.sqrt(nr*nr+nc*nc/h0)/(2*cn**iv.mpf('1.5'))
    e=1-k*k;b=B*B/4
    U0=2*m*(1+m)/e
    U1=2*m*(1+2*m+k*k)/e**2
    T=2*a*a/(e*iv.sqrt(e))+2*a*a/e**2
    CL=U1+U0*T+2*B*a*U0/iv.sqrt(e)
    CG=(1+2*U0)*U1+U0*(1+U0)*T
    A0=CL+b*CG;A1=2*b*(1+U0)
    rho_star=iv.sqrt(cn/b)
    far=(A0+A1*rho_star)/(2*cn**iv.mpf('1.5'))
    return dict(prolate_split_R_over_D=1.,near_boxes=256,
                near_N_R_upper_D=maxima[0],near_N_c_upper=maxima[1],
                minimum_h_over_D2=lower(h0),far_U0_upper=upper(U0),far_U1_upper=upper(U1),
                far_t_gradient_coefficient_upper=upper(T),far_L_gradient_coefficient_upper=upper(CL),
                far_G_gradient_coefficient_upper=upper(CG),far_A0_upper_D=upper(A0),
                far_A1_upper_D2=upper(A1),far_rho_maximizer_upper_over_D=upper(rho_star),
                near_commutator_upper_D=upper(near),far_commutator_upper_D=upper(far),
                global_commutator_upper_D=max(upper(near),upper(far)))


def geometry_controls(bound):
    iv=MPIntervalContext();iv.dps=60
    k,m,a,B,_=[float(x.mid) for x in calibration(iv)]
    actual=symbolic_fields()[0];direct=lapse_expression(np)
    errors=[];rows=[]
    for R in (.5001,.51,.7,1.,2.,10.,100.,10000.):
        for c in (-.99,-.8,0.,.8,.99):
            if R-k*abs(c)<.025:continue
            N,Nr,Nc=map(float,direct(R,c,m,a,B))
            supplied=np.asarray(actual(R+m,c,m,a,B),float)
            errors.append(float(max(abs(N-supplied[0]),abs(Nr-supplied[3]),abs(Nc-supplied[6]))/max(1.,abs(N),abs(Nr),abs(Nc))))
            delta=R*R-k*k;s=1-c*c;h=delta+k*k*s;sigma=(R+m)**2-a*a*c*c
            rho=math.sqrt(delta*s)
            gradient=math.sqrt(delta/h*Nr*Nr+s/h*Nc*Nc)
            comm=(h/sigma)**2*gradient/(2*N**1.5)
            far_enclosed=None
            if R>=1:far_enclosed=bool(gradient<=bound['far_A0_upper_D']+bound['far_A1_upper_D2']*rho)
            rows.append(dict(R_over_D=R,c=c,commutator_norm_D=comm,
                             enclosed=bool(comm<=bound['global_commutator_upper_D']),
                             far_gradient_enclosed=far_enclosed))
    # Orthogonal prolate scale factors and the radius contraction are analytic,
    # not inferred from the sampled points.
    u,v,w,dw=sp.symbols('u v w dw',real=True)
    psi=sp.Function('psi')(u);F=sp.Function('F')(u);C=sp.Function('C')(u)
    action=lambda x:C*sp.diff(x,u)+sp.diff(C,u)*x/2
    paired=(F*action(psi)+action(F*psi))/2
    factor=F*C*sp.diff(psi,u)+sp.diff(F*C,u)*psi/2
    return dict(lapse_jet_relative_errors=errors,sampled_global_points=rows,
                radial_form_factorization=sp.simplify(paired-factor)==0)


def positive_root(Q):
    values,vectors=np.linalg.eigh((Q+Q.conj().T)/2)
    if min(values)<-1e-12:raise ArithmeticError('Covariance is not positive.')
    return (vectors*np.sqrt(np.maximum(values,0)))@vectors.conj().T


def matrix_controls():
    rng=np.random.default_rng(174)
    V=np.linalg.qr(rng.normal(size=(6,3))+1j*rng.normal(size=(6,3)))[0]
    window=np.diag([.35,1.,.7]);h=np.diag([-1.3,.2,2.])
    eye=np.eye(6);small=np.eye(3);Pi=V@V.conj().T
    norm=lambda x:float(np.linalg.norm(x,2))
    trace_norm=lambda x:float(np.linalg.svd(x,compute_uv=False).sum())
    rows=[];weights=[];partition_rows=[]
    for largest in (2.,10.,100.):
        H,beta,N,dH,dN,A0,_,_=benchmark(largest,0.)
        sqrtN=np.diag(np.sqrt(np.diag(N)));invroot=np.diag(1/np.sqrt(np.diag(N)))
        lapseF=np.linalg.solve(N,dN);D=beta@invroot@dH@invroot-(A0@lapseF-lapseF@A0)/2
        ca=norm(lapseF);da=norm(D);dot_slope=-1j*lapseF@beta@np.linalg.inv(N)
        def operators(kappa,M,eta):
            _,_,_,_,_,A,P,dP=benchmark(largest,kappa)
            dot=D-lapseF@A
            masses=np.sqrt(eta+np.arange(4)*M*M);alpha=math.sqrt(6)*M**3
            R=[np.linalg.inv(A+m*eye) for m in masses]
            K=alpha*R[1]@R[2]@R[3];F=R[0]@K
            B=beta@sqrtN@K.conj().T@dot@K@sqrtN
            return A,P,dP,dot,masses,R,K,F,B
        for kappa in (-.3,.3,1.):
            for M in (.7,2.):
                for eta_ratio in (0.,.01):
                    eta=eta_ratio*M*M
                    A,P,dP,dot,masses,R,K,F,B=operators(kappa,M,eta)
                    g=np.linalg.inv(H+masses[0]*N@beta+1j*kappa*eye)
                    Bsrc=-F.conj().T@dP@F/2
                    full=float(np.trace(Bsrc).real)
                    cyclic=-float(np.trace(g@B).real)
                    z=K@sqrtN
                    Bbound=norm(dot)*np.linalg.norm(z,'fro')**2
                    generator=1j*beta@np.linalg.inv(N)
                    Kprime=-math.sqrt(6)*M**3*sum(
                        math.prod([])*np.linalg.multi_dot(
                            [*R[1:j+1],R[j],generator,R[j],*R[j+1:]])
                        for j in (1,2,3))
                    # Prefix above excludes R[j]; append it on both sides of generator.
                    Bprime=beta@sqrtN@(Kprime.conj().T@dot@K+
                        K.conj().T@dot_slope@K+K.conj().T@dot@Kprime)@sqrtN
                    zp=Kprime@sqrtN
                    derivative_bound=norm(dot_slope)*np.linalg.norm(z,'fro')**2+2*norm(dot)*np.linalg.norm(z,'fro')*np.linalg.norm(zp,'fro')
                    step=1e-5
                    finite=(operators(kappa+step,M,eta)[-1]-operators(kappa-step,M,eta)[-1])/(2*step)
                    for lam in (.5,2.):
                        E=np.diag(np.exp(-np.diag(h)**2/lam**2))
                        Y=V@window@np.diag(np.exp(-np.diag(h)**2/(2*lam**2)))
                        Q=Y@Y.conj().T
                        parts=[eye-Pi,V@(small-window@window)@V.conj().T,
                               V@window@(small-E)@window@V.conj().T]
                        remainder=eye-Q;root=positive_root(remainder)
                        e0=float(np.linalg.norm(F@root,'fro')**2)
                        e1=float(np.linalg.norm(A@F@root,'fro')**2)
                        exact0=float(np.trace(F@F.conj().T).real-np.linalg.norm(F@Y,'fro')**2)
                        exact1=float(np.trace(P@F@F.conj().T).real-np.linalg.norm(A@F@Y,'fro')**2)
                        selected=float(np.trace(Q@Bsrc).real)
                        missing=float(np.trace(remainder@Bsrc).real)
                        direct=-float(np.vdot(A@F@root,dot@F@root).real)
                        wrong=-float(np.trace(Q@dP@F@F.conj().T).real)/2
                        rows.append(dict(largest_lapse=largest,kappa=kappa,M=M,eta_over_M_squared=eta_ratio,Lambda=lam,
                            positive_partition_min_eigenvalue=min(float(np.linalg.eigvalsh(p).min()) for p in parts),
                            partition_error=norm(sum(parts)-remainder),source_additivity_error=abs(full-selected-missing),
                            omitted_form_relative_error=abs(missing-direct)/max(1,abs(missing)),
                            omitted_gram_relative_error=max(abs(e0-exact0)/max(1,e0),abs(e1-exact1)/max(1,e1)),
                            remainder_bound=ca*e1+da*math.sqrt(e0*e1),remainder_value=missing,
                            wrong_projector_placement_difference=abs(selected-wrong),
                            unsplit_full_source_relative_error=abs(full-cyclic)/max(1,abs(full)),
                            full_trace_norm_bound_ratio=trace_norm(B)/Bbound,
                            derivative_trace_norm_bound_ratio=trace_norm(Bprime)/derivative_bound,
                            derivative_relative_error=norm(Bprime-finite)/max(1,norm(Bprime))))
                        if largest==10 and kappa==.3 and M==.7 and eta_ratio==.01:
                            covariances=[F@Q@F.conj().T]+[F@p@F.conj().T for p in parts]
                            total=float(np.trace(F@F.conj().T).real)
                            partition_rows.append(dict(Lambda=lam,fractions=[float(np.trace(c).real)/total for c in covariances]))
                    mass=1.;res=np.linalg.inv(A+mass*eye);chi=np.diag([1.,1.,0.,1.,1.,0.])
                    C=sqrtN@A-A@sqrtN
                    identity=res@sqrtN-sqrtN@res-res@C@res
                    weights.append(dict(largest_lapse=largest,kappa=kappa,
                        weighted_resolvent_identity_error=norm(identity),
                        compact_weighted_resolvent_norm=norm(chi@res@sqrtN),
                        bound=math.sqrt(2)/mass+norm(C)/mass**2))
    return rows,weights,partition_rows


def heat_controls():
    rows=[]
    for x in (.003,.4,2.):
        eta=.02;M=.7
        exact=6*M**6/math.prod(x+eta+j*M*M for j in range(4))
        value,error=quad(lambda t:math.exp(-(x+eta)*t)*(-math.expm1(-M*M*t))**3,0,np.inf,epsabs=1e-11,epsrel=1e-11)
        rows.append(dict(eigenvalue=x,rational=exact,heat_integral=value,relative_error=abs(value-exact)/exact,quadrature_error=error))
    return rows


def run():
    path=project_path('dihole_paired_low_frequency_results.json');prior=json.loads(path.read_text())
    inputs=dict(prior['input_sha256'])
    inputs.update({path.name:sha(path),'dihole_paired_low_frequency.py':prior['source_sha256']})
    if not prior['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing prior inputs required.')
    iv=MPIntervalContext();iv.dps=60
    geometry=weight_bound(iv,prior['global_lapse']['lapse_lower'])
    controls=geometry_controls(geometry);matrices,weights,parts=matrix_controls();heat=heat_controls()
    lapse=json.loads(project_path('dihole_massive_channel_results.json').read_text())['records'][0]['lapse']
    channel=json.loads(project_path('dihole_channel_bounds_results.json').read_text())
    n=iv.mpf(lapse['lapse_lower']);w=iv.mpf(channel['uniform_bound']['constants']['W_lower'][0]);length=iv.mpf(1)/100
    insertion=dict(buffer_support_over_D=[-.005,.055],buffer_one_on_over_D=[0.,.05],
                   auxiliary_length_over_D=.01,J_zero_upper_D=upper(4/(2*length*n*w)),
                   J_frequency_coefficient_upper=upper(4/n))
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        original_geometry_jets=max(controls['lapse_jet_relative_errors'])<1e-10,
        global_weight_samples=all(r['enclosed'] for r in controls['sampled_global_points']),
        far_gradient_samples=all(r['far_gradient_enclosed'] is not False for r in controls['sampled_global_points']),
        radial_form_factorization=controls['radial_form_factorization'],
        positive_completion=min(r['positive_partition_min_eigenvalue'] for r in matrices)>-1e-12,
        exact_completion=max(r['partition_error'] for r in matrices)<1e-12,
        source_additivity=max(r['source_additivity_error'] for r in matrices)<1e-9,
        omitted_form=max(r['omitted_form_relative_error'] for r in matrices)<1e-9,
        omitted_full_gram=max(r['omitted_gram_relative_error'] for r in matrices)<1e-9,
        omitted_source_bound=all(abs(r['remainder_value'])<=r['remainder_bound']*(1+1e-12) for r in matrices),
        wrong_placement_detected=max(r['wrong_projector_placement_difference'] for r in matrices)>1e-4,
        unsplit_full_source=max(r['unsplit_full_source_relative_error'] for r in matrices)<1e-9,
        weighted_factor_trace_bound=max(r['full_trace_norm_bound_ratio'] for r in matrices)<=1+1e-12,
        weighted_derivative_trace_bound=max(r['derivative_trace_norm_bound_ratio'] for r in matrices)<=1+1e-12,
        differentiated_full_source=max(r['derivative_relative_error'] for r in matrices)<2e-7,
        weighted_resolvent_identity=max(r['weighted_resolvent_identity_error'] for r in weights)<1e-10,
        compact_weighted_resolvent=all(r['compact_weighted_resolvent_norm']<=r['bound']*(1+1e-12) for r in weights),
        common_positive_heat=max(r['relative_error'] for r in heat)<1e-10,
        covariance_partition=sum(abs(sum(r['fractions'])-1) for r in parts)<1e-12)
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,global_weight=geometry,insertion_factor=insertion,
                geometry_controls=controls,matrix_controls=matrices,weighted_resolvent_controls=weights,
                common_heat_controls=heat,finite_covariance_partition=parts,
                source_sha256=sha(__file__),input_sha256=inputs,
                full_weighted_HS_squared=None,full_weighted_HS_frequency_derivative=None,
                physical_source=None,source_sign=None,atom_transit=None,
                scope=['The lapse commutator and insertion factor are bounded on the actual sewn geometry, for all angular states.',
                       'The completion and unsplit cyclic identities retain the full proper operator and common regulator.',
                       'Finite matrices validate identities and Schatten inequalities; they are not the atomic geometry source.',
                       'Weighted full-field Hilbert-Schmidt smoothing norms and the full high-frequency remainder remain unevaluated.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Source completion checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=result['passed'],tests=tests,global_weight=geometry,insertion_factor=insertion),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=result['finite_covariance_partition']
    fig,ax=plt.subplots(figsize=(8.5,4.6),layout='constrained')
    labels=['Selected','Transverse / outside tube','Window','Gaussian']
    bottom=np.zeros(len(rows));x=np.arange(len(rows))
    for j,label in enumerate(labels):
        values=[r['fractions'][j] for r in rows]
        ax.bar(x,values,bottom=bottom,label=label);bottom+=values
    ax.set(xticks=x,xticklabels=[str(r['Lambda']) for r in rows],
           xlabel='Gaussian cutoff in the finite algebra control',ylabel='Fraction of positive filtered covariance',
           title='Every omitted field contribution remains in the accounting',ylim=(0,1))
    ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1,1))
    fig.supxlabel('Six-dimensional algebra control. Covariance fractions are not stress fractions or atom-transmission probabilities.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
