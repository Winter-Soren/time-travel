"""Enclosed Dirichlet spatial-complement gap and numerical coupling diagnostics.

The gap is a continuum form bound. Ritz residual integrals are floating
quadrature pilots, not certified errors or pointwise heat bounds.
"""
import hashlib
import inspect
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
from scipy.linalg import eigh
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from dihole_continuum_heat_tail import target
from dihole_finite_heat import sine_table, form_matrices, physical_setup
from dihole_optical_heat import optical_fields
from dihole_quadrature_remainder import dyadic, outward
from null_saddle_leads import AU, AV, SY
from project_paths import project_path


def evaluators(iv):
    rho,z=sp.symbols('rho z',real=True)
    adapter=SimpleNamespace(mpf=sp.Rational,sqrt=sp.sqrt,log=sp.log)
    a,b,g=optical_fields(adapter)(rho,z)
    ar,az,br,bz=[sp.diff(f,x) for f,x in ((a,rho),(a,z),(b,rho),(b,z))]
    first=[sp.exp(a),ar,az,br,bz]
    enclosed=sp.lambdify((rho,z),first,
        modules=[dict(sqrt=iv.sqrt,log=iv.log,exp=iv.exp,mpf=iv.mpf),'mpmath'],cse=True)
    assert set(enclosed.__code__.co_names)<={'sqrt','mpf','log','exp'}
    # Floating differentiation is symbolic; only its evaluation/integration is numerical.
    full=[sp.exp(-a),ar,az,sp.diff(ar,rho)+sp.diff(az,z),sp.exp(-b),br,bz,
          g,sp.diff(g,rho),sp.diff(g,z)]
    floating=sp.lambdify((rho,z),full,'numpy',cse=True)
    return enclosed,floating,dict(interval=inspect.getsource(enclosed),floating=inspect.getsource(floating))


def cover(iv,function,constants,divisions):
    L,rs,vs=(constants[k] for k in ('L','rho_star','vstar'))
    assert divisions>0 and divisions&(divisions-1)==0
    wlo,whi,Bhi=np.inf,-np.inf,-np.inf
    width=8/divisions
    for i in range(divisions):
        for j in range(divisions):
            u=iv.mpf([-4+i*width,-4+(i+1)*width])
            v=iv.mpf([-4+j*width,-4+(j+1)*width])
            rho,z=rs+L*(u+v)/iv.sqrt(2),L*(u-v)/iv.sqrt(2)
            alpha,ar,az,br,bz=function(rho,z)
            assert alpha>0 and rho>0
            ww=outward(1/(vs*alpha)**2)
            BB=outward(L**2*((2*ar+br)**2+(2*az+bz)**2+ar**2+az**2)/4)
            wlo=min(wlo,ww['lower']);whi=max(whi,ww['upper']);Bhi=max(Bhi,BB['upper'])
    return dict(divisions=divisions,cells=divisions**2,radius=4,
        w_squared_lower=wlo,w_squared_upper=whi,connection_squared_upper=Bhi)


def complement_gap(iv,bounds,C,count):
    X=dyadic(iv,bounds['w_squared_lower'])*iv.pi**2*((count+1)**2+1)/64
    Y=dyadic(iv,bounds['w_squared_upper'])*dyadic(iv,bounds['connection_squared_upper'])
    assert X>Y and Y>0
    # This exact dyadic epsilon is admissible; optimality is not needed by the proof.
    epsilon=float(np.sqrt(outward(Y)['upper']/outward(X)['lower']))
    assert 0<epsilon<1
    e=dyadic(iv,epsilon)
    gamma=(1-e)*X-(1/e-1)*Y-dyadic(iv,C)
    return dict(count=count,epsilon=epsilon,X=outward(X),Y=outward(Y),gap=outward(gamma))


def coefficients(function,constants,u,v,offset):
    L,rs,vs=(constants[k] for k in ('L','rho_star','vstar'))
    E=vs/L
    rho,z=rs+L*(u+v)/np.sqrt(2),L*(u-v)/np.sqrt(2)
    speed,ar,az,lap_a,beta_inv,br,bz,g,gr,gz=function(rho,z)
    w=speed/vs
    wr,wz=-w*ar,-w*az
    wu,wv=L*(wr+wz)/np.sqrt(2),L*(wr-wz)/np.sqrt(2)
    lap_w=L*L*w*(ar*ar+az*az-lap_a)
    d=398855.5+offset-512000*g
    mass=beta_inv*d/E
    mr,mz=-beta_inv*(br*d+512000*gr)/E,-beta_inv*(bz*d+512000*gz)/E
    mu,mv=L*(mr+mz)/np.sqrt(2),L*(mr-mz)/np.sqrt(2)
    return w,wu,wv,lap_w,mass,mu,mv


def square_action(phi,du,dv,lap,coefficient):
    """A=-K^2 in the actual canonical frame; arrays have (point,spin,column)."""
    w,wu,wv,lapw,m,mu,mv=[np.asarray(x)[...,None,None] for x in coefficient]
    J=SY.imag
    potential=m*m-w*lapw/2-(wu*wu+wv*wv)/4
    return (-w*w*lap-2*w*(wu*du+wv*dv)+potential*phi
        +w*np.einsum('ab,pbk->pak',J,wu*dv-wv*du)
        +w*(mv*np.einsum('ab,pbk->pak',AU.real,phi)
            -mu*np.einsum('ab,pbk->pak',AV.real,phi)))


def controls(iv,enclosed,constants):
    mp=MPContext();mp.dps=85;fields=optical_fields(mp)
    checks=[]
    for rho,z in (('.6672354660415234','0'),('.61','.03'),('.73','-.04')):
        r,s=mp.mpf(rho),mp.mpf(z)
        a=lambda x,y:fields(x,y)[0]
        b=lambda x,y:fields(x,y)[1]
        truth=[mp.exp(a(r,s)),mp.diff(a,(r,s),(1,0)),mp.diff(a,(r,s),(0,1)),
               mp.diff(b,(r,s),(1,0)),mp.diff(b,(r,s),(0,1))]
        values=enclosed(iv.mpf(rho),iv.mpf(z))
        endpoints=[tuple(mp.make_mpf(x) for x in v._mpi_) for v in values]
        assert all(lo<=t<=hi for (lo,hi),t in zip(endpoints,truth))
        checks.append(float(max(abs((lo+hi)/2-t)/(1+abs(t)) for (lo,hi),t in zip(endpoints,truth))))
    # Independent polynomial K(K(psi)) tests all variable-coefficient terms/signs.
    u,v=sp.symbols('u v',real=True)
    U,V,J=sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2),sp.Matrix([[-1,1],[1,1]])/sp.sqrt(2),sp.Matrix([[0,-1],[1,0]])
    w=1+u/10+v*v/20;m=(u*v+u)/3
    psi=sp.Matrix([1+u*u*v,u*v*v+u])
    K=lambda f:-U*(w*f.diff(u)+sp.diff(w,u)*f/2)-V*(w*f.diff(v)+sp.diff(w,v)*f/2)+J*m*f
    exact=-K(K(psi))
    funcs=sp.lambdify((u,v),[psi,psi.diff(u),psi.diff(v),psi.diff(u,2)+psi.diff(v,2),
        w,sp.diff(w,u),sp.diff(w,v),sp.diff(w,u,2)+sp.diff(w,v,2),m,sp.diff(m,u),sp.diff(m,v),exact],'numpy')
    errors=[]
    for x,y in ((.3,-.2),(-.7,.8),(1.,2.)):
        got=funcs(x,y)
        action=square_action(*(np.asarray(f)[None,:,:] for f in got[:4]),got[4:11])[0]
        errors.append(float(np.max(abs(action-got[11]))))
    assert max(errors)<2e-14,errors
    return dict(passed=True,independent_derivative_errors=checks,nonconstant_K_squared_errors=errors)


def basis_grid(count,order):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    nodes,weights=4*nodes,4*weights
    basis,derivative=sine_table(nodes,count,4)
    product=lambda a,b:np.einsum('pi,qj->pqij',a,b).reshape(order*order,count*count)
    phi,du,dv=product(basis,basis),product(derivative,basis),product(basis,derivative)
    freq=np.arange(1,count+1)*np.pi/8
    lap=-phi*(freq[:,None]**2+freq[None,:]**2).ravel()
    spin=lambda b:np.einsum('pk,ab->pakb',b,np.eye(2)).reshape(order*order,2,2*count*count)
    u,v=np.meshgrid(nodes,nodes,indexing='ij')
    return (u.ravel(),v.ravel(),np.outer(weights,weights).ravel()),tuple(spin(b) for b in (phi,du,dv,lap))


def coupling_pilot(function,constants,gaps):
    old=project_path('dihole_heat_angular_state.npz')
    assert hashlib.sha256(old.read_bytes()).hexdigest()=='218545a22a2f1d9d757f7d7535c5534c82a2233e1cdac0b938dd885ceed01bdd'
    coefficient,_=physical_setup();rows=[];saved={}
    with np.load(old) as state:
        frozen=tuple(state[f'M{i}'].copy() for i in range(3))
    for count,order in ((14,48),(18,60)):
        matrices=frozen if count==18 else form_matrices(count,4,order,coefficient)[0]
        grid,basis=basis_grid(count,order)
        u,v,weight=grid;phi,du,dv,lap=basis
        gap=next(g['gap']['lower'] for g in gaps if g['count']==count)
        for offset in (-128,0,128):
            M=matrices[0]+offset*matrices[1]+offset**2*matrices[2]
            values,vectors=eigh(M,subset_by_index=(0,7),driver='evr',check_finite=False)
            action=square_action(phi,du,dv,lap,coefficients(function,constants,u,v,offset))
            # Green's identity checks the second-order action against the existing form.
            direct=phi.reshape(-1,2*count*count).T@(np.repeat(weight,2)[:,None]*action.reshape(-1,2*count*count))
            mismatch=float(np.max(abs(direct-M)))
            assert mismatch<2e-7,mismatch
            residual=action@vectors-(phi@vectors)*values[None,None,:]
            norms=np.sqrt(np.einsum('p,pak,pak->k',weight,residual,residual))
            grid2,basis2=basis_grid(count,order+12)
            u2,v2,weight2=grid2
            fields2=[b@vectors for b in basis2]
            action2=square_action(*fields2,coefficients(function,constants,u2,v2,offset))
            residual2=action2-fields2[0]*values[None,None,:]
            refined=np.sqrt(np.einsum('p,pak,pak->k',weight2,residual2,residual2))
            row=dict(count=count,offset=offset,order=order,refined_order=order+12,
                ritz_values=values.tolist(),residual_norms=norms.tolist(),refined_residual_norms=refined.tolist(),
                quadrature_change=float(np.max(abs(norms-refined))),form_action_max_difference=mismatch,
                complement_gap_lower=gap)
            if count==18 and offset==0:
                coupling=action-phi@M
                weighted=np.repeat(np.sqrt(weight),2)[:,None]*coupling.reshape(-1,2*count*count)
                gram=weighted.T@weighted
                norm2=float(eigh(gram,subset_by_index=(len(M)-1,len(M)-1),check_finite=False)[0][0])
                row['full_coupling_norm_estimate']=np.sqrt(norm2)
                row['zero_shift_schur_norm_upper_estimate']=norm2/gap
            if count==18:
                saved[f'j{offset}_vectors']=vectors;saved[f'j{offset}_values']=values
            rows.append(row)
            print(json.dumps(dict(stage='spatial_coupling_pilot',count=count,offset=offset,
                first_value=values[0],first_residual=refined[0],max_quadrature_change=row['quadrature_change'])),flush=True)
    state_path=project_path('dihole_spatial_complement_state.npz')
    np.savez_compressed(state_path,**saved)
    return rows,hashlib.sha256(state_path.read_bytes()).hexdigest()


def run():
    started=time.perf_counter();iv=MPIntervalContext();iv.dps=60
    constants=target(iv);enclosed,floating,generated=evaluators(iv)
    checks=controls(iv,enclosed,constants)
    old=project_path('dihole_continuum_heat_tail_results.json')
    assert hashlib.sha256(old.read_bytes()).hexdigest()=='4b3f5bb7dfcdebeabcaf8bd79336307432ebdcbfc27e9d3caed3340670bd108a'
    C=json.loads(old.read_text())['covers'][-1]['maxima']['C']
    covers=[cover(iv,enclosed,constants,n) for n in (16,32,64)]
    gaps=[complement_gap(iv,covers[-1],C,n) for n in (10,14,18,22)]
    floats={k:float(v.a) for k,v in constants.items()}
    pilots,state_hash=coupling_pilot(floating,floats,gaps)
    inputs=['dihole_continuum_heat_tail.py','dihole_continuum_heat_tail_results.json',
        'dihole_finite_heat.py','dihole_heat_angular_state.npz','dihole_optical_heat.py',
        'dihole_probe_certificate_results.json','dihole_quadrature_remainder.py','null_saddle_leads.py','project_paths.py']
    result=dict(passed=True,controls=checks,covers=covers,endomorphism_lower=-C,gaps=gaps,pilots=pilots,
        target='Exact A=4 Dirichlet form square H_min*H_min, q=512000, angular center398855.5.',
        certified='Whole-box interval geometry and the stated Young/Poincare form inequality certify the complementary-space gaps, uniformly in angular offset.',
        numerical_only='Ritz values, continuum differential residual integrals and coupling norms use floating quadrature/eigensolves. Refinement changes are not certified errors.',
        limitations=['The complementary subspace is not invariant; its gap alone cannot bound omitted heat.',
            'The Schur formula needs coupling bounds; individual residuals do not control the whole coupling operator.',
            'No pointwise heat certificate, natural-exterior boundary bound, renormalized stress, stable medium or backward passage.'],
        formulas=dict(gap='(1-epsilon)*wmin^2*pi^2*((N+1)^2+1)/64-C-(1/epsilon-1)*wmax^2*Bmax^2',
            connection='Bmax^2 >= L^2/4 sup(|grad log(alpha^2 beta)|^2+|grad log alpha|^2)',
            schur='For real z<gamma: 0 <= R*(A_Q-z)^-1 R <= R*R/(gamma-z).',
            residual='For normalized v in D(A): dist(lambda,spectrum(A)) <= ||(A-lambda)v||; this does not identify or count eigenvalues.'),
        generated_evaluators=generated,generated_sha256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in generated.items()},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),state_sha256=state_hash,
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in inputs},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_spatial_complement_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,seconds=result['runtime_seconds'],covers=covers,gaps=gaps),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,4.8))
    left.plot([r['count'] for r in result['gaps']],[r['gap']['lower'] for r in result['gaps']],'o-')
    left.axhline(0,color='#888888',lw=1)
    left.set(xlabel='Retained sine modes per coordinate',ylabel='Certified lower spectral bound',title='Energy bound for the omitted subspace')
    for offset in (-128,0,128):
        rows=[r for r in result['pilots'] if r['offset']==offset]
        right.semilogy([r['count'] for r in rows],[max(r['refined_residual_norms']) for r in rows],'o-',label=f'Angular offset {offset:+d}')
    right.set(xlabel='Retained sine modes per coordinate',ylabel='Largest residual among eight lowest Ritz states',title='Continuum residual pilot (not certified)')
    right.legend(fontsize=8)
    for ax in (left,right):ax.grid(alpha=.2)
    fig.suptitle('Spatial truncation: a proven gap and a measured coupling test',fontsize=14)
    fig.text(.5,.02,'Dirichlet box only. Residual integration, pointwise heat, exterior matching and physical support remain open.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.94))
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_spatial_complement_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
