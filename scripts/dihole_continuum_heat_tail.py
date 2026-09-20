"""Whole-box continuum angular heat tail, with interval geometry coefficients.

The target is the exact balanced Dirichlet-box Friedrichs square, not a
Galerkin matrix. Retained-sector spatial and natural-exterior errors remain.
"""
import hashlib
import inspect
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from dihole_optical_heat import optical_fields, local_invariants
from dihole_quadrature_remainder import dyadic, outward
from project_paths import project_path


def symbolic_evaluator(context):
    """Reuse exact geometry, differentiate symbolically, preserve all rationals."""
    rho,z=sp.symbols('rho z',real=True)
    adapter=SimpleNamespace(mpf=sp.Rational,sqrt=sp.sqrt,log=sp.log)
    alpha_log,beta_log,gauge=optical_fields(adapter)(rho,z)
    lapse_log=(sp.log(rho)-beta_log)/2
    expressions=[sp.exp(alpha_log),sp.exp(beta_log),gauge,
                 sp.diff(gauge,rho),sp.diff(gauge,z),sp.diff(lapse_log,rho),sp.diff(lapse_log,z)]
    function=sp.lambdify((rho,z),expressions,
        modules=[dict(sqrt=context.sqrt,log=context.log,exp=context.exp,mpf=context.mpf),'mpmath'],cse=True)
    # The mpmath printer emits mpf(n)/mpf(d), never pre-rounded rational constants.
    assert set(function.__code__.co_names)<={'sqrt','mpf','log','exp'}
    return function,inspect.getsource(function)


def target(iv):
    path=project_path('dihole_probe_certificate_results.json')
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='d3c3e55a125c539baa8e9b5c41f97429e2759850a419c23280b70522e26f2843'
    saved=json.loads(path.read_text())['target']['coordinate_constants']
    constants={key:dyadic(iv,float.fromhex(value)) for key,value in saved.items()}
    constants['E']=constants['vstar']/constants['L']
    return constants


def intersect(iv,left,right):
    lo=max(left.a,right.a);hi=min(left.b,right.b)
    if not lo<=hi:raise ArithmeticError('Disjoint purported enclosures.')
    return iv.mpf([lo,hi])


def cover_box(iv,function,constants,divisions):
    """Tiles cover the entire rotated A=4 box, including every edge."""
    started=time.perf_counter();A=4.;q=512000;mu0=iv.mpf(797711)/2
    assert divisions>0 and divisions&(divisions-1)==0
    minima=dict(m=np.inf,beta=np.inf,rho=np.inf)
    maxima=dict(m=-np.inf,beta=-np.inf,D=-np.inf,C=-np.inf,direct_D=-np.inf)
    limiting={}
    E,L,rs=(constants[k] for k in ('E','L','rho_star'))
    width=2*A/divisions
    for i in range(divisions):
        for j in range(divisions):
            ulo,uhi=-A+i*width,-A+(i+1)*width
            vlo,vhi=-A+j*width,-A+(j+1)*width
            u=iv.mpf([dyadic(iv,ulo).a,dyadic(iv,uhi).b])
            v=iv.mpf([dyadic(iv,vlo).a,dyadic(iv,vhi).b])
            rho=rs+L*(u+v)/iv.sqrt(2);z=L*(u-v)/iv.sqrt(2)
            assert rho>0
            alpha,beta,gauge,gr,gz,nr,nz=function(rho,z)
            assert alpha>0 and beta>0
            uc,vc=dyadic(iv,(ulo+uhi)/2),dyadic(iv,(vlo+vhi)/2)
            rc,zc=rs+L*(uc+vc)/iv.sqrt(2),L*(uc-vc)/iv.sqrt(2)
            center_gauge=function(rc,zc)[2]
            # Mean-value enclosure on the containing physical rectangle.
            centered=center_gauge+gr*(rho-rc)+gz*(z-zc)
            gauge_tight=intersect(iv,gauge,centered)
            grad_gauge=iv.sqrt(gr**2+gz**2)
            grad_lapse_squared=nr**2+nz**2
            magnetic=q*grad_gauge/(alpha*beta)
            # Exact static electrovac identity implies R >= -6 alpha^-2 |grad log N|^2.
            potential=(magnetic+iv.mpf(3)*grad_lapse_squared/(2*alpha**2))/E**2
            quantities=dict(m=alpha**2*beta,beta=beta,rho=rho,D=abs(mu0-q*gauge_tight),
                            C=potential,direct_D=abs(mu0-q*gauge))
            for key,value in quantities.items():
                enclosure=outward(value)
                if key in minima:minima[key]=min(minima[key],enclosure['lower'])
                if key in maxima and enclosure['upper']>maxima[key]:
                    maxima[key]=enclosure['upper'];limiting[key]=[ulo,uhi,vlo,vhi]
        if i%16==0:
            print(json.dumps(dict(stage='continuum_geometry_cover',divisions=divisions,
                completed_rows=i+1,seconds=time.perf_counter()-started)),flush=True)
    mlo,mhi,blo,bhi=(dyadic(iv,x) for x in (minima['m'],maxima['m'],minima['beta'],maxima['beta']))
    nash=constants['vstar']**2*(mhi/mlo)**2/(iv.pi*blo)
    cmin=1/(E*bhi)
    return dict(divisions=divisions,closed_cells=divisions**2,radius=A,minima=minima,maxima=maxima,
        nash_constant=outward(nash),angular_slope_lower=outward(cmin)['lower'],limiting_cells=limiting,
        seconds=time.perf_counter()-started)


def continuum_tail(iv,constants,cover,cutoff,tau,method='half_time'):
    D=dyadic(iv,cover['maxima']['D']);n0=iv.mpf(cutoff+1)
    if not n0>D:return dict(cutoff=cutoff,tau=tau,method=method,applicable=False,reason='The whole-box mass-sign threshold is not passed.')
    t=dyadic(iv,tau);c=dyadic(iv,cover['angular_slope_lower'])
    C=dyadic(iv,cover['maxima']['C']);CN=dyadic(iv,cover['nash_constant']['upper'])
    distance=n0-D
    smoothing=tau/2
    if method=='chosen_time':
        approximate_gap=cover['angular_slope_lower']**2*(cutoff+1-cover['maxima']['D'])**2
        smoothing=min(smoothing,1/(cover['maxima']['C']+approximate_gap))
    elif method!='half_time':raise ValueError('Unknown time split.')
    assert 0<smoothing<tau
    # Any common split is valid: the floating choice is then an exact dyadic input.
    smooth=dyadic(iv,smoothing);a=(t-smooth)*c*c
    log_raw=iv.log(4*CN/(iv.pi*constants['L']**2*smooth))+C*smooth-a*distance**2-iv.log(1-iv.exp(-a*(2*distance+1)))
    log_normalized=log_raw+iv.mpf(3)*iv.log(4*iv.pi*t/constants['E']**2)/2
    upper=outward(iv.exp(log_normalized))['upper']
    return dict(cutoff=cutoff,tau=tau,method=method,smoothing_time=smoothing,applicable=True,normalized_tail_upper=upper,
                log10_normalized_tail_upper=outward(log_normalized/iv.log(10))['upper'])


def controls(iv,function,constants):
    mp=MPContext();mp.dps=85;fields=optical_fields(mp)
    reports=[]
    for rho,z in (('.6672354660415234','0'),('.61','.03'),('.73','-.04'),('.64','-.06')):
        r,y=mp.mpf(rho),mp.mpf(z)
        a,b,gauge=fields(r,y)
        gf=lambda x,z:fields(x,z)[2]
        nf=lambda x,z:(mp.log(x)-fields(x,z)[1])/2
        expected=[mp.exp(a),mp.exp(b),gauge,mp.diff(gf,(r,y),(1,0)),mp.diff(gf,(r,y),(0,1)),
                  mp.diff(nf,(r,y),(1,0)),mp.diff(nf,(r,y),(0,1))]
        enclosed=function(iv.mpf(rho),iv.mpf(z))
        # Independent high-precision differentiation checks the symbolic formulas.
        endpoints=[tuple(mp.make_mpf(x) for x in got._mpi_) for got in enclosed]
        assert all(lo<=want<=hi for (lo,hi),want in zip(endpoints,expected))
        errors=[abs((lo+hi)/2-want)/(1+abs(want)) for (lo,hi),want in zip(endpoints,expected)]
        assert max(errors)<mp.mpf('1e-14')
        inv=local_invariants(mp,fields,r,y)
        N=mp.sqrt(r/mp.exp(b));Q=mp.mpf(1)/5
        R_static=6*mp.exp(-2*a)*((2*Q*N/r)**2*(expected[3]**2+expected[4]**2)-expected[5]**2-expected[6]**2)
        error=abs(R_static-inv['R'])/(1+abs(inv['R']))
        assert error<mp.mpf('1e-70')
        reports.append(dict(rho=rho,z=z,symbolic_derivative_max_error=float(max(errors)),
                            independent_curvature_identity_error=mp.nstr(error,12)))
    # Every interval cell contains off-center points, not just its midpoint.
    rbox,zbox=iv.mpf(['.66','.68']),iv.mpf(['-.01','.02'])
    enclosure=function(rbox,zbox)
    rc,zc=iv.mpf('.67'),iv.mpf('.005')
    tight_gauge=intersect(iv,enclosure[2],function(rc,zc)[2]+enclosure[3]*(rbox-rc)+enclosure[4]*(zbox-zc))
    for r,y in (('.661','-.009'),('.679','.019'),('.670','.005')):
        point_values=function(iv.mpf(r),iv.mpf(y))
        assert all(a.a<=b.a and b.b<=a.b for a,b in zip(enclosure,point_values))
        assert tight_gauge.a<=point_values[2].a and point_values[2].b<=tight_gauge.b
    # Geometric Gaussian tail is valid in a constant-coefficient continuum control.
    a=iv.mpf(1)/4;d=iv.mpf(3)/2
    direct=sum((iv.exp(-a*(d+n)**2) for n in range(100)),iv.mpf(0))
    bound=iv.exp(-a*d*d)/(1-iv.exp(-a*(2*d+1)))
    assert direct<bound
    return dict(passed=True,point_checks=reports,whole_interval_containment=True,
        gaussian_tail_control=dict(direct=outward(direct),bound=outward(bound)))


def run():
    started=time.perf_counter();iv=MPIntervalContext();iv.dps=60
    constants=target(iv);function,generated_source=symbolic_evaluator(iv)
    checks=controls(iv,function,constants)
    covers=[cover_box(iv,function,constants,n) for n in (16,32,64)]
    chosen=covers[-1]
    tails=[continuum_tail(iv,constants,chosen,J,tau,method) for method in ('half_time','chosen_time')
           for J in (2047,2559,3071,4095) for tau in (.5,1.,2.)]
    names=['dihole_optical_heat.py','dihole_probe_certificate_results.json','dihole_quadrature_remainder.py',
           'dihole_heat_angular_results.json','project_paths.py']
    result=dict(passed=True,controls=checks,exact_target={k:outward(v) for k,v in constants.items()},
        angular_center=398855.5,charge=512000,covers=covers,tail_bounds=tails,
        geometric_derivation='Ropt=6 alpha^-2[(2Q N/rho)^2 |grad a|^2-|grad logN|^2]; drop the nonnegative electromagnetic term for a curvature lower bound.',
        scope=['Each closed-cell interval cover bounds every point of the exact rotated A4 box, with exact balanced radicals and dyadic coordinates.',
               'The omitted-angular heat bound includes every spatial mode of the Dirichlet-box Friedrichs square.',
               'It does not enclose the retained-sector spatial error, full heat sum, or difference from the natural exterior.',
               'No proper-time energy integral, renormalized stress, supporting material, stability or backward passage is established.'],
        assumptions='Smooth exact off-axis electrovac geometry, the stated Dirichlet form theorem, and correct mpmath interval arithmetic; no formal verification of dependencies.',
        sources=['https://arxiv.org/abs/hep-th/9906160','https://arxiv.org/html/1711.07225'],
        generated_evaluator_source=generated_source,
        generated_evaluator_sha256=hashlib.sha256(generated_source.encode()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names},
        versions=dict(sympy=sp.__version__,interval_digits=iv.dps),runtime_seconds=time.perf_counter()-started)
    project_path('dihole_continuum_heat_tail_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,seconds=result['runtime_seconds'],cover=chosen,tail_bounds=tails),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,5))
    resolutions=[c['divisions'] for c in result['covers']]
    left.plot(resolutions,[c['maxima']['direct_D'] for c in result['covers']],'o--',label='Direct gauge enclosure')
    left.plot(resolutions,[c['maxima']['D'] for c in result['covers']],'o-',label='Mean-value gauge enclosure')
    left.set(xlabel='Cover divisions per coordinate',ylabel='Whole-box bound on |mu0 - q a|',title='Certified geometric bounds')
    for method,label in (('half_time','Equal time split'),('chosen_time','Chosen time split')):
        rows=[r for r in result['tail_bounds'] if r['tau']==.5 and r['method']==method and r['applicable']]
        right.plot([r['cutoff'] for r in rows],[r['log10_normalized_tail_upper'] for r in rows],'o-',label=label)
    right.axhline(-10,color='#777777',ls='--',label='1e-10 reference')
    right.set(xlabel='Angular offset cutoff',ylabel='log10 of certified normalized upper bound',title='Continuum-box remainder at heat parameter 0.5')
    for ax in (left,right):ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('All spatial modes included in the omitted-angular bound',fontsize=14)
    fig.text(.5,.02,'Dirichlet box only. Retained-sector spatial error, exterior boundary effects and physical stress remain unresolved.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.94))
    for ext in ('png','pdf'):fig.savefig(project_path('dihole_continuum_heat_tail_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
