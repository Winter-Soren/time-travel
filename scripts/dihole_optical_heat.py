"""Local ultraviolet subtraction coefficients for the full optical Dirac field.

Rank four, one massless charged Dirac species, D=1. The heat diagonal trace is
(4*pi*s)^(-3/2)*(a0+s*a2+s^2*a4+...). These local coefficients include all
spatial/angular ultraviolet degrees of freedom, unlike one packet resolvent.
They do not determine the finite heat kernel or a renormalized physical stress.
Numerical fourth derivatives are precision-checked, not interval-certified.
"""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from mpmath.ctx_mp import MPContext

from ori_disk_constraints import metric_geometry
from project_paths import project_path


def optical_fields(mp):
    """Exact balanced geometry, returning log(alpha), log(beta), unit gauge."""
    k,Q = mp.mpf(1)/2,mp.mpf(1)/5
    M = (mp.sqrt(65)-5)/20
    ag = mp.sqrt(k*k-M*M)
    magnetic = M/((k+ag)*(M+k))
    def fields(rho,z):
        plus,minus = mp.sqrt(rho*rho+(z+k)**2),mp.sqrt(rho*rho+(z-k)**2)
        r,c = M+(plus+minus)/2,2*z/(plus+minus)
        delta,sigma = r*r-2*M*r-ag*ag,r*r-ag*ag*c*c
        sin2 = 1-c*c
        W = (r*r-ag*ag)**2+delta*ag*ag*sin2
        F = delta+ag*ag*sin2+2*magnetic*M*r*ag*sin2+magnetic**2*W*sin2/4
        N = F/sigma
        alpha = sigma**2/(plus*minus)**2
        beta = rho/N**2
        unit_gauge = (2*M*r*ag+magnetic*W/2)*sin2/(2*Q*F)
        return mp.log(alpha),mp.log(beta),unit_gauge
    return fields


def local_invariants(mp,fields,rho,z,charge=512000):
    """Warped circle over a conformal two-dimensional base; sphere-positive R."""
    alpha_log = lambda x,y:fields(x,y)[0]
    beta_log = lambda x,y:fields(x,y)[1]
    def base(x,y):
        a,b,_ = fields(x,y)
        ar,az = (mp.diff(alpha_log,(x,y),order) for order in ((1,0),(0,1)))
        br,bz = (mp.diff(beta_log,(x,y),order) for order in ((1,0),(0,1)))
        al = mp.diff(alpha_log,(x,y),(2,0))+mp.diff(alpha_log,(x,y),(0,2))
        brr,brz,bzz = (mp.diff(beta_log,(x,y),order) for order in ((2,0),(1,1),(0,2)))
        scale = mp.exp(-2*a)
        K = -scale*al
        Qrr = scale*(brr+br*br-ar*br+az*bz)
        Qzz = scale*(bzz+bz*bz+ar*br-az*bz)
        Qrz = scale*(brz+br*bz-az*br-ar*bz)
        S = Qrr+Qzz
        R = 2*(K-S)
        ricci2 = (K-Qrr)**2+(K-Qzz)**2+2*Qrz**2+S*S
        return R,ricci2,scale,br,bz,a,b,(K-Qrr,-Qrz,K-Qzz,-S)
    R,ricci2,scale,br,bz,a,b,ricci = base(rho,z)
    scalar = lambda x,y:base(x,y)[0]
    Rr,Rz = (mp.diff(scalar,(rho,z),order) for order in ((1,0),(0,1)))
    lapR = scale*(mp.diff(scalar,(rho,z),(2,0))+mp.diff(scalar,(rho,z),(0,2))+br*Rr+bz*Rz)
    gauge = lambda x,y:fields(x,y)[2]
    gr,gz = (mp.diff(gauge,(rho,z),order) for order in ((1,0),(0,1)))
    F2 = 2*charge**2*mp.exp(-2*a-2*b)*(gr*gr+gz*gz)
    riemann2 = 4*ricci2-R*R
    a4gravity = (R*R-3*ricci2-lapR)/30
    a4gauge = 2*F2/3
    return dict(R=R,Ricci_squared=ricci2,Riemann_squared=riemann2,Laplacian_R=lapR,
                F_squared=F2,a0=mp.mpf(4),a2=-R/3,a4_gravity=a4gravity,
                a4_gauge=a4gauge,a4=a4gravity+a4gauge,alpha=mp.exp(a),beta=mp.exp(b),
                R_rho=Rr,R_z=Rz,Ricci_rr_frame=ricci[0],Ricci_rz_frame=ricci[1],
                Ricci_zz_frame=ricci[2],Ricci_phi_frame=ricci[3])


def lapse_counterterm_control(mp,fields,rho,z,center,length):
    """Direct metric perturbation versus the local conformal density identities.

The Gaussian equals120's central bump within its untapered patch. All points
used here are inside that patch; the global bump remains smooth and compact.
"""
    f=lambda x,y:mp.exp(-((x-center)**2+y*y)/(8*length*length))
    a=lambda x,y:fields(x,y)[0]
    b=lambda x,y:fields(x,y)[1]
    def laplacian(function,x,y):
        br,bz=mp.diff(b,(x,y),(1,0)),mp.diff(b,(x,y),(0,1))
        return mp.exp(-2*a(x,y))*(mp.diff(function,(x,y),(2,0))+mp.diff(function,(x,y),(0,2))
            +br*mp.diff(function,(x,y),(1,0))+bz*mp.diff(function,(x,y),(0,1)))
    baseline=local_invariants(mp,fields,rho,z)
    value=f(rho,z); scale=mp.exp(-2*a(rho,z))
    ar,az=mp.diff(a,(rho,z),(1,0)),mp.diff(a,(rho,z),(0,1))
    br,bz=mp.diff(b,(rho,z),(1,0)),mp.diff(b,(rho,z),(0,1))
    fr,fz=mp.diff(f,(rho,z),(1,0)),mp.diff(f,(rho,z),(0,1))
    frr,frz,fzz=(mp.diff(f,(rho,z),order) for order in ((2,0),(1,1),(0,2)))
    hrr=scale*(frr-ar*fr+az*fz); hzz=scale*(fzz+ar*fr-az*fz)
    hrz=scale*(frz-az*fr-ar*fz); hphi=scale*(br*fr+bz*fz)
    ric_hess=(baseline['Ricci_rr_frame']*hrr+2*baseline['Ricci_rz_frame']*hrz
              +baseline['Ricci_zz_frame']*hzz+baseline['Ricci_phi_frame']*hphi)
    grad=scale*(fr*baseline['R_rho']+fz*baseline['R_z'])
    lapf=lambda x,y:laplacian(f,x,y)
    expected=dict(a0=-3*value*baseline['a0'],a2=-value*baseline['a2']-4*lapf(rho,z)/3,
                  a4=value*baseline['a4']-(6*ric_hess+3*grad+4*laplacian(lapf,rho,z))/30)
    def perturbed(e):
        def newfields(x,y):
            aa,bb,A=fields(x,y); shift=mp.log(1+e*f(x,y))
            return aa-shift,bb-shift,A
        got=local_invariants(mp,newfields,rho,z)
        return {key:got[key]/(1+e*value)**3 for key in expected}
    errors=[]
    for step in (mp.mpf('1e-5'),mp.mpf('5e-6')):
        values=[perturbed(j*step) for j in (-2,-1,1,2)]
        derivatives={key:(values[0][key]-8*values[1][key]+8*values[2][key]-values[3][key])/(12*step) for key in expected}
        error={key:abs(derivatives[key]-expected[key])/(1+abs(expected[key])) for key in expected}
        assert max(error.values())<mp.mpf('1e-16'),error
        errors.append({key:mp.nstr(value,12) for key,value in error.items()})
    return dict(rho=mp.nstr(rho,30),z=mp.nstr(z,30),local_bump_value=mp.nstr(value,30),
        expected_density_derivatives={key:mp.nstr(value,40) for key,value in expected.items()},
        five_point_errors=errors,steps=['1e-5','5e-6'],
        convention='Derivative of dV_e*a_j(e) divided by baseline dV; source/gauge fixed. Local total derivatives retained.')


def tensor_control(mp,fields,rho,z,reference):
    """Independent existing general Christoffel implementation, no warp formulas."""
    metric = lambda x,y:np.array([float(mp.exp(2*fields(x,y)[i])) for i in (0,0,1)])
    g=np.diag(metric(rho,z)); dg=np.zeros((3,3,3)); ddg=np.zeros((3,3,3,3))
    for axis,log_index in enumerate((0,0,1)):
        component=lambda x,y:mp.exp(2*fields(x,y)[log_index])
        for i,order in enumerate(((1,0),(0,1))):
            dg[i,axis,axis]=float(mp.diff(component,(rho,z),order))
        for i,j,order in ((0,0,(2,0)),(0,1,(1,1)),(1,0,(1,1)),(1,1,(0,2))):
            ddg[i,j,axis,axis]=float(mp.diff(component,(rho,z),order))
    inverse,_,connection,derivative,scalar=metric_geometry(g,dg,ddg)
    ricci=np.array([[sum(derivative[k,k,i,j]-derivative[j,k,i,k]+sum(
        connection[k,k,l]*connection[l,i,j]-connection[k,j,l]*connection[l,i,k]
        for l in range(3)) for k in range(3)) for j in range(3)] for i in range(3)])
    norm=float(np.einsum('ij,kl,ik,jl',ricci,ricci,inverse,inverse))
    riemann=np.empty((3,3,3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    riemann[i,j,k,l]=g[i,i]*(derivative[k,i,l,j]-derivative[l,i,k,j]+sum(
                        connection[i,k,m]*connection[m,l,j]-connection[i,l,m]*connection[m,k,j] for m in range(3)))
    iv=np.diag(inverse)
    rnorm=float(np.einsum('abcd,a,b,c,d->',riemann**2,iv,iv,iv,iv))
    errors={key:abs(value-float(reference[key]))/(1+abs(float(reference[key])))
            for key,value in dict(R=float(scalar),Ricci_squared=norm,Riemann_squared=rnorm).items()}
    assert max(errors.values())<2e-12,errors
    return errors


def checks(mp):
    cases=[('flat_uniform_B',lambda r,z:(mp.mpf(0),mp.log(r),r*r/2),1,0,
            dict(R=0,Ricci_squared=0,Laplacian_R=0,F_squared=2,a2=0,a4=mp.mpf(4)/3)),
           ('unit_sphere_times_line',lambda r,z:(mp.mpf(0),mp.log(mp.sin(r)),mp.mpf(0)),mp.mpf('.7'),0,
            dict(R=2,Ricci_squared=2,Laplacian_R=0,F_squared=0,a2=-mp.mpf(2)/3,a4=-mp.mpf(1)/15)),
           ('nonconstant_warp',lambda r,z:(mp.mpf(0),r*r/2,mp.mpf(0)),mp.mpf('.5'),0,
            dict(R=-mp.mpf(5)/2,Ricci_squared=mp.mpf(25)/8,Laplacian_R=-5,F_squared=0,a2=mp.mpf(5)/6,a4=mp.mpf(1)/16))]
    reports=[]
    for name,fields,r,z,expected in cases:
        got=local_invariants(mp,fields,mp.mpf(r),mp.mpf(z),charge=1)
        err=max(abs(got[key]-value) for key,value in expected.items())
        assert err<mp.mpf('1e-45'),(name,err)
        reports.append(dict(name=name,maximum_analytic_error=mp.nstr(err,10),
                            tensor_errors=tensor_control(mp,fields,mp.mpf(r),mp.mpf(z),got)))
    # Independent exact Landau heat diagonal: rank4 * sB*coth(sB).
    magnetic=[]
    for s in map(mp.mpf,('.1','.03','.01','.003')):
        exact=4*s/mp.tanh(s)
        remainder=exact-4-4*s*s/3
        coefficient=remainder/s**4
        assert abs(coefficient+mp.mpf(4)/45)<s*s/100
        magnetic.append(dict(s=float(s),flat_subtracted_integrand=float((exact-4)/s**3),
                             fully_subtracted_integrand=float(remainder/s**3),
                             next_coefficient=float(coefficient)))
    # Full-rank S2 x R heat diagonal from the exact Dirac spectrum, radius1.
    sphere=[]
    for s in map(mp.mpf,('.02','.01','.005')):
        stop=int(mp.ceil(mp.sqrt(160/s)))
        normalized=8*s*mp.fsum(n*mp.exp(-s*n*n) for n in range(1,stop+1))
        a4=(normalized-4+2*s/3)/s**2
        assert abs(a4+mp.mpf(1)/15)<s/60
        sphere.append(dict(s=float(s),spectrum_terms=stop,inferred_a4=float(a4)))
    return dict(passed=True,geometric_cases=reports,exact_uniform_magnetic_heat=magnetic,
                exact_sphere_product_heat=sphere)


def run():
    started=time.perf_counter(); mp=MPContext();mp.dps=60
    controls=checks(mp)
    certificate=json.loads(project_path('dihole_probe_certificate_results.json').read_text())
    assert hashlib.sha256(project_path('dihole_probe_certificate_results.json').read_bytes()).hexdigest()=='d3c3e55a125c539baa8e9b5c41f97429e2759850a419c23280b70522e26f2843'
    def exact_float(value,context):
        n,d=float.fromhex(value).as_integer_ratio();return context.mpf(n)/d
    constants=certificate['target']['coordinate_constants']
    rs,L=exact_float(constants['rho_star'],mp),exact_float(constants['L'],mp)
    field=optical_fields(mp)
    rows=[]
    for c in range(-6,7):
        cc=mp.mpf(c)/2
        rho,z=rs+mp.sqrt(2)*L*cc,mp.mpf(0)
        invariants=local_invariants(mp,field,rho,z)
        rows.append(dict(c=float(cc),rho=mp.nstr(rho,42),z='0',
                         values={key:mp.nstr(value,62) for key,value in invariants.items()}))
    for rho,z in ((mp.mpf('.64'),mp.mpf('.03')),(mp.mpf('.7'),mp.mpf('.08'))):
        invariants=local_invariants(mp,field,rho,z)
        rows.append(dict(c=None,rho=mp.nstr(rho,42),z=mp.nstr(z,42),
                         values={key:mp.nstr(value,62) for key,value in invariants.items()}))
    print(json.dumps(dict(stage='optical_invariants',points=len(rows),seconds=time.perf_counter()-started)),flush=True)
    high=MPContext();high.dps=90;highfield=optical_fields(high)
    precision=[];tensor=[]
    for index in (0,6,12,14):
        row=rows[index]
        if row['c'] is not None:
            rho=exact_float(constants['rho_star'],high)+high.sqrt(2)*exact_float(constants['L'],high)*high.mpf(row['c'])
            z=high.mpf(0)
        else:
            rho,z=high.mpf('.7'),high.mpf('.08')
        got=local_invariants(high,highfield,rho,z)
        errors={key:abs(high.mpf(row['values'][key])-value)/(1+abs(value)) for key,value in got.items()}
        assert max(errors.values())<high.mpf('1e-50'),errors
        precision.append(dict(index=index,maximum_normalized_change=high.nstr(max(errors.values()),12)))
        tensor.append(dict(index=index,errors=tensor_control(high,highfield,rho,z,got)))
    controls['lapse_counterterm_variations']=[lapse_counterterm_control(mp,field,rho,z,rs,L)
        for rho,z in ((rs,mp.mpf(0)),(mp.mpf('.64'),mp.mpf('.03')))]
    names=['dihole_probe_certificate_results.json','ori_disk_constraints.py','project_paths.py']
    result=dict(passed=True,controls=controls,points=rows,precision_comparisons=precision,tensor_comparisons=tensor,
        conventions=dict(rank=4,charge=512000,geometry='Exact balanced Q=1/5,D=1, classical uncut exterior.',
            optical_metric='alpha^2*(d_rho^2+d_z^2)+beta^2*d_phi^2; alpha=P/N,beta=rho/N^2.',
            gauge='A_phi=q*aunit, F_squared=2*q^2*|grad aunit|^2/(alpha^2*beta^2).',
            curvature='R sphere-positive; Delta=div grad, nonpositive on oscillatory flat modes.',
            heat_diagonal='tr K(s;x,x)~(4*pi*s)^(-3/2)*(a0+s*a2+s^2*a4+...), measured per optical volume.',
            a0='4',a2='-R/3',a4='(R^2-3*Ricci_squared-Delta R)/30+(2/3)*F_squared',
            compact_lapse_integrated_variations='delta integral a_(2j) dV=(2j-3)*integral f*a_(2j) dV; weights(-3,-1,+1). This is a local-coefficient variation, not an unregulated infinite-volume trace.',
            dimensions='D=1: R and a2 have dimension D^-2; Ricci_squared, Delta R, F_squared and a4 have dimension D^-4.'),
        scope=['Local ultraviolet subtraction coefficients of the full rank-four operator, not a finite heat-kernel approximation with a remainder bound.',
            'mpmath differentiation and60/90-digit comparisons are numerical checks, not interval certificates.',
            'Exact geometry and charge normalization match119; no single angular-sector coefficient is substituted for the full diagonal.',
            'Optical subtraction is not a full locally covariant physical stress prescription: finite counterterms, conformal anomaly, nonlocal kernel and metric variation remain.',
            'No negative supporting stress, stable throat, new physical law, material recipe or backward passage.'],
        sources=['https://arxiv.org/html/hep-th/0306138v3','https://arxiv.org/pdf/hep-th/9906160',
                 'https://arxiv.org/pdf/1210.4031'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_optical_heat_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,central=rows[6],seconds=result['runtime_seconds'],precision=precision),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[row for row in result['points'] if row['c'] is not None]
    x=np.array([row['c'] for row in rows])
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,5.2))
    for key,label,color in (('a4_gauge','Magnetic contribution','#235b9f'),('a4_gravity','Curvature contribution','#b06a16')):
        left.plot(x,[float(row['values'][key]) for row in rows],'o-',label=label,color=color)
    left.set_yscale('symlog',linthresh=.01)
    left.set(xlabel='Radial probe position c (u = v = c)',ylabel=r'Local coefficient $a_4$ ($D^{-4}$)',
             title='The ultraviolet subtraction varies in space')
    left.legend(fontsize=9);left.grid(alpha=.15)
    check=result['controls']['exact_uniform_magnetic_heat']
    s=np.array([row['s'] for row in check])
    right.loglog(s,[row['flat_subtracted_integrand'] for row in check],'o-',label='Only flat subtraction',color='#a34537')
    right.loglog(s,[-row['fully_subtracted_integrand'] for row in check],'o-',label='Also subtract magnetic term',color='#247b62')
    right.set(xlabel='Proper time s (toward zero = shorter distance)',ylabel='Absolute regulated integrand (common factor omitted)',
              title='Exact uniform-field check of cancellation')
    right.legend(fontsize=9);right.grid(alpha=.15)
    fig.suptitle('Preparing the short-distance subtraction for quantum stress',fontsize=14)
    fig.text(.5,.018,'These are subtraction coefficients and an exact benchmark, not a computed supporting energy or a time-travel medium.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.93))
    for extension in ('png','pdf'):
        fig.savefig(project_path('dihole_optical_heat_progress.'+extension),dpi=170)
    plt.close(fig)


if __name__=='__main__':
    run()
