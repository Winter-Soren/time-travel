"""Radial Galerkin evolution of the nonlinear action144 candidate.

Analytic Bessel derivatives and oversampled products avoid the mismatched
grid operators rejected in the preserved finite-difference/spline pilots.
Finite basis/domain/time results do not prove nonlinear continuum stability.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn

from nonlinear_packet_constraints import background, evaluate, solve
from project_paths import project_path


class RadialModel:
    def __init__(self, modes=32, xmax=12., quadrature_factor=4):
        if not (isinstance(modes,int) and modes>=12 and xmax>=10 and quadrature_factor>=4):
            raise ValueError('Require modes>=12, xmax>=10 and quadrature_factor>=4.')
        self.base = background()
        self.modes, self.xmax = modes, xmax
        self.k = np.arange(1,modes+1)*np.pi/xmax
        nodes,weights = np.polynomial.legendre.leggauss(quadrature_factor*modes)
        self.x,self.weights = np.r_[0.,xmax*(nodes+1)/2],np.r_[0.,xmax*weights/2]
        self.B = spherical_jn(0,self.x[:,None]*self.k)
        self.Bx = -self.k*spherical_jn(1,self.x[:,None]*self.k)
        self.project = (2*self.k**2/xmax)[:,None]*(self.weights*self.x**2)[None,:]*self.B.T
        order=len(nodes)
        legendre=np.polynomial.legendre.legvander(nodes,order)
        anti=np.empty((order,order));anti[:,0]=nodes+1
        for degree in range(1,order):
            anti[:,degree]=(legendre[:,degree+1]-legendre[:,degree-1])/(2*degree+1)
        self.integrate = np.zeros((order+1,order+1))
        self.integrate[1:,1:] = xmax/2*(anti@(((2*np.arange(order)+1)/2)[:,None]
                                           *weights[None,:]*legendre[:,:order].T))
        local_nodes,local_weights=np.polynomial.legendre.leggauss(2*modes)
        self.v=(local_nodes+1)/2
        self.vweights=local_weights/2*self.v**3
        self.shear_x=np.r_[self.x,xmax]
        q=(self.shear_x[:,None]*self.v).ravel()
        self.local_B=spherical_jn(0,q[:,None]*self.k)
        self.local_Bx=-self.k*spherical_jn(1,q[:,None]*self.k)
        self.aux_k=np.arange(1,2*modes+1)*np.pi/xmax
        self.aux_B=spherical_jn(0,self.x[:,None]*self.aux_k)
        self.aux_Bx=-self.aux_k*spherical_jn(1,self.x[:,None]*self.aux_k)
        self.aux_project=(2*self.aux_k**2/xmax)[:,None]*(self.weights*self.x**2)[None,:]*self.aux_B.T

    def fields(self,tau,values):
        zc,mc=values.reshape(2,-1)
        P,c,X,H,R=[self.base[k] for k in ('P','c','X','H','radius')]
        z,mu=self.B@zc,self.B@mc
        zx,mux=self.Bx@zc,self.Bx@mc
        lap_z=self.B@(-self.k**2*zc)
        R3=np.exp(-2*(tau+z))/R**2*(-4*lap_z-2*zx*zx)
        local_z=(self.local_B@zc).reshape(len(self.shear_x),-1)
        local_mux=(self.local_Bx@mc).reshape(len(self.shear_x),-1)
        shear_z=np.r_[z,0.]
        shear=-self.shear_x*np.exp(-3*shear_z)*((np.exp(3*local_z)*local_mux)@self.vweights)
        all_mu=np.r_[mu,0.];all_R3=np.r_[R3,0.]
        dy=np.zeros_like(all_mu)
        m0=3*H*(P+3*c);dm=P*H*all_mu;m=m0+dm
        for _ in range(7):
            y=1+dy
            dk=(dm-9*H*c*dy)/(P+3*c*y)
            grav=(P+9*c)*(6*H*dk+dk*dk)/3+3*c*dy*(3*H+dk)**2
            residual=(.5*P*all_R3+grav-P*H*H*shear*shear/3
                      +dy*(.5*X-1.5*X*X-.75*X*X*dy))
            slope=c*m*m*(P-9*c*y)/(P+3*c*y)**3-1.5*X*X*y+.5*X
            if np.any(slope>=-.5) or np.any(y<=0):
                raise ArithmeticError(f'Lapse branch failed at Ht={tau:.9g}.')
            dy-=residual/slope
        if np.max(abs(residual))>2e-12:
            raise ArithmeticError(f'Lapse residual failed at Ht={tau:.9g}.')
        sqrt_y=np.sqrt(1+dy)
        n=-dy/(sqrt_y*(1+sqrt_y));N=1+n
        dk=(dm-9*H*c*dy)/(P+3*c*(1+dy))
        integrand=np.divide(N[:-1]*shear[:-1],self.x,out=np.zeros_like(self.x),where=self.x!=0)
        w=self.integrate@integrand-self.weights@integrand-shear[-1]/3
        nc=self.aux_project@(n[:-1]-n[-1])
        nx=self.aux_Bx@nc
        lap_n_flat=self.aux_B@(-self.aux_k**2*nc)
        nx_over_x=np.divide(nx,self.x,out=lap_n_flat/3,where=self.x!=0)
        lap_n=np.exp(-2*(tau+z))/R**2*(lap_n_flat+zx*nx)
        return dict(zeta=z,mu=mu,zx=zx,mux=mux,curvature=R3,
            shear=shear[:-1],lapse=N[:-1],n=n[:-1],dy=dy[:-1],dk=dk[:-1],
            lapse_x=nx,lapse_x_over_x=nx_over_x,lap_n=lap_n,w=w,maximum_lapse_residual=float(np.max(abs(residual))),
            minimum_lapse_slope=float(np.min(-slope)),endpoint_lapse_deviation=float(n[-1]))

    def rhs(self,tau,values):
        f=self.fields(tau,values)
        P,c,X,H=[self.base[k] for k in ('P','c','X','H')]
        N,n,dy,dk,s,w,mu=[f[k] for k in ('lapse','n','dy','dk','shear','w','mu')]
        zt=n+N*dk/(3*H)+self.x*w*f['zx']+w+N*s/3
        dm=P*H*mu;m0=3*H*(P+3*c)
        bracket=(3*H*dm+m0*dk+dm*dk+P*H*H*s*s
                 +.75*(2*(X-1)*X*dy+X*X*dy*dy)+.5*P*f['curvature'])
        mt=self.x*w*f['mux']+f['lap_n']/H**2-N*bracket/(2*P*H*H)
        return np.r_[self.project@zt,self.project@mt]

    def initial(self,epsilon=.001):
        solution,_=solve(epsilon,xmax=max(12.,self.xmax+1.),tolerance=1e-11)
        old=evaluate(solution,self.base,epsilon,np.maximum(self.x,solution.x[0]))
        N,e=old['lapse'],old['e'];P,c=self.base['P'],self.base['c']
        mu=((P+3*c/N**2)*e/N-3*(P+3*c))/P
        z=-epsilon*np.exp(-self.x**2/2)
        return np.r_[self.project@z,self.project@mu],old

    def evolve(self,epsilon=.001,end=.5,rtol=2e-9,atol=2e-12,max_step=.005):
        if not (0<=epsilon<=.001 and 0<end<=1 and rtol>0 and atol>0 and max_step>0):
            raise ValueError('Require tested amplitude/time range and positive integration tolerances.')
        values,_=self.initial(epsilon)
        solution=solve_ivp(self.rhs,(0.,end),values,method='DOP853',
                           t_eval=np.linspace(0.,end,11),rtol=rtol,atol=atol,max_step=max_step)
        if not solution.success:
            raise ArithmeticError(solution.message)
        return dict(tau=solution.t,coefficients=solution.y.T,nfev=solution.nfev)

    def stress(self,tau,values,difference_step=1e-4):
        f=self.fields(tau,values)
        rate=self.rhs(tau,values)
        side=[self.fields(tau+sign*difference_step,values+sign*difference_step*rate) for sign in (1,-1)]
        P,c,X,H,R=[self.base[k] for k in ('P','c','X','H','radius')]
        N,dy,dk,s,z,zx,nx,mu=[f[k] for k in ('lapse','dy','dk','shear','zeta','zx','lapse_x','mu')]
        chi=lambda v:c*(v['dk']+v['dy']*(3*H+v['dk']))
        hd=lambda v:(v['dk']+H*v['shear'])/3
        chi_tau=(chi(side[0])-chi(side[1]))/(2*difference_step)
        h_tau=(hd(side[0])-hd(side[1]))/(2*difference_step)
        yx=-2*nx/N**3
        dkx=(P*H*f['mux']-3*c*(3*H+dk)*yx)/(P+3*c*(1+dy))
        s_over_x=np.divide(s,self.x,out=np.zeros_like(self.x),where=self.x!=0)
        sx=-f['mux']-3*zx*s-3*s_over_x
        hx=(dkx+H*sx)/3
        chi_x=c*((1+dy)*dkx+(3*H+dk)*yx)
        inverse_A=np.exp(-tau-z)/R
        q=6*H*dk+dk*dk+dy*(3*H+dk)**2
        rho=X*dy*(1.5*X-.5)+.75*X*X*dy*dy-3*c*q
        mean=X*dy*(2*X-1)+X*X*dy*dy-2*c*q+2*H*(chi_tau-self.x*f['w']*chi_x)/N
        flux=2*inverse_A*chi_x
        zc=values[:self.modes]
        lap_z=self.B@(-self.k**2*zc)
        zx_over_x=np.divide(zx,self.x,out=lap_z/3,where=self.x!=0)
        radial=lap_z-zx_over_x-(nx*zx+f['lapse_x_over_x'])/N
        hdelta=hd(f)
        geometric_mean=-2*P*(H*(h_tau-self.x*f['w']*hx)/N
                              +(H+hdelta)*H*s+inverse_A**2*radial)
        geometric_flux=-2*P*inverse_A*(hx+H*s*zx+H*s_over_x)
        areal=R*self.x*np.exp(tau+z);u=self.x*zx
        compactness=areal**2*(2*H*hdelta+hdelta*hdelta)-2*u-u*u
        return dict(density_excess=rho,null_plus=mean+2*flux,null_minus=mean-2*flux,
            geometric_null_plus=geometric_mean+2*geometric_flux,
            geometric_null_minus=geometric_mean-2*geometric_flux,
            flux=flux,compactness_excess=compactness,areal_derivative_factor=1+u,
            maximum_lapse_residual=f['maximum_lapse_residual'],minimum_lapse=float(min(N)),
            minimum_lapse_slope=f['minimum_lapse_slope'])


def symbolic_checks():
    import sympy as s
    P,c,y,m,X=s.symbols('P c y m X',positive=True)
    D=P+3*c*y;k=m/D;N=1/s.sqrt(y)
    F=(P+9*c*y)*m*m/(3*D*D)-s.Rational(3,4)*X*X*y*y+X*y/2
    Fy=s.diff(F,y);C=k*(P+9*c*y)/(3*D)
    am=N**3*s.diff(F,m)/(2*Fy);b=-P*N**3/Fy
    B=-N*c*y/(P*D)
    symbol=s.Matrix([[C*b,B+C*am],[P*b,P*am]])
    return dict(principal_trace_cancellation=s.simplify(symbol.trace())==0,
        principal_quartic_coefficient=s.simplify(symbol.det()+P*c*N*N/(D*Fy))==0,
        lapse_derivative=s.simplify(Fy-(c*m*m*(P-9*c*y)/D**3-s.Rational(3,2)*X*X*y+X/2))==0)


def compare_states(model,state,other_model,other_state):
    x=np.linspace(0,8,801)
    B=spherical_jn(0,x[:,None]*model.k)
    other_B=spherical_jn(0,x[:,None]*other_model.k)
    values=state['coefficients'].reshape(-1,2,model.modes)@B.T
    other_values=other_state['coefficients'].reshape(-1,2,other_model.modes)@other_B.T
    return dict(zeta=float(np.max(abs(values[:,0]-other_values[:,0]))),
                mu=float(np.max(abs(values[:,1]-other_values[:,1]))))


def plot(result,state):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12.,4.3),layout='constrained')
    for i in (0,5,10):
        axes[0].plot(state['x'],state['null_plus'][i]*1e6,label=f"Ht={state['tau'][i]:.2f}")
    axes[0].axhline(0,color='black',lw=.7)
    axes[0].set(xlim=(0,5),xlabel='Comoving radius / initial width',
                ylabel='Null stress / cutoff^4 (x 1,000,000)',title='Outward negative null stress persists')
    axes[0].legend(fontsize=8)
    axes[1].plot(state['tau'],state['null_plus'][:,0]*1e6,color='#205e80')
    axes[1].set(xlabel='Expansion time Ht',ylabel='Central null stress (x 1,000,000)',
                title='The core deepens, then eases')
    axes[2].plot(state['tau'],np.min(state['areal_derivative_factor'],axis=1),color='#b85432')
    axes[2].axhline(0,color='black',ls='--',label='Necessary for a spatial neck')
    axes[2].set(xlabel='Expansion time Ht',ylabel='Minimum areal-radius derivative factor',ylim=(-.05,1.1),
                title='No spatial neck at saved times')
    axes[2].legend(fontsize=7)
    fig.suptitle('Short nonlinear evolution of the candidate field',fontsize=13)
    fig.supxlabel('Finite basis and domain, Ht up to 0.5. No long-term stability, available material, traversable bridge or backward passage established.',fontsize=8)
    for extension in ('png','pdf'):
        fig.savefig(project_path('nonlinear_packet_evolution_progress.'+extension),dpi=160)
    plt.close(fig)


def run():
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    baseline_path=project_path('nonlinear_packet_constraints_results.json')
    baseline=json.loads(baseline_path.read_text(encoding='utf-8'))
    if (not baseline['passed'] or baseline['source_sha256']!=sha(project_path('nonlinear_packet_constraints.py'))
            or baseline['state_sha256']!=sha(project_path('nonlinear_packet_constraints_state.npz'))):
        raise ArithmeticError('Nonlinear initial-state provenance failed.')
    model=RadialModel(48);state=model.evolve()
    controls={}
    for name,modes,radius,factor,kwargs in [
        ('coarse_basis',32,12.,4,{}),('finer_basis',64,12.,4,{}),
        ('larger_domain',72,18.,4,{}),('quadrature',48,12.,6,{}),
        ('time_accuracy',48,12.,4,dict(rtol=2e-11,atol=2e-14,max_step=.0025))]:
        other=RadialModel(modes,radius,factor);trial=other.evolve(**kwargs)
        controls[name]=compare_states(model,state,other,trial)
        controls[name]['rhs_evaluations']=trial['nfev']
        print(json.dumps(dict(completed=name,changes=controls[name])),flush=True)
    rows=[model.stress(t,v) for t,v in zip(state['tau'],state['coefficients'])]
    arrays={name:np.stack([r[name] for r in rows]) for name in rows[0] if isinstance(rows[0][name],np.ndarray)}
    scales=np.maximum(np.max(abs(arrays['null_plus']),axis=1),np.max(abs(arrays['null_minus']),axis=1))
    geometry_error=max(float(np.max(abs(arrays[name]-arrays['geometric_'+name])/scales[:,None]))
                       for name in ('null_plus','null_minus'))
    derivative_error=0.
    for i in (0,5,10):
        finer=model.stress(state['tau'][i],state['coefficients'][i],difference_step=5e-5)
        derivative_error=max(derivative_error,float(np.max(abs(finer['null_plus']-rows[i]['null_plus']))/scales[i]))
    # Independent known linear evolution, using the established mode solver.
    from localized_expanding_packet import packet
    from expanding_field_completion import coefficients
    weak=model.evolve(epsilon=1e-7)
    reference=packet(model.base['H'],epsilon=1e-7)
    x=np.linspace(0,8,401);j=spherical_jn(0,reference['q'][:,None]*x)
    mode=-1e-7*reference['weights']*reference['mode_amplitude'][10]
    velocity=-1e-7*reference['weights']*model.base['H']*reference['mode_velocity'][10]
    d=-(reference['q']/model.base['radius'])**2*np.exp(-1)*mode
    K=coefficients(.01,model.base['H'])['kinetic_over_cutoff2']
    linear_z=mode@j
    linear_mu=((-K*velocity+model.base['J']*d)/(model.base['P']*model.base['H']))@j
    B=spherical_jn(0,x[:,None]*model.k)
    weak_z=B@weak['coefficients'][-1,:model.modes]
    weak_mu=B@weak['coefficients'][-1,model.modes:]
    linear_errors=dict(zeta=float(np.max(abs(weak_z-linear_z))/np.max(abs(linear_z))),
                       mu=float(np.max(abs(weak_mu-linear_mu))/np.max(abs(linear_mu))))
    initial_values,old=model.initial()
    reconstructed=model.fields(0.,initial_values)
    initial_errors=dict(lapse=float(np.max(abs(reconstructed['lapse']-old['lapse']))),
                        shift=float(np.max(abs(reconstructed['w']-old['w']))),
                        zeta_tau=float(np.max(abs(model.B@model.rhs(0.,initial_values)[:model.modes]))))
    symbolic=symbolic_checks()
    checks=dict(**symbolic,
        homogeneous_control=bool(np.max(abs(model.rhs(0.,np.zeros(2*model.modes))))==0),
        initial_slice_reconstructed=max(initial_errors.values())<3e-10,
        basis_refinement=controls['finer_basis']['zeta']<2e-10 and controls['finer_basis']['mu']<2e-8,
        domain_refinement=controls['larger_domain']['zeta']<2e-10 and controls['larger_domain']['mu']<2e-8,
        nonlinear_quadrature=controls['quadrature']['zeta']<2e-11 and controls['quadrature']['mu']<2e-9,
        time_refinement=controls['time_accuracy']['zeta']<2e-11 and controls['time_accuracy']['mu']<2e-9,
        linear_limit=max(linear_errors.values())<1e-4,
        independent_geometric_null_stress=geometry_error<2e-6,
        stress_time_derivative=derivative_error<2e-7,
        lapse_constraint=max(r['maximum_lapse_residual'] for r in rows)<2e-12,
        continued_lapse_branch=min(r['minimum_lapse'] for r in rows)>.99
                              and min(r['minimum_lapse_slope'] for r in rows)>.9,
        regular_center_flux=bool(np.max(abs(arrays['flux'][:,0]))==0),
        negative_central_null_stress=bool(np.all(arrays['null_plus'][:,0]<0)),
        no_sampled_neck=bool(np.min(arrays['areal_derivative_factor'])>.999))
    summary=dict(initial_central_null=float(arrays['null_plus'][0,0]),final_central_null=float(arrays['null_plus'][-1,0]),
        minimum_saved_central_null=float(np.min(arrays['null_plus'][:,0])),
        time_of_saved_central_minimum=float(state['tau'][np.argmin(arrays['null_plus'][:,0])]),
        final_central_zeta=float(sum(state['coefficients'][-1,:model.modes])),
        largest_absolute_compactness=float(np.max(abs(arrays['compactness_excess']))),
        smallest_sampled_areal_derivative=float(np.min(arrays['areal_derivative_factor'])),
        minimum_lapse=min(r['minimum_lapse'] for r in rows),maximum_geometric_null_relative_error=geometry_error,
        maximum_stress_time_derivative_relative_change=derivative_error,
        initial_reconstruction_errors=initial_errors,linear_limit_relative_errors=linear_errors,
        evolution_rhs_evaluations=state['nfev'])
    # |v*j0'(v)|<=2 gives an all-space bound for every saved coefficient state.
    derivative_bound=2*np.sum(abs(state['coefficients'][:,:model.modes]),axis=1)
    summary['all_space_neck_exclusion_bound_at_saved_times']=float(max(derivative_bound))
    checks['no_neck_bound_at_saved_times']=bool(np.max(derivative_bound)<1)
    saved=dict(x=model.x,tau=state['tau'],k=model.k,coefficients=state['coefficients'],**arrays)
    state_path=project_path('nonlinear_packet_evolution_state.npz');np.savez_compressed(state_path,**saved)
    result=dict(complete=True,passed=all(checks.values()),checks=checks,summary=summary,refinement_changes=controls,
        parameters=dict(epsilon=.001,modes=48,xmax=12.,time_Ht=[0.,.5],quadrature_factor=4,
                        largest_retained_initial_wavenumber_over_cutoff=float(model.k[-1]/model.base['radius'])),
        source_sha256=sha(__file__),state_sha256=sha(state_path),
        inputs_sha256={name:sha(project_path(name)) for name in ('nonlinear_packet_constraints.py',
            'nonlinear_packet_constraints_results.json','expanding_field_completion.py','expanding_field_completion_results.json',
            'localized_expanding_packet.py','project_paths.py')},
        scope=['Finite radial Bessel Galerkin evolution with numerical basis, quadrature, domain and time refinement. No continuum error enclosure.',
               'Scalars vanish at a finite outer radius; auxiliary lapse and shift use an approximate asymptotic tail. Domain refinement is a check, not an exact infinite-domain construction.',
               'The full nonlinear matter null stress is compared with independent geometric projections of the evolving fields.',
               'The negative-null-energy region is an assumed EFT source. No material production, UV-validity bound, long-term stability, throat, safe transit or chronology change is established.',
               'Failed nodal-grid pilots are preserved; their grid-sensitive branch failures were not treated as physical singularities.',
               'An all-space derivative bound excludes a spatial neck for the finite-basis states at saved times; unsaved times and the exact nonlinear continuum are not certified.'])
    project_path('nonlinear_packet_evolution_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],checks=checks,summary=summary),indent=2),flush=True)
    if not result['passed']:
        raise ArithmeticError('Nonlinear evolution checks failed; saved evidence retained.')
    plot(result,saved)
    return result


if __name__=='__main__':
    run()
