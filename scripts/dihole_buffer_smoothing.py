"""Full-field weighted smoothing through buffered Dirichlet estimates.

The auxiliary boundary is a factor in an estimate, not a boundary imposed on
the physical field. The finite-regulator, low-frequency bound is not stress.
"""
import json
import math

import numpy as np
from mpmath import mp

from dihole_buffer_geometry import buffer_bounds
from dihole_channel_bounds import MPIntervalContext, lower, upper
from dihole_slice_channel import sha
from dirac_green_mass_limit import benchmark
from project_paths import project_path


NAME = 'dihole_buffer_smoothing'


def weighted_graph(iv, M, eta_ratio, radius, nu, nmax, weight):
    """Local ||1_Omega A_z^r K_z sqrtN||, r=0..3, on a complex disk."""
    if M.a <= 0 or eta_ratio.a < 0 or radius.a < 0 or (M-radius*nu).a <= 0:
        raise ValueError('Positive masses and disk radius below M/nu required.')
    t = [iv.sqrt(eta_ratio+j) for j in (1, 2, 3)]
    graph = []
    for r in range(4):
        total = iv.sqrt(6)*M**3*iv.sqrt(nmax) if r == 3 else iv.mpf(0)
        for j in range(3):
            # Rationalized mass gaps avoid subtracting overlapping intervals.
            coefficient = iv.sqrt(6)*M**(r+1)*t[j]**r
            for k in range(3):
                if k != j:
                    coefficient *= (t[k]+t[j])/abs(k-j)
            gap = M*t[j]-radius*nu
            total += coefficient*(iv.sqrt(nmax)/gap+weight/gap**2)
        graph.append(iv.mpf(upper(total)))
    return graph


def scalar_comparison(iv, geometry):
    get = lambda key: iv.mpf(geometry[key])
    nlo,nhi,plo,phi,rlo,rhi = map(get, ('N_lower','N_upper','P_lower','P_upper','R_lower','R_upper'))
    r0 = iv.mpf(1)/40
    density_lo = get('scalar_volume_ratio_lower')
    density_hi = get('scalar_volume_ratio_upper')
    metric = min(lower(x) for x in (1/phi**2,r0**2/(phi**2*rhi**2),r0**2*nlo**2/rhi**2))
    c = lower(density_lo/density_hi*iv.mpf(metric))
    if c <= 0:
        raise ArithmeticError('Scalar reference must remain elliptic.')
    return dict(reference_radius_over_D=lower(r0), interval_length_over_D=.09,
                density_lower=lower(density_lo),density_upper=upper(density_hi),
                scalar_eigenvalue_factor_lower=c,collars=2,spin_components=4)


def dirichlet_trace(iv, c, b):
    """Two collars, rank four; b=a^2-Cq>0; output has units D^4."""
    if c.a <= 0 or b.a <= 0:
        raise ValueError('Positive comparison factor and shifted mass required.')
    length,r0 = iv.mpf(9)/100,iv.mpf(1)/40
    return 2*(length/(iv.sqrt(c)*b**iv.mpf('1.5'))
              +2*length*r0/(c*b)+2*length*r0*r0/(c**iv.mpf('1.5')*iv.sqrt(b)))


def full_bound(iv, geometry, comparison, prior, selection, eta_cap):
    get = lambda key: iv.mpf(geometry[key])
    q = selection['q']; M = iv.mpf(selection['fixed_regulator_M_D'])
    split = iv.mpf(selection['selected_split']['split_K_D'])
    nu = iv.mpf(prior['inverse_lapse_upper'])
    cw = iv.mpf(prior['weight_upper_D'])
    cq = get('scalar_curvature_negative_upper_D2')/4+q*get('magnetic_unit_flux_upper_D2')+get('log_lapse_gradient_squared_upper_D2')/4
    cq = iv.mpf(upper(cq)); rootc = iv.sqrt(cq)
    g = iv.mpf(1600)/get('P_lower')
    lap = (iv.mpf(4800000)+1600*get('radial_laplacian_drift_upper_D'))/get('P_lower')**2
    c = iv.mpf(comparison['scalar_eigenvalue_factor_lower'])
    eta = iv.mpf([0,eta_cap]); radius = M/(2*nu)
    cn = 2*cw*nu**iv.mpf('1.5')
    perturbation = radius*cn+(2*split*radius+radius**2)*nu**2
    connection = get('radial_connection_upper_D'); pmax = get('P_upper')
    ell = iv.mpf(1)/100

    def estimate(disk, b):
        rho = radius if disk else iv.mpf(0)
        G = weighted_graph(iv,M,eta,rho,nu,get('N_upper'),cw)
        a2 = cq+b; difference = perturbation if disk else iv.mpf(0)
        S = iv.mpf(upper(dirichlet_trace(iv,c,b)))
        H = [G[r+2]+(a2+difference+lap+2*g*(2*g+rootc+rho*nu))*G[r]+2*g*G[r+1]
             for r in (0,1)]
        B0,B1 = (iv.sqrt(S)*h for h in H)
        derivative = pmax*(B1+(2*g+rootc+rho*nu)*B0)+connection*B0
        Z2 = B0**2+ell**2*derivative**2
        return dict(Z_squared_upper=upper(Z2),local_graph_upper=[upper(x) for x in G],
                    localized_HS_upper=[upper(B0),upper(B1)],
                    dirichlet_resolvent_HS_squared_upper_D4=upper(S),
                    auxiliary_b_D2=float(b.mid),auxiliary_b_D2_interval=[lower(b),upper(b)],
                    auxiliary_a_squared_upper_D2=upper(a2))

    # ponytail: a finite grid chooses a valid auxiliary estimate, not an optimum.
    scale = max(upper(M*M),upper(cq))
    choices = [iv.mpf(scale)*iv.mpf(f) for f in ('0.1','1','10')]
    real = min((estimate(False,b) for b in choices),key=lambda row:row['Z_squared_upper'])
    disk = min((estimate(True,b) for b in choices),key=lambda row:row['Z_squared_upper'])
    z2 = iv.mpf(real['Z_squared_upper']); z = iv.sqrt(z2)
    zp = iv.sqrt(iv.mpf(disk['Z_squared_upper']))/radius
    j0 = iv.mpf(prior['J_zero_upper_D']); j1 = iv.mpf(prior['J_frequency_coefficient_upper'])
    b0 = j0*z2
    bp = j1*z2+2*(j0+j1*split)*z*zp
    low = b0/2+split*bp/iv.pi
    return dict(q=q,regulator_M_D=float(M.mid),split_K_D=float(split.mid),
                eta_ratio_interval=[0,eta_cap],Cq_upper_D2=upper(cq),
                cutoff_gradient_upper_D=upper(g),cutoff_laplacian_upper_D2=upper(lap),
                cauchy_radius_lower_D=lower(radius),complex_square_perturbation_upper_D2=upper(perturbation),
                real=real,complex_disk=disk,full_Z_HS_squared_upper=upper(z2),
                full_Z_derivative_HS_upper_over_D=upper(zp),
                full_B_zero_trace_upper_D=upper(b0),full_B_derivative_trace_upper=upper(bp),
                full_paired_low_integral_absolute_upper_D=upper(low))


def algebra_controls():
    rows=[]; eye=np.eye(6); chi=np.diag([1.,1.,0.,1.,1.,0.])
    for largest in (10.,1000.):
        for kappa in (0.,.3,2.):
            _,beta,N,_,_,A,_,_=benchmark(largest,kappa)
            W=np.diag(np.sqrt(np.diag(N)));nu=np.linalg.norm(np.linalg.inv(N),2)
            cw=np.linalg.norm(A@W-W@A,2);gen=1j*beta@np.linalg.inv(N)
            for M in (.7,2.):
                eta=.01*M*M; masses=np.sqrt(eta+M*M*np.arange(1,4));alpha=np.sqrt(6)*M**3
                radius=M/(2*nu)
                for angle in (0.,.7,2.1):
                    delta=radius*np.exp(1j*angle);Az=A+delta*gen
                    resolvents=[np.linalg.inv(Az+m*eye) for m in masses]
                    K=alpha*resolvents[0]@resolvents[1]@resolvents[2]
                    iv=MPIntervalContext();iv.dps=50
                    bounds=weighted_graph(iv,iv.mpf(M),iv.mpf('.01'),iv.mpf(radius),iv.mpf(nu),iv.mpf(2),iv.mpf(cw))
                    for r in range(4):
                        action=np.linalg.matrix_power(Az,r)@K
                        partial=alpha*eye if r==3 else np.zeros_like(Az)
                        for j,m in enumerate(masses):
                            residue=alpha*(-m)**r/math.prod(masses[k]-m for k in range(3) if k!=j)
                            partial=partial+residue*resolvents[j]
                        rows.append(dict(largest_lapse=largest,kappa=kappa,M=M,angle=angle,power=r,
                            partial_fraction_error=float(np.linalg.norm(action-partial)/max(1,np.linalg.norm(action))),
                            weighted_graph_ratio=float(np.linalg.norm(chi@action@W,2)/float(bounds[r].b))))
    return rows


def cauchy_controls():
    rows=[];eye=np.eye(6);chi=np.diag([1.,1.,0.,1.,1.,0.])
    for largest in (10.,1000.):
        for kappa in (0.,.3):
            _,beta,N,_,_,A,_,_=benchmark(largest,kappa)
            W=np.diag(np.sqrt(np.diag(N)));gen=1j*beta@np.linalg.inv(N)
            nu=np.linalg.norm(gen,2);cw=np.linalg.norm(A@W-W@A,2)
            M=.7;radius=M/(2*nu);masses=M*np.sqrt(np.arange(1,4)+.01)
            alpha=np.sqrt(6)*M**3
            R=[np.linalg.inv(A+m*eye) for m in masses]
            product=np.linalg.multi_dot
            exact=-alpha*sum(product(R[:j]+[R[j],gen]+R[j:]) for j in range(3))
            contour=np.zeros_like(A);difference_errors=[];ratios=[]
            for angle in 2*np.pi*np.arange(64)/64:
                delta=radius*np.exp(1j*angle);Az=A+delta*gen
                value=alpha*product([np.linalg.inv(Az+m*eye) for m in masses])
                contour+=value/delta/64
                perturbation=-Az@Az+A@A
                analytic=-delta*(A@gen+gen@A)-delta**2*gen@gen
                difference_errors.append(np.linalg.norm(perturbation-analytic))
                ceiling=radius*2*cw*nu**1.5+(2*abs(kappa)*radius+radius**2)*nu**2
                ratios.append(np.linalg.norm(perturbation,2)/ceiling)
            iv=MPIntervalContext();iv.dps=50
            g0=weighted_graph(iv,iv.mpf(M),iv.mpf('.01'),iv.mpf(radius),iv.mpf(nu),iv.mpf(2),iv.mpf(cw))[0]
            rows.append(dict(largest_lapse=largest,kappa=kappa,
                derivative_cauchy_relative_error=float(np.linalg.norm(exact-contour)/max(1,np.linalg.norm(exact))),
                complex_square_error=float(max(difference_errors)),complex_square_bound_ratio=float(max(ratios)),
                cauchy_weighted_HS_ratio=float(np.linalg.norm(chi@exact@W,'fro')/(np.sqrt(6)*float(g0.b)/radius))))
    return rows


def heat_controls():
    mp.dps=65; rows=[]
    for value in ('0.01','0.1','1','10'):
        a=mp.mpf(value)
        cutoff=int(mp.ceil(mp.sqrt(200/a)))+5
        angular=mp.fsum((2*l+1)*mp.exp(-a*l*(l+1)) for l in range(cutoff))
        ceiling=1+1/a+mp.sqrt(mp.pi/a)
        rows.append(dict(t=value,angular_sum=float(angular),angular_ceiling=float(ceiling),enclosed=angular<ceiling))
    integrals=[]
    for c,b in (('0.2','0.7'),('0.0001','100')):
        c,b=mp.mpf(c),mp.mpf(b);L=mp.mpf(9)/100;r=mp.mpf(1)/40
        # t=u^2 removes the integrable origin singularity.
        numeric=mp.quad(lambda u:8*L/mp.sqrt(mp.pi*c)*mp.exp(-b*u*u)*(u*u+mp.sqrt(mp.pi/c)*r*u+r*r/c),[0,mp.inf])
        analytic=2*(L/(mp.sqrt(c)*b**mp.mpf('1.5'))+2*L*r/(c*b)+2*L*r*r/(c**mp.mpf('1.5')*mp.sqrt(b)))
        integrals.append(float(abs(numeric-analytic)/analytic))
    return dict(angular=rows,heat_integral_relative_errors=integrals)


def run():
    path=project_path('dihole_source_completion_results.json');source=json.loads(path.read_text())
    inputs=dict(source['input_sha256']);inputs.update({path.name:sha(path),'dihole_source_completion.py':source['source_sha256'],
        'dihole_buffer_geometry.py':sha(project_path('dihole_buffer_geometry.py'))})
    trial_name='dihole_buffer_smoothing_prerate_trial.json'
    trial=json.loads(project_path(trial_name).read_text())
    trial_sources=['provenance/dihole_buffer_geometry_prerate.py','provenance/dihole_buffer_smoothing_prerate.py']
    if (trial['source_sha256']!=sha(project_path(trial_sources[1])) or
            trial['input_sha256']['dihole_buffer_geometry.py']!=sha(project_path(trial_sources[0]))):
        raise ArithmeticError('Exact pre-correction reporting diagnostic required.')
    inputs.update({name:sha(project_path(name)) for name in [trial_name,*trial_sources]})
    if not source['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing source-completion inputs required.')
    prior=json.loads(project_path('dihole_paired_low_frequency_results.json').read_text())
    iv=MPIntervalContext();iv.dps=60
    geometry=buffer_bounds(iv);comparison=scalar_comparison(iv,geometry)
    constants={**prior['global_lapse'],**source['insertion_factor'],
               'weight_upper_D':source['global_weight']['global_commutator_upper_D']}
    records=[full_bound(iv,geometry,comparison,constants,selection,cap)
             for selection in (prior['records'][0],prior['records'][-1]) for cap in (0.,.01)]
    algebra=algebra_controls();cauchy=cauchy_controls();heat=heat_controls()
    tests=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        prerate_diagnostic_retained=trial['geometry']['intervals']['log_lapse_gradient_squared'][0]>0,
        constant_cylinder_in_derived_intervals=all(geometry['intervals'][key][0]<=0<=geometry['intervals'][key][1]
            for key in ('log_lapse_gradient_squared','radial_volume_connection','radial_laplacian_drift')),
        buffered_geometry=geometry['geometry_controls']['passed'],
        warped_curvature=max(geometry['geometry_controls']['warped_product_relative_errors'])<1e-10,
        flat_curvature=geometry['geometry_controls']['flat_spherical_scalar_error']<1e-10,
        pole_regularity=max(geometry['geometry_controls']['balanced_polar_value_derivative_errors'])<1e-9,
        positive_scalar_comparison=comparison['scalar_eigenvalue_factor_lower']>0,
        partial_fraction_algebra=max(r['partial_fraction_error'] for r in algebra)<1e-10,
        local_weighted_graph=max(r['weighted_graph_ratio'] for r in algebra)<=1+1e-12,
        cauchy_derivative=max(r['derivative_cauchy_relative_error'] for r in cauchy)<1e-10,
        complex_square_identity=max(r['complex_square_error'] for r in cauchy)<1e-10,
        complex_square_perturbation=max(r['complex_square_bound_ratio'] for r in cauchy)<=1+1e-12,
        cauchy_weighted_HS=max(r['cauchy_weighted_HS_ratio'] for r in cauchy)<=1+1e-12,
        angular_heat_envelope=all(r['enclosed'] for r in heat['angular']),
        heat_integral=max(heat['heat_integral_relative_errors'])<1e-55,
        positive_finite_full_bounds=all(math.isfinite(r['full_paired_low_integral_absolute_upper_D']) and r['full_paired_low_integral_absolute_upper_D']>0 for r in records),
        uniform_shift_bounds=all(records[i+1]['full_Z_HS_squared_upper']>=records[i]['full_Z_HS_squared_upper'] for i in (0,2)))
    tests={k:bool(v) for k,v in tests.items()}
    result=dict(passed=all(tests.values()),tests=tests,geometry=geometry,scalar_comparison=comparison,
        records=records,algebra_controls=algebra,cauchy_controls=cauchy,heat_controls=heat,source_sha256=sha(__file__),input_sha256=inputs,
        full_high_frequency_integral=None,renormalized_physical_source=None,source_sign=None,atom_transit=None,
        references=['https://arxiv.org/pdf/hep-th/0306138','https://arxiv.org/pdf/1105.0532','https://arxiv.org/pdf/1608.06798'],
        scope=['All angular states, both collars and exterior propagation retained in the weighted source estimate.',
               'Auxiliary Dirichlet inverse factors localized fields; it is not a substitute physical exterior.',
               'Finite positive-mass product and bounded paired frequency interval only. Common physical regulator removal is unresolved.',
               'No supporting stress sign, stability, atomic transmission, clock gain or new physical law is established.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Buffered smoothing checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,tests=tests,atomic_records=records[-2:]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    for lo,hi,label,y in ((0,.05,'Physical variation',4),(-.005,.055,'Jet cutoff',3),(-.01,.06,'Localization cutoff',2),(-.015,.065,'Gradient cutoff',1),(-.02,.07,'Auxiliary Dirichlet domain',0)):
        axes[0].plot([lo,hi],[y,y],lw=9,solid_capstyle='butt',label=label)
    axes[0].set(yticks=range(5),yticklabels=['Auxiliary domain','Gradient','Localization','Jet','Physical variation'],xlabel='Collar coordinate s / D',title='Buffers around each collar\nFull exterior remains in the operator')
    rows=result['records'][-2:]
    for row,label in zip(rows,('Zero artificial shift','Uniform shift interval')):
        axes[1].semilogy(range(4),[v/row['regulator_M_D']**r for r,v in enumerate(row['real']['local_graph_upper'])],'o-',label=label)
    axes[1].set(xticks=range(4),xlabel='Number of Dirac-operator factors r',ylabel='Upper bound divided by M^r',title='Weighted local graph bounds\nAtomic target')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.supxlabel('Mathematical estimate on a prescribed geometry. No opening, stress sign or atom passage demonstrated.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
