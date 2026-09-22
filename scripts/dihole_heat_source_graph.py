"""One-derivative HS heat error and the actual positive-time source insertion.

All norms are on the unchanged full physical operator. HS quantities use
the static trace per Euclidean time; no derivative crosses a time window.
"""
import json
import math

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext
from dihole_full_frequency import multiply
from dihole_heat_defect import unpack
from dihole_local_transport import packed
from dihole_localized_heat_error import enclosure
from dihole_regulator_heat_tail import heat_multiplier
from dihole_slice_channel import sha
from project_paths import project_path


NAME = 'dihole_heat_source_graph'


def filter_HS(iv, row):
    tail = row['tail']
    H0, S = ({int(k):iv.mpf(v) for k,v in tail[name].items()}
             for name in ('H0_polynomial','heat_polynomial'))
    poly = multiply(S,multiply(H0,H0))
    if min(poly)<=1: raise ArithmeticError('Nonintegrable frequency tail.')
    K = iv.mpf(row['selected']['cutoff_K_D'])
    low = K*iv.mpf(row['selected_low_detail']['scalar_HS_upper'])**2/iv.pi
    high = sum(v*K**(1-p)/(p-1) for p,v in poly.items())/iv.pi
    return low+high, dict(low=packed(iv,low),high=packed(iv,high),
        total=packed(iv,low+high),minimum_tail_power=min(poly))


def heat_HS(iv, V, nu, mass, t, order):
    """||D^order H(t) Pi_buffer||HS,/T via the saved comparison filter."""
    spectral = iv.mpf(0)
    for k,c in enumerate((iv.mpf(1),iv.mpf(11)/6,iv.mpf(1),iv.mpf(1)/6)):
        p = order+k
        spectral += c/mass**(2*k)*(1 if p==0 else (iv.mpf(p)/(2*iv.e*t))**p)
    return iv.sqrt(nu*V*spectral)


def geometry(iv, frozen, prior, frame, local, buffer, q):
    c = {k:unpack(iv,v) for k,v in frozen['constants'].items()}
    H = iv.mpf(local['real_geometry']['inverse_metric_norm_upper'])
    gamma0,gamma1 = (unpack(iv,v) for v in local['normal_geometry']['christoffel_component_upper'][:2])
    Rabs = 128*H*gamma1+512*H*gamma0**2
    Omega = 48*H*iv.mpf(frame['frame_bounds']['coframe'][1])
    F,F1 = map(iv.mpf,prior['jets']['mask_derivative_upper'][:2])
    gradF = F1/iv.mpf(buffer['P_lower'])
    C4 = Rabs/4+abs(q)*c['magnetic_per_q']
    a = iv.sqrt(2)*F
    b = 2*iv.sqrt(2)*gradF+a*iv.sqrt(C4)+gradF/2+4*F*Omega
    exported = {k:packed(iv,v) for k,v in dict(R_E_absolute=Rabs,
        orthonormal_connection=Omega,F=F,grad_F=gradF,C4=C4,
        insertion_graph_a=a,insertion_graph_b=b).items()}
    return c,a,b,exported


def gaussian_terms(iv,c,row,local):
    """Coefficients times r^k t^p before L2 moments; cutoff kept separate."""
    a,alpha,b0 = (c[k] for k in ('a','alpha','b0'))
    v = {k:unpack(iv,x) for k,x in row['coefficients'].items()}
    u = iv.mpf(1)/100; cc = 1+u*u/3; ell = cc+u
    G = c['metric_difference']; g,h = c['cutoff_gradient'],c['cutoff_hessian']
    interior = [(3,-2,G*(cc*cc+u*u)/4),(1,-1,2*G*cc+a*ell*(2*alpha+b0)),
        (0,0,v['V0']),(3,-1,a*ell*v['b2']),(1,0,v['V1']),
        (2,0,v['V2']),(4,0,v['V4'])]
    cutoff = [(1,-1,a*a*g*ell),(0,0,4*a*a*h+2*a*g*(2*alpha+b0)),(2,0,2*a*g*v['b2'])]
    common = 2*iv.sqrt(iv.pi/5*iv.mpf(local['real_geometry']['volume_upper']))*iv.exp(u)
    moment = lambda k:iv.mpf(2)**(iv.mpf(k)/2-3)*iv.sqrt(math.factorial(k+1))/iv.pi
    def convert(rows,annulus):
        return [(iv.mpf(k)/2-1+p,common*x*moment(k)*(iv.mpf(2)**(iv.mpf(k)/2+1) if annulus else 1))
                for k,p,x in rows]
    return convert(interior,False),convert(cutoff,True),common,moment,ell,v


def calculate(iv,frozen,row,heat,filter_row,prior,frame,local,buffer,fraction):
    if not 0<fraction<=1: raise ValueError('Heat-time fraction must lie in(0,1].')
    q = row['q']; s0 = iv.mpf(heat['cutoff_s0_over_D_squared']); T = s0*fraction; t = T/2
    c,a,b,geo = geometry(iv,frozen,prior,frame,local,buffer,q)
    V,Vrecord = filter_HS(iv,filter_row)
    mass,nu = iv.mpf(heat['comparison_filter_Lambda_D']),iv.mpf(heat['inverse_lapse_upper'])
    inside,cut,common,moment,ell,v = gaussian_terms(iv,c,row,local)
    H = iv.mpf(local['real_geometry']['inverse_metric_norm_upper'])
    beta = iv.mpf(2)**-64/(32*H)
    peak = lambda p:(t/2)**p if p.a<0 else t**p
    late_HS = sum(z*peak(p) for p,z in inside)+iv.exp(-beta/t)*sum(z*peak(p) for p,z in cut)
    early_op = unpack(iv,enclosure(iv,frozen,row,fraction/4)['localized_global_heat_error_upper'])
    F0early = heat_HS(iv,V,nu,mass,t/2,0)*early_op
    F1early = heat_HS(iv,V,nu,mass,t/2,1)*early_op
    F0late,F1late = t*late_HS/2,iv.sqrt(t/iv.e)*late_HS
    F0,F1 = F0early+F0late,F1early+F1late
    Q0 = common*moment(0)/t
    Q1 = common*(c['a']*ell*moment(1)/(2*t**iv.mpf('1.5'))
        +c['b0']*moment(0)/t+v['b2']*moment(2)
        +2*c['a']*c['cutoff_gradient']*moment(0)/t*iv.exp(-beta/t))
    direct = [heat_HS(iv,V,nu,mass,t,j) for j in (0,1)]
    corrected = [Q0+F0,Q1+F1]
    X = [x if x.b<=y.b else y for x,y in zip(direct,corrected)]
    source = 2*T*X[1]*(a*X[1]+b*X[0])
    old = 2*T*nu*iv.mpf(heat['W']['complete_upper_D2'])*heat_multiplier(iv,mass,T)
    selected = source if source.b<old.b else old
    # Independent one-Dirac operator bound, using the uniform R4 majorant.
    eps = c['epsilon']; rootR4 = ((1+eps)/(1-eps))**2
    weights = dict(leading=iv.pi/2,constant=iv.mpf(2),varying_field=3*iv.pi/4,
        quadratic=iv.mpf(8)/3,quartic=iv.mpf(16)/5)
    powers = dict(leading=1,constant=2,varying_field=3,quadratic=4,quartic=6)
    f = fraction/2
    grad_sum = sum(unpack(iv,x)*f**(iv.mpf(powers[k])/2)*weights[k] for k,x in row['interior_terms'].items())
    cutweights = dict(gradient=iv.pi/2,constant=iv.mpf(2),quadratic=iv.mpf(8)/3)
    cutpowers = dict(gradient=1,constant=2,quadratic=4)
    grad_sum += iv.exp(-beta/s0*(1/f-1))*sum(unpack(iv,x)*f**(iv.mpf(cutpowers[k])/2)*cutweights[k]
        for k,x in row['cutoff_terms'].items())
    normalized_graph = rootR4/iv.sqrt(iv.e)*grad_sum
    return dict(q=q,heat_time_fraction=float(fraction.a),leg_time_fraction=float(f.a),
        geometry=geo,comparison_filter_HS_squared=Vrecord,
        normalized_Dirac_operator_error_upper=packed(iv,normalized_graph),
        late_defect_HS_sup_upper=packed(iv,late_HS),early_defect_operator_integral_upper=packed(iv,early_op),
        HS_error_terms={k:packed(iv,x) for k,x in dict(value_early=F0early,value_late=F0late,
            Dirac_early=F1early,Dirac_late=F1late).items()},
        value_HS_error_upper=packed(iv,F0),Dirac_HS_error_upper=packed(iv,F1),
        frozen_value_HS_upper=packed(iv,Q0),frozen_Dirac_HS_upper=packed(iv,Q1),
        direct_heat_HS_upper=[packed(iv,x) for x in direct],
        corrected_heat_HS_upper=[packed(iv,x) for x in corrected],
        selected_heat_HS_upper=[packed(iv,x) for x in X],
        physical_heat_variation_absolute_upper=packed(iv,source),
        previous_heat_variation_upper=packed(iv,old),selected_heat_variation_upper=packed(iv,selected),
        source_ceiling_improves=bool(source.b<old.a),
        source_ceiling_log10_change_upper=packed(iv,iv.log(old/selected,10)) if selected.b<old.a else packed(iv,0),
        signed_source=None,subtracted_short_time_remainder=None,renormalized_source=None,atom_transit=None)


def controls():
    mp.dps = 80
    moment_errors=[]; beta_errors=[]
    for k in range(5):
        exact=2**(mp.mpf(k)/2-3)*mp.sqrt(mp.factorial(k+1))/mp.pi
        integral=2*mp.pi**2*mp.quad(lambda r:r**(2*k+3)*mp.exp(-r*r/2)/(4*mp.pi)**4,[0,1,4,mp.inf])
        moment_errors.append(abs(mp.sqrt(integral)/exact-1))
    for p,w in zip(map(mp.mpf,('-.5','0','.5','1','2')),
                   (mp.pi/2,2,3*mp.pi/4,mp.mpf(8)/3,mp.mpf(16)/5)):
        val=2*(p+1)*mp.quad(lambda z:mp.sin(z)**(2*p+1),[0,mp.pi/2])
        beta_errors.append(abs(val-w))
    P=mp.matrix([[3,mp.mpf('.3'),mp.mpf('.4')],[mp.mpf('.3'),4,mp.mpf('.2')],[mp.mpf('.4'),mp.mpf('.2'),5]])
    lam,U=mp.eigsy(P); S=U*mp.diag([mp.sqrt(x) for x in lam])*U.T
    A=mp.diag([2,3]);J=mp.matrix([[1,0],[0,1],[0,0]]);Pi=J*J.T
    V=mp.matrix([[1,mp.mpf('.2'),0],[mp.mpf('.2'),-1,0],[0,0,0]])
    B=P*J-J*A;ratios=[];source_ratios=[];leakage=[]
    for t in map(mp.mpf,('.01','.1','.3')):
        H=lambda s:U*mp.diag([mp.exp(-s*x) for x in lam])*U.T
        R=H(t)*J-J*mp.expm(-t*A)
        norms=[mp.norm(R),mp.norm(S*R)]
        ceilings=[mp.norm(H(t/2))*t*mp.norm(B)/2+t*mp.norm(B)/2,
            mp.norm(S*H(t/2))*t*mp.norm(B)/2+mp.sqrt(t/mp.e)*mp.norm(B)]
        ratios.extend(float(x/y) for x,y in zip(norms,ceilings))
        HH=H(t/2);variation=-t*mp.fsum(((S*V+V*S)*H(t))[i,i] for i in range(3))
        bound=2*t*mp.norm(Pi*S*HH)*mp.norm(V)*mp.norm(Pi*HH)
        source_ratios.append(float(abs(variation)/bound));leakage.append(float(mp.norm((H(t)*J)[2,:])))
    return dict(Gaussian_L2_moment_relative_error=float(max(moment_errors)),
        Beta_weight_error=float(max(beta_errors)),split_HS_max_ratio=max(ratios),
        physical_insertion_max_ratio=max(source_ratios),exterior_leakage_min=min(leakage))


def run():
    iv=MPIntervalContext();iv.dps=100
    path=project_path('dihole_localized_heat_error_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_localized_heat_error.py':previous['source_sha256']}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):
        raise ArithmeticError('Exact passing188 inputs required.')
    frozen,prior,frame,local,buffer,heat,frequency=(json.loads(project_path(n+'_results.json').read_text()) for n in
        ('dihole_frozen_magnetic','dihole_regular_jets','dihole_positive_frame','dihole_local_transport',
         'dihole_buffer_smoothing','dihole_regulator_heat_tail','dihole_full_frequency'))
    records=[]
    for row,hr in zip(frozen['records'],heat['records']):
        if row['q']!=hr['q']:raise ArithmeticError('Flux mismatch.')
        fr=next(x for x in frequency['records'] if x['q']==row['q'] and x['eta_ratio_interval']==[0.,0.])
        for n in (1,2,4,16): records.append(calculate(iv,frozen,row,hr,fr,prior,frame,local,buffer['geometry'],iv.mpf(1)/n))
    control=controls()
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        Gaussian_L2_moments=control['Gaussian_L2_moment_relative_error']<1e-65,
        Beta_weights=control['Beta_weight_error']<1e-65,
        split_HS_control=control['split_HS_max_ratio']<1,
        actual_insertion_control=control['physical_insertion_max_ratio']<1,
        exterior_coupling_present=control['exterior_leakage_min']>0,
        frequency_tail_integrable=all(r['comparison_filter_HS_squared']['minimum_tail_power']>1 for r in records),
        selected_source_not_worse=all(unpack(iv,r['selected_heat_variation_upper']).b<=unpack(iv,r['previous_heat_variation_upper']).b for r in records),
        positive_time_only=all(r['subtracted_short_time_remainder'] is None and r['renormalized_source'] is None for r in records),
        no_signed_source_or_transport=all(r['signed_source'] is None and r['atom_transit'] is None for r in records))
    result=dict(passed=all(checks.values()),checks=checks,controls=control,records=records,
        input_sha256=inputs,source_sha256=sha(__file__),bound_encoding=previous['bound_encoding'],
        scope=['Global one-Dirac graph and per-Euclidean-time Hilbert-Schmidt heat errors; no exterior truncation.',
            'The saved compact metric probe and fixed gauge/volume identification are retained in dotD and dotP.',
            'An absolute bound on the actual positive-time heat variation is evaluated, not its sign or a signed frozen-reference error.',
            'The subtracted short-time remainder, physical regulator matching and renormalized supporting source remain unresolved.',
            'One hydrogen atom, trial1nm opening. No stability, atom transit, survival or clock gain is established.'],
        sources=['https://arxiv.org/abs/math/0303095','https://web.math.ku.dk/grubb/notes/heat.pdf'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Heat-source graph checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,controls=control,records=[dict(q=r['q'],fraction=r['heat_time_fraction'],
        source_log=r['physical_heat_variation_absolute_upper']['log10_upper'],old_log=r['previous_heat_variation_upper']['log10_upper'],
        improves=r['source_ceiling_improves']) for r in records]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=result['records'][4]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.bar(['Previous full-field ceiling','Derivative + HS estimate'],
        [r['previous_heat_variation_upper']['log10_upper'],r['physical_heat_variation_absolute_upper']['log10_upper']],color=['#8b5897','#386e94'])
    ax.set(ylabel='log10 of absolute heat-variation ceiling (D = 1)',title='Atomic parameters: bounding the actual field insertion')
    fig.supxlabel('Positive-time absolute ceilings, not a signed or renormalized supporting stress.\nNo stable opening, atom passage or clock gain is established.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
