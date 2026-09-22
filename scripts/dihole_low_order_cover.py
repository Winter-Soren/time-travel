"""Directed low-order regular-frame cover for the same frozen heat reference.

All source/output buffers used by192 lie inside u in [-.0025,.0525].
R in [1/40,21/400], R' in [0,1], and the saved R'' ceiling cover that region.
Rotational invariance reduces tensor Frobenius norms to y=(sqrt(t),0).
"""
import json

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext,calibration,field_function,lower,upper
from dihole_frozen_magnetic import first_defect
from dihole_heat_defect import unpack
from dihole_heat_source_graph import calculate
from dihole_local_transport import packed
from dihole_nested_heat_bound import DEPTHS,comparison
from dihole_pointwise_dirac import PointwiseDirac
from dihole_regular_low_fields import build_evaluator
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_low_order_cover'


def scalar_jets(c,raw,power_f,power_w):
    f,w=raw[:2]
    fr,wr,fx,wx=raw[3],raw[4],raw[6],raw[7]
    rr=power_f*fr/f+power_w*wr/w;xx=power_f*fx/f+power_w*wx/w
    value=f**power_f*w**power_w
    return (value,value*rr,value*xx,
        value*(rr**2+power_f*(raw[9]/f-(fr/f)**2)+power_w*(raw[10]/w-(wr/w)**2)),
        value*(rr*xx+power_f*(raw[12]/f-fr*fx/f**2)+power_w*(raw[13]/w-wr*wx/w**2)),
        value*(xx**2+power_f*(raw[15]/f-(fx/f)**2)+power_w*(raw[16]/w-(wx/w)**2)))


def tensor_norms(c,R,t,pole,rate,second,raw,regular):
    """Norms in the actual frame, before source/output comparison factors."""
    half=c.mpf(1)/2
    N,NR,Nx,_,_,_=scalar_jets(c,raw,half,-half)
    J,JR,Jx,JRR,JRx,Jxx=scalar_jets(c,raw,-half,-half)
    Ni,NiR,Nix,_,_,_=scalar_jets(c,raw,-half,half)
    h=c.sqrt(raw[0]);hR=raw[3]/(2*h);hx=raw[6]/(2*h)
    z,zR,zx=(regular[k] for k in ('zeta','zeta_R','zeta_x'))
    Z=z*(h+1)*J
    ZR=zR*(h+1)*J+z*hR*J+z*(h+1)*JR
    Zx=zx*(h+1)*J+z*hx*J+z*(h+1)*Jx
    d=-4*Z/(1+t)**2;dR=-4*ZR/(1+t)**2
    dt=8*Z/(1+t)**3+8*pole*Zx/(1+t)**4
    xt=-2*pole/(1+t)**2;y=c.sqrt(t)
    a=(1+t)*N/(2*R);b=(1+t)*d/(2*R);rad=(1+t)*J/(2*R)
    aR=(1+t)*(NR-N/R)/(2*R);at=(N+(1+t)*Nx*xt)/(2*R)
    radR=(1+t)*(JR-J/R)/(2*R);radt=(J+(1+t)*Jx*xt)/(2*R)
    diagonal_u=[rate*NiR/Ni,rate*JR/J,rate*radR/rad,rate*aR/a]
    diagonal_y=[2*y*xt*Nix/Ni,2*y*xt*Jx/J,2*y*radt/rad,2*y*at/a]
    frame=J**2*sum(v**2 for v in diagonal_u)+rad**2*sum(v**2 for v in diagonal_y)
    frame+=a**2*((b*y/a)**2+(b*y/rad)**2)
    numerator=N-2*pole*Jx/(1+t)+(1+3*t)*d/2
    K=numerator/R
    KR=(NR-2*pole*JRx/(1+t)+(1+3*t)*dR/2-numerator/R)/R
    Kt=(Nx*xt-2*pole*(Jxx*xt/(1+t)-Jx/(1+t)**2)
        +(3*d+(1+3*t)*dt)/2)/R
    zero=((rate*JR)**2+t*K**2)/4
    zero_gradient=J**2*((second*JR+rate**2*JRR)**2+(rate*KR*y)**2)/4
    zero_gradient+=rad**2*((rate*JRx*y*xt)**2+(K+2*t*Kt)**2/4)+a**2*K**2/4
    G,GR,Gx=(regular[k] for k in ('G','G_R','G_x'))
    C=4*pole*rate*G/(1+t)**2
    Cu=4*pole*(second*G+rate**2*GR)/(1+t)**2
    Ct=4*pole*rate*(Gx*xt/(1+t)**2-2*G/(1+t)**3)
    Du=-4*rate*raw[14]/(1+t)**2
    Dt=-4*raw[17]*xt/(1+t)**2+8*raw[8]/(1+t)**3
    # Each antisymmetric pair occurs twice in the full three-index norm.
    field=2*((J*J*a*Cu*y)**2+(J*rad*a*Du)**2
        +(rad*J*a*(C+2*t*Ct))**2+(rad*rad*a*2*y*Dt)**2+(a*J*rad*C)**2)
    return dict(frame_derivative=c.sqrt(frame),geometric_zero=c.sqrt(zero),
        geometric_zero_derivative=c.sqrt(zero_gradient),curvature_derivative=c.sqrt(field))


def cover(iv,prior,radial=16,angular=24):
    function=field_function(iv);regular=build_evaluator(iv);params=calibration(iv)
    caps=dict.fromkeys(('frame_derivative','geometric_zero','geometric_zero_derivative','curvature_derivative'),0.)
    guards={'NP':float('inf'),'W':float('inf')}
    second=iv.mpf(prior['step_derivative_upper'][1])*20
    for i in range(radial):
        R=iv.mpf([iv.mpf(1)/40+iv.mpf(11)*i/(400*radial),iv.mpf(1)/40+iv.mpf(11)*(i+1)/(400*radial)])
        for j in range(angular):
            t=iv.mpf([iv.mpf(2)*j/angular,iv.mpf(2)*(j+1)/angular])
            for pole in (-1,1):
                x=pole*(1-t)/(1+t)
                raw=function(R,x,*params)
                if raw[0].a<=0 or raw[1].a<=0:raise ArithmeticError('Metric positivity unresolved on a tile.')
                guards['NP']=min(guards['NP'],lower(raw[0]));guards['W']=min(guards['W'],lower(raw[1]))
                values=tensor_norms(iv,R,t,pole,iv.mpf([0,1]),iv.mpf([0,second]),raw,regular(R,x))
                for k,v in values.items():caps[k]=max(caps[k],upper(v))
    return dict(radial_tiles=radial,angular_tiles=angular,boxes=2*radial*angular,
        radius_interval=['1/40','21/400'],squared_chart_radius_interval=[0,2],
        radius_rate_interval=[0,1],radius_second_upper=upper(second),
        positive_guards=guards,actual_frame_tensor_upper=caps,
        coverage='Whole radial/chart product, both poles; no sampling substitutes for intervals.')


def improved_constants(iv,frozen,enclosed):
    old={k:unpack(iv,v) for k,v in frozen['constants'].items()}
    a,eps=old['a'],old['epsilon']
    v={k:iv.mpf(x) for k,x in enclosed['actual_frame_tensor_upper'].items()}
    alpha=a*iv.sqrt(1+eps)*v['frame_derivative']
    proposed=dict(alpha=alpha,metric_difference=2*a*alpha,b0=v['geometric_zero'],
        b1=iv.sqrt(1+eps)*v['geometric_zero_derivative'],
        curvature_gradient_per_q=(1+eps)**iv.mpf('1.5')*v['curvature_derivative'])
    updated=dict(old)
    for k,value in proposed.items():
        if value.b<old[k].b:updated[k]=value
    return old,updated,proposed


def point_controls():
    """Compare all norm contractions with numerical coordinate derivatives."""
    mp.dps=60;e=PointwiseDirac();field=field_function(mp);reg=build_evaluator(mp);params=calibration(mp)
    errors=[]
    for pole,y in ((1,mp.mpf('.4')),(-1,mp.mpf('.7'))):
        coords=[mp.mpf('.021'),y,mp.mpf(0)]
        p=e.point(*coords,pole=pole);R,rate,second=e.radial(coords[0]);x=pole*(1-y*y)/(1+y*y)
        predicted=tensor_norms(mp,R,y*y,pole,rate,second,field(R,x,*params),reg(R,x))
        derivatives=[];bg=[];curves=[]
        for nu in range(3):
            def at(z):
                v=list(coords);v[nu]=z
                return e.point(*v,pole=pole)
            derivatives.append(mp.matrix([[mp.diff(lambda z:at(z)['I'][a,b],coords[nu]) for b in range(4)] for a in range(4)]))
            bg.append(mp.matrix([mp.diff(lambda z:at(z)['divI'][a]/2,coords[nu]) for a in range(4)]))
            curves.append(mp.matrix([[mp.diff(lambda z:at(z)['F'][a,b],coords[nu]) for b in range(4)] for a in range(4)]))
        I,E=p['I'],p['E'];f=mp.mpf(0);g=mp.mpf(0);h=mp.mpf(0)
        for c in range(4):
            f+=mp.norm(sum((I[nu+1,c]*derivatives[nu]*E for nu in range(3)),mp.zeros(4)))**2
            g+=mp.norm(sum((I[nu+1,c]*bg[nu] for nu in range(3)),mp.zeros(4,1)))**2
            h+=mp.norm(I.T*sum((I[nu+1,c]*curves[nu] for nu in range(3)),mp.zeros(4))*I)**2
        expected=dict(frame_derivative=mp.sqrt(f),geometric_zero=mp.norm(p['divI'])/2,
            geometric_zero_derivative=mp.sqrt(g),curvature_derivative=mp.sqrt(h))
        errors.extend(abs(predicted[k]-v)/max(1,abs(v)) for k,v in expected.items())
    return dict(independent_tensor_relative_error=str(max(errors)))


def run():
    iv=MPIntervalContext();iv.dps=90
    path=project_path('dihole_nested_heat_bound_results.json');previous=json.loads(path.read_text())
    inputs={**previous['input_sha256'],path.name:sha(path),'dihole_nested_heat_bound.py':previous['source_sha256'],
        'dihole_regular_low_fields.py':sha(project_path('dihole_regular_low_fields.py'))}
    if not previous['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing192 inputs required.')
    names=('dihole_frozen_magnetic','dihole_regular_jets','dihole_positive_frame','dihole_local_transport',
        'dihole_buffer_smoothing','dihole_regulator_heat_tail','dihole_full_frequency')
    frozen,prior,frame,local,buffer,heat,frequency=(json.loads(project_path(n+'_results.json').read_text()) for n in names)
    enclosed=cover(iv,prior)
    old,new,proposed=improved_constants(iv,frozen,enclosed)
    adjusted=dict(frozen,constants={k:packed(iv,v) for k,v in new.items()})
    records=[]
    for hr in heat['records']:
        row=first_defect(iv,new,hr['q'],iv.mpf(hr['cutoff_s0_over_D_squared']))
        fr=next(x for x in frequency['records'] if x['q']==row['q'] and x['eta_ratio_interval']==[0.,0.])
        steps=[calculate(iv,adjusted,row,hr,fr,prior,frame,local,buffer['geometry'],iv.mpf(2)**-j) for j in range(13)]
        for exponent in (0,1,2,4):
            T=iv.mpf(hr['cutoff_s0_over_D_squared'])*iv.mpf(2)**-exponent
            trials=[comparison(iv,steps[exponent:],T,depth) for depth in DEPTHS]
            best=min(trials,key=lambda v:unpack(iv,v['signed_reference_error_per_q2_upper']).b)
            prior_row=next(v for v in previous['records'] if v['q']==hr['q'] and v['heat_time_fraction']==float(2**-exponent))
            records.append(dict(q=hr['q'],heat_time_fraction=float(2**-exponent),selected=best,
                previous_reference_error_per_q2_upper=prior_row['selected_reference_error_per_q2_upper'],
                first_defect_column_integral_upper=row['total_column_time_integral_upper'],
                signed_full_heat_interval=None,renormalized_source=None,atom_transit=None))
    control=point_controls()
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        independent_tensor_derivatives=mp.mpf(control['independent_tensor_relative_error'])<mp.mpf('1e-45'),
        whole_buffer_positive=min(enclosed['positive_guards'].values())>0,
        full_chart_and_both_poles=enclosed['boxes']==2*enclosed['radial_tiles']*enclosed['angular_tiles'],
        retained_radius_and_metric_comparison=new['epsilon']==old['epsilon'] and new['radius']==old['radius'],
        every_replacement_tighter=all(new[k].b<old[k].a for k in proposed),
        every_atomic_comparison_tighter=all(unpack(iv,r['selected']['signed_reference_error_per_q2_upper']).b
            <unpack(iv,r['previous_reference_error_per_q2_upper']).a for r in records if r['q']!=1000),
        physical_source_still_missing=all(r['signed_full_heat_interval'] is None and r['renormalized_source'] is None for r in records))
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,controls=control,cover=enclosed,
        previous_constants={k:packed(iv,old[k]) for k in proposed},
        selected_constants={k:packed(iv,new[k]) for k in proposed},records=records,
        source_sha256=sha(__file__),input_sha256=inputs,
        scope=['Directed real cover of low-order regular-frame tensors; axis rotation preserves full tensor norms.',
            'Only proved tighter coefficients replace the old defect constants; geometry, cutoff and physical field are unchanged.',
            'The191 signed reference remains numerical, not enclosed; no renormalized support or atom transit.'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Low-order cover checks failed; result retained.')
    plot(result,previous)
    print(json.dumps(dict(passed=True,checks=checks,cover=enclosed,atomic=[dict(fraction=r['heat_time_fraction'],
        radius_log=r['selected']['signed_reference_error_per_q2_upper']['log10_upper']) for r in records if r['q']!=1000]),indent=2))


def plot(result,previous):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    keys=list(result['selected_constants']);x=list(range(len(keys)))
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    old=[result['previous_constants'][k]['log10_upper'] for k in keys]
    new=[result['selected_constants'][k]['log10_upper'] for k in keys]
    axes[0].plot(x,old,'x--',color='#777777',label='Previous ceiling')
    axes[0].plot(x,new,'o-',color='#247886',label='Directed tensor cover')
    axes[0].set(xticks=x,xticklabels=['Frame','Metric','Zero term','Zero derivative','Field derivative'],
        ylabel='log10 of coefficient ceiling',title='Bound the actual low-order field tensors')
    axes[0].tick_params(axis='x',labelsize=8);axes[0].legend(fontsize=8)
    row=max((r for r in result['records'] if r['heat_time_fraction']==1),key=lambda v:v['q'])
    axes[1].bar(['Previous error','New error','Reference magnitude'],
        [row['previous_reference_error_per_q2_upper']['log10_upper'],
         row['selected']['signed_reference_error_per_q2_upper']['log10_upper'],
         float(mp.log10(abs(mp.mpf(previous['prior_numerical_reference_per_q2']))))],
        color=['#777777','#247886','#73488e'])
    axes[1].axhline(0,color='#888888',linewidth=.7)
    axes[1].set(ylabel='log10 of normalized magnitude',title='Comparison at the same atomic heat time')
    fig.supxlabel('Whole-region interval bounds; the physical field and opening are unchanged.\nThe reference quadrature is still uncertified. No supported opening or atom time travel is demonstrated.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
