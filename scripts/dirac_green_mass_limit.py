"""Finite algebra controls for the localized Green identity, not bridge stress.

The six-dimensional benchmark respects the exact factorization. A compressed
physical form need not do so; no continuum limit is inferred from this test.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from project_paths import project_path


def benchmark(largest_lapse, kappa):
    K=np.array([[0.,1.,0.],[1.,0.,1.],[0.,1.,0.]])
    Z=np.zeros_like(K);I=np.eye(3)
    H=np.block([[K,Z],[Z,-K]]);beta=np.block([[Z,I],[I,Z]])
    N=np.diag([.2,2.,largest_lapse]*2)
    F=np.diag([.2,-.1,.3]*2);B=np.diag(1/np.sqrt(np.diag(N)))
    V0=np.array([[.3,.1,-.2],[.1,-.4,.05],[-.2,.05,.2]])
    V=np.block([[V0,Z],[Z,-V0]]);dN=N@F;dB=-F@B/2
    A=H+1j*kappa*np.eye(6)
    D=beta@B@A@B
    dD=beta@(dB@A@B+B@V@B+B@A@dB)
    P=D.conj().T@D;dP=dD.conj().T@D+D.conj().T@dD
    return H,beta,N,V,dN,D,P,dP


def run():
    eye=np.eye(6);errors={k:[] for k in ('inverse','proper_square','dirac_weighted_inverse','variation','pv_variation')}
    gaps=[];wrong=[]
    relative=lambda a,b:float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
    for largest in (2.,10.,100.):
        for kappa in (.05,.7,3.):
            H,beta,N,V,dN,D,P,dP=benchmark(largest,kappa)
            sqrtN=np.diag(np.sqrt(np.diag(N)))
            for eta in (.0009,.04,2.89):
                pv_P=0.;pv_H=0.
                for j,cj in enumerate((1,-3,3,-1)):
                    m=np.sqrt(eta+j)
                    Hm=H+m*N@beta;dHm=V+m*dN@beta
                    resolvent=np.linalg.inv(Hm+1j*kappa*eye)
                    first=np.linalg.inv(D+m*eye)
                    mapped=sqrtN@resolvent@sqrtN@beta
                    errors['inverse'].append(relative(first,mapped))
                    squared=np.linalg.inv(P+m*m*eye)
                    signed=(first-np.linalg.inv(D-m*eye))/(2*m)
                    errors['proper_square'].append(relative(squared,signed))
                    averaged=(first+np.linalg.inv(D-m*eye))/2
                    errors['dirac_weighted_inverse'].append(relative(-D@squared,averaged))
                    proper=float(np.trace(squared@dP).real)
                    canonical=2*float(np.trace(resolvent@dHm).real)
                    lapse=2*float(np.trace(np.linalg.solve(N,dN)).real)
                    errors['variation'].append(abs(proper-canonical+lapse)/max(1.,abs(proper),abs(canonical)))
                    pv_P+=cj*proper;pv_H+=cj*canonical
                    gaps.append(float(np.min(abs(np.linalg.eigvalsh(Hm)))/(m*.2)))
                errors['pv_variation'].append(abs(pv_P-pv_H)/max(1.,abs(pv_P),abs(pv_H)))
                _,_,_,_,_,_,wrongP,_=benchmark(largest,np.sqrt(kappa*kappa+eta))
                wrong.append(dict(largest_lapse=largest,kappa=kappa,eta=eta,
                    relative_wrong_frequency_shift=relative(np.linalg.inv(P+eta*eye),np.linalg.inv(wrongP))))
    # An exact zero-energy block demonstrates a state-choice mismatch.
    beta0=np.diag([1.,-1.]);Cempty=np.zeros((2,2));Cpositive=(np.eye(2)-beta0)/2
    Cnegative=(np.eye(2)+beta0)/2;Caverage=(Cpositive+Cnegative)/2
    probe=np.array([[0.,1.],[1.,0.]])
    zero=dict(positive_mass_occupation=Cpositive.tolist(),negative_mass_occupation=Cnegative.tolist(),
        declared_empty_occupation=Cempty.tolist(),average_occupation=Caverage.tolist(),
        positive_mass_distance_from_empty=float(np.linalg.norm(Cpositive-Cempty,2)),
        averaged_distance_from_empty=float(np.linalg.norm(Caverage-Cempty,2)),
        metric_probe_zero_correction=float(np.trace((Caverage-Cempty)@probe)))
    maxima={key:max(value) for key,value in errors.items()}
    checks={key:bool(value<1e-7) for key,value in maxima.items()}
    checks.update(mass_gap_bound=min(gaps)>=1-1e-10,
        wrong_frequency_shift_detected=max(row['relative_wrong_frequency_shift'] for row in wrong)>.1,
        positive_mass_not_empty_zero=zero['positive_mass_distance_from_empty']==1.,
        paired_mass_not_empty_zero=zero['averaged_distance_from_empty']==.5,
        specific_zero_probe_cancels=zero['metric_probe_zero_correction']==0.)
    result=dict(passed=all(checks.values()),checks=checks,maximum_errors=maxima,
        minimum_gap_to_proved_bound_ratio=min(gaps),wrong_frequency_shift=wrong,zero_mode_control=zero,
        parameter_grid=dict(largest_lapse=[2.,10.,100.],kappa=[.05,.7,3.],eta=[.0009,.04,2.89],M=1.,matrix_dimension=6),
        scope='Finite exact-factorization algebra controls only. Not the sewn bridge, a stress value, a continuum remainder bound, or evidence for a time-travel medium.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs_sha256={'project_paths.py':hashlib.sha256(project_path('project_paths.py').read_bytes()).hexdigest()})
    project_path('dirac_green_mass_limit_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if not result['passed']: raise ArithmeticError('Green-identity controls failed; results retained.')
    plot(result)
    print(json.dumps({key:result[key] for key in ('passed','checks','maximum_errors','minimum_gap_to_proved_bound_ratio','zero_mode_control')},indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained')
    rows=[r for r in result['wrong_frequency_shift'] if r['kappa']==.7]
    for largest in (2.,10.,100.):
        selected=[r for r in rows if r['largest_lapse']==largest]
        axes[0].semilogx([r['eta'] for r in selected],[r['relative_wrong_frequency_shift'] for r in selected],
            'o-',label=f'max lapse {largest:g}')
    axes[0].set(xlabel='Proper squared-mass shift eta',ylabel='Normalized inverse difference',
        title='A proper mass is not a frequency shift')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[1].axis('off')
    steps=[('ESTABLISHED','Local Green identity at positive mass','#e0f0e9'),
        ('ESTABLISHED','Massless limit at fixed nonzero frequency','#e0f0e9'),
        ('OPEN','Joint infrared limit and stress coincidence','#fff0d3'),
        ('OPEN','Physical stress, supported bridge, stability','#fff0d3')]
    for i,(status,label,color) in enumerate(steps):
        y=.85-i*.22
        axes[1].text(.5,y,status+'\n'+label,ha='center',va='center',fontsize=9,
            transform=axes[1].transAxes,bbox=dict(boxstyle='round,pad=.55',fc=color,ec='#b4b4b4'))
        if i<3: axes[1].annotate('',xy=(.5,y-.14),xytext=(.5,y-.09),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
    fig.suptitle('Quantum-source calculation: one correspondence proved, physical support still open',fontsize=12)
    fig.supxlabel('Left: six-dimensional algebra control only. Right: analytic progress on the existing static model. No backward passage established.',fontsize=8)
    for ext in ('png','pdf'): fig.savefig(project_path('dirac_green_mass_limit_progress.'+ext),dpi=160)
    plt.close(fig)


if __name__=='__main__': run()
