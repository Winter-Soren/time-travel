"""One-sector zero-frequency return enclosure on the frozen finite DG model.

This is a bounded threshold experiment, not a frequency scan.  The actual
q=1 magnetic radial tail is enclosed to infinity at the retained outer
angular cutoff.  An open-DG accretivity identity supplies a bulk Lipschitz
bound before the physical mouth and facet tail maps are sewn.

The analytic inequalities below are evaluated with ordinary floating-point
arithmetic.  They are not interval-certified roundoff or a continuum PDE
certificate.  A reported failed Neumann gate is retained in the results.
Only mu=1/2 is tested; no nine-sector or infinite-angular conclusion follows.

Run: .venv\\Scripts\\python.exe -B -W error threshold_return_bound.py
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import block_diag, expm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from complete_angular_source import CompleteSector, selected_map
from exterior_boundary_map import ID, _complex_json
from magnetic_tail_map import _AngularCoupling
from throat_null_response import bump, mode_variation


INPUTS = ('complete_angular_source.py', 'complete_angular_modes.py',
    'dirac_sewing_audit.py', 'exterior_boundary_map.py',
    'exterior_tail_transport.py', 'magnetic_tail_map.py',
    'imaginary_weyl_leads.py', 'spherical_weyl_leads.py',
    'two_mouth_mesh.py', 'throat_null_response.py',
    'long_throat_integral_results.json')


def upper_norm(a):
    """Frobenius upper bound on the spectral norm, not an interval routine."""
    return float(np.linalg.norm(a, 'fro'))


def hermitian_upper(a):
    """Gershgorin upper bound on the largest Hermitian-part eigenvalue."""
    h = (a+a.conj().T)/2
    return float(np.max(h.diagonal().real+np.sum(abs(h), axis=1)-abs(h.diagonal())))


def far_enclosure(L=40., nmin=1, J=9, q=1., d=1., K=.37):
    a = q*d/(2*(1-d/(2*L))**3)
    ball = 2*a/(2*nmin+1)+(2*a*K+(J+a/L)*K*K)/((2*nmin+2)*L)
    contraction = (2*a+2*(J+a/L)*K)/((2*nmin+2)*L)
    c = nmin-(a+J*K)/L-a*K/L**2
    return dict(L=L,nmin=nmin,J=J,q=q,d=d,a=a,K=K,
        ball_image_radius=ball,contraction=contraction,
        Z_terminal_error=K/L,radial_decay_exponent_lower=c,
        energy_terminal_upper=2*L/(2*c-1),
        passed=bool(ball<K and contraction<1 and c>.5),
        formula='C=(D+E)/r, ||E||<=a/r; sup r||Z0+I||<=K; E0(L)<=2L/(2c-1) I',
        scope='Actual point-pair tail for r>=L at fixed finite angular projection J, not the complete angular operator.')


def coupling_with_derivative(angular, radius):
    r=float(radius);a=angular.d/2;c=angular.cosine;s=angular.sine
    rp=np.sqrt(r*r+2*a*r*c+a*a);rm=np.sqrt(r*r-2*a*r*c+a*a)
    potential=angular.q/2*((r*c+a)/rp-(r*c-a)/rm)
    derivative=-angular.q*angular.d*r*s*s/4*(rp**-3+rm**-3)
    J=np.einsum('p,pmn->mn',potential,angular.kernel)
    Jp=np.einsum('p,pmn->mn',derivative,angular.kernel)
    C=(angular.free+1j*J)/r
    Cp=-(angular.free+1j*J)/r**2+1j*Jp/r
    return C,Cp


def coupling_bounds(lo,hi,J=9,q=1.,d=1.):
    """Interval bounds on ||C||,||C'||,||C''|| from the actual pair field."""
    a=d/2;den=lo-a
    # Projection is contractive in exact angular arithmetic. The retained
    # harmonic Grams are GL-exact polynomials; this tiny allowance is stated.
    g=1+1e-12
    b0=g*q*d*hi*hi/(2*den**3)
    b1=g*q*d*hi/(2*den**3)
    b2=g*q*d/2*(den**-3+3*hi*(hi+a)/den**5)
    return (J/lo+b0/lo, J/lo**2+b0/lo**2+b1/lo,
            2*J/lo**3+2*b0/lo**3+2*b1/lo**2+b2/lo)


def scalar_riccati_step(value,linear,quadratic,forcing,length):
    """Exact constant-coefficient comparison flow, barring floating rounding."""
    matrix=np.array([[linear,forcing],[-quadratic,0.]])
    packed=expm(length*matrix)@np.array([value,1.])
    if packed[1]<=0:return math.inf
    return max(0.,float(packed[0]/packed[1]))


def propagate_enclosure(angular,solution,far,R=2.,log_step=.001):
    """A posteriori residual enclosure of a piecewise-linear reference curve.

    This does not treat solve_ivp's requested tolerance as an error bound.
    The actual polynomial curve is differentiated and its Riccati defect is
    enclosed with explicit C,C',C'' bounds on every radial subinterval.
    """
    L=far['L'];count=len(angular.channels)
    intervals=int(math.ceil(math.log(L/R)/log_step))
    nodes=np.geomspace(L,R,intervals+1)
    values=solution.sol(nodes).T.reshape(-1,count,count)
    error=far['Z_terminal_error'];energy=far['energy_terminal_upper']
    maximum_defect=maximum_error=0.;largest_log_norm=-math.inf
    floor=1e-13
    for j in range(intervals):
        hi,lo=nodes[j:j+2];length=hi-lo;mid=(hi+lo)/2
        zh,zl=values[j:j+2];z=(zh+zl)/2;zp=(zh-zl)/length
        C,Cp=coupling_with_derivative(angular,mid)
        c0,c1,c2=coupling_bounds(lo,hi,angular.nmax,angular.q,angular.d)
        # Spectral norms are upper-bounded via exact Hermitian row sums.
        zn=math.sqrt(max(0.,hermitian_upper(zh.conj().T@zh),hermitian_upper(zl.conj().T@zl)))
        zpn=upper_norm(zp)
        residual=zp-(C.conj().T-z@C@z)
        residual_derivative=-(Cp.conj().T-zp@C@z-z@Cp@z-z@C@zp)
        second=(1+zn*zn)*c2+4*zpn*zn*c1+2*zpn*zpn*c0
        defect=upper_norm(residual)+length/2*upper_norm(residual_derivative)+length**2/8*second
        defect+=floor*(1+c0+zpn)
        half_variation=length/2*(zpn*c0+zn*c1)
        linear=hermitian_upper(z@C)+hermitian_upper(C@z)+2*half_variation
        previous=error
        error=scalar_riccati_step(error,linear,c0,defect,length)
        if not math.isfinite(error):
            return dict(passed=False,reason='Riccati comparison envelope developed a pole',at_radius=lo)
        emax=max(previous,error)
        logarithmic=hermitian_upper(C@z)+half_variation+c0*emax
        multiplier=math.exp(2*logarithmic*length)
        contribution=math.expm1(2*logarithmic*length)/logarithmic if logarithmic else 2*length
        energy=multiplier*energy+contribution
        maximum_defect=max(maximum_defect,defect)
        maximum_error=max(maximum_error,error)
        largest_log_norm=max(largest_log_norm,logarithmic)
    Z=values[-1]
    return dict(passed=True,R=R,L=L,log_step=log_step,intervals=intervals,
        Z_error_upper=error,energy_upper=energy,
        maximum_local_defect_upper=maximum_defect,
        maximum_propagated_error=maximum_error,
        maximum_logarithmic_norm_upper=largest_log_norm,
        nominal_unitarity_error=upper_norm(Z.conj().T@Z-np.eye(count)),
        arithmetic_allowance_per_residual=floor,
        rounding_scope='Analytic defect and comparison formulas evaluated in float64; the allowance is not an interval-arithmetic proof.',
        Z_nominal=_complex_json(Z))


def zero_tail(R=2.,L=40.,mu=.5,J=9):
    angular=_AngularCoupling(mu=mu,nmax=J)
    count=len(angular.channels)
    def rhs(r,state):
        Z=state.reshape(count,count);C=angular.matrix(r)
        return (C.conj().T-Z@C@Z).ravel()
    solved=solve_ivp(rhs,(L,R),(-np.eye(count,dtype=complex)).ravel(),
        method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True,max_step=.1)
    if not solved.success:raise RuntimeError(solved.message)
    far=far_enclosure(L,int(angular.channels[0]),J)
    coarse=propagate_enclosure(angular,solved,far,R,.002)
    fine=propagate_enclosure(angular,solved,far,R,.001)
    Z=solved.y[:,-1].reshape(count,count)
    return Z,dict(far=far,coarse=coarse,fine=fine,nfev=solved.nfev,
        angular_gram_error=angular.gram_error,
        terminal_nominal='-I is only the center of the certified far-tail ball; it is not asserted to equal the magnetic zero-frequency map at L.',
        exact_identity='E0-Ek=integral (Fk-F0)^dagger(Fk-F0) dr+(Zk-Z0)^dagger(Zk-Z0)/k; ||Zk-Z0||<=k||E0||.',
        scope='Finite J actual pair projected on r>=R, including an analytic enclosure of r>L. No infinite-angular or full continuum DG certificate.')


def tail_green_control():
    """Independent nonnormal constant-C check, without an exterior solve."""
    rng=np.random.default_rng(104)
    def unitary():
        z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        return np.linalg.qr(z)[0]
    U,V=unitary(),unitary();m=np.array([.4,1.2,3.]);k=.23
    omega=np.hypot(m,k);beta=m/(omega+k)
    Z0=-V@U.conj().T;Z=-V@np.diag(beta)@U.conj().T
    E0=U@np.diag(1/m)@U.conj().T
    Ek=U@np.diag((1+beta**2)/(2*omega))@U.conj().T
    cross=U@np.diag((1+beta)/(m+omega))@U.conj().T
    delta=Z-Z0;difference_energy=E0+Ek-cross-cross.conj().T
    derivative=V@np.diag(1/m)@U.conj().T
    return dict(seed=104,kappa=k,
        cross_identity_error=upper_norm(np.eye(3)-Z0.conj().T@Z-k*cross),
        current_identity_error=upper_norm(np.eye(3)-Z.conj().T@Z-2*k*Ek),
        strong_energy_identity_error=upper_norm(E0-Ek-difference_energy-delta.conj().T@delta/k),
        derivative_order_error=upper_norm(derivative+Z0@E0),
        wrong_order_negative_control=upper_norm(derivative+E0@Z0),
        lipschitz_ratio=float(np.linalg.norm(delta,2)/(k*np.linalg.norm(E0,2))),
        free_radial_normalization='C=n/r gives E0(R)=2R/(2n-1), including both radial characteristics.')


def open_gram_check(sector):
    model=sector.model;rows=[];cols=[];values=[]
    def add(i,j,block):
        ii=6*i+np.arange(6);jj=6*j+np.arange(6)
        rows.extend(np.repeat(ii,6));cols.extend(np.tile(jj,6));values.extend(block.ravel())
    for edge in model.boundary:
        mass=np.einsum('p,pa,pb->ab',edge['weights'],edge['basis'],edge['basis'])
        add(edge['triangle'],edge['triangle'],np.kron(mass,ID))
    for edge in model.interior:
        a,b=edge['basis'],edge['other_basis'];w=edge['weights'];i,j=edge['triangle'],edge['other']
        aa=np.einsum('p,pa,pb->ab',w,a,a);bb=np.einsum('p,pa,pb->ab',w,b,b)
        ab=np.einsum('p,pa,pb->ab',w,a,b)
        add(i,i,np.kron(aa,ID));add(j,j,np.kron(bb,ID))
        add(i,j,-np.kron(ab,ID));add(j,i,-np.kron(ab.T,ID))
    gram=coo_matrix((values,(rows,cols)),shape=sector.unshifted.shape).tocsc()
    difference=sector.unshifted+sector.unshifted.conj().T-gram
    absolute=float(np.max(abs(difference.data))) if difference.nnz else 0.
    return dict(maximum_absolute_residual=absolute,
        relative_residual=absolute/float(np.max(abs(sector.unshifted.data))),
        formula='A0+A0_dagger = sum boundary rho|u|^2 + sum interior rho|jump u|^2, with the same face quadrature and active axis domain.',
        minimum_cell_mass_eigenvalue=sector.assembly['minimum_cell_mass_eigenvalue'])


def open_response(sector,left=False):
    active=sector.active
    zero=sector.left_zero if left else sector.unshifted
    mouthG,mouthC=(sector.GL,sector.CL) if left else (sector.Gm,sector.Cm)
    outer=sector.left_outer if left else sector.outer
    G=np.column_stack((mouthG,outer['G']))[active]
    C=np.vstack((mouthC,outer['C']))[:,active]
    A=zero[active][:,active].tocsc();M=sector.mass[active][:,active]
    lu=splu(A);X=lu.solve(G);Y=lu.solve(C.conj().T,trans='H')
    S=C@X;EU=X.conj().T@(M@X);EV=Y.conj().T@(M@Y)
    u2=max(0.,hermitian_upper(EU));v2=max(0.,hermitian_upper(EV))
    return S,dict(energy_U=_complex_json(EU),energy_V=_complex_json(EV),
        norm_U_upper=math.sqrt(u2),norm_V_upper=math.sqrt(v2),
        response_lipschitz_upper=math.sqrt(u2*v2),
        response_norm=float(np.linalg.norm(S,2)),
        primal_relative_residual=float(np.linalg.norm(A@X-G)/np.linalg.norm(G)),
        dual_relative_residual=float(np.linalg.norm(A.conj().T@Y-C.conj().T)/np.linalg.norm(C)),
        formula='S(k)-S(0)=-k V_dagger (I+k T^-1)^-1 U; T=M^-1/2 A0 M^-1/2 accretive.',
        solve_scope='Numerically computed U,V with residuals; no validated interval inverse or roundoff enclosure.')


def zero_core(nmax=4,rho_c=32.):
    ell=.2*math.hypot(1.,rho_c);cut=math.atan(rho_c);count=2*nmax+1
    B=np.zeros((count,count),complex);P=np.zeros((2,count,count),complex)
    B[0,nmax]=1.;records=[];slope=np.zeros(2)
    length=2*ell*cut
    for n in range(1,nmax+1):
        lam=math.sqrt(n*(n+1));action=2*lam*math.asinh(rho_c)
        t=1/math.cosh(action);r=-1j*math.tanh(action)
        j=n-1;upper=n;lower=nmax+n
        B[upper,j]=B[lower,lower]=r;B[upper,lower]=B[lower,j]=t
        local=[]
        for probe in range(2):
            Jf=quad(lambda s:lam/math.cos(s)*bump(s,cut,probe),-cut,cut,epsabs=2e-13,epsrel=2e-12)[0]
            dr=1j*Jf*t*t;dt=Jf*t*math.tanh(action)
            P[probe,upper,j]=P[probe,lower,lower]=dr
            P[probe,upper,lower]=P[probe,lower,j]=dt
            definition=sector_probe(probe,cut)
            vmax=lam/(ell*min(math.cos(definition[0]),math.cos(definition[1])))
            # Conservative continuous open-core inverse bound 2*length;
            # incoming/adjoint zero response energies <=2*length.
            derivative=8*length*length*vmax+4*length
            slope[probe]+=4*derivative  # pair trace norm, both actual Weyls
            local.append(dict(J_f=Jf,zero_derivative_norm=Jf*t,
                mass_perturbation_sup_upper=vmax,mixed_frequency_metric_norm_upper=derivative))
        records.append(dict(n=n,mass_action=action,transmission=t,probes=local))
    return B,P,dict(ell=ell,optical_length=length,records=records,
        frequency_lipschitz=2*length,pressure_trace_norm_slope=slope.tolist(),
        derivation='Open zero-frequency core Green kernel has pointwise norm<=2, hence inverse L2 norm<=2L; incoming and adjoint response energy<=2L. Accretivity gives ||P_n(k)-P_n(0)||<=k*(8L^2 sup(m f)+4L).',
        zero_scope='The chiral zero has P(0)=0 and exact ||P(k)||=2k I_f exp(-kL); it is added separately.')


def sector_probe(index,cut):
    center,width=((0.,.65),(.60,.25))[index]
    return cut*(center-width),cut*(center+width)


def convert(Fin,value,Fout):
    return np.linalg.solve(Fout.T,(Fin@value).T).T


def _hashes():
    return {p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in INPUTS}


def run_checks():
    started=time.perf_counter();inputs=_hashes()
    sector=CompleteSector(.5,h=.05,nmax=4,outer_nmax=9,rho_c=32.)
    N=4;c=2*N+1;p=len(sector.config['outer_channels']);ell=.2*math.hypot(1.,32.)
    kmax=1e-6/ell
    print('104: m=0 open DG assembled',flush=True)
    gram=open_gram_check(sector);SR,dr=open_response(sector);SL,dl=open_response(sector,True)
    print('104: both open zero-frequency responses solved',flush=True)
    Z,tail=zero_tail();enclosure=tail['fine']
    raw,P,core=zero_core();raw_norms=[float(np.linalg.norm(block_diag(v,v.conj().T),ord='nuc')) for v in P]
    outRF=sector.right_factors['outer'];outLF=sector.left_factors['outer']
    ZR=convert(outRF[1],Z,outRF[0]);ZL=convert(outLF[1],Z.conj().T,outLF[0])
    gamma_core=max(float(np.linalg.norm(sector.Fin,2)*np.linalg.norm(np.linalg.inv(sector.Fout),2)),
                   float(np.linalg.norm(sector.Lin,2)*np.linalg.norm(np.linalg.inv(sector.Lout),2)))
    gamma_outer=max(float(np.linalg.norm(outRF[1],2)*np.linalg.norm(np.linalg.inv(outRF[0]),2)),
                    float(np.linalg.norm(outLF[1],2)*np.linalg.norm(np.linalg.inv(outLF[0]),2)))
    openS=block_diag(SR,SL);lopen=max(dr['response_lipschitz_upper'],dl['response_lipschitz_upper'])
    eye=np.eye(len(openS));mouth_indices=np.r_[np.arange(c),c+p+np.arange(c)]
    outer_indices=np.r_[c+np.arange(p),2*c+p+np.arange(p)]
    phases=[]
    for angle in (0.,math.pi):
        b=selected_map(raw,.5,angle,N)
        br=convert(sector.Fin,b,sector.Fout);bl=convert(sector.Lin,b.conj().T,sector.Lout)
        B=block_diag(br,ZR,bl,ZL);M=eye-openS@B
        inverse=np.linalg.solve(M,eye);inverse_residual=upper_norm(eye-inverse@M)
        rinv=upper_norm(inverse);RS=inverse@openS
        tail_error=gamma_outer*(enclosure['Z_error_upper']+kmax*enclosure['energy_upper'])
        core_error=kmax*gamma_core*core['frequency_lipschitz']
        total_B_error=max(tail_error,core_error)
        Bsup=float(np.linalg.norm(B,2))+total_B_error
        neumann=(inverse_residual+rinv*kmax*lopen*Bsup
            +upper_norm(RS[:,mouth_indices])*core_error
            +upper_norm(RS[:,outer_indices])*tail_error)
        inverse_upper=rinv/(1-neumann) if neumann<1 else None
        pressures=[];Dfull=[]
        for v in P:
            vv=selected_map(v,.5,angle,N)
            d=block_diag(convert(sector.Fin,vv,sector.Fout),np.zeros((p,p)),
                convert(sector.Lin,vv.conj().T,sector.Lout),np.zeros((p,p)))
            pressures.append(float(np.trace(RS@d).real/math.pi));Dfull.append(d)
        ref_gap=math.sqrt(2)-kmax*(core['frequency_lipschitz']+2*ell*(math.pi/2-math.atan(32.)))
        source_bounds=[]
        for j,probe in enumerate(sector.probes):
            pnorm=raw_norms[j]+kmax*(core['pressure_trace_norm_slope'][j]+4*probe['I_f'])
            snorm=kmax*probe['I_f']*2*c
            factor=((inverse_upper*(1+kmax*lopen)*gamma_core+1/ref_gap)/math.pi
                    if inverse_upper is not None else None)
            source_bounds.append(dict(pressure_derivative_trace_norm_upper=pnorm,
                shift_derivative_trace_norm_upper=snorm,
                pressure_integral_absolute_upper=kmax*factor*pnorm if factor is not None else None,
                shift_integral_absolute_upper=kmax*factor*snorm if factor is not None else None,
                either_null_integral_absolute_upper=kmax*factor*(pnorm+2*snorm) if factor is not None else None))
        phases.append(dict(phase=angle,minimum_nominal_sewing_singular_value=float(np.linalg.svd(M,compute_uv=False)[-1]),
            nominal_inverse_frobenius_upper=rinv,inverse_residual=inverse_residual,
            core_map_error_upper=core_error,tail_map_error_upper=tail_error,
            bulk_neumann_term=rinv*kmax*lopen*Bsup,
            core_neumann_term=upper_norm(RS[:,mouth_indices])*core_error,
            tail_neumann_term=upper_norm(RS[:,outer_indices])*tail_error,
            neumann_epsilon=neumann,neumann_passed=bool(neumann<1),inverse_upper=inverse_upper,
            reference_inverse_gap_lower=ref_gap,
            nominal_connected_pressure_at_zero=pressures,
            conditional_source_remainder=source_bounds,
            B_nominal=_complex_json(B),pressure_derivatives=_complex_json(Dfull)))
    # A separate exact zero-frequency core comparison uses the frozen radial ODE.
    core_errors=[]
    for n in (0,1,4):
        S,dS,_=mode_variation(0.,n,ell,32.,rtol=3e-13,atol=3e-15)
        action=2*math.sqrt(n*(n+1))*math.asinh(32.)
        t=1/math.cosh(action);expected=np.array([[-1j*math.tanh(action),t],[t,-1j*math.tanh(action)]])
        dp=[]
        for j in range(2):
            Jf=0. if n==0 else core['records'][n-1]['probes'][j]['J_f']
            dp.append(Jf*np.array([[1j*t*t,t*math.tanh(action)],
                                   [t*math.tanh(action),1j*t*t]]))
        core_errors.append(dict(n=n,error=float(np.max(abs(S-expected))),
            derivative_error=float(np.max(abs(dS-np.array(dp))))))
    green=tail_green_control()
    tests=dict(far_volterra=tail['far']['passed'],inward_comparison=enclosure['passed'],
        open_dg_gram_identity=gram['relative_residual']<1e-12,
        positive_mass=gram['minimum_cell_mass_eigenvalue']>0,
        open_solves=max(dr['primal_relative_residual'],dr['dual_relative_residual'],dl['primal_relative_residual'],dl['dual_relative_residual'])<1e-10,
        opposite_open_adjoint=float(np.max(abs(SL-SR.conj().T)))<1e-10,
        nominal_tail_unitarity=enclosure['nominal_unitarity_error']<1e-9,
        exact_zero_core=max(max(v['error'],v['derivative_error']) for v in core_errors)<1e-10,
        nonnormal_tail_green_identity=max(green[k] for k in ('cross_identity_error','current_identity_error',
            'strong_energy_identity_error','derivative_order_error'))<1e-12,
        wrong_order_control=green['wrong_order_negative_control']>1e-2,
        inputs_preserved=_hashes()==inputs)
    return dict(passed=all(tests.values()),tests=tests,
        certificate_gate_passed=all(v['neumann_passed'] for v in phases),
        parameters=dict(mu=.5,m=0,q=1,b=.2,d=1,rho_c=32.,ell=ell,h=.05,R=2.,L=40.,nmax=4,outer_nmax=9,
            kappa_interval=[0.,kmax],scaled_x_interval=[0.,1e-6]),
        open_gram=gram,open_right=dr,open_left=dl,tail=tail,core=core,
        tail_green_control=green,probes=sector.probes,
        core_controls=core_errors,gram_conversion=dict(core_norm_factor=gamma_core,outer_norm_factor=gamma_outer),
        phases=phases,matrices=dict(S_open_right=_complex_json(SR),S_open_left=_complex_json(SL),
            Z_zero_nominal=_complex_json(Z),F_mouth_right_in=_complex_json(sector.Fin),
            F_mouth_right_out=_complex_json(sector.Fout),F_mouth_left_in=_complex_json(sector.Lin),
            F_mouth_left_out=_complex_json(sector.Lout),F_outer_right=[_complex_json(v) for v in outRF],
            F_outer_left=[_complex_json(v) for v in outLF]),
        scope=[
            'Only the actual retained m=0, N4, J9, h=.05 model is assessed; no other m sector is inferred.',
            'The far tail is enclosed to infinity in the finite J projected equation; this does not retrofit the L40 truncation of the saved nonzero-frequency integral.',
            'Conditional analytic finite-DG bound: open inverse responses and matrix norms use floating arithmetic with recorded residuals, not verified interval roundoff.',
            'The source trace bounds include the AdS reference, both actual Weyls and the two named compact probes. They are not pointwise stress bounds.',
            'All Neumann outcomes are saved. A failed gate is not repaired by sign tests or endpoint sampling.',
            'No continuum/infinite-angular/UV completion, smooth junction, Einstein solution, absolute support or chronology claim.'],
        input_sha256=inputs,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.perf_counter()-started)


if __name__=='__main__':
    result=run_checks()
    path=Path(__file__).with_name('threshold_return_bound_results.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert json.loads(path.read_text(encoding='utf-8'))==result
    print(json.dumps(dict(passed=result['passed'],certificate_gate_passed=result['certificate_gate_passed'],
        seconds=result['seconds'],phases=[dict(phase=p['phase'],epsilon=p['neumann_epsilon'],
        gap=p['minimum_nominal_sewing_singular_value'],bounds=p['conditional_source_remainder']) for p in result['phases']]),indent=2))
    if not result['passed']:raise AssertionError('Threshold analytic/numerical control failed; diagnostics saved')
