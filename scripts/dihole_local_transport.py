"""Explicit normal balls, local heat-transport jets and graph evaluation.

All constants use D=1 and the actual static collar. No heat-defect norm,
renormalized source, stable opening or atom transport is asserted.
"""
import json
import math

from mpmath import mp

from dihole_channel_bounds import MPIntervalContext, calibration, field_function, lower, upper
from dihole_slice_channel import sha
from project_paths import project_path


NAME='dihole_local_transport'


def dprod(a,b):
    """Ordinary derivative bounds, with the full binomial product rule."""
    return [sum(math.comb(n,j)*a[j]*b[n-j] for j in range(n+1))
            for n in range(min(len(a),len(b)))]


def pmul(a,b,n):
    return [sum(a[j]*b[k-j] for j in range(max(0,k-len(b)+1),min(k+1,len(a)))) for k in range(n+1)]


def compose(a,z,n):
    """Normalized Frechet jets; the increment must have zero constant term."""
    if z[0] != 0:raise ValueError('Composition requires a zero-constant increment.')
    result=[0]*(n+1);power=[1]+[0]*n
    for coefficient in a[:n+1]:
        result=[x+coefficient*y for x,y in zip(result,power)]
        power=pmul(power,z,n)
    return result


def packed(iv,x):
    """Exact upward dyadic endpoint, including values beyond binary64 range."""
    x=iv.mpf(x)
    sign,man,exponent,bits=x._mpi_[1]
    if sign or bits<0:raise ArithmeticError('Finite nonnegative ceiling required.')
    return dict(mantissa=str(man),binary_exponent=exponent,
                log10_upper=upper(iv.log(x,10)) if man else None)


def real_metric_bounds(iv,g,radial=24,angular=32):
    fun,par=field_function(iv),calibration(iv)
    mag=lambda z:max(abs(lower(z)),abs(upper(z)))
    b=dict.fromkeys(('T_R','V_R','U_R','T_x','V_x','U_x','F_x','angular_R'),0.)
    rlo,rhi=iv.mpf(1)/40,iv.mpf(7)/100
    for i in range(radial):
        R=iv.mpf([rlo+(rhi-rlo)*i/radial,rlo+(rhi-rlo)*(i+1)/radial])
        for j in range(angular):
            x=iv.mpf([-1+iv.mpf(2)*j/angular,-1+iv.mpf(2)*(j+1)/angular])
            z=fun(R,x,*par);f,W=z[:2];fR,WR=z[3:5];fx,Wx=z[6:8]
            if f.a<=0 or W.a<=0:raise ArithmeticError('Positive real factors required.')
            U,V=W/f,f*W
            TR=fR/W-f*WR/W**2;VR=fR*W+f*WR;UR=WR/f-W*fR/f**2
            Tx=fx/W-f*Wx/W**2;Vx=fx*W+f*Wx;Ux=Wx/f-W*fx/f**2
            Fx=Wx*(f-1/f)+W*fx*(1+1/f**2)
            for key,val in dict(T_R=TR,V_R=VR,U_R=UR,T_x=Tx,V_x=Vx,U_x=Ux,F_x=Fx).items():
                b[key]=max(b[key],mag(val))
            b['angular_R']=max(b['angular_R'],mag(4*(2*R*U+R**2*UR)),mag(4*(2*R*V+R**2*VR)))
    nl,nh,pl,ph=(iv.mpf(g[k]) for k in ('N_lower','N_upper','P_lower','P_upper'))
    C2=3*iv.sqrt(3)/16;C3=25*iv.sqrt(5)/216;C4=343*iv.sqrt(7)/4096
    M0=max(ph**2,1/nl**2);Mx=iv.mpf(max(b['U_x'],b['V_x']))
    angular_y=16*rhi**2*C3*M0+16*rhi**2*C4*Mx+8*rhi**2*C3*iv.mpf(b['F_x'])
    scalar_y=4*C2*iv.mpf(max(b['T_x'],b['V_x']))
    Ms=max(iv.mpf(b[k]) for k in ('T_R','V_R','angular_R'));My=max(angular_y,scalar_y)
    M=max(nh**2,ph**2,4*rhi**2*max(ph**2,1/nl**2))
    H=max(1/nl**2,1/pl**2,9/(4*rlo**2)*max(1/pl**2,nh**2))
    return dict(boxes=radial*angular,jets=b,metric_norm_upper=upper(M),inverse_metric_norm_upper=upper(H),
        metric_lower=lower(1/H),first_component_upper=upper(max(Ms,My)),radial_first_upper=upper(Ms),
        angular_first_upper=upper(My),directional_first_upper=upper(iv.sqrt(Ms**2+2*My**2)),
        volume_lower=lower(4*pl**2*rlo**2/9),volume_upper=upper(4*ph**2*rhi**2))


def normal_geometry(iv,prior,real):
    fields=prior['jets']['total_order_derivative_upper']
    g,h=(list(map(iv.mpf,fields[k])) for k in ('metric','inverse_metric'))
    g[0]=iv.mpf(real['metric_norm_upper']);g[1]=iv.mpf(real['first_component_upper'])
    h[0]=iv.mpf(real['inverse_metric_norm_upper'])
    # Differentiate g^{-1}g=I. The factor4 turns component bounds into
    # matrix operator bounds; using the minimum retains both certificates.
    inverse_op=[h[0]]
    for n in range(1,20):
        value=h[0]*sum(math.comb(n,j)*4*g[j]*inverse_op[n-j] for j in range(1,n+1))
        inverse_op.append(value)
        h[n]=iv.mpf(min(upper(h[n]),upper(value)))
    gamma=[6*x for x in dprod(h,g[1:])]
    C=[2**(n+3)*gamma[n] for n in range(18)]
    A=iv.log(iv.mpf(3)/2)
    root=A/(2*C[0]+iv.sqrt(4*C[0]**2+2*C[1]*A))
    radius=iv.mpf(min(.005,lower(root)))
    m=iv.mpf(real['metric_lower']);M=iv.mpf(real['metric_norm_upper'])
    rho=iv.mpf(lower(iv.sqrt(m)*radius))
    exponent=4*C[0]*radius+2*C[1]*radius**2
    if exponent.b>A.a or radius.b>iv.mpf(1/200).a:raise ArithmeticError('Normal-ball inequality failed.')
    H=[]
    for k in range(18):
        H.append(radius*sum(math.comb(k,j)*math.factorial(2)/math.factorial(2-j)*2**(2-j)
                           *radius**(k-j)*C[k-j] for j in range(min(2,k)+1)))
    if H[1].b>=1:raise ArithmeticError('Scaled phase Lipschitz constant exceeds1.')
    Z=[iv.mpf(0),iv.mpf(3)/2]
    for n in range(2,18):
        forcing=compose([0,0]+[H[k]/math.factorial(k) for k in range(2,n+1)],Z+[0],n)[n]
        Z.append(iv.exp(1)*forcing)
    E=[iv.mpf(0)]+[radius**(1-n)*Z[n] for n in range(1,18)]
    Y=[iv.mpf(0),iv.mpf(2)]
    for n in range(2,18):Y.append(2*compose([0,0]+E[2:n+1],Y+[0],n)[n])
    return dict(radius=radius,rho=rho,gamma=gamma,C=C,Z=Z,exp=E,inverse=Y,
        exported=dict(coordinate_velocity_radius_lower=lower(radius),normal_metric_radius_lower=lower(rho),
            closed_cutoff_radius_lower=lower(rho/2),coordinate_clearance=.02,
            variational_exponent_upper=upper(exponent),phase_lipschitz_upper=upper(H[1]),
            jacobian_singular_values=[.5,1.5],inverse_jacobian_norm_upper=2.,
            christoffel_component_upper=[packed(iv,x) for x in gamma[:18]],
            exp_normalized_frechet_upper=[packed(iv,x) for x in E[1:]],
            inverse_normalized_frechet_upper=[packed(iv,x) for x in Y[1:]],
            scaled_phase_normalized_frechet_upper=[packed(iv,x) for x in Z[1:]],
            orthonormal_volume_density_lower=lower(iv.mpf(real['volume_lower'])/(16*iv.mpf(real['volume_upper']))),
            orthonormal_volume_density_upper=upper(81*iv.mpf(real['volume_upper'])/(16*iv.mpf(real['volume_lower'])))))


def dirac_sequences(iv,prior,q):
    fields=prior['jets']['total_order_derivative_upper']
    E,I,A,V,Vi=(list(map(iv.mpf,fields[k])) for k in
        ('coframe','inverse_coframe','gauge_per_abs_q','half_density','inverse_half_density'))
    J=[abs(q)*x+y for x,y in zip(A,dprod(V[1:],Vi))]
    b=[768*x+16*y for x,y in zip(dprod(dprod(I,I),E[1:]),dprod(I,J))]
    return [4*x for x in I],b,(E,I,A,V,Vi)


def point_bound(iv,prior,record,real):
    a,b,_=dirac_sequences(iv,prior,record['q'])
    sigma=iv.mpf(1)/10;r=iv.mpf(2)**-36;h=r/8
    if iv.mpf(real['metric_norm_upper']).b>=100 or (16*a[1]*r).b>sigma.a:
        raise ArithmeticError('Declared frozen-principal-symbol cube is not certified.')
    L=2/sigma;A=L*(b[0]+40*a[0]/h)+20/h
    U=[iv.mpf(1)]
    for ell in range(1,7):
        U.append((L+A)*U[ell-1]+L*sum(math.comb(ell-1,j)*
            (2**(j+1)*a[j]*U[ell-j]+2**j*b[j]*U[ell-1-j]) for j in range(1,ell)))
    W=[sum(math.comb(j,k)*2**k*math.factorial(k)*(20/h)**k*U[j-k] for k in range(j+1)) for j in range(7)]
    sobolev=W[0]+3*W[2]+3*W[4]+W[6]
    coeff=record['coefficient_operator_norm_upper']
    C=sobolev*(iv.mpf(coeff['deltaP_second'][0])/12+iv.mpf(coeff['deltaP_first'][0])/30+iv.mpf(coeff['deltaP_zero'][0])/48)
    return dict(C=C,exported=dict(q=record['q'],cube_half_width=lower(r),sigma=lower(sigma),
        freeze_ratio_upper=upper(16*a[1]*r/sigma),local_derivative_upper=[packed(iv,x) for x in U],
        H6_graph_factor_upper=packed(iv,sobolev),canonical_point_evaluation_upper=packed(iv,C),
        convention='||dotP_c u(x)|| <= Cev ||(1+P_c)^3 u||_L2(dx); fixed-space canonical operator. Physical-frame evaluation additionally has factor1/v(x).'))


def transport_bounds(iv,prior,record,normal,real):
    """Physical-volume coefficients u0..u4, not the heat-error defect."""
    q=record['q'];_,_,(frame,invframe,gauge,V,Vi)=dirac_sequences(iv,prior,q)
    R=normal['radius'];Z=normal['Z'];E=normal['exp'];Y=normal['inverse'];n=16
    # Physical spin+U(1) connection; unitary along every real geodesic.
    spin=dprod(dprod(dprod(frame,invframe),invframe),frame[1:])
    omega=[768*x+abs(q)*y for x,y in zip(spin,gauge)]
    omega_taylor=[2**(k+1)*omega[k]/math.factorial(k) for k in range(n+1)]
    path=[0]+[R*Z[k] for k in range(1,n+1)]
    velocity=[2*R]+path[1:]
    K=pmul(compose(omega_taylor,path,n),velocity,n)
    U=[iv.mpf(1)]
    for k in range(1,n+1):U.append(sum(K[j]*U[k-j] for j in range(1,k+1)))
    # J is differentiated in scaled initial velocity; retain its true
    # determinant lower bound, not the crude constant of 24*Jseries^4.
    J=[iv.mpf(3)/2]+[(k+1)*Z[k+1] for k in range(1,n+1)]
    determinant=[1]+[0]*n
    for _ in range(4):determinant=pmul(determinant,J,n)
    delta_det=[0]+[24*x*16 for x in determinant[1:]]
    inverse_root=compose([4*iv.mpf(math.comb(2*k,k))/4**k for k in range(n+1)],delta_det,n)
    vi_taylor=[2**k*Vi[k]/math.factorial(k) for k in range(n+1)]
    vicomposed=compose(vi_taylor,path,n)
    u0w=pmul(pmul(vicomposed,inverse_root,n),U,n)
    vmax=iv.sqrt(iv.mpf(real['volume_upper']))
    u0w=[vmax*x for x in u0w]
    u0=compose(u0w,[0]+[Y[k]/R for k in range(1,n+1)],n)
    # An independent zero-order inverse follows from unitary transport.
    inverse0=iv.sqrt(iv.mpf(real['volume_upper'])/iv.mpf(real['volume_lower']))*iv.mpf(9)/4
    u0inv=[inverse0]
    for k in range(1,n+1):u0inv.append(inverse0*sum(u0[j]*u0inv[k-j] for j in range(1,k+1)))
    # Convert the canonical P coefficients to physical volume before
    # applying the heat-transport recursion; ordered coordinate pairs.
    saved=record['coefficient_operator_norm_upper']
    p2,p1,p0=(list(map(iv.mpf,saved[key])) for key in ('P_second','P_first','P_zero'))
    v1=dprod(V[1:],Vi);v2=dprod(V[2:],Vi)
    p1g=[a+8*b for a,b in zip(p1,dprod(p2,v1))]
    p0g=[a+4*b+16*c for a,b,c in zip(p0,dprod(p1,v1),dprod(p2,v2))]
    p2,p1,p0=([2**k*x/math.factorial(k) for k,x in enumerate(row)] for row in (p2,p1g,p0g))
    X=compose(E,Y,n);X[0]=0 # exp_y(tau exp_y^-1(x)), with0<=tau<=1.
    coefficients=[u0]
    for j in range(1,5):
        order=n-2*j;previous=coefficients[-1]
        derivative=[(k+1)*previous[k+1] for k in range(order+1)]
        second=[(k+1)*(k+2)*previous[k+2] for k in range(order+1)]
        applied=[16*a+4*b+c for a,b,c in zip(pmul(p2,second,order),pmul(p1,derivative,order),pmul(p0,previous,order))]
        integrand=compose(pmul(u0inv,applied,order),X,order)
        coefficients.append([x/j for x in pmul(u0,integrand,order)])
    return dict(coefficients=coefficients,exported=dict(q=q,
        convention='Physical-volume kernel (4*pi*t)^-2 exp(-d(x,y)^2/(4t)) sum_j t^j uj(x,y); coefficients below bound ||D_x^n uj||/n! in coordinate Euclidean Frechet norm, source y fixed.',
        u_normalized_frechet_upper=[[packed(iv,x) for x in row] for row in coefficients],
        derivative_orders=[len(row)-1 for row in coefficients],
        u0_inverse_normalized_frechet_upper=[packed(iv,x) for x in u0inv],
        parallel_transport_scaled_normalized_frechet_upper=[packed(iv,x) for x in U],
        unitary_transport_used=True,heat_defect_B4=None))


def algebra_controls(iv):
    mp.dps=90
    product=dprod([1,2,6],[3,4,10])
    polynomial=compose([1,2,3],[0,4,5],4)
    expected=[mp.diff(lambda z:1+2*(4*z+5*z*z)+3*(4*z+5*z*z)**2,0,k)/math.factorial(k) for k in range(5)]
    # Inverse F(v)=v+c*v^2 has exact signed Catalan coefficients;
    # our absolute-series recursion uses the same positive magnitudes.
    c=mp.mpf(1)/10;y=[mp.mpf(0),mp.mpf(1)]
    for k in range(2,12):y.append(compose([0,0,c],y+[0],k)[k])
    inverse_error=max(abs(y[k]-mp.mpf(math.comb(2*k-2,k-1))/k*c**(k-1)) for k in range(1,12))
    # Constant unitary connection exp(i*q*v): normalized derivatives are
    # |q|^n/n!, below the unitarity-based recurrence |q|^n.
    unitary=max(float(abs(mp.diff(lambda v:mp.exp(7j*v),mp.mpf('.2'),k))/math.factorial(k)/7**k) for k in range(1,17))
    # An independent exact geodesic: g=exp(2*x) dx^2+dy^2+dz^2+dt^2,
    # exp_0(v)_x=log(1+v_x), inverse exp_x=expm1(x).
    gm=iv.exp(iv.mpf(1)/25)
    metric=[upper(2**k*gm) for k in range(21)]
    toy=dict(jets=dict(total_order_derivative_upper=dict(metric=metric,inverse_metric=metric)))
    real=dict(metric_norm_upper=upper(gm),inverse_metric_norm_upper=upper(gm),
        metric_lower=lower(1/gm),first_component_upper=upper(2*gm),
        volume_lower=lower(iv.exp(-iv.mpf(1)/50)),volume_upper=upper(iv.exp(iv.mpf(1)/50)))
    geometry=normal_geometry(iv,toy,real);r=mp.mpf(geometry['exported']['coordinate_velocity_radius_lower'])
    flow_ratios=[];inverse_ratios=[]
    for k in range(1,18):
        flow_ratios.append(float(1/(k*(1-r)**k)/mp.mpf(geometry['exp'][k]._mpi_[0])))
        inverse_ratios.append(float((1+r)/math.factorial(k)/mp.mpf(geometry['inverse'][k]._mpi_[0])))
    fourier=[]
    for k,ceiling in ((0,mp.mpf(1)/48),(1,mp.mpf(1)/60),(2,mp.mpf(1)/48)):
        value=mp.sqrt(mp.quad(lambda t:t**(3+2*k)/(1+t*t)**6,[0,1,mp.inf])/(8*mp.pi**2))
        fourier.append(float(value/ceiling))
    return dict(product=product==[3,10,44],composition_error=float(max(abs(x-y) for x,y in zip(polynomial,expected))),
                inverse_series_error=float(inverse_error),unitary_transport_ratio=unitary,
                exact_geodesic_max_ratio=max(flow_ratios),exact_inverse_geodesic_max_ratio=max(inverse_ratios),
                fourier_evaluation_max_ratio=max(fourier))


def run():
    iv=MPIntervalContext();iv.dps=80
    path=project_path('dihole_regular_jets_results.json');prior=json.loads(path.read_text())
    inputs={**prior['input_sha256'],path.name:sha(path),'dihole_regular_jets.py':prior['source_sha256']}
    if not prior['passed'] or not all(sha(project_path(n))==h for n,h in inputs.items()):raise ArithmeticError('Exact passing inputs required.')
    geometry=json.loads(project_path('dihole_buffer_smoothing_results.json').read_text())['geometry']
    real=real_metric_bounds(iv,geometry);normal=normal_geometry(iv,prior,real)
    points=[point_bound(iv,prior,row,real) for row in prior['operator_bounds']]
    transport=[transport_bounds(iv,prior,row,normal,real) for row in prior['operator_bounds']]
    controls=algebra_controls(iv)
    checks=dict(exact_inputs=all(sha(project_path(n))==h for n,h in inputs.items()),
        full_real_cover=real['boxes']==768,real_ellipticity=real['metric_lower']>0,
        improved_first_derivative=real['first_component_upper']<prior['jets']['total_order_derivative_upper']['metric'][1],
        normal_ball_positive=normal['exported']['normal_metric_radius_lower']>0,
        normal_jacobian=normal['exported']['variational_exponent_upper']<math.log(1.5)+1e-15,
        scaled_flow_lipschitz=normal['exported']['phase_lipschitz_upper']<1,
        frozen_principal_symbol=all(p['exported']['freeze_ratio_upper']<1 for p in points),
        product_control=controls['product'],composition_control=controls['composition_error']<1e-80,
        inverse_series_control=controls['inverse_series_error']<1e-80,
        unitary_control=controls['unitary_transport_ratio']<=1,
        exact_geodesic_control=controls['exact_geodesic_max_ratio']<1 and controls['exact_inverse_geodesic_max_ratio']<1,
        fourier_evaluation_control=controls['fourier_evaluation_max_ratio']<1,
        heat_transport_schedule=all(t['exported']['derivative_orders']==[16,14,12,10,8] for t in transport),
        physical_source_not_claimed=all(t['exported']['heat_defect_B4'] is None for t in transport))
    result=dict(passed=all(checks.values()),checks=checks,controls=controls,real_geometry=real,
        normal_geometry=normal['exported'],point_evaluation=[p['exported'] for p in points],
        heat_transport=[t['exported'] for t in transport],source_sha256=sha(__file__),input_sha256=inputs,
        bound_encoding='Each packed upper endpoint is exactly int(mantissa)*2**binary_exponent. log10_upper is an upward rounded display bound, not the stored constant.',
        missing=dict(heat_defect_B4=None,curved_heat_remainder=None,finite_target_source=None,stability=None,atom_transit=None),
        scope=['Actual static collar, full regular angular charts, exact atomic flux; no reduced angular sector or small-flux replacement.',
               'Local normal radius, graph-to-point bound and heat-transport coefficient derivative bounds are evaluated. No full heat-defect norm or renormalized supporting stress is evaluated.',
               'Physical heat transport and canonical coordinate graph estimates are kept distinct using the physical half-density.',
               'The mathematical radii and cutoffs are estimates, not opening sizes or physical manufacturing specifications.'],
        sources=['https://math.uchicago.edu/~dannyc/courses/riem_geo_2013/riem_geo_notes.pdf','https://arxiv.org/abs/1607.05152','https://dlmf.nist.gov/5.12.E1'])
    project_path(NAME+'_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not result['passed']:raise ArithmeticError('Local transport checks failed; evidence retained.')
    plot(result)
    print(json.dumps(dict(passed=True,checks=checks,real=real,normal_radius=normal['exported']['normal_metric_radius_lower'],
        point_log10=[p['exported']['canonical_point_evaluation_upper']['log10_upper'] for p in points],
        u4_log10=[[x['log10_upper'] for x in t['exported']['u_normalized_frechet_upper'][4]] for t in transport]),indent=2))


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.5,4.7),layout='constrained')
    for j,row in enumerate(result['heat_transport'][-1]['u_normalized_frechet_upper']):
        ax.plot(range(len(row)),[v['log10_upper'] for v in row],marker='.',label=f'u{j}')
    ax.set(xlabel='Endpoint derivative order',ylabel='log10 of normalized derivative ceiling, D = 1',
           title='Local heat-transport coefficients now have explicit bounds')
    ax.legend(ncol=5);ax.grid(alpha=.2)
    fig.supxlabel('Actual atomic-flux reference. Very conservative bounds, not actual coefficient sizes.\nThe heat defect, supporting source, stability and atom passage remain unproved.',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(project_path(NAME+'_progress.'+ext),dpi=170)
    plt.close(fig)


if __name__=='__main__':run()
