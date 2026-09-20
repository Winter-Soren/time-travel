"""A posteriori compact-probe bounds in the EXACT balanced dihole exterior.

The finite trial is the saved117 primal; the dual is constructed exactly by
reflection: b(u,v)=-sigma_z*a(v,u). No new PDE solve is performed. Quadrature
nodes, basis functions, contractions, geometry and sums are enclosed before
adding the analytic118 quadrature remainder. This is one compact resolvent
probe at positive imaginary frequency, not a quantum stress or a throat.
"""
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import mpmath
from mpmath.ctx_iv import MPIntervalContext

from dihole_interval_arithmetic import point, add, sub, mul, div, square, sqrt, sum_interval, matmul
from dihole_interval_arithmetic import run_checks as arithmetic_checks
from dihole_interval_quadrature import interval_grid, self_check as quadrature_checks
from dihole_quadrature_remainder import dyadic, outward, run_checks as remainder_checks
from project_paths import project_path, source_matches


def trans(a):
    return a[0].T, a[1].T


def outer(a, b):
    return mul((a[0][:, None], a[1][:, None]), (b[0][None, :], b[1][None, :]))


def export(a):
    return [float(a[0]), float(a[1])]


def exact_target(calibration, previous):
    """Only the coordinate map uses dyadic constants; gravity is balanced exactly."""
    iv = MPIntervalContext(); iv.dps = 75
    k, Q = iv.mpf(1)/2, iv.mpf(1)/5
    M = (iv.sqrt(k*k+4*Q*k)-k)/2
    ag = iv.sqrt(k*k-M*M)
    B = M/((k+ag)*(M+k))
    saved = previous['refined_corrected']['wide_q512000']
    L = dyadic(iv, saved['length_scale'])
    vstar = dyadic(iv, calibration['coefficients']['v_star'])
    rho_star = dyadic(iv, calibration['coefficients']['rho_star'])
    E = vstar/L
    assert iv.mpf('.153') < M < iv.mpf('.154') and 0 < ag < iv.mpf('.476')
    assert 0 < B < iv.mpf('.241') and L < iv.mpf('.011762')
    assert E > iv.mpf('43.3') and vstar > iv.mpf('.509')
    extent = L*iv.mpf('18.25')/iv.sqrt(2)
    assert rho_star-extent > iv.mpf('.515') and rho_star+extent < iv.mpf('.820')
    assert extent < iv.mpf('.152') and L*iv.mpf('.75')/iv.sqrt(2) < iv.mpf('.0063')
    support_rho = rho_star-9*iv.sqrt(2)*L
    assert support_rho > 0
    raw = dict(k=k, Q=Q, M=M, ag=ag, B=B, L=L, E=E, vstar=vstar,
               rho_star=rho_star, sqrt2=iv.sqrt(2))
    bounds = {name: outward(value) for name, value in raw.items()}
    target = {name: (np.asarray(value['lower']), np.asarray(value['upper'])) for name,value in bounds.items()}
    # The enlarged118 modulus3e25 also bounds these two additional integrands.
    cut_field = iv.mpf(18)*520000
    residual = 14+8*(2*(18*16+25)*iv.mpf(520000))+iv.mpf('2.2')*cut_field+500001*cut_field
    center_bound, source_bound = cut_field*(residual+14), iv.mpf(14)**2
    assert center_bound < iv.mpf('3e25') and source_bound < iv.mpf('3e25')
    return target, dict(exact_parameters='D=1,Q=1/5,k=1/2,M=(sqrt(k*k+4*Q*k)-k)/2,ag=sqrt(k*k-M*M),B=M/((k+ag)*(M+k))',
        coordinate_constants={name: float(value).hex() for name,value in
            (('L', saved['length_scale']), ('vstar', calibration['coefficients']['v_star']),
             ('rho_star', calibration['coefficients']['rho_star']))},
        exact_energy='E=vstar/L as an exact ratio, not the rounded117 energy.',
        parameter_enclosures=bounds, minimum_support_rho=outward(support_rho),
        center_integrand_modulus_bound=outward(center_bound),
        source_norm_integrand_modulus_bound=outward(source_bound))


def cutoff_grid(grid, inner, outer_radius):
    """Panel-aware real polynomial; all changes occur at exact panel edges."""
    x, ids, edges = grid['x'], grid['panel_ids'], grid['panel_edges']
    n = len(x[0])
    # These arrays are filled in place; point() deliberately aliases its ends.
    values = np.zeros(n), np.zeros(n)
    derivatives = np.zeros(n), np.zeros(n)
    for panel in range(len(edges)-1):
        low, high = edges[panel:panel+2]
        mask = ids == panel
        if -inner <= low and high <= inner:
            values[0][mask] = values[1][mask] = 1.
        elif high <= -outer_radius or low >= outer_radius:
            continue
        else:
            assert (inner <= low < high <= outer_radius or -outer_radius <= low < high <= -inner)
            sign = 1 if low >= inner else -1
            xx = x[0][mask], x[1][mask]
            t = div(sub(mul(point(sign), xx), point(inner)), point(outer_radius-inner))
            one_minus = sub(point(1), t)
            value = mul(mul(square(one_minus), one_minus),
                        add(add(point(1), mul(point(3),t)), mul(point(6),square(t))))
            derivative = mul(point(-30*sign/(outer_radius-inner)), mul(square(t), square(one_minus)))
            values[0][mask], values[1][mask] = value
            derivatives[0][mask], derivatives[1][mask] = derivative
    return values, derivatives


def geometry(grid, target):
    """Exact-balanced Weyl coefficients, including interval cancellation in mass."""
    t = target
    u = grid['x'][0][:,None], grid['x'][1][:,None]
    v = grid['x'][0][None,:], grid['x'][1][None,:]
    rho = add(t['rho_star'], div(mul(t['L'],add(u,v)),t['sqrt2']))
    z = div(mul(t['L'],sub(u,v)),t['sqrt2'])
    zp, zm = add(z,t['k']), sub(z,t['k'])
    dp, dm = sqrt(add(square(rho),square(zp))), sqrt(add(square(rho),square(zm)))
    total = add(dp,dm)
    r, c = add(t['M'],div(total,point(2))), div(mul(point(2),z),total)
    ag2 = square(t['ag'])
    sigma = sub(square(r),mul(ag2,square(c)))
    h = mul(dp,dm)
    speed = div(square(h),square(sigma))
    drho = div(add(div(rho,dp),div(rho,dm)),point(2))
    dz = div(add(div(zp,dp),div(zm,dm)),point(2))
    crho = mul(point(-2),div(mul(c,drho),total))
    cz = div(sub(point(2),mul(point(2),mul(c,dz))),total)
    derivatives = []
    for rj, cj, numerator_plus, numerator_minus in ((drho,crho,rho,rho),(dz,cz,zp,zm)):
        logh = add(div(numerator_plus,square(dp)),div(numerator_minus,square(dm)))
        sigma_j = mul(point(2),sub(mul(r,rj),mul(ag2,mul(c,cj))))
        derivatives.append(mul(point(2),mul(speed,sub(logh,div(sigma_j,sigma)))))
    vrho, vz = derivatives
    delta = sub(sub(square(r),mul(point(2),mul(t['M'],r))),ag2)
    sin2 = sub(point(1),square(c))
    W = add(square(sub(square(r),ag2)),mul(delta,mul(ag2,sin2)))
    F = add(sub(sigma,mul(point(2),mul(t['M'],r))),
            add(mul(point(2),mul(t['B'],mul(t['M'],mul(r,mul(t['ag'],sin2))))),
                div(mul(square(t['B']),mul(W,sin2)),point(4))))
    Pnum = mul(add(mul(point(2),mul(t['M'],mul(r,t['ag']))),div(mul(t['B'],W),point(2))),sin2)
    angular = div(mul(F,sub(mul(point(398908.5),F),mul(point(1280000),Pnum))),
                  mul(rho,mul(t['E'],square(sigma))))
    connection_factor = div(t['L'],mul(point(2),mul(t['sqrt2'],t['vstar'])))
    return dict(speed=div(speed,t['vstar']), cu=mul(connection_factor,add(vrho,vz)),
                cv=mul(connection_factor,sub(vrho,vz)), mass=angular)


def rotate(field, axis, sqrt2):
    a, b = field
    if axis == 'u':
        return [div(add(a,b),sqrt2), div(sub(a,b),sqrt2)]
    return [div(sub(b,a),sqrt2), div(add(a,b),sqrt2)]


def trial_fields(grid, a):
    basis, derivative = grid['basis'], grid['derivative']
    values, du, dv = [], [], []
    for spin in (0,1):
        projected = matmul(basis,point(a[:,:,spin]))
        dprojected = matmul(derivative,point(a[:,:,spin]))
        values.append(matmul(projected,trans(basis)))
        du.append(matmul(dprojected,trans(basis)))
        dv.append(matmul(projected,trans(derivative)))
        print(json.dumps(dict(stage='interval_contractions', spin=spin)),flush=True)
    chi1, dchi1 = cutoff_grid(grid,7,9)
    chi, chiu, chiv = outer(chi1,chi1),outer(dchi1,chi1),outer(chi1,dchi1)
    cut = [mul(chi,f) for f in values]
    cut_u = [add(mul(chi,f),mul(chiu,g)) for f,g in zip(du,values)]
    cut_v = [add(mul(chi,f),mul(chiv,g)) for f,g in zip(dv,values)]
    return cut,cut_u,cut_v


def independent_geometry_control(target, metadata):
    """Pointwise implementation control from the original separate N/P and gauge."""
    from mpmath.ctx_mp import MPContext
    mp = MPContext(); mp.dps = 100
    def exact_float(x):
        n,d = float(x).as_integer_ratio()
        return mp.mpf(n)/d
    Q,k = mp.mpf(1)/5,mp.mpf(1)/2
    M = (mp.sqrt(65)-5)/20
    ag = mp.sqrt(mp.mpf(1)/4-M*M)
    B = (k-ag)/(M*(M+k))
    constants = {name:exact_float(float.fromhex(value)) for name,value in metadata['coordinate_constants'].items()}
    L,vs,rs = (constants[name] for name in ('L','vstar','rho_star'))
    E = vs/L
    parameters = dict(k=k,Q=Q,M=M,ag=ag,B=B,L=L,E=E,vstar=vs,rho_star=rs,sqrt2=mp.sqrt(2))
    for name,value in parameters.items():
        assert exact_float(target[name][0]) <= value <= exact_float(target[name][1]),name
    def metric(rho,z):
        dp,dm = mp.sqrt(rho*rho+(z+k)**2),mp.sqrt(rho*rho+(z-k)**2)
        r,c = M+(dp+dm)/2,(dp-dm)/(2*k)
        delta = r*r-2*M*r-ag*ag
        sigma,s2 = r*r-ag*ag*c*c,1-c*c
        W = (r*r-ag*ag)**2+delta*ag*ag*s2
        N = (delta+ag*ag*s2+2*B*M*r*ag*s2+B*B*W*s2/4)/sigma
        P = N*sigma*sigma/(delta+k*k*s2)**2
        original_A = -(2*M*r*ag+B*W/2)*s2/(N*sigma)
        return N,P,-original_A/(2*Q)
    def velocity(rho,z):
        N,P,_ = metric(rho,z)
        return N/P
    widths = {name:mp.mpf(0) for name in ('speed','cu','cv','mass','vrho','vz')}
    points = [(0.,0.),(7.5,-8.25),(-6.25,5.75),(8.5,8.)]
    for uu,vv in points:
        got = geometry({'x':point(np.array([uu,vv]))},target)
        factor = div(mul(target['sqrt2'],target['vstar']),target['L'])
        got['vrho'],got['vz'] = mul(factor,add(got['cu'],got['cv'])),mul(factor,sub(got['cu'],got['cv']))
        rho,z = rs+L*(exact_float(uu)+exact_float(vv))/mp.sqrt(2),L*(exact_float(uu)-exact_float(vv))/mp.sqrt(2)
        N,P,A = metric(rho,z)
        vrho,vz = mp.diff(lambda x:velocity(x,z),rho),mp.diff(lambda x:velocity(rho,x),z)
        truth = dict(speed=N/P/vs,cu=L*(vrho+vz)/(2*mp.sqrt(2)*vs),
                     cv=L*(vrho-vz)/(2*mp.sqrt(2)*vs),
                     mass=N*N*(mp.mpf(797817)/2-512000*A)/(rho*E),vrho=vrho,vz=vz)
        for name,value in truth.items():
            lo,hi = exact_float(got[name][0][0,1]),exact_float(got[name][1][0,1])
            assert lo <= value <= hi,(uu,vv,name)
            widths[name] = max(widths[name],hi-lo)
    sx,sz = mp.matrix([[0,1],[1,0]]),mp.matrix([[1,0],[0,-1]])
    vector = mp.matrix([mp.mpf(3)/8,-mp.mpf(5)/4])
    for axis,matrix in (('u',(sx+sz)/mp.sqrt(2)),('v',(sx-sz)/mp.sqrt(2))):
        got,truth = rotate([point(.375),point(-1.25)],axis,target['sqrt2']),matrix*vector
        for i in range(2):
            assert exact_float(got[i][0]) <= truth[i] <= exact_float(got[i][1])
    return dict(passed=True,precision_digits=100,points=points,parameter_inclusions=10,
                coefficient_inclusions=24,rotation_inclusions=4,
                maximum_interval_widths={name:mp.nstr(value,12) for name,value in widths.items()},
                scope='Independent pointwise implementation control; not a global enclosure proof.')


def run():
    started = time.perf_counter()
    controls = dict(arithmetic=arithmetic_checks(),quadrature=quadrature_checks())
    calibration = json.loads(project_path('dihole_null_scaling_results.json').read_text())
    previous = json.loads(project_path('dihole_probe_residual_results.json').read_text())
    assert source_matches('dihole_probe_residual.py',previous['source_sha256'])
    assert source_matches('dihole_null_scaling.py',calibration['source_sha256'])
    assert hashlib.sha256(project_path('dihole_null_scaling_results.json').read_bytes()).hexdigest() == previous['input_sha256']['dihole_null_scaling_results.json']
    assert hashlib.sha256(project_path('dihole_probe_residual_fields.npz').read_bytes()).hexdigest() == previous['field_sha256']
    remainder = remainder_checks()
    epsilon = remainder['remainder_components'][2]['squared_residual_integral_error']['upper']
    with np.load(project_path('dihole_probe_residual_fields.npz')) as fields:
        a = fields['wide_q512000_a']
        dual = -a.transpose(1,0,2)*np.array([1.,-1.])
        old_dual_difference = float(np.linalg.norm(dual-fields['wide_q512000_b']))
    assert a.shape == (128,128,2) and np.isfinite(a).all()
    target, target_metadata = exact_target(calibration,previous)
    controls['geometry'] = independent_geometry_control(target,target_metadata)
    grid = interval_grid(order=72,count=128,digits=75)
    print(json.dumps(dict(stage='validated_grid',seconds=time.perf_counter()-started)),flush=True)
    fields, du, dv = trial_fields(grid,a)
    coeff = geometry(grid,target)
    chi_source,_ = cutoff_grid(grid,5,7)
    phi0 = grid['basis'][0][:,0],grid['basis'][1][:,0]
    source = outer(mul(chi_source,phi0),mul(chi_source,phi0))
    au, av = rotate(fields,'u',target['sqrt2']),rotate(fields,'v',target['sqrt2'])
    kinetic_u, kinetic_v = rotate(du,'u',target['sqrt2']),rotate(dv,'v',target['sqrt2'])
    jfield = [sub(point(0),fields[1]),fields[0]]
    residual = []
    for spin in (0,1):
        kinetic = mul(coeff['speed'],add(kinetic_u[spin],kinetic_v[spin]))
        connection = add(mul(coeff['cu'],au[spin]),mul(coeff['cv'],av[spin]))
        r = sub(add(sub(sub(point(0),kinetic),connection),mul(coeff['mass'],jfield[spin])),fields[spin])
        residual.append(add(source,r) if spin==0 else r)
    # S=sigma_z Pz anti-commutes with K on compact fields; S h=h.
    # With b=-Sa, the adjoint residual is S*r and has identical norm.
    b = [sub(point(0),trans(fields[0])),trans(fields[1])]
    center_integrand = sub(mul(source,fields[0]),add(mul(b[0],residual[0]),mul(b[1],residual[1])))
    weight = outer(grid['weights'],grid['weights'])
    integrate = lambda f: sum_interval(mul(weight,f))
    integral = dict(residual_squared=integrate(add(square(residual[0]),square(residual[1]))),
                    corrected_center=integrate(center_integrand),source_norm_squared=integrate(square(source)))
    errors = point(epsilon)
    total = {name: (sub(value,errors)[0],add(value,errors)[1]) for name,value in integral.items()}
    assert total['residual_squared'][0]>0 and total['source_norm_squared'][0]>0
    r2upper = float(total['residual_squared'][1])
    center = total['corrected_center']
    imaginary_G = sub(center,point(r2upper))[0],add(center,point(r2upper))[1]
    physical_G = div(imaginary_G,target['E'])
    normalized_G = div(imaginary_G,total['source_norm_squared'])
    reference = previous['refined_corrected']['wide_q512000']
    tests = dict(input_provenance=True, exact_balanced_target_inside_analytic_boxes=True,
        exact_coordinate_scaling=True, compact_support_inside_physical_domain=True,
        validated_distinct_Gauss_nodes=True, interval_Hermite_contractions=True,
        positive_residual_and_source_integrals=True,
        saved_estimate_inside_new_response_projection=bool(imaginary_G[0] < reference['corrected_imaginary_G'] < imaginary_G[1]),
        useful_radius=bool(r2upper < 3e-6),
        source_normalization_resolved=bool(total['source_norm_squared'][0]>.999999999))
    assert all(tests.values()),tests
    names=['dihole_interval_arithmetic.py','dihole_interval_quadrature.py','dihole_quadrature_remainder.py',
           'dihole_null_scaling_results.json','dihole_probe_residual_results.json',
           'dihole_probe_residual_fields.npz','dihole_probe_residual.py','project_paths.py']
    result=dict(passed=True,tests=tests,target=target_metadata,controls=controls,
        runtime_versions=dict(numpy=np.__version__,mpmath=mpmath.__version__),
        quadrature=dict(order_per_width2_panel=72,total_nodes_per_axis=648,
            analytic_absolute_remainder_per_integral=epsilon,metadata=grid['metadata']),
        trial=dict(primal='Exact saved wide_q512000_a dyadic Hermite coefficients, cutoff7/9.',
            dual='NEW exact b(u,v)=-sigma_z*a(v,u), then identical cutoff7/9.',
            old_dual_coefficient_difference_norm=old_dual_difference,
            source='h=chi_source*f0, source cutoff5/7, unnormalized.',
            symmetry='S K S=-K, S h=h, b=-Sa => adjoint residual=S*primal residual.'),
        quadrature_sum_enclosures={key:export(value) for key,value in integral.items()},
        integral_enclosures={key:export(value) for key,value in total.items()},
        response=dict(center_imaginary_interval=export(center),dual_disk_radius_upper=r2upper,
            imaginary_G_interval=export(imaginary_G),physical_imaginary_G_interval=export(physical_G),
            unit_source_imaginary_G_interval=export(normalized_G),
            solution_error_upper=float(sqrt(total['residual_squared'])[1]),
            interpretation='Exact response lies within r2upper of some i*C with C in center interval; imaginary_G is its imaginary-axis projection, not an asserted real-axis scattering probability.'),
        scope=['One compact source in one azimuthal sector q512000,mu398908.5 at scaled kappa1.',
            'Natural uncut exact balanced exterior domain from docs/DIHOLE_DIRAC_DOMAIN.md.',
            'Elementary arithmetic uses documented IEEE754 interval assumptions; no unbounded BLAS rounding estimate.',
            'No vacuum state, ultraviolet sum, renormalized stress, backreaction, throat replacement or backward travel.'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={name:hashlib.sha256(project_path(name).read_bytes()).hexdigest() for name in names},
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_probe_certificate_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,response=result['response'],seconds=result['runtime_seconds']),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,5.4),gridspec_kw={'width_ratios':[1,1.25]})
    left.axis('off')
    stages=[('Exact balanced geometry',.87,'#247b62'),
            ('Validated nodes and quantum trial',.66,'#247b62'),
            ('Rounding + integration bounds',.45,'#247b62'),
            ('One response enclosed',.24,'#235b9f')]
    for text,y,color in stages:
        left.text(.5,y,text,ha='center',va='center',transform=left.transAxes,
                  bbox=dict(boxstyle='round,pad=.6',facecolor='white',edgecolor=color,lw=2),color=color)
        if y>.24:
            left.annotate('',xy=(.5,y-.14),xytext=(.5,y-.055),xycoords='axes fraction',
                          arrowprops=dict(arrowstyle='->',color='#7b828a'))
    left.text(.5,.045,'Still open: quantum stress and stable throat',ha='center',
              transform=left.transAxes,color='#6a4a37',fontsize=10)
    r=result['response']; lower,upper=r['imaginary_G_interval']; cl,ch=r['center_imaginary_interval']
    offset=.4190046
    scale=lambda x:(x-offset)*1e6
    right.plot([scale(lower),scale(upper)],[.75,.75],color='#235b9f',lw=10,solid_capstyle='butt')
    right.plot([scale(lower),scale(upper)],[.75,.75],'|',color='#123a6c',ms=19)
    right.plot([scale(cl),scale(ch)],[.30,.30],color='#247b62',lw=7)
    right.plot([scale((cl+ch)/2)],[.30],'o',color='#247b62',ms=6)
    right.set(yticks=[.30,.75],yticklabels=['Corrected trial center','Certified response'],
              ylim=(-.10,1.13),xlim=(-3.5,3.5),
              xlabel=r'Imaginary response minus 0.4190046 ($\times10^{-6}$)',
              title='Full bound for the compact probe')
    right.text(.5,.98,f'[{lower:.10f}, {upper:.10f}]',ha='center',va='top',transform=right.transAxes,
               color='#235b9f',fontsize=13)
    right.grid(axis='x',alpha=.16)
    fig.suptitle('Curved exterior: the response now has a complete numerical enclosure',fontsize=14,y=.98)
    fig.text(.5,.015,'One source, one angular mode, one imaginary frequency. No traversable medium or backward passage established.',
             ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.07,1,.93))
    for extension in ('png','pdf'):
        fig.savefig(project_path('dihole_probe_certificate_progress.'+extension),dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()
