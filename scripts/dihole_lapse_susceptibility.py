"""Certified local lapse derivatives of the frozen119 compact Dirac response.

The spatial metric, magnetic connection, source, Hilbert map and energy scale
stay fixed. N_e=N(1+e*f), with smooth compact 0<=f<=1, gives the first-order
Dirac identity H_e=H+e*{f,H}/2. Bounded overlaps and the existing solution-error
certificate enclose G'(0); no new PDE solve or supporting stress is asserted.
See docs/EXTERIOR_STRESS_BRIDGE.md for the domain and quantum-state distinction.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import mpmath
from mpmath.ctx_iv import MPIntervalContext

from dihole_interval_arithmetic import point, add, sub, mul, div, square, sqrt, sum_interval, matmul
from dihole_interval_quadrature import interval_grid
from dihole_probe_certificate import trans, outer, export, cutoff_grid
from dihole_quadrature_remainder import dyadic, outward
from project_paths import project_path


def finite_control():
    """Independent finite-matrix congruence derivative, including noncommuting f."""
    H = np.array([[0., -1j], [1j, 0.]])
    h, eye = np.array([1., 0.]), np.eye(2)
    cases = []
    for kappa in (1., 2.):
        psi = np.linalg.solve(H-1j*kappa*eye, h)
        chi = np.linalg.solve(H+1j*kappa*eye, h)
        for diagonal in ([1., 0.], [1., 1.]):
            f = np.diag(diagonal)
            V = (f@H+H@f)/2
            overlap = -(np.vdot(chi, f@h)+np.vdot(h, f@psi))/2-1j*kappa*np.vdot(chi, f@psi)
            derivative = -np.vdot(chi, V@psi)
            factor = sum(diagonal)
            exact = -1j*factor*kappa/(1+kappa*kappa)**2
            A, B = (psi/1j).real, (chi/1j).real
            center = -(h@f@(A-B)/2+kappa*B@f@A)
            assert max(abs(overlap-exact), abs(derivative-exact), abs(1j*center-exact)) < 1e-15
            def response(e):
                root = np.diag(np.sqrt(1+e*np.asarray(diagonal)))
                return np.vdot(h, np.linalg.solve(root@H@root-1j*kappa*eye, h))
            errors = [abs((response(e)-response(-e))/(2*e)-exact) for e in (.001, .0005)]
            if kappa == 1 and diagonal == [1., 1.]:
                # The leading central-difference error cancels in this case.
                assert max(errors) < 1e-12
            else:
                assert errors[1] < .26*errors[0] and errors[1] < 1e-7
            cases.append(dict(kappa=kappa, mask=diagonal, exact_imaginary_derivative=exact.imag,
                              finite_difference_errors=errors))
    return dict(passed=True, cases=cases,
                scope='Checks derivative signs and pairings. Exact affine dependence is a differential-operator identity, not a finite-matrix identity.')


def gaussian_1d(x, center):
    """Interval exp at the enclosed real nodes, using an isolated context."""
    iv = MPIntervalContext(); iv.dps = 75
    lo, hi = np.empty_like(x[0]), np.empty_like(x[1])
    for j, (left, right) in enumerate(zip(*x)):
        xx = iv.mpf([dyadic(iv, left).a, dyadic(iv, right).b])
        value = outward(iv.exp(-(xx-center)**2/8))
        lo[j], hi[j] = value['lower'], value['upper']
    return lo, hi


def run():
    started = time.perf_counter()
    name = 'dihole_probe_certificate_results.json'
    assert hashlib.sha256(project_path(name).read_bytes()).hexdigest() == 'd3c3e55a125c539baa8e9b5c41f97429e2759850a419c23280b70522e26f2843'
    saved = json.loads(project_path(name).read_text())
    inputs = dict(saved['input_sha256'])
    inputs['dihole_probe_certificate.py'] = saved['source_sha256']
    inputs[name] = hashlib.sha256(project_path(name).read_bytes()).hexdigest()
    for filename, expected in inputs.items():
        assert hashlib.sha256(project_path(filename).read_bytes()).hexdigest() == expected, filename
    assert saved['passed'] and all(saved['tests'].values())
    controls = dict(finite_matrix=finite_control())

    iv = MPIntervalContext(); iv.dps = 75
    constants = {key:dyadic(iv, float.fromhex(value))
                 for key,value in saved['target']['coordinate_constants'].items()}
    minimum_rho = constants['rho_star']-11*iv.sqrt(2)*constants['L']
    assert minimum_rho > iv.mpf('.484')
    F, field, source_bound = iv.exp(iv.mpf(9)/128), iv.mpf(9360000), iv.mpf(14)
    budgets = dict(derivative_center=F*(source_bound*field+field**2),
                   weighted_source_squared=F**2*source_bound**2,
                   weighted_trial_squared=F**2*field**2)
    assert F < iv.mpf('1.073') and all(value < iv.mpf('3e25') for value in budgets.values())
    remainder = saved['quadrature']['analytic_absolute_remainder_per_integral']
    assert 8*324*iv.mpf('3e25')*iv.mpf(2)**(1-2*72) <= dyadic(iv, remainder)
    energy_saved = saved['target']['parameter_enclosures']['E']
    energy = point(energy_saved['lower'])[0], point(energy_saved['upper'])[1]
    source_norm2 = tuple(np.asarray(v) for v in saved['integral_enclosures']['source_norm_squared'])
    error = point(saved['response']['solution_error_upper'])

    grid = interval_grid(order=72, count=128, digits=75)
    print(json.dumps(dict(stage='validated_grid', seconds=time.perf_counter()-started)), flush=True)
    with np.load(project_path('dihole_probe_residual_fields.npz')) as data:
        coefficients = data['wide_q512000_a']
    assert coefficients.shape == (128,128,2) and np.isfinite(coefficients).all()
    cutoff = cutoff_grid(grid,7,9)[0]
    cut = outer(cutoff,cutoff)
    basis = grid['basis']
    A = [mul(cut,matmul(matmul(basis,point(coefficients[:,:,spin])),trans(basis))) for spin in (0,1)]
    B = [sub(point(0),trans(A[0])), trans(A[1])]
    source_cutoff = cutoff_grid(grid,5,7)[0]
    phi0 = basis[0][:,0], basis[1][:,0]
    h = outer(mul(source_cutoff,phi0),mul(source_cutoff,phi0))
    weight = outer(grid['weights'],grid['weights'])
    def integrate(field):
        value = sum_interval(mul(weight,field))
        return sub(value,point(remainder))[0], add(value,point(remainder))[1]
    computed_source_norm = integrate(square(h))
    assert computed_source_norm[0] <= source_norm2[1] and source_norm2[0] <= computed_source_norm[1]
    # kappa=1. The exact dual is reflected, so S-even masks have ||fB||=||fA||.
    density = sub(point(0),add(mul(point(.5),mul(h,sub(A[0],B[0]))),
                              add(mul(B[0],A[0]),mul(B[1],A[1]))))
    rows = []
    for center in (-2,0,2):
        gaussian = gaussian_1d(grid['x'],center)
        f = outer(gaussian,gaussian)
        D = integrate(mul(f,density))
        fh2 = integrate(square(mul(f,h)))
        fA2 = integrate(add(square(mul(f,A[0])),square(mul(f,A[1]))))
        assert fh2[0] > 0 and fA2[0] > 0
        radius = add(mul(sqrt(fh2),error),add(mul(point(2),mul(sqrt(fA2),error)),square(error)))
        radius_upper = float(radius[1])
        enclosure = sub(D,point(radius_upper))[0], add(D,point(radius_upper))[1]
        location = constants['rho_star']+iv.sqrt(2)*constants['L']*center
        rows.append(dict(center_uv=[center,center], center_rho=outward(location),
            imaginary_center_interval=export(D), weighted_source_norm_squared=export(fh2),
            weighted_primal_and_dual_norm_squared=export(fA2), disk_radius_upper=radius_upper,
            imaginary_derivative_interval=export(enclosure),
            physical_imaginary_derivative_interval=export(div(enclosure,energy)),
            unit_source_imaginary_derivative_interval=export(div(enclosure,source_norm2))))
    tests = dict(exact_input_provenance=True, derivative_pairing_control=True,
        smooth_mask_support_inside_exterior=True, analytic_overlap_budgets=True,
        source_norm_agrees_with_certificate=True, positive_weighted_norm_integrals=True,
        three_finite_derivative_enclosures=all(np.isfinite(row['imaginary_derivative_interval']).all() for row in rows))
    assert all(tests.values()),tests
    result = dict(passed=True, tests=tests, controls=controls, responses=rows,
        variation=dict(lapse='N_e=N*(1+e*f), |e|<1; fixed spatial metric, connection, source and baseline energy.',
            mask='f_c=zeta(u)zeta(v)exp(-((u-c)^2+(v-c)^2)/8), c=-2,0,2.',
            smooth_taper='zeta=1 on [-9,9], 0 outside [-11,11], 0<=zeta<=1; explicit flat smoothstep in EXTERIOR_STRESS_BRIDGE.md.',
            minimum_support_rho=outward(minimum_rho), operator='H_e=H+e*(fH+Hf)/2 exactly for this first-order differential operator.',
            derivative_center='D=-integral f*(h^T*(A-B)/2+B^T*A), kappa=1.',
            error_radius='||fh||*e_field+2*||fA||*e_field+e_field^2; ||fB||=||fA|| by reflection.',
            solution_error_upper=saved['response']['solution_error_upper'],
            interpretation='G_prime lies within disk_radius of some i*D in the center interval; displayed bounds project that disk onto the imaginary axis.'),
        quadrature=dict(order_per_panel=72, nodes_per_axis=648, analytic_remainder_per_integral=remainder,
            complex_gaussian_modulus=outward(F), integrand_modulus_bounds={k:outward(v) for k,v in budgets.items()},
            source_norm_control=export(computed_source_norm)),
        scope=['One compact source, one chirality/azimuthal sector, scaled imaginary frequency kappa=1.',
            'A metric derivative of a packet response, not a ground-state expectation, energy density, null stress or return probability.',
            'No supporting throat, new physical law, material recipe or backward passage.',
            'Reuses119 IEEE754 and mpmath interval arithmetic assumptions and full solution-error bound.'],
        runtime_versions=dict(numpy=np.__version__,mpmath=mpmath.__version__),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), input_sha256=inputs,
        runtime_seconds=time.perf_counter()-started)
    project_path('dihole_lapse_susceptibility_results.json').write_text(json.dumps(result,indent=2)+'\n')
    draw(result)
    print(json.dumps(dict(passed=True,responses=rows,seconds=result['runtime_seconds']),indent=2),flush=True)
    return result


def draw(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,(left,right)=plt.subplots(1,2,figsize=(11.5,5.4),gridspec_kw={'width_ratios':[1,1.35]})
    colors=['#276c9d','#c78019','#357b55']
    axis=np.linspace(-7,7,201); u,v=np.meshgrid(axis,axis)
    for row,color in zip(result['responses'],colors):
        c=row['center_uv'][0]
        g=np.exp(-((u-c)**2+(v-c)**2)/8)
        left.contour(u,v,g,levels=[.3,.7],colors=[color],linewidths=[1,2])
        left.plot(c,c,'o',color=color)
        left.text(c+.2,c+.2,f'c={c}',color=color)
    left.set(aspect='equal',xlabel='Local coordinate u',ylabel='Local coordinate v',
             title='Three local clock-rate variations')
    for j,(row,color) in enumerate(zip(result['responses'],colors)):
        lo,hi=row['imaginary_derivative_interval']; cl,ch=row['imaginary_center_interval']
        right.plot([lo,hi],[j,j],color=color,lw=8,solid_capstyle='butt')
        right.plot([lo,hi],[j,j],'|',color=color,ms=17)
        right.plot([(cl+ch)/2],[j],'o',color='white',mec=color,ms=6)
    right.axvline(0,color='#888888',lw=1)
    right.set(yticks=[0,1,2],yticklabels=['c = -2','c = 0','c = 2'],ylim=(-.6,2.6),
              xlabel=r'Imaginary part of $dG/d\epsilon$ at $\epsilon=0$',title='Response derivative: enclosed intervals')
    right.grid(axis='x',alpha=.15)
    fig.suptitle('Testing how a quantum probe responds to a local change in lapse',fontsize=14,y=.98)
    fig.text(.5,.025,'Contours show mathematical perturbations, not material. Derivative signs do not establish negative energy or a time-travel medium.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.93))
    for extension in ('png','pdf'):
        fig.savefig(project_path('dihole_lapse_susceptibility_progress.'+extension),dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    run()
