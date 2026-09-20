"""A massive radial mode on the saved smooth large-flux collar.

This compares two auxiliary self-adjoint circle holonomies on identical
local geometry.  It does not solve the physical two-mouth return region.
The full frequency bound is analytic within that stated radial model;
its background integrals are evaluated numerically, not by interval arithmetic.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from numpy.polynomial.legendre import leggauss

from dirac_return_closure import collar


OLD_COLLAR_SHA = '595b0d4578bc13f828020b956a2d5ce41e0e11b12de47df4e0d5a19119829140'
W = 1.6194435263649365
ALPHA = 1e-9
Q = 100000
WIDTH = 1.


def bump(l):
    """Peak-one smooth compact proper-distance weight, |l|<1 r_e."""
    u = (l/WIDTH)**2
    if u >= 1.:
        return 0., 0.
    value = math.exp(1.-1./(1.-u))
    return value, -2.*l*value/(WIDTH**2*(1.-u)**2)


def build_profile(refinement=1):
    record, raw = collar(W, alpha=ALPHA, refinement=refinement, return_profile=True)
    sol = raw['solution']
    end = float(sol.t[-1])

    def fields(l):
        state = sol.sol(abs(l))
        sign = -1. if l < 0 else 1.
        return (1.+raw['central_radius_offset']+ALPHA*state[0],
                math.exp(state[2]), sign*ALPHA*state[1], sign*state[3])

    return dict(record=record, fields=fields, end=end, refinement=refinement, solution=sol)


def profile_integrals(profile, n=1, q=Q, eta=0.):
    lam = math.sqrt(n*(n+q))
    end = profile['end']
    fields = profile['fields']
    opts = dict(epsabs=3e-11, epsrel=3e-12, limit=250)
    def integrate(fun, cutoff=end):
        return 2.*quad(fun, 0., cutoff, points=[WIDTH] if cutoff>WIDTH else None, **opts)[0]
    length = integrate(lambda l: math.exp(-2*eta*bump(l)[0])/fields(l)[1])
    barrier = integrate(lambda l: lam*math.exp(-eta*bump(l)[0])/fields(l)[0])
    I = integrate(lambda l: bump(l)[0]*math.exp(-2*eta*bump(l)[0])/fields(l)[1], WIDTH)
    J = integrate(lambda l: lam*bump(l)[0]*math.exp(-eta*bump(l)[0])/fields(l)[0], WIDTH)
    G = 4.*math.pi*ALPHA/Q
    def geom(l):
        R, N, B, A = fields(l)
        f, fp = bump(l)
        return N*B*((2*A*R+B)*f+R*fp)/G
    def source(l):
        R, N, B, A = fields(l)
        f, fp = bump(l)
        z = A-B/R
        return Q/(12*math.pi)*N*(-fp*z-f*z*z-f*W/(N*N))
    return dict(optical_length=length, mass_barrier=barrier,
                probe_optical_integral=I, probe_mass_integral=J,
                geometric_flare_probe=integrate(geom, WIDTH),
                retained_action_flare_probe=integrate(source, WIDTH),
                lambda_n=lam, G_over_re_squared=G)


def radial_mode(kappa, profile, integrals, n=1, q=Q, eta=0., rtol=2e-9):
    """Positive-imaginary-frequency scattering, with log amplitude retained.

    z=1-y is stable near reflection -i.  h=log(t)+integral(m dl)
    avoids subtracting the large mass barrier in the source integrand.
    The last two components differentiate z,h under the compact null probe.
    """
    if not math.isfinite(kappa) or kappa < 0 or n < 0 or q <= 0:
        raise ValueError('Finite kappa>=0, n>=0, q>0 required')
    lam = math.sqrt(n*(n+q))
    fields = profile['fields']
    end = profile['end']
    def coefficients(l):
        R, N, _, _ = fields(l)
        f = bump(l)[0]
        return kappa/N*math.exp(-2*eta*f), lam/R*math.exp(-eta*f), f
    def rhs(l, state):
        z, h, dz, dh = state
        k, m, f = coefficients(l)
        a = -2*k-2*m*(1-z)
        return [2*k*(1-z)-m*z*(2-z), -k+m*z,
                a*dz+f*(-4*k*(1-z)+m*z*(2-z)),
                2*f*k-f*m*z+m*dz]
    def jac(l, state):
        z, h, dz, dh = state
        k, m, f = coefficients(l)
        a = -2*k-2*m*(1-z)
        return np.array([[a,0,0,0], [m,0,0,0],
                         [2*m*dz+f*(4*k+2*m*(1-z)),0,a,0],
                         [-f*m,0,m,0]])
    # Split at the compact support, avoiding a large accepted step skipping it.
    state = np.array([1., 0., 0., 0.])
    nfev = 0
    for left, right in ((-end,-WIDTH),(-WIDTH,0.),(0.,WIDTH),(WIDTH,end)):
        solved = solve_ivp(rhs,(left,right),state,method='Radau',jac=jac,
                           rtol=rtol,atol=np.array([rtol*1e-5,rtol*.05,rtol*1e-4,rtol*.05]))
        if not solved.success:
            raise RuntimeError('Smooth-collar radial solve: '+solved.message)
        state = solved.y[:,-1]
        nfev += solved.nfev
    z, h, dz, dh = state
    y = 1.-z
    log_t = -integrals['mass_barrier']+h
    # At nonzero eta this tangent is not used for the production derivative.
    dlog_t = dh+integrals['probe_mass_integral']
    dy = -dz
    with mp.workdps(60):
        lt = mp.mpf(str(log_t)); t2 = mp.exp(2*lt)
        C = 1+mp.mpf(str(y))**2+t2
        log_a = mp.log(2)+lt-mp.log(C)  # a=2/Tr(M)
        dlog_X = (2*mp.mpf(str(y))*dy+2*t2*dlog_t)/C-dlog_t
        a2 = mp.exp(2*log_a)
        log_contrast = mp.log(2*mp.atanh(mp.exp(log_a)))
        if dlog_X:
            log_source = mp.log(2/mp.pi)+log_a+mp.log(abs(dlog_X))-mp.log1p(-a2)
            source_sign = int(mp.sign(dlog_X))
        else:
            log_source = mp.ninf
            source_sign = 0
    return dict(kappa=kappa, reflection_imaginary=-y, one_minus_y=z,
                log_transmission=log_t, log_transmission_plus_barrier=h,
                pressure_dlog_transmission=dlog_t,
                pressure_dreflection_imaginary=-dy,
                log_abs_logdet_phase_contrast=float(log_contrast),
                pressure_dlog_half_transfer_trace=float(dlog_X),
                source_integrand_sign=source_sign,
                log_abs_source_integrand=float(log_source) if source_sign else None,
                scaled_source_integrand=(source_sign*math.exp(float(log_source)+integrals['mass_barrier'])
                                         if source_sign else 0.),
                transmission_is_below_float_range=bool(log_t<math.log(np.finfo(float).tiny)),
                nfev=nfev,rtol=rtol,eta=eta)


def frequency_bound(integrals, multiplicity):
    """Whole kappa integral bound from two positive transfer-matrix bases."""
    with mp.workdps(60):
        L, barrier, I, J = [mp.mpf(str(integrals[key])) for key in
            ('optical_length','mass_barrier','probe_optical_integral','probe_mass_integral')]
        log_sech = mp.log(2)-barrier-mp.log1p(mp.exp(-2*barrier))
        K = (barrier+mp.log1p(mp.exp(-2*barrier)))/L
        bracket = I*K*K+J*K+2*I*K/L+2*I/L**2+J/L
        log_bound = mp.log(2*multiplicity/mp.pi)+log_sech-2*mp.log(mp.tanh(barrier))+mp.log(bracket)
        ratio = log_bound-mp.log(abs(mp.mpf(str(integrals['geometric_flare_probe']))))
        return dict(frequency_split=float(K),log_abs_bound=float(log_bound),
                    log10_abs_bound=float(log_bound/mp.log(10)),
                    log10_bound_over_geometric_probe=float(ratio/mp.log(10)),
                    finite_angular_multiplicity=multiplicity,
                    bound_as_decimal=mp.nstr(mp.exp(log_bound),18))


def pilot():
    started = time.perf_counter()
    profile = build_profile(1)
    ints = profile_integrals(profile)
    print(json.dumps(dict(stage='profile',seconds=time.perf_counter()-started,integrals=ints)),flush=True)
    rows = []
    for kappa in (0.,1.):
        begin = time.perf_counter()
        rows.append(radial_mode(kappa,profile,ints))
        print(json.dumps(dict(stage='radial',seconds=time.perf_counter()-begin,**rows[-1])),flush=True)
    print(json.dumps(dict(stage='bound',result=frequency_bound(ints,2*(Q+2)))),flush=True)
    return dict(seconds=time.perf_counter()-started,rows=rows)


def _direct_transfer_control(profile):
    """Independent unfactored linear transfer at a smaller, finite barrier.

    q=1 here is solely an ODE control on the same geometry, not a new
    semiclassical background or a physical continuation of the q=1 model.
    """
    kappa, lam = 1., math.sqrt(2.)
    def rhs(l, flat):
        R,N,_,_ = profile['fields'](l)
        generator=np.array([[-kappa/N,lam/R],[lam/R,kappa/N]])
        return (generator@flat.reshape(2,2)).ravel()
    direct=solve_ivp(rhs,(-profile['end'],profile['end']),np.eye(2).ravel(),
                     method='DOP853',rtol=2e-11,atol=1e-12)
    if not direct.success:
        raise RuntimeError(direct.message)
    matrix=direct.y[:,-1].reshape(2,2)
    ints=profile_integrals(profile,q=1)
    scalar=radial_mode(kappa,profile,ints,q=1,rtol=2e-10)
    return dict(operator_q=1,geometry_q=Q,kappa=kappa,
                log_transmission_error=abs(scalar['log_transmission']+math.log(matrix[1,1])),
                reflection_error=abs(-scalar['reflection_imaginary']-matrix[0,1]/matrix[1,1]),
                log_half_trace_error=abs(-scalar['log_abs_logdet_phase_contrast']+math.log(2.)
                                         -math.log(np.trace(matrix)/2.)),
                nfev=direct.nfev,scope='Independent linear transfer ODE control, not an additional physical candidate.')


def _bound_proof_control():
    """Independent finite transfer-products test of the positive-basis lemma."""
    from scipy.linalg import expm
    masses=np.array([.2,1.7,.6,2.3]); clocks=np.array([1.4,.3,1.2,.7])
    lengths=np.array([.3,.2,.5,.4]); weights=np.array([0.,.8,1.,.1])
    def transfer(eta):
        result=np.eye(2)
        for m,k,d,f in zip(masses,clocks,lengths,weights):
            m*=math.exp(-eta*f);k*=math.exp(-2*eta*f)
            result=expm(d*np.array([[-k,m],[m,k]]))@result
        return result
    M=transfer(0.);X=np.trace(M)/2.
    derivative=(math.log(np.trace(transfer(1e-5)))-math.log(np.trace(transfer(-1e-5))))/2e-5
    rhs=float(np.sum(weights*(masses+2*clocks)*lengths))
    return dict(half_trace=float(X),mass_lower_bound=math.cosh(float(masses@lengths)),
                clock_lower_bound=math.cosh(float(clocks@lengths)),
                absolute_log_trace_derivative=abs(derivative),derivative_bound=rhs,
                determinant_error=abs(float(np.linalg.det(M))-1.))


def checks():
    started=time.perf_counter()
    coarse=build_profile(1); fine=build_profile(2)
    ci=profile_integrals(coarse); fi=profile_integrals(fine)
    saved=json.loads(Path('dirac_return_closure_results.json').read_text(encoding='utf-8'))
    expected=saved['records'][-1]['roots'][-1]['background']
    default=collar(W,alpha=ALPHA,refinement=2)
    background_equal=default==fine['record']==expected
    old=Path('dirac_return_closure.py').read_bytes()
    reverted=old.replace(b', *, return_profile=False):',b'):').replace(
        b"    record={'clock':clock",b"    return {'clock':clock").replace(
        b"    if return_profile:\n        return record, {'solution':sol,'central_radius_offset':offset,'alpha':alpha}\n    return record\n",b'')
    preservation=hashlib.sha256(reverted).hexdigest()==OLD_COLLAR_SHA
    print(json.dumps(dict(stage='background',default_saved_exact=background_equal,
                          old_source_restored_exact=preservation,seconds=time.perf_counter()-started)),flush=True)

    scale=math.sqrt(fi['lambda_n'])
    frequencies=(0.,1.,scale,4*scale)
    pairs=[]
    for kappa in frequencies:
        a=radial_mode(kappa,coarse,ci,rtol=2e-9)
        b=radial_mode(kappa,fine,fi,rtol=2e-10)
        pairs.append(dict(kappa=kappa,coarse=a,fine=b,
                          log_transmission_difference=abs(a['log_transmission']-b['log_transmission']),
                          source_log_difference=abs(a['log_abs_source_integrand']-b['log_abs_source_integrand'])))
    zero=pairs[0]['fine']
    with mp.workdps(60):
        barrier=mp.mpf(str(fi['mass_barrier']))
        zero_log=-barrier+mp.log(2)-mp.log1p(mp.exp(-2*barrier))
        zero_derivative=fi['probe_mass_integral']*mp.tanh(barrier)
    zero_errors=dict(log_transmission=abs(zero['log_transmission']-float(zero_log)),
                     pressure_log_transmission=abs(zero['pressure_dlog_transmission']-float(zero_derivative)))

    zi=profile_integrals(fine,n=0)
    massless=radial_mode(.4,fine,zi,n=0,rtol=2e-10)
    massless_errors=dict(log_transmission=abs(massless['log_transmission']+.4*zi['optical_length']),
                        pressure_log_transmission=abs(massless['pressure_dlog_transmission']-.8*zi['probe_optical_integral']),
                        reflection=abs(massless['reflection_imaginary']))
    L=zi['optical_length'];I=zi['probe_optical_integral']
    def zero_integrand(k):
        z=k*L
        if z>700:return 0.
        return -4*k*I/(math.pi*math.sinh(z)) if z else -4*I/(math.pi*L)
    zero_integral=quad(zero_integrand,0.,np.inf,epsabs=2e-12,epsrel=2e-12)[0]
    massless_errors['integrated_pressure']=abs(zero_integral+math.pi*I/L**2)

    derivative_rows=[]
    central=pairs[2]['fine']
    for step in (1e-4,5e-5):
        values=[]
        for eta in (-step,step):
            ii=profile_integrals(fine,eta=eta)
            values.append(radial_mode(scale,fine,ii,eta=eta,rtol=2e-10)['log_abs_logdet_phase_contrast'])
        actual=(values[1]-values[0])/(2*step)
        predicted=-central['pressure_dlog_half_transfer_trace']
        derivative_rows.append(dict(step=step,finite_variation=actual,tangent=predicted,
                                    relative_error=abs(actual-predicted)/max(1.,abs(predicted))))

    finite=[]; maximum=4*scale
    for order in (8,12):
        nodes,weights=leggauss(order); ks=(nodes+1)*maximum/2; ws=weights*maximum/2
        rows=[radial_mode(float(k),fine,fi,rtol=2e-9) for k in ks]
        value=float(np.dot(ws,[r['scaled_source_integrand'] for r in rows]))
        multiplicity=2*(Q+2)
        log_abs=math.log(abs(value))+math.log(multiplicity)-fi['mass_barrier']
        finite.append(dict(order=order,kappa_interval=[0.,maximum],
                           barrier_scaled_integral_per_paired_weyl=value,
                           full_first_multiplet_sign=int(np.sign(value)),
                           full_first_multiplet_log10_abs=log_abs/math.log(10),rows=rows))
    finite_error=abs(finite[0]['barrier_scaled_integral_per_paired_weyl']
                     -finite[1]['barrier_scaled_integral_per_paired_weyl'])/abs(finite[1]['barrier_scaled_integral_per_paired_weyl'])
    bound=frequency_bound(fi,2*(Q+2)); coarse_bound=frequency_bound(ci,2*(Q+2))
    direct=_direct_transfer_control(fine)
    lemma=_bound_proof_control()
    quadrature_barrier=abs(fi['mass_barrier']-2*fi['lambda_n']*fine['solution'].y[5,-1])
    quadrature_optical=abs(fi['optical_length']-2*fine['solution'].y[4,-1])
    flare_relative=abs(fi['geometric_flare_probe']-fi['retained_action_flare_probe'])/abs(fi['geometric_flare_probe'])
    tests=dict(default_saved84_profile_exact=background_equal,
               only_authorized_historical_api_extension=preservation,
               barrier_independent_ode_accumulator=quadrature_barrier<2e-6,
               optical_independent_ode_accumulator=quadrature_optical<2e-10,
               same_probe_geometric_source=flare_relative<2e-7,
               exact_zero_frequency_massive=max(zero_errors.values())<2e-7,
               exact_massless_controls=max(massless_errors.values())<2e-7,
               radial_background_refinement=max(p['source_log_difference'] for p in pairs)<2e-5,
               compact_finite_variation=max(r['relative_error'] for r in derivative_rows)<2e-6,
               finite_frequency_quadrature=finite_error<2e-3,
               independent_linear_transfer=max(direct[k] for k in ('log_transmission_error','reflection_error','log_half_trace_error'))<2e-7,
               transfer_positivity_lemma=(lemma['half_trace']>=max(lemma['mass_lower_bound'],lemma['clock_lower_bound'])
                    and lemma['absolute_log_trace_derivative']<=lemma['derivative_bound'] and lemma['determinant_error']<2e-12))
    tests={key:bool(value) for key,value in tests.items()}
    paths=('dirac_return_closure.py','dirac_return_closure_results.json','dirac_matched_source.py','dirac_matched_source_results.json')
    result=dict(passed=all(tests.values()),tests=tests,seconds=time.perf_counter()-started,
        conclusions=dict(first_mode_bound_below_geometric_budget=bound['log10_bound_over_geometric_probe']<0,
                         scope='Outcome of the specified auxiliary phase-contrast model, not a verification requirement or evidence for absolute physical support.'),
        parameters=dict(q=Q,n=1,alpha=ALPHA,w=W,outer_areal_radius=100.,
                        proper_probe_half_width=WIDTH,geometry_refinements=[1,2],
                        angular_degeneracy=Q+2,weyl_copies=2,
                        units='r_e=hbar=c=1; center clock N(0)=1; first-multiplet source includes both Weyl operators.'),
        profile=fi,coarse_profile=ci,radial_samples=pairs,
        zero_frequency_errors=zero_errors,massless_control_errors=massless_errors,
        massless_integrated_contrast=zero_integral,
        compact_finite_variations=derivative_rows,finite_frequency_integrals=finite,
        finite_frequency_relative_difference=finite_error,whole_frequency_bound=bound,
        bound_background_log10_difference=abs(bound['log10_abs_bound']-coarse_bound['log10_abs_bound']),
        independent_transfer_control=direct,positive_basis_control=lemma,
        independent_barrier_accumulator_error=quadrature_barrier,
        independent_optical_accumulator_error=quadrature_optical,
        geometric_retained_source_relative_difference=flare_relative,
        equations=dict(radial='H=-i sigma3 d/dx + sqrt(n(n+q))*N/R sigma1; dx=dl/N',
            boundary='f(+x_b)=exp(i theta) f(-x_b), theta=0 or pi, in the same transported angular/spin/gauge frame.',
            null_variation='N_eta=N exp(eta f), A_eta=exp(-eta f), R fixed; k_eta=k exp(-2 eta f), m_eta=m exp(-eta f).',
            phase_determinant='D_theta=C-2t cos(theta), C=1+t^2-rL*rR; log(D_pi/D_0)=2atanh(2t/C).',
            source='Delta Q=-(multiplicity/pi) integral_0^infinity d_kappa d_eta log(D_pi/D_0).',
            geometric_probe='Q_geom=-G^-1 integral N^2 R f (R_prime/N)_prime dl; evaluated by compact-support integration by parts.',
            retained_probe='Q_action=q/(12pi) integral N[-f_prime*z-f*z^2-f*w/N^2]dl, z=N_prime/N-R_prime/R.'),
        bound_proof=dict(variables='Lambda=int m dl, L=int dl/N, J=int m f dl, I=int f dl/N, X=Tr M/2, a=1/X.',
            positivity='Real cooperative generator [[-k,m],[m,k]] and its constant Hadamard/sign rotation each have nonnegative off-diagonal entries. Dyson positivity gives X>=cosh Lambda and X>=cosh(kappa L). Cyclic Duhamel traces bound |delta log X|<=2 kappa I+J.',
            envelope='a<=min(sech Lambda,2 exp(-kappa L)); |delta log(D_pi/D_0)|<=2a(2kappa I+J)/tanh(Lambda)^2.',
            integrated='K=log(2cosh Lambda)/L; bound=(2 multiplicity sech Lambda)/(pi tanh(Lambda)^2)*(I K^2+J K+2I K/L+2I/L^2+J/L).',
            interpretation='Analytic all-frequency bound for this finite first multiplet and these two specified unitary circle domains. Numerical background integration/rounding is checked by refinement, not interval-certified.'),
        limitations=[
            'The identical auxiliary S2 x S1 completion has continuous Lipschitz metric coefficients and a derivative jump at its seam outside the compact probe. It is not a physical smooth exterior or an Einstein-Maxwell junction; no absorbing boundary is used.',
            'The periodic and antiperiodic self-adjoint radial domains preserve the boundary current. They are specific holonomies, not a bound over arbitrary exterior scattering operators.',
            'Only n=1 is included. Its genuine q+2 degeneracy and two Weyl copies are retained; this is not the old q=1 answer multiplied by q.',
            'Same-geometry state differences cancel local counterterms/anomaly contributions. No invariant-product reference is used to assert absolute null stress on the varying-radius collar.',
            'The collar retains its assumed w and truncated83 bulk action. The physical zero-mode return, other angular channels, absolute completed source, formation and stability remain unsolved.',
            'The finite quadrature is only a finite-frequency contrast. Its omitted interval is not inferred from quadrature agreement; the separate analytic bound covers the complete integral.',
            'Exponentially small quantities are kept as logarithms. A floating-point zero is never reported as an exact vanishing response.'],
        provenance=dict(original_collar_source_sha256=OLD_COLLAR_SHA,
            original_default_bytes_recovered=preservation,saved_fixed_w_background_exact=background_equal,
            input_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}),
        sources=['https://arxiv.org/html/1807.04726v3#S2.SS3','https://arxiv.org/html/1807.04726v3#S5.SS2'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if not result['passed']:
        raise AssertionError(json.dumps(dict(failed=[k for k,v in tests.items() if not v],
            finite_frequency_relative_difference=finite_error,direct=direct,zero=zero_errors,
            massless=massless_errors,barrier=quadrature_barrier,flare=flare_relative,
            derivatives=derivative_rows,refinement=pairs),indent=2))
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pilot',action='store_true')
    args = parser.parse_args()
    if args.pilot:
        pilot()
    else:
        result=checks()
        destination=Path(__file__).with_name('smooth_collar_mode_results.json')
        destination.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(dict(passed=result['passed'],seconds=result['seconds'],tests=result['tests'],
                              whole_frequency_bound=result['whole_frequency_bound'])),flush=True)
