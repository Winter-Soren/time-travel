"""Finite-channel imaginary-frequency axial tail in the curved dihole.

This uses the actual Weyl metric and gauge field, projected on a small fixed
Melvin radial basis. It is not a two-dimensional exterior or throat map.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import eval_genlaguerre, roots_laguerre

from balanced_dihole import symbolic_fields
from melvin_leads import coefficients, spectrum, zero_profile


def radial_basis(p, pairs=3, size=32, order=96, omega=.8791333609701094):
    """Regular lower Ritz modes, differentiated upper partners, exact zero."""
    if (not isinstance(pairs, int) or pairs < 1 or not isinstance(size, int)
            or size < 4 or size <= pairs or not isinstance(order, int)
            or order < size+3 or order > 144 or not math.isfinite(omega) or omega <= 0):
        raise ValueError('Require pairs>=1, size>pairs with size>=4, size+3<=order<=144, and finite omega>0; exp-weight quadrature is capped at144.')
    mu = .5
    beta, flux = p['B']**2/4, p['q']/(p['Q']*p['B'])
    linear, cubic = coefficients(mu, beta, flux)
    values, vectors, residuals = spectrum(mu, -1, size, omega, beta, flux, count=pairs)
    # Fix the arbitrary real eigenvector signs across sizes/channel cutoffs.
    for j in range(pairs):
        pivot = int(np.argmax(abs(vectors[:, j])))
        if vectors[pivot, j] < 0:
            vectors[:, j] *= -1
    masses = np.sqrt(values)
    x, wx = roots_laguerre(order)
    rho = np.sqrt(x/omega)
    weights = wx*np.exp(x)/(2*np.sqrt(omega*x))
    normalizations = 1/np.sqrt(np.arange(size)+1.)

    def profiles(r):
        xx = omega*np.asarray(r)**2
        polynomials = np.column_stack([eval_genlaguerre(j, 1., xx) for j in range(size)])
        derivatives = np.column_stack([np.zeros_like(xx) if j == 0 else -eval_genlaguerre(j-1, 2., xx) for j in range(size)])
        poly = (polynomials*normalizations)@vectors
        dp = (derivatives*normalizations)@vectors
        base = math.sqrt(2*omega**2)*np.asarray(r)**1.5*np.exp(-xx/2)
        lower = base[:, None]*poly
        lower_prime = base[:, None]*((1.5/np.asarray(r)-omega*np.asarray(r))[:, None]*poly+2*omega*np.asarray(r)[:, None]*dp)
        w = mu/np.asarray(r)+linear*np.asarray(r)+cubic*np.asarray(r)**3
        upper = (lower_prime+w[:, None]*lower)/masses
        return lower, lower_prime, upper

    lower, lower_prime, upper_massive = profiles(rho)
    zero = zero_profile(mu, beta, flux)
    upper_zero = np.exp(mu*np.log(rho)+linear*rho*rho/2+cubic*rho**4/4-zero['log_norm']/2)
    upper = np.column_stack([upper_zero, upper_massive])
    W = mu/rho+linear*rho+cubic*rho**3
    gram_lower = lower.T@(weights[:, None]*lower)
    gram_upper = upper.T@(weights[:, None]*upper)
    C0 = upper.T@(weights[:, None]*(lower_prime+W[:, None]*lower))
    target = np.vstack([np.zeros(pairs), np.diag(masses)])
    # Canonical profiles are already normalized analytically. Do not conceal
    # their finite Ritz/quadrature errors by replacing the measured Gram data.
    h = 1e-5
    sample = np.array([.3, .8, 1.5, 2.4])
    a = profiles(sample)
    near = [profiles(sample+j*h)[0] for j in (-2, -1, 1, 2)]
    fd = (near[0]-8*near[1]+8*near[2]-near[3])/(12*h)
    derivative_error = float(np.max(abs(fd-a[1])))
    return dict(pairs=pairs, ritz_size=size, radial_order=order, omega=omega,
                mu=mu, beta=beta, flux=flux, masses=masses,
                lower=lower, lower_prime=lower_prime, upper=upper,
                rho=rho, weights=weights, W_Melvin=W,
                gram_lower=gram_lower, gram_upper=gram_upper,
                C_Melvin_measured=C0, C_Melvin_target=target,
                squared_operator_ritz_residuals=residuals,
                zero_profile=zero, lower_derivative_control=derivative_error)


def curved_coefficients(z, basis, p, field_function, melvin=False):
    rho = basis['rho']
    if melvin:
        vinverse = np.ones_like(rho)
        W = basis['W_Melvin']
    else:
        plus = np.hypot(rho, z+p['k'])
        minus = np.hypot(rho, z-p['k'])
        total = plus+minus
        r = p['M']+total/2
        cosine = 2*z/total  # Rationalized (dplus-dminus)/(2k).
        raw = field_function(r, cosine, p['M'], p['a'], p['B'])
        N, P, Ahat = [np.broadcast_to(np.asarray(v, float), rho.shape) for v in raw[:3]]
        A = p['physical_A_factor']*Ahat
        vinverse = P/N
        W = N*P*(basis['mu']-A)/rho
        if np.min(N) <= 0 or np.min(P) <= 0:
            raise ValueError('Axial projection left the positive-lapse exterior.')
    weights = basis['weights']
    fp, fm = basis['upper'], basis['lower']
    plus_matrix = fp.T@((weights*vinverse)[:, None]*fp)
    minus_matrix = fm.T@((weights*vinverse)[:, None]*fm)
    coupling = fp.T@(weights[:, None]*(basis['lower_prime']+W[:, None]*fm))
    return plus_matrix, minus_matrix, coupling


def pure_start(kappa, masses):
    if not math.isfinite(kappa) or kappa <= 0 or not np.all(np.isfinite(masses)) or np.min(masses) <= 0:
        raise ValueError('Require finite kappa>0 and finite positive massive eigenvalues.')
    return np.column_stack([np.zeros(len(masses)), np.diag(masses/(kappa+np.hypot(kappa, masses)))])


def leading_start_coefficient(kappa, basis, p):
    """Order1/z graph correction in this fixed, finite projected basis."""
    linear, cubic = coefficients(basis['mu'], basis['beta'], basis['flux'])
    rho = basis['rho']
    W1 = 4*p['M']*(linear*rho+2*cubic*rho**3)
    C1 = basis['upper'].T@((basis['weights']*W1)[:, None]*basis['lower'])
    Y0 = pure_start(kappa, basis['masses'])
    h = np.hypot(kappa, basis['masses'])
    right = np.r_[kappa, h]
    forcing = C1.T-Y0@C1@Y0-8*kappa*p['M']*Y0
    Y1 = forcing/(h[:, None]+right[None, :])
    residual = h[:, None]*Y1+Y1*right[None, :]-forcing
    return C1, Y1, float(np.max(abs(residual)))


def integrate_tail(start, basis, p, fields, kappa=.2, target=4., melvin=False, tighter=False, leading=False):
    """Z=-iY; integrating the real stable-graph Riccati equation backwards."""
    if not np.all(np.isfinite([start, target, kappa])) or not start > target > p['k'] or kappa <= 0:
        raise ValueError('Require finite start>target>k and kappa>0 for this positive axial lead.')
    if melvin and leading:
        raise ValueError('The curved leading correction does not apply to the exact Melvin control.')
    shape = (basis['pairs'], basis['pairs']+1)
    pure_initial = pure_start(kappa, basis['masses'])
    initial = pure_initial.copy()
    if leading:
        initial += leading_start_coefficient(kappa, basis, p)[1]/start

    def rhs(z, flat):
        Y = flat.reshape(shape)
        Bp, Bm, C = curved_coefficients(z, basis, p, fields, melvin=melvin)
        return (kappa*(Bm@Y+Y@Bp)-C.T+Y@C@Y).reshape(-1)

    solution = solve_ivp(rhs, (start, target), initial.reshape(-1), method='DOP853',
                         rtol=2e-12 if tighter else 2e-10,
                         atol=2e-13 if tighter else 2e-11, dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    records = []
    for z in (32., 16., 8., 4.):
        if z < target or z > start:
            continue
        Y = solution.sol(z).reshape(shape)
        Bp, Bm, C = curved_coefficients(z, basis, p, fields, melvin=melvin)
        records.append(dict(z=z, Y=Y.tolist(), Z_real=np.zeros_like(Y).tolist(), Z_imag=(-Y).tolist(),
            spectral_norm=float(np.linalg.norm(Y, 2)),
            departure_from_pure_Melvin=float(np.linalg.norm(Y-pure_initial, 2)),
            zero_column_norm=float(np.linalg.norm(Y[:, 0])),
            Bplus=Bp.tolist(), Bminus=Bm.tolist(), C=C.tolist(),
            Bplus_minimum_eigenvalue=float(np.min(np.linalg.eigvalsh(Bp))),
            Bminus_minimum_eigenvalue=float(np.min(np.linalg.eigvalsh(Bm)))))
    return dict(pairs=basis['pairs'], ritz_size=basis['ritz_size'], radial_order=basis['radial_order'],
                start_z=start, kappa=kappa, target_z=target, initial_Y=initial.tolist(),
                initialization=('Finite projected order1/L graph imposed at L; remaining far tail is not certified.' if leading else
                    'Pure Melvin asymptotic graph imposed at finite L; omitted L-to-infinity curved tail is not certified.'),
                leading_start=leading,
                pure_Melvin_control=melvin, tighter_ODE=tighter,
                ode_evaluations=solution.nfev, records=records)


def basis_summary(basis):
    p = basis['pairs']
    return dict(pairs=p, ritz_size=basis['ritz_size'], radial_order=basis['radial_order'],
        omega=basis['omega'], masses=basis['masses'].tolist(),
        squared_operator_ritz_residuals=basis['squared_operator_ritz_residuals'].tolist(),
        upper_Gram_error=float(np.linalg.norm(basis['gram_upper']-np.eye(p+1), 2)),
        lower_Gram_error=float(np.linalg.norm(basis['gram_lower']-np.eye(p), 2)),
        Melvin_coupling_error=float(np.linalg.norm(basis['C_Melvin_measured']-basis['C_Melvin_target'], 2)),
        lower_derivative_control=basis['lower_derivative_control'], zero_profile=basis['zero_profile'])


def run():
    inputs = ('balanced_dihole.py', 'balanced_dihole_results.json', 'melvin_leads.py')
    hashes = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in inputs}
    saved = json.loads(Path('balanced_dihole_results.json').read_text(encoding='utf-8'))
    if not saved['passed'] or saved['source_sha256'] != hashes['balanced_dihole.py']:
        raise ValueError('Successful matching dihole source/results are required.')
    p = saved['parameters']
    fields = symbolic_fields()[0]
    bases = [radial_basis(p, pairs=n) for n in (2, 3)]
    refined = radial_basis(p, pairs=3, size=48, order=144)
    summaries = [basis_summary(b) for b in bases+[refined]]
    if max(max(s['upper_Gram_error'], s['lower_Gram_error'], s['Melvin_coupling_error']) for s in summaries) > 2e-8:
        raise AssertionError('Radial basis normalization/pairing failed before axial integration.')
    records = [integrate_tail(L, b, p, fields) for b in bases for L in (32., 64., 128.)]
    refined_record = integrate_tail(128., refined, p, fields)
    tight_record = integrate_tail(128., refined, p, fields, tighter=True)
    pure = integrate_tail(32., refined, p, fields, melvin=True)
    leading_record = integrate_tail(32., refined, p, fields, leading=True)
    control = next(r for r in records if r['pairs'] == 3 and r['start_z'] == 128.)
    def difference(a, b, z):
        aa = next(r for r in a['records'] if r['z'] == z)
        bb = next(r for r in b['records'] if r['z'] == z)
        return float(np.linalg.norm(np.asarray(aa['Y'])-bb['Y'], 2))
    length_controls = []
    for n in (2, 3):
        selected = [r for r in records if r['pairs'] == n]
        for a, b in zip(selected[:-1], selected[1:]):
            length_controls.append(dict(pairs=n, starts=[a['start_z'], b['start_z']],
                difference_at4=difference(a,b,4.), difference_at8=difference(a,b,8.)))
    refinement = dict(radial_Ritz_and_quadrature_at4=difference(control,refined_record,4.),
        radial_Ritz_and_quadrature_at8=difference(control,refined_record,8.),
        tighter_ODE_at4=difference(refined_record,tight_record,4.),
        tighter_ODE_at8=difference(refined_record,tight_record,8.))
    C1, Y1, leading_residual = leading_start_coefficient(.2, refined, p)
    asymptotic_controls = []
    for z in (32., 64., 128.):
        Bp, Bm, C = curved_coefficients(z, refined, p, fields)
        asymptotic_controls.append(dict(z=z,
            Bplus_remainder_norm=float(np.linalg.norm(Bp-(1+4*p['M']/z)*np.eye(4), 2)),
            Bminus_remainder_norm=float(np.linalg.norm(Bm-(1+4*p['M']/z)*np.eye(3), 2)),
            C_remainder_norm=float(np.linalg.norm(C-refined['C_Melvin_target']-C1/z, 2))))
    leading_comparison = dict(
        first_order_Riccati_residual=leading_residual, C1=C1.tolist(), Y1=Y1.tolist(),
        leading_L32_vs_pure_L128_at4=difference(leading_record,refined_record,4.),
        leading_L32_vs_pure_L128_at8=difference(leading_record,refined_record,8.),
        leading_L32_vs_pure_L128_at32=difference(leading_record,refined_record,32.),
        projected_coefficient_remainders=asymptotic_controls,
        derivation='Bplus/minus=I+4M*I/z+O(z^-2); Wtrace=W+2M*(rho*W_prime+W)/z. C1=4M<plus,(a*rho+2c*rho^3)minus>; (h_i+h_right,j)*Y1_ij=(C1^T-Y0*C1*Y0-8*kappa*M*Y0)_ij, h_right=(kappa,h_j).',
        scope='The correction solves the asymptotic equation for the finite projected basis. Remainder samples and start comparison do not enclose the whole radial tail or omitted channels.')
    two = next(r for r in records if r['pairs'] == 2 and r['start_z'] == 128.)
    channel_sensitivity = []
    for z in (4., 8., 16., 32.):
        low = np.asarray(next(r for r in two['records'] if r['z'] == z)['Y'])
        larger = np.asarray(next(r for r in control['records'] if r['z'] == z)['Y'])
        channel_sensitivity.append(dict(z=z, pairs=[2,3],
            common_block_difference=float(np.linalg.norm(low-larger[:2,:3], 2)),
            interpretation='Measured finite-channel change, not a bound on omitted channels or a convergence certificate.'))
    reflection = []
    for z in (4., 8., 32.):
        a, b = curved_coefficients(z, refined, p, fields), curved_coefficients(-z, refined, p, fields)
        reflection.append(dict(z=z, maximum_projected_reflection_error=max(float(np.max(abs(v-w))) for v,w in zip(a,b))))
    mass_zero = dict(p, M=0., k=p['a'])
    actual_melvin_limit = curved_coefficients(4., refined, mass_zero, fields)
    analytic_melvin_limit = curved_coefficients(4., refined, p, fields, melvin=True)
    melvin_limit_error = max(float(np.max(abs(a-b))) for a,b in zip(actual_melvin_limit,analytic_melvin_limit))
    tests = dict(radial_normalization_and_pairing=max(max(s['upper_Gram_error'],s['lower_Gram_error'],s['Melvin_coupling_error']) for s in summaries)<2e-8,
        analytic_radial_derivative=max(s['lower_derivative_control'] for s in summaries)<2e-8,
        Ritz_squared_operator_residual=max(max(s['squared_operator_ritz_residuals']) for s in summaries)<2e-8,
        exact_zero_normalization=max(s['zero_profile']['quadrature_relative_error'] for s in summaries)<2e-9,
        pure_Melvin_constant_graph=max(r['departure_from_pure_Melvin'] for r in pure['records'])<2e-8,
        positive_frequency_weight=min(min(r['Bplus_minimum_eigenvalue'],r['Bminus_minimum_eigenvalue']) for r in refined_record['records'])>0,
        passive_finite_graph=max(r['spectral_norm'] for a in records+[refined_record,tight_record] for r in a['records'])<1+2e-8,
        reflected_projected_coefficients=max(r['maximum_projected_reflection_error'] for r in reflection)<2e-10,
        actual_Mzero_projected_Melvin_limit=melvin_limit_error<2e-10,
        leading_start_Riccati_identity=leading_residual<2e-13,
        radial_refinement=refinement['radial_Ritz_and_quadrature_at4']<2e-7,
        axial_ODE_refinement=refinement['tighter_ODE_at4']<2e-7,
        finite_start_stability=max(r['difference_at4'] for r in length_controls)<2e-7,
        historical_inputs_preserved=all(hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest for name,digest in hashes.items()))
    tests = {k:bool(v) for k,v in tests.items()}
    result = dict(passed=all(tests.values()), tests=tests, parameters=p,
        pilot=dict(mu=.5,kappa=.2, pairs=[2,3], starts=[32.,64.,128.], targets=[4.,8.,16.,32.]),
        bases=summaries, records=records, refined_record=refined_record, tight_ODE_record=tight_record,
        leading_start_record=leading_record, leading_start_control=leading_comparison,
        pure_Melvin_control=pure, finite_start_sensitivity=length_controls,
        refinement=refinement, finite_channel_sensitivity=channel_sensitivity,
        projected_reflection_controls=reflection,
        actual_Mzero_projected_Melvin_limit_error=melvin_limit_error,
        conventions=dict(
            trace='w=sqrt(v)*u with v=N/P. The axial current is integral w_dagger*sigma_z*w drho. Radial profiles use Lebesgue drho, with upper real and lower i times real.',
            operator='[-i*sigma_x*partial_rho-i*sigma_z*partial_z+sigma_y*N*P*(mu-Aphysical)/rho]w=(E/v)w.',
            projected_equation="A_prime=-kappa*Bplus*A-i*C*B; B_prime=kappa*Bminus*B+i*C_transpose*A.",
            Riccati="Z=B/A=-iY; Y_prime=kappa*(Bminus*Y+Y*Bplus)-C_transpose+Y*C*Y.",
            spaces='A has exact zero then n upper massive components; B has n lower massive components. Y is n by(n+1). No second artificial zero channel is inserted.',
            massive_mode_phases='Each lower Ritz vector has its largest-magnitude Laguerre coefficient positive; its upper partner is (partial_rho+W_Melvin)f_minus/lambda. This fixes phases before comparing different cutoffs.',
            infinity='Pure Melvin Z_paired=-i*lambda/(kappa+sqrt(kappa^2+lambda^2)), zero column0.',
            opposite_end='At the reflected negative-z end, the opposite graph is Zminus(z)=Zplus(-z)^T in these paired real phases. This Weyl has one physical incoming upper zero channel there and no outgoing zero channel. The rectangular zero row may mix with massive columns; it is not an extra channel.',
            gauge='Use actual project Aphysical=-q*Ahat/(2Q), regular outer-axis Aphysical=0; no removal of asymptotic flux without shifting mu.',
            units='D=1 Weyl coordinate separation, kappa and masses inverse D. z is Weyl axial coordinate, not proper distance.'),
        scope=[
            'A fixed finite-channel axial projection of the actual curved dihole at |z|>=4; not a full meridional PDE, mouth map, throat sewing or quantum-source calculation.',
            'Finite-L pure Melvin initial data omit the actual exterior beyond L. L refinement measures numerical sensitivity only; it is not a rigorous tail enclosure or proof of the infinite-domain map.',
            'One leading1/L corrected start includes the fixed-basis dilation mixing through the projected coefficients. Its remaining far-tail error is not certified.',
            'Radial Ritz residuals, quadrature, paired-basis and ODE controls are floating-point tests. Omitted massive channels and other azimuthal sectors remain uncontrolled.',
            'This helper caps its exp-weight Laguerre quadrature at144 nodes; larger orders require a different stable weight evaluation and are explicitly rejected.',
            'The chiral zero mode is exact in the Melvin limit. Curved finite-z projections can mix it with retained massive channels.',
            'The background retains its classical horizons and Melvin environment. No supporting quantum state, stability, traversability or time travel is established.'],
        input_sha256=hashes, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=saved.get('sources', []))
    Path(__file__).with_name('dihole_axial_tail_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'],tests=tests,bases=summaries,
        refined_profile=[{k:r[k] for k in ('z','spectral_norm','departure_from_pure_Melvin','zero_column_norm')} for r in refined_record['records']],
        finite_start_sensitivity=length_controls,finite_channel_sensitivity=channel_sensitivity,
        refinement=refinement),indent=2))
    if not result['passed']:
        raise AssertionError(tests)
    return result


if __name__ == '__main__':
    run()
