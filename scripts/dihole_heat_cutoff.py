"""Exact heat-cutoff bookkeeping in checkpoint136's fixed physical scheme.

Only frozen local coefficients are evaluated. The actual ultraviolet remainder
and infrared heat trace are unknown here; this module does not produce a total
source or a physical error bound. Its independent kernel is an algebra control.
"""
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

from mpmath import mp

from dihole_pv_counterterms import MU_H_D0, coefficients
from project_paths import project_path, source_matches


NAME = 'dihole_heat_cutoff'
FROZEN_INPUTS = {
    'dihole_local_source_results.json': 'cc219789568fda40584506837197cb7aae5a564d60b71356e365d026faab69a6',
    'dihole_pv_counterterms_results.json': 'c62c5a7853592e7cd91ac990fd857fee6b0e1c269d1a6bd2aa168c64401e1e8e',
}
EULER_GAMMA = float(mp.euler)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_number(value, name, *, positive=False, nonnegative=False):
    if value is None or isinstance(value, bool):
        raise ValueError(name+' must be a supplied finite real number.')
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(name+' is outside its finite permitted range.')
    return value


def local_cutoff_term(s0, data, mu_D0=MU_H_D0):
    """Return L(s0;mu), with s0 measured in D0^2 and response in hbar*c/D0.

    data has the three keys returned by dihole_pv_counterterms.coefficients.
    This evaluates stored numerical coefficients, not coefficient enclosures.
    The optional nonzero a0 is used only in the independent control below.
    """
    s0 = finite_number(s0, 's0', positive=True)
    mu = finite_number(mu_D0, 'mu_D0', positive=True)
    a0, a2, a4 = [finite_number(data[key], key) for key in
                  ('a0_variation', 'a2_action_derivative', 'a4_action_derivative')]
    factor = 1/(32*math.pi**2)
    pieces = dict(a0_term=-factor*a0/(2*s0)/s0, a2_term=-factor*a2/s0,
                  a4_term=factor*a4*(math.log(s0)+2*math.log(mu)+EULER_GAMMA))
    total = math.fsum(pieces.values())
    if not all(math.isfinite(x) for x in (*pieces.values(), total)):
        raise ValueError('The local cutoff expression exceeds binary64 range.')
    return dict(s0_over_D0_squared=s0, mu_times_D0=mu, **pieces, local_term=total)


def conditional_absolute_bound(s0, data, *, uv_integral_bound, ir_integral_bound,
                               mu_D0=MU_H_D0):
    """Triangle expression only, requiring BOTH supplied unhalved integral bounds.

    UV bounds integral_0^s0 |deltaK-H| ds/s; IR bounds integral_s0^infinity
    |deltaK| ds/s. This binary64 evaluation does not certify supplied bounds,
    local coefficients, or rounding. An Abel value alone is not an IR bound.
    """
    uv = finite_number(uv_integral_bound, 'UV integral bound', nonnegative=True)
    ir = finite_number(ir_integral_bound, 'absolute IR integral bound', nonnegative=True)
    value = abs(local_cutoff_term(s0, data, mu_D0)['local_term'])+uv/2+ir/2
    return finite_number(value, 'conditional triangle expression', nonnegative=True)


def authenticated_coefficients():
    inputs, evidence, records = {}, {}, {}
    for name, digest in FROZEN_INPUTS.items():
        path = project_path(name)
        if sha(path) != digest:
            raise ArithmeticError('Frozen local input mismatch: '+name)
        record = json.loads(path.read_text(encoding='utf-8'))
        if record.get('passed') is not True:
            raise ArithmeticError('Frozen local calculation did not pass: '+name)
        inputs[name] = digest
        source = name.replace('_results.json', '.py')
        dependencies = {source: record['source_sha256'], **record.get('input_sha256', {})}
        for old_path, expected in dependencies.items():
            filename = Path(old_path).name
            actual = sha(project_path(filename))
            if not source_matches(filename, expected, actual):
                raise ArithmeticError('Historical dependency mismatch: '+filename)
            inputs[filename] = actual
            evidence[filename] = dict(recorded_sha256=expected, current_sha256=actual,
                accepted_as='exact' if expected == actual else 'recorded layout migration')
        records[name] = record
    data = coefficients(records['dihole_local_source_results.json'])
    saved = records['dihole_pv_counterterms_results.json']
    if data != saved['coefficients'] or saved['mu_H_times_D0'] != MU_H_D0:
        raise ArithmeticError('Frozen counterterms and regenerated local coefficients differ.')
    return data, inputs, evidence


def kernel_controls():
    """A decaying auxiliary kernel, checked independently by quadrature/Mellin.

    deltaK=exp(-s)*(2/s^2+3/s+4)/(16*pi^2) has (a0,a2,a4)=(2,1,2).
    Its Mellin transform is [2 Gamma(z-2)+3 Gamma(z-1)+4 Gamma(z)]/(16*pi^2).
    Finite part plus this scheme's gamma_E term gives the exact target below.
    """
    with mp.workdps(70):
        mu, factor = mp.mpf(7)/5, 1/(32*mp.pi**2)
        expected = factor*(-mp.mpf(3)/2+2*mp.log(mu**2))
        # The cancellation-safe Taylor expansion is used only close to zero.
        series = [2*(-1)**n/mp.factorial(n)+3*(-1)**(n-1)/mp.factorial(n-1)
                  +4*(-1)**(n-2)/mp.factorial(n-2) for n in range(3, 64)]

        def heat(s):
            return mp.exp(-s)*(2/s**2+3/s+4)/(16*mp.pi**2)

        def subtracted_over_s(s):
            if abs(s) < mp.mpf('.01'):
                return mp.polyval(series[::-1], s)/(16*mp.pi**2)
            return (mp.exp(-s)*(2+3*s+4*s*s)-(2+s+2*s*s))/(16*mp.pi**2*s**3)

        def local(s):
            return factor*(-1/s**2-1/s+2*(mp.log(mu**2*s)+mp.euler))

        cut_rows, public_errors = [], []
        test_data = dict(a0_variation=2., a2_action_derivative=1., a4_action_derivative=2.)
        for s0 in map(mp.mpf, ('.03', '.3', '3')):
            uv = mp.quad(subtracted_over_s, [0, mp.mpf('.01'), s0])/2
            ir = mp.quad(lambda s: heat(s)/s, [s0, max(1, s0), mp.inf])/2
            value = uv+ir+local(s0)
            cut_rows.append(dict(s0=float(s0), value=mp.nstr(value, 62),
                                 absolute_error=float(abs(value-expected))))
            public = local_cutoff_term(float(s0), test_data, float(mu))['local_term']
            public_errors.append(float(abs(mp.mpf(public)-local(s0))/max(1, abs(local(s0)))))
        pv_rows = []
        beta0, beta1 = (9*mp.log(3)-12*mp.log(2))/2, 6*mp.log(2)-3*mp.log(3)
        for mass in map(mp.mpf, ('.7', '2')):
            def pv_integrand(s):
                if not s:
                    return 2*mass**6/(16*mp.pi**2)
                return (-mp.expm1(-mass**2*s))**3*heat(s)/s
            def correction(s):
                e = mp.exp(-mass**2*s)
                return (3*e-3*e**2+e**3)*subtracted_over_s(s)
            anchor = mp.quad(pv_integrand, [0, mp.mpf('.01'), 1, mp.inf])/2
            remainder = mp.quad(correction, [0, mp.mpf('.01'), 1, mp.inf])/2
            ct = factor*(2*beta0*mass**4+beta1*mass**2
                         +2*(mp.log(mass**2/mu**2)+mp.log(mp.mpf(3)/8)))
            recovered = anchor-ct+remainder
            pv_rows.append(dict(mass=float(mass), recovered=mp.nstr(recovered, 62),
                                absolute_error=float(abs(recovered-expected))))
        # Ordinary truncations of integral sin(s) ds alternate between 0 and 2;
        # its Abel value is lim_eta 1/(1+eta^2)=1. This is NOT a geometry kernel.
        oscillatory_error = max(abs(mp.quad(mp.sin, [0, 2*mp.pi])),
                                abs(mp.quad(mp.sin, [0, mp.pi])-2))
        return dict(kernel='exp(-s)*(2/s^2+3/s+4)/(16*pi^2), auxiliary only',
            coefficients=[2, 1, 2], mu=mp.nstr(mu, 62), expected=mp.nstr(expected, 62),
            expected_formula='[-3/2+2 log(mu^2)]/(32 pi^2)', cutoff_rows=cut_rows,
            finite_PV_rows=pv_rows, public_local_max_scaled_error=max(public_errors),
            IR_counterexample=dict(deltaK='s*sin(s)', bound='|deltaK(s)| <= s',
                short_time='O(s^2), hence valid O(s) remainder with H=0',
                ordinary_truncations='Integral_0^T deltaK(s) ds/s = 1-cos(T); no limit.',
                abel_integral='Integral_0^infinity exp(-eta*s)*sin(s) ds = 1/(1+eta^2) -> 1.',
                absolute_integral='Integral_0^infinity |sin(s)| ds diverges.',
                finite_quadrature_control_error=float(oscillatory_error)))


def run():
    start, own_hash = perf_counter(), sha(__file__)
    data, inputs, provenance = authenticated_coefficients()
    controls = kernel_controls()
    rows = [local_cutoff_term(s, data) for s in (1e-4, 1e-3, .01, .1, 1.)]
    rejected = []
    for uv, ir in ((None, 1.), (1., None), (math.inf, 1.), (1., -1.)):
        try:
            conditional_absolute_bound(1., data, uv_integral_bound=uv, ir_integral_bound=ir)
        except ValueError:
            rejected.append(True)
        else:
            rejected.append(False)
    scale_difference = (local_cutoff_term(1., data, 2.)['local_term']
                        -local_cutoff_term(1., data, 1.)['local_term'])
    scale_error = abs(scale_difference-data['a4_action_derivative']*math.log(4)/(32*math.pi**2))
    tests = dict(authenticated_inputs_unchanged=all(sha(project_path(n)) == h for n, h in inputs.items()),
        cutoff_invariance_with_gamma_E=max(r['absolute_error'] for r in controls['cutoff_rows']) < 1e-50,
        independent_finite_PV_identity=max(r['absolute_error'] for r in controls['finite_PV_rows']) < 1e-50,
        numerical_local_formula=controls['public_local_max_scaled_error'] < 2e-15,
        fixed_scheme_scale_sign=scale_error < 1e-16,
        unknown_or_invalid_bounds_rejected=all(rejected),
        Abel_is_not_ordinary_IR_control=controls['IR_counterexample']['finite_quadrature_control_error'] < 1e-50)
    result = dict(complete=True, passed=all(tests.values()), tests=tests, coefficients=data,
        mu_H_times_D0=MU_H_D0, local_rows=rows, controls=controls, scale_identity_error=scale_error,
        equations=dict(heat='deltaK(s)=-s Tr_per_Killing_time(dotP exp(-sP)); P is the common physical4D operator.',
            subtraction='H(s)=(a0/s^2+a2/s+a4)/(16 pi^2), aj=delta integrated physical4D Aj including volume variation.',
            combined='Q_HK=-C_M0+(1/2) Abel Integral_0^infinity [deltaK-(1-w_M0)H] ds/s, w_M=(1-exp(-M^2 s))^3.',
            cutoff='Q_HK=(1/2)Integral_0^s0(deltaK-H)ds/s+(1/2)Abel Integral_s0^infinity deltaK ds/s+L(s0;mu_H).',
            local='L=[-a0/(2 s0^2)-a2/s0+a4*(log(mu_H^2 s0)+gamma_E)]/(32 pi^2).',
            cancellation='dL/ds0=H(s0)/(2 s0), cancelling the derivative of the two split integrals.',
            PV_log_moment='Integral_0^s0 w_M ds/s=log(M^2 s0)+gamma_E+log(3/8)+3E1(M^2s0)-3E1(2M^2s0)+E1(3M^2s0).',
            conditional_bound='If B_UV >= Integral_0^s0 |deltaK-H|ds/s and B_IR >= Integral_s0^infinity |deltaK|ds/s, then |Q_HK| <= |L|+(B_UV+B_IR)/2.'),
        actual_unknowns=dict(UV_remainder_integral=None, UV_absolute_integral_bound=None,
            IR_heat_integral=None, IR_absolute_integral_bound=None, total_response=None,
            total_absolute_response_bound=None),
        units='s0 is in D0^2, mu is in 1/D0; numerical responses are in hbar*c/D0 with D0=1.',
        prescription='Checkpoint136 physical heat/zeta scheme: mu_H=sqrt(2)*exp(-gamma_E)/D0, fixed couplings, no added finite R2/Ric2/F2 actions; same state and ordered IR prescription.',
        scope=['The physical anchor cancels algebraically before estimates; neither it nor a global C_R is evaluated here.',
            'The O(s) ultraviolet remainder makes the subtracted short-time integral locally integrable; its magnitude is not bounded by the three local coefficients alone.',
            'Abel means eta down to zero after multiplying the heat integrand by exp(-eta*s), using the established ordered physical prescription. It is not an absolute IR bound.',
            'The bound |deltaK|<=C*s by itself permits an Abel limit without an ordinary improper integral, as the auxiliary oscillatory control shows.',
            'Real local coefficient inputs are saved numerical quadratures, not interval enclosures. Conditional bound arithmetic also needs rounding and coefficient error control for a certificate.',
            'No actual heat kernel, geometry, angular/frequency completion, or physical source magnitude/sign is calculated; a compressed finite model cannot replace this physical4D kernel.'],
        input_sha256=inputs, historical_source_authentication=provenance,
        source_sha256=own_hash, elapsed_seconds=perf_counter()-start)
    project_path(NAME+'_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(passed=result['passed'], tests=tests, local_rows=rows,
        elapsed_seconds=result['elapsed_seconds'], source_sha256=own_hash,
        results_sha256=sha(project_path(NAME+'_results.json'))), indent=2))
    if not result['passed']:
        raise ArithmeticError('Cutoff controls failed; saved evidence retained.')
    return result


if __name__ == '__main__':
    run()
