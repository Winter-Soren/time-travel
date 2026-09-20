"""Local PV subtraction for the explicit physical4D heat/zeta reference.

Reuse133's saved compact action variations; do not rerun geometry or fit a
coupling. These terms are the large-proper-mass subtraction, not a vacuum
source. They must not be subtracted from a fixed angular/absorbing pilot and
called a renormalized physical answer.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from mpmath import mp

from project_paths import project_path


MU_H_D0 = math.sqrt(2)*math.exp(-np.euler_gamma)
BETA1 = 6*math.log(2)-3*math.log(3)


def coefficients(local):
    row = next(r for r in local['records'] if r['order'] == 128)
    action = row['action_derivative']
    charged_f2_action = -row['charge_included_F2_energy_derivative']
    # Static Euclidean invariant integrals equal the Lorentzian action
    # integrals here; their energy/action sign conventions differ.
    a2 = -action[0]/3
    a4 = action[1]/30-action[2]/10+2*charged_f2_action/3
    factor = 1/(32*math.pi**2)
    return dict(a0_variation=0., a2_action_derivative=a2, a4_action_derivative=a4,
                quadratic_mass_response_coefficient=BETA1*a2*factor,
                logarithmic_response_coefficient=a4*factor,
                scale_only_response=(math.log(2)-2*np.euler_gamma)*a4*factor)


def counterterm_response(mass_D0, data, mu_D0=MU_H_D0):
    if not math.isfinite(mass_D0) or mass_D0 <= 0 or not math.isfinite(mu_D0) or mu_D0 <= 0:
        raise ValueError('Positive finite proper mass and fixed renormalization scale required.')
    bracket = 2*(math.log(mass_D0)-math.log(mu_D0))+math.log(3/8)
    quadratic = data['quadratic_mass_response_coefficient']*mass_D0**2
    logarithmic = data['logarithmic_response_coefficient']*bracket
    return dict(mass_D0=mass_D0, quadratic=quadratic, logarithmic=logarithmic,
                response_to_subtract=quadratic+logarithmic)


def flat_scale_control():
    """Independent auxiliary Green-kernel limit, not a charged stress proof."""
    with mp.workdps(70):
        mass, radius, ell = mp.mpf('0.7'), mp.mpf('1e-12'), mp.mpf('1.3')
        sigma = radius**2/2
        exact = mass*mp.besselk(1, mass*radius)/(4*mp.pi**2*radius)
        parametrix = 1/(8*mp.pi**2*sigma)+mass**2/(16*mp.pi**2)*mp.log(sigma/ell**2)
        mu = mp.sqrt(2)*mp.exp(-mp.euler)/ell
        expected = mass**2/(16*mp.pi**2)*(mp.log(mass**2/mu**2)-1)
        return dict(mass=float(mass), radius=float(radius), ell=float(ell),
                    hadamard_finite_part=float(exact-parametrix), zeta_finite_part=float(expected),
                    finite_radius_error=float(abs(exact-parametrix-expected)))


def run():
    path = project_path('dihole_local_source_results.json')
    local = json.loads(path.read_text())
    if not local['passed']: raise RuntimeError('The frozen local-action source did not pass its checks.')
    data = coefficients(local)
    flat = flat_scale_control()
    weights = (1, -3, 3, -1); a = Fraction(3, 2)
    rational = all(sum(Fraction(c)/(x+j*a) for j, c in enumerate(weights))
                   == 6*a**3/(x*(x+a)*(x+2*a)*(x+3*a))
                   for x in (Fraction(1, 100), Fraction(1), Fraction(10000)))
    rows = [counterterm_response(float(m), data) for m in (1, 2, 4, 8, 16, 32, 64, 128)]
    scale_error = max(abs(counterterm_response(m, data)['response_to_subtract']
                         -counterterm_response(m, data, 1.)['response_to_subtract']
                         +data['scale_only_response']) for m in (1., 4., 16.))
    # Euler and total derivatives drop out only after the compact bulk
    # variation. No actual physical determinant/IR limit is evaluated here.
    tests = dict(auxiliary_scale_matching=flat['finite_radius_error'] < 1e-22,
                 rational_PV_identity=rational, scale_response_sign=scale_error < 3e-16,
                 finite_counterterms=all(math.isfinite(r['response_to_subtract']) for r in rows))
    result = dict(complete=True, passed=all(tests.values()), tests=tests,
                  mu_H_times_D0=MU_H_D0, beta1=BETA1, coefficients=data,
                  flat_control=flat, scale_response_identity_error=scale_error, counterterm_rows=rows,
                  prescription='Q_HK=lim_M [Q_PV(M)-Q_local(M;mu_H)] after physical fixed-M limits; no extra finite R2/Ric2/F2 terms.',
                  units='D0=1 in numerical values; all responses are in hbar*c/D0.',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest()},
                  scope=['Paired variation preserves N*sqrt(h) exactly; a0 contribution vanishes.',
                         'Physical4D a2=-R/3 and compact variation of a4=R2/30-Ric2/10+(2/3)F2 reused from133.',
                         'Large-M local asymptotic subtraction, not an approximation of the full response at MD0=1.',
                         'No finite-M pilot result is modified. Gravitational/gauge couplings remain fixed in the declared scheme.',
                         'Chosen action reference is explicit; equivalence to133 charged point-split reference is not claimed.',
                         'No physical-exterior, angular, frequency, infrared or regulator-removal limit computed.'])
    project_path('dihole_pv_counterterms_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
    if not result['passed']: raise RuntimeError('Counterterm normalization controls failed; evidence retained.')
    return result


if __name__ == '__main__': run()
