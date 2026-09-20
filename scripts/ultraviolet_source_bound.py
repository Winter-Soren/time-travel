"""Analytic all-angular ultraviolet source bound for the rho_c=32 core.

Every frequency integral is a finite positive polynomial/exponential formula.
The saved coefficient enclosures use mpmath interval arithmetic. No exterior
solve, angular truncation, fitted tail, or quadrature supplies the bound.
Canonical physical-flux passivity and the declared continuum core are assumed.
"""
from project_paths import project_path
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import time

import mpmath as mp


PROBE_FRACTIONS = (("central", Fraction(0), Fraction(65, 100)),
                   ("upper_side", Fraction(60, 100), Fraction(25, 100)))
MASS_GRID = tuple(map(Fraction, (0,))) + (
    Fraction(1, 8), Fraction(1, 4), Fraction(1, 2),
    Fraction(3, 4), Fraction(7, 8), Fraction(1))


def rational(ctx, value):
    value = Fraction(value)
    return ctx.mpf(value.numerator) / value.denominator


def add(*polys):
    result = [0] * max(map(len, polys))
    for poly in polys:
        for degree, coefficient in enumerate(poly):
            result[degree] += coefficient
    return result


def scale(poly, coefficient):
    return [coefficient * value for value in poly]


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def power(poly, exponent):
    result = [1]
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def moment_polynomial(ctx, order, pole, shift):
    """S_j(pole/(k+shift)) <= this polynomial in k, for n>=1."""
    result = [0]
    for degree in range(order + 1):
        coefficient = (math.comb(order, degree) * rational(ctx, Fraction(3, 2))
                       ** (order - degree) * math.factorial(degree)
                       / pole ** (degree + 1))
        result = add(result, scale(power([shift, 1], degree + 1), coefficient))
    return result


def exponential_moments(ctx, beta, start, degree):
    """Exact integral k^j exp(-beta*k), k=start..infinity; positive recurrence."""
    endpoint = ctx.exp(-beta * start)
    result = [endpoint / beta]
    for order in range(1, degree + 1):
        result.append((endpoint * start ** order + order * result[-1]) / beta)
    return result


def integrate_polynomial(ctx, coefficients, beta, start):
    moments = exponential_moments(ctx, beta, start, len(coefficients) - 1)
    return sum(a * b for a, b in zip(coefficients, moments))


def evaluate_polynomial(coefficients, value):
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def build_majorant(ctx, probe, bins=32, start=48, rho=32):
    if probe not in (0, 1) or not isinstance(bins, int) or bins < 1:
        raise ValueError("Use a known probe and positive integer force-bin count")
    if start <= 0 or rho <= 0:
        raise ValueError("Positive frequency threshold and throat cut required")
    start, rho = ctx.mpf(start), ctx.mpf(rho)
    cut = ctx.atan2(rho, ctx.mpf(1))
    length, ccut, lam = 2 * cut, ctx.sqrt(1 + rho * rho), ctx.sqrt(2)
    ell, shift = ccut / 5, 1 / (2 * start)

    def tanh_positive(value):
        small = ctx.exp(-2 * value)
        return (1 - small) / (1 + small)

    def primitive(fraction):
        if fraction == 0:
            return ctx.mpf(0)
        angle = rational(ctx, abs(fraction)) * cut
        tangent = ctx.sin(angle) / ctx.cos(angle)
        value = ctx.log(tangent + ctx.sqrt(1 + tangent * tangent))
        return value if fraction > 0 else -value

    t_length = tanh_positive(start * length)
    c_delta = ccut * (1 + shift / start) / (lam * t_length)
    barrier_pole = sum(
        lam / 2 * tanh_positive(start * rational(ctx, 2 * left) * cut)
        * (primitive(2 * right - 1) - primitive(2 * left - 1))
        for left, right in zip(MASS_GRID[:-1], MASS_GRID[1:]))
    transmission_s2 = moment_polynomial(ctx, 2, barrier_pole, shift)
    transmission_s3 = moment_polynomial(ctx, 3, barrier_pole, shift)
    factor = 8 * c_delta / (ctx.pi * ell)
    name, center, width = PROBE_FRACTIONS[probe]
    edges = [center - width + 2 * width * i / bins for i in range(bins + 1)]
    terms, force_weights = [], []
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        nearest = min(max(center, left), right)
        u = ((nearest - center) / width) ** 2
        peak = ctx.exp(1 - 1 / (1 - rational(ctx, u)))
        weight = rational(ctx, right - left) * cut * peak
        force_weights.append(weight)
        a, distance_right = rational(ctx, 1 + left) * cut, rational(ctx, 1 - right) * cut
        ta, tr = tanh_positive(start * a), tanh_positive(start * distance_right)
        max_sec = 1 / ctx.cos(rational(ctx, max(abs(left), abs(right))) * cut)
        al = lam * tr * (primitive(left) - primitive(Fraction(-1)))
        ar = lam * ta * (primitive(Fraction(1)) - primitive(right))
        for side, pole, distance in (("left", al, a), ("right", ar, distance_right)):
            s2 = moment_polynomial(ctx, 2, pole, shift)
            s3 = moment_polynomial(ctx, 3, pole, shift)
            poly = add(multiply([0, 4], s2), scale(s3, max_sec))
            terms.append(dict(kind="massive_reflection", bin=index, side=side,
                              beta=2 * distance, polynomial=scale(poly, factor * weight)))
        a0 = add([0, 2], scale(multiply([0, 1], [shift, 1]), 4 / (lam * ta)))
        a1 = scale(add([1], scale([shift, 1], 1 / (lam * ta))), max_sec)
        poly = scale(add(multiply(a0, transmission_s2),
                         multiply(a1, transmission_s3)), 2)
        terms.append(dict(kind="massive_transmission_pressure", bin=index,
                          beta=length, polynomial=scale(poly, factor * weight)))
        terms.append(dict(kind="massive_null_shift", bin=index, beta=length,
                          polynomial=scale(multiply([0, 4], transmission_s2),
                                           factor * weight)))
    force_integral_upper = sum(force_weights)
    zero_factor = 8 * force_integral_upper / (
        ctx.pi * ell * (1 - ctx.exp(-2 * start * length)))
    terms.append(dict(kind="zero_block", beta=length, polynomial=[0, zero_factor]))
    return dict(probe=name, start=start, rho_c=rho, ell=ell, length=length,
                shift=shift, delta_inverse_coefficient=c_delta,
                transmission_barrier_pole=barrier_pole,
                force_integral_upper=force_integral_upper, terms=terms)


def integrated_parts(ctx, data):
    parts = {name: ctx.mpf(0) for name in (
        "massive_reflection", "massive_transmission_pressure", "massive_null_shift", "zero_block")}
    for term in data["terms"]:
        parts[term["kind"]] += integrate_polynomial(
            ctx, term["polynomial"], term["beta"], data["start"])
    parts["total_pressure"] = (parts["massive_reflection"]
                               + parts["massive_transmission_pressure"] + parts["zero_block"])
    parts["total_null"] = parts["total_pressure"] + parts["massive_null_shift"]
    return parts


def interval_record(value):
    """One extra binary64 ulp encloses conversion of each interval endpoint."""
    lo, hi = float(value.a), float(value.b)
    return dict(lower=math.nextafter(lo, -math.inf), upper=math.nextafter(hi, math.inf))


def controls():
    from radial_angular_envelope import coefficients, envelope, y_lower
    from throat_null_response import PROBES

    tests = {"probe_definitions_match": all(
        old["name"] == name and Fraction(str(old["center_fraction"])) == center
        and Fraction(str(old["width_fraction"])) == width
        for old, (name, center, width) in zip(PROBES, PROBE_FRACTIONS))}
    moment_ratios, integral_errors, comparisons = [], [], []
    with mp.workdps(80):
        for alpha in map(mp.mpf, (".001", ".03", ".2", "1", "5", "20")):
            q, gap = mp.exp(-alpha), -mp.expm1(-alpha)
            exact = (q * (mp.mpf("2.25") / gap + 3*q/gap**2 + q*(1+q)/gap**3),
                     q * (mp.mpf("3.375") / gap + mp.mpf("6.75")*q/gap**2
                          + mp.mpf("4.5")*q*(1+q)/gap**3 + q*(1+4*q+q*q)/gap**4))
            for order, actual in zip((2, 3), exact):
                bound = sum(math.comb(order, j)*mp.mpf("1.5")**(order-j)
                            * math.factorial(j)/alpha**(j+1) for j in range(order+1))
                moment_ratios.append(float(actual/bound))
        for beta in map(mp.mpf, (".4", "1", "3.1")):
            values = exponential_moments(mp, beta, mp.mpf(48), 5)
            for degree, actual in enumerate(values):
                independent = mp.gammainc(degree+1, 48*beta, mp.inf)/beta**(degree+1)
                integral_errors.append(float(abs(actual/independent-1)))
        for probe in (0, 1):
            majorant = build_majorant(mp, probe)
            for frequency in (48., 64., 96., 192., 1024.):
                old = coefficients(frequency, 32., probe, cutoff=0)
                sharper_cdelta = math.hypot(1., 32.) / (
                    2*frequency*y_lower(math.sqrt(2), frequency, old["length"]))
                radial = envelope(old)["scaled_null_bound"] * (
                    sharper_cdelta/old["delta_inverse_coefficient"])
                ours = sum(evaluate_polynomial(t["polynomial"], mp.mpf(frequency))
                           * mp.exp(-t["beta"]*frequency) for t in majorant["terms"]
                           if t["kind"] != "zero_block")
                comparisons.append(dict(probe=probe, x=frequency,
                    log10_analytic_majorant=float(mp.log10(ours)),
                    radial_to_majorant_ratio=float(radial/ours)))
    tests.update(moment_cell_majorant=max(moment_ratios) <= 1,
                 positive_gamma_recurrence=max(integral_errors) < 1e-70,
                 pointwise_radial_formula_comparison=max(r["radial_to_majorant_ratio"] for r in comparisons) <= 1+1e-11)
    return tests, dict(moment_actual_over_bound=moment_ratios,
                       maximum_gamma_relative_error=max(integral_errors), pointwise_controls=comparisons,
                       role="Independent implementation checks; sampled controls do not establish the analytic inequalities or integral bound.")


def run():
    start_time = time.perf_counter()
    root = Path(__file__).resolve().parent
    inputs = ("radial_angular_envelope.py", "throat_null_response.py", "long_throat_integral_results.json")
    hashes = {name: hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs}
    frozen = json.loads((project_path("long_throat_integral_results.json")).read_text(encoding="utf-8"))
    tests, details = controls()
    tests["frozen_integral_passed"] = frozen["passed"]
    signals = frozen["integration"]["cutoffs"]["4"]["partial_sums"]["fine"]["totals"]["pi"]["pressure"]
    records, enclosure_checks = [], []
    for probe in (0, 1):
        mp.iv.dps = 70
        data = build_majorant(mp.iv, probe)
        coarse = integrated_parts(mp.iv, data)
        mp.iv.dps = 100
        fine_data = build_majorant(mp.iv, probe)
        fine = integrated_parts(mp.iv, fine_data)
        parts = {name: interval_record(value) for name, value in fine.items()}
        with mp.workdps(120):
            point = integrated_parts(mp, build_majorant(mp, probe))
            for name, value in point.items():
                lo, hi = (mp.make_mpf(endpoint) for endpoint in fine[name]._mpi_)
                enclosure_checks.append(bool(lo <= value <= hi))
        for name, value in fine.items():
            enclosure_checks.append(bool(coarse[name].a <= value.a and value.b <= coarse[name].b))
        records.append(dict(probe=PROBE_FRACTIONS[probe][0], integrals=parts,
            primary_finite_window_pressure=signals[probe],
            null_bound_over_absolute_primary_signal=parts["total_null"]["upper"]/abs(signals[probe]),
            constants={name: interval_record(fine_data[name]) for name in (
                "ell", "length", "shift", "delta_inverse_coefficient", "transmission_barrier_pole", "force_integral_upper")},
            exponential_polynomials=[dict(kind=t["kind"], bin=t.get("bin"), side=t.get("side"),
                beta=interval_record(t["beta"]), coefficients=[interval_record(mp.iv.mpf(c)) for c in t["polynomial"]])
                for t in fine_data["terms"]]))
    tests["interval_precision_and_independent_point_enclosure"] = all(enclosure_checks)
    tests["positive_integral_enclosures"] = all(
        value["lower"] > 0 for row in records for value in row["integrals"].values())
    tests["uv_bound_below_primary_finite_window_signal"] = all(
        row["null_bound_over_absolute_primary_signal"] < 1 for row in records)
    tests["frozen_inputs_preserved"] = hashes == {
        name: hashlib.sha256((project_path(name)).read_bytes()).hexdigest() for name in inputs}
    assert all(tests.values()), tests
    result = dict(passed=True, tests=tests, seconds=time.perf_counter()-start_time,
        parameters=dict(q=1, b=.2, rho_c=32, x_min=48, force_bins=32,
                        mass_grid=[str(x) for x in MASS_GRID], angular_range="zero plus every massive n>=1",
                        interval_decimal_precision=100, integration_measure="dx/ell=dκ"),
        records=records, controls=details, input_sha256=hashes,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        equations=dict(
            comparison="λ*=sqrt2, h=1/(2x0), T(s)=tanh(x0*s): Y(s)>=λ*T(s)/(2(x+h)).",
            defect="g=δ-(2x/M)y obeys g'=-2(my+x)g+2x(1+y²)(1-m/M); m<=M=λ_n sec(c), hence δ>=2xy/M.",
            uniform_defect="δ_n^-1<=cUV λ_n; cUV=sec(c)*(1+h/x0)/(sqrt2*T(L)).",
            moment="S_j(α)<=sum_(l=0)^j binom(j,l)*(3/2)^(j-l)*l!/α^(l+1), from positive unit-cell integrals.",
            gamma="I0=exp(-βx0)/β; Ij=[exp(-βx0)*x0^j+j*I_(j-1)]/β.",
            zero="Integral <=[8 F/(πell(1-exp(-2x0L)))]*exp(-Lx0)*(x0/L+1/L²)."),
        scope=[
            "Conditional absolute upper bound on the x>=48 integral of each full-Dirac relative pressure/null source; both Weyls and every angular mode counted once.",
            "Assumes the declared product core, nonnegative compact probes, canonical physical-flux contraction of connected/reference exterior maps, and their well-posed continuum limit.",
            "All angular blocks are included, so no omitted-channel feedback term is excluded from this UV bound.",
            "mpmath interval evaluation encloses the analytic majorant coefficients and positive integrals; saved binary64 endpoints are expanded outward by one ulp.",
            "This does not bound the low-frequency threshold, finite-frequency angular feedback, discretization errors, or supporting geometry; it does not certify the full source sign."])
    target = project_path("ultraviolet_source_bound_results.json")
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(dict(passed=True, tests=tests, seconds=result["seconds"],
                         records=[{k:v for k,v in row.items() if k != "exponential_polynomials"} for row in records]), indent=2))
    return result


if __name__ == "__main__":
    run()
