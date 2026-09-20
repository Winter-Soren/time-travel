"""A specified unitary loop/reservoir vacuum benchmark, not the Ori/MMP exterior.

Units are hbar=c=1.  A complex chiral fermion propagates around a loop of
length L with NS phase -1.  Its vertex is S(E)=exp(i E tau) S0, where
S0=[[-r*, sqrt(1-|r|^2)], [sqrt(1-|r|^2), r]] and r=R exp(i beta).
The nonnegative tau is an explicit causal vertex delay.  All incoming
reservoir levels with E<0 are occupied, with chemical potential zero.
A Dirac channel here consists of two chiral copies with phases +/-beta.

The reference retains the same one-pass delay, but has R=0.  Only the
return-dependent interaction energy is calculated.  No real exterior
amplitude is inferred from a probability, and no general analytic-
continuation API or physical state-preparation claim is supplied.

Primary foundations (the prefactors below are independently checked):
  Harrison, Weyand, Kirsten, J. Math. Phys. 57, 102301 (2016), Eqs. 18-21:
  https://arxiv.org/pdf/1606.07834
  Wirzba, J. Phys. A 41, 164003 (2008), Eqs. 8-13, relative spectral density:
  https://arxiv.org/pdf/0711.2395

Run only this benchmark: .venv\\Scripts\\python.exe -B -W error scattering_vacuum.py
"""

from __future__ import annotations
from project_paths import project_path

import cmath
from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
import time

import mpmath as mp


def _finite_real(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _phase(angle: float) -> complex:
    if not math.isfinite(angle):
        raise ValueError("combined propagation phase must be finite")
    return cmath.exp(1j * angle)


@dataclass(frozen=True)
class LoopReservoir:
    """This particular scale-independent vertex and its causal delay only."""

    modulus: float
    length: float
    phase: float = 0.0
    delay: float = 0.0

    def __post_init__(self) -> None:
        for name in ("modulus", "length", "phase", "delay"):
            object.__setattr__(self, name, _finite_real(name, getattr(self, name)))
        if not 0.0 <= self.modulus <= 1.0:
            raise ValueError("modulus must lie in [0, 1]")
        if self.length <= 0.0 or self.delay < 0.0:
            raise ValueError("length must be positive and causal delay nonnegative")
        if not math.isfinite(self.length + 2.0 * self.delay):
            raise ValueError("combined reference delay must remain finite")

    @property
    def round_trip_length(self) -> float:
        return self.length + self.delay

    def scattering(self, energy: float) -> dict:
        """Real-energy, flux-normalized amplitudes; includes the reservoir."""
        energy = _finite_real("energy", energy)
        r = self.modulus * _phase(self.phase)
        t = math.sqrt((1.0 - self.modulus) * (1.0 + self.modulus))
        delay_phase = _phase(energy * self.delay)
        vertex = ((-delay_phase * r.conjugate(), delay_phase * t),
                  (delay_phase * t, delay_phase * r))
        z = -_phase(energy * self.length)
        determinant = 1.0 - vertex[1][1] * z
        # At R=1 the reservoir decouples. Its scattering alone omits the
        # isolated loop eigenstates, which the vacuum determinant retains.
        effective = vertex[0][0]
        if self.modulus < 1.0:
            effective += vertex[0][1] * z * vertex[1][0] / determinant
        reference = delay_phase**2 * z
        return {"vertex": vertex, "loop_phase": z,
                "determinant": determinant, "effective": effective,
                "reference": reference}

    def densities(self, energy: float) -> dict[str, float]:
        """Ordinary continuum densities exist here only for R<1."""
        energy = _finite_real("energy", energy)
        if self.modulus == 1.0:
            raise ValueError("R=1 has a discrete loop spectrum, not this density")
        amplitude = self.modulus * _phase(
            self.phase + energy * self.round_trip_length)
        denominator = abs(1.0 + amplitude)**2
        poisson = (1.0 - self.modulus**2) / denominator
        interaction = -self.round_trip_length / math.pi * (
            amplitude / (1.0 + amplitude)).real
        loop_reference = self.length / (2.0 * math.pi)
        reference = (self.length + 2.0 * self.delay) / (2.0 * math.pi)
        return {"interaction": interaction,
                "loop": loop_reference * poisson,
                "loop_reference": loop_reference,
                "reservoir_total": reference + interaction,
                "reservoir_reference": reference}

    def vacuum_energy(self, kind: str = "dirac", quadrature: bool = True) -> mp.mpf:
        """Relative filled-sea energy, from this model's explicit determinant.

        D(i*kappa)=1+R*exp(i*beta-kappa*(L+tau)).  The second Dirac
        chirality has conjugate beta and the same real energy. This is
        not a routine for continuing arbitrary sampled scattering data.
        """
        if kind not in ("complex_chiral", "dirac"):
            raise ValueError("kind must be 'complex_chiral' or 'dirac'")
        if not isinstance(quadrature, bool):
            raise TypeError("quadrature must be a boolean")
        copies = 1 if kind == "complex_chiral" else 2
        with mp.workdps(60):
            r = mp.mpf(self.modulus) * mp.exp(1j * mp.mpf(self.phase))
            total = mp.mpf(self.length) + mp.mpf(self.delay)
            if quadrature:
                # y=kappa*(L+tau). A return-dependent contour integral;
                # a general scattering system requires a separate justification.
                integral = mp.quad(
                    lambda y: mp.re(mp.log1p(r * mp.exp(-y))),
                    [0, 1, mp.inf])
                return -copies * integral / (mp.pi * total)
            return copies * mp.re(mp.polylog(2, -r)) / (mp.pi * total)

    def chiral_stress(self) -> mp.mpf:
        """Uniform physical-loop T_uu (or conjugate T_vv) in this toy.

        The global energy also includes the vertex delay: when tau>0,
        dividing that energy by L does not give the loop density.
        """
        with mp.workdps(60):
            r = mp.mpf(self.modulus) * mp.exp(1j * mp.mpf(self.phase))
            total = mp.mpf(self.length) + mp.mpf(self.delay)
            return mp.re(mp.polylog(2, -r)) / (mp.pi * total**2)


def _mp_scattering_ratio(model: LoopReservoir, energy: mp.mpf) -> mp.mpc:
    """Independent vertex elimination for the spectral-phase derivative check."""
    r = mp.mpf(model.modulus) * mp.exp(1j * mp.mpf(model.phase))
    t = mp.sqrt(1 - abs(r)**2)
    a = mp.exp(1j * energy * mp.mpf(model.delay))
    z = -mp.exp(1j * energy * mp.mpf(model.length))
    s00, s01, s10, s11 = -a * mp.conj(r), a*t, a*t, a*r
    effective = s00 + s01 * z * s10 / (1 - s11*z)
    return effective / (a*a*z)


def run_checks() -> dict:
    start = time.perf_counter()
    maxima = {key: 0.0 for key in (
        "vertex_unitarity", "reservoir_unitarity", "relative_determinant",
        "spectral_phase_derivative", "local_mode_density", "delay_density_relation",
        "winding_series", "imaginary_integral_vs_dilogarithm", "zeta_normalization",
        "lower_half_plane_pole_residual")}
    models = [LoopReservoir(r, 1.3, phase, delay)
              for r, phase, delay in ((.2, 0., 0.), (.5, .45, 0.),
                                     (.9, -.7, .7), (.5, 2*math.pi/3, .4))]
    with mp.workdps(60):
        for model in models:
            for winding in (-1, 0, 1):
                total = mp.mpf(model.length)+mp.mpf(model.delay)
                pole = ((2*winding+1)*mp.pi-mp.mpf(model.phase) +
                        1j*mp.log(mp.mpf(model.modulus)))/total
                assert mp.im(pole) < 0
                residual = abs(1+mp.mpf(model.modulus)*mp.exp(
                    1j*mp.mpf(model.phase)+1j*pole*total))
                maxima["lower_half_plane_pole_residual"] = max(
                    maxima["lower_half_plane_pole_residual"], float(residual))
            for energy in (-3.2, -.7, 0., .63, 2.9):
                data = model.scattering(energy)
                vertex = data["vertex"]
                unitary_error = max(abs(sum(vertex[k][i].conjugate()*vertex[k][j]
                                            for k in range(2)) - (i == j))
                                    for i in range(2) for j in range(2))
                maxima["vertex_unitarity"] = max(maxima["vertex_unitarity"], unitary_error)
                maxima["reservoir_unitarity"] = max(maxima["reservoir_unitarity"],
                                                    abs(abs(data["effective"]) - 1))
                D = data["determinant"]
                maxima["relative_determinant"] = max(maxima["relative_determinant"],
                    abs(data["effective"] / data["reference"] - D.conjugate()/D))
                density = model.densities(energy)
                ratio = lambda e: _mp_scattering_ratio(model, e)
                independently_differentiated = mp.im(mp.diff(ratio, energy)/ratio(energy))/(2*mp.pi)
                maxima["spectral_phase_derivative"] = max(maxima["spectral_phase_derivative"],
                    float(abs(independently_differentiated-density["interaction"])))
                loop_amplitude = vertex[1][0] / D
                direct_local_density = model.length * abs(loop_amplitude)**2/(2*math.pi)
                maxima["local_mode_density"] = max(maxima["local_mode_density"],
                    abs(direct_local_density - density["loop"]))
                maxima["delay_density_relation"] = max(maxima["delay_density_relation"],
                    abs(density["loop"]-density["loop_reference"]-
                        model.length/model.round_trip_length*density["interaction"]))
                amplitude = -mp.mpf(model.modulus)*mp.exp(
                    1j*(mp.mpf(model.phase)+mp.mpf(energy)*model.round_trip_length))
                series = model.round_trip_length/mp.pi * mp.re(
                    mp.fsum(amplitude**n for n in range(1, 601)))
                tail_bound = (model.round_trip_length/mp.pi *
                              model.modulus**601/(1-model.modulus))
                error = float(abs(series-density["interaction"]))
                assert error <= float(tail_bound)+2e-12
                maxima["winding_series"] = max(maxima["winding_series"], error)
                assert density["loop"] > 0 and density["reservoir_total"] > 0

        examples = []
        for name, model in (
            ("zero_return", LoopReservoir(0, 1)),
            ("partial_NS", LoopReservoir(.5, 1)),
            ("same_probability_positive_phase", LoopReservoir(.5, 1, 2*math.pi/3)),
            ("same_probability_negative_phase", LoopReservoir(.5, 1, -2*math.pi/3)),
            ("same_probability_causal_delay", LoopReservoir(.5, 1, 0, .8)),
            ("closed_NS", LoopReservoir(1, 1))):
            integral = model.vacuum_energy()
            dilog = model.vacuum_energy(quadrature=False)
            difference = float(abs(integral-dilog))
            maxima["imaginary_integral_vs_dilogarithm"] = max(
                maxima["imaginary_integral_vs_dilogarithm"], difference)
            assert difference < 1e-50
            chiral = model.vacuum_energy("complex_chiral", quadrature=False)
            assert abs(integral-2*chiral) < mp.mpf("1e-50")
            stress = model.chiral_stress()
            assert abs(integral-2*model.round_trip_length*stress) < mp.mpf("2e-16")
            examples.append({"name": name, "return_modulus": model.modulus,
                "return_probability": model.modulus**2, "phase_radians": model.phase,
                "physical_loop_length": model.length, "causal_vertex_delay": model.delay,
                "round_trip_length": model.round_trip_length,
                "dirac_interaction_energy": float(integral),
                "complex_chiral_interaction_energy": float(chiral),
                "physical_loop_Tuu_equals_Tvv": float(stress),
                "physical_loop_interaction_energy": float(2*model.length*stress)})
        values = {row["name"]: row for row in examples}
        assert values["zero_return"]["dirac_interaction_energy"] == 0
        negative = values["partial_NS"]["dirac_interaction_energy"]
        assert negative < 0 < values["same_probability_positive_phase"]["dirac_interaction_energy"]
        assert values["same_probability_positive_phase"]["dirac_interaction_energy"] == \
               values["same_probability_negative_phase"]["dirac_interaction_energy"]
        assert abs(values["same_probability_causal_delay"]["dirac_interaction_energy"] /
                   negative - 1/1.8) < 1e-15
        assert len({values[name]["return_probability"] for name in values
                    if name.startswith("same_probability") or name == "partial_NS"}) == 1
        zeta_energy = -4*mp.pi*mp.zeta(-1, mp.mpf("0.5"))
        exact_closed = LoopReservoir(1, 1).vacuum_energy()
        maxima["zeta_normalization"] = float(abs(zeta_energy-exact_closed))
        assert abs(exact_closed+mp.pi/6) < mp.mpf("1e-50")
        assert maxima["zeta_normalization"] < 1e-50

        # Conditional mathematical ceiling, not a demonstrated property of
        # the physical exterior map. Here |R(i*kappa)|=R exp(-kappa*tau)<=1
        # is explicit. The variation changes L with this return map fixed.
        energy_gaps, stress_gaps = [], []
        for modulus in (0., .5, 1.):
            for phase in (0., -2*math.pi/3, 2*math.pi/3):
                for delay in (0., .8):
                    model = LoopReservoir(modulus, 1., phase, delay)
                    for kappa in (.01, .4, 3.):
                        assert modulus*mp.exp(-mp.mpf(kappa)*delay) <= 1
                    energy_gaps.append(model.vacuum_energy(quadrature=False)+mp.pi/6)
                    stress_gaps.append(model.chiral_stress()+mp.pi/12)
        assert min(energy_gaps) >= -mp.mpf("1e-50")
        assert min(stress_gaps) >= -mp.mpf("1e-50")
        conditional_bounds = {
            "hypotheses": ["same filled negative-energy occupation and relative subtraction as this model",
                "scalar return R(i*kappa) is analytic as required for the contour and has modulus at most one",
                "local Tuu=-one_half*dE/dL is applicable, with the external return map fixed under the length variation"],
            "dirac_energy_lower_bound": "E >= -pi/(6L)",
            "chiral_stress_lower_bound": "Tuu >= -pi/(12L^2)",
            "proof": "For x=exp(-kappa L), |1+z|<=1+x and Re[z/(1+z)]<=x/(1+x) whenever |z|<=x<1; integrate these inequalities.",
            "q_channel_conditional_generalization": "A q-by-q contraction obeys |det(I+xR)|<=(1+x)^q, so the same Dirac normalization gives E>=-q*pi/(6L).",
            "physical_status": "real-energy unitarity alone has not established the required imaginary-axis contractivity of the actual exterior return operator",
            "toy_parameters_checked": len(energy_gaps),
            "minimum_energy_above_bound": float(min(energy_gaps)),
            "minimum_stress_above_bound": float(min(stress_gaps)),
            "arithmetic_tolerance": 1e-50}

    invalid = [dict(modulus=-.01, length=1), dict(modulus=1.01, length=1),
               dict(modulus=.5, length=0), dict(modulus=.5, length=-1),
               dict(modulus=.5, length=1, delay=-1),
               dict(modulus=float("nan"), length=1),
               dict(modulus=.5, length=float("inf")),
               dict(modulus=.5, length=1, phase=float("inf")),
               dict(modulus=.5, length=1, delay=float("nan")),
               dict(modulus=True, length=1), dict(modulus=.5, length=1, phase=1j),
               dict(modulus=.5, length=1e308, delay=1e308)]
    for arguments in invalid:
        try:
            LoopReservoir(**arguments)
        except (TypeError, ValueError):
            pass
        else:
            raise AssertionError(f"invalid inputs accepted: {arguments}")
    for action in (lambda: models[0].scattering(float("nan")),
                   lambda: models[0].densities(float("inf")),
                   lambda: LoopReservoir(1, 1).densities(0),
                   lambda: models[0].vacuum_energy("majorana")):
        try:
            action()
        except (TypeError, ValueError):
            pass
        else:
            raise AssertionError("invalid operation accepted")
    assert max(maxima[key] for key in maxima if key not in
               ("imaginary_integral_vs_dilogarithm", "zeta_normalization")) < 2e-12

    return {"schema_version": 1, "status": "passed",
        "source_file": Path(__file__).name,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "units": "hbar=c=1; E and all optical delays use the same time coordinate",
        "scope": {
            "model": "specified unitary loop plus reservoir with constant vertex mixing and causal delay",
            "occupation": "zero-temperature incoming reservoir sea: every E<0 level occupied, chemical potential zero",
            "dirac_counting": "two independent complex chiral copies with conjugate phases; no additional spin factor",
            "reference": "same one-pass delay with return modulus zero; local and vertex reference terms subtracted",
            "closed_limit": "isolated loop levels retained through D or the R->1- limit; reservoir scattering alone omits them",
            "local_stress": "uniform physical-loop stress only in this specified toy; exterior delay carries part of global energy",
            "limitations": ["no actual exterior scattering amplitude or returning state has been computed here",
                "no inference of a complex amplitude or its continuation from return probabilities",
                "not the complete curved-space determinant, local anomaly, or semiclassical Einstein source",
                "phase-dependent occupied ground state is specified; arbitrary fixed charge sectors are not constructed",
                "no claim of physical preparation, stability, traversal, or chronology change"]},
        "equations": {
            "return_amplitude": "r(E)=R exp(i beta+i E tau)",
            "determinant": "D(E)=1+R exp(i beta+i E(L+tau))",
            "relative_density": "Delta rho=-(1/pi) Im d_E log D",
            "dirac_energy": "2 Re Li2(-R exp(i beta))/(pi(L+tau))",
            "local_chiral_stress": "Tuu=Tvv=Re Li2(-R exp(i beta))/(pi(L+tau)^2)",
            "future_interface": "D=det[I-P_throat(E) R_exterior(E)] after a physical unitary completion, state choice and reference subtraction"},
        "checks": {"max_absolute_errors": maxima, "real_energy_test_points": 20,
                   "energy_integrals": len(examples), "invalid_inputs_rejected": len(invalid)+4,
                   "quadrature_decimal_precision": 60,
                   "rounding_scope": "high-precision and floating-point checks; not interval-certified bounds"},
        "conditional_bounds": conditional_bounds,
        "examples": examples,
        "sources": [
            {"url": "https://arxiv.org/pdf/1606.07834", "context": "Dirac graph unitary transition and secular determinant, Eqs18-21"},
            {"url": "https://arxiv.org/pdf/0711.2395", "context": "reference-relative scattering spectral density, Eqs8-13; fermionic factor checked independently here"}],
        "runtime_seconds": time.perf_counter()-start}


if __name__ == "__main__":
    result = run_checks()
    destination = project_path("scattering_vacuum_results.json")
    destination.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    assert json.loads(destination.read_text(encoding="utf-8")) == result
    print(json.dumps({"status": result["status"], "results": destination.name,
                      "source_sha256": result["source_sha256"],
                      "runtime_seconds": result["runtime_seconds"]}, sort_keys=True))
