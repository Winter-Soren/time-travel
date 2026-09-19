"""Verify that the reported outcome agrees with all saved computations."""
import json
from pathlib import Path
p=Path(__file__).resolve().parent
load=lambda n: json.loads((p/n).read_text())
s=load('symbolic_results.json'); n=load('numerical_results.json')
c=load('mode_cases.json'); q=load('kinetic_repair_results.json')
qc=load('kinetic_repair_cases.json'); ini=load('initial_energy_results.json')
assert s['all_16_Ricci_components_zero'] and s['massless_wave_equation_all_terms_zero']
assert len(c)==n['count']==144
assert sum(i['mass']>0 for i in c)==108
assert all(i['success'] and i['energy_amplification']>1e20 for i in c)
assert n['max_independent_exact_error']<2e-9
assert n['independent_solver_relative_error']<1e-9
assert n['regular_branch_control']['bounded']
assert len(qc)==q['algebraic_grid_cases']==2916
assert all(i['has_nonconstant_divergent_branch'] for i in qc)
assert ini['energy_in_dimensionless_model_units']>0
assert ini['quadrature_refinement_difference']<1e-9
for name in ('research_note.md','verify_geometry.py','run_modes.py',
             'test_kinetic_repair.py','check_initial_energy.py'):
    assert (p/name).is_file() and (p/name).stat().st_size>0
print('REPORT CHECKS PASS. The tested candidate FAILS the bounded-response diagnostic.')
print('No stable time-travel medium or nonlinear outcome has been established.')
