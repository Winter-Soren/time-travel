# Directory guide

Research files keep their existing paths: scripts load neighboring results and
data, reports link to figures, and saved audits name checkpoints explicitly.
Use this index to navigate the project. Run existing research commands from the
project root; consult the reports before starting potentially expensive runs.

## Start here

| Material | Location |
| --- | --- |
| Current work and continuity | [TASK_STATE.md](TASK_STATE.md) |
| Research results and reproduction notes | [RESEARCH.md](RESEARCH.md) |
| Long-throat finite-window integral | [long_throat_integral.py](long_throat_integral.py), [results](long_throat_integral_results.json), [figure](long_throat_integral_progress.png) |
| Infinite-angular direct-trace bound | [radial_angular_envelope.py](radial_angular_envelope.py), [results](radial_angular_envelope_results.json); conditional collar-feedback and threshold arguments in [RESEARCH.md](RESEARCH.md) |
| First massive mode on the saved smooth collar | [smooth_collar_mode.py](smooth_collar_mode.py), [results](smooth_collar_mode_results.json); reuses [dirac_return_closure.py](dirac_return_closure.py) through an optional profile argument |
| Longer-throat sampled source | [long_throat_source.py](long_throat_source.py), [results](long_throat_source_results.json), [figure](long_throat_source_progress.png) |
| Defect-weighted angular remainder | [angular_defect_bound.py](angular_defect_bound.py), [results](angular_defect_bound_results.json); proof and assumptions in [RESEARCH.md](RESEARCH.md) |
| Complete finite angular-source calculation | [complete_angular_source.py](complete_angular_source.py), [pilot and integration evidence](complete_angular_source_results.json), [figure](complete_angular_source_progress.png); continuum-state and angular-cutoff arguments in [RESEARCH.md](RESEARCH.md) |
| Three-sector reference sum and complete angular bases | [angular_reference_sum.py](angular_reference_sum.py), [saved results](angular_reference_sum_results.json), [complete_angular_modes.py](complete_angular_modes.py), [figure](angular_reference_sum_progress.png) |
| Invariant AdS reference and interior source | [ads_reference_source.py](ads_reference_source.py), [saved results](ads_reference_source_results.json), [figure](ads_reference_source_progress.png) |
| Geometric angular connection and missing sectors | [azimuthal_modes.py](azimuthal_modes.py), [connected source calculation](azimuthal_source.py), [saved results](azimuthal_source_results.json), [figure](azimuthal_source_progress.png) |
| Opposite-chirality boundary matching and combined response | [dirac_sewing_audit.py](dirac_sewing_audit.py), [saved physical and adjoint controls](dirac_sewing_audit_results.json), [figure](dirac_sewing_progress.png) |
| Throat pressure and energy-flow diagnostics | [throat_null_response.py](throat_null_response.py), [saved metric derivatives and checks](throat_null_response_results.json), [figure](throat_null_progress.png) |
| Compact local source diagnostic | [local_lapse_response.py](local_lapse_response.py), [saved responses and direct perturbation checks](local_lapse_response_results.json), [figure](local_lapse_progress.png) |
| Computed phase-dependent spectral diagnostic | [phase_vacuum_contrast.py](phase_vacuum_contrast.py), [saved determinants and integrals](phase_vacuum_contrast_results.json), [figure](phase_vacuum_progress.png) |
| Curved modes and junction support | [curved_throat_modes.py](curved_throat_modes.py), [curved_exterior_closure.py](curved_exterior_closure.py), [closure results](curved_exterior_closure_results.json); [throat_junction_audit.py](throat_junction_audit.py), [required-stress results](throat_junction_audit_results.json) |
| Imaginary response and throat connection | [imaginary_weyl_leads.py](imaginary_weyl_leads.py), [imaginary_exterior.py](imaginary_exterior.py), [connected results](imaginary_exterior_results.json); [exterior_loop_closure.py](exterior_loop_closure.py), [closure results](exterior_loop_closure_results.json), [curved-mode derivation](CURVED_THROAT_MATCHING.md) |
| Magnetic exterior continuation | [magnetic_tail_map.py](magnetic_tail_map.py), [tail checks](magnetic_tail_map_results.json); [exterior_tail_transport.py](exterior_tail_transport.py), [connected results](exterior_tail_transport_results.json), [figure](magnetic_tail_progress.png) |
| Connected two-mouth exterior | [exterior_boundary_map.py](exterior_boundary_map.py), [saved results](exterior_boundary_map_results.json); [lead modes](spherical_weyl_leads.py), [mesh](two_mouth_mesh.py) |
| Incoming spinor modes | [INCOMING_MODES.md](INCOMING_MODES.md), [incoming_modes.py](incoming_modes.py), [saved results](incoming_modes_results.json) |
| Exact angular collar modes | [dipole_collar_modes.py](dipole_collar_modes.py), [saved results](dipole_collar_modes_results.json) |
| Local radial correction bound | [collar_dressing_bound.py](collar_dressing_bound.py), [saved results](collar_dressing_bound_results.json); derivation in [RESEARCH.md](RESEARCH.md) |
| Open-boundary scattering | [OPEN_BOUNDARY_SCATTERING.md](OPEN_BOUNDARY_SCATTERING.md), [open_boundary_spinor.py](open_boundary_spinor.py), [saved results](open_boundary_spinor_results.json) |
| Quantum-source analysis | [QUANTUM_COMPLETION.md](QUANTUM_COMPLETION.md) |
| Source-energy analysis | [SOURCE_FEASIBILITY.md](SOURCE_FEASIBILITY.md) |
| Fiction-inspired diagnostics | [FICTION_PHYSICS.md](FICTION_PHYSICS.md) |

## Research families

| Family | Source and retained artifacts at the root |
| --- | --- |
| Coupled models | `coupled_misner.py`, `coupled_ori.py`, corresponding `*_results.json` |
| Ori constraints and evolution | `disk_basis.py`, `ori_disk_constraints.py`, `ori_evolution.py`, `ori_large_analysis.py`; `ori_*`, `smooth_weak_*`, `polynomial_weak_*` data and comparisons |
| Casimir and quantum source | `casimir_cavity.py`, `ori_quantum_source.py`, corresponding `*_results.json` |
| Pulses and transport | `mmp_pulse.py`, `future_time.py`, corresponding `*_results.json` |
| Response and moving exterior | `polyakov_response.py`, `moving_exterior_maps.py`, `moving_exterior_nonlinear.py`; `polyakov_*`, `distributed_*`, `moving_exterior_*` results |
| Magnetic and Dirac calculations | `magnetic_*.py`, `dirac_*.py`, corresponding `*_results.json` |
| Spinful return | [spinful_return.py](spinful_return.py), [saved results](spinful_return_results.json), [figure](spinful_return_progress.png) |
| Figures | Root `*_progress.png` / `*_progress.pdf` pairs and Ori evolution figures; generated by `quantum_horizon_plot.py`, `casimir_source_plot.py`, and `ori_large_analysis.py` |

The PNG and PDF versions serve different viewing/export purposes; both are
retained. JSON results, NPZ states, coefficient snapshots, failed-case records,
and checkpoint metadata are research evidence, not disposable temporary files.

## Existing directories

- [chronology_onset_audit/](chronology_onset_audit/): original audit sources,
  saved computations, report, run log, and SHA-256 manifest; kept together.
- [provenance/](provenance/): historical evolution source used for provenance,
  plus `ori_initial_data.cpython-312.pyc`, an unidentified historic
  **CPython 3.12-only** artifact without a matching source. It was archived
  from `__pycache__/ori_initial_data.cpython-312.pyc` with its SHA-256 unchanged;
  retain it for provenance, not as a supported executable entry point.
- [maintenance/](maintenance/): cleanup records. The
  [20 September 2026 manifest](maintenance/cleanup-2026-09-20.json) records exact
  deleted paths, byte counts, hashes, retained-file inventory, and verification.
- `.venv/`: local Python environment; retained untouched.
- `.cache/` and `__pycache__/`: runtime caches may reappear automatically.

## Cleanup boundary

On 20 September 2026, cleanup removed 16 CPython bytecode files with retained
source counterparts and one generated Matplotlib font cache: **17 files,
2,284,220 bytes**. One orphan bytecode file was moved into `provenance/` as
described above. No scientific source, result, checkpoint, report, or figure
was moved or deleted. No scientific calculations were rerun for cleanup.

Accepted coefficient snapshots and the audit run log were retained even where
content is redundant: their filenames participate in provenance or references.
Checkpoint files remain necessary to saved audits and resumption workflows.
Future relocation requires updating and verifying those consumers together.

Initial verification matched all 369 retained files to their pre-cleanup
SHA-256 hashes. The subsequent archival check confirmed the bytecode's hash;
active `incoming_modes.py` changed concurrently and was left untouched by
cleanup. The final manifest records observed concurrent changes separately.
These seven files already differed from `chronology_onset_audit/sha256_manifest.json`:

- `chronology_onset_audit/symbolic_results.json`
- `chronology_onset_audit/mode_cases.json`
- `chronology_onset_audit/sample_trajectories.json`
- `chronology_onset_audit/kinetic_repair_cases.json`
- `chronology_onset_audit/initial_energy_results.json`
- `chronology_onset_audit/numerical_results.json`
- `chronology_onset_audit/kinetic_repair_results.json`

The cause was not investigated during cleanup. Both the files and historical
manifest were preserved without alteration; no integrity repair is claimed.
