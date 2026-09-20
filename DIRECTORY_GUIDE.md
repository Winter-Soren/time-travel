# Directory guide

Files are grouped by purpose. Run commands from the project root and consult
the reports before starting potentially expensive calculations. Every chat
session must follow [AGENTS.md](AGENTS.md) when creating or relocating files.

| Folder | Contents |
| --- | --- |
| [scripts/](scripts/) | Python calculations, shared modules, and plot generators |
| [results/](results/) | Saved JSON results and numerical reports |
| [data/](data/) | NPZ states, initial coefficients, and resumable checkpoints |
| [figures/](figures/) | PNG previews and PDF exports |
| [docs/](docs/) | Research notes, derivations, and companion reports |
| [maintenance/](maintenance/) | Verification and cleanup/migration records |
| [provenance/](provenance/) | Preserved historical source artifacts |

Use `project_path` from `scripts/project_paths.py` to locate project assets in
Python. Direct commands use `python scripts/name.py`. For `python -c` imports,
set `$env:PYTHONPATH = (Join-Path (Get-Location) 'scripts')` in PowerShell.

The [layout migration](maintenance/layout-migration.json) records every old/new
path and hash. Saved JSON, NPZ, and figure bytes were preserved. Historical JSON
recipes still describe their original layout; translate their paths using this
mapping. Source-only path edits are recorded explicitly, and new runs record
actual new hashes. Existing Ori checkpoints retain strict source checks:
use their original source revision for exact historical resumption.

Run `python -B maintenance/verify_layout.py` for the lightweight layout check.

## Start here

| Material | Location |
| --- | --- |
| Current work and continuity | [TASK_STATE.md](TASK_STATE.md) |
| Research results and reproduction notes | [docs/RESEARCH.md](docs/RESEARCH.md) |
| Quantum channels and curved axial response | [scripts/melvin_leads.py](scripts/melvin_leads.py), [results](results/melvin_leads_results.json), [figure](figures/melvin_channel_progress.png); [scripts/dihole_axial_tail.py](scripts/dihole_axial_tail.py) and [results](results/dihole_axial_tail_results.json) |
| Gravitating exterior and attachment data | [scripts/balanced_dihole.py](scripts/balanced_dihole.py), [results](results/balanced_dihole_results.json), [figure](figures/balanced_exterior_progress.png); [RN compatibility gate](scripts/rn_exterior_gate.py) and [results](results/rn_exterior_gate_results.json) |
| Smooth collar, magnetic response and same-geometry support target | [scripts/smooth_magnetic_collar.py](scripts/smooth_magnetic_collar.py), [results](results/smooth_magnetic_collar_results.json), [figure](figures/smooth_collar_progress.png); conditional RN exterior gate and coupled Dirac operator in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Finite magnetic matching on the actual mouth data | [scripts/magnetic_mouth_matching.py](scripts/magnetic_mouth_matching.py), [results](results/magnetic_mouth_matching_results.json), [figure](figures/magnetic_matching_progress.png); source definition and the remaining smooth-collar requirements in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Actual conforming field and failed residual target | [scripts/conforming_exterior_trial.py](scripts/conforming_exterior_trial.py), [saved field/results](results/conforming_exterior_trial_results.json), [figure](figures/conforming_trial_progress.png); [round-domain quadrature](scripts/physical_exterior_quadrature.py) and [controls](results/physical_exterior_quadrature_results.json) |
| Two-leg angular feedback and continuum-residual target | [scripts/feedback_residual_bound.py](scripts/feedback_residual_bound.py), [results](results/feedback_residual_bound_results.json), [diagram](figures/tail_completion_progress.png); finite mixing witness plus requirements for conforming physical trial fields |
| Complete retained low-frequency enclosure | [scripts/complete_threshold_bound.py](scripts/complete_threshold_bound.py), [results](results/complete_threshold_bound_results.json); all9 retained sectors, conditional finite N4/J9/DG result |
| Magnetic-tail angular quadrature error | [scripts/tail_quadrature_bound.py](scripts/tail_quadrature_bound.py), [results](results/tail_quadrature_bound_results.json); exact-GL truncation versus exact finite-J projection, separate from roundoff |
| All-angular high-frequency source bound | [scripts/ultraviolet_source_bound.py](scripts/ultraviolet_source_bound.py), [results](results/ultraviolet_source_bound_results.json), [frequency ledger](figures/source_remainder_progress.png) |
| One-sector zero-frequency enclosure | [scripts/threshold_return_bound.py](scripts/threshold_return_bound.py), [results](results/threshold_return_bound_results.json); finite N4/J9/DG model and non-interval floating evaluation |
| Buffered angular-feedback estimate | [scripts/buffered_angular_bound.py](scripts/buffered_angular_bound.py), [results](results/buffered_angular_bound_results.json); explicit decreasing bound, numerically too loose |
| Long-throat finite-window integral | [scripts/long_throat_integral.py](scripts/long_throat_integral.py), [results](results/long_throat_integral_results.json), [figure](figures/long_throat_integral_progress.png) |
| Infinite-angular direct-trace bound | [scripts/radial_angular_envelope.py](scripts/radial_angular_envelope.py), [results](results/radial_angular_envelope_results.json); conditional collar-feedback and threshold arguments in [docs/RESEARCH.md](docs/RESEARCH.md) |
| First massive mode on the saved smooth collar | [scripts/smooth_collar_mode.py](scripts/smooth_collar_mode.py), [results](results/smooth_collar_mode_results.json); reuses [scripts/dirac_return_closure.py](scripts/dirac_return_closure.py) through an optional profile argument |
| Longer-throat sampled source | [scripts/long_throat_source.py](scripts/long_throat_source.py), [results](results/long_throat_source_results.json), [figure](figures/long_throat_source_progress.png) |
| Defect-weighted angular remainder | [scripts/angular_defect_bound.py](scripts/angular_defect_bound.py), [results](results/angular_defect_bound_results.json); proof and assumptions in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Complete finite angular-source calculation | [scripts/complete_angular_source.py](scripts/complete_angular_source.py), [pilot and integration evidence](results/complete_angular_source_results.json), [figure](figures/complete_angular_source_progress.png); continuum-state and angular-cutoff arguments in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Three-sector reference sum and complete angular bases | [scripts/angular_reference_sum.py](scripts/angular_reference_sum.py), [saved results](results/angular_reference_sum_results.json), [scripts/complete_angular_modes.py](scripts/complete_angular_modes.py), [figure](figures/angular_reference_sum_progress.png) |
| Invariant AdS reference and interior source | [scripts/ads_reference_source.py](scripts/ads_reference_source.py), [saved results](results/ads_reference_source_results.json), [figure](figures/ads_reference_source_progress.png) |
| Geometric angular connection and missing sectors | [scripts/azimuthal_modes.py](scripts/azimuthal_modes.py), [connected source calculation](scripts/azimuthal_source.py), [saved results](results/azimuthal_source_results.json), [figure](figures/azimuthal_source_progress.png) |
| Opposite-chirality boundary matching and combined response | [scripts/dirac_sewing_audit.py](scripts/dirac_sewing_audit.py), [saved physical and adjoint controls](results/dirac_sewing_audit_results.json), [figure](figures/dirac_sewing_progress.png) |
| Throat pressure and energy-flow diagnostics | [scripts/throat_null_response.py](scripts/throat_null_response.py), [saved metric derivatives and checks](results/throat_null_response_results.json), [figure](figures/throat_null_progress.png) |
| Compact local source diagnostic | [scripts/local_lapse_response.py](scripts/local_lapse_response.py), [saved responses and direct perturbation checks](results/local_lapse_response_results.json), [figure](figures/local_lapse_progress.png) |
| Computed phase-dependent spectral diagnostic | [scripts/phase_vacuum_contrast.py](scripts/phase_vacuum_contrast.py), [saved determinants and integrals](results/phase_vacuum_contrast_results.json), [figure](figures/phase_vacuum_progress.png) |
| Curved modes and junction support | [scripts/curved_throat_modes.py](scripts/curved_throat_modes.py), [scripts/curved_exterior_closure.py](scripts/curved_exterior_closure.py), [closure results](results/curved_exterior_closure_results.json); [scripts/throat_junction_audit.py](scripts/throat_junction_audit.py), [required-stress results](results/throat_junction_audit_results.json) |
| Imaginary response and throat connection | [scripts/imaginary_weyl_leads.py](scripts/imaginary_weyl_leads.py), [scripts/imaginary_exterior.py](scripts/imaginary_exterior.py), [connected results](results/imaginary_exterior_results.json); [scripts/exterior_loop_closure.py](scripts/exterior_loop_closure.py), [closure results](results/exterior_loop_closure_results.json), [curved-mode derivation](docs/CURVED_THROAT_MATCHING.md) |
| Magnetic exterior continuation | [scripts/magnetic_tail_map.py](scripts/magnetic_tail_map.py), [tail checks](results/magnetic_tail_map_results.json); [scripts/exterior_tail_transport.py](scripts/exterior_tail_transport.py), [connected results](results/exterior_tail_transport_results.json), [figure](figures/magnetic_tail_progress.png) |
| Connected two-mouth exterior | [scripts/exterior_boundary_map.py](scripts/exterior_boundary_map.py), [saved results](results/exterior_boundary_map_results.json); [lead modes](scripts/spherical_weyl_leads.py), [mesh](scripts/two_mouth_mesh.py) |
| Incoming spinor modes | [docs/INCOMING_MODES.md](docs/INCOMING_MODES.md), [scripts/incoming_modes.py](scripts/incoming_modes.py), [saved results](results/incoming_modes_results.json) |
| Exact angular collar modes | [scripts/dipole_collar_modes.py](scripts/dipole_collar_modes.py), [saved results](results/dipole_collar_modes_results.json) |
| Local radial correction bound | [scripts/collar_dressing_bound.py](scripts/collar_dressing_bound.py), [saved results](results/collar_dressing_bound_results.json); derivation in [docs/RESEARCH.md](docs/RESEARCH.md) |
| Open-boundary scattering | [docs/OPEN_BOUNDARY_SCATTERING.md](docs/OPEN_BOUNDARY_SCATTERING.md), [scripts/open_boundary_spinor.py](scripts/open_boundary_spinor.py), [saved results](results/open_boundary_spinor_results.json) |
| Quantum-source analysis | [docs/QUANTUM_COMPLETION.md](docs/QUANTUM_COMPLETION.md) |
| Source-energy analysis | [docs/SOURCE_FEASIBILITY.md](docs/SOURCE_FEASIBILITY.md) |
| Fiction-inspired diagnostics | [docs/FICTION_PHYSICS.md](docs/FICTION_PHYSICS.md) |

## Research families

| Family | Source and retained artifacts |
| --- | --- |
| Coupled models | `scripts/coupled_misner.py`, `scripts/coupled_ori.py`, corresponding `*_results.json` |
| Ori constraints and evolution | `scripts/disk_basis.py`, `scripts/ori_disk_constraints.py`, `scripts/ori_evolution.py`, `scripts/ori_large_analysis.py`; `ori_*`, `smooth_weak_*`, `polynomial_weak_*` data and comparisons |
| Casimir and quantum source | `scripts/casimir_cavity.py`, `scripts/ori_quantum_source.py`, corresponding `*_results.json` |
| Pulses and transport | `scripts/mmp_pulse.py`, `scripts/future_time.py`, corresponding `*_results.json` |
| Response and moving exterior | `scripts/polyakov_response.py`, `scripts/moving_exterior_maps.py`, `scripts/moving_exterior_nonlinear.py`; `polyakov_*`, `distributed_*`, `moving_exterior_*` results |
| Magnetic and Dirac calculations | `magnetic_*.py`, `dirac_*.py`, corresponding `*_results.json` |
| Spinful return | [scripts/spinful_return.py](scripts/spinful_return.py), [saved results](results/spinful_return_results.json), [figure](figures/spinful_return_progress.png) |
| Figures | Root `*_progress.png` / `*_progress.pdf` pairs and Ori evolution figures; generated by `scripts/quantum_horizon_plot.py`, `scripts/casimir_source_plot.py`, and `scripts/ori_large_analysis.py` |

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

## Earlier cleanup record

On 20 September 2026, cleanup removed 16 CPython bytecode files with retained
source counterparts and one generated Matplotlib font cache: **17 files,
2,284,220 bytes**. One orphan bytecode file was moved into `provenance/` as
described above. No scientific source, result, checkpoint, report, or figure
was moved or deleted during that earlier cache cleanup. The later layout
migration above relocates these files without deleting scientific artifacts.
No scientific calculations were rerun for cache cleanup.

Accepted coefficient snapshots and the audit run log were retained even where
content is redundant: their filenames participate in provenance or references.
Checkpoint files remain necessary to saved audits and resumption workflows.
The later migration updates and verifies those consumers together.

Initial verification matched all 369 retained files to their pre-cleanup
SHA-256 hashes. The subsequent archival check confirmed the bytecode's hash;
active `scripts/incoming_modes.py` changed concurrently and was left untouched by
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
