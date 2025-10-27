# Details

## New Inlist Options
Below some explanation for the more opaque inlist options.

### Adaptive Damping
When using ` gentle_mixing_use_adaptive_damping = .true.`, `gentle_mixing` tries to find the ideal damping parameters for your model within the range you provide. As with most of the `gentle_mixing` code, this applies to CPM and the solver separately.

#### CPM
For CPM, the algorithm first expands up to `gentle_mixing_critical_msd_cpm` using no damping at all. In case that the solver does not find a valid solution (i.e., does a retry), the `gentle_mixing_critical_msd_cpm` is multiplied by `gentle_mixing_adaptive_damping_factor` and another try is made. This continues until either a solution is found or `gentle_mixing_critical_msd_cpm` drops below `gentle_mixing_adaptive_damping_minimum_msd_cpm`, in which case CPM is skipped for this timestep. In case of a successful CPM step, `gentle_mixing_critical_msd_cpm` is increased again by 1/`gentle_mixing_adaptive_damping_factor` for the next timestep, up to the maximum value provided by the user.

#### Solver
For the solver, the algorithm is similar. First, the solver tries to converge with no damping up to `gentle_mixing_critical_msd_solver`. If it fails, `s$ mix_factor` is initially reduced by `1d-4` and afterwords by `1d-1` on each subsequent retry. There is again a lower boundary beneath which the solver will skip mixing entirely for this timestep, and `mix_factor` is increased again by a factor of `10` (or `1d4` at the end) on successful convergence. The initial greater reduction is to quickly approach a reasonable damping value. The exact parameters can be found in the `struct_burn_mix.f90` file of this repo. They are not yet implemented as an inlist option, but you can modify them there (and recompile afterwards) if needed.


## Files
- `controls.defaults`: Default controls/documentation for the new code.
- `ctrls_io.f90`: File containing the assignment of the new star_info variables.
- `eos_support.f90`: Defining new `solve_eos_given_PgasS` function.
- `evolve.f90`: Adding `gentle_mixing` to the main evolution loop.
- `gentle_mixing.f90`: The main Fortran module implementing gentle mixing based on CPM.
- `init.f90`: Sets some state variables for gentle mixing at initialization.
- `makefile_base`: Add `gentle_mixing.f90` to the list of source files to compile.
- `micro.f90`: Add new CPM modes to `do_eos_for_cell`.
- `star_controls.inc`: Contains the type definitions for the new controls used by the element sedimentation module.
- `star_data_def.f90`: defines new variables for `gentle_mixing`.
- `star_data_step_work.inc`: Other star_type variables that are used by `gentle_mixing` but are not a user input.
- `struct_burn_mix.f90`: Contains the adaptive damping routine for the solver.
- `timestep.f90`: Adding `check_Z_MSD` to keep track of how much mixing has occurred and trigger a redo/retry if necessary.


## New Variables

Below is a compact overview of the new variables and parameters introduced by gentle_mixing.

### Controls (star_info; inlist-controlled)

- `s% do_gentle_mixing` – Enables/disables gentle mixing (CPM with limits).
- `s% conv_premix_fix_vars` – Which variables are held fixed during CPM, one of `'PT'`, `'PD'`, `'PS'`, `'DT'`.
- `s% gentle_mixing_iso` – Isotope (e.g., `'o16'`) whose abundance profile/change is monitored.
- `s% gentle_mixing_critical_msd_cpm` – Critical MSD threshold for CPM.
- `s% gentle_mixing_use_adaptive_damping` – Enables adaptive damping (reduces CPM MSD limit and solver mixing when needed).
- `s% gentle_mixing_cpm_damping_factor` – Damping factor (<1) for reducing the CPM MSD limit on retries.
- `s% gentle_mixing_min_msd_cpm` – Lower bound for the adaptively damped CPM MSD threshold.
- `s% gentle_mixing_critical_msd_solver` – Critical MSD threshold during the solver step.
- `s% gentle_mixing_critical_msd_solver_at_high_Z` – Stricter MSD threshold when Z in the convective region is high.
- `s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z` – Z threshold above which the “high-Z” MSD limit applies (<0 disables).
- `s% gentle_mixing_critical_dXi_cpm` – Max allowed absolute abundance change per cell in CPM.
- `s% gentle_mixing_critical_dXi_solver` – Max allowed absolute abundance change per cell in the solver.
- `s% gentle_mixing_timpestep_reduction_threshold_yr` – Time threshold (years) below which no further dt reduction is applied.
- `s% gentle_mixing_maximum_cpm_extend` – Maximum mass fraction `q` to which CPM may extend a zone.
- `s% skip_cpm_if_initial_model_too_unstable` – Skip CPM if the initial model would already exceed the MSD limit.
- `s% gentle_mixing_skip_boundaries_if_any_unstable` – Abort CPM for all boundaries once any boundary exceeds the limit.
- `s% gentle_mixing_reduce_mixing_length` – Damp the mixing length `mixing_length_alpha` instead of `mix_factor` in the solver.

### Runtime/state variables (star_info)

- `s% fix_eos_result` – Mode switch for special EOS solves in CPM (default/`PD`/`PS`).
- `s% gentle_mixing_timestep_reduction_type` – Flag for dt reduction due to gentle_mixing (0: none, 1: reduce).
- `s% mixing_length_alpha_reference` – Reference mixing length to restore after damping.
- `s% mix_factor_reference` – Reference `mix_factor` to restore after damping.
- `s% gentle_mixing_critical_msd_cpm_reference` – Reference CPM MSD threshold used by adaptive damping.

### Enums/Constants

- `mixing_without_dt_reduction (=0)` – No dt reduction.
- `mixing_with_dt_reduction (=1)` – dt reduction active.
- `gentle_mixing_default_modes (=0)` – Default EOS mode for CPM.
- `gentle_mixing_PD_MODE (=1)` – EOS mode: pressure + density fixed (CPM).
- `gentle_mixing_PS_MODE (=2)` – EOS mode: pressure + entropy fixed (CPM).
- `Tlim_Z_MSD` – New timestep limiter index for the gentle_mixing MSD limit.

### Timestep logic

- `check_Z_MSD` – New dt limiter based on MSD and `dXi`; triggers retry/redo and dt reduction if needed.
- `gentle_mixing_reduce_timestep` – Reduces `dt` via `dt_limit_ratio` when required (controlled by `s% gentle_mixing_timpestep_reduction_threshold_yr`).

### Solver damping (struct_burn_mix.f90)

- `adaptive_damping` – Generic damping scheme for `mixing_length_alpha` or `mix_factor` in the solver; uses:
  - `mixing_threshold = 1d-4`, `mlt_threshold = 1d-3`, `no_mixing_threshold = 1d-25`,
  - `increasing_factor = 10d0`, `decreasing_factor = 0.1d0`,
  - `set_to_zero_if_too_low = .true.`

### EOS/CPM extensions

- Extended EOS paths in `do_eos_for_cell` via `s% fix_eos_result` (`PD`/`PS`).
- `solve_eos_given_PgasS` – New EOS routine for given gas pressure `Pgas` and entropy `S` (supports CPM `'PS'`).

### Module-internal additions (for completeness)

- `TRACE_GENTLE_MIXING` – Debug output for gentle_mixing.
- `pre_cpm_xa(:,:)` – Abundances before CPM for comparison (stability/MSD checks).
- `was_split`, `too_unstable` – State flags used during stability checks/splitting.
- `CMP_MODE` – Selected CPM update mode (one of `FIXED_PD`, `FIXED_PS`, `FIXED_PT`, `FIXED_DT`).


