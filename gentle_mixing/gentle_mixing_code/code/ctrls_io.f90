! anchor = 'conv_premix_fix_pgas, conv_premix_dump_snapshots, do_premix_heating, &'
    do_gentle_mixing, conv_premix_fix_vars, gentle_mixing_iso, gentle_mixing_critical_msd_cpm, gentle_mixing_use_adaptive_damping, gentle_mixing_cpm_damping_factor, gentle_mixing_min_msd_cpm, gentle_mixing_critical_msd_solver,gentle_mixing_critical_msd_solver_at_high_Z, gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z, gentle_mixing_critical_dXi_cpm, gentle_mixing_critical_dXi_solver, gentle_mixing_timestep_reduction_threshold_yr, gentle_mixing_maximum_cpm_extend, skip_cpm_if_initial_model_too_unstable, gentle_mixing_skip_boundaries_if_any_unstable, gentle_mixing_reduce_mixing_length, &
! anchor = 's% do_premix_heating = do_premix_heating'
 ! gentle mixing
 s% do_gentle_mixing = do_gentle_mixing
 s% conv_premix_fix_vars = conv_premix_fix_vars
 s% gentle_mixing_iso = gentle_mixing_iso
 s% gentle_mixing_critical_msd_cpm = gentle_mixing_critical_msd_cpm
 s% gentle_mixing_use_adaptive_damping = gentle_mixing_use_adaptive_damping
 s% gentle_mixing_cpm_damping_factor = gentle_mixing_cpm_damping_factor
 s% gentle_mixing_min_msd_cpm = gentle_mixing_min_msd_cpm
 s% gentle_mixing_critical_msd_solver = gentle_mixing_critical_msd_solver
 s% gentle_mixing_critical_msd_solver_at_high_Z = gentle_mixing_critical_msd_solver_at_high_Z
 s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z = gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z
 s% gentle_mixing_critical_dXi_cpm = gentle_mixing_critical_dXi_cpm
 s% gentle_mixing_critical_dXi_solver = gentle_mixing_critical_dXi_solver
 s% gentle_mixing_timestep_reduction_threshold_yr = gentle_mixing_timestep_reduction_threshold_yr
 s% gentle_mixing_maximum_cpm_extend = gentle_mixing_maximum_cpm_extend
 s% skip_cpm_if_initial_model_too_unstable = skip_cpm_if_initial_model_too_unstable
 s% gentle_mixing_skip_boundaries_if_any_unstable = gentle_mixing_skip_boundaries_if_any_unstable
 s% gentle_mixing_reduce_mixing_length = gentle_mixing_reduce_mixing_length
! anchor = 'do_premix_heating = s% do_premix_heating'
 ! gentle mixing
 do_gentle_mixing = s% do_gentle_mixing
 conv_premix_fix_vars = s% conv_premix_fix_vars
 gentle_mixing_iso = s% gentle_mixing_iso
 gentle_mixing_critical_msd_cpm = s% gentle_mixing_critical_msd_cpm
 gentle_mixing_use_adaptive_damping = s% gentle_mixing_use_adaptive_damping
 gentle_mixing_cpm_damping_factor = s% gentle_mixing_cpm_damping_factor
 gentle_mixing_min_msd_cpm = s% gentle_mixing_min_msd_cpm
 gentle_mixing_critical_msd_solver = s% gentle_mixing_critical_msd_solver
 gentle_mixing_critical_msd_solver_at_high_Z = s% gentle_mixing_critical_msd_solver_at_high_Z
 gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z = s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z
 gentle_mixing_critical_dXi_cpm = s% gentle_mixing_critical_dXi_cpm
 gentle_mixing_critical_dXi_solver = s% gentle_mixing_critical_dXi_solver
 gentle_mixing_timestep_reduction_threshold_yr = s% gentle_mixing_timestep_reduction_threshold_yr
 gentle_mixing_maximum_cpm_extend = s% gentle_mixing_maximum_cpm_extend
 skip_cpm_if_initial_model_too_unstable = s% skip_cpm_if_initial_model_too_unstable
 gentle_mixing_skip_boundaries_if_any_unstable = s% gentle_mixing_skip_boundaries_if_any_unstable
 gentle_mixing_reduce_mixing_length = s% gentle_mixing_reduce_mixing_length

! anchor = 'call check_controls(s, ierr)'
 ! reference values for gentle mixing
 s% mixing_length_alpha_reference = s% mixing_length_alpha
 s% mix_factor_reference = s% mix_factor
 s% gentle_mixing_critical_msd_cpm_reference = s%gentle_mixing_critical_msd_cpm