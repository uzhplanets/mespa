! anchor = 'logical        :: do_premix_heating ! include eps_pre_mix in energy equation'
         ! gentle mixing params
         logical        :: do_gentle_mixing
         character(2)   :: conv_premix_fix_vars
         character(len=strlen)    :: gentle_mixing_iso
         real(dp)       :: gentle_mixing_critical_msd_cpm
         logical        :: gentle_mixing_use_adaptive_damping
         real(dp)       :: gentle_mixing_cpm_damping_factor
         real(dp)       :: gentle_mixing_min_msd_cpm
         real(dp)       :: gentle_mixing_critical_msd_solver
         real(dp)       :: gentle_mixing_critical_msd_solver_at_high_Z
         real(dp)       :: gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z
         real(dp)       :: gentle_mixing_critical_dXi_cpm
         real(dp)       :: gentle_mixing_critical_dXi_solver
         real(dp)       :: gentle_mixing_timestep_reduction_threshold_yr
         real(dp)       :: gentle_mixing_maximum_cpm_extend
         logical        :: skip_cpm_if_initial_model_too_unstable
         logical        :: gentle_mixing_skip_boundaries_if_any_unstable
         logical        :: gentle_mixing_reduce_mixing_length