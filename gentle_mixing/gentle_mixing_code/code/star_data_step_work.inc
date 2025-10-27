! anchor = 'logical :: fix_Pgas ! Flag to keep Pgas fixed during cell updates'
      ! used for improved CPM.
      integer  :: fix_eos_result ! Flag to keep eos results fixed during cell updates
      integer  :: gentle_mixing_timestep_reduction_type ! 0 for no reduction, 1 for reduction, (integers to allow for future expansion)
      real(dp) :: mixing_length_alpha_reference
      real(dp) :: mix_factor_reference
      real(dp) :: gentle_mixing_critical_msd_cpm_reference