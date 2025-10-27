! anchor = 'logical :: do_chem'
         real(dp) :: adaptive_damping_factor, ratio_to_reference
         real(dp), parameter :: mixing_threshold = 1d-4  ! value at which damping starts to become significant
         real(dp), parameter :: mlt_threshold = 1d-3  ! > 1d0 means off
         real(dp), parameter :: no_mixing_threshold = 1d-25  ! 1d-12; value after which we turn off mixing
         real(dp), parameter :: increasing_factor = 10d0
         real(dp), parameter :: decreasing_factor = 0.1d0
         logical,  parameter :: set_to_zero_if_too_low = .true.
         integer :: n_repeat
! anchor = 'use hydro_vars, only: set_vars_if_needed'
         use gentle_mixing, only: adaptive_damping
! anchor = 'do_struct_burn_mix = retry'

         ! gentle mixing - damp mixing if we are in trouble
         if (s% do_gentle_mixing .and. s% mixing_length_alpha_reference /= 0d0 .and. s% mix_factor_reference /= 0d0) then

               ! determine the current damping factor
               if ( s% gentle_mixing_reduce_mixing_length ) then
                  call adaptive_damping(s,  s%mixing_length_alpha , s%mixing_length_alpha_reference, no_mixing_threshold, s%mixing_length_alpha_reference, increasing_factor, decreasing_factor, set_to_zero_if_too_low, mlt_threshold, ierr)
                  if (s%mixing_length_alpha /= s% mixing_length_alpha_reference) write(*,*) '  gentle_mixing: mixing_length_alpha = ', s% mixing_length_alpha, "mixing_length_alpha / mixing_length_alpha_reference = ", s% mixing_length_alpha / s% mixing_length_alpha_reference
               else
                  call adaptive_damping(s,  s%mix_factor , s%mix_factor_reference, no_mixing_threshold, s%mix_factor_reference, increasing_factor, decreasing_factor, set_to_zero_if_too_low, mixing_threshold, ierr)
                  if (s%mix_factor /= s% mix_factor_reference) write(*,*) '  gentle_mixing: mix_factor = ', s% mix_factor, "mix_factor / mix_factor_reference = ", s% mix_factor / s% mix_factor_reference
               end if 
               
         end if