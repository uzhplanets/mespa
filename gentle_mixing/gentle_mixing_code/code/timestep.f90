! anchor = 'if (return_now(Tlim_dX_nuc_drop)) return'

            if (s% do_gentle_mixing) then
            do_timestep_limits = check_Z_MSD( &
               s, skip_hard_limit, dt, dt_limit_ratio(Tlim_Z_MSD))
            if (return_now(Tlim_Z_MSD)) return
            end if
! anchor = 'end function check_dt_div_min_dr_div_cs'
            
      ! gentle mixing timestep control
      integer function check_Z_MSD(s, skip_hard_limit, dt, dt_limit_ratio)
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         logical, intent(in) :: skip_hard_limit
         real(dp), intent(in) :: dt
         real(dp), intent(inout) :: dt_limit_ratio

         ! local
         real(dp) :: Z_change_MSD, Z_change_MSD_crit, dXi
         integer  :: i_Xi

         logical :: verbose = .true.
         integer :: ierr

         include 'formats'

         ierr = 0 ! 0 means AOK
         check_Z_MSD = keep_going

         if (s% mix_factor == 0d0 .or. .not. s% do_mix) then
            if (verbose) write(*,*) "  gentle_mixing: no mixing due to mix_factor = 0d0 or do_mix = .false."
            return
         end if


         ! critical value for Z_change_MSD that shall not be exceeded
         i_Xi = s%net_iso(chem_get_iso_id(s%gentle_mixing_iso))
         call get_msd_crit(s, i_Xi, Z_change_MSD_crit,  ierr)
         call Xi_change_MSD(i_Xi, Z_change_MSD, ierr)

         dXi = MAXVAL(ABS(s% xa_sub_xa_start(i_Xi,1:s%nz)))
         if (verbose) write(*, *) "  gentle_mixing: maximum change in dXi = ", dXi
         

         if (Z_change_MSD > Z_change_MSD_crit) then
            
            check_Z_MSD = check_Z_MSD_retry_or_redo(s, dt)

            s% gentle_mixing_timestep_reduction_type = mixing_with_dt_reduction ! reduce timestep

            ! restore mix_factor
         else if (dXi > s% gentle_mixing_critical_dXi_solver) then
            if (verbose) write(*, *) "  gentle_mixing: maximum change in dXi exceeds gentle_mixing_critical_dXi_solver = ", s% gentle_mixing_critical_dXi_solver
            
            check_Z_MSD = check_Z_MSD_retry_or_redo(s, dt)

            s% gentle_mixing_timestep_reduction_type = mixing_with_dt_reduction ! reduce timestep

         end if

         ! At the moment, we only distinguish between no reduction and reduction. However, we could
         ! distinguish between reduction in cpm, reduction in solver, and reduction in both.
         select case(s% gentle_mixing_timestep_reduction_type)
         case(mixing_without_dt_reduction)
            ! no reduction
         case(mixing_with_dt_reduction) ! reduction
            call gentle_mixing_reduce_timestep(dt_limit_ratio)
         case default
            ! unknown integer
            write(*,*) "  gentle_mixing: unknown integer for gentle_mixing_timestep_reduction_type"
            ierr = 1
            return
         end select

         ! restore the default
         s% gentle_mixing_timestep_reduction_type = mixing_without_dt_reduction

         if (verbose) then
            write(*,*) "  gentle_mixing: Z_change_MSD = ", Z_change_MSD
            write(*,*) "  gentle_mixing: Z_change_MSD_crit = ", Z_change_MSD_crit
            write(*,*) "  gentle_mixing: retry with smaller mix_factor: ", check_Z_MSD == retry, "redo: ", check_Z_MSD == redo
            write(*,*) ""
         end if

      contains

         function check_Z_MSD_retry_or_redo(s, dt) result(redo_or_retry)
            type (star_info), pointer :: s
            real(dp), intent(in) :: dt
            integer :: redo_or_retry
            logical :: verbose = .false.

            if ( dt <  secyer * s% gentle_mixing_timestep_reduction_threshold_yr) then
               if (verbose) write(*, *) "  gentle_mixing: dt < dt_retry_threshold. dt = ", dt, "dt/gentle_mixing_timestep_reduction_threshold_yr = ", dt/(secyer *s% gentle_mixing_timestep_reduction_threshold_yr)
               redo_or_retry = redo ! redo with smaller mix_factor
            else 
               redo_or_retry = retry ! retry with smaller mix_factor
            end if
            
         end function check_Z_MSD_retry_or_redo
         
         subroutine get_msd_crit(s, i_Xi, Z_change_MSD_crit,  ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: i_Xi
            real(dp), intent(out) :: Z_change_MSD_crit
            integer, intent(out) :: ierr

            ! local
            real(dp) :: Z_conv_max
            logical :: verbose = .false.

            ierr = 0 ! 0 means AOK.

            if (s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z <= 0d0) then
               Z_change_MSD_crit = s% gentle_mixing_critical_msd_solver
            else
               Z_conv_max = MAXVAL(s% xa(i_Xi, 1:s%nz), mask=(s%mixing_type(1:s%nz) == 1))
               if ( Z_conv_max > s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z) then
                  if (verbose) write(*, *) "  gentle_mixing: Z_conv_max >  s% gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z = ", Z_conv_max
                  Z_change_MSD_crit = s% gentle_mixing_critical_msd_solver_at_high_Z
               else
                  Z_change_MSD_crit = s% gentle_mixing_critical_msd_solver
               end if
            end if
            
            ! test that Z_change_MSD_crit is positive
            if (Z_change_MSD_crit < 0d0) then
               ierr = 1
               write(*,*) "get_msd_crit: Z_change_MSD_crit is negative:", Z_change_MSD_crit
            end if
            
         end subroutine get_msd_crit

         subroutine Xi_change_MSD(i_Xi, dXi_MSD, ierr)
            integer, intent(in) :: i_Xi
            real(dp), intent(out) :: dXi_MSD
            integer, intent(out) :: ierr

            ierr = 0
            dXi_MSD = dot_product(pow2(s% xa_sub_xa_start(i_Xi,1:s%nz)), s% dm(1:s%nz)) / s% mstar
            dXi_MSD = sqrt(dXi_MSD)

            if (dXi_MSD < 0d0) then
               ierr = 1
               write(*,*) "Xi_change_MSD: dXi_MSD is negative:", dXi_MSD
            end if
            
         end subroutine Xi_change_MSD

         subroutine gentle_mixing_reduce_timestep(dt_limit_ratio)
            ! timestep reduction for gentle mixing
            real(dp), intent(out) :: dt_limit_ratio

            ! local
            real(dp) :: f
            real(dp) :: log_dt_reduction_limit
            log_dt_reduction_limit = log10(s% gentle_mixing_timestep_reduction_threshold_yr)

            ierr = 0 ! 0 means AOK.

            if (safe_log10(s%dt_years) > log_dt_reduction_limit) then
               f = 0.5
               dt_limit_ratio = 1d0 / f
            end if

            if (dt_limit_ratio < 0d0) then
               ierr = 1
               write(*,*) "gentle_mixing_reduce_timestep: dt_limit_ratio is negative:", dt_limit_ratio
            end if

         end subroutine gentle_mixing_reduce_timestep


      end function check_Z_MSD
