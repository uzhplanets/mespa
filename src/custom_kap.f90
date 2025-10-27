module custom_kap
   ! This module modifies the opacity calculation in the star module
   ! to add support for grain opacities, cloud opacities, and opacity windows

   ! kap control parameters:
   ! kap_logical_ctrl(1): Enable/disable grain opacity
   ! kap_logical_ctrl(2): Enable/disable cloud opacity
   ! kap_logical_ctrl(3): Enable/disable the opacity window
   ! kap_ctrl(1): Scaling factor for the grain opacity
   ! kap_ctrl(2): Ramp-up time for the grain opacity
   ! kap_ctrl(3): Cloud-opacity normalisation
   ! kap_ctrl(4): Cloud-deck width
   ! kap_ctrl(5): Cloud-deck location
   ! kap_ctrl(6): Ramp-up time for the cloud opacity
   ! kap_ctrl(7): Scaling factor for the opacity window
   ! kap_ctrl(8): Scale opacity after this amount of time (in years)


   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   private
   public custom_kap_get, custom_opacity_factor

contains


   subroutine custom_opacity_factor(id, ierr)
      use kap_def, only: Kap_General_Info, kap_handles
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer :: k
      real(dp) :: logT
      real(dp), parameter :: w_loc = 3.3  ! in logT
      real(dp), parameter :: w_scale = 0.15  ! in logT
      
      type(star_info), pointer :: s
      type(Kap_General_Info), pointer :: kap_rq
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      kap_rq => kap_handles(s% kap_handle)

      if (.not. kap_rq% kap_logical_ctrl(3)) then
         return
      end if

      if (s% star_age < kap_rq% kap_ctrl(8)) then
         ! do not apply the opacity window before this time
         s% extra_opacity_factor = 1d0
      else
         do k = 1, s% nz
            logT = log10(s% T(k))
            s% extra_opacity_factor(k) = 1d0 &
            - kap_rq% kap_ctrl(7) &
            * exp(-0.5d0 * pow2((logT - w_loc) / w_scale))
         end do
      end if

   end subroutine custom_opacity_factor


   subroutine custom_kap_get( &
      id, k, handle, species, chem_id, net_iso, xa, &
      logRho, logT, &
      lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
      eta, d_eta_dlnRho, d_eta_dlnT, &
      kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr &
      )

      use chem_def, only: chem_isos
      use chem_lib, only: basic_composition_info
      use eos_lib, only: eosDT_get
      use eos_def, only: i_lnPgas, num_eos_basic_results
      use kap_lib, only: kap_ptr, kap_get, kap_get_elect_cond_opacity, &   
         kap_get_radiative_opacity
      use kap_def, only: kap_is_initialized, Kap_General_Info, num_kap_fracs, &
         i_frac_Type2

      ! INPUT
      integer, intent(in) :: id ! star id if available; 0 otherwise
      integer, intent(in) :: k ! cell number or 0 if not for a particular cell
      integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
      integer, intent(in) :: species
      integer, pointer :: chem_id(:) ! maps species to chem id
      ! index from 1 to species
      ! value is between 1 and num_chem_isos
      integer, pointer :: net_iso(:) ! maps chem id to species number
      ! index from 1 to num_chem_isos (defined in chem_def)
      ! value is 0 if the iso is not in the current net
      ! else is value between 1 and number of species in current net
      real(dp), intent(in) :: xa(:) ! mass fractions
      real(dp), intent(in) :: logRho ! density
      real(dp), intent(in) :: logT ! temperature
      real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
      ! free_e := total combined number per nucleon of free electrons and positrons
      real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
      ! eta := electron degeneracy parameter

      ! OUTPUT
      real(dp), intent(out) :: kap_fracs(num_kap_fracs)
      real(dp), intent(out) :: kap ! total opacity
      real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
      real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
      real(dp), intent(out) :: dlnkap_dxa(:) ! partial derivative w.r.t. to species
      integer, intent(out) :: ierr ! 0 means AOK.

      ! variables for the (total) radiative opacity
      real(dp) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
      real(dp) :: kap_rad_tot, dlnkap_rad_tot_dlnRho, dlnkap_rad_tot_dlnT
      real(dp) :: frac_lowT, frac_highT, frac_Type2

      ! variables for the electron conduction opacity
      real(dp) :: kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT

      ! support variables for the opacity calculation
      real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      real(dp) :: XC, XN, XO, XNe
      integer :: i, iz

      ! variables for the grain opacity
      real(dp) :: kap_grains, dlnkap_grains_dlnRho, dlnkap_grains_dlnT
      real(dp) :: f_grains, T_6, logR_bar, logT1_star, logT2_star
      
      ! variables for the cloud opacity
      real(dp) :: P_gas
      real(dp) :: kap_clouds  ! cloud opacity
      real(dp) :: dlnkap_clouds_dlnRho, dlnkap_clouds_dlnT
      real(dp) :: kap_c ! cloud opacity normalisation
      real(dp) :: delta_c  ! cloud thickness parameter
      real(dp) :: P_c  ! cloud deck location

      ! variables for the equation of state call
      real(dp) :: eos_res(num_eos_basic_results)
      real(dp) :: eos_d_dlnd(num_eos_basic_results)
      real(dp) :: eos_d_dlnT(num_eos_basic_results)
      real(dp) :: eos_d_dxa(num_eos_basic_results, species)
      
      type(star_info), pointer :: s
      type(Kap_General_Info), pointer :: kap_rq

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (.not. kap_is_initialized) then
         ierr = -1
         return
      end if

      call kap_ptr(handle, kap_rq, ierr)
      if (ierr /= 0) return

      ! not doing grains or clouds; simply call the default kap_get and return
      if (.not. kap_rq% kap_logical_ctrl(1) .and. .not. kap_rq% kap_logical_ctrl(2)) then
         call kap_get( &
            handle, species, chem_id, net_iso, xa, &
            logRho, logT, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr &
         )
         return
      end if

      call basic_composition_info( &
         species, chem_id, xa, X, Y, Z, &
         abar, zbar, z2bar, z53bar, ye, mass_correction, sumx &
         )

      xc = 0; xn = 0; xo = 0; xne = 0
      do i=1, species
         iz = chem_isos% Z(chem_id(i))
         select case(iz)
         case (6)
            xc = xc + xa(i)
         case (7)
            xn = xn + xa(i)
         case (8)
            xo = xo + xa(i)
         case (10)
            xne = xne + xa(i)
         end select
      end do

      call kap_get_radiative_opacity( &
         handle, X, Z, XC, XN, XO, XNe, logRho, logT, &
         frac_lowT, frac_highT, frac_Type2, kap_rad, & 
         dlnkap_rad_dlnRho, dlnkap_rad_dlnT, ierr &
         )

      call kap_get_elect_cond_opacity( &
         handle, zbar, logRho, logT, &
         kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr &
         )

      if (kap_rq% kap_logical_ctrl(1)) then
         ! include grain opacity
         T_6 = exp10(logT) / pow6(10d0)
         logR_bar = log10(exp10(logRho) / pow3(T_6))
         logT1_star = 0.0245d0 * logR_bar + 3.096d0
         logT2_star = 0.0245d0 * logR_bar + 3.221d0

         if (s% star_age < kap_rq% kap_ctrl(2)) then
            f_grains = kap_rq% kap_ctrl(1) &
               * (max(s% star_age, 1d-3) / kap_rq% kap_ctrl(2))
         else
            f_grains = kap_rq% kap_ctrl(1)
         end if

         if (logT < logT1_star) then
            call get_kap_grains( &
               logT, f_grains, kap_grains, &
               dlnkap_grains_dlnT, dlnkap_grains_dlnRho &
               )
         else if (logT > logT2_star) then
            ! no grains if the temperature is high
            kap_grains = 0d0
            dlnkap_grains_dlnT = 0d0
            dlnkap_grains_dlnRho = 0d0
         else
            ! linearly interpolate in the intermediate region
            call get_kap_grains_linear( &
               logT, logRho, f_grains, logT1_star, logT2_star, &
               kap_grains, dlnkap_grains_dlnT, &
               dlnkap_grains_dlnRho &
               )
         end if
      else
         ! no grain opacity
         kap_grains = 0d0
         dlnkap_grains_dlnT = 0d0
         dlnkap_grains_dlnRho = 0d0
      end if

      if (kap_rq% kap_logical_ctrl(2)) then
         ! include cloud opacity
         call eosDT_get(s% eos_handle, species, chem_id, net_iso, xa, &
            exp10(logRho), logRho, exp10(logT), logT, &
            eos_res, eos_d_dlnd, eos_d_dlnT, eos_d_dxa, ierr)
         if (ierr /= 0) then
            write(*, *) "eosDT_get failed in custom_kap_get"
            return
         end if
         if (s% star_age < kap_rq% kap_ctrl(6)) then
            kap_c = kap_rq% kap_ctrl(3) &
               * max(s% star_age, 1d-3) / kap_rq% kap_ctrl(6)
         else
            kap_c = kap_rq% kap_ctrl(3)
         end if
         delta_c = kap_rq% kap_ctrl(4)
         P_c = kap_rq% kap_ctrl(5)
         P_gas = exp(eos_res(i_lnPgas))
         kap_clouds = get_kap_clouds(P_gas, kap_c, delta_c, P_c) 
         dlnkap_clouds_dlnRho = 2d0 * delta_c * (1 - P_gas / P_c) &
            * P_gas / P_c * eos_d_dlnd(i_lnPgas)
         dlnkap_clouds_dlnT = 2d0 * delta_c * (1 - P_gas / P_c) &
            * P_gas / P_c * eos_d_dlnT(i_lnPgas)
      else
         kap_clouds = 0d0
         dlnkap_clouds_dlnRho = 0d0
         dlnkap_clouds_dlnT = 0d0
      end if

      ! add the grain and cloud opacities to the radiative opacity
      kap_rad_tot = kap_rad + kap_grains + kap_clouds
      dlnkap_rad_tot_dlnRho = 1d0 / kap_rad_tot &
         * (kap_rad * dlnkap_rad_dlnRho &
            + kap_grains * dlnkap_grains_dlnRho &
            + kap_clouds * dlnkap_clouds_dlnRho &
            )
      dlnkap_rad_tot_dlnT = 1d0 / kap_rad_tot &
         * (kap_rad * dlnkap_rad_dlnT &
            + kap_grains * dlnkap_grains_dlnT &
            + kap_clouds * dlnkap_clouds_dlnT &
            )

      ! combine the total radiative opacity (gas, grains and clouds)
      ! with the electron conductive opacity
      call combine_rad_with_conduction( &
         kap_rq, logRho, logT, zbar, &
         kap_rad_tot, dlnkap_rad_tot_dlnRho, dlnkap_rad_tot_dlnT, &
         kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr &
      )
      dlnkap_dxa = 0  ! composition derivatives not implemented

      ! if (kap_clouds > 0) then
      !    write(*, *) "----------------"
      !    ! write(*, *) "kap_tot", kap
      !    ! write(*, *) "dln_kap_tot_dlnRho", dlnkap_dlnRho
      !    ! write(*, *) "dln_kap_tot_dlnT", dlnkap_dlnT
      !    ! write(*, *) "kap_rad_tot", kap_rad_tot
      !    ! write(*, *) "dln_kap_rad_dlnRho", dlnkap_rad_dlnRho
      !    ! write(*, *) "dln_kap_rad_dlnT", dlnkap_rad_dlnT
      !    write(*, *) "kap_ec", kap_ec
      !    write(*, *) "dln_kap_ec_dlnRho", dlnkap_ec_dlnRho
      !    write(*, *) "dln_kap_ec_dlnT", dlnkap_ec_dlnT
      !    write(*, *) "kap_rad", kap_rad
      !    write(*, *) "dln_kap_rad_dlnRho", dlnkap_rad_dlnRho
      !    write(*, *) "dln_kap_rad_dlnT", dlnkap_rad_dlnT
      !    ! write(*, *) "kap_grains", kap_grains
      !    ! write(*, *) "dln_kap_grains_dlnRho", dlnkap_grains_dlnRho
      !    ! write(*, *) "dln_kap_grains_dlnT", dlnkap_grains_dlnT
      !    write(*, *) "kap_clouds", kap_clouds
      !    write(*, *) "dln_kap_clouds_dlnRho", dlnkap_clouds_dlnRho
      !    write(*, *) "dln_kap_clouds_dlnT", dlnkap_clouds_dlnT
      !    write(*, *) "----------------"
      ! end if

   end subroutine custom_kap_get


   real(dp) function get_only_kap_grains(logT_in, f_grains) result(kap_out)
      real(dp), intent(in) :: logT_in, f_grains
      kap_out = f_grains * exp10(0.430d0 + 1.3143d0 * (logT_in - 2.85d0))
   end function get_only_kap_grains


   subroutine get_kap_grains(logT_in, f_grains, kap_grains, &
      dlnkap_grains_dlnT, dlnkap_grains_dlnRho &
      )
      real(dp), intent(in) :: logT_in, f_grains
      real(dp), intent(out) :: kap_grains
      real(dp), intent(out) :: dlnkap_grains_dlnT, dlnkap_grains_dlnRho
      kap_grains = get_only_kap_grains(logT_in, f_grains)
      dlnkap_grains_dlnT = 1.3143d0
      dlnkap_grains_dlnRho = 0d0
   end subroutine get_kap_grains


   subroutine get_kap_grains_linear( &
      logT, logRho, f_grains, logT1, logT2, &
      kap_grains, dlnkap_grains_dlnT, dlnkap_grains_dlnRho &
      )

      real(dp), intent(in) :: logT, logRho, f_grains, logT1, logT2
      real(dp), intent(out) :: kap_grains
      real(dp), intent(out) :: dlnkap_grains_dlnT, dlnkap_grains_dlnRho
      real(dp) :: kap_grains_T1, kap_grains_T2
      real(dp) :: dlogT1_dlogT, dlogT2_dlogT
      real(dp) :: dkap_grains_T1_dlogT, dkap_grains_T2_dlogT
      real(dp) :: dlogT_dlogT
      real(dp) :: dkap_grains_dlogT
      real(dp) :: dlogT1_dlogRho, dlogT2_dlogRho
      real(dp) :: dkap_grains_T1_dlogRho, dkap_grains_T2_dlogRho
      real(dp) :: dlogT_dlogRho
      real(dp) :: dkap_grains_dlogRho
      real(dp), parameter :: c0 = 1.3143d0
      real(dp), parameter :: c1 = 0.0245d0

      kap_grains_T1 = get_only_kap_grains(logT1, f_grains)
      kap_grains_T2 = 0d0  ! get_only_kap_grains(logT2_star, f_grains)
      kap_grains = kap_grains_T1 + (kap_grains_T2 - kap_grains_T1) / (logT2 - logT1) * (logT - logT1)

      dlogT1_dlogT = -3d0 * c1
      dlogT2_dlogT = -3d0 * c1
      dkap_grains_T1_dlogT = kap_grains_T1 * ln10 * c0 * dlogT1_dlogT
      dkap_grains_T2_dlogT = 0d0
      dlogT_dlogT = 1d0
      dkap_grains_dlogT = dkap_grains_T1_dlogT &
                                 + (kap_grains - kap_grains_T1) &
                                 * ((dkap_grains_T2_dlogT - dkap_grains_T1_dlogT) / (kap_grains_T2 - kap_grains_T1) &
                                    - (dlogT2_dlogT - dlogT1_dlogT) / (logT2 - logT1) &
                                    + (dlogT_dlogT - dlogT1_dlogT) / (logT - logT1))
      dlnkap_grains_dlnT = 1d0 / (kap_grains * ln10) * dkap_grains_dlogT

      dlogT1_dlogRho = c1
      dlogT2_dlogRho = c1
      dkap_grains_T1_dlogRho = kap_grains_T1 * ln10 * c0 * dlogT1_dlogRho
      dkap_grains_T2_dlogRho = 0d0  ! kap_grains_T2 * ln10 * c0 * d_log10_T2_d_log10_Rho
      dlogT_dlogRho = 0d0
      dkap_grains_dlogRho = dkap_grains_T1_dlogRho &
                                 + (kap_grains - kap_grains_T1) &
                                 * ((dkap_grains_T2_dlogRho - dkap_grains_T1_dlogRho) / (kap_grains_T2 - kap_grains_T1) &
                                    - (dlogT2_dlogRho - dlogT1_dlogRho) / (logT2 - logT1) &
                                    + (dlogT_dlogRho - dlogT1_dlogRho) / (logT - logT1))
      dlnkap_grains_dlnRho = 1d0 / (kap_grains * ln10) * dkap_grains_dlogRho
   end subroutine get_kap_grains_linear


   real(dp) function get_kap_clouds(P, kap_c, delta_c, P_c) result(kap_clouds)
      real(dp), intent(in) :: P, kap_c, delta_c, P_c
      kap_clouds = kap_c &
         * exp(-delta_c * pow2((1 - P / P_c)))
   end function get_kap_clouds


   subroutine combine_rad_with_conduction( &
         rq, logRho, logT, zbar, &
         kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
         kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr &
      )
      use utils_lib, only: is_bad, mesa_error
      
      type (Kap_General_Info), pointer :: rq
      real(dp), intent(in) :: logRho, logT, zbar
      real(dp), intent(in) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
      real(dp), intent(in) :: kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT
      real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
      integer, intent(out) :: ierr  ! 0 means AOK.
      logical, parameter :: dbg = .false.

      include 'formats'

      ierr = 0

      if (is_bad(kap_rad)) then
         write(*,*) 'kap_rad', kap_rad
         call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
      end if
      if (dbg) write(*, 1) 'kap_rad', kap_rad
      if (dbg) write(*, 1) 'dlnkap_rad_dlnRho', dlnkap_rad_dlnRho
      if (dbg) write(*, 1) 'dlnkap_rad_dlnT', dlnkap_rad_dlnT
      if (dbg) write(*, *)

      if (.not. rq% include_electron_conduction) then
         kap = kap_rad
         dlnkap_dlnRho = dlnkap_rad_dlnRho
         dlnkap_dlnT = dlnkap_rad_dlnT
         return
      end if

      if (is_bad(kap_ec)) then
         write(*,*) 'kap_ec', kap_ec
         call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
      end if
      if (dbg) write(*, 1) 'kap_ec', kap_ec
      if (dbg) write(*, 1) 'dlnkap_ec_dlnRho', dlnkap_ec_dlnRho
      if (dbg) write(*, 1) 'dlnkap_ec_dlnT', dlnkap_ec_dlnT
      if (dbg) write(*, *)

      kap = 1d0 / (1d0 / kap_rad + 1d0 / kap_ec)
      if (dbg) write(*, 1) 'kap_rad', kap_rad
      if (dbg) write(*, 1) 'kap', kap
      if (dbg) write(*, 1) 'log10(kap)', log10(kap)

      if (is_bad(kap)) then
         ierr = -1; return
         write(*, 1) 'kap', kap
         call mesa_error(__FILE__,__LINE__,'Get_kap_Results')
      end if

      dlnkap_dlnRho = (kap / kap_rad) * dlnkap_rad_dlnRho &
         + (kap / kap_ec) * dlnkap_ec_dlnRho

      if (is_bad(dlnkap_dlnRho)) then
         ierr = -1; return
         write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
         write(*,1) 'kap', kap
         write(*,1) 'dkap_dlnRho', kap * dlnkap_dlnRho
         write(*,1) 'dkap_ec_dlnRho', kap_ec * dlnkap_ec_dlnRho
         write(*,1) 'dkap_rad_dlnRho', kap_rad * dlnkap_rad_dlnRho
         write(*,1) 'kap_rad', kap_rad
         write(*,1) 'kap_ec', kap_ec
         call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
      end if

      dlnkap_dlnT = (kap / kap_rad) * dlnkap_rad_dlnT &
         + (kap / kap_ec) * dlnkap_ec_dlnT

      if (is_bad(dlnkap_dlnT)) then
         ierr = -1; return
         write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
         write(*,1) 'kap', kap
         write(*,1) 'dkap_dlnT', kap * dlnkap_dlnT
         write(*,1) 'dkap_ec_dlnT', kap_ec * dlnkap_ec_dlnT
         write(*,1) 'dkap_rad_dlnT', kap_rad * dlnkap_rad_dlnT
         write(*,1) 'kap_rad', kap_rad
         write(*,1) 'kap_ec', kap_ec
         call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
      end if

      if (dbg) write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
      if (dbg) write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
      if (dbg) call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')

   end subroutine combine_rad_with_conduction

end module custom_kap
