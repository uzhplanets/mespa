! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   ! Globale Variable für das Regenmodell

   real(kind=dp), parameter :: r_jupiter_avg = 6.995d9 ! cm
contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      s% extras_startup => extras_startup
      s% extras_start_step => extras_start_step
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

      s% how_many_extra_history_header_items => how_many_extra_history_header_items
      s% data_for_extra_history_header_items => data_for_extra_history_header_items
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

   end subroutine extras_controls


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

   end subroutine extras_startup


   integer function extras_start_step(id)
      use chem_def, only: chem_isos
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      extras_start_step = 0

   end function extras_start_step


   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going

      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 1
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ! error handling
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      names(1) = 'radius_Jup'
      vals(1) = exp10(s% log_surface_radius) / r_jupiter_avg * Rsun

   end subroutine data_for_extra_history_columns


   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 2
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      integer :: i_mass_Jup, i_radius_Jup

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      i_mass_Jup = 1
      i_radius_Jup = i_mass_Jup + 1

      names(i_mass_Jup) = 'mass_Jup'
      vals(1:s%nz,i_mass_Jup) = s%m(1:s%nz) / m_jupiter

      names(i_radius_Jup) = 'radius_Jup'
      vals(1:s%nz,i_radius_Jup) = s% r(1:s%nz) / r_jupiter_avg


   end subroutine data_for_extra_profile_columns


   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 0
   end function how_many_extra_history_header_items


   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

   end subroutine data_for_extra_history_header_items


   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 0
   end function how_many_extra_profile_header_items


   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

   end subroutine data_for_extra_profile_header_items


   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_check_model can do that.
   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer :: order_of_magnitude
      real :: time_floor, time_10_percent_floor_old, time_10_percent_floor_new

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going

      ! is there a star_age already in log10?
      s% xtra(1) = log10(s% time_years)

      if ((floor(s% xtra_old(1)) - floor(s% xtra(1)) .ne. 0)) then

         s% need_to_save_profiles_now = .true.
         s% need_to_update_history_now = .true.

         ! Check for half order of magnitude (5 * 10^n) increases
      else if (mod(s% xtra_old(1), 1d0) < 0.7d0 .and. mod(s% xtra(1), 1d0) >= 0.7d0) then

         ! we don't need to save the early stages of the star
         order_of_magnitude = floor(s% xtra_old(1))
         if (order_of_magnitude < 5) return

         ! Time has increased by half an order of magnitude
         s% need_to_save_profiles_now = .true.
         s% need_to_update_history_now = .true.

         ! the following code saves profiles every 10% of the age above 1e9 years
      else
         ! Check for 10% increases within the same order of magnitude
         order_of_magnitude = floor(s% xtra_old(1))
         ! exit if the order of magnitude is too small
         if (order_of_magnitude < 9) return

         time_floor = 10.0 ** real(order_of_magnitude)
         time_10_percent_floor_old = floor(10**(s% xtra_old(1))/time_floor)
         time_10_percent_floor_new = floor(10**(s% xtra(1))/time_floor)

         if (time_10_percent_floor_old /= time_10_percent_floor_new) then
            ! Time has increased by 10% within the same order of magnitude
            s% need_to_save_profiles_now = .true.
            s% need_to_update_history_now = .true.
         endif
      endif

      ! see extras_check_model for information about custom termination codes
      ! by default, indicate where (in the code) MESA terminated
      if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
   end function extras_finish_step


   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine extras_after_evolve

end module run_star_extras


