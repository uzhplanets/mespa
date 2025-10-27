! anchor = 'use phase_separation, only: do_phase_separation'
         use gentle_mixing, only: do_gentle_mixing
! anchor = 'call do_conv_premix(s, ierr)'; replace
                  if (s% do_gentle_mixing) then
                     call do_gentle_mixing(s, ierr)
                  else
                     call do_conv_premix(s, ierr)
                  end if  