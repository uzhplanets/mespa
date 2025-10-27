! anchor = 'use paquette_coeffs, only: free_collision_integrals'
         use element_sedimentation, only: shutdown_phase_diagram
! anchor = 'call free_collision_integrals()'
         call shutdown_phase_diagram()