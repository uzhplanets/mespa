! anchor = 'character (len=32) :: phase_separation_option'
         
         ! element sedimentation parameters
         logical        :: do_element_sedimentation
         character(len=strlen)   :: element_sedimentation_method
         character(len=strlen)   :: element_sedimentation_velocity_method
         real(dp)       :: element_sedimentation_velocity
         real(dp)       :: element_sedimentation_advection_diffusion_time_step_factor
         real(dp)       :: element_sedimentation_maximum_dXi_per_cell
         character(len=strlen)   :: element_sedimentation_phase_diagram
         real(dp)       :: element_sedimentation_phase_diagram_temperature_shift
         logical        :: element_sedimentation_extrapolate_phase_diagram
         character(len=strlen)   :: element_sedimentation_chemical_species
