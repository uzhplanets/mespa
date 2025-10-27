! anchor = '! gfortran seems to require "save" here.  at least it did once upon a time.'
      ! gentle mixing
      integer, parameter :: mixing_without_dt_reduction = 0
      integer, parameter :: mixing_with_dt_reduction = 1
      integer, parameter :: gentle_mixing_default_modes = 0
      integer, parameter :: gentle_mixing_PD_MODE = 1
      integer, parameter :: gentle_mixing_PS_MODE = 2