# Modules for Experiments in Stellar and Planetary Astrophysics (MESPA)

Welcome! Modules for Experiments in Stellar and Planetary Astrophysics (MESPA) implements various modifications to the [MESA](https://docs.mesastar.org/) stellar evolution code to make it suitable for sophisticated planetary models. It adds a custom equation of state (`custom_eos`), opacity modifications (`custom_kap`), a new mixing scheme (`gentle_mixing`) and support for element sedimentation (`element_sedimentation`)

## Getting started with the custom equation of state and opacity modifications
This quickstart guide assumes you have installed MESA and have some familiarity with it.

To get started, following the following steps:
1. Unzip the files into a directory of your choice. This is where your models will live.
2. Set the `MESPA_DIR` environment variable to point towards your installation directory; `export MESPA_DIR=PATH_TO_MESPA`.
3. Run `./mk` in `$MESPA_DIR`. This creates the planet executable that you will use to run your evolution models. You may have to give executable permissions to the `mk` and `clean` files.

In the `test_suite` folder, there are two test suite examples to get you started with using the custom equation of state and opacity subroutines. For example, you can run the custom_eos_test test suite by navigating into its folder and running `./rn`.

### Custom equation of state (`custom_eos`)
The custom equation of state (eos) is enabled by setting the following pointers in the run_star_extras.f90 file:

```Fortran
    s% eos_rq% other_eos_frac => custom_eos_frac
    s% eos_rq% other_eos_component => custom_eos_component
    s% eos_rq% other_eos_results => custom_eos_results
```

You also have to set the following options in the inlist of your project:
```Fortran
    ! use_other_eos_component
    ! ~~~~~~~~~~~~~~

    ! Set to .true. to enable the custom eos.
    ! ::
    
    use_other_eos_component = .true.
    
    ! use_other_eos_results
    ! ~~~~~~~~~~~~~~

    ! Set to .true. to enable modifying the eos results.
    ! ::
    use_other_eos_results = .true.

    ! eos_integer_ctrl(1)
    ! ~~~~~~~~~~~~~~

    ! Which equation of state to use. Options are (depending on the availability of tables):
    !   1: A low-temperature extension of SCvH for hydrogen and helium and the QEOS heavy-element equation of state for water.
    !   2: A low-temperature extension of SCvH for hydrogen and helium and the QEOS heavy-element equation of state for a mixture 50-50 mixture of water and rocks.
    !   3: The Chabrier & Debras (2021) hydrogen-helium equation of state and the QEOS heavy-element equation of state for water.
    !   4: The Chabrier & Debras (2021) hydrogen-helium equation of state and the QEOS heavy-element equation of state for a 50-50 mixture of water and rocks.
    !   5: A user-generated table with prefix “eosdt_extra_tables”.

    ! ::
    eos_integer_ctrl(1) = 4

    ! eos_integer_ctrl(2)
    ! ~~~~~~~~~~~~~~

    ! Which interpolation method to use to interpolate between the heavy-element tables. Options are 1 for linear, and 2 for cubic.

    ! ::
    eos_integer_ctrl(2) = 2
    
    ! eos_integer_ctrl(3)
    ! ~~~~~~~~~~~~~~

    ! In-table interpolation method. Options are 1 for the MESA default (bi-cubic), 2 for bi-linear if the specific heats (Cp or Cv) are negative, and 3 for always bi-linear.

    ! ::
    eos_integer_ctrl(3) = 1
```

### Custom opacity module (`custom_kap`)
To-do: Add description.


## `gentle_mixing`

Welcome! `gentle_mixing` is a module for [MESA](https://docs.mesastar.org/) that limits the impact of mixing within a single timestep. It is based on the study [Knierim & Helled 2024](https://arxiv.org/abs/2407.09341). If `gentle_mixing` is useful for your research, please consider citing this paper. :)

### Quickstart

Before installing, make sure `$MESA_DIR` is set to your main MESA directory. Then, run the following commands:

1.  First, clone the repository
2.  `cd gentle_mixing`
3.  Run `python setup.py install`
4.  `cd` into your MESA work directory and run `./clean && ./mk`.

That's it! You can now use `gentle_mixing` in your MESA projects. You can find the new options that `gentle_mixing` introduces either below or in the `controls.defaults` file inside the `$MESA_DIR/star/default` directory (after installation).

If you get annoyed by the extra output by `gentle_mixing`, set `TRACE_GENTLE_MIXING = .false.` in `gentle_mixing.f90` and recompile MESA and your work directory. (And yes, I will make this an inlist option in the future. :))

### Introduction

Especially in planetary evolution, mixing heavy elements too rapidly can lead to poor convergence (or prevent it entirely). `gentle_mixing` addresses this by checking the potential impact of mixing during convective premixing (CPM) and solver iterations, limiting that impact to a user-defined threshold.

The module evaluates the mixing impact using the mean squared displacement (MSD) between two subsequent models:

$$\sigma_Z^2 = \int_{0}^{M} (Z(t+\Delta t) - Z(t))^2 \frac{dm}{M}$$

where $Z$ is the mass fraction of the heavy element, $M$ is the total mass of the star, and $\Delta t$ is the timestep. Note that during CPM, no time elapses; instead, we compare $Z$ after CPM to $Z$ before CPM.

You can set two critical MSD values: one for CPM and one for the solver iterations. If the calculated MSD exceeds this critical value, the module reduces the mixing efficiency to ensure the change does not surpass this limit. If the efficiency is reduced, the algorithm also reduces the timestep according to a damping function.

#### Adaptive Mode

For added stability, `gentle_mixing` also features an **adaptive mode** (controlled by `gentle_mixing_use_adaptive_damping`). When enabled, the code attempts to mix *up to* the specified critical MSD tolerance. If the model fails to converge (crashes), it automatically retries the timestep, but this time using tighter, more restrictive tolerances to find a stable solution.

#### CPM Modes

In addition, this module introduces two new modes during CPM to hold the pressure and density ('PD') or the pressure and total entropy ('PS') of the convective region constant. These modes help mitigate the temperature spikes often encountered during CPM for planets. You can find more details in the paper [Knierim & Helled 2024](https://arxiv.org/abs/2407.09341).

### Common Issues and Tips

Ideally, you'll want to set the damping tolerances as high as possible and use both CPM and MESA's normal mixing modes. Finding a converging set of parameters can be tricky, though. Here are some tips that might help:

* At least in the planetary context, `conv_premix_fix_vars = 'PS'` is often the most stable option. Keep in mind though this option requires two root-finding steps for the pressure and entropy, which can slow down the code.
* Make sure your timesteps are sufficiently large before **turning** on mixing. "Sufficiently large" depends on your problem, but in the planetary context, timesteps below $10^2$ years often cause convergence issues.
* If you have convergence issues during CPM, try **disabling** it first and see if the solver alone can converge. If it still fails, try running your model without any mixing. Does it converge then? If not, the problem is likely not related to `gentle_mixing`. If it *does* converge without mixing, first try enabling mixing *without* CPM by slowly varying the solver parameters until you reach your desired values. Then, repeat the process for CPM.
* The default MESA tolerances (e.g., gold tolerances) are **often** too tight for planetary problems with mixing. Try loosening them, but be careful: If they are too loose, the model might not be accurate anymore and you may encounter other convergence issues.
* If you use adaptive mode, start with relatively loose tolerances for both CPM and the solver. In principle, the code will tighten them automatically if needed. However, if the initial tolerances are too tight, the code might not be able to find a stable solution at all. So experiment a little bit!

### Options
The Options for `gentle_mixing` are set in the `controls` namelist. They are:
```Fortran

      ! do_gentle_mixing
      ! ~~~~~~~~~~~~~~

      ! Set to .true. to limit how much the abundace gradient can change per time step.
      ! This is useful for avoiding numerical instabilities in the convective premixing scheme.

      ! ::

    do_gentle_mixing = .false.

      ! conv_premix_fix_vars
      ! ~~~~~~~~~~~~~~
      
      ! Set to the variables that should be kept constant during the convective premixing iterations.
      ! Options are 'PT', 'PD', 'PS', and 'DT'.
      ! 'PT' keeps pressure and temperature constant.
      ! 'PD' keeps pressure and density constant.
      ! 'PS' keeps pressure and total entropy over the mixed region constant.
      ! 'DT' keeps density and temperature constant.
      ! This option replaces the old ``conv_premix_fix_pgas`` control.
      
      ! ::
    
    conv_premix_fix_vars = 'PT'
    
      ! gentle_mixing_iso
      ! ~~~~~~~~~~~~~~
      
      ! Sets the chemical species that is investigated to determine the maximum change in the abundance profile.
      
      ! ::
    
    gentle_mixing_iso = 'o16'
    
      ! gentle_mixing_critical_msd_cpm
      ! ~~~~~~~~~~~~~~
      
      ! Maximum allowed change in the mean square deviation of the abundance profile during convective premixing iterations.
      
      ! ::
    
    gentle_mixing_critical_msd_cpm = 1d-3
    
      ! gentle_mixing_use_adaptive_damping
      ! ~~~~~~~~~~~~~~

      ! If true, the solver will reduce gentle_mixing_critical_msd_cpm and mix_factor as a function of the number of retries.
      ! If false, the gentle_mixing_critical_msd_cpm in CPM will be constant while mix_factor will be dampened with the timestep.

      ! ::
  
    gentle_mixing_use_adaptive_damping = .true.

      ! gentle_mixing_cpm_damping_factor
      ! ~~~~~~~~~~~~~~

      ! Damping factor for the convective premixing iterations if gentle_mixing_use_adaptive_damping is set to .true.

      ! ::

    gentle_mixing_cpm_damping_factor = 0.5d0

      ! gentle_mixing_min_msd_cpm
      ! ~~~~~~~~~~~~~~
      
      ! Minimum allowed change in the mean square deviation of the abundance profile during the solver iterations.
      ! If gentle_mixing_use_adaptive_damping is set to .true., the solver won't damp the mixing if the MSD is below this value.

      ! ::

    gentle_mixing_min_msd_cpm = 1d-4


      ! gentle_mixing_critical_msd_solver
      ! ~~~~~~~~~~~~~~

      ! Maximum allowed change in the mean square deviation of the abundance profile during the solver iterations.

      ! ::

    gentle_mixing_critical_msd_solver = 1d-3


      ! gentle_mixing_critical_msd_solver_at_high_Z
      ! ~~~~~~~~~~~~~~

      ! Maximum allowed change in the mean square deviation of the abundance profile if Z in any cell that convects is above a threshold.

      ! ::

    gentle_mixing_critical_msd_solver_at_high_Z = 1d-3

    
      ! gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z
      ! ~~~~~~~~~~~~~~

      ! Threshold for the metallicity above which gentle_mixing_critical_msd_solver_at_high_Z is used. If <0 then it is not used.

      ! ::

    gentle_mixing_min_Z_for_critical_msd_solver_at_high_Z = -1d0


      ! gentle_mixing_critical_dXi_cpm
      ! ~~~~~~~~~~~~~~

      ! Maximum allowed change in the absolute abundance of a cell during the convective premixing iterations.

      ! ::

    gentle_mixing_critical_dXi_cpm = 1d0

      ! gentle_mixing_critical_dXi_solver
      ! ~~~~~~~~~~~~~~

      ! Maximum allowed change in the absolute abundance of a cell during the solver iterations.

      ! ::

    gentle_mixing_critical_dXi_solver = 1d0

      ! gentle_mixing_timestep_reduction_threshold_yr
      ! ~~~~~~~~~~~~~~

      ! No timestep reduction will be done if the timestep is smaller than this value.
      ! The luminosity tends to become unsteady if the timestep is too small.
      ! This option helps to avoid this by not reducing the timestep further if it is already small.
      ! The right value depends on the evolution timescale of your problem.

      ! ::

    gentle_mixing_timestep_reduction_threshold_yr = 1d4

      ! gentle_mixing_maximum_cpm_extend
      ! ~~~~~~~~~~~~~~

      ! Maximum extend of the convective premixing region.
      ! If a convective region expands outside of this limit, the convective premixing is stopped.
      ! The solver, however, will still be able to mix the region.
      ! This helps when there is an abundace discontinuity between the envelope and the atmosphere.
      ! Sometimes, convective premixing will lead to a negative surface luminoisty, which is unphysical.
      ! To avoid this, this option allows to ignore the convective premixing in the outer layers.

      ! ::

    gentle_mixing_maximum_cpm_extend = 1d0

      ! skip_cpm_if_initial_model_too_unstable
      ! ~~~~~~~~~~~~~~

      ! If true, checks the MSD of `gentle_mixing_iso` from mixing the existing convective regions fully (i.e., without advancing any boundary).
      ! This is also done during CPM while advancing a convective boundary.
      ! If `MSD > gentle_mixing_critical_msd_cpm` CPM will be skipped. 
      ! Typically, you want to set this to `.false.` 
      ! However, sometimes the model converges better if you mix most of the material first just using the solver and only later CPM.

      ! ::

    skip_cpm_if_initial_model_too_unstable = .false.

      ! gentle_mixing_skip_boundaries_if_any_unstable
      ! ~~~~~~~~~~~~~~

      ! If true, exits convective premixing if the last boundary that was extended exceeded the critical MSD.
      ! This saves time if the model is already close to the mixing limit and expending any of the other boundaries would also exceed the limit.
      ! Sometimes, however, the model tries to mix a large convective layer early on that already exceeds the limit.
      ! In this case, you want to set this to `.false.`, because you otherwise skip the mixing of the other boundaries.
      ! If unsure, set this to `.false.`.
      
      ! ::

    gentle_mixing_skip_boundaries_if_any_unstable = .false.

      ! gentle_mixing_reduce_mixing_length
      ! ~~~~~~~~~~~~~~

      ! If true, the mixing length - instead of the mixing factor - will be reduced by gentle_mixing.

      ! ::

      gentle_mixing_reduce_mixing_length = .false.
```
For more details on the code and options, see the documentation in `gentle_mixing/docs/index.md`.


