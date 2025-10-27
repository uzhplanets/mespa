# Modules for Experiments in Stellar and Planetary Astrophysics (MESPA)

Welcome! Modules for Experiments in Stellar and Planetary Astrophysics (MESPA) implements various modifications to the [MESA](https://docs.mesastar.org/) stellar evolution code to make it suitable for sophisticated planetary models. It adds a custom equation of state (`custom_eos`), opacity modifications (`custom_kap`), a new mixing scheme (`gentle_mixing`) and support for element sedimentation (`element_sedimentation`)

## Getting started with the custom equation of state and opacity modifications
This quickstart guide assumes you have installed MESA and have some familiarity with it.

To get started, following the following steps:
1. Unzip the files into a directory of your choice. This is where your models will live.
2. Set the `MESPA_DIR` environment variable to point towards your installation directory; `export MESPA_DIR=PATH_TO_MESPA`.
3. Run `./mk` in `$MESPA_DIR`. This creates the planet executable that you will use to run your evolution models. You may have to give executable permissions to the `mk` and `clean` files.

In the `test_suite` folder, there are two test suite examples to get you started with using the custom equation of state and opacity subroutines. For example, you can run the custom_eos_test test suite by navigating into its folder and running `./rn`.

## Custom equation of state (`custom_eos`)
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

## Custom opacity module (`custom_kap`)
To-do: Add description.
