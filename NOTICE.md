Unless otherwise stated the source code and scripts in this project are licensed under the GNU General Public License Version 3, see [LICENSE.md](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/LICENSE.md).

This project includes a modified version of several functions of the `MixSIAR` R package:

- Original Authors: Brian Stock, Brice Semmens, Eric Ward, Andrew Parnell, Andrew Jackson, Donald Phillips
- Project Repository: https://github.com/brianstock/MixSIAR
- Original License: GPL-3 License
- Modified by: Mattia Ghilardi
- Modified on: 2026
- Modifications: 
  - function `run_model()` has been modified to run chains in parallel and renamed `run_model_parallel()`
  - function `compare_models()` has been modified to compare models fitted with parallel chains and renamed `compare_models_parallel()`
  - function `output_JAGS()` has been modified to allow optional suppression of output plots and renamed `output_JAGS_custom()`
  - function `plot_data_two_iso()` has been modified to allow further plot customisation and renamed `plot_data_two_iso_custom()`
  - function `plot_continuous_var()` has been modified to allow further plot customisation and renamed `plot_MixSIAR_continuous()`
  - function `combine_sources()` has been modified to allow optional suppression of output plot and renamed `combine_sources_custom()`

This project also includes a function, `mixing_polygon_simulation()`, whose code has been adapted by Mattia Ghilardi from DataS1 in Smith et al. (2013) https://doi.org/10.1111/2041-210X.12048. Changes include converting the original code to a function and improvements to the figures.

All modified code resides in the [R/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/R) directory.
