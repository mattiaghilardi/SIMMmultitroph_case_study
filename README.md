
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SIMMmultitroph_case_study

<!-- badges: start -->

<!-- badges: end -->

The goal of this project is to reproduce all results and figures of the
case study in the paper

> TITLE AND AUTHORS
> <!-- > Title: Best practices for estimating consumer reliance on basal food sources with bulk stable isotopes -->
> <!-- > -->
> <!-- > Authors: Ghilardi M, Morais RA, Brandl SJ, Casey JM, Mercière A, Morat F, Schiettekatte NMD, Kayal M, Letourneur Y, Parravicini V -->

## Instructions

The project uses `targets` as pipeline tool to manage the entire project
workflow, `renv` for package management and `crew` for parallel
computing. It depends on R version 4.4.2 (2024-10-31), JAGS version
4.3.2, which is used to fit Bayesian mixing models, and CmdStan version
2.35.0, which is used to fit Bayesian regression models with `brms`.

> [!NOTE] 
> If JAGS is not installed on the machine, you can install it
> from <https://mcmc-jags.sourceforge.io/>.

> [!CAUTION] 
> The project takes several hours to run using parallel
> computing. It was run on a machine with 20 cores and a minimum of 12
> cores are required.
>
> To test whether the project is running correctly, it is recommended to
> modify `.targets.R` by reducing the number of iterations of Bayesian
> mixing models by replacing
> `run = list(chainLength = 200000, ..., calcDIC = TRUE)` with
> `run = "test"` on lines 78 and 90.

To reproduce the project:

1.  Open the R project in RStudio or open an R session with working
    directory set to the root of the project.
2.  Install the required R packages by calling:

``` r
renv::restore()
```

3.  Check the installed version of CmdStan by calling:

``` r
cmdstanr::cmdstan_version()
```

> [!NOTE] 
> If CmdStan is not installed on the machine, you can install
> it by calling:
>
> ``` r
> CmdStanR::install_cmdstan(version = "`r cmdstanr::cmdstan_version()`")
> ```
>
> as explained at <https://mc-stan.org/cmdstanr/articles/cmdstanr.html>.

4.  Run the pipeline by calling:

``` r
targets::tar_make()
```

## Content

This repository is structured as follow:

- [:file_folder:
  data/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/data):
  contains all data used in the analyses.

- [:file_folder:
  derived_data/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/derived_data):
  contains intermediary data required by `MixSIAR`.

- [:file_folder:
  mixing_models/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/mixing_models):
  contains Bayesian mixing models written in JAGS by `MixSIAR`.

- [:file_folder:
  output/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/output):
  contains all the outputs created by the pipeline, including figures
  and tables.

- [:file_folder:
  R/](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/tree/main/R):
  contains R functions.

- [:page_facing_up:
  \_targets.R](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/blob/main/_targets.R):
  project pipeline.

- [:page_facing_up:
  \_targets_packages.R](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/blob/main/_targets_packages.R):
  list of package dependencies created by `targets::tar_renv()` for
  compatibility with `renv`.

- [:page_facing_up:
  .Rprofile](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/blob/main/.Rprofile):
  script to activate `renv` that is automatically executed every time an
  R session starts.

- [:page_facing_up:
  renv.lock](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/blob/main/renv.lock):
  file that records the library used to run the project and makes it
  easier to reinstall it in the future and on different machines.

- [:page_facing_up:
  SIMMmultitroph_simulations.Rproj](https://github.com/mattiaghilardi/SIMMmultitroph_case_study/blob/main/SIMMmultitroph_case_study.Rproj):
  project file used by RStudio to define and manage the R project. Can
  be used as a shortcut for opening the project directly from the
  filesystem.

## Dependency network

The following graph shows the dependency network of the project:

``` mermaid
graph LR
  style Legend fill:#FFFFFF00,stroke:#000000;
  style Graph fill:#FFFFFF00,stroke:#000000;
  subgraph Legend
    direction LR
    x7420bd9270f8d27d([""Up to date""]):::uptodate --- xbf4603d6c2c2ad6b([""Stem""]):::none
    xbf4603d6c2c2ad6b([""Stem""]):::none --- xf0bce276fe2b9d3e>""Function""]:::none
    xf0bce276fe2b9d3e>""Function""]:::none --- x5bffbffeae195fc9{{""Object""}}:::none
  end
  subgraph Graph
    direction LR
    x52b646ad5661aeea>"plot_MixSIAR_continuous"]:::uptodate --> xaeeacd478c4fad9c>"plot_rel_contribution_vs_tl"]:::uptodate
    x33f1329a60edfce2>"mixing_polygon_simulation"]:::uptodate --> x0b99d0ca42f92a63>"run_mixing_polygon_simulation"]:::uptodate
    xf48d3b6327a7cf5e>"run_model_parallel"]:::uptodate --> xd9dbadeaf839ef00>"check_baselines_diet"]:::uptodate
    xf48d3b6327a7cf5e>"run_model_parallel"]:::uptodate --> x68113335f6347c2b>"run_MixSIAR_models"]:::uptodate
    x9af3cebb4b257c6d>"fit_TP_model"]:::uptodate --> xff752a073fe9c83f>"estimate_TP"]:::uptodate
    x2d6d3d16b184a8dd>"set_output_options"]:::uptodate --> x0753d7762c879e81>"make_MixSIAR_stats"]:::uptodate
    x2d6d3d16b184a8dd>"set_output_options"]:::uptodate --> x566596e258c49171>"save_MixSIAR_stats_diag"]:::uptodate
    x530a48c297c1222b>"diag_custom"]:::uptodate --> x728931985b13eb6b>"plot_mixing_model_comparison"]:::uptodate
    xd2c2ad5099df2bb2>"point_custom"]:::uptodate --> x728931985b13eb6b>"plot_mixing_model_comparison"]:::uptodate
    x2099c459998438b9>"cor_custom"]:::uptodate --> x728931985b13eb6b>"plot_mixing_model_comparison"]:::uptodate
    xb473ac5dfee6f55e>"compare_models_parallel"]:::uptodate --> xc98c9990b6cc6a30>"select_best_models"]:::uptodate
    xca670edab52a76d5>"output_JAGS_custom"]:::uptodate --> x566596e258c49171>"save_MixSIAR_stats_diag"]:::uptodate
    xc825fbce21649751>"plot_data_two_iso_custom"]:::uptodate --> x54e9877832c0f69f>"plot_isospace_mixsiar"]:::uptodate
    x421a64d350c50cc6>"check_isotopes_across_years"]:::uptodate --> x7341aad9693c7473(["fish_isotope_check"]):::uptodate
    x8cb65925a21b734a(["sia_fish"]):::uptodate --> x7341aad9693c7473(["fish_isotope_check"]):::uptodate
    xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate --> xc859af69afe1465e(["MixSIAR_stats_diag_Post"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> xc859af69afe1465e(["MixSIAR_stats_diag_Post"]):::uptodate
    x566596e258c49171>"save_MixSIAR_stats_diag"]:::uptodate --> xc859af69afe1465e(["MixSIAR_stats_diag_Post"]):::uptodate
    x8d837fef00b2e817>"prepare_TDF_data"]:::uptodate --> xc6672540d1e7777e(["TDF"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> xc6672540d1e7777e(["TDF"]):::uptodate
    x6531b1ee1ea59353(["file_sia_sources"]):::uptodate --> xee0490802e3f9c6f(["sia_sources"]):::uptodate
    x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate --> x0fde4fd6ad542f8b(["MixSIAR_rel_contributions_vs_tl_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> x0fde4fd6ad542f8b(["MixSIAR_rel_contributions_vs_tl_McCutchan"]):::uptodate
    xaeeacd478c4fad9c>"plot_rel_contribution_vs_tl"]:::uptodate --> x0fde4fd6ad542f8b(["MixSIAR_rel_contributions_vs_tl_McCutchan"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> xed03bbb4665843f3(["MixSIAR_isospace_and_rel_contributions_Post"]):::uptodate
    xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate --> xed03bbb4665843f3(["MixSIAR_isospace_and_rel_contributions_Post"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> xed03bbb4665843f3(["MixSIAR_isospace_and_rel_contributions_Post"]):::uptodate
    x2452010f7f119092>"plot_isospace_and_rel_contributions"]:::uptodate --> xed03bbb4665843f3(["MixSIAR_isospace_and_rel_contributions_Post"]):::uptodate
    xd2d3f5375a95ccfa(["baselines_diet"]):::uptodate --> x5a91308517f5ebfd(["TP"]):::uptodate
    xff752a073fe9c83f>"estimate_TP"]:::uptodate --> x5a91308517f5ebfd(["TP"]):::uptodate
    x8cb65925a21b734a(["sia_fish"]):::uptodate --> x5a91308517f5ebfd(["TP"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> x5a91308517f5ebfd(["TP"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> x7bac5e9ae53a7775(["MixSIAR_summary_McCutchan"]):::uptodate
    x17612f01167adf3d>"make_MixSIAR_summary_table"]:::uptodate --> x7bac5e9ae53a7775(["MixSIAR_summary_McCutchan"]):::uptodate
    x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate --> x7bac5e9ae53a7775(["MixSIAR_summary_McCutchan"]):::uptodate
    xed8121b94a5eae48(["file_sia_fish"]):::uptodate --> x8cb65925a21b734a(["sia_fish"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate
    xc98c9990b6cc6a30>"select_best_models"]:::uptodate --> xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate
    x0753d7762c879e81>"make_MixSIAR_stats"]:::uptodate --> x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate
    x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate --> x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate
    x1113c56f2c0e3d44>"make_model_comparison_table"]:::uptodate --> x43ae57e26738935a(["MixSIAR_comparison_Post"]):::uptodate
    xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate --> x43ae57e26738935a(["MixSIAR_comparison_Post"]):::uptodate
    x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate --> x9342d817833f1295(["MixSIAR_rel_contributions_taxonomy_McCutchan"]):::uptodate
    xc1948ba637fb6d45>"plot_all_rel_contribution"]:::uptodate --> x9342d817833f1295(["MixSIAR_rel_contributions_taxonomy_McCutchan"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate
    x68113335f6347c2b>"run_MixSIAR_models"]:::uptodate --> xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate
    xc6672540d1e7777e(["TDF"]):::uptodate --> xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> x7a97f41ca1b5aad8(["MixSIAR_summary_Post"]):::uptodate
    x17612f01167adf3d>"make_MixSIAR_summary_table"]:::uptodate --> x7a97f41ca1b5aad8(["MixSIAR_summary_Post"]):::uptodate
    xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate --> x7a97f41ca1b5aad8(["MixSIAR_summary_Post"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate
    x68113335f6347c2b>"run_MixSIAR_models"]:::uptodate --> xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate
    xc6672540d1e7777e(["TDF"]):::uptodate --> xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate
    x0738df6c8eb6329b(["file_guilds"]):::uptodate --> x01b7ba2d65a46f54(["trophic_guilds"]):::uptodate
    x2d21ff1824673907(["consumers"]):::uptodate --> xd2dc6df2f1239860(["simulation_mixing_polygon"]):::uptodate
    x0b99d0ca42f92a63>"run_mixing_polygon_simulation"]:::uptodate --> xd2dc6df2f1239860(["simulation_mixing_polygon"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> xd2dc6df2f1239860(["simulation_mixing_polygon"]):::uptodate
    xc6672540d1e7777e(["TDF"]):::uptodate --> xd2dc6df2f1239860(["simulation_mixing_polygon"]):::uptodate
    xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate --> xffb370e30be3607c(["MixSIAR_rel_contributions_vs_tl_Post"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> xffb370e30be3607c(["MixSIAR_rel_contributions_vs_tl_Post"]):::uptodate
    xaeeacd478c4fad9c>"plot_rel_contribution_vs_tl"]:::uptodate --> xffb370e30be3607c(["MixSIAR_rel_contributions_vs_tl_Post"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> x0b5b3adaf1fdf7ce(["MixSIAR_isospace_Post"]):::uptodate
    x54e9877832c0f69f>"plot_isospace_mixsiar"]:::uptodate --> x0b5b3adaf1fdf7ce(["MixSIAR_isospace_Post"]):::uptodate
    x951bfa52b1ec58c7>"prepare_consumer_data"]:::uptodate --> x2d21ff1824673907(["consumers"]):::uptodate
    x8cb65925a21b734a(["sia_fish"]):::uptodate --> x2d21ff1824673907(["consumers"]):::uptodate
    x5a91308517f5ebfd(["TP"]):::uptodate --> x2d21ff1824673907(["consumers"]):::uptodate
    x01b7ba2d65a46f54(["trophic_guilds"]):::uptodate --> x2d21ff1824673907(["consumers"]):::uptodate
    x0753d7762c879e81>"make_MixSIAR_stats"]:::uptodate --> xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate
    xe872beb09fd20a3b(["MixSIAR_best_models_Post"]):::uptodate --> xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate
    xd65f8845d077b9bc(["MixSIAR_models_Post"]):::uptodate --> xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate
    xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate --> xc47165df53baa141(["MixSIAR_rel_contributions_taxonomy_Post"]):::uptodate
    xc1948ba637fb6d45>"plot_all_rel_contribution"]:::uptodate --> xc47165df53baa141(["MixSIAR_rel_contributions_taxonomy_Post"]):::uptodate
    x2d21ff1824673907(["consumers"]):::uptodate --> xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate
    xd464358f22aa0b2a>"remove_outliers"]:::uptodate --> xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate
    xd2dc6df2f1239860(["simulation_mixing_polygon"]):::uptodate --> xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate
    x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate --> x7ff40d374e80fa37(["MixSIAR_stats_diag_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> x7ff40d374e80fa37(["MixSIAR_stats_diag_McCutchan"]):::uptodate
    x566596e258c49171>"save_MixSIAR_stats_diag"]:::uptodate --> x7ff40d374e80fa37(["MixSIAR_stats_diag_McCutchan"]):::uptodate
    x1113c56f2c0e3d44>"make_model_comparison_table"]:::uptodate --> xad22377feb7cecb1(["MixSIAR_comparison_McCutchan"]):::uptodate
    x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate --> xad22377feb7cecb1(["MixSIAR_comparison_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> x0552aff7858b1c08(["MixSIAR_isospace_McCutchan"]):::uptodate
    x54e9877832c0f69f>"plot_isospace_mixsiar"]:::uptodate --> x0552aff7858b1c08(["MixSIAR_isospace_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate
    xc98c9990b6cc6a30>"select_best_models"]:::uptodate --> x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate
    xc6c78c2be3de1bfd(["consumers_clean"]):::uptodate --> xa4e6c98d947ea246(["MixSIAR_isospace_and_rel_contributions_McCutchan"]):::uptodate
    x650cc0c834d289d4(["MixSIAR_best_models_McCutchan"]):::uptodate --> xa4e6c98d947ea246(["MixSIAR_isospace_and_rel_contributions_McCutchan"]):::uptodate
    xb2c63c1dbcf3dd6c(["MixSIAR_models_McCutchan"]):::uptodate --> xa4e6c98d947ea246(["MixSIAR_isospace_and_rel_contributions_McCutchan"]):::uptodate
    x2452010f7f119092>"plot_isospace_and_rel_contributions"]:::uptodate --> xa4e6c98d947ea246(["MixSIAR_isospace_and_rel_contributions_McCutchan"]):::uptodate
    x31d5311b21dbef23(["MixSIAR_stats_McCutchan"]):::uptodate --> x48c566c6d3dabe92(["MixSIAR_comparison_rel_contributions"]):::uptodate
    xc20c67c0eabb23a2(["MixSIAR_stats_Post"]):::uptodate --> x48c566c6d3dabe92(["MixSIAR_comparison_rel_contributions"]):::uptodate
    x728931985b13eb6b>"plot_mixing_model_comparison"]:::uptodate --> x48c566c6d3dabe92(["MixSIAR_comparison_rel_contributions"]):::uptodate
    x631c79307bb6a1f1>"prepare_source_data"]:::uptodate --> xf053c63c658d7219(["sources"]):::uptodate
    xee0490802e3f9c6f(["sia_sources"]):::uptodate --> xf053c63c658d7219(["sources"]):::uptodate
    xd9dbadeaf839ef00>"check_baselines_diet"]:::uptodate --> xd2d3f5375a95ccfa(["baselines_diet"]):::uptodate
    xf053c63c658d7219(["sources"]):::uptodate --> xd2d3f5375a95ccfa(["baselines_diet"]):::uptodate
    xc89be4ed5763b132{{"controller"}}:::uptodate --> xc89be4ed5763b132{{"controller"}}:::uptodate
  end
  classDef uptodate stroke:#000000,color:#ffffff,fill:#354823;
  classDef none stroke:#000000,color:#000000,fill:#94a4ac;
  linkStyle 0 stroke-width:0px;
  linkStyle 1 stroke-width:0px;
  linkStyle 2 stroke-width:0px;
  linkStyle 101 stroke-width:0px;
```

## Computational environment

    #> R version 4.4.2 (2024-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 22.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    #>  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    #>  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: Atlantic/Canary
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices datasets  utils     methods   base     
    #> 
    #> other attached packages:
    #> [1] targets_1.6.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] base64url_1.4        compiler_4.4.2       renv_1.0.5          
    #>  [4] rjags_4-15           tidyselect_1.2.1     callr_3.7.6         
    #>  [7] yaml_2.3.8           fastmap_1.2.0        lattice_0.22-6      
    #> [10] coda_0.19-4.1        R6_2.5.1             generics_0.1.3      
    #> [13] igraph_2.0.3         distributional_0.4.0 knitr_1.45          
    #> [16] backports_1.4.1      checkmate_2.3.1      tibble_3.2.1        
    #> [19] pillar_1.9.0         posterior_1.5.0      rlang_1.1.4         
    #> [22] utf8_1.2.4           xfun_0.43            cli_3.6.3           
    #> [25] magrittr_2.0.3       ps_1.7.6             digest_0.6.35       
    #> [28] grid_4.4.2           processx_3.8.4       rstudioapi_0.16.0   
    #> [31] cmdstanr_0.8.1       secretbase_0.4.0     lifecycle_1.0.4     
    #> [34] vctrs_0.6.5          data.table_1.15.4    evaluate_0.23       
    #> [37] glue_1.7.0           tensorA_0.36.2.1     codetools_0.2-20    
    #> [40] abind_1.4-5          fansi_1.0.6          rmarkdown_2.26      
    #> [43] tools_4.4.2          pkgconfig_2.0.3      htmltools_0.5.8.1
