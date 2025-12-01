library(targets)
library(crew)

# Set up controller to launch 2 workers on local processes
controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 2,
  seconds_idle = 10,
  crashes_error = 10
)

# Set targets options and required packages
tar_option_set(packages = c("here", "cli", "readr", "tibble", "dplyr", "tidyr", "purrr", "rlang",
                            "ggplot2", "ggtext", "patchwork", "ggdist", "ggblend", "lemon", 
                            "ggplotify", "cowplot", "GGally", "tRophicPosition", "MixSIAR", "R2jags",
                            "brms", "tidybayes", "parallel", "cmdstanr", "loo", "splancs", "sp"),
               controller = controller,
               memory = "transient", 
               garbage_collection = TRUE)

# Load functions
tar_source()

# Create required folders
if (!dir.exists("derived_data")) dir.create("derived_data")
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/figures")) dir.create("output/figures")
if (!dir.exists("output/tables")) dir.create("output/tables")
if (!dir.exists("output/MixSIAR")) dir.create("output/MixSIAR")
if (!dir.exists("output/MixSIAR/summary_statistics")) dir.create("output/MixSIAR/summary_statistics")
if (!dir.exists("output/MixSIAR/diagnostics")) dir.create("output/MixSIAR/diagnostics")
if (!dir.exists("mixing_models")) dir.create("mixing_models")

# Pipeline
list(
  # Prepare data MixSIAR ----
  
  # Load isotope data and trophic guilds
  tar_target(file_sia_fish, "data/sia_fish.csv", format = "file"),
  tar_target(sia_fish, read_csv(file_sia_fish)),
  tar_target(file_sia_sources, "data/sia_sources.csv", format = "file"),
  tar_target(sia_sources, read_csv(file_sia_sources)),
  tar_target(file_guilds, "data/trophic_guilds.csv", format = "file"),
  tar_target(trophic_guilds, read_csv(file_guilds)),
  # Check fish isotopic signatures across years
  tar_target(fish_isotope_check, check_isotopes_across_years(sia_fish)),
  # Prepare source and TDF data for MixSIAR
  tar_target(sources, prepare_source_data(sia_sources)),
  tar_target(TDF, prepare_TDF_data(sources)),
  # Check source contribution to potential baselines
  tar_target(baselines_diet, check_baselines_diet(sources, 
                                                  run = "short")),
  # Estimate trophic position
  # NOTE: results of this step and consequently of the rest of the pipeline
  # will differ slightly each time the code runs
  # Modifications to the source code of `tRophicPosition::TPmodel()` 
  # are required to obtain reproducible results
  tar_target(TP, estimate_TP(sia_fish, sources, baselines_diet)),
  # Prepare consumer data for MixSIAR
  tar_target(consumers, prepare_consumer_data(sia_fish, trophic_guilds, TP)),
  # Mixing polygon simulations
  tar_target(simulation_mixing_polygon, run_mixing_polygon_simulation(consumers, 
                                                                      sources, 
                                                                      TDF, 
                                                                      filename = "mixing_polygon_simulation", 
                                                                      its = 5000,
                                                                      threshold = 0.01)),
  # Remove few outliers
  tar_target(consumers_clean, remove_outliers(consumers, simulation_mixing_polygon)),
  
  # MixSIAR models ----
  
  # Run separately for each TDF source on different workers
  # "resid_error = FALSE" because we have some species with only one observation
  tar_target(MixSIAR_models_Post, run_MixSIAR_models(consumers_clean,
                                                     sources,
                                                     TDF,
                                                     run = list(chainLength = 200000, 
                                                                burn = 100000, 
                                                                thin = 100, 
                                                                chains = 3,
                                                                calcDIC = TRUE),
                                                     TDF_source = "Post",
                                                     resid_err = FALSE,
                                                     process_err = TRUE,
                                                     guild = "all")),
  tar_target(MixSIAR_models_McCutchan, run_MixSIAR_models(consumers_clean,
                                                          sources,
                                                          TDF,
                                                          run = list(chainLength = 200000, 
                                                                     burn = 100000, 
                                                                     thin = 100, 
                                                                     chains = 3,
                                                                     calcDIC = TRUE),
                                                          TDF_source = "McCutchan",
                                                          resid_err = FALSE,
                                                          process_err = TRUE,
                                                          guild = "all")),
  
  # Models output ----
  
  # Isospace
  tar_target(MixSIAR_isospace_Post, plot_isospace_mixsiar(MixSIAR_models_Post$mix,
                                                          MixSIAR_models_Post$source,
                                                          MixSIAR_models_Post$discr, 
                                                          filename = "isospace",
                                                          TDF_source = "Post", 
                                                          type = "png")),
  tar_target(MixSIAR_isospace_McCutchan, plot_isospace_mixsiar(MixSIAR_models_McCutchan$mix,
                                                               MixSIAR_models_McCutchan$source,
                                                               MixSIAR_models_McCutchan$discr, 
                                                               filename = "isospace",
                                                               TDF_source = "McCutchan", 
                                                               type = "png")),
  # Best models
  tar_target(MixSIAR_best_models_Post, select_best_models(MixSIAR_models_Post$models)),
  tar_target(MixSIAR_best_models_McCutchan, select_best_models(MixSIAR_models_McCutchan$models)),
  # Model comparison table
  tar_target(MixSIAR_comparison_Post, 
             make_model_comparison_table(MixSIAR_best_models_Post, 
                                         filename = "summary_MixSIAR_comparison",
                                         TDF_source = "Post")),
  tar_target(MixSIAR_comparison_McCutchan, 
             make_model_comparison_table(MixSIAR_best_models_McCutchan, 
                                         filename = "summary_MixSIAR_comparison",
                                         TDF_source = "McCutchan")),
  # Save diagnostics and statistics
  tar_target(MixSIAR_stats_diag_Post, save_MixSIAR_stats_diag(MixSIAR_best_models_Post, 
                                                              MixSIAR_models_Post$mix, 
                                                              MixSIAR_models_Post$source,
                                                              TDF_source = "Post")),
  tar_target(MixSIAR_stats_diag_McCutchan, save_MixSIAR_stats_diag(MixSIAR_best_models_McCutchan,
                                                                   MixSIAR_models_McCutchan$mix,
                                                                   MixSIAR_models_McCutchan$source,
                                                                   TDF_source = "McCutchan")),
  # Combine summary statistics
  tar_target(MixSIAR_stats_Post, make_MixSIAR_stats(MixSIAR_best_models_Post,
                                                    MixSIAR_models_Post$mix,
                                                    MixSIAR_models_Post$source)),
  tar_target(MixSIAR_stats_McCutchan, make_MixSIAR_stats(MixSIAR_best_models_McCutchan,
                                                         MixSIAR_models_McCutchan$mix,
                                                         MixSIAR_models_McCutchan$source)),
  # Save summary statistics table
  tar_target(MixSIAR_summary_Post, 
             make_MixSIAR_summary_table(MixSIAR_stats_Post, 
                                        consumers_clean, 
                                        filename = "summary_MixSIAR_proportions", 
                                        TDF_source = "Post")),
  tar_target(MixSIAR_summary_McCutchan, 
             make_MixSIAR_summary_table(MixSIAR_stats_McCutchan, 
                                        consumers_clean, 
                                        filename = "summary_MixSIAR_proportions", 
                                        TDF_source = "McCutchan")),
  # Figures
  tar_target(MixSIAR_comparison_rel_contributions,
             plot_mixing_model_comparison(stats_list = list(Post = MixSIAR_stats_Post,
                                                            McCutchan = MixSIAR_stats_McCutchan),
                                          type = "png")),
  tar_target(MixSIAR_isospace_and_rel_contributions_Post, 
             plot_isospace_and_rel_contributions(consumers_clean,
                                                 MixSIAR_models_Post, 
                                                 MixSIAR_best_models_Post,
                                                 TDF_source = "Post",
                                                 filename = "source_contribution_by_trophic_guild",
                                                 type = "png")),
  tar_target(MixSIAR_isospace_and_rel_contributions_McCutchan, 
             plot_isospace_and_rel_contributions(consumers_clean,
                                                 MixSIAR_models_McCutchan, 
                                                 MixSIAR_best_models_McCutchan,
                                                 TDF_source = "McCutchan",
                                                 filename = "source_contribution_by_trophic_guild",
                                                 type = "png")),
  tar_target(MixSIAR_rel_contributions_taxonomy_Post, 
             plot_all_rel_contribution(MixSIAR_stats_Post,
                                       TDF_source = "Post",
                                       type = "png")),
  tar_target(MixSIAR_rel_contributions_taxonomy_McCutchan, 
             plot_all_rel_contribution(MixSIAR_stats_McCutchan,
                                       TDF_source = "McCutchan",
                                       type = "png")),
  tar_target(MixSIAR_rel_contributions_vs_tl_Post, 
             plot_rel_contribution_vs_tl(MixSIAR_models_Post, 
                                         MixSIAR_best_models_Post, 
                                         TDF_source = "Post",
                                         type = "png",
                                         plot_type = "lineribbon_gradient", 
                                         .width = ppoints(40),
                                         add_line_CE_center = TRUE,
                                         add_text_CE_center = TRUE)),
  tar_target(MixSIAR_rel_contributions_vs_tl_McCutchan, 
             plot_rel_contribution_vs_tl(MixSIAR_models_McCutchan, 
                                         MixSIAR_best_models_McCutchan, 
                                         TDF_source = "McCutchan",
                                         type = "png",
                                         plot_type = "lineribbon_gradient", 
                                         .width = ppoints(40),
                                         add_line_CE_center = TRUE,
                                         add_text_CE_center = TRUE))
)
