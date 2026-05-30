library(targets)
library(crew)

# Set up controller to launch 2 workers on local processes
controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 2
)

# Set targets options and required packages
tar_option_set(packages = c("readr", "tibble", "dplyr", "tidyr", "stringr", "purrr", "rlang",
                            "ggplot2", "ggtext", "patchwork", "ggdist", "ggblend", "lemon", 
                            "ggplotify", "cowplot", "GGally", "ggthemes", "splancs", "sp",
                            "tRophicPosition", "MixSIAR", "R2jags", "brms", "tidybayes", 
                            "parallel", "cmdstanr", "loo", "cli"),
               controller = controller,
               memory = "transient", 
               garbage_collection = TRUE,
               seed = 25,
               storage = "worker", 
               retrieval = "worker")

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

  # Load data ---------------------------------------------------------------
  
  # Consumers
  tar_target(file_sia_fish, 
             "data/sia_fish.csv", 
             format = "file"),
  tar_target(sia_fish, 
             read_csv(file_sia_fish)),
  # Sources
  tar_target(file_sia_sources, 
             "data/sia_sources.csv", 
             format = "file"),
  tar_target(sia_sources, 
             read_csv(file_sia_sources)),
  # Baselines
  tar_target(file_sia_baselines, 
             "data/sia_baselines.csv", 
             format = "file"),
  tar_target(sia_baselines, 
             read_csv(file_sia_baselines)),
  # Trophic guilds
  tar_target(file_guilds, 
             "data/trophic_guilds.csv", 
             format = "file"),
  tar_target(trophic_guilds, 
             read_csv(file_guilds)),
  # Priors
  tar_target(file_priors,
             "data/priors.csv",
             format = "file"),
  tar_target(priors, 
             read_csv(file_priors)),

  # Lipid correction --------------------------------------------------------
  
  # Perform lipid correction for baselines 
  tar_target(sia_baselines_corrected,
             baselines_lipid_correction(sia_baselines)),
  
  # Perform lipid correction for consumers
  tar_target(sia_fish_corrected,
             fish_lipid_correction(sia_fish)),

  # Plot isospace -----------------------------------------------------------
  
  # Plot initial isospace
  tar_target(initial_isospace, 
             plot_initial_isospace(sia_sources,
                                   sia_baselines_corrected,
                                   sia_fish_corrected)),
  
  # Initial check of source and consumer data -------------------------------
  
  # Check overlap between source isotopic signatures
  tar_target(source_overlap, 
             check_source_overlap(sia_sources)),
  
  # Check fish isotopic signatures across years
  tar_target(fish_isotope_check, 
             check_isotopes_across_years(sia_fish_corrected)),
  
  # Prepare source and TDF data for MixSIAR ---------------------------------
  
  # Prepare source data for MixSIAR
  tar_target(sources, 
             prepare_source_data(sia_sources),
             format = "file"),
  
  # Prepare TDF data for MixSIAR
  tar_target(TDF, 
             prepare_TDF_data(sources),
             format = "file"),
  
  # Source groups and colours -----------------------------------------------
  
  tar_target(source_groups, 
             list("Algae" = c("Red algae", "Green-Brown algae"),
                  "Cyanobacteria" = "Cyanobacteria",
                  "POM" = "POM")),
  
  tar_target(source_colours,
             c("#117733", "#000000", "#56B4E9")),
  
  # Check baselines diet ----------------------------------------------------
  
  # Check source contribution to baselines
  tar_target(baselines_diet, 
             check_baselines_diet(sia_baselines_corrected,
                                  sources, 
                                  run = "normal",
                                  combine_sources = TRUE,
                                  groups = source_groups,
                                  colours = source_colours)),
  
  # Estimate trophic position -----------------------------------------------
  
  # NOTE: results of this step and consequently of the rest of the pipeline
  # will differ slightly each time the code runs
  # Modifications to the source code of `tRophicPosition::TPmodel()` 
  # are required to obtain reproducible results
  tar_target(TP, 
             estimate_TP(sia_fish_corrected, 
                         sia_baselines_corrected)),
  
  # Plot correlation between TPs
  tar_target(TP_corr_plot,
             plot_TP_correlation(TP)),
  
  # Plot species-specific TPs
  tar_target(TP_plots,
             plot_TP(TP)),
  
  # Prepare consumer data for MixSIAR ---------------------------------------

  # Add trophic guilds and rescale isotopic values to TP = 1
  tar_target(consumers, 
             prepare_consumer_data(sia_fish_corrected, 
                                   trophic_guilds, 
                                   TP)),
  
  # Plot rescaled isospace
  tar_target(rescaled_isospace,
             plot_rescaled_isospace(sources,
                                    consumers)),
  
  # Mixing polygon simulations
  tar_target(simulation_mixing_polygon, 
             run_mixing_polygon_simulation(consumers, 
                                           sources, 
                                           TDF, 
                                           its = 5000,
                                           threshold = 0.01)),
  
  # Plot probability in mixing polygon
  tar_target(prob_in_mixing_polygon,
             plot_prob_in_mixing_polygon(simulation_mixing_polygon)),
  
  # Identify consumer outliers
  tar_target(consumer_outliers, 
             identify_outliers(simulation_mixing_polygon)),
  
  # Plot consumer outliers
  tar_target(consumer_outliers_plot, 
             plot_outliers(consumers,
                           consumer_outliers)),
  
  # Remove outliers
  tar_target(consumers_clean, 
             remove_outliers(consumers, 
                             consumer_outliers)),
  
  # Fit MixSIAR models ------------------------------------------------------
  
  # Prepare list of Dirichlet priors on p.global
  tar_target(prior_list,
             prepare_prior_list(priors)),
  
  # Plot Dirichlet priors
  tar_target(prior_plot,
             plot_priors(sources,
                         prior_list)),
  
  # Run separately for each TDF source on different workers
  # "resid_error = FALSE" because we have some species with only one observation
  # and a nested hierarchical structure
  tar_target(MixSIAR_models_Post, 
             run_MixSIAR_models(consumers_clean,
                                sources,
                                TDF,
                                run = "long",
                                alpha.prior = prior_list,
                                TDF_source = "Post",
                                resid_err = FALSE,
                                process_err = TRUE,
                                guild = "all")),
  tar_target(MixSIAR_models_McCutchan, 
             run_MixSIAR_models(consumers_clean,
                                sources,
                                TDF,
                                run = "long",
                                alpha.prior = prior_list,
                                TDF_source = "McCutchan",
                                resid_err = FALSE,
                                process_err = TRUE,
                                guild = "all")),
  # ~15h

  # Compute polygon area ----------------------------------------------------
  
  # Get standardised surface area of mixing polygon
  # This can be computed once as it is identical for all models
  tar_target(polygon_area,
             MixSIAR::calc_area(source = MixSIAR_models_Post$source, 
                                mix = MixSIAR_models_Post$mix$`herbivores-detritivores`$full, 
                                discr = MixSIAR_models_Post$discr)),
  
  # Models output ----
  
  # Isospace
  tar_target(MixSIAR_isospace_Post, 
             plot_isospace_mixsiar(MixSIAR_models_Post$mix,
                                   MixSIAR_models_Post$source,
                                   MixSIAR_models_Post$discr,
                                   TDF_source = "Post")),
  tar_target(MixSIAR_isospace_McCutchan, 
             plot_isospace_mixsiar(MixSIAR_models_McCutchan$mix,
                                   MixSIAR_models_McCutchan$source,
                                   MixSIAR_models_McCutchan$discr, 
                                   TDF_source = "McCutchan")),
  
  # Best models
  tar_target(MixSIAR_best_models_Post, 
             select_best_models(MixSIAR_models_Post$models)),
  tar_target(MixSIAR_best_models_McCutchan, 
             select_best_models(MixSIAR_models_McCutchan$models)),
  
  # Model comparison table
  tar_target(MixSIAR_comparison_Post, 
             make_model_comparison_table(MixSIAR_best_models_Post, 
                                         filename = "summary_MixSIAR_comparison",
                                         TDF_source = "Post"),
             format = "file"),
  tar_target(MixSIAR_comparison_McCutchan, 
             make_model_comparison_table(MixSIAR_best_models_McCutchan, 
                                         filename = "summary_MixSIAR_comparison",
                                         TDF_source = "McCutchan"),
             format = "file"),
  
  # Save diagnostics and statistics
  tar_target(MixSIAR_stats_diag_Post, 
             save_MixSIAR_stats_diag(MixSIAR_best_models_Post, 
                                     MixSIAR_models_Post$mix, 
                                     MixSIAR_models_Post$source,
                                     TDF_source = "Post")),
  tar_target(MixSIAR_stats_diag_McCutchan, 
             save_MixSIAR_stats_diag(MixSIAR_best_models_McCutchan,
                                     MixSIAR_models_McCutchan$mix,
                                     MixSIAR_models_McCutchan$source,
                                     TDF_source = "McCutchan")),
  
  # Combine summary statistics
  tar_target(MixSIAR_stats_Post, 
             make_MixSIAR_stats(MixSIAR_best_models_Post,
                                MixSIAR_models_Post$mix,
                                MixSIAR_models_Post$source)),
  tar_target(MixSIAR_stats_McCutchan, 
             make_MixSIAR_stats(MixSIAR_best_models_McCutchan,
                                MixSIAR_models_McCutchan$mix,
                                MixSIAR_models_McCutchan$source)),
  
  # Save summary statistics table
  tar_target(MixSIAR_summary_Post, 
             make_MixSIAR_summary_table(MixSIAR_stats_Post, 
                                        consumers_clean, 
                                        filename = "summary_MixSIAR_proportions", 
                                        TDF_source = "Post"),
             format = "file"),
  tar_target(MixSIAR_summary_McCutchan, 
             make_MixSIAR_summary_table(MixSIAR_stats_McCutchan, 
                                        consumers_clean, 
                                        filename = "summary_MixSIAR_proportions", 
                                        TDF_source = "McCutchan"),
             format = "file"),
  
  # Figures
  tar_target(MixSIAR_comparison_rel_contributions,
             plot_mixing_model_comparison(stats_list = list(Post = MixSIAR_stats_Post,
                                                            McCutchan = MixSIAR_stats_McCutchan),
                                          combine_sources = TRUE,
                                          groups = source_groups,
                                          colours = source_colours)),
  tar_target(MixSIAR_isospace_and_rel_contributions_Post, 
             plot_isospace_and_rel_contributions(MixSIAR_best_models_Post,
                                                 MixSIAR_models_Post$mix, 
                                                 MixSIAR_models_Post$source, 
                                                 MixSIAR_models_Post$discr, 
                                                 combine_sources = TRUE,
                                                 prior_list = prior_list,
                                                 groups = source_groups,
                                                 colours = source_colours)),
  tar_target(MixSIAR_isospace_and_rel_contributions_McCutchan, 
             plot_isospace_and_rel_contributions(MixSIAR_best_models_McCutchan,
                                                 MixSIAR_models_McCutchan$mix, 
                                                 MixSIAR_models_McCutchan$source, 
                                                 MixSIAR_models_McCutchan$discr,
                                                 combine_sources = TRUE,
                                                 prior_list = prior_list,
                                                 groups = source_groups,
                                                 colours = source_colours)),
  tar_target(MixSIAR_rel_contributions_taxonomy_Post, 
             plot_all_rel_contribution(MixSIAR_best_models_Post,
                                       MixSIAR_models_Post$mix,
                                       MixSIAR_models_Post$source,
                                       combine_sources = TRUE,
                                       prior_list = prior_list,
                                       groups = source_groups,
                                       colours = source_colours, 
                                       TDF_source = "Post", 
                                       filename = "rel_contributions",
                                       filetype = "png")),
  tar_target(MixSIAR_rel_contributions_taxonomy_McCutchan, 
             plot_all_rel_contribution(MixSIAR_best_models_McCutchan,
                                       MixSIAR_models_McCutchan$mix,
                                       MixSIAR_models_McCutchan$source,
                                       combine_sources = TRUE,
                                       prior_list = prior_list,
                                       groups = source_groups,
                                       colours = source_colours, 
                                       TDF_source = "McCutchan", 
                                       filename = "rel_contributions",
                                       filetype = "png")),
  tar_target(MixSIAR_rel_contributions_vs_tl_Post,
             plot_rel_contribution_vs_tl(MixSIAR_models_Post,
                                         MixSIAR_best_models_Post,
                                         TDF_source = "Post", 
                                         filename = "rel_contributions_vs_tl",
                                         filetype = "png",
                                         plot_type = "lineribbon",
                                         .width = c(0.5, 0.95),
                                         alpha = 0.6,
                                         add_line_CE_center = TRUE,
                                         add_text_CE_center = FALSE,
                                         resolution = 50,
                                         combine_sources = TRUE,
                                         groups = source_groups,
                                         colours = source_colours),
             format = "file"),
  tar_target(MixSIAR_rel_contributions_vs_tl_McCutchan,
             plot_rel_contribution_vs_tl(MixSIAR_models_McCutchan,
                                         MixSIAR_best_models_McCutchan,
                                         TDF_source = "McCutchan", 
                                         filename = "rel_contributions_vs_tl",
                                         filetype = "png",
                                         plot_type = "lineribbon",
                                         .width = c(0.5, 0.95),
                                         alpha = 0.6,
                                         add_line_CE_center = TRUE,
                                         add_text_CE_center = FALSE,
                                         resolution = 50,
                                         combine_sources = TRUE,
                                         groups = source_groups,
                                         colours = source_colours),
             format = "file")
)
