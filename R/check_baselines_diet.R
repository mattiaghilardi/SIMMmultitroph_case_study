#' Check sources contributions to potential baselines using two TDFs (Post and McCutchan)
#'
#' @param sources Output of [prepare_source_data()]
#' @param run The `run` argument passed to [run_model_parallel()]
#'
#' @return A list with 5 elements:
#'  - "models": mixing models output
#'  - "diagnostics": model diagnostics
#'  - "stats": model summary statistics
#'  - "plot": patchwork of density plots of source proportions
#'  - "baselines": names of selected baselines
check_baselines_diet <- function(sources, run = "test") {
  
  # Inverts data
  sia_inverts <- sources$invertebrates$data
  
  # Summarise data
  summary_inverts <- sia_inverts |>
    group_by(taxon) |>
    summarise(n = n(),
              mean_d13C = mean(d13C, na.rm = TRUE), 
              sd_d13C = sd(d13C, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE)) |>
    ungroup()
  
  # Potential baselines (primary consumers)
  potential_baselines <- c("Bivalvia", "Porifera", "Gastropoda")
  consumers <- sia_inverts |> 
    filter(taxon %in% potential_baselines)
  
  # Prepare data for MixSIAR
  consumerspath <- tempfile(pattern = "consumers", fileext = ".csv")
  readr::write_csv(consumers, 
                   file = consumerspath)
  
  sourcespath <- tempfile(pattern = "sources", fileext = ".csv")
  readr::write_csv(sources$sources$data, 
                   file = sourcespath)
  # TDFs from Post 2002
  TDF_Post <- data.frame(source = c("Macroalgae", "Cyanobacteria", "POM"),
                         Meand13C = 0.39,
                         SDd13C = 1.3, 
                         Meand15N = 3.4, 
                         SDd15N = 0.98)
  # TDFs from McCutchan 2003
  TDF_McCutchan <- data.frame(source = c("Macroalgae", "Cyanobacteria", "POM"),
                              Meand13C = 1.3,
                              SDd13C = 0.3,
                              Meand15N = 2.9, 
                              SDd15N = 0.32)
  tdfpath_Post <- tempfile(pattern = "tdf_Post", fileext = ".csv")
  readr::write_csv(TDF_Post,
                   file = tdfpath_Post)
  tdfpath_McCutchan <- tempfile(pattern = "tdf_McCutchan", fileext = ".csv")
  readr::write_csv(TDF_McCutchan,
                   file = tdfpath_McCutchan)
  
  # Load mixture data
  mix <- MixSIAR::load_mix_data(consumerspath, 
                                iso_names = c("d13C", "d15N"), 
                                factors = "taxon",
                                fac_random = FALSE,
                                fac_nested = NULL,
                                cont_effects = NULL)
  # Load source data
  source <- MixSIAR::load_source_data(sourcespath, 
                                      data_type = "raw", 
                                      mix = mix,
                                      source_factors = NULL,
                                      conc_dep = FALSE)
  # Load discrimination data
  tdf_Post <- MixSIAR::load_discr_data(tdfpath_Post,
                                       mix = mix)
  tdf_McCutchan <- MixSIAR::load_discr_data(tdfpath_McCutchan,
                                            mix = mix)
  # Write model
  # Residual only error because filter feeders don't sample individual prey (Stock et al. 2018)
  modelpath <- tempfile(pattern = "model", fileext = ".txt")
  MixSIAR::write_JAGS_model(modelpath, 
                            resid_err = TRUE,
                            process_err = FALSE, 
                            mix = mix, 
                            source = source)
  
  # MixSIAR::plot_data_two_iso(isotopes = 1:2,
  #                            mix = mix,
  #                            source = source, 
  #                            discr = tdf_Post, 
  #                            plot_save_pdf = FALSE,
  #                            plot_save_png = FALSE,
  #                            return_obj = TRUE)
  
  # Fit models with TDFs from Post and McCutchan
  fit_inverts <- purrr::map(
    list(Post = tdf_Post,
         McCutchan = tdf_McCutchan),
    ~ run_model_parallel(run = run, 
                         mix = mix, 
                         source = source, 
                         discr = .x, 
                         model_filename = modelpath, 
                         alpha.prior = 1,
                         seed = 123))
  
  # Set output options to avoid saving summary, diagnostics and plots
  mixsiar_options <- list(summary_save = FALSE,
                          summary_name = "summary_statistics",
                          sup_post = TRUE,
                          plot_post_save_pdf = FALSE,
                          plot_post_name = "posterior_density",
                          sup_pairs = TRUE,
                          plot_pairs_save_pdf = FALSE,
                          plot_pairs_name = "pairs_plot",
                          sup_xy = TRUE,
                          plot_xy_save_pdf = FALSE,
                          plot_xy_name = "xy_plot",
                          gelman = TRUE,
                          heidel = FALSE,
                          geweke = TRUE,
                          diag_save = FALSE,
                          diag_name = "diagnostics",
                          indiv_effect = FALSE,       
                          plot_post_save_png = FALSE, 
                          plot_pairs_save_png = FALSE,
                          plot_xy_save_png = FALSE,
                          diag_save_ggmcmc = FALSE,
                          return_obj = TRUE)
  
  # Get diagnostics
  diagnostics <- purrr::map(fit_inverts,
                            ~ MixSIAR::output_diagnostics(.x,
                                                          mix = mix,
                                                          source = source,
                                                          output_options = mixsiar_options))
  
  # Get statistics
  stats <- purrr::map(fit_inverts,
                      ~ MixSIAR::output_stats(.x,
                                              mix = mix,
                                              source = source,
                                              output_options = mixsiar_options) |>
                        as.data.frame() |>
                        tibble::rownames_to_column("parameter"))
  
  # Plot posteriors
  source_names <- source$source_names
  p <- purrr::map(
    c("Post", "McCutchan"),
    ~ apply(fit_inverts[[.x]][["BUGSoutput"]][["sims.list"]][["p.fac1"]], 
            3, 
            function(x) as.data.frame(x) |> 
              rlang::set_names(mix$FAC[[1]]$labels)) |> 
      bind_rows(.id = "source") |> 
      tidyr::pivot_longer(cols = -source, names_to = "taxon") |> 
      mutate(source = factor(source, labels = source_names),
             tdf = factor(.x))
  ) |> 
    bind_rows() |> 
    ggplot(aes(x = value, color = source, fill = source)) +
    geom_density(aes(y = after_stat(scaled)), alpha = 0.3) +
    facet_grid(rows = vars(tdf),
               cols = vars(taxon)) +
    xlab("Proportion") +
    ylab("Scaled Posterior Density") +
    scale_fill_viridis_d(aesthetics = c("color", "fill")) +
    theme_bw() +
    theme(legend.title = element_blank())
  
  ggsave(here::here("output", "figures", "mixsiar_baselines.png"), 
         p,
         width = 18, height = 12, units = "cm")
  
  # Use Bivalvia (more pelagic) and Gastropoda (benthic)
  selected_baselines <- c("Bivalvia", "Gastropoda")
  
  return(list(models = fit_inverts,
              diagnostics = diagnostics,
              stats = stats,
              plot = p,
              baselines = selected_baselines))
}
