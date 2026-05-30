#' Check sources contributions to baselines using TDFs from two meta-analyses (Post and McCutchan)
#'
#' @param sia_baselines_corrected Output of [baselines_lipid_correction()]
#' @param sources Output of [prepare_source_data()]
#' @param run The `run` argument passed to [run_model_parallel()]
#' @param combine_sources Logical; if sources should be combined
#' @param groups Only if `combine_sources = TRUE`. Named list; which sources to combine, 
#' and what names to give the new combined sources
#' @param colours Vector of colours to use for sources, plotted in alphabetic order.
#' If NULL (default), the default ggplot2 colour palette will be used
#'
#' @return A list with 5 elements:
#'  - "models": mixing models output
#'  - "diagnostics": model diagnostics
#'  - "stats": model summary statistics
#'  - "plot": patchwork of density plots of source proportions
check_baselines_diet <- function(sia_baselines_corrected, 
                                 sources,
                                 run = "test",
                                 combine_sources = FALSE,
                                 groups = NULL,
                                 colours = NULL) {
  
  # Read sources data and extract names
  sources <- readr::read_csv(sources)
  source_names <- levels(as.factor(sources$source))
  n_sources <- length(source_names)
  
  # Check combined sources
  if (combine_sources) {
    is.list(groups) ||
      cli::cli_abort("If {.arg combine_sources} is TRUE, {.arg groups} must be a named list.")
    
    # Stack list
    groups_df <- stack(groups)
    all(source_names %in% groups_df$values) ||
      cli::cli_abort(c("{.arg groups} does not include all initial sources.",
                       "Please correct source groups."))
    
    # New source names
    new_source_names <- unique(as.character(groups_df$ind))
    n_sources <- length(new_source_names)
  }
  
  # Check colours
  if (!is.null(colours)) {
    length(colours) == n_sources ||
      cli::cli_abort(
        c("{.arg colours} length incorrect.", 
          "You provided {.value {length(colours)}} colours, but there are {.value {n_sources}} sources.")
      )
  }
  
  # Prepare data for MixSIAR
  # Consumers
  consumerspath <- tempfile(pattern = "consumers", fileext = ".csv")
  readr::write_csv(sia_baselines_corrected |> 
                     select(baseline, d13C = d13C_corrected, d15N), 
                   file = consumerspath)
  # Sources
  sourcespath <- tempfile(pattern = "sources", fileext = ".csv")
  readr::write_csv(sources |> 
                     select(source, d13C, d15N), 
                   file = sourcespath)
  # TDFs from Post 2002
  TDF_Post <- data.frame(source = unique(sources$source),
                         Meand13C = 0.39,
                         SDd13C = 1.3, 
                         Meand15N = 3.4, 
                         SDd15N = 0.98)
  # TDFs from McCutchan 2003
  TDF_McCutchan <- data.frame(source = unique(sources$source),
                              Meand13C = 1.3,
                              SDd13C = 0.30,
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
                                factors = "baseline",
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
  
  # # Plot isospace Post
  # MixSIAR::plot_data_two_iso(isotopes = 1:2,
  #                            mix = mix,
  #                            source = source,
  #                            discr = tdf_Post,
  #                            plot_save_pdf = FALSE,
  #                            plot_save_png = FALSE,
  #                            return_obj = TRUE) +
  #   theme(legend.position = "inside",
  #         legend.position.inside = c(0, 0.1),
  #         legend.background = element_rect(fill = NA))
  
  # Fit models with TDFs from Post and McCutchan
  fit_baselines <- purrr::map(
    list(Post = tdf_Post,
         McCutchan = tdf_McCutchan),
    ~ run_model_parallel(run = run, 
                         mix = mix, 
                         source = source, 
                         discr = .x, 
                         model_filename = modelpath, 
                         alpha.prior = 1,
                         jags.seed = 123) # default seed
    )
  
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
  diagnostics <- purrr::map(
    fit_baselines,
    ~ MixSIAR::output_diagnostics(.x,
                                  mix = mix,
                                  source = source,
                                  output_options = mixsiar_options)
  )
  
  # Get statistics
  stats <- purrr::map(
    fit_baselines,
    ~ MixSIAR::output_stats(.x,
                            mix = mix,
                            source = source,
                            output_options = mixsiar_options) |>
      as.data.frame() |>
      tibble::rownames_to_column("parameter")
  )
  
  # Get posteriors
  df <- purrr::map(
    c("Post", "McCutchan"),
    ~ apply(fit_baselines[[.x]][["BUGSoutput"]][["sims.list"]][["p.fac1"]], 
            3, 
            function(x) as.data.frame(x) |> 
              rlang::set_names(mix$FAC[[1]]$labels) |> 
              mutate(draw = 1:3000) |> 
              tidyr::pivot_longer(cols = -draw, names_to = "taxon")) |> 
      bind_rows(.id = "source") |> 
      mutate(source = factor(source, labels = source_names),
             tdf = factor(.x))
  ) |> 
    bind_rows()
  
  # Combine sources
  if (combine_sources) {
    df <- df |>
      dplyr::left_join(groups_df,
                       by = c("source" = "values")) |> 
      dplyr::group_by(taxon, tdf, draw, ind) |> 
      dplyr::summarise(value = sum(value), 
                       .groups = "drop") |> 
      dplyr::rename("source" = "ind") |> 
      dplyr::mutate(source = as.character(source))
  }
  
  # Plot posterior
  p <- df |> 
    ggplot(aes(x = value, color = source, fill = source)) +
    geom_density(aes(y = after_stat(scaled)), alpha = 0.5) +
    facet_grid(rows = vars(tdf),
               cols = vars(taxon)) +
    labs(x = "Proportion",
         y = "Scaled Posterior Density",
         color = "Source",
         fill = "Source") +
    theme_bw()
  
  if (!is.null(colours)) {
    p <- p +
      scale_colour_manual(values = colours, 
                          aesthetics = c("colour", "fill")) 
  }
  
  return(list(models = fit_baselines,
              diagnostics = diagnostics,
              stats = stats,
              plot = p))
}
