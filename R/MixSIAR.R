#' Fitting mixing models with MixSIAR
#'
#' @param consumers_clean List of paths to consumer csv files
#' @param sources Path to sources csv file
#' @param TDF Path to TDF csv file
#' @param run The `run` argument passed to [run_model_parallel()]
#' @param alpha.prior Dirichlet prior on p.global (i.e. global source proportions).
#' Must be a numeric vector of length 1 or equal to the number of sources,
#' or a named list of numeric vectors, one for each guild (default = 1, uninformative)
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#' @param guild Trophic guilds for which models should be fitted. 
#' Either "all" to fit models to all trophic guilds, or one or more of 
#' "macrocarnivores", "microphages", "herbivores-detritivores", "planktivores",
#' "invertivores-benthic", "invertivores-sessile", "omnivores-benthic", "omnivores-pelagic"
#' 
#' @return A list of models, four for each trophic guild
run_MixSIAR_models <- function(consumers_clean,
                               sources,
                               TDF,
                               run,
                               alpha.prior = 1,
                               TDF_source = c("Post", "McCutchan"),
                               resid_err = TRUE,
                               process_err = TRUE,
                               guild = "all") {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  # Trophic guilds
  guild_names <- names(consumers_clean[[TDF_source]]) |>
    purrr::set_names()
  
  # Check guild names
  is.character(guild) ||
    cli::cli_abort(c(
      "{.var guild} must be a character vector",
      "x" = "You've supplied a {.cls {class(guild)}} vector."
    ))
  
  if (! "all" %in% guild) {
    guild_error <- guild[! guild %in% guild_names]
    length(guild_error) == 0 ||
      cli::cli_abort(c(
        "{.var guild} must be one or more of {.val {guild_names}}, 
      or {.val all} to model all trophic guilds.", 
        "x" = "{length(guild_error)} {.var guild} {?name/names} {?is/are} 
      incorrect: {.val {guild_error}}."
      ))
    guild_names <- guild_names[guild_names %in% guild] # keep only supplied guilds
  }
  
  cli::cli_inform("Setting up models")
  
  # Set up 4 models, no family random for microphages
  model_types <- purrr::map(
    guild_names,
    function(x) {
      if (x != "microphages") {
        list("full" = list(name = "full",
                           factors = c("family", "species"), 
                           fac_random = c(TRUE, TRUE), 
                           fac_nested = c(FALSE, TRUE), 
                           cont_effects = "log_tl"),
             "tl" = list(name = "tl",
                         factors = NULL, 
                         fac_random = NULL, 
                         fac_nested = NULL, 
                         cont_effects = "log_tl"),
             "taxonomy" = list(name = "taxonomy",
                               factors = c("family", "species"), 
                               fac_random = c(TRUE, TRUE), 
                               fac_nested = c(FALSE, TRUE), 
                               cont_effects = NULL),
             "null" = list(name = "null",
                           factors = NULL, 
                           fac_random = NULL, 
                           fac_nested = NULL, 
                           cont_effects = NULL))
      } else {
        list("full" = list(name = "full",
                           factors = "species", 
                           fac_random = TRUE, 
                           fac_nested = FALSE,
                           cont_effects = "log_tl"),
             "tl" = list(name = "tl",
                         factors = NULL, 
                         fac_random = NULL, 
                         fac_nested = NULL, 
                         cont_effects = "log_tl"),
             "taxonomy" = list(name = "taxonomy",
                               factors = "species", 
                               fac_random = TRUE, 
                               fac_nested = FALSE, 
                               cont_effects = NULL),
             "null" = list(name = "null",
                           factors = NULL, 
                           fac_random = NULL, 
                           fac_nested = NULL, 
                           cont_effects = NULL))
      }
    }
  )
  
  # Load mixture data: one for each model type
  mix <- purrr::map(
    guild_names,
    function(i) 
      purrr::map(
        model_types[[i]], 
        ~MixSIAR::load_mix_data(filename = consumers_clean[[TDF_source]][[i]], 
                                iso_names = c("d13C", "d15N"), 
                                factors = .x[["factors"]], 
                                fac_random = .x[["fac_random"]], 
                                fac_nested = .x[["fac_nested"]], 
                                cont_effects = .x[["cont_effects"]])
      )
  )
  
  # Load source data - same for all trophic guilds
  source <- MixSIAR::load_source_data(filename = sources, 
                                      source_factors = NULL,
                                      conc_dep = FALSE, 
                                      data_type = "raw", 
                                      mix[[1]][["full"]]) # mix only used to check isotope names
  
  # Load discrimination data - same for all trophic guilds
  discr <- MixSIAR::load_discr_data(filename = TDF,
                                    mix[[1]][["full"]]) # mix only used to check isotope names
  
  # Write models
  cli::cli_inform("Writing models")
  
  # Changes by model type, but same for all trophic guilds and TDF sources
  model_names <- purrr::map(
    guild_names,
    function(i) 
      purrr::map(
        model_types[[i]],
        ~ paste0("mixing_models/MixSIAR_", i, "_", .x[["name"]], ".txt")
      )
  )
  resid_err <- resid_err
  process_err <- process_err
  purrr::map(
    guild_names,
    function(i) 
      purrr::map(
        model_types[[i]],
        ~ MixSIAR::write_JAGS_model(model_names[[i]][[.x[["name"]]]], 
                                    resid_err, 
                                    process_err, 
                                    mix[[i]][[.x[["name"]]]], 
                                    source)
      )
  )
  
  # Priors
  
  # Check priors
  # Class must be numeric or list
  prior_class <- class(alpha.prior)
  prior_class %in% c("numeric", "list") ||
    cli::cli_abort(c(
      "{.var alpha.prior} must be a numeric vector or a list of numeric vectors."
    ))
  # Numeric must be of length 1 or same as the sources
  if (prior_class == "numeric") {
    length(alpha.prior) %in% c(1, source$n.sources) ||
      cli::cli_abort(c(
        "{.var alpha.prior} must be a numeric vector of length 1 or equal to the number of sources."
      ))
  } else {
    # If list, must be the same length as the guilds
    length(alpha.prior) == length(guild_names) ||
      cli::cli_abort(c(
        "{.var alpha.prior} must be a list of numeric vectors of length equal to the number of guilds."
      ))
    # Names must be the same as the guilds
    all(names(alpha.prior) == guild_names) ||
      cli::cli_abort(c(
        "{.var alpha.prior} must be a named list with names equal to the guilds."
      ))
    # Each element must be a numeric vector of length 1 or same as the sources
    all(sapply(alpha.prior, function(i) is.numeric(i))) & 
      all(sapply(alpha.prior, function(i) length(i) %in% c(1, source$n.sources))) ||
      cli::cli_abort(c(
        "{.var alpha.prior} must be a list of numeric vectors each of length 1 or equal to the number of sources."
      ))
  }
  
  cli::cli_inform("Setting priors")
  
  # Make prior list if a single vector is supplied
  if (is.numeric(alpha.prior)) {
    # Same prior for all guilds
    priors <- purrr::map(guild_names, ~ alpha.prior)
  } else {
    priors <- alpha.prior
  }
  
  # Fit models
  cli::cli_inform("Fitting models")
  
  models <- purrr::map(
    guild_names,
    function(i) 
      purrr::map(
        model_types[[i]],
        function(j) {
          cli::cli_progress_step('Fitting model: {.val {i}} - {.val {j[["name"]]}}', spinner = TRUE)
          fit <- run_model_parallel(run = run, 
                                    mix[[i]][[j[["name"]]]], 
                                    source, 
                                    discr, 
                                    model_names[[i]][[j[["name"]]]], 
                                    alpha.prior = priors[[i]],
                                    jags.seed = 123) # default seed
          cli::cli_progress_update()
          fit
        }
      )
  )

  list(models = models,
       mix = mix,
       source = source,
       discr = discr)
}

#' Make combined isospace plot
#'
#' @param mix Mixture data. List returned by [run_MixSIAR_models()]
#' @param source Source data. List returned by [run_MixSIAR_models()]
#' @param discr Discrimination data. List returned by [run_MixSIAR_models()]
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#' 
#' @return A patchwork object
plot_isospace_mixsiar <- function(mix,
                                  source,
                                  discr, 
                                  TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  isospace <- purrr::map(
    mix,
    ~ plot_data_two_iso_custom(isotopes = 1:2,
                               plot_save_pdf = FALSE,
                               plot_save_png = FALSE,
                               mix = .x[["full"]],
                               source = source, 
                               discr = discr, 
                               return_obj = TRUE,
                               point_size = 0.7, 
                               text_size = 2.5, 
                               linewidth = 0.4,
                               color_sources = "firebrick")
  )
  
  final_plot <- patchwork::wrap_plots(
    purrr::map(
      1:length(isospace), 
      function(.x) {
        guild_name <- names(isospace)[.x]
        isospace[[.x]]$data <- isospace[[.x]]$data |> 
          mutate(trophic_guild = guild_name)
        
        isospace[[.x]] +
          facet_wrap(~trophic_guild) +
          guides(color = guide_legend(label.theme = element_text(size = 6),
                                      keywidth = 0.1,
                                      keyheight = 0.1,
                                      direction = "vertical",
                                      default.unit = "line",
                                      ncol = 1)) +
          theme(strip.text = element_text(face = "bold"),
                legend.position = "right", 
                legend.box.spacing = unit(0, "lines"),
                axis.title = element_text(size = 8),
                axis.text = element_text(size = 8)) +
          ggthemes::scale_color_tableau(palette = "Tableau 20",
                                        breaks = levels(factor(mix[[.x]][["full"]]$FAC[[1]]$values)),
                                        labels = mix[[.x]][["full"]]$FAC[[1]]$labels)
        }
      ), 
    ncol = 2) + 
    patchwork::plot_layout(axes = "collect")
  
  return(final_plot)
}

#' Best MixSIAR model selection with LOO
#'
#' @param models List of MixSIAR models returned by [run_MixSIAR_models()]
#' 
#' @return A list
select_best_models <- function(models) {
  
  purrr::map(
    models,
    function(x) {
      loo_summary <- compare_models_parallel(x)
      name_best_model <- loo_summary |> 
        filter(weight > 0.5) |> 
        pull(Model)
      if (length(name_best_model) == 0) {
        competitive_models <- loo_summary |>
          filter(weight > 0.01 & (dLOOic < 3 | se_dLOOic > dLOOic))
        complexity <- tibble(Model = c("null", "tl", "taxonomy", "full")) |> 
          rowwise() |> 
          mutate(n_params = dim(x[[Model]]$BUGSoutput$sims.matrix)[2])
        name_best_model <- competitive_models |>
          left_join(complexity) |>
          filter(n_params == min(n_params)) |>
          pull(Model)
      }
      best_model <- x[[name_best_model]]
      list(loo_summary = loo_summary,
           name_best_model = name_best_model,
           best_model = best_model)
    }
  )
}

#' Summary table of MixSIAR model comparison with LOO
#'
#' @param best_models Output from [select_best_models()]
#' @param filename Name of csv file which will be saved in `/output/tables`
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#'
#' @return Path to the csv file
make_model_comparison_table <- function(best_models, 
                                        filename = "summary_MixSIAR_comparison", 
                                        TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  # Make table
  out <- purrr::map(best_models, 
                    ~ .x$loo_summary |> 
                      mutate(across(where(is.double), function(i) round(i, 2)))) |> 
    bind_rows(.id = "trophic_guild")
  
  # Save table
  path <- paste0("output/tables/", filename, "_", TDF_source, ".csv")
  readr::write_csv(out, path)
  
  return(path)
}

#' Set output options for MixSIAR
#'
#' Default only return object. See [MixSIAR::output_JAGS()] for details.
#'
#' @return A list of output options
set_output_options <- function(summary_save = FALSE,                 
                               summary_name = "output/MixSIAR/summary_statistics", 
                               sup_post = TRUE,                    
                               plot_post_save_pdf = FALSE,           
                               plot_post_name = "output/MixSIAR/posterior_density",
                               sup_pairs = TRUE,             
                               plot_pairs_save_pdf = FALSE,    
                               plot_pairs_name = "output/MixSIAR/pairs_plot",
                               sup_xy = TRUE,           
                               plot_xy_save_pdf = FALSE,
                               plot_xy_name = "output/MixSIAR/xy_plot",
                               gelman = TRUE,
                               heidel = FALSE,  
                               geweke = TRUE,   
                               diag_save = FALSE,
                               diag_name = "output/MixSIAR/diagnostics",
                               indiv_effect = FALSE,       
                               plot_post_save_png = FALSE, 
                               plot_pairs_save_png = FALSE,
                               plot_xy_save_png = FALSE,
                               diag_save_ggmcmc = FALSE,
                               sup_post_resid = TRUE,
                               sup_cont_eff = TRUE,
                               return_obj = TRUE) {
  
  list(summary_save = summary_save,                 
       summary_name = summary_name, 
       sup_post = sup_post,                    
       plot_post_save_pdf = plot_post_save_pdf,           
       plot_post_name = plot_post_name,
       sup_pairs = sup_pairs,             
       plot_pairs_save_pdf = plot_pairs_save_pdf,    
       plot_pairs_name = plot_pairs_name,
       sup_xy = sup_xy,           
       plot_xy_save_pdf = plot_xy_save_pdf,
       plot_xy_name = plot_xy_name,
       gelman = gelman,
       heidel = heidel,  
       geweke = geweke,   
       diag_save = diag_save,
       diag_name = diag_name,
       indiv_effect = indiv_effect,       
       plot_post_save_png = plot_post_save_png, 
       plot_pairs_save_png = plot_pairs_save_png,
       plot_xy_save_png = plot_xy_save_png,
       diag_save_ggmcmc = diag_save_ggmcmc,
       sup_post_resid = sup_post_resid,
       sup_cont_eff = sup_cont_eff,
       return_obj = return_obj)
}

#' Save statistics and diagnostics of best MixSIAR models in `output/MixSIAR`
#'
#' @param best_models Output from [select_best_models()]
#' @param mix Mixture data used for the best models. List returned by [run_MixSIAR_models()]
#' @param source Source data used for the best models. List returned by [run_MixSIAR_models()]
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#'
#' @return Paths to the saved txt files
save_MixSIAR_stats_diag <- function(best_models, 
                                    mix, 
                                    source,
                                    TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  # Main output folder
  folder <- "output/MixSIAR"
  
  # Save files
  out <- purrr::map(
    names(best_models),
    function(.x) {
      # Get best model name
      model_type <- best_models[[.x]]$name_best_model
      # Make paths
      path_stats <- paste0(folder, "/summary_statistics/", .x, "_", TDF_source)
      path_diag <- paste0(folder, "/diagnostics/", .x, "_", TDF_source)
      # Save
      output_JAGS_custom(best_models[[.x]]$best_model, 
                         mix[[.x]][[model_type]], 
                         source,
                         output_options = set_output_options(
                           summary_save = TRUE,
                           summary_name = path_stats,
                           diag_save = TRUE,
                           diag_name = path_diag,
                           sup_post_resid = TRUE))
      # Return paths
      return(c(paste0(path_stats, ".txt"), 
               paste0(path_diag, ".txt")))
    }
  )
  
  return(unlist(out))
}

# REQUIRES GITHUB VERSION OF MixSIAR
#' Combine statistics from the best MixSIAR models
#'
#' @param best_models Output from [select_best_models()]
#' @param mix Mixture data used for the best models. List returned by [run_MixSIAR_models()]
#' @param source Source data used for the best models. List returned by [run_MixSIAR_models()]
#'
#' @return A data frame
make_MixSIAR_stats <- function(best_models, 
                               mix, 
                               source) {
  
  # Get statistics for each trophic guild
  stats_list <- purrr::map(
    names(best_models) |> 
      purrr::set_names(),
    function(.x) {
      model_type <- best_models[[.x]]$name_best_model
      capture.output(
        stats <- MixSIAR::output_stats(best_models[[.x]]$best_model, 
                                       mix[[.x]][[model_type]], 
                                       source, 
                                       output_options = set_output_options(
                                         summary_save = FALSE,
                                         return_obj = TRUE))
        )
      stats |>
        as.data.frame() |>
        tibble::rownames_to_column("parameter")
      })
  
  # Merge into a single data frame
  stats_combined <- stats_list |> 
    bind_rows(.id = "trophic_guilds")
  
  return(stats_combined)
}

#' Save summary table with proportions for each trophic guild, family and species
#'
#' @param MixSIAR_stats Output from [make_MixSIAR_stats()]
#' @param consumers_clean List of paths to consumer csv files
#' @param filename Name of csv file which will be saved in `/output/tables`
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#'
#' @return Path to the saved csv file
make_MixSIAR_summary_table <- function(MixSIAR_stats,
                                       consumers_clean,
                                       filename = "summary_MixSIAR_proportions", 
                                       TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  # Consumers
  consumers <- purrr::map(
    consumers_clean[[TDF_source]],
    ~ readr::read_csv(.x)
  ) |> 
    bind_rows()
  
  # Make table
  out <- MixSIAR_stats |>
    dplyr::filter(startsWith(parameter, "p.")) |>
    dplyr::mutate(parameter = gsub("p\\.", "", parameter)) |>
    tidyr::separate_wider_delim(parameter, delim = ".", names = c("var", "source")) |>
    dplyr::left_join(consumers |> 
                       dplyr::select(species, family) |>
                       dplyr::distinct(),
                     by = c("var" = "species")) |>
    dplyr::mutate(species = ifelse(is.na(family), "-", var),
                  var = ifelse(var == "global", "-", var),
                  family = ifelse(is.na(family), var, family),
                  value = paste0(`50%`, "\n[", `2.5%`, ", ", `97.5%`, "]")) |>
    dplyr::arrange(trophic_guilds, family, species) |>
    dplyr::select("Trophic\nguild" = trophic_guilds, Family = family, Species = species, source, value) |>
    tidyr::pivot_wider(names_from = source, values_from = value)
  
  # Save table
  path <- paste0("output/tables/", filename, "_", TDF_source, ".csv")
  readr::write_csv(out, path)
  
  return(path)
}

#' Make combined plot with isospace and relative source contributions
#' to each trophic guild
#'
#' @param best_models Output from [select_best_models()]
#' @param mix Mixture data. List returned by [run_MixSIAR_models()]
#' @param source Source data. List returned by [run_MixSIAR_models()]
#' @param discr Discrimination data. List returned by [run_MixSIAR_models()]
#' @param combine_sources Logical; if sources should be combined
#' @param prior_list Only if `combine_sources = TRUE`. Output from [prepare_prior_list()]
#' @param groups Only if `combine_sources = TRUE`. Named list; which sources to combine, 
#' and what names to give the new combined sources
#' @param colours Vector of colours to use for sources, plotted in alphabetic order.
#' If NULL (default), the default ggplot2 colour palette will be used
#'
#' @return A patchwork object
plot_isospace_and_rel_contributions <- function(best_models,
                                                mix, 
                                                source,
                                                discr,
                                                combine_sources = FALSE,
                                                prior_list = NULL,
                                                groups = NULL,
                                                colours = NULL) {
  
  ## Checks
  
  # Get sources names and number
  n_sources <- source$n.sources
  source_names <- source$source_names
  
  # Check combined sources
  if (combine_sources) {
    is.list(groups) ||
      cli::cli_abort("If {.arg combine_sources} is TRUE, {.arg groups} must be a named list.")
    
    all(source_names %in% unlist(groups)) ||
      cli::cli_abort(c("{.arg groups} does not include all initial sources.",
                       "Please correct source groups."))
    
    # New source names
    source_names <- names(groups)
    n_sources <- length(groups)
  }
  
  # Check colours
  if (!is.null(colours)) {
    length(colours) == n_sources ||
      cli::cli_abort(
        c("{.arg colours} length incorrect.", 
          "You provided {.value {length(colours)}} colours, but there are {.value {n_sources}} sources.")
      )
  }
  
  ## Prepare data
  
  # Consumers
  consumers <- purrr::map(
    mix,
    ~ .x[[1]][["data"]]
  ) |> 
    bind_rows()
  
  # Sources
  sources <- data.frame(source = source$source_names) |>
    bind_cols(source$S_MU |> 
                as.data.frame() |> 
                dplyr::rename_with(~ paste("mean", .x, sep = "_"), starts_with("d"))) |>
    bind_cols(source$S_SIG |> 
                as.data.frame() |> 
                dplyr::rename_with(~ paste("sd", .x, sep = "_"), starts_with("d")))
  
  # TDFs
  tdf <- discr$mu |> 
    as.data.frame() |> 
    tibble::rownames_to_column("source") |> 
    bind_cols(discr$sig2)
  
  # Source + TDF
  source_tdf_corrected <- sources |>
    left_join(tdf, by = "source") |>
    mutate(mean_d13C = mean_d13C + Meand13C,
           mean_d15N = mean_d15N + Meand15N,
           sd_d13C = sqrt(sd_d13C^2 + SDd13C), # TDF SDs are variances, already squared by MixSIAR
           sd_d15N = sqrt(sd_d15N^2 + SDd15N)) |>
    select(1:5)
  
  if (combine_sources) {
    source_tdf_corrected <- source_tdf_corrected |> 
      bind_rows(data.frame(source = "Algae (combined)"))
  }
  
  ## Plot isospace + source contributions for each trophic guild
  
  # Colours for trophic guilds
  guild_colours = RColorBrewer::brewer.pal(length(mix), "Dark2")
  
  p <- purrr::map(
    1:length(mix),
    function(.x) {
      
      # Guild name
      guild <- names(mix)[.x]
      
      # Isospace
      p1 <- source_tdf_corrected |>
        ggplot(aes(x = mean_d13C, 
                   y = mean_d15N)) +
        # Add all consumers
        geom_point(data = consumers |>
                     select(mean_d13C = d13C,
                            mean_d15N = d15N),
                   colour = "grey80",
                   alpha = 0.2) +
        facet_wrap(~trophic_guild) +
        theme_bw() +
        labs(x = "&delta;<sup>13</sup>C (&permil;)",
             y = "&delta;<sup>15</sup>N (&permil;)") +
        theme(axis.title.x = ggtext::element_markdown(size = 10),
              axis.title.y = ggtext::element_markdown(size = 10), 
              strip.text = element_text(face = "bold", size = 10),
              legend.position = "bottom", 
              legend.title.position = "left",
              legend.title = element_text(hjust = 0.5),plot.margin = margin(1, 0, 1, 1, "pt"))
      
      # Adjust sources and consumers colours and shape  if combined sources
      if (combine_sources) {
        p1 <- p1 + 
          # Highlight consumers of this guild
          geom_point(data = consumers |>
                       filter(trophic_guild == guild) |>
                       rename(mean_d13C = d13C,
                              mean_d15N = d15N),
                     fill = "white",
                     alpha = 0.6,
                     shape = 21) +
          # Add sources
          geom_linerange(aes(ymin = mean_d15N - sd_d15N,
                             ymax = mean_d15N + sd_d15N,
                             colour = source),
                         na.rm = TRUE) +
          geom_linerange(aes(xmin = mean_d13C - sd_d13C,
                             xmax = mean_d13C + sd_d13C,
                             colour = source),
                         na.rm = TRUE) +
          geom_point(aes(shape = source, fill = source), 
                     size = 2.5,
                     colour = "black",
                     na.rm = TRUE) +
          scale_shape_manual(values = 21:25) + 
          scale_fill_manual(values = c("#117733", "#000000", "#999933", "#56B4E9", "#D55E00"),
                            aesthetics = c("colour", "fill")) +
          labs(fill = "Source",
               colour = "Source",
               shape = "Source")
      } else {
        p1 <- p1 + 
          # Highlight consumers of this guild
          geom_point(data = consumers |>
                       filter(trophic_guild == guild) |>
                       rename(mean_d13C = d13C,
                              mean_d15N = d15N),
                     fill = guild_colours[.x],
                     alpha = 0.6,
                     shape = 21) +
          # Add sources
          geom_linerange(aes(ymin = mean_d15N - sd_d15N,
                             ymax = mean_d15N + sd_d15N)) +
          geom_linerange(aes(xmin = mean_d13C - sd_d13C,
                             xmax = mean_d13C + sd_d13C)) +
          geom_point(aes(shape = source), 
                     size = 2.5,
                     fill = "white") +
          scale_shape_manual(values = 22:25) +
          labs(shape = "Source")
      }
      
      # Source relative contributions
      # Get posteriors
      model_type <- best_models[[.x]]$name_best_model
      if (combine_sources) {
        combined <- combine_sources_custom(best_models[[.x]]$best_model, 
                                           mix[[.x]][[model_type]], 
                                           source, 
                                           alpha.prior = prior_list[[.x]], 
                                           groups = groups,
                                           plot_prior = FALSE)
        post <- combined$post
      } else {
        post <- best_models[[.x]]$best_model$BUGSoutput$sims.matrix
      }
      
      # Make data frame with source names
      source_labels <- data.frame(source_index = 1:n_sources,
                                  source = source_names)
      
      # Wrangle posterior
      p2_data <- post |> 
        as.data.frame() |>
        dplyr::select(dplyr::starts_with("p.global")) |> 
        tidyr::pivot_longer(dplyr::everything(), names_to = "parameter") |> 
        mutate(source_index = stringr::str_extract(parameter, "(?<=\\[).*?(?=\\])"), 
               source_index = as.numeric(source_index)) |> 
        left_join(source_labels, by = "source_index")
      
      # Plot
      p2 <- p2_data |> 
        ggplot(aes(x = source, y = value)) +
        labs(y = "Relative contribution") +
        theme_classic() +
        theme(axis.title.y = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = margin(1, 1, 1, 0, "pt")) + 
        coord_cartesian(clip = "off")
      
      # Add sources
      if (combine_sources) {
        p2 <- p2 +
          ggdist::stat_pointinterval(aes(colour = source, fill = source, shape = source),
                                     .width = c(0.5, 0.95), 
                                     point_colour = "black",
                                     stroke = 0.5,
                                     show.legend = FALSE) +
          scale_shape_manual(values = c(21, 22, 24))
      } else {
        p2 <- p2 +
          ggdist::stat_pointinterval(aes(shape = source),
                                     .width = c(0.5, 0.95), 
                                     fill = "white",
                                     show.legend = FALSE) +
          scale_shape_manual(values = 22:25)
      }
      
      if (combine_sources & !is.null(colours)) {
        p2 <- p2 +
          scale_colour_manual(values = colours, aesthetics = c("colour", "fill"))
      }
      
      # Remove title x axis p1
      if (.x %in% 1:6) p1 <- p1 + theme(axis.title.x = element_blank())
      
      # Adjust y axes position and limits and combine plots
      if (.x %in% seq(1, length(mix), by = 2)) {
        p1 <- p1 + scale_y_continuous(position = "right")
        p2 <- p2 + scale_y_continuous(limits = c(0, 1), 
                                      expand = c(0, 0), 
                                      position = "left")
        p2 + p1 + patchwork::plot_layout(widths = c(0.2, 0.8))
      } else {
        p1 <- p1 + theme(axis.title.y = element_blank())
        p2 <- p2 + scale_y_continuous(limits = c(0, 1), 
                                      expand = c(0, 0), 
                                      position = "right")
        p1 + p2 + patchwork::plot_layout(widths = c(0.8, 0.2))
      }
    }
  )
  
  ## Combine all plots
  p <- patchwork::wrap_plots(p, ncol = 2) + 
    patchwork::plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom')
  
  return(p)
}

#' Plot relative source contributions for each species, family and trophic guild
#'
#' @param best_models Output from [select_best_models()]
#' @param mix Mixture data. List returned by [run_MixSIAR_models()]
#' @param source Source data. List returned by [run_MixSIAR_models()]
#' @param combine_sources Logical; if sources should be combined
#' @param prior_list Only if `combine_sources = TRUE`. Output from [prepare_prior_list()]
#' @param groups Only if `combine_sources = TRUE`. Named list; which sources to combine, 
#' and what names to give the new combined sources
#' @param colours Vector of colours to use for sources, plotted in alphabetic order.
#' If NULL (default), the default ggplot2 colour palette will be used
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#' @param filename Name of the file to save
#' @param filetype Type of the file to save (default to png)
#' 
#' @return A list of paths to saved plots
plot_all_rel_contribution <- function(best_models,
                                      mix,
                                      source,
                                      combine_sources = FALSE,
                                      prior_list = NULL,
                                      groups = NULL,
                                      colours = NULL,
                                      TDF_source = c("Post", "McCutchan"),
                                      filename = "rel_contributions_",
                                      filetype = "png") {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  # Get sources names and number
  n_sources <- source$n.sources
  source_names <- source$source_names
  
  # Check combined sources
  if (combine_sources) {
    is.list(groups) ||
      cli::cli_abort("If {.arg combine_sources} is TRUE, {.arg groups} must be a named list.")
    
    all(source_names %in% unlist(groups)) ||
      cli::cli_abort(c("{.arg groups} does not include all initial sources.",
                       "Please correct source groups."))
    
    # New source names
    source_names <- names(groups)
    n_sources <- length(groups)
  }
  
  # Check colours
  if (!is.null(colours)) {
    length(colours) == n_sources ||
      cli::cli_abort(
        c("{.arg colours} length incorrect.", 
          "You provided {.value {length(colours)}} colours, but there are {.value {n_sources}} sources.")
      )
  }
  
  # Plot point interval for each species, family and trophic guild
  plots <- purrr::map(
    names(mix) |>
      rlang::set_names(),
    function(.x) {
      # Get posteriors and number and names of sources
      model_type <- best_models[[.x]]$name_best_model
      if (combine_sources) {
        combined <- combine_sources_custom(best_models[[.x]]$best_model, 
                                           mix[[.x]][[model_type]], 
                                           source, 
                                           alpha.prior = prior_list[[.x]], 
                                           groups = groups,
                                           plot_prior = FALSE)
        post <- combined$post
      } else {
        post <- best_models[[.x]]$best_model$BUGSoutput$sims.matrix
      }
      
      # Make data frame with labels of each factor variable
      fac1_labels <- mix[[.x]][[model_type]]$FAC[[1]]$labels
      fac_labels <- data.frame(var = rep("fac1", length(fac1_labels)),
                               fac_name = mix[[.x]][[model_type]]$FAC[[1]]$name,
                               fac_index = 1:length(fac1_labels),
                               fac_label = fac1_labels)
      if (length(mix[[.x]][[model_type]]$factors) == 2) {
        fac2_labels <- mix[[.x]][[model_type]]$FAC[[2]]$labels
        fac_labels <- fac_labels |> 
          dplyr::bind_rows(data.frame(var = rep("fac2", length(fac2_labels)),
                                      fac_name = mix[[.x]][[model_type]]$FAC[[2]]$name,
                                      fac_index = 1:length(fac2_labels),
                                      fac_label = fac2_labels))
      }
      
      # Make data frame with source names
      source_labels <- data.frame(source_index = 1:n_sources,
                                  source = source_names)
      
      # Wrangle posteriors
      plot_data <- post |> 
        as.data.frame() |>
        dplyr::select(dplyr::starts_with("p.") & !starts_with("p.ind")) |> 
        tidyr::pivot_longer(dplyr::everything(), names_to = "parameter") |> 
        dplyr::mutate(parameter = gsub("p\\.", "", parameter)) |>
        tidyr::separate_wider_delim(parameter, delim = "[", names = c("var", "index")) |> 
        dplyr::mutate(index = gsub("]", "", index),
                      index = ifelse(!stringr::str_detect(index, ","), paste0("NA,", index), index)) |> 
        tidyr::separate_wider_delim(index, delim = ",", names = c("fac_index", "source_index")) |> 
        dplyr::mutate(fac_index = suppressWarnings(as.numeric(fac_index)),
                      source_index = suppressWarnings(as.numeric(source_index))) |> 
        dplyr::left_join(fac_labels, by = c("var", "fac_index")) |> 
        dplyr::left_join(source_labels, by = "source_index") |> 
        dplyr::mutate(fac_name = ifelse(is.na(fac_name), "", fac_name),
                      var_type = factor(fac_name, 
                                        levels = c("", "family", "species")),
                      label = dplyr::case_when(fac_name == "species" ~ glue::glue("<i>{fac_label}</i>"), 
                                               var == "global" ~ var,
                                               .default = fac_label))
      
      # Plot
      p <- plot_data |> 
        ggplot(aes(x = value, y = label, colour = source)) +
        ggdist::stat_pointinterval(.width = c(0.5, 0.95),
                                   interval_size_range = c(0.4, 0.8),
                                   point_size = 1.5,
                                   position = position_dodge(width = 0.7)) +
        facet_grid(rows = vars(var_type),
                   scales = "free_y",
                   space = "free") +
        labs(x = "Relative contribution",
             title = .x,
             colour = "Source") + 
        theme_bw() +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_text(size = 9),
              axis.text.y = ggtext::element_markdown(size = 7),
              axis.text.x = element_text(size = 7),
              plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
              legend.position = "bottom", 
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 8),
              legend.box.spacing = unit(0, "lines"))
      
      if (!is.null(colours)) {
        p <- p +
          scale_colour_manual(values = colours)
      }
      
      return(p)
    }
  )
  
  # Save plots and return paths
  paths <- names(plots) |>
    purrr::set_names() |>
    purrr::map(~ paste0("output/figures/", filename, "_", .x, "_", TDF_source, ".", filetype))
  
  purrr::map(names(plots),
             ~ ggsave(paths[[.x]], 
                      plots[[.x]], 
                      width = 12, 
                      height = ifelse(.x %in% c("invertivores-benthic", "macrocarnivores"), 24, 18), 
                      units = "cm"))
  
  return(paths)
}

#' Plot total length effect on relative source contributions by trophic guild
#'
#' @param models Output from [run_MixSIAR_models()]
#' @param best_models Output from [select_best_models()]
#' @param plot_type One of "lineribbon" (Default), "spaghetti", or "lineribbon_gradient"
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used for trophic position estimation
#' @param filename Name of the file to save
#' @param filetype Type of the file to save (default to png)
#' @param colours Vector of colours to use for sources, plotted in alphabetic order.
#' If NULL (default), the default ggplot2 colour palette will be used
#' @param ... Arguments passed on to [plot_MixSIAR_continuous()]
#'
#' @return Path to the saved plot
plot_rel_contribution_vs_tl <- function(models, 
                                        best_models, 
                                        plot_type = "lineribbon",
                                        TDF_source = c("Post", "McCutchan"), 
                                        filename = "rel_contributions_vs_tl",
                                        filetype = "png",
                                        colours = NULL,
                                        ...) {
  
  # Plot total length effect for guilds with best model "full" or "tl"
  plots <- purrr::map(
    names(best_models) |> 
      rlang::set_names(),
    function(.x) {
      name_best_model <- best_models[[.x]]$name_best_model
      if (name_best_model %in% c("full", "tl")) {
        p <- plot_MixSIAR_continuous(jags.1 = best_models[[.x]]$best_model, 
                                     mix = models$mix[[.x]][[name_best_model]], 
                                     source = models$source, 
                                     plot_type = plot_type,
                                     ...) + 
          labs(x = "Total length (cm) [natural-log]") +
          theme(panel.grid = element_blank())
        
        p$data <- p$data |> 
          mutate(trophic_guild = .x)
        
        p <- p +
          facet_wrap(~trophic_guild) +
          theme(strip.text = element_text(face = "bold"))
        
        if (!is.null(colours)) {
          p <- p +
            scale_colour_manual(name = "Source",
                                values = colours,
                                aesthetics = c("colour", "fill"))
        }
        
        return(p)
      }
    }
  )
  
  # Combine plots
  final_plot <- plots |> 
    purrr::keep(purrr::negate(is.null)) |> 
    patchwork::wrap_plots(ncol = 2) + 
    patchwork::plot_layout(guides = "collect", 
                           axis_titles = "collect") & 
    theme(legend.position = "bottom")
  
  if (plot_type == "lineribbon_gradient") {
    final_plot <- final_plot &
      guides(fill_ramp = ggdist::guide_rampbar(title = "CI", 
                                               theme = theme(legend.title = element_text(vjust = 0.8))))
    
  }
  
  # Save plot
  path <- paste0("output/figures/", filename, "_", TDF_source, ".", filetype)
  ggsave(path, final_plot, width = 18, height = 22, units = "cm")
  
  return(path)
}
