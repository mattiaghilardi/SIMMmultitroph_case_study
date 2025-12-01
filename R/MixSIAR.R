#' Fitting mixing models with MixSIAR
#'
#' @param consumers_clean List of paths to consumer csv files
#' @param sources Paths to sources csv file
#' @param TDF Paths to TDF csv file
#' @param run The `run` argument passed to [run_model_parallel()]
#' @param alpha.prior Dirichlet prior on p.global (i.e. global source proportions) (default = 1, uninformative)
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#' @param guild Trophic guilds for which models should be fitted. Either "all"
#' to fit models to all trophic guilds, or one or more of 
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
  guild_names <- names(consumers_clean$path[[TDF_source]]) |>
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
        ~MixSIAR::load_mix_data(filename = consumers_clean$path[[TDF_source]][[i]], 
                                iso_names = c("d13C", "d15N"), 
                                factors = .x[["factors"]], 
                                fac_random = .x[["fac_random"]], 
                                fac_nested = .x[["fac_nested"]], 
                                cont_effects = .x[["cont_effects"]])
      )
  )
  
  # Load source data - same for all trophic guilds
  source <- MixSIAR::load_source_data(filename = sources$sources$path, 
                                      source_factors = NULL,
                                      conc_dep = FALSE, 
                                      data_type = "raw", 
                                      mix[[1]][["full"]]) # mix only used to check isotope names
  
  # Load discrimination data - same for all trophic guilds
  discr <- MixSIAR::load_discr_data(filename = TDF$path,
                                    mix[[1]][["full"]]) # mix only used to check isotope names
  
  # Write models
  cli::cli_inform("Writing models")
  
  # Changes by model type, but same for all trophic guilds and TDF sources
  model_names <- purrr::map(
    guild_names,
    function(i) 
      purrr::map(
        model_types[[i]],
        ~ here::here("mixing_models", paste0("MixSIAR_", i, "_", .x[["name"]], ".txt"))
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
                                    alpha.prior = alpha.prior,
                                    seed = 123)
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

#' Save combined isospace plot
#'
#' @param mix Mixture data. List returned by [run_MixSIAR_models()]
#' @param source Source data. List returned by [run_MixSIAR_models()]
#' @param discr Discrimination data. List returned by [run_MixSIAR_models()]
#' @param filename Name of file to save
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#' @param height Plot height
#' @param width Plot width
#' @param units Units for height and width
#' @param type File type (default to png)
#' 
#' @return Path to the saved file
plot_isospace_mixsiar <- function(mix,
                                  source,
                                  discr,
                                  filename = "isospace", 
                                  TDF_source = c("Post", "McCutchan"),
                                  width = 18, 
                                  height = 20, 
                                  units = "cm", 
                                  type = "png") {
  
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
  
  patchwork::wrap_plots(
    purrr::map(
      1:length(isospace), 
      ~ isospace[[.x]] +
        ggtitle(names(isospace)[.x]) +
        guides(color = guide_legend(label.theme = element_text(size = 6),
                                    keywidth = 0.1,
                                    keyheight = 0.1,
                                    direction = "vertical",
                                    default.unit = "line",
                                    ncol = 1)) +
        theme(plot.title = element_text(size = 8,
                                        face = "bold",
                                        hjust = 0.5),
              legend.position = "right", 
              legend.box.spacing = unit(0, "lines"),
              axis.title = element_text(size = 8),
              axis.text = element_text(size = 8)) +
      scale_color_viridis_d(breaks = levels(factor(mix[[.x]][["full"]]$FAC[[1]]$values)),
                            labels = mix[[.x]][["full"]]$FAC[[1]]$labels)
    ), 
    ncol = 2)
  
  path <- here::here("output", "figures", paste0(filename, "_", TDF_source, ".", type))
  
  ggsave(path, width = width, height = height, units = units)
  
  return(path)
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
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#'
#' @return A data frame
make_model_comparison_table <- function(best_models, 
                                        filename = "summary_MixSIAR_comparison", 
                                        TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  out <- purrr::map(best_models, 
                    ~ .x$loo_summary |> 
                      mutate(across(where(is.double), function(i) round(i, 2)))) |> 
    bind_rows(.id = "trophic_guild")
  
  path <- here::here("output", "tables", paste0(filename, "_", TDF_source, ".csv"))
  
  write_csv(out, path)
  
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
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#'
#' @return Paths to the saved txt files
save_MixSIAR_stats_diag <- function(best_models, 
                                    mix, 
                                    source,
                                    TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  folder <- "output/MixSIAR"
  out <- purrr::map(
    names(best_models),
    function(x) {
      model_type <- best_models[[x]]$name_best_model
      output_JAGS_custom(best_models[[x]]$best_model, 
                         mix[[x]][[model_type]], 
                         source,
                         output_options = set_output_options(
                           summary_save = TRUE,
                           summary_name = paste0(folder, "/summary_statistics/", x, "_", TDF_source),
                           diag_save = TRUE,
                           diag_name = paste0(folder, "/diagnostics/", x, "_", TDF_source),
                           sup_post_resid = TRUE))
      c(here::here(folder, "summary_statistics", paste0(x, "_", TDF_source, ".txt")),
        here::here(folder, "diagnostics", paste0(x, "_", TDF_source, ".txt")))
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
  
  purrr::map(
    names(best_models) |> 
      purrr::set_names(),
    function(x) {
      model_type <- best_models[[x]]$name_best_model
      capture.output(
        stats <- MixSIAR::output_stats(best_models[[x]]$best_model, 
                                       mix[[x]][[model_type]], 
                                       source, 
                                       output_options = set_output_options(
                                         summary_save = FALSE,
                                         return_obj = TRUE))
        )
      stats |>
        as.data.frame() |>
        tibble::rownames_to_column("parameter")
      }) |>
    bind_rows(.id = "trophic_guilds")
}

#' Save summary table with proportions for each trophic guild, family and species
#'
#' @param MixSIAR_stats Output from [make_MixSIAR_stats()]
#' @param consumers_clean List of paths to consumer csv files
#' @param filename Name of csv file which will be saved in `/output/tables`
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#'
#' @return Path to the saved csv file
make_MixSIAR_summary_table <- function(MixSIAR_stats,
                                       consumers_clean,
                                       filename = "summary_MixSIAR_proportions", 
                                       TDF_source = c("Post", "McCutchan")) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  consumers <- consumers_clean$data[[TDF_source]]
  
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
  
  path <- here::here("output", "tables", paste0(filename, "_", TDF_source, ".csv"))
  write_csv(out, path)
  
  return(path)
}

#' Make combined plot with isospace and relative source contributions
#' to each trophic guild
#'
#' @param consumers_clean List of paths to consumer csv files
#' @param models Output from [run_MixSIAR_models()]
#' @param best_models Output from [select_best_models()]
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#' @param colors Vector of colors with length equal to the number of trophic guilds
#' @param filename Name of the file to save
#' @param type File type (default to png)
#'
#' @return Path to the saved file
plot_isospace_and_rel_contributions <- function(consumers_clean, 
                                                models, 
                                                best_models, 
                                                TDF_source = c("Post", "McCutchan"),
                                                colors = RColorBrewer::brewer.pal(length(consumers_clean$path[[TDF_source]]), "Dark2"), 
                                                filename = "source_contribution_by_trophic_guild",
                                                type = "png") {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  mix <- consumers_clean$data[[TDF_source]]
  
  source <- data.frame(source = models$source$source_names) |>
    bind_cols(models$source$S_MU |> 
                as.data.frame() |> 
                dplyr::rename_with(~ paste("mean", .x, sep = "_"), starts_with("d"))) |>
    bind_cols(models$source$S_SIG |> 
                as.data.frame() |> 
                dplyr::rename_with(~ paste("sd", .x, sep = "_"), starts_with("d")))
  
  tdf <- models$discr$mu |> 
    as.data.frame() |> 
    tibble::rownames_to_column("source") |> 
    bind_cols(models$discr$sig2)
  
  source_tdf_corrected <- source |>
    left_join(tdf, by = "source") |>
    mutate(mean_d13C = mean_d13C + Meand13C,
           mean_d15N = mean_d15N + Meand15N,
           sd_d13C = sqrt(sd_d13C^2 + SDd13C), # TDF SDs are variances, already squared by MixSIAR
           sd_d15N = sqrt(sd_d15N^2 + SDd15N)) |>
    select(1:5)
  
  p <- purrr::map2(
    .x = 1:length(names(models$mix)),
    .y = colors,
    function(.x, .y) {
      guild <- names(models$mix)[.x]
      
      p1 <- source_tdf_corrected |>
        ggplot(aes(x = mean_d13C, 
                   y = mean_d15N)) +
        geom_point(data = mix |>
                     select(mean_d13C = d13C,
                            mean_d15N = d15N),
                   color = "grey",
                   alpha = 0.2) +
        geom_point(data = mix |>
                     filter(trophic_guild == guild) |>
                     rename(mean_d13C = d13C,
                            mean_d15N = d15N),
                   fill = .y,
                   alpha = 0.6,
                   shape = 21) +
        geom_linerange(aes(ymin = mean_d15N - sd_d15N,
                           ymax = mean_d15N + sd_d15N)) +
        geom_linerange(aes(xmin = mean_d13C - sd_d13C,
                           xmax = mean_d13C + sd_d13C)) +
        geom_point(aes(shape = source), size = 2.5, show.legend = FALSE) +
        facet_wrap(~trophic_guild) +
        theme_bw() +
        labs(x = "&delta;<sup>13</sup>C (&permil;)",
             y = "&delta;<sup>15</sup>N (&permil;)") +
        theme(axis.title.x = ggtext::element_markdown(size = 10),
              axis.title.y = ggtext::element_markdown(size = 10), 
              strip.text = element_text(face = "bold", size = 10),
              plot.margin = margin(1, 0, 1, 1, "pt"))
      
      R2jags::attach.jags(best_models[[guild]]$best_model)
      
      p2 <- p.global |> 
        as.data.frame() |> 
        rlang::set_names(source_tdf_corrected$source) |>
        tidyr::pivot_longer(cols = everything(), 
                            names_to = "source", 
                            values_to = "value") |>
        ggplot(aes(x = source, y = value, shape = source)) +
        ggdist::stat_gradientinterval(aes(fill_ramp = after_stat(y/max(y))),
                                      point_interval = "median_qi",
                                      .width = c(0.5, 0.95),
                                      point_size = 2.5,
                                      fill_type = "gradient", 
                                      scale = 0.7) + 
        ggdist::scale_fill_ramp_continuous(guide = ggdist::guide_rampbar(title = "Normalised density")) +
        labs(y = "Relative contribution",
             shape = "Source") +
        theme_classic() +
        theme(axis.title.y = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "bottom", 
              legend.title.position = "top",
              legend.title = element_text(hjust = 0.5),
              plot.margin = margin(1, 1, 1, 0, "pt"))
      
      if (.x %in% 1:6) p1 <- p1 + theme(axis.title.x = element_blank())
      
      if (.x %in% seq(1, length(unique(mix$trophic_guild)), by = 2)) {
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
  
  p <- patchwork::wrap_plots(p, ncol = 2) + 
    patchwork::plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom')
  
  path <- here::here("output", 
                     "figures", 
                     paste0(filename, "_", TDF_source, ".", type))
  
  ggsave(path, p, width = 7, height = 9)
  
  path
}

#' Plot relative source contributions for each species, family and trophic guild
#'
#' @param stats_df Output from [make_MixSIAR_stats()]
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#' @param filename Name of the file to save
#' @param type File type (default to png)
#' 
#' @return Path to the saved files
plot_all_rel_contribution <- function(stats_df,
                                      TDF_source = c("Post", "McCutchan"),
                                      filename = "rel_contributions_",
                                      type = "png") {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  p <- purrr::map(
    unique(stats_df$trophic_guilds) |>
      rlang::set_names(),
    ~ stats_df |>
      dplyr::filter(trophic_guilds == .x & startsWith(parameter, "p.")) |>
      dplyr::mutate(parameter = gsub("p\\.", "", parameter)) |>
      tidyr::separate_wider_delim(parameter, delim = ".", names = c("var", "source")) |>
      mutate(var_type = case_when(var == "global" ~ "",
                                  stringr::str_detect(var, " ") ~ "species",
                                  .default = "family"),
             var_type = factor(var_type, levels = c("", "family", "species")),
             var = ifelse(var_type == "species", glue::glue("<i>{var}</i>"), var)) |>
      ggplot(aes(x = `50%`, xmin = `2.5%`, xmax = `97.5%`, y = var, color = source)) +
      geom_pointrange(position = position_dodge(width = 0.7), fatten = 2) +
      scale_color_viridis_d() +
      facet_grid(rows = vars(var_type),
                 scales = "free_y",
                 space = "free") +
      labs(x = "Relative contribution",
           title = .x) + 
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
  )
  
  paths <- names(p) |>
    purrr::set_names() |>
    purrr::map(~ here::here("output", "figures", paste0(filename, .x, "_", TDF_source, ".", type)))
  
  purrr::map(1:length(p),
             ~ ggsave(paths[[.x]], p[[.x]], 
                      width = 12, 
                      height = ifelse(names(p)[.x] %in% c("invertivores-benthic", 
                                                          "macrocarnivores"), 
                                      24, 18), 
                      units = "cm"))
  
  paths
}

#' Plot total length effect on relative source contributions by trophic guild
#'
#' @param models Output from [run_MixSIAR_models()]
#' @param best_models Output from [select_best_models()]
#' @param TDF_source A string; "Post" or "McCutchan", the TDF used of trophic position estimation
#' @param filename Name of the file to save
#' @param type File type (default to png)
#' @param ... Arguments passed on to [plot_MixSIAR_continuous()]
#'
#' @return Path to the saved file
plot_rel_contribution_vs_tl <- function(models, 
                                        best_models, 
                                        TDF_source = c("Post", "McCutchan"), 
                                        filename = "rel_contributions_vs_tl_",
                                        type = "png", 
                                        ...) {
  
  TDF_source <- rlang::arg_match(TDF_source)
  
  p <- purrr::map(
    names(best_models) |> 
      rlang::set_names(),
    function(x) {
      name_best_model <- best_models[[x]]$name_best_model
      if (name_best_model %in% c("full", "tl")) {
        mod <- best_models[[x]]$best_model
        source <- models$source
        mix <- models$mix[[x]][[name_best_model]]
        plot_MixSIAR_continuous(mod, mix, source, ...) +
          scale_fill_viridis_d(aesthetics = c("color", "fill")) +
          theme(panel.grid = element_blank()) + 
          labs(title = x, 
               x = "Total length (cm) [natural-log]")
      }
    }
  )
  
  path <- here::here("output", "figures", paste0(filename, TDF_source, ".png"))
  
  p <- p |> 
    purrr::keep(purrr::negate(is.null)) |> 
    patchwork::wrap_plots(ncol = 2) + 
    patchwork::plot_layout(guides = "collect", 
                           axes = "collect") & 
    theme(legend.position = "bottom") &
    guides(fill_ramp = ggdist::guide_rampbar(title = "CI", 
                                             theme = theme(legend.title = element_text(vjust = 0.8))))
  
  ggsave(path, p, width = 7, height = 9)
  
  path
}
