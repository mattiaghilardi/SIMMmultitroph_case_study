#' Plot correlation of mixing model estimates across methods
#'
#' @param stats_list Named list of outputs from [make_MixSIAR_stats()] for different analyses
#' @param combine_sources Logical; if sources should be combined
#' @param groups Only if `combine_sources = TRUE`. Named list; which sources to combine, 
#' and what names to give the new combined sources
#' @param colours Vector of colours to use for sources, plotted in alphabetic order.
#' If NULL (default), the default ggplot2 colour palette will be used
#' 
#' @return A patchwork object
plot_mixing_model_comparison <- function(stats_list = list(Post = MixSIAR_stats_Post,
                                                           McCutchan = MixSIAR_stats_McCutchan),
                                         combine_sources = FALSE,
                                         groups = NULL,
                                         colours = NULL) {
  
  # Get names of TDF sources
  TDF_sources <- names(stats_list)
  
  # Merge data
  df <- stats_list |> 
    dplyr::bind_rows(.id = "TDF_source") |> 
    # Keep proportions for each species
    dplyr::filter(startsWith(parameter, "p.")) |> 
    dplyr::select(trophic_guilds, parameter, `50%`, TDF_source) |> 
    dplyr::mutate(parameter = gsub("p\\.", "", parameter),
                  TDF_source = factor(TDF_source, levels = TDF_sources)) |>
    tidyr::separate_wider_delim(parameter, delim = ".", names = c("var", "source")) |>  
    dplyr::filter(stringr::str_detect(var, " "))
  
  # Get sources names and number
  source_names <- unique(df$source)
  n_sources <- length(source_names)
  
  # Combine sources
  if (combine_sources) {
    # Check combined sources
    is.list(groups) ||
      cli::cli_abort("If {.arg combine_sources} is TRUE, {.arg groups} must be a named list.")
    
    # Stack list and check names
    groups_df <- stack(groups)
    all(source_names %in% groups_df$values) ||
      cli::cli_abort(c("{.arg groups} does not include all initial sources.",
                       "Please correct source groups."))
    
    # Combine medians
    df <- df |>
      dplyr::left_join(groups_df,
                       by = c("source" = "values")) |> 
      dplyr::group_by(TDF_source, trophic_guilds, var, ind) |> 
      dplyr::summarise(`50%` = sum(`50%`), 
                       .groups = "drop") |> 
      dplyr::rename("source" = "ind")
    
    # New source names
    source_names <- names(groups)
    n_sources <- length(groups)
  }
  
  # Convert to wide format: one column per TDF source
  df2 <- df |> 
    tidyr::pivot_wider(names_from = TDF_source, values_from = `50%`)
  
  # Check colours
  if (!is.null(colours)) {
    length(colours) == n_sources ||
      cli::cli_abort(
        c("{.arg colours} length incorrect.", 
          "You provided {.value {length(colours)}} colours, but there are {.value {n_sources}} sources.")
      )
  }
  
  # Model formula
  formula <- paste0("mvbind(", paste(TDF_sources, collapse = ", "), ") ~ 1")
  
  # Model all data
  fit_all <- brms::brm(brms::bf(formula) + 
                         brms::set_rescor(TRUE),
                       data = df2,
                       cores = 4, 
                       backend = "cmdstan",
                       seed = NA) # Seed set through targets
  
  # Model by individual source
  fit_ind <- purrr::map(
    rlang::set_names(source_names),
    ~ brms::brm(brms::bf(formula) + 
                  brms::set_rescor(TRUE),
                data = df2 |> filter(source == .x),
                cores = 4, 
                backend = "cmdstan",
                seed = NA) # Seed set through targets
  )
  
  # Plot correlations
  p1 <- df2 |> 
    GGally::ggpairs(
      aes(colour = source), columns = 4:ncol(df2),
      upper = list(continuous = GGally::wrap(cor_custom, main_model = fit_all,
                                             sub_models = tibble(source = source_names,
                                                                 model = fit_ind),
                                             digits = 2, wrap_interval = FALSE)),
      lower = list(continuous = GGally::wrap(point_custom,
                                             shape = 19, alpha = 0.3, linecolor = "firebrick")),
      diag = list(continuous = GGally::wrap(diag_custom,
                                            xtext = "centre", ytext = Inf, vjust = 1.5, size = 4, fontface = "bold"))
    ) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(), 
          panel.border = element_rect(color = "black"))
  
  if (!is.null(colours)) {
    p1 <- p1 +
      scale_colour_manual(values = colours)
  }
  
  # Plot contrasts
  p2 <- purrr::map2(
    .x = fit_ind,
    .y = source_names,
    ~ brms::as_draws_df(.x) |> 
      select(ends_with("Intercept")) |> 
      tibble::rownames_to_column(".draw") |> 
      tidyr::pivot_longer(cols = -.draw) |> 
      mutate(name = gsub("b_",  "", name),
             name = gsub("_Intercept",  "", name)) |> 
      tidybayes::compare_levels(value, name) |> 
      mutate(source = .y)
  ) |> 
    bind_rows() |> 
    ggplot(aes(x = source, y = value, color = source)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggdist::stat_pointinterval(.width = c(0.5, 0.95)) +
    theme_bw() +
    ylab(paste(TDF_sources, collapse = "-")) +
    theme(axis.title.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  if (!is.null(colours)) {
    p2 <- p2 +
      scale_colour_manual(values = colours)
  }
  
  # Scatterplot by trophic guild
  p3 <- point_custom(df2, 
                     aes(x = Post, y = McCutchan, color = source), 
                     shape = 19, alpha = 0.3, linecolor = "firebrick") +
    facet_wrap(~trophic_guilds, ncol = 2) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))
  
  if (!is.null(colours)) {
    p3 <- p3 +
      scale_colour_manual(values = colours)
  }
  
  # Combine plots
  design <- "AAAACCC
             AAAACCC
             AAAACCC
             AAAACCC
             BBBBCCC
             BBBBCCC"
  
  p <- patchwork::wrap_plots(GGally::ggmatrix_gtable(p1),
                             p2,
                             p3) +
    patchwork::plot_layout(design = design) +
    patchwork::plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(0.025, 0.975),
          plot.tag = element_text(hjust = 1, vjust = 0, face = "bold")) 
  
  return(p)
}
