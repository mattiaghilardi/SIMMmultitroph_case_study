#' Plot correlation of mixing model estimates across methods
#'
#' @param stats_list Named list of outputs from [make_MixSIAR_stats()] for different analyses
#' @param type File type (default to png)
#'
#' @return Path to the saved file
plot_mixing_model_comparison <- function(stats_list = list(Post = MixSIAR_stats_Post,
                                                           McCutchan = MixSIAR_stats_McCutchan),
                                         type = "png") {
  
  # get names
  TDF_source <- names(stats_list)
  
  # merge data
  df <- dplyr::bind_rows(stats_list, .id = "TDF_source")
  
  # keep proportions for each species
  df2 <- df |>
    dplyr::filter(startsWith(parameter, "p.")) |> 
    dplyr::select(trophic_guilds, parameter, `50%`, TDF_source) |> 
    tidyr::pivot_wider(names_from = TDF_source, values_from = `50%`) |> 
    dplyr::mutate(parameter = gsub("p\\.", "", parameter)) |>
    tidyr::separate_wider_delim(parameter, delim = ".", names = c("var", "source")) |>  
    dplyr::filter(stringr::str_detect(var, " ")) # keep only species
  
  formula <- paste0("mvbind(", paste(TDF_source, collapse = ", "), ") ~ 1")
  
  # model all data
  fit_all <- brms::brm(brms::bf(formula) + 
                         brms::set_rescor(TRUE),
                       data = df2,
                       cores = 4, 
                       backend = "cmdstan",
                       seed = 123)
  
  # model cyanobacteria
  fit_cyano <- brms::brm(brms::bf(formula) + 
                           brms::set_rescor(TRUE),
                         data = df2 |> filter(source == "Cyanobacteria"),
                         cores = 4, 
                         backend = "cmdstan",
                         seed = 123)
  
  # model macroalgae
  fit_algae <- brms::brm(brms::bf(formula) + 
                           brms::set_rescor(TRUE),
                         data = df2 |> filter(source == "Macroalgae"),
                         cores = 4, 
                         backend = "cmdstan",
                         seed = 123)
  
  # model pom
  fit_pom <- brms::brm(brms::bf(formula) + 
                         brms::set_rescor(TRUE),
                       data = df2 |> filter(source == "POM"),
                       cores = 4, 
                       backend = "cmdstan",
                       seed = 123)
  
  # plot correlations
  p1 <- df2 |> 
    GGally::ggpairs(
      aes(color = source), columns = 4:ncol(df2),
      upper = list(continuous = GGally::wrap(cor_custom, main_model = fit_all,
                                             sub_models = tibble(source = c("Cyanobacteria", "Macroalgae", "POM"), 
                                                                 model = list(fit_cyano, fit_algae, fit_pom)),
                                             digits = 2, wrap_interval = FALSE)),
      lower = list(continuous = GGally::wrap(point_custom,
                                             shape = 19, alpha = 0.3, linecolor = "firebrick")),
      diag = list(continuous = GGally::wrap(diag_custom,
                                            xtext = "centre", ytext = Inf, vjust = 1.5, size = 4, fontface = "bold"))
    ) + 
    scale_color_viridis_d() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(), 
          panel.border = element_rect(color = "black"))
  
  # plot contrasts
  p2 <- purrr::map2(
    .x = list(fit_cyano, fit_algae, fit_pom),
    .y = c("Cyanobacteria", "Macroalgae", "POM"),
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
    ggdist::stat_pointinterval(.width = c(0.5, 0.95),
                               position = position_dodge(width = 0.5)) +
    scale_color_viridis_d(aesthetics = c("color", "fill")) +
    theme_bw() +
    ylab(#paste("Difference in estimated\nrelative contribution",
               paste(TDF_source, collapse = "-")#,
               #sep = "\n")
  ) +
    theme(axis.title.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # scatterplot by trophic guild
  p3 <- point_custom(df2, 
                     aes(x = Post, y = McCutchan, color = source), 
                     shape = 19, alpha = 0.3, linecolor = "firebrick") +
    facet_wrap(~trophic_guilds, ncol = 2) +
    scale_color_viridis_d() +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(strip.background = element_blank(), 
          legend.position = "none")
  
  # combine plots
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
  
  # save plot
  path <- here::here("output", "figures", paste0("mixing_model_comparison", ".", type))
  ggsave(path, p, width = 18, height = 18, units = "cm")
  
  path
}
