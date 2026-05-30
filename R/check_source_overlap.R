#' Check overlap between source isotopic signatures
#'
#' @param sia_sources Raw sia_sources data
#' 
#' @return A list with two elements:
#'  - "models": a list of two brmsfit objects
#'    - "fit_n": d15N model
#'    - "fit_c": d13C model
#'  - "plot": a patchwork object
check_source_overlap <- function(sia_sources) {
  
  # Summary of isotopic signatures
  summary_sources <- sia_sources |>
    group_by(source) |>
    summarise(n = n(),
              mean_d13C = mean(d13C, na.rm = TRUE), 
              sd_d13C = sd(d13C, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE),
              .groups = "drop")
  
  # Plot isospace of sources
  p1 <- summary_sources |> 
    ggplot(aes(x = mean_d13C, 
               y = mean_d15N,
               color = source)) +
    geom_linerange(aes(xmin = mean_d13C - sd_d13C,
                       xmax = mean_d13C + sd_d13C)) +
    geom_pointrange(aes(ymin = mean_d15N - sd_d15N,
                        ymax = mean_d15N + sd_d15N)) +
    geom_point(data = sia_sources,
               aes(x = d13C, y = d15N),
               alpha = 0.4) + 
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)",
         color = "Source") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  # Bayesian regression models
  fit_n <- brms::brm(d15N ~ 0 + source, 
                     data = sia_sources,
                     cores = 4,
                     backend = "cmdstan",
                     seed = NA) # seed set through targets
  fit_c <- brms::brm(d13C ~ 0 + source, 
                     data = sia_sources,
                     cores = 4,
                     backend = "cmdstan",
                     seed = NA) # seed set through targets
  
  brms::pp_check(fit_n)
  brms::pp_check(fit_c)
  brms::bayes_R2(fit_n) # 0.81
  brms::bayes_R2(fit_c) # 0.94
  
  # Plot contrasts
  p2 <- purrr::map2(
    .x = list(fit_n, fit_c),
    .y = c("&delta;<sup>15</sup>N", "&delta;<sup>13</sup>C"),
    ~ sia_sources |> 
      select(source) |>
      distinct() |> 
      tidybayes::add_epred_draws(.x) |> 
      tidybayes::compare_levels(.epred, by = source) |> 
      mutate(isotope = .y, 
             above = ifelse(.epred > 0, 1, 0), 
             prob = sum(above) / n(), 
             prob = ifelse(prob >= 0.5, prob, 1 - prob))) |> 
    bind_rows() |> 
    ggplot(aes(x = .epred, y = source, fill = prob)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_halfeye(normalize = "panels",
                         # slab_alpha = 0.5,
                         .width = c(0.5, 0.95),
                         scale = 0.8,
                         point_size = 2) +
    facet_grid(cols = vars(isotope), scales = "free") +
    scale_fill_gradient2("Probability of difference",
                         midpoint = 0.75, 
                         low = "grey99", mid = "grey", high = "deepskyblue") +
    theme_bw() + 
    xlab("Difference in isotopic signature (&permil;)") +
    scale_y_discrete(position = "right") +
    theme(strip.text.x = ggtext::element_markdown(face = "bold"), 
          axis.title.y = element_blank(),
          axis.title.x = ggtext::element_markdown(),
          legend.position = "bottom")
  
  # Final plot
  p <- p1 + p2 + 
    patchwork::plot_layout(
      design = "
      AAAAAAAAAA#
      AAAAAAAAAA#
      BBBBBBBBBBB
      BBBBBBBBBBB
      BBBBBBBBBBB
      ") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(hjust = 1, vjust = 0, face = "bold"))
  
  return(list(models = list(fit_n = fit_n, 
                            fit_c = fit_c), 
              plot = p)
  )
}
