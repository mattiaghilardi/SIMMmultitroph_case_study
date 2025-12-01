#' Check variation in fish d13C and d15N across years
#'
#' @param sia_fish Raw sia_fish data
#'
#' @return A list with two brmsfit objects:
#'  - "fit_n": d15N model
#'  - "fit_c": d13C model
check_isotopes_across_years <- function(sia_fish) {
  
  # keep species in all 3 years
  sp_year <- sia_fish |> 
    select(year, species) |> 
    distinct() |>
    group_by(species) |> count()
  sia_fish_all_years <- sia_fish |> 
    filter(species %in% (sp_year |> filter(n > 2) |> pull(species))) |> 
    mutate(year = as.factor(year))
  
  # models
  fit_n <- brms::brm(d15N ~ 0 + year:species, 
                     data = sia_fish_all_years, 
                     cores = 4,
                     backend = "cmdstan",
                     seed = 123)
  fit_c <- brms::brm(d13C ~ 0 + year:species, 
                     data = sia_fish_all_years, 
                     cores = 4,
                     backend = "cmdstan",
                     seed = 123)
  
  brms::pp_check(fit_n)
  brms::pp_check(fit_c)
  brms::bayes_R2(fit_n) # 0.88
  brms::bayes_R2(fit_c) # 0.64
  
  # species-level predictions
  nd <- sia_fish_all_years |> 
    select(year, species) |> 
    distinct() |> 
    arrange(species)
  preds <- purrr::map2(
    .x = list(fit_n, fit_c),
    .y = c("&delta;<sup>15</sup>N (&permil;)", "&delta;<sup>13</sup>C (&permil;)"),
    ~ nd |> 
      tidybayes::add_epred_draws(.x) |> 
      mutate(isotope = .y)) |> 
    bind_rows()
  
  p1 <- sia_fish_all_years |> 
    tidyr::pivot_longer(cols = c(d15N, d13C), names_to = "isotope") |> 
    mutate(isotope = ifelse(isotope == "d15N",
                            "&delta;<sup>15</sup>N (&permil;)",
                            "&delta;<sup>13</sup>C (&permil;)")) |> 
    ggplot(aes(x = species, color = year)) + 
    geom_point(aes(y = value), alpha = 0.2, 
               position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) + 
    ggdist::stat_pointinterval(data = preds,
                               aes(y = .epred),
                               position = position_dodge(width = 0.8),
                               .width = c(0.5, 0.95),
                               point_size = 2) + 
    facet_grid(rows = vars(isotope), scales = "free", switch = "y") +
    scale_colour_viridis_d("Year", option = "G", end = 0.8) +
    theme_bw() + 
    theme(strip.text.y = ggtext::element_markdown(), 
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title = element_blank(), 
          axis.text.x = element_text(angle = 30, hjust = 1, face = "italic", size = 7),
          axis.text.y = element_text(size = 8),
          legend.box.spacing = unit(0.1, "cm"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 7))
  
  # contrasts
  p2 <- purrr::map2(
    .x = list(fit_n, fit_c),
    .y = c("&delta;<sup>15</sup>N", "&delta;<sup>13</sup>C"),
    ~ nd |> 
      tidybayes::add_epred_draws(.x) |> 
      ungroup() |> 
      group_by(species) |> 
      tidybayes::compare_levels(.epred, by = year) |> 
      mutate(isotope = .y) |> 
      ungroup()) |> 
    bind_rows() |> 
    ggplot(aes(x = .epred, y = year)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_gradientinterval(geom = "slab") +
    ggdist::stat_pointinterval(aes(color = species), 
                               .width = c(0.5, 0.95),
                               position = position_dodge(width = 0.85), 
                               point_size = 1.5,
                               alpha = 0.5) +
    ggdist::stat_pointinterval(.width = c(0.5, 0.95)) + 
    ggdist::scale_slab_alpha_continuous(guide = "none") + 
    facet_grid(cols = vars(isotope), scales = "free") +
    scale_color_viridis_d(option = "H") +
    labs(x = "Difference in isotopic signature (&permil;)",
         y = "Contrast",
         color = "Species") +
    theme_bw() +
    theme(strip.text.x = ggtext::element_markdown(),
          legend.text = element_text(size = 6, face = "italic"),
          legend.title = element_text(size = 9),
          axis.title.x = ggtext::element_markdown(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(size = 8))
  
  p1 / p2 + 
    patchwork::plot_layout(
      heights = c(0.7, 1),
      design = "
  AAAAAAAAAAAAAAAAAAAAAAAAAA
  #BBBBBBBBBBBBBBBBBBBBBBBBB
  ", 
      guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") & 
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(hjust = 1, vjust = 0, face = "bold"))
  
  ggsave(here::here("output", "figures", "isotope_variation_across_years.png"),
         width = 20, height = 18, units = "cm")
  
  list(fit_n = fit_n, fit_c = fit_c)
}
