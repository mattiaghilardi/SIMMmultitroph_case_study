#' Plot initial isospace
#'
#' @param sia_sources Raw sia_sources data
#' @param sia_baselines_corrected Output of [baselines_lipid_correction()]
#' @param sia_fish_corrected Output of [fish_lipid_correction()]
#'
#' @returns A ggplot object
plot_initial_isospace <- function(sia_sources, 
                                  sia_baselines_corrected, 
                                  sia_fish_corrected) {
  
  # Summary of source isotopic signatures
  summary_sources <- sia_sources |>
    group_by(source) |>
    summarise(n = n(),
              mean_d13C = mean(d13C, na.rm = TRUE), 
              sd_d13C = sd(d13C, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE),
              .groups = "drop")
  
  # Summary of baseline isotopic signatures
  summary_baselines <- sia_baselines_corrected |>
    group_by(baseline) |>
    summarise(n = n(),
              mean_d13C = mean(d13C_corrected, na.rm = TRUE), 
              sd_d13C = sd(d13C_corrected, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE),
              .groups = "drop")
  
  # Plot
  p <- summary_sources |> 
    ggplot(aes(x = mean_d13C, y = mean_d15N)) + 
    # add mixing polygon
    geom_polygon(data = summary_sources |> 
                   slice(chull(mean_d13C, mean_d15N)), 
                 alpha = 0.1) +
    # add source means and SDs
    geom_linerange(aes(ymin = mean_d15N - sd_d15N, 
                       ymax = mean_d15N + sd_d15N,
                       color = source)) +
    geom_pointrange(aes(xmin = mean_d13C - sd_d13C, 
                        xmax = mean_d13C + sd_d13C,
                        color = source)) + 
    # add baseline means and SDs
    geom_linerange(data = summary_baselines,
                   aes(ymin = mean_d15N - sd_d15N, 
                       ymax = mean_d15N + sd_d15N)) +
    geom_pointrange(data = summary_baselines,
                    aes(xmin = mean_d13C - sd_d13C, 
                        xmax = mean_d13C + sd_d13C,
                        shape = baseline)) +
    # add source individual samples
    geom_point(data = sia_sources,
               aes(x = d13C, y = d15N, color = source),
               alpha = 0.4) +
    # add baseline individual samples
    geom_point(data = sia_baselines_corrected,
               aes(x = d13C_corrected, y = d15N, shape = baseline)) +
    # add fish
    geom_point(data = sia_fish_corrected, 
               aes(x = d13C_corrected, y = d15N),
               color = "grey",
               alpha = 0.4) +
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)",
         color = "Source",
         shape = "Baseline") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  return(p)
}

# ggsave("output/figures/initial_isospace.png",
#        width = 18, height = 16, units = "cm")
