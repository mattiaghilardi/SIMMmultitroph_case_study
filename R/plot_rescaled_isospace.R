#' Plot isospace with rescaled fish isotopic signatures
#'
#' @param sources Output of [prepare_source_data()]
#' @param consumers Output of [prepare_consumer_data()]
#'
#' @returns A ggplot object
plot_rescaled_isospace <- function(sources, 
                                   consumers) {
  
  
  # Summary of source isotopic signatures
  summary_sources <- readr::read_csv(sources) |>
    group_by(source) |>
    summarise(n = n(),
              mean_d13C = mean(d13C, na.rm = TRUE), 
              sd_d13C = sd(d13C, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE),
              .groups = "drop")
  
  # Plot
  p <- consumers |>
    bind_rows(.id = "TDF_source") |> 
    mutate(TDF_source = factor(TDF_source, levels = c("Post", "McCutchan"))) |> 
    ggplot() +
    # add mixing polygon
    geom_polygon(data = summary_sources |>
                   slice(chull(mean_d13C, mean_d15N)), aes(x = mean_d13C, y = mean_d15N),
                 alpha = 0.1) +
    # add fish
    geom_point(aes(x = d13C, y = d15N),
               alpha = 0.3) + 
    # add source means and SDs
    geom_linerange(data = summary_sources, 
                   aes(x = mean_d13C, y = mean_d15N,
                       ymin = mean_d15N - sd_d15N,
                       ymax = mean_d15N + sd_d15N,
                       colour = source)) +
    geom_linerange(data = summary_sources, 
                   aes(x = mean_d13C, y = mean_d15N,
                       xmin = mean_d13C - sd_d13C,
                       xmax = mean_d13C + sd_d13C,
                       colour = source)) +
    geom_point(data = summary_sources, 
               aes(x = mean_d13C, y = mean_d15N, shape = source, fill = source), 
               size = 2.5,
               colour = "black") +
    scale_shape_manual(values = 22:25) + 
    scale_fill_manual(values = c("#000000", "#999933", "#56B4E9", "#D55E00"),
                      aesthetics = c("colour", "fill")) +
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)",
         shape = "Source",
         colour = "Source",
         fill = "Source") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown()) + 
    facet_grid(~TDF_source)
  
  return(p)
}
