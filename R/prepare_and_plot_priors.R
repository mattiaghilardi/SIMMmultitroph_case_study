#' Make list of custom Dirichlet priors on p.global
#' 
#' % algae divided equally between algal groups
#'
#' @param priors Data frame with percentage for algae, cyanobacteria and POM 
#' for each trophic guild
#'
#' @return A list of tibbles, one per trophic guild
prepare_prior_list <- function(priors) {
  
  # Number of sources
  n_sources <- 4
  
  # Make list of Dirichlet priors on p.global
  prior_list <- priors |> 
    # Divide the percentage of algae equally between the groups
    dplyr::mutate(green_brown_algae = algae * 1/2,
                  red_algae = algae * 1/2) |> 
    # Remove algae
    dplyr::select(-algae) |> 
    # Convert percentage in proportion and 
    # multiply by number of source as suggested by Stock et al 2018
    dplyr::mutate(dplyr::across(dplyr::where(is.double), ~ .x / 100 * n_sources)) |> 
    # Convert to long format
    tidyr::pivot_longer(cols = -trophic_guild, 
                        names_to = "source", 
                        values_to = "prior") |> 
    # Arrange sources alphabetically as sources in MixSIAR
    dplyr::arrange(source) |> 
    # Split by trophic guild
    split(~trophic_guild) |> 
    # Extract prior vector
    purrr::map(~ dplyr::pull(.x, prior))
  
  return(prior_list)
}

#' Plot Dirichlet priors
#'
#' @param sources Output of [prepare_source_data()]
#' @param prior_list Output of [prepare_prior_list()]
#'
#' @return A ggplot object
plot_priors <- function(sources,
                        prior_list) {
  
  # source names
  source_names <- read.csv(sources)$source |> 
    as.factor() |> 
    levels()
  
  # trophic guild labels
  labels <- prior_list |> 
    purrr::map(~ paste0("&alpha; = (", 
                        paste(.x, collapse = ", "), 
                        ")") |> 
                 as.data.frame() |> 
                 rlang::set_names("prior")) |> 
    bind_rows(.id = "trophic_guild") |> 
    mutate(label = paste(trophic_guild, prior, sep = "<br><br>"))
  
  # plot
  p <- prior_list |> 
    purrr::map(~ brms::rdirichlet(10000, .x) |> 
                 as.data.frame() |> 
                 rlang::set_names(source_names) |> 
                 tidyr::pivot_longer(cols = everything(), 
                                     names_to = "source")) |> 
    bind_rows(.id = "trophic_guild") |> 
    left_join(labels) |> 
    ggplot(aes(x = value, fill = source)) + 
    geom_density() + 
    facet_grid(rows = vars(label), 
               cols = vars(source),
               scales = "free_y") + 
    labs(x = "Proportion", 
         y = "Density") + 
    scale_fill_manual(values = c("#000000", "#999933", "#56B4E9", "#D55E00")) + 
    scale_x_continuous(limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2), 
                       labels = round(seq(0, 1, by = 0.2), 1)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = "none",
          strip.text.y = ggtext::element_markdown(angle = 0))
  
  return(p)
}
