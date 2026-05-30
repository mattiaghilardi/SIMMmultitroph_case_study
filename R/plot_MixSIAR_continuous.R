#' Plot proportions by a continuous covariate in a MixSIAR model
#'
#' This function uses code from [MixSIAR::plot_continuous_var()]. 
#' It creates a plot of how the mixture proportions change according to a continuous covariate.
#' It allows to plot the global proportions as well as proportions of specific factor levels,
#' including nested factors.
#' 
#' @param jags.1 Output from [MixSIAR::run_model()]
#' @param mix Output from [MixSIAR::load_mix_data()]
#' @param source Output from [MixSIAR::load_source_data()]
#' @param fac1 A string, one of the levels of factor 1. 
#' If NULL (Default) plots the global proportions
#' @param fac2 A string, one of the levels of factor 2. 
#' If specified, then also `fac1` must be specified
#' @param exclude_sources_below Don't plot sources with median proportion below this 
#' level for the entire range of continuous effect variable (default = 0.1)
#' @param plot_type One of "lineribbon" (Default), "spaghetti", or "lineribbon_gradient"
#' @param ndraws Only if `plot_type = "spaghetti"`. The number of posterior draws to plot
#' @param .width Only if `plot_type = "lineribbon"` or `plot_type = "lineribbon_gradient`.
#' The `.width` argument passed to `ggdist::point_interval()`
#' @param add_line_CE_center Logical; if a vertical dashed line at the mean of the 
#' continuous effect variable should be added. Values of sources crossing this line  
#' are those reported by MixSIAR. Default to FALSE
#' @param add_text_CE_center Logical; if a label with the mean of the continuous effect 
#' variable should be added on top of the vertical line. Default to FALSE
#' @param resolution Number of support points used to generate the plot. 
#' Higher resolution leads to smoother plots. Defaults to 100. 
#' It might be necessary to reduce resolution when only few RAM is available.
#' @param combine_sources Logical; if sources should be combined
#' @param groups Only if `combine_sources = TRUE`. Named list; which sources to combine, 
#' and what names to give the new combined sources
#' @param ... Other arguments passed to [ggdist::stat_ribbon()] if 
#' `plot_type = "lineribbon"` or `plot_type = "lineribbon_gradient"`, or
#' to [ggplot2::geom_line()] if `plot_type = "spaghetti"`
#'
#' @return A ggplot object
#' 
plot_MixSIAR_continuous <- function(jags.1, 
                                    mix, 
                                    source, 
                                    fac1 = NULL, 
                                    fac2 = NULL, 
                                    exclude_sources_below = 0.1, 
                                    plot_type = c("lineribbon", "spaghetti", "lineribbon_gradient"), 
                                    ndraws = 100,
                                    .width = if (plot_type == "lineribbon") c(0.5, 0.8, 0.95) 
                                      else if (plot_type == "lineribbon_gradient") ppoints(30),
                                    add_line_CE_center = FALSE,
                                    add_text_CE_center = FALSE,
                                    resolution = 100,
                                    combine_sources = FALSE,
                                    groups = NULL,
                                    ...) {
  
  # Check plot type
  plot_type <- rlang::arg_match(plot_type)
  
  # Get factor labels
  f1 <- f2 <- fac.lab <- NULL
  if (!is.null(fac1)) {
    fac1 %in% mix$FAC[[1]]$labels ||
      cli::cli_abort("{.val {fac1}} is not a valid level of factor {.val {mix$FAC[[1]]$name}}")
    f1 <- which(mix$FAC[[1]]$labels == fac1)
    fac.lab <- mix$FAC[[1]]$labels[f1]
  }
  if (!is.null(fac2)) {
    !is.null(fac1) ||
      cli::cli_abort("If {.arg fac2} is not NULL, {.arg fac1} must also be specified.")
    fac2 %in% mix$FAC[[2]]$labels ||
      cli::cli_abort("{.val {fac2}} is not a valid level of factor {.val {mix$FAC[[2]]$name}}")
    f2 <- which(mix$FAC[[2]]$labels == fac2)
    fac.lab <- paste(fac.lab, mix$FAC[[2]]$labels[f2], sep = " - ")
  }
  
  # Prepare data
  df <- prepare_cont_data(jags.1 = jags.1, 
                          mix = mix, 
                          source = source, 
                          f1 = f1, 
                          f2 = f2, 
                          exclude_sources_below = exclude_sources_below, 
                          resolution = resolution,
                          combine_sources = combine_sources,
                          groups = groups)
  
  # Median values
  df_median <- df |>
    dplyr::group_by(source, x) |>
    dplyr::summarise(.value = median(.value),
                     .groups = "drop")
  
  # Remove heavy object
  rm(jags.1)
  
  # Plot
  if (plot_type == "spaghetti") {
    p <- df |>
      dplyr::filter(.draw %in% sample(1:chain.len, ndraws)) |>
      dplyr::group_by(source) |>
      ggplot2::ggplot(ggplot2::aes(x = x, y = .value, color = source)) +
      ggplot2::geom_line(ggplot2::aes(y = .value, group = paste(source, .draw)), 
                         ...) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") +
      ggplot2::geom_line(data = df_median,
                         linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply")
  }
    
  if (plot_type == "lineribbon") {
    p <- df |>
      ggplot2::ggplot(ggplot2::aes(x = x, y = .value, color = source, fill = source)) +
      ggdist::stat_ribbon(aes(fill_ramp = ggplot2::after_stat(level)), 
                          .width = .width,
                          ...) |> 
      ggblend::blend("multiply") +
      ggplot2::geom_line(data = df_median,
                         linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      ggdist::scale_fill_ramp_discrete(name = "CI")
  }
  
  if (plot_type == "lineribbon_gradient") {
    p <- df |>
      ggplot2::ggplot(ggplot2::aes(x = x, y = .value, color = source, fill = source)) +
      ggdist::stat_ribbon(aes(fill_ramp = ggplot2::after_stat(.width)), 
                          .width = .width,
                          ...) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      ggplot2::geom_line(data = df_median,
                         linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      ggdist::scale_fill_ramp_continuous(range = c(1, 0), 
                                         guide = ggdist::guide_rampbar(title = "CI"))
  }
  
  if (add_line_CE_center) {
    p <- p + 
      ggplot2::geom_vline(xintercept = mix$CE_center, linetype = "dashed")
    if (add_text_CE_center)
      p <- p + 
        ggplot2::annotate("label", x = mix$CE_center, y = Inf, 
                          label = paste("mean =", round(mix$CE_center, 2)), 
                          vjust = 1)
  }
  
  p <- p +
    ggplot2::labs(title = fac.lab,
                  y = "Relative contribution",
                  x = mix$cont_effects[1]) +
    ggplot2::theme_bw()
  
  return(p)
}

#' Prepare data for plotting continuous effect
#' 
#' @param f1 The index of the level of factor 1. 
#' If NULL (Default) returns the global proportions
#' @param f2 The index of the level of factor 2. 
#' If specified, then also `f1` must be specified
prepare_cont_data <- function(jags.1, 
                              mix, 
                              source, 
                              f1 = NULL, 
                              f2 = NULL, 
                              exclude_sources_below = 0.1, 
                              resolution = 100,
                              combine_sources = FALSE,
                              groups = NULL) {
  
  # Attach model
  R2jags::attach.jags(jags.1)
  
  # Get number and names of sources
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  # Check combined sources
  if (combine_sources) {
    is.list(groups) ||
      cli::cli_abort("If {.arg combine_sources} is TRUE, {.arg groups} must be a named list.")
    
    # Stack list
    groups_df <- stack(groups)
    all(source_names %in% groups_df$values) ||
      cli::cli_abort(c("{.arg groups} does not include all initial sources.",
                       "Please correct source groups."))
  }
  
  # Extract continuous effect
  label <- mix$cont_effects[1]
  cont <- mix$CE[[1]]
  ilr.cont <- get(paste("ilr.cont", 1, sep = ""))
  
  # Compute proportions in ILR-space
  chain.len = dim(p.global)[1]
  Cont1.plot <- seq(from = round(min(cont), 1), to = round(max(cont), 1), length.out = resolution)
  ilr.plot <- array(NA, dim = c(resolution, n.sources - 1, chain.len))
  for(src in 1:n.sources - 1) {
    for(i in 1:resolution){
      if (!is.null(f1) & !is.null(f2))
        ilr.plot[i, src,] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i] + ilr.fac1[, f1, src] + ilr.fac2[, f2, src]
      else if (!is.null(f1) & is.null(f2))
        ilr.plot[i, src,] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i] + ilr.fac1[, f1, src]
      else
        ilr.plot[i, src,] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i]
    }
  }
  
  # Transform every draw from ILR-space to p-space
  e <- matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, ncol = (n.sources - 1))
  for(i in 1:(n.sources-1)){
    e[, i] <- exp(c(rep(sqrt(1 / (i * (i + 1))), i), -sqrt(i / (i + 1)), rep(0, n.sources - i - 1)))
    e[, i] <- e[, i] / sum(e[, i])
  }
  
  # Dummy variables for inverse ILR calculation
  cross <- array(data = NA, dim = c(resolution, chain.len, n.sources, n.sources - 1))  
  tmp <- array(data = NA, dim = c(resolution, chain.len, n.sources))  
  p.plot <- array(data = NA, dim = c(resolution, chain.len, n.sources))  
  for(i in 1:resolution){
    for(d in 1:chain.len){
      for(j in 1:(n.sources - 1)){
        cross[i, d, , j] <- (e[, j]^ilr.plot[i, j, d]) / sum(e[, j]^ilr.plot[i, j, d]);
      }
      for(src in 1:n.sources){
        tmp[i, d, src] <- prod(cross[i, d, src, ]);
      }
      for(src in 1:n.sources){
        p.plot[i, d, src] <- tmp[i, d, src] / sum(tmp[i, d, ]);
      }
    }
  }
  
  # Transform Cont1.plot (x-axis) back to the original scale
  Cont1.plot <- Cont1.plot * mix$CE_scale + mix$CE_center
  
  # Make data frame
  df <- p.plot |>
    apply(3, 
          function(x) {
            x |> 
              as.data.frame() |> 
              dplyr::mutate(x = Cont1.plot) |> 
              tidyr::pivot_longer(cols = 1:all_of(chain.len), 
                                  names_to = NULL, 
                                  values_to = ".value") |>
              dplyr::mutate(.draw = rep(1:chain.len, resolution))
          }
    ) |>
    purrr::set_names(source_names) |>
    dplyr::bind_rows(.id = "source")
  
  # Combine sources
  if (combine_sources) {
    df <- df |>
      dplyr::left_join(groups_df,
                       by = c("source" = "values")) |> 
      dplyr::group_by(ind, x, .draw) |> 
      dplyr::summarise(.value = sum(.value), 
                       .groups = "drop") |> 
      dplyr::rename("source" = "ind")
    
    # New source names
    source_names <- groups_df$ind
  }
  
  # Remove sources with very low proportions
  df_median <- df |>
    dplyr::group_by(source, x) |>
    dplyr::summarise(.value = median(.value),
                     .groups = "drop")
  rm.srcs <- df_median |>
    tidyr::pivot_wider(names_from = source, values_from = .value) |>
    dplyr::select(-x) |>
    apply(2, function(x) all(x < exclude_sources_below))
  df <- dplyr::filter(df, source %in% source_names[!rm.srcs])
  
  return(df)
}

# # Test
# jags1 <- MixSIAR_models_Post$models$`invertivores-benthic`$full
# source <- MixSIAR_models_Post$source
# mix <- MixSIAR_models_Post$mix$`invertivores-benthic`$full
# 
# #spaghetti
# plot_MixSIAR_continuous(jags1, mix, source, fac1 = "Holocentridae", fac2 = "Myripristis berndti", plot_type = "spaghetti", alpha = 0.25) +
#   scale_color_viridis_d()
# # lineribbon
# plot_MixSIAR_continuous(jags1, mix, source, fac1 = "Holocentridae", .width = c(0.5, 0.95), alpha = 0.5) +
#   scale_fill_viridis_d(aesthetics = c("color", "fill"))
# # lineribbon gradients
# plot_MixSIAR_continuous(jags1, mix, source, fac1 = "Holocentridae", plot_type = "lineribbon_gradient") +
#   scale_fill_viridis_d(aesthetics = c("color", "fill")) +
#   theme(panel.grid = element_blank())
