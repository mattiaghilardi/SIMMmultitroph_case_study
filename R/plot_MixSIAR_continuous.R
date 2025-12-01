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
#' @param ... Other arguments passed to [ggdist::stat_lineribbon()] if `plot_type = "lineribbon`,
#' to [ggplot2::geom_line()] if `plot_type = "spaghetti`, 
#' or to [ggdist::stat_ribbon()] if `plot_type = "lineribbon_gradient`
#'
#' @return A ggplot object
#' 
plot_MixSIAR_continuous <- function(jags.1, mix, source, 
                                    fac1 = NULL, fac2 = NULL, 
                                    exclude_sources_below = 0.1, 
                                    plot_type = c("lineribbon", "spaghetti", "lineribbon_gradient"), 
                                    ndraws = 100,
                                    .width = if (plot_type == "lineribbon") c(0.5, 0.8, 0.95) 
                                      else if (plot_type == "lineribbon_gradient") ppoints(30),
                                    add_line_CE_center = FALSE,
                                    add_text_CE_center = FALSE,
                                    ...) {
  
  plot_type <- rlang::arg_match(plot_type)
  
  R2jags::attach.jags(jags.1)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
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
  
  label <- mix$cont_effects[1]
  cont <- mix$CE[[1]]
  ilr.cont <- get(paste("ilr.cont", 1, sep = ""))
  
  n.plot = 200
  chain.len = dim(p.global)[1]
  Cont1.plot <- seq(from = round(min(cont), 1), to = round(max(cont), 1), length.out = n.plot)
  ilr.plot <- array(NA, dim = c(n.plot, n.sources - 1, chain.len))
  for(src in 1:n.sources - 1) {
    for(i in 1:n.plot){
      if (!is.null(fac1) & !is.null(fac2))
        ilr.plot[i, src,] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i] + ilr.fac1[, f1, src] + ilr.fac2[, f2, src]
      else if (!is.null(fac1) & is.null(fac2))
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
  
  # dummy variables for inverse ILR calculation
  cross <- array(data = NA, dim = c(n.plot, chain.len, n.sources, n.sources - 1))  
  tmp <- array(data = NA, dim = c(n.plot, chain.len, n.sources))  
  p.plot <- array(data = NA, dim = c(n.plot, chain.len, n.sources))  
  for(i in 1:n.plot){
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
  
  # transform Cont1.plot (x-axis) back to the original scale
  Cont1.plot <- Cont1.plot * mix$CE_scale + mix$CE_center
  
  # make data frame
  df <- p.plot |>
    apply(3, 
          function(x) {
            x |> 
              as.data.frame() |> 
              dplyr::mutate(x = Cont1.plot) |> 
              tidyr::pivot_longer(cols = 1:all_of(chain.len), 
                                  names_to = NULL, 
                                  values_to = ".value") |>
              dplyr::mutate(.draw = rep(1:chain.len, n.plot))
          }
    ) |>
    purrr::set_names(source_names) |>
    dplyr::bind_rows(.id = "source")
  
  # remove sources from plot with very low proportions
  df_median <- df |>
    dplyr::group_by(source, x) |>
    dplyr::summarise(.value = median(.value))
  rm.srcs <- df_median |>
    tidyr::pivot_wider(names_from = source, values_from = .value) |>
    dplyr::select(-x) |>
    apply(2, function(x) all(x < exclude_sources_below))
  df <- dplyr::filter(df, source %in% source_names[!rm.srcs])
  df_median <- dplyr::filter(df_median, source %in% source_names[!rm.srcs])
  
  # plot
  if (plot_type == "spaghetti") {
    p <- df |>
      dplyr::filter(.draw %in% sample(1:chain.len, ndraws)) |>
      dplyr::group_by(source) |>
      ggplot(aes(x = x, y = .value, color = source)) +
      geom_line(aes(y = .value, group = paste(source, .draw)), ...) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") +
      geom_line(data = df_median,
                linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2.5) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply")
  }
    
  if (plot_type == "lineribbon") {
    p <- df |>
      ggplot(aes(x = x, y = .value, color = source, fill = source)) +
      ggdist::stat_ribbon(aes(fill_ramp = after_stat(level)), 
                          .width = .width,
                          ...) |> 
      ggblend::blend("multiply") +
      geom_line(data = df_median,
                linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2.5) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      ggdist::scale_fill_ramp_discrete(name = "CI")
  }
  
  if (plot_type == "lineribbon_gradient") {
    p <- df |>
      ggplot(aes(x = x, y = .value, color = source, fill = source)) +
      ggdist::stat_ribbon(aes(fill_ramp = after_stat(.width)), 
                          .width = .width,
                          ...) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      geom_line(data = df_median,
                linewidth = 1) |> 
      ggblend::copy_under(color = "white", linewidth = 2.5) |> 
      ggblend::partition(vars(source)) |> 
      ggblend::blend("multiply") + 
      ggdist::scale_fill_ramp_continuous(range = c(1, 0), 
                                         guide = ggdist::guide_rampbar(title = "CI"))
  }
  
  if (add_line_CE_center) {
    p <- p + 
      geom_vline(xintercept = mix$CE_center, linetype = "dashed")
    if (add_text_CE_center)
      p <- p + 
        annotate("label", x = mix$CE_center, y = Inf, 
                 label = paste("mean =", round(mix$CE_center, 2)), 
                 vjust = 1)
  }
  
  p +
    labs(title = fac.lab,
         y = "Relative contribution",
         x = label) +
    theme_bw()
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
