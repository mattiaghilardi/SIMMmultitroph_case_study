#' Fit a Bayesian model to estimate TP
#'
#' @param data a list containing the data
#' @param model.string model string containing the model description
#' @param n.adapt number of iterations for adaptation (initial sampling phase)
#' @param n.iter integer defining the number of iterations
#' @param burnin number of iterations discarded as burn in
#' @param thin thinning interval to get posterior samples
#' @param n.chains number of parallel chains for the model
#' @param quiet logical value to indicate whether messages 
#' generated during compilation will be suppressed, 
#' as well as the progress bar during adaptation
#' @param .point Point summary function, which takes a vector and 
#' returns a single value, e.g. mean(), median(), or `ggdist::Mode()`
#' @param .width vector of probabilities to use that determine the 
#' widths of the resulting intervals
#' @param .interval Interval function, which takes a vector and a 
#' probability (.width) and returns a two-element vector representing 
#' the lower and upper bound of an interval; e.g. `ggist::qi()`, `ggist::hdi()`
#' @param ... additional arguments passed to `rjags::jags.model`
#'
#' @return A list with 4 elements:
#'  - "summary": results summary
#'  - "TP": all posterior draws for TP
#'  - "muDeltaN": all posterior draws for muDeltaN
#'  - "alpha": all posterior draws for alpha
#'  - "samples": output of `tRophicPosition::posteriorTP`
fit_TP_model <- function(data,
                         model.string,
                         n.adapt = 10000,
                         n.iter = 10000,
                         burnin = NULL,
                         thin = 10,
                         n.chains = 2,
                         quiet = FALSE,
                         .point = ggdist::Mode, 
                         .width = 0.95, 
                         .interval = ggdist::hdci,
                         ...) {
  
  if (class(data) != "isotopeData") stop("`data` must be an object of class 'isotopeData'")
  model_name <- class(model.string)
  if (!any(c("oneBaseline", "twoBaselines", "twoBaselinesFull") %in% model_name)) {
    stop("`model.string` must be a string containing the description of a 
         'oneBaseline', 'twoBaselines', or 'twoBaselinesFull' model")
  }
  
  model_name <- model_name[!model_name == "character"]
  
  
  # Note: not possible to set seed, need modification of `TPmodel()` source code
  fit <- tRophicPosition::TPmodel(data = data, 
                                  model.string = model.string, 
                                  n.adapt = n.adapt, 
                                  n.chains = n.chains, 
                                  quiet = quiet)
  
  variables <- if (model_name == "oneBaseline") {
    c("TP", "muDeltaN")
  } else {
    c("TP", "muDeltaN", "alpha") 
  }
  
  post <- tRophicPosition::posteriorTP(fit,
                                       variable.names = variables,
                                       n.iter = n.iter,
                                       burnin = burnin,
                                       thin = thin,
                                       quiet = quiet,
                                       ...)
  
  post_combined <- dplyr::bind_rows(lapply(post, as.data.frame))
  consumer <- ifelse(!is.null(attributes(data)$consumer),
                     attributes(data)$consumer,
                     NA)
  point_name <- gsub(".*::", "", tolower(quo_name(enquo(.point))))
  interval_name <- gsub(".*::", "", tolower(quo_name(enquo(.interval))))
  summary <- data.frame(consumer = consumer,
                        TP = post_combined$TP) |> 
    group_by(consumer) |> 
    ggdist::point_interval(TP, .point = .point, .width = .width, .interval = .interval) |> 
    mutate(.point = point_name, 
           .interval = interval_name)
  
  list(summary = summary,
       TP = post_combined$TP,
       alpha = if (model_name == "oneBaseline") NULL else post_combined$alpha,
       muDeltaN = post_combined$muDeltaN,
       samples = post)
}


#' Estimate trophic position
#'
#' @param sia_fish_corrected Output of [fish_lipid_correction()]
#' @param sia_baselines_corrected Output of [baselines_lipid_correction()]
#'
#' @return A list with 2 elements:
#'  - "models": a nested list of TP models output (one per species and TDF)
#'  - "TP_summary": tibble with summary of estimated trophic positions
estimate_TP <- function(sia_fish_corrected, sia_baselines_corrected) {
  
  # TDFs
  TDF_values <- list(
    "Post" = tRophicPosition::TDF(author = "Post",
                                  element = "both",
                                  seed = 123),
    "McCutchan" = tRophicPosition::TDF(author = "McCutchan", 
                                       element = "both", 
                                       type = "muscle",
                                       seed = 123)
  )
  
  # Prepare data set
  consumer_df <- sia_baselines_corrected |> 
    bind_rows(sia_fish_corrected) |> 
    as.data.frame()
  
  # Extract data by species
  consumer_list <- purrr::map(
    TDF_values,
    ~ tRophicPosition::extractIsotopeData(consumer_df, 
                                          b1 = "Bivalvia",
                                          b2 = "Gastropoda",
                                          baselineColumn = "baseline", 
                                          consumersColumn = "species",
                                          groupsColumn = NULL,
                                          d13C = "d13C_corrected", 
                                          d15N = "d15N",
                                          deltaC = .x$deltaC,
                                          deltaN = .x$deltaN)
  )
  
  # Run models
  model.string <- tRophicPosition::jagsBayesianModel(model = "twoBaselinesFull", 
                                                     TP = "dunif(2, 6)", # constrain TP between 2 and 6
                                                     lambda = 2)
  
  consumer_models <- purrr::map(
    consumer_list,
    ~ parallel::mclapply(.x,
                         fit_TP_model,
                         model.string = model.string,
                         n.adapt = 2000, n.iter = 2000,
                         burnin = 2000, thin = 1, n.chains = 4,
                         quiet = FALSE, 
                         mc.cores = 10)
  )
  
  # Extract TP
  df <- expand.grid(i = c("Post", "McCutchan"),
                    j = 1:length(consumer_models[[1]]))
  
  TP <- purrr::map(
    1:nrow(df),
    ~ consumer_models[[df$i[.x]]][[df$j[.x]]][["summary"]] |> 
      mutate(TDF_source = df$i[.x])) |> 
    bind_rows()
  
  return(list(models = consumer_models,
              TP_summary = TP)
         )
  
}

#' Plot TP correlation
#'
#' @param TP Output of [estimate_TP()]
#'
#' @returns A ggplot object
plot_TP_correlation <- function(TP) {
  
  plot <- TP$TP_summary |> 
    select(consumer, TP, TDF_source) |> 
    tidyr::pivot_wider(names_from = TDF_source,
                       values_from = TP) |> 
    ggplot(aes(x = Post, y = McCutchan)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(shape = 21, fill = "grey", alpha = 0.5) + 
    scale_x_continuous(breaks = seq(2, 4, by = 1)) +
    scale_y_continuous(breaks = seq(2, 4, by = 1)) + 
    ggtitle("Estimated TP") + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  
  return(plot)
}

#' Plot TP mode and HDI for both TDF sources
#'
#' @param TP Output of [estimate_TP()]
#'
#' @returns A list of ggplot objects
plot_TP <- function(TP) {
  
  # Split plots in two pages
  TP_long <- purrr::map(
    TP$models,
    function(i) {
      purrr::map(
        1:length(i),
        function(j) 
          data.frame(i[[j]][["TP"]]) |> 
          rlang::set_names(names(i)[j])
      ) |> 
        bind_cols() |> 
        tidyr::pivot_longer(cols = everything()) |> 
        mutate(id = as.numeric(as.factor(name)),
               page = ifelse(id < max(id)/2, 1, 2))
    }
  )
  
  plots <- purrr::map(
    TP_long,
    function(i) {
      purrr::map(
        1:2,
        function(j) {
          p <- i |> 
            filter(page == j) |> 
            ggplot(aes(y = name, x = value)) +
            ggdist::stat_pointinterval(point_interval = "mode_hdci", 
                                       .width = c(0.5, 0.95),
                                       fatten_point = 1.2,
                                       interval_size_range = c(0.4, 0.8)) +
            xlab("Trophic position") +
            theme_bw() + 
            scale_y_discrete(limits = function(x) rev(x)) +
            theme(axis.text.y = element_text(face = "italic", size = 7),
                  axis.title.y = element_blank())
          
          return(p)
          }
        )
      }
    )
  
  return(plots)
}
