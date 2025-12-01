# Functions to prepare data as required by MixSIAR

#' Prepare source data for TP and MixSIAR and save csv files in "derived_data"
#'
#' @param sia_sources Raw sia_sources data
#' 
#' @return A list with data and paths to csv files
prepare_source_data <- function(sia_sources) {
  
  # Remove 3 observations with unmeasured d15N
  sia_sources <- sia_sources |> 
    filter(d15N != 0)
  
  # check outliers
  sia_sources |> ggplot(aes(x = names, y = d13C)) + geom_boxplot()
  sia_sources |> ggplot(aes(x = names, y = d15N)) + geom_boxplot()
  # Few outliers but most are close to the data
  # Check on biplot
  df <- sia_sources |>
    group_by(names) |>
    mutate(n_low = quantile(d15N, probs = 0.25) - 1.5*IQR(d15N),
           n_upp = quantile(d15N, probs = 0.75) + 1.5*IQR(d15N),
           c_low = quantile(d13C, probs = 0.25) - 1.5*IQR(d13C),
           c_upp = quantile(d13C, probs = 0.75) + 1.5*IQR(d13C),
           outlier = factor(ifelse(d13C > c_upp |
                                     d13C < c_low |
                                     d15N > n_upp |
                                     d15N < n_low,
                                   1, 0))) |>
    ungroup()
  
  df |> 
    ggplot(aes(x = d13C, y = d15N, color = outlier)) + 
    geom_point() + 
    facet_wrap(~names) + 
    theme_bw()
  
  # Some outliers are still close to other data points and may represent natural variation
  # But four outliers are quite extreme, one for POM, one for Bivalvia with relatively low d15N,
  # and two Gasteropoda with relatively high d15N
  df <- df |>
    mutate(true_outlier = factor(case_when(names == "Bivalvia" & d15N < 4 ~ 1,
                                           names == "POM" & d15N < 0 ~ 1,
                                           names == "Gastropoda" & d15N > 8 ~ 1,
                                           .default = 0)))
  
  df |> 
    ggplot(aes(x = d13C, y = d15N, fill = true_outlier)) + 
    geom_point(shape = 21, alpha = 0.8) + 
    facet_wrap(~names) + 
    labs(x = "&delta;<sup>13</sup>C (&permil;)",
         y = "&delta;<sup>15</sup>N (&permil;)",
         fill = "Outlier") +
    scale_fill_manual(values = c("white", "red")) +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  ggsave(here::here("output", "figures", "source_inverts_outliers.png"), 
         width = 18, height = 16, units = "cm")
  
  # Remove them
  df2 <- df |>
    filter(true_outlier == 0)
  
  # Check source overlap
  summary_sources <- df2 |>
    group_by(type, names) |>
    summarise(n = n(),
              mean_d13C = mean(d13C, na.rm = TRUE), 
              sd_d13C = sd(d13C, na.rm = TRUE),
              mean_d15N = mean(d15N, na.rm = TRUE),
              sd_d15N = sd(d15N, na.rm = TRUE)) |>
    ungroup()
  
  p1 <- ggplot(summary_sources |>
                 filter(type != "Invertebrates"),
               aes(x = mean_d13C, 
                   y = mean_d15N,
                   xmin = mean_d13C - sd_d13C,
                   xmax = mean_d13C + sd_d13C,
                   ymin = mean_d15N - sd_d15N,
                   ymax = mean_d15N + sd_d15N,
                   color = names)) +
    geom_linerange(orientation = "x") +
    geom_pointrange(orientation = "y") +
    xlab("&delta;<sup>13</sup>C (&permil;)") +
    ylab("&delta;<sup>15</sup>N (&permil;)") +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          legend.title = element_blank())
  
  # models
  fit_n <- brms::brm(d15N ~ 0 + names, 
                     data = df2 |> filter(type != "Invertebrates"),
                     cores = 4,
                     backend = "cmdstan",
                     seed = 123)
  fit_c <- brms::brm(d13C ~ 0 + names, 
                     data = df2 |> filter(type != "Invertebrates"),
                     cores = 4,
                     backend = "cmdstan",
                     seed = 123)
  brms::pp_check(fit_n)
  brms::pp_check(fit_c)
  
  p2 <- purrr::map2(
    .x = list(fit_n, fit_c),
    .y = c("&delta;<sup>15</sup>N", "&delta;<sup>13</sup>C"),
    ~ brms::posterior_epred(.x, 
                            newdata = data.frame(names = df2 |> 
                                                   filter(type != "Invertebrates") |> 
                                                   pull(names) |> 
                                                   unique()),
                            re_formula = NA) |> 
      as.data.frame() |> 
      rlang::set_names(df2 |> 
                         filter(type != "Invertebrates") |> 
                         pull(names) |> 
                         unique()) |> 
      tibble::rownames_to_column(".draw") |> 
      tidyr::pivot_longer(cols = -.draw, names_to = "source") |> 
      tidybayes::compare_levels(value, by = source) |> 
      mutate(isotope = .y, 
             # for setting threshold at 0.75
             # lower = quantile(value, probs = 0.25), 
             # upper = quantile(value, probs = 0.75), 
             # overlap = ifelse(lower > 0 | upper < 0, 0, 1),
             # overlap = factor(overlap,
             #                  levels = c(0, 1),
             #                  labels = c(">25%", "<25%")),
             above = ifelse(value > 0, 1, 0), 
             prob = sum(above) / n(), 
             prob = ifelse(prob >= 0.5, prob, 1 - prob))) |> 
    bind_rows() |> 
    ggplot(aes(x = value, y = source, fill = prob)) +
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
    # scale_fill_manual("Posterior distribution overlapping 0", 
    #                   values = c("deepskyblue", "grey")) +
    theme_bw() + 
    xlab("Difference in isotopic signature (&permil;)") +
    scale_y_discrete(position = "right") +
    theme(strip.text.x = ggtext::element_markdown(face = "bold"), 
          axis.title.y = element_blank(),
          axis.title.x = ggtext::element_markdown(),
          legend.position = "bottom")
  
  # final plot
  p1 + p2 + 
    patchwork::plot_layout(
      design = "
      #AAAAAAAAAAAAAAA#
      #AAAAAAAAAAAAAAA#
      BBBBBBBBBBBBBBBBB
      BBBBBBBBBBBBBBBBB
      BBBBBBBBBBBBBBBBB
      ")
  
  ggsave(here::here("output", "figures", "source_isotopic_overlap.png"), 
         width = 18, height = 20, units = "cm")
  
  # All algal groups have similar isotopic composition
  # We can exclude turf as it has only one data point and 
  # combine green and brown algae as "Macroalgae"
  df3 <- df2 |> 
    filter(names != "Turf") |>
    mutate(source = ifelse(type == "Algae", "Macroalgae", names))
  
  # save cleaned files for sources and invertebrates
  paths <- list(sources = here::here("derived_data", "sources.csv"),
                invertebrates = here::here("derived_data", "invertebrates.csv"))
  
  source_data <- df3 |>
    filter(type != "Invertebrates") |>
    select(source, names, id, d13C, d15N)
  
  write_csv(source_data, paths$sources)
  
  inverts_data <- df3 |>
    filter(type == "Invertebrates") |>
    select(id, type, taxon = names, d13C, d15N)
  
  write_csv(inverts_data, paths$invertebrates)
  
  list(sources = list(data = source_data,
                      path = paths$sources),
       invertebrates = list(data = inverts_data,
                            path = paths$invertebrates))
}

#' Prepare consumer data 
#'
#' @param sia_fish Raw sia_fish data
#' @param trophic_guilds Raw trophic_guilds data
#' @param TP Data frame with trophic positions.
#' Output of [estimate_TP()]
#'
#' @return A list of consumers data frames
prepare_consumer_data <- function(sia_fish,
                                  trophic_guilds,
                                  TP) {
  
  # Add trophic guilds
  consumers <- sia_fish |>
    left_join(trophic_guilds, by = "species")
  
  # Add column with log(TL) for MixSIAR
  consumers <- consumers |> mutate(log_tl = log(tl))
  
  # Check number of observations per guild
  consumers |> count(trophic_guild)
  consumers |> filter(trophic_guild == "invertivores-pelagic") |> pull(species) |> unique()
  # only "Odonus niger" and "Chaetodon trichrous"
  # assign as planktivores
  consumers <- consumers |>
    mutate(trophic_guild = ifelse(trophic_guild == "invertivores-pelagic",
                                  "planktivores", 
                                  trophic_guild)) |>
    select(id, year, species, family, 
           tl, log_tl, sl, weight, trophic_guild, 
           d13C, d15N)
  
  # TDF_summary <- read_csv(here::here("output", "tables", "TDF_values.csv"))
  TDF_summary <- data.frame(TDF_source = c("Post", "McCutchan"),
                            deltaC_mean = c(0.39, 1.3),
                            deltaN_mean = c(3.4, 2.9))
  
  consumers_list <- purrr::map(
    TP$TDF_source |> 
      unique() |> 
      rlang::set_names(),
    ~ consumers |>
      rename(d13C_raw = d13C,
             d15N_raw = d15N) |> 
      left_join(TP |> 
                  filter(TDF_source == .x) |> 
                  select(species = consumer, TP, TDF_source),
                by = "species") |>
      left_join(TDF_summary,
                by = "TDF_source") |>
      mutate(d13C = d13C_raw - deltaC_mean * (TP - 1),
             d15N = d15N_raw - deltaN_mean * (TP - 1))
  )
  
  consumers_list
}

#' Prepare TDF data for MixSIAR and save csv files in "derived_data"
#'
#' @param sources Output of [prepare_source_data()]
#'
#' @return Path to csv file
prepare_TDF_data <- function(sources) {
  
  source_names <- sources$sources$data |> 
    pull(source) |> 
    unique()
  
  tdf <- data.frame(source = source_names,
                    Meand13C = 0,
                    SDd13C = 0,
                    Meand15N = 0,
                    SDd15N = 0)
  
  path <- here::here("derived_data", "TDF.csv")
  write_csv(tdf, path)
  
  list(data = tdf, 
       path = path)
}

#' Remove consumer outliers and save csv files in "derived_data"
#' 
#' @param consumers Output of [prepare_consumer_data()]
#' @param simulation_mixing_polygon Output of [run_mixing_polygon_simulation()]
#'
#' @return Paths to csv files
remove_outliers <- function(consumers, simulation_mixing_polygon) {
  
  # Get ID of outliers (probability < 0.01)
  id_outliers <- purrr::map(
    simulation_mixing_polygon,
    ~ .x[["consumer_probabilities"]] |>
      tibble::rownames_to_column("id") |>
      filter(probability < 0.01) |>
      select(id, probability))
  
  # # Which has less outliers
  # n_outliers <- purrr::map(id_outliers, ~ nrow(.x))
  # min_outliers <- which.min(n_outliers |> purrr::as_vector())
  
  purrr::map(
    names(id_outliers),
    function(i) {
      # Data for plot
      df_outliers <- consumers[[i]] |> 
        left_join(id_outliers[[i]]) |> 
        mutate(outlier = ifelse(is.na(probability), 0, 1)) |> 
        group_by(species) |> 
        filter(max(outlier) == 1)
      
      # Compute number of facet columns
      n_cols <- df_outliers$species |> unique() |> length() |> sqrt() |> round()
      
      # Plot
      ggplot(df_outliers,
             aes(x = d13C_raw, y = d15N_raw, fill = factor(outlier))) + 
        geom_point(shape = 21) + 
        facet_wrap(~species, ncol = n_cols) + 
        scale_fill_manual(values = c("white", "red")) + 
        labs(x = "&delta;<sup>13</sup>C (&permil;)",
             y = "&delta;<sup>15</sup>N (&permil;)", 
             fill = "Outside mixing polygon") + 
        theme_bw() + 
        theme(panel.grid.minor = element_blank(), 
              strip.text = element_text(face = "italic", size = 8), 
              axis.title.x = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              legend.position = "top")
      
      ggsave(here::here("output", "figures", paste0("mixing_polygon_outliers_", i, ".png")), 
             width = 18, height = 18, units = "cm")
    }
  )
  
  # Remove outliers
  consumers_clean <- purrr::map(
    consumers,
    ~ .x |> 
      filter(!id %in% (id_outliers[[unique(.x$TDF_source)]]$id)))
  
  # Save consumer data by trophic guild
  paths <- purrr::map(
    consumers_clean,
    function(.x) {
      
      sia_fish_by_guild <- split(.x,
                                 as.factor(.x$trophic_guild))
      TDF_source <- unique(.x$TDF_source)
      paths <- names(sia_fish_by_guild) |>
        rlang::set_names() |>
        purrr::map(~ here::here("derived_data", paste0("consumers_", .x, "_", TDF_source, ".csv")))
      
      purrr::map(names(sia_fish_by_guild),
                 ~ write_csv(sia_fish_by_guild[[.x]], paths[[.x]]))
      
      return(paths)
    }
  )
  
  return(list(data = consumers_clean,
              path = paths,
              outliers = id_outliers))
}
