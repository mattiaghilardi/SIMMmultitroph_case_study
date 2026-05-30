# Functions to prepare data as required by MixSIAR

#' Prepare source data for TP and MixSIAR and save csv file in "derived_data"
#'
#' @param sia_sources Raw sia_sources data
#' 
#' @return Path to csv file
prepare_source_data <- function(sia_sources) {
  
  # Check carbon content
  max(sia_sources$C_percent)
  # 37.47349 -> below 40%, no need for lipid-correction
  
  # Combine green and brown algae 
  sources <- sia_sources |>
    mutate(source = ifelse(source %in% c("Green algae", "Brown algae"),
                           "Green-Brown algae",
                           source))
  
  # Reorder columns for MixSIAR, first column must be the source names
  sources <- sources |> 
    select(source, sample_id, sample_type, species, 
           d13C, d15N, C_percent, N_percent, CN_ratio)
  
  # Save csv file in "derived_data"
  path <- "derived_data/sources.csv"
  readr::write_csv(sources, path)
  
  return(path)
}

#' Prepare TDF data for MixSIAR and save csv files in "derived_data"
#'
#' @param sources Output of [prepare_source_data()]
#'
#' @return Path to csv file
prepare_TDF_data <- function(sources) {
  
  # Extract source names
  source_names <- sources |> 
    read_csv() |> 
    arrange(source) |> 
    pull(source) |> 
    unique()
  
  # Create data frame
  tdf <- data.frame(source = source_names,
                    Meand13C = 0,
                    SDd13C = 0,
                    Meand15N = 0,
                    SDd15N = 0)
  
  # Save csv file in "derived_data"
  path <- "derived_data/TDF.csv"
  readr::write_csv(tdf, path)
  
  return(path)
}

#' Prepare consumer data
#' 
#' Add trophic guilds and rescale isotopic values to `TP = 1`
#'
#' @param sia_fish_corrected Output of [fish_lipid_correction()]
#' @param trophic_guilds Raw trophic_guilds data
#' @param TP Output of [estimate_TP()]
#'
#' @return A list of consumers data frames
prepare_consumer_data <- function(sia_fish_corrected,
                                  trophic_guilds,
                                  TP) {
  
  # Add trophic guilds
  consumers <- sia_fish_corrected |>
    left_join(trophic_guilds, by = "species") |> 
    # Add column with log(TL) for MixSIAR
    mutate(log_tl = log(tl))
  
  # Check number of observations per guild
  consumers |> count(trophic_guild)
  # Only 18 observations for "invertivores-pelagic" -> check number of species
  consumers |> 
    filter(trophic_guild == "invertivores-pelagic") |> 
    pull(species) |> 
    unique()
  # Only two: "Odonus niger" and "Chaetodon trichrous" -> assign as "planktivores"
  consumers <- consumers |>
    mutate(trophic_guild = ifelse(trophic_guild == "invertivores-pelagic",
                                  "planktivores", 
                                  trophic_guild))
  
  # Mean TDFs
  TDF_summary <- data.frame(TDF_source = c("Post", "McCutchan"),
                            deltaC_mean = c(0.39, 1.3),
                            deltaN_mean = c(3.4, 2.9))
  
  # Rescale consumers to TP = 1 for each TDF source used
  consumers_list <- purrr::map(
    TDF_summary$TDF_source |> 
      rlang::set_names(),
    ~ consumers |>
      rename(d13C_raw = d13C,
             d15N_raw = d15N) |> 
      left_join(TP$TP_summary |> 
                  filter(TDF_source == .x) |> 
                  select(consumer, TP, TDF_source),
                by = c("species" = "consumer")) |>
      left_join(TDF_summary,
                by = "TDF_source") |>
      mutate(d13C = d13C_corrected - deltaC_mean * (TP - 1),
             d15N = d15N_raw - deltaN_mean * (TP - 1))
  )
  
  consumers_list
}

#' Remove consumer outliers and save csv files in "derived_data"
#' 
#' @param consumers Output of [prepare_consumer_data()]
#' @param consumer_outliers Output of [identify_outliers()]
#'
#' @return Paths to csv files
remove_outliers <- function(consumers, consumer_outliers) {
  
  # Remove outliers
  consumers_clean <- purrr::map(
    consumers,
    ~ .x |> 
      filter(! sample_id %in% (consumer_outliers[[unique(.x$TDF_source)]]$sample_id))
    )
  
  # Save consumer data by trophic guild
  paths <- purrr::map(
    consumers_clean,
    function(.x) {
      # Split by trophic guild
      sia_fish_by_guild <- split(.x,
                                 as.factor(.x$trophic_guild))
      # Make paths
      TDF_source <- unique(.x$TDF_source)
      paths <- names(sia_fish_by_guild) |>
        rlang::set_names() |>
        purrr::map(~ paste0("derived_data/consumers_", .x, "_", TDF_source, ".csv"))
      # Save csv files
      purrr::map(names(sia_fish_by_guild),
                 ~ readr::write_csv(sia_fish_by_guild[[.x]], paths[[.x]]))
      
      return(paths)
    }
  )
  
  return(paths)
}
