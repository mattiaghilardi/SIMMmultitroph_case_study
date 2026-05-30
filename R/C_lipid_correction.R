#' Lipid-normalisation of d13C for aquatic organisms
#'
#' @param d13C Measured \u03b4^13^C value
#' @param CN_ratio C:N ratio of the sample
#' @param model String; the lipid-normalisation model to use. 
#' There are three options:
#' - `Post`: linear model from Post et al 2007
#' - `KP`: non-linear model from Kiljunen et al 2006, 
#'         with lipid percent model from Post et al 2007
#' - `KMM`: non-linear model from Kiljunen et al 2006, 
#'          with lipid percent model from McConnaughey and McRoy 1979
#'
#' @returns Lipid-normalised value
C_lipid_correction <- function(d13C, CN_ratio, model = c("Post", "KP", "KMM")) {
  
  model <- match.arg(model)
  
  if (model == "Post") {
    # Model from Post et al 2007 for aquatic organisms
    d13C_corrected <- d13C + 0.99 * CN_ratio - 3.32
  } else {
    
    # Model from Kiljunen et al 2006 with different % lipid model
    if (model == "KP") {
      # % lipid model from Post et al 2007 for aquatic organisms
      L <- -20.54 + 7.24 * CN_ratio
    } else { # model == "KMM"
      # % lipid model from McConnaughey and McRoy 1979
      L <- 93 / (1 + (0.246 * CN_ratio - 0.775)^-1)
    }
    
    # D and I constants from Kiljunen et al 2006: 
    # D = 7.018 ± 0.263
    # I = 0.048 ± 0.013
    D <- 7.018
    I <- 0.048
    
    # Model from Kiljunen et al 2006, originally from McConnaughey anf McRoy 1979
    d13C_corrected <- d13C + D * (I + (3.90 / (1 + 287 / L)))
  }
  
  return(d13C_corrected)
}

#' Lipid-normalisation of baselines d13C
#'
#' @param sia_baselines Raw sia_baselines data
#'
#' @returns The input data with an additional column "d13C_corrected"
baselines_lipid_correction <- function(sia_baselines) {
  
  # Lipid correction from Post et al 2007 for aquatic organisms 
  data_corrected <- sia_baselines |> 
    dplyr::mutate(
      d13C_corrected = ifelse(!is.na(CN_ratio) & CN_ratio >= 3.5,
                              C_lipid_correction(d13C, CN_ratio, model = "Post"),
                              d13C)
      )
  
  return(data_corrected)
}

#' Lipid-normalisation of fish d13C
#'
#' @param sia_fish Raw sia_fish data
#'
#' @returns The input data with an additional column "d13C_corrected"
fish_lipid_correction <- function(sia_fish) {
  
  # Lipid correction recommended by Skinner et al 2016:
  # Lipid normalisation model from Kiljunen et al 2006
  # with % lipid model from Post et al 2007 for aquatic organisms 
  data_corrected <- sia_fish |> 
    dplyr::mutate(
      d13C_corrected = ifelse(!is.na(CN_ratio) & CN_ratio >= 3.5,
                              C_lipid_correction(d13C, CN_ratio, model = "KP"),
                              d13C)
    )

  return(data_corrected)
}
