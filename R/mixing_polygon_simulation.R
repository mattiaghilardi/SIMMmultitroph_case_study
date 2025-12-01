#' Mixing Polygon Simulation (2 isotopes)
#' 
#' Code adapted from Smith et al. 2013
#' Original code (v1.2) available at: http://www.famer.unsw.edu.au/software/polygon.html
#'
#' @param sources Data frame with four columns, mean and SD of delta13C and delta15N of sources. 
#' Always put 13C before 15N
#' @param mixture Data frame with two columns, delta13C and delta15N of consumers. 
#' Always put 13C before 15N
#' @param TDF Data frame with four columns, mean and SD TDFs of sources. 
#' Always put 13C before 15N
#' @param its Number of iterations
#' @param min_C Minimum delta13C
#' @param max_C Maximum delta13C
#' @param min_N Minimum delta15N
#' @param max_N Maximum delta15N
#' @param res Resolution
#' @param threshold Probability threshold used to determine point in polygon (default to 0.05)
#' @param seed Seed for reproducibility
#' @param facet Column name to be used for faceting isospace plots
#'
#' @return A list with six elements:
#'  - `variance_plot`: ggplot showing variance of polygon area during simulation
#'  - `probability_plot`: ggplot showing the proportion of iterations that each 
#'  consumer was inside mixing polygon
#'  - `isospace_plot_filled`: ggplot showing the isospace with mixing region, 
#'  consumers, and average enriched source signatures
#'  - `isospace_plot_bw`: ggplot showing isospace in black and white
#'  - `consumer_probabilities`: data frame of consumer probabilities. 
#'  Consumers on rows and iterations on columns. Last column is the average probability
#'  - `parameter_values`: data frame of parameter values for each iteration
#'
#' @importFrom sp point.in.polygon
#' @importFrom splancs areapl
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang set_names
#' @importFrom dplyr mutate
#' @importFrom ggtext element_markdown
#'
#' @references Smith, J.A., D. Mazumder, I. M. Suthers, M. D. Taylor (2013) 
#' To fit or not to fit: evaluating stable isotope mixing models using simulated mixing polygons. 
#' Methods in Ecology and Evolution 4 612-618, DOI: 10.1111/2041-210X.12048
#' 
#' @examples dontrun{
#' # Using example data from: http://www.famer.unsw.edu.au/software/polygon.html
#' sources <- read.table("Sources_example.csv",header=T,sep=",") #always put 13C(x) before 15N(y)
#' mixture <- read.table("Mixture_example.csv",header=T,sep=",") #some error is required for every value
#' TDF <-  read.table("TEF_example.csv", header=T,sep=",")
#' sim_mixpol <- mixing_polygon_simulation(sources, mixture, TDF)
#' sim_mixpol$isospace_plot_filled
#' }
mixing_polygon_simulation <- function(sources,
                                      mixture,
                                      TDF,
                                      its = 1500,  #specify the number of iterations ("its")
                                      min_C = -50,  #specify the dimensions and resolution for the mixing region figure
                                      max_C = -20,    #choose values outside the 95% mixing region
                                      min_N = -2,  
                                      max_N = 10,
                                      res = 250, #resolution of the mixing region figure; reducing this improves performance
                                      threshold = 0.05,
                                      seed = NULL,
                                      facet = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ##Now RUN the simulation 
  step_C <- (max_C - min_C) / (res - 1)
  step_N <- (max_N - min_N) / (res - 1)   
  C_g <- seq(min_C, max_C, by = step_C) #values must be in ascending order
  N_g <- seq(min_N, max_N, by = step_N) #values must be in ascending order
  mgrid <- function(a, b) {  #create a grid of values to test for P-I-P
    list(
      x = outer(b * 0, a, FUN = "+"),
      y = outer(b, a * 0, FUN = "+")
    )
  }
  m <- mgrid(C_g, N_g)
  Par_values <- array(0, c(its, nrow(sources) * 4 + 3))  #create files to store data
  p <- array(0, c(its, nrow(mixture)))
  mix_reg <- array(0, c(res, res))
  
  cli::cli_progress_bar("Iterations", total = its)
  for (i in 1:its) {    #run loops to generate source isotopic signatures, for iterations = 'its'
    v <- array(0, c(nrow(sources), 2))
    f <- array(0, c(nrow(TDF), 2))
    for (j in 1:nrow(sources)) {
      v[j, 1] <- rnorm(1, mean = sources[j, 1], sd = sources[j, 2])  #generate values from norm. dist. for d13C
      v[j, 2] <- rnorm(1, mean = sources[j, 3], sd = sources[j, 4])  #generate values from norm. dist. for d15N
      f[j, 1] <- rnorm(1, mean = TDF[j, 1], sd = TDF[j, 2])  #generate values from norm. dist. for d13C enrichment
      f[j, 2] <- rnorm(1, mean = TDF[j, 3], sd = TDF[j, 4])  #generate values from norm. dist. for d15N enrichment
    }
    V <- v + f
    hull <- chull(V)  #create a 2D convex hull from the enriched sources
    hull_a <- append(hull, hull[1])  #closes the polygon
    P <- sp::point.in.polygon(mixture[, 1], mixture[, 2], V[hull_a, 1], V[hull_a, 2]) #calculate P_I_P 
    P_n <- as.numeric(P)
    p[i,] <- P_n
    poly_a <- splancs::areapl(V[hull_a,])  #calculate polygon area, for evaluating the quantity of iterations
    m$y_f <- m$y[res:1,]  #flip y grid data to resemble axes (d13C=x, d15N=y) 
    m_r <- sp::point.in.polygon(m$x, m$y_f, V[hull_a, 1], V[hull_a, 2]) #calculate P-I-P for the mixing region
    m_r_s <- matrix(m_r, nrow = res, byrow = F)  #convert vector into square matrix
    m_r_s[m_r_s > 1] <- 1  #point.in.polygon can return '2' or '3'
    mix_reg <- mix_reg + m_r_s
    vals <- c(v[, 1], v[, 2], f[, 1], f[, 2], 0, 0, 0)  #concatenate values for this iteration
    Par_values[i,] <- vals  #store values
    Par_values[i, ncol(Par_values) - 2] <- poly_a
    Par_values[i, ncol(Par_values) - 1] <- i
    Par_values[i, ncol(Par_values)] <- var(Par_values[1:i, ncol(Par_values) - 2])
    cli::cli_progress_update()
  }
  
  # FIGURE 1: variance of polygon area during simulation
  p1 <- data.frame(Iterations = 1:its, 
                   Variance = Par_values[, ncol(Par_values)]) |> 
    ggplot(aes(x = Iterations, y = Variance)) + 
    geom_line() + 
    theme_bw()
  
  # FIGURE 2: proportion of iterations that each consumer was inside mixing polygon
  p2 <- data.frame(Consumer = as.factor(1:ncol(p)), 
                   prob = colSums(p) / its) |> 
    dplyr::mutate(inside = factor(ifelse(prob > threshold, "TRUE", "FALSE"),
                                  levels = c("TRUE", "FALSE"))) |>
    ggplot(aes(x = Consumer, y = prob, fill = inside)) +
    geom_col() +
    geom_hline(yintercept = threshold, linetype = "dashed") +
    scale_fill_manual(paste0("Probability > ", threshold),
                      values = c("darkgrey", "firebrick")) +
    ylab("Probability consumer inside mixing polygon") +
    ylim(0, 1) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  # FIGURE 3: mixing region, consumers, average enriched source signatures
  mix_reg <- mix_reg / its  #convert to 0-1 scale
  mix_reg[mix_reg == 0] <- NA #make the zeros white
  mix_regt <- t(mix_reg[ncol(mix_reg):1,])  #transpose matrix
  mix_regt <- mix_regt |> 
    as.data.frame(row.names = C_g) |> 
    rlang::set_names(N_g) |>  
    tibble::rownames_to_column("d13C") |> 
    tidyr::pivot_longer(cols = where(is.double), 
                        names_to = "d15N", 
                        values_to = "density") |> 
    dplyr::mutate(d13C = as.numeric(d13C), 
                  d15N = as.numeric(d15N))
  sources_TEF <- sources + TDF
  names(sources_TEF)[1] <- "d13C"
  names(sources_TEF)[3] <- "d15N"
  names(mixture)[1:2] <- c("d13C", "d15N")
  
  make_plot <- function(fill = TRUE, 
                        source_color = "white",
                        mixture_color = "black",
                        mixture_fill = "grey",
                        facet = NULL) {
    p <- mix_regt |> 
      ggplot(aes(x = d13C, y = d15N))
    
    if (fill) p <- p + geom_raster(aes(fill = density))
    
    breaks <- if(threshold < 0.05) 
      c(threshold, 0.05, seq(0.1, 1, by = 0.1))
    else 
      c(threshold, seq(0.1, 1, by = 0.1))
    
    p <- p +
      geom_contour(aes(z = density), 
                   breaks = breaks, 
                   color = "black") + 
      geom_point(data = sources_TEF,
                 color = source_color, 
                 shape = 4, size = 4) + 
      geom_point(data = mixture,
                 color = mixture_color,
                 fill = mixture_fill,
                 shape = 21,
                 alpha = 0.5) + 
      labs(x = "&delta;<sup>13</sup>C (&permil;)",
           y = "&delta;<sup>15</sup>N (&permil;)") +
      theme_bw() +
      theme(axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            panel.grid = element_blank())
    
    if (!is.null(facet)) p <- p + facet_wrap(vars(.data[[facet]]))
    
    p
  }
  
  p3 <- make_plot(facet = facet) +
    scale_fill_viridis_c()
  
  # FIGURE 4: Same as Fig. 3, but black and white
  p4 <- make_plot(fill = FALSE,
                  source_color = "black",
                  mixture_color = "darkgrey",
                  mixture_fill = "darkgrey",
                  facet = facet)
  
  # Table 1: Consumer probabilities
  p_df <- p |>
    rbind(colSums(p) / its) |>
    t() |>
    as.data.frame(row.names = rownames(mixture)) |>
    rlang::set_names(c(paste0("iteration.", 1:its), "probability"))
  
  # Table 1: Parameter values
  par_values_df <- Par_values |>
    as.data.frame() |>
    rlang::set_names(
      c(paste(c(rep("d13C", nrow(sources)),
                rep("d15N", nrow(sources)),
                rep("13C_TEF", nrow(sources)),
                rep("15N_TEF", nrow(sources))
      ), 
      rownames(sources), 
      sep = "_"),
      "Poly_Area", "Iteration", "Variance")
    )
  
  list(variance_plot = p1, 
       probability_plot = p2, 
       isospace_plot_filled = p3, 
       isospace_plot_bw = p4, 
       consumer_probabilities = p_df, 
       parameter_values = par_values_df)
  
}

#' Run mixing polygon simulations
#'
#' @param consumers Output of `prepare_consumer_data()`
#' @param sources Output of `prepare_source_data()`
#' @param TDF Output of `prepare_TDF_data()`
#' @param its Number of iterations
#' @param threshold Probability threshold used to determine point in polygon
#' @param filename Name of the file to save
#' @param width Width of the file to save
#' @param height Height of the file to save
#' @param units Unit of "width" and "height"
#' @param type File type (default to png)
#'
#' @return A list
run_mixing_polygon_simulation <- function(consumers, 
                                          sources, 
                                          TDF, 
                                          its = 5000,
                                          threshold = 0.05,
                                          filename,
                                          width = 18, 
                                          height = 18, 
                                          units = "cm", 
                                          type = "png") {
  
  # Prepare data
  source <- sources$sources$data |> 
    dplyr::group_by(source) |>
    dplyr::summarise(mean_d13C = mean(d13C, na.rm = TRUE), 
                     sd_d13C = sd(d13C, na.rm = TRUE),
                     mean_d15N = mean(d15N, na.rm = TRUE),
                     sd_d15N = sd(d15N, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    tibble::column_to_rownames("source")
  
  tdf <- TDF$data |> 
    dplyr::arrange(source) |>
    tibble::column_to_rownames("source")
  
  mixture <- consumers |> 
    purrr::map(~ .x |>
                 dplyr::select(id, d13C, d15N, trophic_guild) |>
                 tibble::column_to_rownames("id"))
  
  # Run simulations
  sim <- purrr::map(
    mixture,
    ~ mixing_polygon_simulation(sources = source, 
                                mixture = .x, 
                                TDF = tdf,
                                its = its, 
                                min_C = -25, max_C = -1,
                                min_N = -2, max_N = 8, 
                                res = 250,
                                threshold = threshold,
                                seed = 123,
                                facet = "trophic_guild")
    )
  
  # Add plot titles
  sim <- purrr::map(
    names(sim) |> 
      rlang::set_names(),
    function(.x) {
      sim[[.x]]$isospace_plot_filled <- sim[[.x]]$isospace_plot_filled +
        ggtitle(.x) +
        theme(plot.title = element_text(hjust = 0.5))
      
      sim[[.x]]$isospace_plot_bw <- sim[[.x]]$isospace_plot_bw +
        ggtitle(.x)+
        theme(plot.title = element_text(hjust = 0.5))
      
      return(sim[[.x]])
    }
  )
  
  # Move legend to empty panel
  sim <- purrr::map(
    sim,
    function(.x) {
      empty_panels <- .x$isospace_plot_filled |>
        cowplot::plot_to_gtable() |>
        gtable::gtable_filter("panel") |>
        with(setNames(grobs, layout$name)) |>
        purrr::keep(~identical(.x, zeroGrob())) |>
        names()
      
      if (length(empty_panels) > 0) {
        .x$isospace_plot_filled <- .x$isospace_plot_filled |>
          lemon::reposition_legend(position = "center", 
                                   panel = empty_panels) |> 
          ggplotify::as.ggplot()
      }
      
      return(.x)
    }
  )
  
  # Save plots
  purrr::map(
    names(sim),
    ~ ggsave(here::here("output", "figures", paste0(filename, "_", .x , ".", type)),
             sim[[.x]]$isospace_plot_filled, width = width, height = height, units = units)
  )
  
  # Plot histogram of probabilities across methods
  purrr::map(sim, 
             ~ .x$consumer_probabilities |> 
               select(probability)) |> 
    bind_rows(.id = "TDF_source") |> 
    ggplot(aes(x = probability)) + 
    geom_histogram(bins = 30, color = "black", fill = "grey") + 
    labs(x = "Probability consumer inside mixing polygon", 
         y = "Number of individual consumers") + 
    facet_grid(~TDF_source) + 
    theme_bw()
  
  ggsave(here::here("output", "figures", paste0("probability_inside_mixing_polygon.", type)),
         width = 18, height = 10, units = units)
  
  return(sim)
}
