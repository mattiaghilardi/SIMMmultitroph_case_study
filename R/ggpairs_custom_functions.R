# custom functions for GGally::ggpairs

# diagonal: density and label
diag_custom <- function(data, mapping, xtext = -Inf, ytext = Inf, ...) {
  
  xtext <- ifelse(xtext %in% c("center", "centre"),
                  mean(range(data[[as_label(mapping$x)]])),
                  xtext)
  ggplot(data = data, mapping = mapping, ...) +
    geom_density() +
    annotate("text", label = as_label(mapping$x), 
             x = xtext, y = ytext, ...) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
}

# points and identity line
point_custom <- function(data, mapping, linecolor = "black", ...) {
  
  ggplot(data = data, mapping = mapping, ...) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = linecolor) +
    geom_point(...) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

# bayesian correlation
cor_custom <- function(data, mapping, main_model, sub_models = NULL,
                       digits = 3, wrap_interval = TRUE, ...) {
  
  a <- gsub("\n.*", "", as_label(mapping$x))
  b <- gsub("\n.*", "", as_label(mapping$y))
  
  extract_cor <- function(model) {
    brms::as_draws_df(model) |> 
      select(starts_with("rescor")) |> 
      tibble::rownames_to_column(".draw") |> 
      tidyr::pivot_longer(cols = -.draw) |> 
      mutate(name = gsub("rescor__",  "", name)) |>
      tidyr::separate_wider_delim(cols = name, delim = "__", names = c("TPa", "TPb")) |> 
      filter((TPa == a & TPb == b) | (TPa == b & TPb == a)) |> 
      group_by(TPa, TPb) |> 
      ggdist::median_qi(.width = 0.95) |> 
      mutate(label = paste0("\u03C1 : ",
                            round(value, digits), 
                            ifelse(wrap_interval, "\n[", " ["),
                            round(.lower, digits), 
                            "-", 
                            round(.upper, digits), 
                            "]")) |> 
      pull(label)
  }
  
  label <- extract_cor(main_model)
  
  if (is.null(sub_models)) {
    GGally::ggally_text(label = label, ...) +
      theme_bw() +
      theme(panel.grid = element_blank())
  } else {
    sub_labels <- purrr::map(1:nrow(sub_models),
                             ~ paste0(sub_models[.x, 1, drop = TRUE], 
                                      "\n", 
                                      extract_cor(sub_models[.x, 2, drop = TRUE][[1]]))) |> 
      purrr::list_simplify()
    
    yPos <- seq(from = 0.9, to = 0.1, length.out = nrow(sub_models) + 1)
    yPos <- yPos[-1]
    
    sub_models$label <- sub_labels
    sub_models$yPos <- yPos
    group_var <- names(sub_models)[1]
    
    GGally::ggally_text(label = label, yP = 0.9, ...) +
      geom_text(data = sub_models, 
                aes(x = 0.5, y = yPos, label = label, color = .data[[group_var]]),
                ...) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
    
}
