library(targets)
library(ggplot2)

# Initial isospace
tar_read(initial_isospace)
ggsave("output/figures/initial_isospace.png",
       width = 14, height = 12, units = "cm")

# Overlap between source isotopic signatures
tar_read(source_overlap)$plot
ggsave("output/figures/source_isotopic_overlap.png",
       width = 18, height = 20, units = "cm")

# Fish isotopic signatures across years
tar_read(fish_isotope_check)$plot
ggsave("output/figures/isotope_variation_across_years.png",
       width = 18, height = 20, units = "cm")

# Source contribution to baselines
tar_read(baselines_diet)$plot
ggsave("output/figures/mixsiar_baselines.png",
       width = 18, height = 12, units = "cm")

# Correlation between trophic positions
tar_read(TP_corr_plot)
ggsave("output/figures/TP_comparison.png",
       width = 10, height = 10, units = "cm")

# Species-specific trophic position
tar_load(TP_plots)
purrr::map2(c("Post", "Post", "McCutchan", "McCutchan"),
            c(1, 2, 1, 2),
            ~ ggsave(paste0("output/figures/", paste0("TP_", .x, "_", .y, ".png")),
                     TP_plots[[.x]][[.y]],
                     height = 22, width = 14, units = "cm")
)

# Rescaled isospace
tar_read(rescaled_isospace)
ggsave("output/figures/rescaled_isospace.png",
       width = 18, height = 12, units = "cm")

# Probability in mixing polygon
tar_load(simulation_mixing_polygon)
purrr::map(
  names(simulation_mixing_polygon),
  ~ ggsave(paste0("output/figures/mixing_polygon_simulation_", .x , ".png"),
           simulation_mixing_polygon[[.x]]$isospace_plot_filled,
           width = 18, height = 18, units = "cm")
)

# Probability in mixing polygon
tar_read(prob_in_mixing_polygon)
ggsave("output/figures/probability_inside_mixing_polygon.png",
       width = 18, height = 10, units = "cm")

# Plot consumer outliers
tar_load(consumer_outliers_plot)
purrr::map(names(consumer_outliers_plot),
           ~ ggsave(paste0("output/figures/mixing_polygon_outliers_", .x, ".png"),
                    consumer_outliers_plot[[.x]],
                    width = 18, height = 18, units = "cm")
)

# Dirichlet priors
tar_read(prior_plot)
ggsave("output/figures/alpha_priors.png",
       width = 18, height = 20, units = "cm")

# Final isospace by trophic guild
tar_read(MixSIAR_isospace_Post)
ggsave("output/figures/isospace_Post.png",
       width = 18, height = 20, units = "cm")
tar_read(MixSIAR_isospace_McCutchan)
ggsave("output/figures/isospace_McCutchan.png",
       width = 18, height = 20, units = "cm")

# Correlation between mixing model estimates
tar_read(MixSIAR_comparison_rel_contributions)
ggsave("output/figures/mixing_model_comparison.png",
       width = 18, height = 18, units = "cm")

# Isospace and relative source contributions
tar_read(MixSIAR_isospace_and_rel_contributions_Post)
ggsave("output/figures/source_contribution_by_trophic_guild_Post.png",
       width = 18, height = 22, units = "cm")
tar_read(MixSIAR_isospace_and_rel_contributions_McCutchan)
ggsave("output/figures/source_contribution_by_trophic_guild_McCutchan.png",
       width = 18, height = 22, units = "cm")
