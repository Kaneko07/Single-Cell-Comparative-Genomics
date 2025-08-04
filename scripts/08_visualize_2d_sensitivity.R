suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(viridis)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript 08_visualize_2d_sensitivity_enhanced.R <2d_results.csv> <output_dir> <snp_changepoint> <sag_changepoint> <cluster_template.csv>", call. = FALSE)
}

sensitivity_file <- args[1] 
output_dir <- args[2]
snp_changepoint <- as.numeric(args[3])
sag_changepoint <- as.numeric(args[4])
cluster_file <- args[5]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("INFO: Starting enhanced 2D visualization for dual-axis validation...\n")

plot_data <- read_csv(cluster_file, show_col_types = FALSE)

numeric_cols <- c("snp_threshold", "sag_threshold", "sag_count", "final_snp_count", "n_content_mean", "cluster_count")
for (col in numeric_cols) {
  if (col %in% colnames(plot_data)) {
    plot_data[[col]] <- as.numeric(plot_data[[col]])
  }
}

plot_data <- plot_data %>%
  mutate(
    is_snp_changepoint = (snp_threshold == snp_changepoint),
    is_sag_changepoint = (sag_threshold == sag_changepoint),
    is_both_changepoint = (is_snp_changepoint & is_sag_changepoint),
    changepoint_type = case_when(
      is_both_changepoint ~ "Both",
      is_snp_changepoint ~ "SNP Only",
      is_sag_changepoint ~ "SAG Only",
      TRUE ~ "Neither"
    )
  )

cat(paste("Data loaded:", nrow(plot_data), "combinations\n"))

theme_publication <- theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95")
  )

cat("INFO: Creating cluster stability analysis plot...\n")

snp1_data <- plot_data %>% 
  filter(snp_threshold == snp_changepoint) %>%
  arrange(sag_threshold)

if (nrow(snp1_data) > 0) {
  p_stability <- ggplot(snp1_data, aes(x = sag_threshold, y = cluster_count)) +
    geom_line(color = "steelblue", size = 1.2, alpha = 0.8) +
    geom_point(size = 4, alpha = 0.9, aes(color = factor(cluster_count))) +
    
    geom_vline(xintercept = sag_changepoint, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = sag_changepoint, y = max(snp1_data$cluster_count) * 0.9,
             label = paste("Statistical\nChangepoint\n(SAG =", sag_changepoint, ")"), 
             color = "red", size = 3.5, hjust = 0.5) +
    
    geom_rect(data = snp1_data %>% filter(cluster_count == 12) %>% 
                summarise(xmin = min(sag_threshold) - 0.4, xmax = max(sag_threshold) + 0.4,
                          ymin = 11.5, ymax = 12.5),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "lightgreen", alpha = 0.3, inherit.aes = FALSE) +
    
    geom_point(data = snp1_data %>% filter(sag_threshold == 14),
               aes(x = sag_threshold, y = cluster_count), 
               color = "darkgreen", size = 6, shape = 17) +
    annotate("text", x = 14, y = 12.5,
             label = "Proposed\nOptimal\n(SAG = 14)", 
             color = "darkgreen", size = 3.5, hjust = 0.5) +
    
    scale_color_viridis_d(name = "Cluster\nCount") +
    labs(
      title = "Cluster Stability Analysis",
      subtitle = paste("SNP Threshold =", snp_changepoint, "(Statistical Changepoint)"),
      x = "SAG Threshold",
      y = "Number of Phylogenetic Clusters"
    ) +
    theme_publication +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "01_cluster_stability_analysis.png"), p_stability, 
         width = 12, height = 8, dpi = 300)
  
  p_data_cluster <- ggplot(snp1_data, aes(x = sag_count + final_snp_count, y = cluster_count)) +
    geom_smooth(method = "loess", se = TRUE, color = "grey70", alpha = 0.3) +
    geom_point(size = 4, aes(color = factor(sag_threshold))) +
    
    geom_point(data = snp1_data %>% filter(sag_threshold %in% c(13, 14)),
               size = 6, shape = 21, stroke = 2, fill = "white") +
    
    geom_text(data = snp1_data %>% filter(sag_threshold == 13),
              aes(label = paste("CP (SAG=13)\nClusters:", cluster_count)),
              vjust = -1.5, hjust = 0.5, size = 3.5, color = "red") +
    
    geom_text(data = snp1_data %>% filter(sag_threshold == 14),
              aes(label = paste("Proposed (SAG=14)\nClusters:", cluster_count)),
              vjust = 1.8, hjust = 0.5, size = 3.5, color = "darkgreen") +
    
    scale_color_viridis_d(name = "SAG\nThreshold") +
    labs(
      title = "Data Volume vs Cluster Count Relationship",
      subtitle = "Trade-off between data quantity and phylogenetic resolution",
      x = "Total Data Volume (SAG + SNP Count)",
      y = "Number of Phylogenetic Clusters"
    ) +
    theme_publication
  
  ggsave(file.path(output_dir, "02_data_volume_vs_clusters.png"), p_data_cluster, 
         width = 10, height = 8, dpi = 300)
}

cat("INFO: Creating dual-axis validation summary...\n")

validation_data <- plot_data %>% 
  filter(snp_threshold == snp_changepoint, sag_threshold >= 10, sag_threshold <= 16) %>%
  mutate(
    analysis_type = case_when(
      sag_threshold == sag_changepoint ~ "Statistical Changepoint",
      sag_threshold == 14 ~ "Proposed Optimal",
      TRUE ~ "Alternative"
    ),
    quality_score = 1 / n_content_mean,
    data_score = (sag_count + final_snp_count) / 1000,
    cluster_score = cluster_count / max(cluster_count)
  )

p_multi_comparison <- validation_data %>%
  select(sag_threshold, analysis_type, sag_count, final_snp_count, n_content_mean, cluster_count) %>%
  pivot_longer(cols = c(sag_count, final_snp_count, n_content_mean, cluster_count),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = case_when(
      metric == "sag_count" ~ "SAG Count",
      metric == "final_snp_count" ~ "SNP Count", 
      metric == "n_content_mean" ~ "N Content (lower is better)",
      metric == "cluster_count" ~ "Cluster Count"
    ),
    metric_label = factor(metric_label, levels = c("SAG Count", "SNP Count", "N Content (lower is better)", "Cluster Count"))
  ) %>%
  ggplot(aes(x = factor(sag_threshold), y = value)) +
  geom_col(aes(fill = analysis_type), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(value, 3), color = analysis_type), 
            vjust = -0.3, size = 3, fontface = "bold") +
  
  facet_wrap(~metric_label, scales = "free_y", ncol = 2) +
  
  scale_fill_manual(values = c("Statistical Changepoint" = "#e31a1c", 
                               "Proposed Optimal" = "#2ca02c", 
                               "Alternative" = "#1f77b4"),
                    name = "Analysis Type") +
  scale_color_manual(values = c("Statistical Changepoint" = "#e31a1c", 
                                "Proposed Optimal" = "#2ca02c", 
                                "Alternative" = "#1f77b4"),
                     guide = "none") +
  
  labs(
    title = "Dual-Axis Validation: Statistical vs Practical Optimization",
    subtitle = "Comparison of metrics around the statistical changepoint",
    x = "SAG Threshold",
    y = "Metric Value"
  ) +
  theme_publication +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "03_dual_axis_validation.png"), p_multi_comparison, 
       width = 14, height = 10, dpi = 300)

cat("INFO: Creating publication-ready summary figure...\n")

p_summary_top <- snp1_data %>%
  ggplot(aes(x = sag_threshold, y = cluster_count)) +
  geom_line(color = "steelblue", size = 1.5, alpha = 0.8) +
  geom_point(size = 3, color = "steelblue", alpha = 0.8) +
  
  geom_point(data = snp1_data %>% filter(sag_threshold == 13), 
             color = "red", size = 5) +
  geom_point(data = snp1_data %>% filter(sag_threshold == 14), 
             color = "darkgreen", size = 5, shape = 17) +
  
  annotate("text", x = 13, y = 10.5, label = "Statistical\nChangepoint", 
           color = "red", size = 3, hjust = 0.5) +
  annotate("text", x = 14, y = 12.5, label = "Proposed\nOptimal", 
           color = "darkgreen", size = 3, hjust = 0.5) +
  
  geom_rect(xmin = 13.5, xmax = 16.5, ymin = 11.5, ymax = 12.5,
            fill = "lightgreen", alpha = 0.3) +
  annotate("text", x = 15, y = 12, label = "Stable Region", 
           color = "darkgreen", size = 3, fontface = "italic") +
  
  labs(title = "A. Cluster Stability Analysis",
       x = "SAG Threshold", y = "Cluster Count") +
  theme_publication +
  theme(plot.title = element_text(size = 14))

comparison_data <- validation_data %>%
  filter(sag_threshold %in% c(13, 14)) %>%
  select(sag_threshold, analysis_type, sag_count, cluster_count, n_content_mean) %>%
  mutate(
    improvement_sag = ifelse(sag_threshold == 14, 
                             paste("+", round((sag_count - validation_data$sag_count[validation_data$sag_threshold == 13]) / validation_data$sag_count[validation_data$sag_threshold == 13] * 100, 1), "%", sep=""), 
                             "Baseline"),
    improvement_cluster = ifelse(sag_threshold == 14,
                                 paste("+", round((cluster_count - validation_data$cluster_count[validation_data$sag_threshold == 13]) / validation_data$cluster_count[validation_data$sag_threshold == 13] * 100, 1), "%", sep=""),
                                 "Baseline")
  )

p_summary_bottom <- comparison_data %>%
  select(sag_threshold, analysis_type, sag_count, cluster_count) %>%
  pivot_longer(cols = c(sag_count, cluster_count), names_to = "metric", values_to = "value") %>%
  mutate(metric_label = ifelse(metric == "sag_count", "SAG Count", "Cluster Count")) %>%
  ggplot(aes(x = metric_label, y = value, fill = analysis_type)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = value), position = position_dodge(width = 0.8), 
            vjust = -0.3, size = 4, fontface = "bold") +
  
  scale_fill_manual(values = c("Statistical Changepoint" = "#e31a1c", 
                               "Proposed Optimal" = "#2ca02c"),
                    name = "Approach") +
  
  labs(title = "B. Quantitative Comparison",
       x = "Metric", y = "Count") +
  theme_publication +
  theme(plot.title = element_text(size = 14), legend.position = "bottom")

summary_combined <- grid.arrange(p_summary_top, p_summary_bottom, ncol = 1, heights = c(2, 1))

ggsave(file.path(output_dir, "04_publication_summary.png"), summary_combined, 
       width = 12, height = 10, dpi = 300)

cat("INFO: Creating SNP threshold stability analysis...\n")

snp_thresholds <- unique(plot_data$snp_threshold)
snp_stability_data <- plot_data %>%
  arrange(snp_threshold, sag_threshold) %>%
  group_by(snp_threshold) %>%
  mutate(
    cluster_change = cluster_count - lag(cluster_count, default = cluster_count[1]),
    is_stable = cluster_change == 0,
    stability_zone = cumsum(cluster_change != 0 | row_number() == 1)
  ) %>%
  ungroup()

stability_zones <- snp_stability_data %>%
  group_by(snp_threshold, stability_zone, cluster_count) %>%
  summarise(
    zone_start = min(sag_threshold),
    zone_end = max(sag_threshold),
    zone_length = n(),
    .groups = 'drop'
  ) %>%
  filter(zone_length >= 2) %>%
  arrange(snp_threshold, zone_start)

p_snp_stability <- snp_stability_data %>%
  ggplot(aes(x = sag_threshold, y = cluster_count)) +
  geom_line(color = "steelblue", size = 1, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.9, aes(color = factor(cluster_count))) +
  
  geom_rect(data = stability_zones, 
            aes(xmin = zone_start - 0.4, xmax = zone_end + 0.4, 
                ymin = cluster_count - 0.3, ymax = cluster_count + 0.3),
            fill = "lightgreen", alpha = 0.4, inherit.aes = FALSE) +
  
  geom_vline(data = plot_data %>% 
               group_by(snp_threshold) %>% 
               summarise(sag_changepoint = first(calculated_sag_changepoint), .groups = 'drop'),
             aes(xintercept = sag_changepoint), 
             linetype = "dashed", color = "red", alpha = 0.7) +
  
  facet_wrap(~paste("SNP Threshold =", snp_threshold), scales = "free", ncol = 2) +
  
  scale_color_viridis_d(name = "Cluster\nCount") +
  labs(
    x = "SAG Threshold",
    y = "Number of Phylogenetic Clusters"
  ) +
  theme_publication +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "05_snp_threshold_stability_analysis.png"), p_snp_stability, 
       width = 14, height = 10, dpi = 300)

cat("INFO: Creating stability vs resolution trade-off analysis...\n")

snp_metrics <- stability_zones %>%
  group_by(snp_threshold) %>%
  summarise(
    max_stable_cluster = max(cluster_count),
    longest_stable_zone = max(zone_length),
    num_stable_zones = n(),
    total_stability_score = sum(zone_length),
    stability_value_score = max(cluster_count * zone_length),
    .groups = 'drop'
  ) %>%
  mutate(
    is_snp_changepoint = (snp_threshold == snp_changepoint),
    high_resolution_stable = ifelse(max_stable_cluster >= 10, "Yes", "No")
  )

p_stability_resolution <- snp_metrics %>%
  ggplot(aes(x = max_stable_cluster, y = longest_stable_zone)) +
  geom_point(aes(color = factor(snp_threshold), size = stability_value_score), 
             alpha = 0.8) +
  
  geom_point(data = snp_metrics %>% filter(snp_threshold == snp_changepoint),
             color = "red", size = 8, shape = 1, stroke = 2) +
  
  geom_text(aes(label = paste("SNP =", snp_threshold)), 
            vjust = -1.5, hjust = 0.5, size = 4, fontface = "bold") +
  
  annotate("rect", xmin = 10, xmax = Inf, ymin = 3, ymax = Inf,
           fill = "lightblue", alpha = 0.3) +
  annotate("text", x = 11, y = 3.5, 
           label = "High Resolution\n& Stability Zone", 
           color = "darkblue", size = 3.5, fontface = "italic") +
  
  scale_color_viridis_d(name = "SNP\nThreshold") +
  scale_size_continuous(name = "Stability Value\nScore", range = c(4, 12)) +
  labs(
    title = "Stability vs Resolution Trade-off Analysis",
    subtitle = "Only SNP=1 achieves high resolution (≥10 clusters) with adequate stability",
    x = "Maximum Stable Cluster Count",
    y = "Longest Stable Zone (consecutive SAG thresholds)"
  ) +
  theme_publication

ggsave(file.path(output_dir, "06_stability_resolution_tradeoff.png"), p_stability_resolution, 
       width = 12, height = 8, dpi = 300)

cat("INFO: Creating phylogenetic value score comparison...\n")

p_value_comparison <- snp_metrics %>%
  mutate(
    snp_label = paste("SNP =", snp_threshold),
    snp_label = factor(snp_label, levels = paste("SNP =", sort(snp_threshold)))
  ) %>%
  ggplot(aes(x = snp_label, y = stability_value_score, fill = is_snp_changepoint)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = stability_value_score), 
            vjust = -0.5, size = 5, fontface = "bold") +
  
  geom_text(aes(label = paste("Max:", max_stable_cluster, "clusters\nStable:", longest_stable_zone, "zones")),
            vjust = 1.5, size = 3.5, color = "white", fontface = "bold") +
  
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkred"),
                    labels = c("FALSE" = "Alternative", "TRUE" = "Statistical Changepoint"),
                    name = "SNP Threshold Type") +
  
  labs(
    title = "Phylogenetic Value Score Comparison",
    subtitle = "Value Score = Maximum Stable Cluster Count × Longest Stable Zone",
    x = "SNP Threshold",
    y = "Phylogenetic Value Score"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "07_phylogenetic_value_comparison.png"), p_value_comparison, 
       width = 10, height = 8, dpi = 300)

cat("INFO: Creating comprehensive SNP threshold analysis summary...\n")

p_summary_stability <- snp_metrics %>%
  select(snp_threshold, max_stable_cluster, longest_stable_zone) %>%
  pivot_longer(cols = c(max_stable_cluster, longest_stable_zone), 
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = ifelse(metric == "max_stable_cluster", "Max Stable\nClusters", "Longest Stable\nZone"),
    is_changepoint = (snp_threshold == snp_changepoint)
  ) %>%
  ggplot(aes(x = factor(snp_threshold), y = value, fill = is_changepoint)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = value), vjust = -0.3, size = 4, fontface = "bold") +
  
  facet_wrap(~metric_label, scales = "free_y") +
  
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkred"),
                    name = "Changepoint") +
  
  labs(title = "A. SNP Threshold Stability Metrics",
       x = "SNP Threshold", y = "Count") +
  theme_publication +
  theme(plot.title = element_text(size = 14), legend.position = "none")

exclusivity_data <- snp_metrics %>%
  mutate(
    category = case_when(
      max_stable_cluster >= 10 & longest_stable_zone >= 3 ~ "High Resolution\n& High Stability",
      max_stable_cluster >= 10 ~ "High Resolution\nOnly",
      longest_stable_zone >= 3 ~ "High Stability\nOnly",
      TRUE ~ "Low Resolution\n& Low Stability"
    ),
    snp_label = paste("SNP", snp_threshold),
    count_value = 1
  )

p_summary_exclusivity <- exclusivity_data %>%
  ggplot(aes(x = snp_label, y = count_value)) +
  geom_col(aes(fill = category), alpha = 0.8, width = 0.6) +
  geom_text(aes(label = snp_label), y = 1.1, size = 4, fontface = "bold") +
  
  scale_fill_manual(values = c(
    "High Resolution\n& High Stability" = "#2ca02c",
    "High Resolution\nOnly" = "#ff7f0e", 
    "High Stability\nOnly" = "#1f77b4",
    "Low Resolution\n& Low Stability" = "#d62728"
  ), name = "Performance Category") +
  
  labs(title = "B. High-Resolution Stability Exclusivity",
       subtitle = "Only SNP=1 achieves both high resolution (≥10 clusters) and stability (≥3 zones)",
       x = "SNP Threshold", y = "") +
  theme_publication +
  theme(plot.title = element_text(size = 14), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

snp_summary_combined <- grid.arrange(p_summary_stability, p_summary_exclusivity, 
                                     ncol = 1, heights = c(1.2, 1))

ggsave(file.path(output_dir, "08_snp_threshold_analysis_summary.png"), snp_summary_combined, 
       width = 14, height = 10, dpi = 300)

sag_summary_stats <- validation_data %>%
  filter(sag_threshold %in% c(13, 14)) %>%
  select(sag_threshold, analysis_type, sag_count, final_snp_count, n_content_mean, cluster_count) %>%
  mutate(
    sag_improvement = round((sag_count - sag_count[sag_threshold == 13]) / sag_count[sag_threshold == 13] * 100, 1),
    snp_improvement = round((final_snp_count - final_snp_count[sag_threshold == 13]) / final_snp_count[sag_threshold == 13] * 100, 1),
    cluster_improvement = round((cluster_count - cluster_count[sag_threshold == 13]) / cluster_count[sag_threshold == 13] * 100, 1),
    n_content_change = round((n_content_mean - n_content_mean[sag_threshold == 13]) / n_content_mean[sag_threshold == 13] * 100, 1)
  )

snp_summary_stats <- snp_metrics %>%
  mutate(
    resolution_efficiency = max_stable_cluster / snp_threshold,
    stability_efficiency = longest_stable_zone / snp_threshold,
    overall_efficiency = stability_value_score / snp_threshold
  )

write_csv(sag_summary_stats, file.path(output_dir, "sag_threshold_validation_summary.csv"))
write_csv(snp_summary_stats, file.path(output_dir, "snp_threshold_stability_summary.csv"))
write_csv(stability_zones, file.path(output_dir, "stability_zones_detailed.csv"))

cat("INFO: Enhanced visualization with SNP threshold analysis completed successfully.\n")
cat("Generated files:\n")
cat("  - 01_cluster_stability_analysis.png: SAG threshold stability (main finding)\n")
cat("  - 02_data_volume_vs_clusters.png: Data volume trade-off analysis\n") 
cat("  - 03_dual_axis_validation.png: Comprehensive multi-metric comparison\n")
cat("  - 04_publication_summary.png: Publication-ready SAG analysis summary\n")
cat("  - 05_snp_threshold_stability_analysis.png: SNP threshold stability comparison\n")
cat("  - 06_stability_resolution_tradeoff.png: Stability vs resolution trade-off\n")
cat("  - 07_phylogenetic_value_comparison.png: Phylogenetic value score comparison\n")
cat("  - 08_snp_threshold_analysis_summary.png: Comprehensive SNP threshold summary\n")
cat("  - sag_threshold_validation_summary.csv: SAG threshold quantitative summary\n")
cat("  - snp_threshold_stability_summary.csv: SNP threshold stability metrics\n")
cat("  - stability_zones_detailed.csv: Detailed stability zone information\n")
cat(paste("Output directory:", output_dir, "\n"))