suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("5つの引数が必要です: input_file, output_threshold_file, output_plot_file, x_label, y_label", call. = FALSE)
}

input_file <- args[1]
output_threshold_file <- args[2]
output_plot_file <- args[3]
x_label <- args[4]
y_label <- args[5]

dir.create(dirname(output_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_threshold_file), showWarnings = FALSE, recursive = TRUE)

data <- tryCatch({
  read_tsv(input_file, show_col_types = FALSE)
}, error = function(e) {
  stop(paste("入力ファイルの読み込みに失敗しました:", input_file, "\nOriginal error:", e$message), call. = FALSE)
})

if (nrow(data) == 0) {
  stop(paste("入力ファイルが空です:", input_file, ". 変化点解析を実行できません。"), call. = FALSE)
}

colnames(data) <- c("Value", "Count")

cpt <- cpt.mean(data$Count, method = "AMOC", penalty = "MBIC")
changepoint_index <- cpts(cpt)

if (length(changepoint_index) == 0 || changepoint_index == 0) {
  cat("WARNING: No changepoint detected. Defaulting threshold to 2 to remove singletons.\n")
  threshold <- 2
} else {
  threshold <- data$Value[changepoint_index]
}

writeLines(as.character(threshold), output_threshold_file)
cat("INFO: Detected changepoint threshold:", threshold, "\n")

plot <- ggplot(data, aes(x = Value, y = Count)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1.2) +
  labs(x = x_label, y = y_label, title = "Changepoint Analysis") +
  theme_bw(base_size = 16) +
  annotate("text", x = threshold, y = max(data$Count) * 0.9,
           label = paste("Threshold =", threshold), color = "red", hjust = -0.1, size = 5)

ggsave(output_plot_file, plot, width = 8, height = 6, dpi = 300)
cat("INFO: Plot saved to:", output_plot_file, "\n")