library(ggplot2)
library(dplyr)

# Veriyi oku
results <- read.csv("results_final.csv")

# Factor sıralamaları
results$Variance <- factor(results$Variance, levels = c("Homogeneous", "Moderate", "High"))
results$Test     <- factor(results$Test, levels = c("TukeyHSD","Scheffe","FisherLSD",
                                                    "GamesHowell","DunnettT3","TamhaneT2","Sidak"))
results$n_label  <- paste0("n = ", results$n)
results$n_label  <- factor(results$n_label, levels = c("n = 6","n = 8","n = 10"))
results$k_label  <- paste0("k = ", results$k)
results$k_label  <- factor(results$k_label, levels = c("k = 3","k = 4","k = 5"))

# Renk ve şekil paleti (7 test)
test_colors <- c(
  "TukeyHSD"   = "#1f77b4",
  "Scheffe"    = "#ff7f0e",
  "FisherLSD"  = "#2ca02c",
  "GamesHowell"= "#d62728",
  "DunnettT3"  = "#9467bd",
  "TamhaneT2"  = "#8c564b",
  "Sidak"      = "#e377c2"
)
test_shapes <- c(
  "TukeyHSD"   = 16,
  "Scheffe"    = 17,
  "FisherLSD"  = 15,
  "GamesHowell"= 18,
  "DunnettT3"  = 4,
  "TamhaneT2"  = 8,
  "Sidak"      = 3
)

# Bradley sınırları
bradley_low  <- 0.025
bradley_high <- 0.075

# ─────────────────────────────────────────
# FIG 1 — POWER
# ─────────────────────────────────────────
fig1 <- ggplot(results, aes(x = Test, y = Power,
                            color = Test, shape = Test, group = Test)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_grid(n_label ~ Variance + k_label) +
  scale_color_manual(values = test_colors) +
  scale_shape_manual(values = test_shapes) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x     = "Post-hoc test",
    y     = "Power",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    strip.text      = element_text(size = 7),
    legend.position = "bottom",
    legend.text     = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

tiff("Fig1.tiff", width = 174, height = 200,
     units = "mm", res = 600, compression = "lzw")
print(fig1)
dev.off()

# ─────────────────────────────────────────
# FIG 2 — TYPE I ERROR
# ─────────────────────────────────────────
fig2 <- ggplot(results, aes(x = Test, y = TypeI_Error,
                            color = Test, shape = Test, group = Test)) +
  geom_hline(yintercept = bradley_low,  linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_hline(yintercept = bradley_high, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_hline(yintercept = 0.05, linetype = "dotted",
             color = "black", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_grid(n_label ~ Variance + k_label) +
  scale_color_manual(values = test_colors) +
  scale_shape_manual(values = test_shapes) +
  scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)) +
  labs(
    x     = "Post-hoc test",
    y     = "Type I Error Rate",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    strip.text       = element_text(size = 7),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

tiff("Fig2.tiff", width = 174, height = 200,
     units = "mm", res = 600, compression = "lzw")
print(fig2)
dev.off()