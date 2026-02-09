suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

theme_nature <- theme_classic(base_size = 8) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(1.2, "mm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    plot.subtitle = element_text(size = 8, hjust = 0),
    plot.margin = margin(6, 6, 6, 6, "pt")
  )

theme_umap <- theme_nature +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

theme_set(theme_nature)

save_figure <- function(plot, folder, filename, width = 3.5, height = 3.5) {
  ggsave(
    filename = file.path(folder, paste0(filename, ".svg")),
    plot = plot,
    width = width, height = height,
    device = svglite::svglite
  )
  ggsave(
    filename = file.path(folder, paste0(filename, ".png")),
    plot = plot,
    width = width, height = height,
    dpi = 600
  )
}

clean_umap <- function(p, title = NULL, subtitle = NULL) {
  p +
    coord_fixed() +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    theme_umap
}
