# ============================================================
# 1_PCA_climate.R
# PCA of population-level climate/geography predictors
# Outputs: pca objects + PCs joined to individual data + biplot
# ============================================================

# ---- packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)

# (modeling packages can be loaded later in the modeling script;
# keeping this file focused on PCA/figure is usually cleaner)

# ---- theme ----
theme_pub <- function(base_size = 12, base_family = "sans"){
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey85"),
      axis.text        = ggplot2::element_text(colour = "black"),
      axis.title       = ggplot2::element_text(colour = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.25, colour = "grey70"),
      panel.border     = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.4),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0, margin = ggplot2::margin(b = 6)),
      legend.title     = ggplot2::element_text(face = "bold")
    )
}
th_pub <- theme_pub(base_size = 11)

# ---- read data ----
# setwd("./")  # optional; prefer RStudio Project paths if possible
data <- read.csv("./Phenology project_Climate.csv")

# ---- required columns ----
req_cols <- c("Population","State","Latitude","Longitude",
              "MAT","MAP","PrecSeason","GrowPrecip","SolarRad","Elevation_m","CGDD")
stopifnot(all(req_cols %in% names(data)))

# ---- rename once (keep consistent naming downstream) ----
df_src <- data %>%
  dplyr::rename(Elevation = Elevation_m)

# ============================================================
# 1) Build population-level site table (one row per population)
# ============================================================

clim_vars <- c("Latitude","CGDD","MAP","Elevation","MAT","PrecSeason","GrowPrecip","SolarRad")

site_df <- df_src %>%
  dplyr::select(Population, State, Longitude, dplyr::all_of(clim_vars)) %>%
  group_by(Population, State) %>%
  summarise(
    N = dplyr::n(),
    Longitude = mean(Longitude, na.rm = TRUE),
    across(all_of(clim_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::drop_na(all_of(clim_vars))

stopifnot(dplyr::n_distinct(site_df$Population) == nrow(site_df))

# ============================================================
# 2) PCA (population-level)
# ============================================================

clim_mat <- site_df %>%
  dplyr::select(all_of(clim_vars)) %>%
  scale(center = TRUE, scale = TRUE)

pca <- prcomp(clim_mat, center = TRUE, scale. = TRUE)

# variance explained table
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_tbl <- tibble::tibble(
  PC = paste0("PC", seq_along(pca_var)),
  var_explained = pca_var,
  cum_var = cumsum(pca_var)
)

# loadings table
loadings_tbl <- as.data.frame(pca$rotation) %>%
  tibble::rownames_to_column("variable")

print(pca_tbl)
print(loadings_tbl)

# ============================================================
# 3) Extract PCs and join back to individual plants
# ============================================================

pcs <- as.data.frame(pca$x[, 1:3, drop = FALSE])
colnames(pcs) <- c("PC1","PC2","PC3")
pcs$Population <- site_df$Population

data_pca <- df_src %>%
  left_join(pcs %>% dplyr::select(Population, PC1, PC2, PC3), by = "Population")

# ============================================================
# 4) PCA biplot figure (PC1 vs PC2)
# ============================================================

scores <- site_df %>%
  dplyr::select(Population, State, N) %>%
  left_join(pcs %>% dplyr::select(Population, PC1, PC2), by = "Population")

ve <- (pca$sdev^2) / sum(pca$sdev^2)
x_lab <- sprintf("PC1 (%.1f%%)", 100 * ve[1])
y_lab <- sprintf("PC2 (%.1f%%)", 100 * ve[2])

loads <- as.data.frame(pca$rotation[, 1:2, drop = FALSE])
loads$var <- rownames(loads)
colnames(loads)[1:2] <- c("PC1","PC2")

# scale arrows to match score space (simple and stable)
arrow_scale <- 0.85 * min(
  diff(range(scores$PC1, na.rm = TRUE)) / diff(range(loads$PC1, na.rm = TRUE)),
  diff(range(scores$PC2, na.rm = TRUE)) / diff(range(loads$PC2, na.rm = TRUE))
)

loads_plot <- dplyr::transmute(loads,
  var,
  xend = PC1 * arrow_scale,
  yend = PC2 * arrow_scale
)

# fallback palette if state_pal not defined
if (!exists("state_pal")) {
  st <- sort(unique(scores$State))
  state_pal <- setNames(scales::hue_pal()(length(st)), st)
}

p_biplot <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_point(
    data = scores,
    ggplot2::aes(PC1, PC2, color = State),
    size = 2.2, alpha = 0.9
  ) +
  ggplot2::geom_segment(
    data = loads_plot,
    ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = ggplot2::arrow(length = grid::unit(0.02, "npc")),
    linewidth = 0.45,
    color = "grey25"
  ) +
  ggrepel::geom_text_repel(
    data = loads_plot,
    ggplot2::aes(x = xend, y = yend, label = var),
    size = 3.1, color = "grey20", min.segment.length = 0
  ) +
  ggplot2::scale_color_manual(values = state_pal, name = "State") +
  ggplot2::labs(x = x_lab, y = y_lab) +
  ggplot2::coord_equal() +
  th_pub +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

print(p_biplot)

# optional: save
# ggsave("figs/Fig_climate_PCA_biplot.png", p_biplot, width = 6.2, height = 5.2, dpi = 400)
