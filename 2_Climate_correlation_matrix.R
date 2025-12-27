# ============================================================
# 2_Climate_correlation_matrix.R
# Correlation matrix of population-level climate variables
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggcorrplot)

# ============================================================
# 1) Load data (standalone-safe)
# ============================================================

data <- read.csv("./Phenology project_Climate.csv")

req_cols <- c("Population","State","Latitude","Longitude",
              "MAT","MAP","PrecSeason","GrowPrecip","SolarRad","Elevation_m","CGDD")
stopifnot(all(req_cols %in% names(data)))

# standardize elevation name
data <- data %>%
  dplyr::rename(Elevation = Elevation_m)

# ============================================================
# 2) Build population-level site table (if needed)
# ============================================================

clim_vars <- c(
  "Latitude",
  "CGDD",
  "MAP",
  "Elevation",
  "MAT",
  "PrecSeason",
  "GrowPrecip",
  "SolarRad"
)

site_df <- data %>%
  dplyr::select(Population, State, all_of(clim_vars)) %>%
  group_by(Population, State) %>%
  summarise(
    across(all_of(clim_vars), ~ mean(.x, na.rm = TRUE)),
    N = dplyr::n(),
    .groups = "drop"
  ) %>%
  drop_na(all_of(clim_vars))

stopifnot(n_distinct(site_df$Population) == nrow(site_df))

# ============================================================
# 3) Correlation matrix (Pearson)
# ============================================================

clim_df <- site_df %>% dplyr::select(all_of(clim_vars))

cor_mat <- cor(
  clim_df,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# optional p-value matrix (for reporting / masking)
cor_p <- matrix(NA, ncol = ncol(clim_df), nrow = ncol(clim_df))
colnames(cor_p) <- rownames(cor_p) <- colnames(clim_df)

for (i in seq_len(ncol(clim_df))) {
  for (j in seq_len(ncol(clim_df))) {
    cor_p[i, j] <- cor.test(
      clim_df[[i]],
      clim_df[[j]],
      method = "pearson"
    )$p.value
  }
}

# ============================================================
# 4) Save correlation table (optional)
# ============================================================

cor_tbl <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("Variable")

# write.csv(cor_tbl, "tables/Climate_correlation_matrix.csv", row.names = FALSE)

# ============================================================
# 5) Plot correlation matrix
# ============================================================

p_cor <- ggcorrplot(
  cor_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  lab_size = 3.2,
  colors = c("#2166AC", "white", "#B2182B"),
  outline.color = "grey60",
  ggtheme = theme_pub(base_size = 11)
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11)
  )

print(p_cor)

# ggsave("figs/Fig_climate_correlation_matrix.png",
#        p_cor, width = 6.2, height = 5.6, dpi = 400)
