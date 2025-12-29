# ============================================================
# 4_Growth_Repro_Gaussian_GLM.R
# Gaussian GLMs for growth & reproduction traits
# ============================================================

# ---- packages ----
library(dplyr)
library(tidyr)
library(readr)
library(broom)

# ============================================================
# 1) Load data
# ============================================================

infile  <- "./Phenology project_Climate.csv"
out_dir <- "tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(infile)

req_cols <- c("Population","Batch","State","Latitude","Longitude",
              "MAT","MAP","PrecSeason","GrowPrecip","SolarRad","Elevation_m","CGDD")
stopifnot(all(req_cols %in% names(dat)))

# standardize elevation name
dat <- dat %>% dplyr::rename(Elevation = Elevation_m)

# ============================================================
# 2) Build population-level climate PCs
# ============================================================

clim_vars <- c(
  "Latitude","CGDD","MAP","Elevation",
  "MAT","PrecSeason","GrowPrecip","SolarRad"
)

site_df <- dat %>%
  dplyr::select(Population, dplyr::all_of(clim_vars)) %>%
  dplyr::group_by(Population) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(clim_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::drop_na(dplyr::all_of(clim_vars))

stopifnot(dplyr::n_distinct(site_df$Population) == nrow(site_df))

clim_mat <- site_df %>%
  dplyr::select(dplyr::all_of(clim_vars)) %>%
  scale(center = TRUE, scale = TRUE)

pca <- stats::prcomp(clim_mat, center = TRUE, scale. = TRUE)

pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
colnames(pcs) <- c("PC1","PC2")
pcs$Population <- site_df$Population

dat_pca <- dat %>%
  dplyr::left_join(pcs, by = "Population") %>%
  dplyr::mutate(
    Batch = factor(Batch),
    Population = factor(Population)
  )

# ============================================================
# 3) Traits and transformations
#    EDIT names to match your CSV exactly
# ============================================================

traits <- tibble::tibble(
  trait = c(
    "Primary_inflorescence",
    "Total_inflorescence_length",
    "Stem_diameter",
    "First_bud_height",
    "Total_height"
  ),
  transform = c("log10","log10","log10","raw","raw")
)

stopifnot(all(traits$trait %in% names(dat_pca)))

safe_log10 <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) log10(x + 1) else log10(x)
}

# ============================================================
# 4) Fit Gaussian GLMs
# ============================================================

fit_one <- function(trait, transform_flag) {

  df <- dat_pca %>%
    dplyr::select(Batch, PC1, PC2, y = all_of(trait)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      y_model = if (transform_flag == "log10") safe_log10(y) else y
    )

  f <- stats::as.formula(
    "y_model ~ PC1 + PC2 + Batch + Batch:PC1 + Batch:PC2"
  )

  m <- stats::glm(f, data = df, family = gaussian())

  list(
    trait = trait,
    transform = transform_flag,
    n = nrow(df),
    model = m
  )
}

fits <- purrr::pmap(
  list(traits$trait, traits$transform),
  fit_one
)
