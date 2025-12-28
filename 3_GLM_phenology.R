# ============================================================
# 3_GLM_phenology.R
# GLM for phenology traits using climate PCs + Batch
# Masking-safe (explicit dplyr:: calls everywhere)
# ============================================================

# ---- packages ----
library(dplyr)
library(tidyr)
library(DHARMa)
library(MASS)

# ============================================================
# 1) Load data
# ============================================================

dat <- read.csv("./Phenology project_Climate.csv")

req_cols <- c("Population","Batch","State","Latitude","Longitude",
              "MAT","MAP","PrecSeason","GrowPrecip","SolarRad","Elevation_m","CGDD")
stopifnot(all(req_cols %in% names(dat)))

# standardize naming
dat <- dat %>%
  dplyr::rename(Elevation = Elevation_m)

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

pca <- prcomp(clim_mat, center = TRUE, scale. = TRUE)

pcs <- as.data.frame(pca$x[, 1:3, drop = FALSE])
colnames(pcs) <- c("PC1","PC2","PC3")
pcs$Population <- site_df$Population

# join PCs back to individual plants
dat_pca <- dat %>%
  dplyr::left_join(pcs, by = "Population") %>%
  dplyr::mutate(
    Batch = factor(Batch),
    Population = factor(Population)
  )

# ============================================================
# 3) Choose response variable
# ============================================================

response <- "Days_to_bud"   # <-- CHANGE AS NEEDED
stopifnot(response %in% names(dat_pca))

df <- dat_pca %>%
  dplyr::select(dplyr::all_of(c(response, "PC1", "PC2", "Batch"))) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    y = .[[response]],
    y_log10 = log10(y)
  )

# ============================================================
# 4) Model formulas
# ============================================================

f_gaus <- stats::as.formula(
  "y_log10 ~ PC1 + PC2 + PC1:PC2 + Batch + Batch:PC1 + Batch:PC2"
)

f_cnt <- stats::as.formula(
  "y ~ PC1 + PC2 + PC1:PC2 + Batch + Batch:PC1 + Batch:PC2"
)

# ============================================================
# 5) Fit models
# ============================================================

# ---- Gaussian GLM (log10 timing) ----
m_gaus <- stats::glm(f_gaus, data = df, family = gaussian())
cat("\n==== Gaussian GLM (log10) ====\n")
print(summary(m_gaus))
print(stats::anova(m_gaus, test = "F"))

sim_gaus <- DHARMa::simulateResiduals(m_gaus, n = 1000)
plot(sim_gaus)
DHARMa::testUniformity(sim_gaus)

# ---- Poisson GLM ----
m_pois <- stats::glm(f_cnt, data = df, family = poisson())
cat("\n==== Poisson GLM ====\n")
print(summary(m_pois))

dispersion <- sum(stats::residuals(m_pois, type = "pearson")^2) /
              stats::df.residual(m_pois)
cat("\nPoisson dispersion: ", round(dispersion, 3), "\n", sep = "")

sim_pois <- DHARMa::simulateResiduals(m_pois, n = 1000)
plot(sim_pois)
DHARMa::testDispersion(sim_pois)

# ---- Negative binomial GLM ----
m_nb <- MASS::glm.nb(f_cnt, data = df)
cat("\n==== Negative Binomial GLM ====\n")
print(summary(m_nb))

sim_nb <- DHARMa::simulateResiduals(m_nb, n = 1000)
plot(sim_nb)
DHARMa::testDispersion(sim_nb)
