################################################################################
#
#  Running Health Models in R
#
#  Script 02 - Second-stage meta-analysis and BLUP by IGR
#
################################################################################

#------------------------------------------------------------------------------
# Purpose
#------------------------------------------------------------------------------
#
# This script performs the second-stage analysis using the estimates obtained
# from the first-stage IGR-specific time-series models.
#
# The aims are:
#
#   1. Fit a multivariate meta-analysis model with random effects only.
#   2. Obtain the pooled temperature-mortality association.
#   3. Plot the pooled exposure-response association for Rio de Janeiro.
#   4. Predict IGR-specific associations using Best Linear Unbiased Predictions
#      (BLUPs).
#   5. Estimate the Minimum Mortality Temperature (MMT) for each IGR, restricting
#      the search to the 25th to 98th percentiles of the local temperature
#      distribution.
#
# This script assumes that the objects created in Script 01 are already
# available in the R environment:
#
#   - data_final
#   - coefs_all_mat
#   - vcov_all
#   - igr_include
#
################################################################################


#==============================================================================
# 1. Output folders
#==============================================================================

dir_save <- "results"
dir_fig  <- file.path(dir_save, "figures")

if (!dir.exists(dir_save)) {
  dir.create(dir_save)
}

if (!dir.exists(dir_fig)) {
  dir.create(dir_fig)
}


#==============================================================================
# 2. Define temperature variable
#==============================================================================

# This must be the same temperature variable used in the first-stage models.

temp_var <- "daily_pop_weighted_mean_temperature"


#==============================================================================
# 3. Define the temperature basis for prediction
#==============================================================================

# The second-stage model combines the reduced cumulative associations obtained
# in the first stage.
#
# To predict pooled and IGR-specific curves, we recreate a one-dimensional
# temperature basis consistent with the first-stage model specification.

predvar <- quantile(
  data_final[[temp_var]],
  probs = (0:100) / 100,
  na.rm = TRUE
)

bvar <- onebasis(
  predvar,
  fun = "ns",
  knots = quantile(predvar, c(0.10, 0.90), na.rm = TRUE)
)


#==============================================================================
# 4. Fit the second-stage meta-analysis model
#==============================================================================

# The meta-analysis combines the IGR-specific reduced coefficients from the
# first stage, accounting for both within-IGR uncertainty and between-IGR
# heterogeneity.
#
# In this course example, the model includes random effects only, without
# meta-predictors.

meta_model <- mixmeta(
  coefs_all_mat,
  vcov_all,
  method = "reml",
  control = list(showiter = TRUE)
)

summary(meta_model)


#==============================================================================
# 5. Obtain the pooled exposure-response curve
#==============================================================================

# First, we predict the pooled curve without centering to identify the pooled
# Minimum Mortality Temperature.

cp_meta_uncentered <- crosspred(
  bvar,
  coef = coef(meta_model),
  vcov = vcov(meta_model),
  model.link = "log",
  by = 0.1
)

MMT_pooled <- cp_meta_uncentered$predvar[
  which.min(cp_meta_uncentered$allRRfit)
]


# Then, we centre the pooled curve at the pooled MMT.
# After centering, RR = 1 at the pooled MMT.

cp_meta <- crosspred(
  bvar,
  coef = coef(meta_model),
  vcov = vcov(meta_model),
  model.link = "log",
  by = 0.1,
  cen = MMT_pooled
)


#==============================================================================
# 6. Organise pooled exposure-response curve
#==============================================================================

pooled_curve <- data.frame(
  Temperature = cp_meta$predvar,
  RR          = cp_meta$allRRfit,
  RR_low      = cp_meta$allRRlow,
  RR_high     = cp_meta$allRRhigh,
  MMT_pooled  = MMT_pooled
)


#==============================================================================
# 7. Plot pooled exposure-response association
#==============================================================================

# The pooled exposure-response curve summarizes the average temperature-
# mortality association across all IGRs in Rio de Janeiro after accounting for
# between-IGR heterogeneity.

P99 <- quantile(
  data_final[[temp_var]],
  probs = 0.99,
  na.rm = TRUE
)

x_breaks <- sort(unique(round(
  c(
    min(pooled_curve$Temperature),
    MMT_pooled,
    P99,
    max(pooled_curve$Temperature)
  ),
  1
)))

p_pooled <- ggplot(
  pooled_curve,
  aes(x = Temperature, y = RR)
) +
  geom_ribbon(
    aes(ymin = RR_low, ymax = RR_high),
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(
    xintercept = MMT_pooled,
    linetype = 2
  ) +
  geom_vline(
    xintercept = P99,
    linetype = 3
  ) +
  scale_x_continuous(breaks = x_breaks) +
  labs(
    title = "Pooled temperature-mortality association - Rio de Janeiro",
    x = "Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_classic()

p_pooled

ggsave(
  filename = file.path(dir_fig, "pooled_exposure_response_curve.png"),
  plot = p_pooled,
  width = 7,
  height = 5,
  dpi = 300
)


#==============================================================================
# 8. Obtain IGR-specific BLUPs
#==============================================================================

# BLUPs provide stabilized IGR-specific estimates by borrowing information
# across all IGRs included in the meta-analysis.

blup_all <- blup(
  meta_model,
  vcov = TRUE
)

names(blup_all) <- as.character(names(blup_all))


#==============================================================================
# 9. Estimate MMT by IGR using BLUPs
#==============================================================================

# The MMT is estimated separately for each IGR using the BLUP coefficients.
#
# To avoid selecting unrealistic values at the tails of the local temperature
# distribution, the MMT search is restricted to temperatures between the 25th
# and 98th percentiles.

mmt_igr <- data.frame(
  IGR      = igr_include,
  MMT      = NA_real_,
  MMT_perc = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(igr_include)) {
  
  igr_code <- as.character(igr_include[i])
  
  cat("Estimating MMT for IGR", i, "of", length(igr_include),
      "- IGR:", igr_code, "\n")
  
  data_igr <- data_final |>
    filter(IGR == igr_code)
  
  perc_seq <- 25:98
  
  temps <- quantile(
    data_igr[[temp_var]],
    probs = perc_seq / 100,
    na.rm = TRUE
  )
  
  bvar_i <- onebasis(
    temps,
    fun = "ns",
    knots = quantile(data_igr[[temp_var]], c(0.10, 0.90), na.rm = TRUE),
    Bound = range(data_igr[[temp_var]], na.rm = TRUE)
  )
  
  coef_i <- blup_all[[igr_code]]$blup
  
  rr_i <- as.numeric(
    exp(bvar_i %*% coef_i)
  )
  
  mmt_position <- which.min(rr_i)
  
  mmt_igr$MMT_perc[i] <- perc_seq[mmt_position]
  mmt_igr$MMT[i] <- as.numeric(temps[mmt_position])
}


#==============================================================================
# 10. Save second-stage results
#==============================================================================

write.csv(
  pooled_curve,
  file = file.path(dir_save, "pooled_exposure_response_curve.csv"),
  row.names = FALSE
)

write.csv(
  mmt_igr,
  file = file.path(dir_save, "MMT_by_IGR_after_BLUP.csv"),
  row.names = FALSE
)

saveRDS(
  meta_model,
  file = file.path(dir_save, "second_stage_meta_model.rds")
)

saveRDS(
  cp_meta,
  file = file.path(dir_save, "pooled_crosspred.rds")
)

saveRDS(
  blup_all,
  file = file.path(dir_save, "BLUP_by_IGR.rds")
)

################################################################################
# End of script
################################################################################
