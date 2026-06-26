################################################################################
#
#  BREATHE
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
#   3. Predict IGR-specific associations using Best Linear Unbiased Predictions
#      (BLUPs).
#   4. Estimate the Minimum Mortality Temperature (MMT) for each IGR, restricting
#      the search to the 25th to 98th percentiles of the local temperature
#      distribution.
#
# For teaching purposes, we fit a meta-analysis model without meta-predictors.
# In applied studies, geographic, climatic, demographic, or socioeconomic
# variables may be included as meta-predictors.
#
################################################################################


#==============================================================================
# 1. Output folder
#==============================================================================

dir_save <- "results"

if (!dir.exists(dir_save)) {
  dir.create(dir_save)
}


#==============================================================================
# 2. Load first-stage results and temperature data
#==============================================================================

# Load the coefficients and variance-covariance matrices estimated in the
# first-stage models.

load(file.path(dir_save, "first_stage_results.RData"))


# Load the original dataset.
# Here we only need the temperature series to define the exposure basis and
# to estimate the MMT within each IGR.

data_final <- read_parquet("banco_BREATHE_IGR_SIM.parquet")

data_final <- data_final |>
  mutate(Date = as.Date(Date)) |>
  filter(!is.na(IGR))


#==============================================================================
# 3. Identify valid first-stage estimates
#==============================================================================

# Before fitting the meta-analysis, we keep only IGRs with valid coefficients
# and valid variance-covariance matrices.
#
# An IGR is excluded if:
#
#   - any coefficient is missing;
#   - the covariance matrix is missing;
#   - the covariance matrix has the wrong dimension;
#   - the covariance matrix contains missing values;
#   - the covariance matrix is not positive definite.

is_invalid_vcov <- function(mat, expected_dim) {
  
  if (!is.matrix(mat)) return(TRUE)
  if (any(is.na(mat))) return(TRUE)
  if (any(dim(mat) != expected_dim)) return(TRUE)
  
  eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  
  any(eig <= 0)
}

expected_dim <- ncol(coefs_all_mat)

invalid_coef_igrs <- rownames(coefs_all_mat)[
  apply(coefs_all_mat, 1, anyNA)
]

invalid_vcov_igrs <- names(vcov_all)[
  sapply(vcov_all, is_invalid_vcov, expected_dim = expected_dim)
]

valid_igrs <- setdiff(
  igr_include,
  union(invalid_coef_igrs, invalid_vcov_igrs)
)


# Keep only valid IGRs

coefs_all_mat_clean <- coefs_all_mat[
  rownames(coefs_all_mat) %in% valid_igrs,
  ,
  drop = FALSE
]

vcov_all_clean <- vcov_all[
  rownames(coefs_all_mat_clean)
]


#==============================================================================
# 4. Define the temperature basis for prediction
#==============================================================================

# The first-stage coefficients describe the overall cumulative temperature-
# mortality association using a natural cubic spline with knots at the 10th and
# 90th percentiles.
#
# To obtain pooled and IGR-specific curves, we must recreate the same
# one-dimensional basis used after cross-reduction.

temp_var <- "daily_pop_weighted_mean_temperature"

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
# 5. Fit the second-stage meta-analysis model
#==============================================================================

# In this simplified course example, the second-stage model includes only
# random effects. This allows the temperature-mortality association to vary
# across IGRs, without introducing meta-predictors.
#
# The input consists of:
#
#   - coefs_all_mat_clean: reduced coefficients from the first-stage models;
#   - vcov_all_clean: corresponding within-IGR covariance matrices.
#
# The model is fitted using restricted maximum likelihood (REML).

meta_model <- mixmeta(
  coefs_all_mat_clean,
  vcov_all_clean,
  method = "reml",
  control = list(showiter = TRUE)
)

summary(meta_model)


#==============================================================================
# 6. Obtain the pooled exposure-response curve
#==============================================================================

# First, we predict the pooled curve without centering in order to identify
# the pooled Minimum Mortality Temperature (MMT).

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


# Then, we predict the pooled curve again, now centering the relative risks
# at the pooled MMT. After centering, RR = 1 at the MMT.

cp_meta <- crosspred(
  bvar,
  coef = coef(meta_model),
  vcov = vcov(meta_model),
  model.link = "log",
  by = 0.1,
  cen = MMT_pooled
)


#==============================================================================
# 7. Save pooled exposure-response curve
#==============================================================================

pooled_curve <- data.frame(
  Temperature = cp_meta$predvar,
  RR          = cp_meta$allRRfit,
  RR_low      = cp_meta$allRRlow,
  RR_high     = cp_meta$allRRhigh,
  MMT_pooled  = MMT_pooled
)

write.csv(
  pooled_curve,
  file = file.path(dir_save, "pooled_exposure_response_curve.csv"),
  row.names = FALSE
)


#==============================================================================
# 8. Obtain IGR-specific BLUPs
#==============================================================================

# BLUPs borrow strength from the full set of IGRs and provide stabilized
# IGR-specific estimates of the temperature-mortality association.
#
# They are useful especially when some IGRs have small populations, sparse
# mortality counts, or imprecise first-stage estimates.

blup_all <- blup(
  meta_model,
  vcov = TRUE
)

save(
  blup_all,
  file = file.path(dir_save, "blup_results.RData")
)


#==============================================================================
# 9. Estimate MMT by IGR using BLUPs
#==============================================================================

# The MMT is estimated separately for each IGR using the BLUP coefficients.
#
# To avoid selecting unrealistic values at the extreme cold or hot tails of the
# local temperature distribution, we restrict the MMT search to temperatures
# between the 25th and 98th percentiles of each IGR-specific distribution.

names(blup_all) <- as.character(names(blup_all))

mmt_igr <- data.frame(
  IGR      = valid_igrs,
  MMT      = NA_real_,
  MMT_perc = NA_real_,
  status   = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(valid_igrs)) {
  
  igr_code <- as.character(valid_igrs[i])
  
  cat("Estimating MMT for IGR", i, "of", length(valid_igrs),
      "- IGR:", igr_code, "\n")
  
  tryCatch({
    
    #----------------------------------------------------------------------
    # Step 1. Extract the temperature series for the current IGR
    #----------------------------------------------------------------------
    
    data_igr <- data_final |>
      filter(IGR == igr_code)
    
    #----------------------------------------------------------------------
    # Step 2. Define candidate temperatures for the MMT search
    #----------------------------------------------------------------------
    
    # Candidate temperatures are defined as local percentiles from P25 to P98.
    # This keeps the MMT search within a plausible range of observed
    # temperatures for each IGR.
    
    perc_seq <- 25:98
    
    temps <- quantile(
      data_igr[[temp_var]],
      probs = perc_seq / 100,
      na.rm = TRUE
    )
    
    #----------------------------------------------------------------------
    # Step 3. Recreate the temperature basis for this IGR
    #----------------------------------------------------------------------
    
    # The basis must be compatible with the first-stage reduced coefficients.
    # Therefore, we use the same spline function and the same knot placement
    # strategy: internal knots at the local 10th and 90th percentiles.
    
    bvar_i <- onebasis(
      temps,
      fun = "ns",
      knots = quantile(data_igr[[temp_var]], c(0.10, 0.90), na.rm = TRUE),
      Bound = range(data_igr[[temp_var]], na.rm = TRUE)
    )
    
    #----------------------------------------------------------------------
    # Step 4. Predict the BLUP curve and locate the MMT
    #----------------------------------------------------------------------
    
    coef_i <- blup_all[[igr_code]]$blup
    
    rr_i <- as.numeric(
      exp(bvar_i %*% coef_i)
    )
    
    mmt_position <- which.min(rr_i)
    
    mmt_igr$MMT_perc[i] <- perc_seq[mmt_position]
    mmt_igr$MMT[i] <- as.numeric(temps[mmt_position])
    mmt_igr$status[i] <- "success"
    
  }, error = function(e) {
    
    mmt_igr$status[i] <<- paste("error:", e$message)
    
    message("Error estimating MMT for IGR ", igr_code, ": ", e$message)
  })
}


#==============================================================================
# 10. Save second-stage results
#==============================================================================

save(
  meta_model,
  cp_meta,
  pooled_curve,
  MMT_pooled,
  blup_all,
  mmt_igr,
  file = file.path(dir_save, "second_stage_results.RData")
)

write.csv(
  mmt_igr,
  file = file.path(dir_save, "MMT_by_IGR_after_blup.csv"),
  row.names = FALSE
)

################################################################################
# End of script
################################################################################

