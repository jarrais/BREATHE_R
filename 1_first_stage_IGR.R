################################################################################
#
#  BREATHE
#
#  Running Health Models in R
#
#  Script 01 - First-stage time-series models by IGR
#
################################################################################

#------------------------------------------------------------------------------
# Purpose
#------------------------------------------------------------------------------
#
# This script fits the first-stage time-series models separately for each
# Immediate Geographic Region (IGR).
#
# The aim is to estimate the IGR-specific temperature-mortality association
# using a Distributed Lag Non-Linear Model (DLNM).
#
# For teaching purposes, we adopt a fixed model specification:
#
#   - maximum lag: 5 days
#   - temperature-response: natural cubic spline with two internal knots
#     located at the 10th and 90th percentiles
#   - lag-response: natural cubic spline with three log-spaced knots
#   - seasonal and long-term adjustment: natural spline with 4 degrees of
#     freedom per year
#
# In a full applied analysis, alternative parameterizations can be compared
# using the quasi-Akaike Information Criterion (qAIC), together with convergence,
# epidemiological plausibility, spatial coverage, and model stability.
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
# 2. Load and prepare data
#==============================================================================

# Mortality and temperature data have already been aggregated at the
# Immediate Geographic Region (IGR) level. Both datasets contain one
# observation per IGR and calendar day.

# Load mortality data
mort_data <- read_parquet("mort_data_RJ.parquet")

# Load daily mean temperature data
tempdata <- read_parquet("tempdata_RJ.parquet")


# Merge mortality and temperature datasets using the IGR identifier
# and calendar date.
data_final <- inner_join(
  x  = mort_data,
  y  = tempdata,
  by = c("IGR", "Date")
)


# Prepare variables required for the time-series analysis.
#
# The day of the week is included as a categorical variable to adjust
# for systematic weekly patterns in mortality. A sequential time index
# is also created within each IGR and will later be used to control
# for long-term and seasonal trends.

data_final <- data_final |>
  mutate(
    Date = as.Date(Date),
    # Day of the week controls for weekly mortality patterns
    dow = factor(
      weekdays(Date),
      levels = c(
        "Monday", "Tuesday", "Wednesday", "Thursday",
        "Friday", "Saturday", "Sunday"
      )
    )
  ) |>
  arrange(IGR, Date) |>
  group_by(IGR) |>
  mutate(
    # Sequential time index used to control long-term and seasonal trends
    time = row_number()
  ) |>
  ungroup()

#==============================================================================
# 3. Define model specification
#==============================================================================

# Maximum lag in days
lag_num <- 5

# Natural cubic spline for both temperature and lag dimensions
lagfun <- "ns"
varfun <- "ns"

# Internal knots for the temperature-response curve
knots_perc <- c(0.10, 0.90)

# Internal knots for the lag-response curve
knots_lag <- logknots(lag_num, nk = 3)

# Degrees of freedom per year for long-term and seasonal control
dfseas <- 4

# Number of coefficients after reducing the DLNM to the overall association
num_cols <- length(knots_perc) + 1

# IGRs included in the analysis
igr_include <- unique(data_final$IGR)

# Common temperature range used as boundary for the spline basis
bound_temp <- range(
  data_final$daily_pop_weighted_mean_temperature,
  na.rm = TRUE
)

#==============================================================================
# 4. Create objects to store first-stage results
#==============================================================================

coefs_all_mat <- matrix(
  NA,
  nrow = length(igr_include),
  ncol = num_cols,
  dimnames = list(igr_include, NULL)
)

vcov_all <- setNames(
  vector("list", length(igr_include)),
  igr_include
)

MMT_vec <- rep(NA, length(igr_include))

status_vec <- rep(NA_character_, length(igr_include))

#==============================================================================
# 5. Stage 1 - Fit IGR-specific time-series models
#==============================================================================

# The first stage consists of fitting an independent time-series model for
# each Immediate Geographic Region (IGR). The resulting exposure-response
# associations are then used in the second-stage multivariate meta-analysis.
#
# The following steps are repeated for every IGR:
#
#   1. Subset the data for the current IGR.
#   2. Construct the DLNM cross-basis.
#   3. Fit the quasi-Poisson regression model.
#   4. Reduce the bi-dimensional association to the overall cumulative curve.
#   5. Store the estimated coefficients and covariance matrix.
#   6. Estimate the Minimum Mortality Temperature (MMT).

for (i in seq_along(igr_include)) {
  
  igr_code <- igr_include[i]
  
  cat("Fitting IGR", i, "of", length(igr_include),
      "- IGR:", igr_code, "\n")
  
  tryCatch({
    
    #----------------------------------------------------------------------
    # Step 1. Extract data for the current IGR
    #----------------------------------------------------------------------
    
    data_igr <- data_final[data_final$IGR == igr_code, ]
    
    # Number of years available for this IGR.
    # This is used to determine the total degrees of freedom for the
    # seasonal and long-term trend adjustment.
    n_years <- length(unique(format(data_igr$Date, "%Y")))
    
    #----------------------------------------------------------------------
    # Step 2. Construct the DLNM cross-basis
    #----------------------------------------------------------------------
    
    # Internal knots are placed at the 10th and 90th percentiles of the
    # local temperature distribution. This allows the exposure-response
    # relationship to adapt to the climate experienced within each IGR.
    
    knots_temp <- quantile(
      data_igr$daily_pop_weighted_mean_temperature,
      probs = knots_perc,
      na.rm = TRUE
    )
    
    cb <- crossbasis(
      data_igr$daily_pop_weighted_mean_temperature,
      lag = lag_num,
      argvar = list(
        fun   = varfun,
        knots = knots_temp,
        Bound = bound_temp
      ),
      arglag = list(
        fun   = lagfun,
        knots = knots_lag
      )
    )
    
    #----------------------------------------------------------------------
    # Step 3. Fit the first-stage time-series model
    #----------------------------------------------------------------------
    
    # The model assumes a quasi-Poisson distribution for the daily mortality
    # counts and adjusts for day of the week and long-term/seasonal trends.
    
    model <- glm(
      count_all ~ cb + dow + ns(time, df = dfseas * n_years),
      family = quasipoisson(),
      data = data_igr
    )
    
    #----------------------------------------------------------------------
    # Step 4. Reduce the DLNM association
    #----------------------------------------------------------------------
    
    # The cross-basis describes both the exposure-response and lag-response
    # dimensions. For the second-stage meta-analysis, this bi-dimensional
    # association is reduced to the overall cumulative exposure-response
    # relationship across all lag days.
    
    red <- crossreduce(
      cb,
      model,
      by = 0.1,
      type = "overall"
    )
    
    
    #----------------------------------------------------------------------
    # Step 5. Store first-stage estimates
    #----------------------------------------------------------------------
    
    # Save the reduced coefficients and their variance-covariance matrix.
    # These quantities constitute the input for the second-stage
    # multivariate meta-analysis.
    
    coefs_all_mat[i, ] <- coef(red)
    vcov_all[[i]] <- vcov(red)
    
    
    #----------------------------------------------------------------------
    # Step 6. Estimate the Minimum Mortality Temperature (MMT)
    #----------------------------------------------------------------------
    
    # The MMT is defined as the temperature associated with the lowest
    # predicted cumulative relative risk.
    
    pred <- crosspred(
      cb,
      model,
      by = 0.1
    )
    
    MMT_vec[i] <- pred$predvar[
      which.min(pred$allRRfit)
    ]
    
    status_vec[i] <- "success"
    
  }, error = function(e) {
    
    status_vec[i] <<- paste("error:", e$message)
    
    message("Error in IGR ", igr_code, ": ", e$message)
  })
}

#==============================================================================
# 6. Organise first-stage summary
#==============================================================================

results_igr <- data.frame(
  IGR    = igr_include,
  MMT    = MMT_vec,
  status = status_vec,
  stringsAsFactors = FALSE
)

results_igr

#==============================================================================
# 7. Save first-stage results
#==============================================================================

save(
  coefs_all_mat,
  vcov_all,
  igr_include,
  results_igr,
  file = file.path(dir_save, "first_stage_results.RData")
)

################################################################################
# End of script
################################################################################








##############################################################################
################################  B R E A T H E  #############################
##############################################################################
# Climate-change attribution (Brazil) –
# Two-stage time-series model by RGI
##############################################################################

# ========================================================================== #
# 1. Directories and file tags
# ========================================================================== #
#setwd("~/Breathe/Data/RGI/")
dir_save <- "/Users/jonyarrais/Documents/UFF/Projetos/BREATHE/Artigo Attibution LANCET/2 - BREATHE RGI SIM/Curso Bristol/RESULTS"

outcome_tag <- "allpop_"
region_tag  <- "RGI_"
model_tag   <- "ts_SIM_"
file_tag    <- paste0(model_tag, region_tag, outcome_tag)

# ========================================================================== #
# 2. Load data already aggregated by RGI
# ========================================================================== #
#Dados de mortalidade
mort_data <- read_parquet("mort_data_RJ.parquet")

tempdata <- read_parquet("tempdata_RJ.parquet")

data_final <- inner_join(x = mort_data, 
                         y = tempdata,
                         by = c("IGR","Date"))


data_final <- data_final %>%
  mutate(Date = as.Date(Date),
         dow  = factor(weekdays(Date),
                       levels = c("Monday","Tuesday","Wednesday",
                                  "Thursday","Friday","Saturday","Sunday"))) %>%
  arrange(IGR, Date) %>%
  group_by(IGR) %>%
  mutate(time = row_number()) %>%
  ungroup()

# ========================================================================== #
# 3. Model parameters
# ========================================================================== #
lag_num    <- 5
lagfun     <- "ns"
varfun     <- "ns"
knots_perc <- c(0.1, 0.9)
knots_lag  <- logknots(lag_num, nk = 3)
dfseas     <- 4
num_cols   <- length(knots_perc) + 1
igr_include <- unique(data_final$IGR)

# ========================================================================== #
# 4. Stage 1 – RGI-specific models
# ========================================================================== #
coefs_all_mat <- matrix(NA, nrow = length(igr_include), ncol = num_cols,
                        dimnames = list(igr_include, NULL))
vcov_all  <- setNames(vector("list", length(igr_include)), igr_include)

MMT_vec = rep(NA, length(igr_include))

bound_temp <- range(data_final$daily_pop_weighted_mean_temperature, na.rm = TRUE)

for (i in seq_along(igr_include)) tryCatch({
  
  cat(i, " ")
  
  igr_code   <- igr_include[i]
  data_igr   <- data_final[data_final$IGR == igr_code, ]
  n_years    <- length(unique(format(data_igr$Date, "%Y")))
  
  knots_temp <- quantile(data_igr$daily_pop_weighted_mean_temperature, probs = knots_perc, na.rm = TRUE)

  cb <- crossbasis(
    data_igr$daily_pop_weighted_mean_temperature,
    lag    = lag_num,
    argvar = list(fun = varfun, knots = knots_temp, Bound = bound_temp),
    arglag = list(fun = lagfun, knots = knots_lag))
  
  #modelo sem umidade
  model <- glm(count_all ~ cb + dow + ns(time, df = dfseas * n_years),
               family = quasipoisson(),
               data   = data_igr)
  
  red <- crossreduce(cb, model, by = 0.1, type = "overall")
  coefs_all_mat[i, ] <- coef(red)
  vcov_all [[i]]     <- vcov(red)
  
  pred     <- crosspred(cb, model, by = 0.1)
  MMT      <- pred$predvar[ which.min(pred$allRRfit) ]
  MMT_vec[i] <- MMT
  })

# ========================================================================== #
# 5. Save
# ========================================================================== #
results_igr <- data.frame(
  IGR   = igr_include,
  MMT       = MMT_vec,
  stringsAsFactors = FALSE
)

save(coefs_all_mat, vcov_all, igr_include,
     file = file.path(dir_save, paste0(file_tag, "first_stage_lag5_knots_p10_90_df4.RData")))
