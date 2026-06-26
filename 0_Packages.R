################################################################################
#
#  Running Health Models in R
#
#  Script 00 - Install and Load Packages
#
################################################################################

#------------------------------------------------------------------------------
# Purpose
#------------------------------------------------------------------------------
#
# This script installs (if necessary) and loads all R packages required for
# the course "Running Health Models in R".
#
# The script is intended to be executed only once when setting up a new R
# environment. Subsequent executions will simply load the installed packages.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Required packages
#------------------------------------------------------------------------------

packages <- c(
  "arrow",
  "data.table",
  "dlnm",
  "dplyr",
  "ggplot2",
  "lubridate",
  "mixmeta",
  "stringr",
  "tidyr"
)

#------------------------------------------------------------------------------
# Install missing packages
#------------------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

invisible(lapply(packages, install_if_missing))

#------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------

invisible(lapply(packages, library, character.only = TRUE))

#------------------------------------------------------------------------------
# Base R packages
#------------------------------------------------------------------------------

# 'splines' is distributed with R and therefore does not need installation.
library(splines)

#------------------------------------------------------------------------------
# End of script
#------------------------------------------------------------------------------