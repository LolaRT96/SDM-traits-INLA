# -----------------------------
# Title: Project Setup for Traits in Space Pipeline
# Author: M.Grazia Pennino
# Date: 2025-07-08
# Description: Install/load packages, create project folder structure,
#              and define base paths for data, output, and plots.
# -----------------------------

# --- 1. Install and load required packages --------------------------------------

required_pkgs <- c("icesDatras", "dplyr", "purrr", "tidyverse", "INLA",
                   "sp", "fields", "viridis", "pROC", "ggplot2")

installed_pkgs <- rownames(installed.packages())

for (pkg in required_pkgs) {
  if (!pkg %in% installed_pkgs) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Note: It is recommended to install INLA from the official website for the latest version:
# https://www.r-inla.org/download-install

# --- 2. Create folder structure --------------------------------------------------

dir.create("data", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# --- 3. Define base paths --------------------------------------------------------

data_path <- "data"
output_path <- "output"
plots_path <- "plots"

# --- 4. Print message -------------------------------------------------------------

message("Project setup complete.")
message("Folders created: data/, output/, plots/")
message("Paths defined:")
message(" - Data directory: ", data_path)
message(" - Output directory: ", output_path)
message(" - Plots directory: ", plots_path)

