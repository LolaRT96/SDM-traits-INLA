# Traits in Space: Integrating Functional Ecology into Spatial Species Distribution Models
This repository contains the R scripts and data pipeline to simulate and analyze spatial species distribution models (SDMs) integrating functional traits, applied both to simulated and empirical ICES DATRAS data. 
The aim is to demonstrate how including species traits improves model performance and ecological inference.

---
## Repository Structure

- `00_project_setup.R`  
Installs and loads all required packages and creates the project folder structure (`data/`, `output/`, `plots/`).
- `01_download_datrast_data.R`  
  Downloads haul header (HH) and catch-at-age (CA) data from ICES DATRAS for selected surveys and years, harmonizes column types, and saves them as RDS files.
- `02_prepare_empirical_sdm_data_merluccius.R`  
  Loads the downloaded ICES data, filters for the target species (*Merluccius merluccius*), merges with FishBase traits, and prepares the dataset for spatial SDM modeling.
- `03_fit_simulated_case_study.R`  
  Fits a spatial SDM using INLA + SPDE to the simulated data, including environmental covariates and traits, evaluates model performance, and produces diagnostic plots.
- `04_fit_spatial_sdm_merluccius_merluccius.R`  
   Fits a spatial SDM using INLA + SPDE to the case study data, including environmental covariates and traits, evaluates model performance, and produces diagnostic plots.
- `05_prepare_multispecies_datasets.R`  
  Loads the downloaded ICES data for the multispecies models
- `06_fit_multi_species_model.R`  
  Fits a spatial SDM using INLA + SPDE to the multispecies case study data, including environmental covariates and traits, evaluates model performance, and produces diagnostic plots.
- `07_simulated_data_case_multispecies.R`  
  Simulates a multi-species dataset with species-specific random slopes modulated by traits, fits a hierarchical random slopes model using INLA, evaluates results, and visualizes trait effects.
  

## Requirements
- R version >= 4.2  
- Packages: `icesDatras`, `dplyr`, `purrr`, `tidyverse`, `INLA`, `sp`, `fields`, `viridis`, `pROC`, `ggplot2`, `rfishbase`, `sf`

Note: INLA is recommended to be installed via their official website for the latest stable version: [www.r-inla.org](https://www.r-inla.org/download-install)

---

## Contact

For questions or suggestions, please contact:  
Mara Dolores Riesgo mdolores.riesgo@ieo.csic.es

---

## License

This repository is open-source under the MIT License.

