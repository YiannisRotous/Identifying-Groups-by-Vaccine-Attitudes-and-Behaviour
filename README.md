# Identifying-Groups-by-Vaccine-Attitudes-and-Behaviour
This repository presents the findings and methodology from a 2024 community-based survey conducted in London, UK. Moving beyond the "pro-vaccine" vs. "anti-vaccine" binary, this study characterizes the multidimensional nature of vaccine hesitancy and uptake within a diverse urban population.

# Data Availability

The dataset used in this project is not included in this repository as it was not collected by the authors. The data originates from the pan-London Why We Get Vaccinated campaign. The pan-London WWGV campaign is a three-year campaign created with local communities, UKHSA, the NHS, London Borough Councils and the London regional Public Health communications network of the Association of Directors of Public Health (ADPH) London; the organisation acts as a coordinating body, championing public health across all of London’s 32 boroughs and the City of London, and providing a forum for the Directors of Public Health (DsPH) to collaborate.

# Data Requests

We do not have the authority to distribute this data. Any requests for access or technical inquiries regarding the dataset should be directed to the Association of Directors of Public Health (ADPH) London or the relevant campaign coordinators.

For more information, please visit the [ADPH London official website](https://www.adph.org.uk/networks/london/)

# Analysis & Scripts

The statistical analysis for the separate regression models for each question and statement, is contained in the R Markdown file:

* **[separate_regressions_wave1.Rmd](./R_scripts/separate_regressions_wave1.Rmd)**: This R Markdown script contains the code for the separate regression analysis of the survey questions and statements.
* **[preprocess_survey_wave.R](./R_scripts/preprocess_survey_wave.R)**: This script is used to clean and preprocess the raw survey data, bringing it into the required format for statistical modeling.

For the joint Bayesian Nonparametric latent class analysis, we utilized a Stan-based workflow:

* **[BNP_LCA_analysis.R](./R_scripts/BNP_LCA_analysis.R)**: The main R script used to execute the analysis.
* **[lca.stan](./Stan_scripts/lca.stan)**: Contains the Stan code for the latent class analysis model.
* **[two_step_mult.stan](./Stan_scripts/two_step_mult.stan)**: Contains the Stan code for the two-step approach, where predicted classifications are regressed on individual-level covariates.
  
