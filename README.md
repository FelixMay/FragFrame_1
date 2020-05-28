
<!-- README.md is generated from README.Rmd. Please edit that file -->
FragFrame\_1
============

<!-- badges: start -->
<!-- badges: end -->
The FragFrame\_1 repository includes R code and data to reproduce the analyses shown in the article

**Ecosystem decay exacerbates biodiversity loss with habitat loss**

*by Jonathan M. Chase, Shane A. Blowes, Tiffany M. Knight, Katharina Gerstner and Felix May*

Here we give a brief overview on the code and data files in this repository

Data
----

There are two data files in data folder:

`fragSAD_predicts_ewers.csv`: This file contains the species abundance data and fragment areas for all studies

`new_meta_2_merge.csv`: This file contains the meta data for all studies

Extendend data figs tabs
------------------------

This folder contains the figures and tables that appear in the supplementary material of the article. It also provides the model fits and simulation results that are required to reproduce the figures.

Intermediate results
--------------------

This folder contains the estimates of all biodiversity indices for all fragments (*biodiv\_frag*) and of beta-diversity partitioning values for all fragment pairs (*betapart*). The files with different numbers refer to different methods for the imputation of the area of continuous habitat and the treatment of non-integer abundances. The reference scenario is labelled with 2\_, while the other versions provide robustness checks.

The file `7_resultsS2000_N40000_mp1_nrep2000.csv` contains the simulation results for the random sampling hypothesis presented in Ext\_Dat\_Fig1.

Main results
------------

There are multiple data files in this folder:

`fragSize_ref.Rdata`: This file contains the models fit to all the data for each metric.

`fragSize_ref_pool.Rdata`: This file contains the models fit to the data with the studies 'pooled' sampling designs removed.

`fragSize_z_score_ref.Rdata`: This file contains the models fit to the z-score transformed biodiversity indices.

`Jne_zi_fragSize.Rdata`: contains the model fit to the nestedness component of Jaccard dissimilarity as a function of patch size difference.

`Jtu_z1i_fragSize.Rdata`: contains the model fit to the turnover component of Jaccard dissimilarity as a function of patch size difference.

`Rne_zi_fragSize.Rdata`: contains the model fit to the nestedness component of Ruzicka dissimilarity as a function of patch size difference.

`Rth_z1i_fragSize.Rdata`: contains the model fit to the turnover component of Ruzicka dissimilarity as a function of patch size difference.

This folder also contains the figures presented in the main text of the article.

Pre-processing files
--------------------

`0_preproc_wide_to_long_format.R`, `0_preproc_ewers_et_al_2007.R`, `0_preproc_predicts_data.R`: These three files do the pre-processing of raw data. Essentially, these three files create the data file `fragSAD_predicts_ewers.csv`, which is published in this repository. The single study data files that are the input for the pre-processing files are not published here and therefore, this files cannot be re-run as they are.

R scripts
---------

`0_init_dirs_load_packages.R`: This script loads necessary R packages and sets up path and working directories. This set up needs to be adjusted for specific users and R sessions.

`1_analysis_get_biodiv_new.R`: This file implements the first main part of the analysis, which includes the standardization and estimation of biodiversity parameters and the community dissimilarity partitioning. Based on the data in `fragSAD_predicts_ewers.csv`, this R-script calculates all biodiversity parameters for all fragments, as well as the incidence- and abundance-based turnover and nestedness components for all fragment pairs.

`2_`: Preliminary visual inspection of the data.

`3_`: These files do the second main part of the analysis by fitting models to the standardized data. The output of `3a_fragSize_fit_models.R` is save in the main results folder as `fragSize_ref.Rdata` and `fragSize_ref_pool.Rdata`; similarly, the output of `3d_` (models fit to beta-diversity components) is saved to individual files found in the main results folder.

`4_`: These files contain code to conduct visual inspection of model fits (e.g., plots of residuals, Markov chain inspection, posterior predictive checks).

`5_`: These files prepare coefficients and posterior distributions for plotting.

`6_`: These files plot results as they are presented in the main manuscript text.

`7_`: These files plot results as they are presented in the extended data section. These files also provide the simulations and model fits required for the extended data section.

`99_gg_legend`: Function for plotting legends as separate panels.
