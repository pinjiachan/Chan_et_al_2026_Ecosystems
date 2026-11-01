# Code for Chan et al. (2026) Ecosystems

Welcome! This is a repository for the code for the litterfall model in our paper:

> Chan PJ, Hu K, Elliott T, Lai HR, Burns B, Schwendenmann L (2026). Dieback disease increases short-term spatial heterogeneity and litterfall inputs in forests dominated by a long-lived foundation species. *Ecosystems* 29:91. https://doi.org/10.1007/s10021-026-01113-8

## Usage

The scripts are in the `R` folder:

- `01_fit_model.R` fits the Bayesian litterfall model in `greta`
- `02_map_predictions.R` uses the fitted model to predict litterfall on a 1 m grid in each stand, using the observed canopy health scores of kauri

## Data

The data are not included. They were collected on the traditional lands of Te Kawerau a Maki and are held under the CARE Principles for Indigenous Data Governance. As stated in the paper, requests to access the data should be directed to Luitgard Schwendenmann (l.schwendenmann@auckland.ac.nz). The code is shared to document the analysis, so it cannot be run without the data.

## Notes

- The map predictions were scripted plot by plot because of computational limitations. The stand is set in the script, and the script is run once for each stand.
- `greta` needs a TensorFlow backend, which can be installed with `greta::install_greta_deps()`.


