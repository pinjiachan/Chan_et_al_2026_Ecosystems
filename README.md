# Code for Chan et al. (2026) Ecosystems

This is a repository for the code for the litterfall model in our paper:

> Chan PJ, Hu K, Elliott T, Lai HR, Burns B, Schwendenmann L (2026). Dieback disease increases short-term spatial heterogeneity and litterfall inputs in forests dominated by a long-lived foundation species. *Ecosystems* 29:91. https://doi.org/10.1007/s10021-026-01113-8

## Usage

The scripts are in the `R` folder:

- `01_fit_model.R` fits the Bayesian litterfall model in `greta`
- `02_map_predictions.R` uses the fitted model to predict litterfall on a 1 m grid in each stand, using the observed canopy health scores of kauri

## Data

The data are not included. They were collected on the traditional lands of Te Kawerau a Maki and are held under the CARE Principles for Indigenous Data Governance. As stated in the paper, requests to access the data should be directed to Luitgard Schwendenmann (l.schwendenmann@auckland.ac.nz). The code is shared to document the analysis, so it cannot be run without the data.

### Data format
 
The scripts read a single file, `data.csv`. Each row is one species-specific litter observation (the leaf litter of one species collected from one trap in one monthly collection) paired with one tree of that species within 25 m of the trap. The litter columns are repeated on every tree row of the same litter observation, and the tree columns are repeated for every trap and collection that a tree is paired with. The scripts collapse the data to one row per litter observation and one row per tree where needed.

| Column | Description |
|---|---|
| Plot | Forest plot (stand) ID |
| Collection | Collection number (monthly collections) |
| Date | Collection date |
| Trap | Trap number within the plot |
| Plot_Trap | Unique trap ID (plot and trap) |
| X_trap, Y_trap | Trap coordinates within the plot (m) |
| Species_code | Species code of the litter and of the tree |
| Mass | Dry mass of the litter of that species in the trap (g) |
| Interval_days | Days since the previous collection |
| Sieve_size | Sieve mesh size, "1mm" or "4mm" |
| tree_id | Tree ID |
| X_tree, Y_tree | Tree coordinates within the plot (m) |
| DBH | Tree diameter at breast height (cm) |
| Canopy | Canopy health score (1 to 4 for kauri; 1 for all other species) |
| dist | Distance between the tree and the trap (m) |

## Notes

- The map predictions were scripted plot by plot because of computational limitations. 
- `greta` needs a TensorFlow backend, which can be installed with `greta::install_greta_deps()`.


