# Supporting Code for: Is the denitrification potential of floodplain sediments controlled by their organic carbon contents?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16903197.svg)](https://doi.org/10.5281/zenodo.16903197)


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> AmmerBatch

It is authored by Vitor Cantarella, Johann Holdt and coauthors.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "AmmerBatch"
```
which auto-activate the project and enable local path handling from DrWatson.

## Folder Structure

```
.
├── data
│   ├── exp_pro
│   ├── exp_raw
│   └── external
├── notebooks
├── plots
├── scripts
├── src
└── test
```

## File Descriptions

### `/scripts`

- **`reading_processing_plotting_raw_data.jl`**: This script reads the raw experimental data from Excel files (`.xlsx`) located in `data/exp_raw`. It processes the data, converts units, and saves the processed data into individual CSV and JLS files for each sample in `data/exp_pro`. It also generates initial concentration plots for each sample and saves them in the `plots` folder.
- **`fit_linear_regression_v2.jl`**: This script performs a linear regression on the processed data to determine reaction rates. It reads the processed data from `data/exp_pro`, calculates zero-order reaction rates for nitrate, and saves the results in `data/exp_pro/linear_regression_params_v2.csv`. It also generates grid plots summarizing the fits, which are saved in the `plots` folder.
- **`analysing_rates_with_flow_v2.jl`**: This script combines the calculated reaction rates with external data on groundwater flow to estimate the length required for denitrification. It reads the regression results from `data/exp_pro` and external data from `data/external`. The resulting plots are saved in the `plots` folder.
- **`analysing_results_v2sep.jl`**: This script generates plots for the publication, summarizing the results from the linear regression. It reads the regression results from `data/exp_pro/linear_regression_params_v2.csv` and produces plots that are saved in the `plots` folder. It also saves a summary of the results by facies in `data/exp_pro/facies_results.csv`.
- **`sulfur_analysis.jl`**: This script analyzes the relationship between nitrate and sulfate reduction rates. It reads the regression results from `data/exp_pro/linear_regression_params_v2.csv` and generates plots showing the stoichiometric relationship between the two processes. The plots are saved in the `plots` folder.
- **`make_plotexplanation_background.jl`**: This script generates a specific plot for a single sample (A308) to be used as a background for a figure in the manuscript. It reads the data for this sample from `data/exp_pro` and saves the resulting plot in various formats in the `plots` folder.

### `/data/exp_raw`

- **`johann_batch_preprocessed.xlsx`**: Raw data from the batch experiments (Part 1). Contains information about the samples, time points, and measured concentrations of various chemical species.
- **`johann_batch_preprocessed_part2.xlsx`**: Raw data from the batch experiments (Part 2).
- **`MK_Trockengewichte_Johann .xlsx`**: Contains the dry weights of the samples.

### `/data/exp_pro`

This folder contains the processed data, derived from the raw data by the `reading_processing_plotting_raw_data.jl` script.

- **`sample_info.csv`**: Information about the samples from Part 1 of the experiments.
- **`sample_info_part2.csv`**: Information about the samples from Part 2 of the experiments.
- **`$(sample).csv` / `$(sample).jls`**: Processed data for each individual sample, with concentrations converted to mmol/L.
- **`linear_regression_params_v2.csv`**: Results of the linear regression analysis, including calculated reaction rates (`r_no3`, `r_so4`), R² values, and other parameters.
- **`facies_results.csv`**: Summary of the reaction rates grouped by facies.

### `/data/external`

- **`velocity_facies_v2.csv`**: External data containing information about porosity and hydraulic conductivity for each facies. Reference Holdt et al. 2025
- **`q_facies.csv`**: External data containing information about specific discharge for each facies. Reference Holdt et al. 2025

### `/plots`

This folder contains all the plots generated by the scripts.

### `/notebooks`

- **`model_formulation.qmd`**: A Quarto notebook containing the model formulation. This can be rendered to HTML or PDF.

### `/src` and `/test`

These folders contain the source code for the project and the tests, respectively.

## Description of the relevant scripts

Please refer to the following scripts to reproduce the works of this paper:
all scripts locate in the subfolder `/scripts`

- reading_processing_plotting_raw_data: data preparation for analysis from the excel sheet with lab results to mmol/L concentrations for each sample.
- fit_linear_regression_v2: data correction and analysis through linear regression as described in the Methods section
- analysing_results_v2sep: making plot for the publication from the results of the linear regression
- analysing_rates_with_flow_v2: makes the last plot of the paper joining the reaction rates from this study with the specific discharge data from Holdt et al. 2025 (groundwater model results from a previous work).
- make_plotexplanation_background: script to draw the background used in Figure 2 of the manuscript. The figure is further edited in excalidraw (raw file is in the subfolder plots).
