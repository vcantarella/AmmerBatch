# Supporting Code for: Is the denitrification potential of floodplain sediments controlled by their organic carbon contents?


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> AmmerBatch

It is authored by Vitor Cantarella, Johann Holdt.

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

## Description of the relevant scripts

Please refer to the following scripts to reproduce the works of this paper:
all scripts locate in the subfolder `/scripts`

- reading_processing_plotting_raw_data: data preparation for analysis from the excel sheet with lab results to mmol/L concentrations for each sample.
- fit_linear_regression_v2: data correction and analysis through linear regression as described in the Methods section
- analysing_results_v2sep: making plot for the publication from the results of the linear regression
- analysing_rates_with_flow_v2: makes the last plot of the paper joining the reaction rates from this study with the specific discharge data from Holdt et al. 2025 (groundwater model results from a previous work).
- make_plotexplanation_background: script to draw the background used in Figure 2 of the manuscript. The figure is further edited in excalidraw (raw file is in the subfolder plots).
