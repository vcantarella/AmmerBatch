#=
This script analyzes the reaction rates obtained from the `fit_linear_regression_v2.jl`
script and combines them with external data on groundwater flow to estimate the length
required for denitrification within different geological facies.

The main steps of the script are:
1.  **Load Data**:
    - It loads the linear regression results, which include the zero-order nitrate
      reduction rates (`r_no3`), from `data/exp_pro/linear_regression_params_v2.csv`.
    - It loads external data on porosity (`velocity_facies_v2.csv`) and specific
      discharge (`q_facies.csv`) for each facies from the `data/external` directory.
2.  **Calculate Physical Properties**:
    - It calculates the solid density (`ρₛ`) and bulk density (`ρᵦ`) of the sediments
      based on Total Organic Carbon (TOC) content.
3.  **Calculate Denitrification Length**:
    - It calculates the theoretical length (`length_to_reduce`) required to completely
      denitrify a given initial concentration of nitrate (50 mg/L) for each sample.
      This calculation uses the measured reaction rates and physical properties, along with
      a calculated specific discharge (`q`).
    - It also calculates a `length_model` using specific discharge values from a model.
4.  **Group by Facies**:
    - The results are grouped by facies, and summary statistics (mean, std, min, max)
      for the denitrification length are calculated for each facies.
5.  **Generate Plots**:
    - It generates several plots to visualize the results:
        - `length_to_reduce_v2.png`: A bar plot showing the mean denitrification length
          for each facies, with error bars representing the min and max values.
        - `length_to_reduce_model_v2.png` and `length_to_reduce_model_v3.png`:
          Similar plots, but using the model-based specific discharge and with different styling.
    - These plots are saved in the `plots` directory.

The script's output helps to understand how the denitrification potential, represented by
the reaction rates, translates into an effective denitrification capacity at the field scale,
considering the hydrogeological properties of the different facies.
=#

using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = DataFrame(CSV.File(datadir("exp_pro","linear_regression_params_v2.csv")))
int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"

# Read the velocity data
vel_df = CSV.read(datadir("external","velocity_facies_v2.csv"), DataFrame)
q_df = CSV.read(datadir("external","q_facies.csv"), DataFrame)


TOC = int_df[!,:TOC].*1e-2 # %wt (proportion) total organic carbon
C_OM = 37.5*1e-2 #%wt (proportion) carbon content in organic matter
ρₒₘ = 1350 #kg/m3 density of organic matter
ρₘₛ = 2710 #kg/m3 density of calcite
# ρₛ calculated based on expression from Ruehlmann et al. 2006
ρₛ = @. (TOC/C_OM/ρₒₘ + (1-TOC/C_OM)/ρₘₛ)^-1 #kg/m3
# porosity from facies as estimated in Strobel et al 2024
# get the values in vel_df and match with the facies in int_df
ϕ = [vel_df[vel_df.facies .== facies, :porosity][1] for facies in int_df.facies]
# bulk density:
ρᵦ = @. (1-ϕ)*ρₛ
# specific discharge given in m/s from Martin et al. 2020
T = 6.7e-5 # transmissivity m2/s
thick = 5 # thickness m
i = 2 / 750 # hydraulic gradient
q = T * i / thick # specific discharge m/s
q = q * 86400 # m/d
q_model = [q_df[q_df.facies .== facies, :mean][1] for facies in int_df.facies]
q_model .= q_model .* 86400 # m/d
# length for denitrifying a given concentration
c_no3 = 50e-3 / 62.0049 # mol L⁻¹ - given concentration inflow of NO₃⁻
c_no3 = c_no3 * 1e3 # mol m⁻³
length_to_reduce = c_no3 * q ./ (ρᵦ .* int_df[!, :r_no3]) # m
int_df[!, :length] = length_to_reduce
length_model = c_no3 * q_model ./ (ρᵦ .* int_df[!, :r_no3]) # m
int_df[!, :length_model] = length_model
k_model = [q_df[q_df.facies .== facies, :k][1] for facies in int_df.facies]

v = 1:nrow(int_df)
ind = v[.!ismissing.(int_df[!, :TOC])]
facies_result_len = @chain(
    # filter missing TOC values
    int_df[ind, :],
    @group_by(facies),
    @summarize(
        mean_length = mean(length),
        std_length = std(length),
        min_length = minimum(length),
        max_length = maximum(length),
    )
)

facies_result_model = @chain(
    # filter missing TOC values
    int_df[ind, :],
    @group_by(facies),
    @summarize(
        mean_length = mean(length_model),
        std_length = std(length_model),
        min_length = minimum(length_model),
        max_length = maximum(length_model),
    )
)


unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result_len.facies_code = map(facies->facies_code[facies],facies_result_len.facies)
facies_result_model.facies_code = map(facies->facies_code[facies],facies_result_model.facies)
facies_result_len.k_result = [q_df[q_df.facies .== facies, :k][1] for facies in facies_result_len.facies]
facies_result_model.k_result = [q_df[q_df.facies .== facies, :k][1] for facies in facies_result_model.facies]
# sort dfs based on facies_code
facies_result_len = sort(facies_result_len, :facies_code)
facies_result_model = sort(facies_result_model, :facies_code)
label_values = ["Clay (FT-C1)", "Tufa grains (FT-01)", "Calcareous silt (FT-02)", "Tufa & reed (FT-04)", "Silt & moss (FT-06)", "Silt & organic debris (FT-07)",
    "Brown peat (FT-08)", "Black peat (FT-09)"]
labels = Dict(zip(unique_facies, label_values))
fontcolor = "#474747"
f = Figure()
ax = Axis(f[1, 1],
    # xlabel = "Facies",
    ylabel = "length to reduce 50 mg L⁻¹ NO₃⁻ [m]",
    xticks = (1:8, label_values),
    # yticks = 1e-1:2e-1:1.2,
    xticklabelrotation = π/6,
    xgridvisible = false,
    ygridvisible = false,
    # title = "NO₃⁻ Reduction Lengths Across Facies Types",
    # titlefont = "Avenir Book",
    titlesize = 21,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    xlabelfont = "Avenir Book",
    ylabelfont = "Avenir Book",
    xlabelcolor = fontcolor,
    ylabelcolor = fontcolor,
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    )
ax2 = Axis(f[1, 1],
    ylabel = "k [m/s]",
    yaxisposition = :right,
    xgridvisible = false,
    ygridvisible = false,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
)
hidespines!(ax2)
hidedecorations!(ax2)
hidespines!(ax, :t, :r)
barplot!(ax, facies_result_len.facies_code, facies_result_len.mean_length, color = :steelblue)
errorbars!(ax, facies_result_len.facies_code, facies_result_len.mean_length,
    facies_result_len.mean_length-facies_result_len.min_length,
    facies_result_len.max_length - facies_result_len.mean_length;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
scatter!(ax2, facies_result_len.facies_code, facies_result_len.k_result, color = :red, markersize = 8, label = "k")
lines!(ax2, facies_result_len.facies_code, facies_result_len.k_result, color = :red, linewidth = 1.5)
resize_to_layout!(f)
f
save("plots/length_to_reduce_v2.png", f)

f = Figure()
ax = Axis(f[1, 1],
    xlabel = "Facies",
    ylabel = "length to reduce 50 mg L⁻¹ NO₃⁻ [m]",
    xticks = (1:8, label_values),
    # yticks = 1e-1:2e-1:1.2,
    xticklabelrotation = π/6,
    xgridvisible = false,
    ygridvisible = false,
    # title = "NO₃⁻ Reduction Lengths Across Facies Types",
    # titlefont = "Avenir Book",
    titlesize = 21,
    xlabelsize = 18,
    ylabelsize = 16,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    ylabelfont = "Avenir Book",
    yscale = Makie.pseudolog10,
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(5)
    )
ax2 = Axis(f[1, 1],
    ylabel = "k [m/s]",
    yaxisposition = :right,
    yticks = [1e-7,5.0e-6,1e-5,1.5e-5],
    xgridvisible = false,
    ygridvisible = false,
    xticksvisible = false,
    xlabelsize = 0,
    ylabelsize = 18,
    xticklabelsize = 0,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
)
#ylims!(ax, -1e-9, maximum(facies_result_len.max_length) + 1e-9)
ylims!(ax2, -1e-6, 1.6e-5)
hidespines!(ax2, :t, :l, :b)
linkxaxes!(ax, ax2)
#ax2.yreversed = true
#hidedecorations!(ax2)
hidespines!(ax, :t, :r)
barplot!(ax, facies_result_model.facies_code, facies_result_model.mean_length, color = :steelblue)
errorbars!(ax, facies_result_model.facies_code, facies_result_model.mean_length,
    facies_result_model.mean_length-facies_result_model.min_length,
    facies_result_model.max_length - facies_result_model.mean_length;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
scatter!(ax2, facies_result_model.facies_code, facies_result_model.k_result, color = :red, markersize = 8, label = "k")
lines!(ax2, facies_result_model.facies_code, facies_result_model.k_result, color = :red, linewidth = 1.5)
resize_to_layout!(f)
f
save("plots/length_to_reduce_model_v2.png", f)
println("DONE!")





f = Figure()
ax = Axis(f[2, 1],
    xlabel = "Facies",
    ylabel = "length to reduce 50 mg L⁻¹ NO₃⁻ [m]",
    xticks = (1:8, label_values),
    yticks = LogTicks(-2:1),
    # yticks = 1e-1:2e-1:1.2,
    xticklabelrotation = π/6,
    xgridvisible = false,
    ygridvisible = false,
    # title = "NO₃⁻ Reduction Lengths Across Facies Types",
    # titlefont = "Avenir Book",
    titlesize = 21,
    xlabelsize = 18,
    ylabelsize = 14,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    ylabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    yscale = log10,
    yminorticksvisible = true,
    yminorticks = IntervalsBetween(5)
    )
ax2 = Axis(f[1, 1],
    ylabel = "K [m/s]",
    yaxisposition = :right,
    yticks = [1e-7,5.0e-6,1e-5,1.5e-5],
    xgridvisible = false,
    ygridvisible = false,
    xticksvisible = false,
    xlabelsize = 0,
    ylabelsize = 14,
    xticklabelsize = 0,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    ylabelfont = "Avenir Book",
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
)
#ylims!(ax, -1e-9, maximum(facies_result_len.max_length) + 1e-9)
ylims!(ax2, -1e-6, 1.6e-5)
hidespines!(ax2, :l, :b)
linkxaxes!(ax, ax2)
ax2.yreversed = true
#hidedecorations!(ax2)
hidespines!(ax, :t, :r)
barplot!(ax, facies_result_model.facies_code, facies_result_model.mean_length, color = :steelblue)
errorbars!(ax, facies_result_model.facies_code, facies_result_model.mean_length,
    facies_result_model.mean_length-facies_result_model.min_length,
    facies_result_model.max_length - facies_result_model.mean_length;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
barplot!(ax2, facies_result_model.facies_code, facies_result_model.k_result, color = :darkgrey)
resize_to_layout!(f)
f
save("plots/length_to_reduce_model_v3.png", f)