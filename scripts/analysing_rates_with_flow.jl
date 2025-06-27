using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = DataFrame(CSV.File(datadir("exp_pro","integration_results_final.csv")))

# Read the velocity data
vel_df = CSV.read(datadir("external","velocity_facies_v2.csv"), DataFrame)


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

# length for denitrifying a given concentration
c_no3 = 50e-3 / 62.0049 # mol L⁻¹ - given concentration inflow of NO₃⁻
c_no3 = c_no3 * 1e3 # mol m⁻³
length_to_reduce = c_no3 * q ./ (ρᵦ .* int_df[!, :r_no3_weight]) # m
int_df[!, :length] = length_to_reduce


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


unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result_len.facies_code = map(facies->facies_code[facies],facies_result_len.facies)
label_values = ["Clay", "Tufa grains", "Calcareous silt", "Tufa & reed", "Silt & moss", "Silt & organic debris",
    "Brown peat", "Black peat"]
labels = Dict(zip(unique_facies, label_values))
fontcolor = "#474747"
# Labeler = label_scientific()
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
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    )
hidespines!(ax, :t, :r)
barplot!(ax, facies_result_len.facies_code, facies_result_len.mean_length, color = :steelblue)
errorbars!(ax, facies_result_len.facies_code, facies_result_len.mean_length,
    facies_result_len.mean_length-facies_result_len.min_length,
    facies_result_len.max_length - facies_result_len.mean_length;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
resize_to_layout!(f)
f
save("plots/length_to_reduce.png", f)