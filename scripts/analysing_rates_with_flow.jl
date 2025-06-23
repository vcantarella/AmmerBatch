using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = DataFrame(XLSX.readtable(datadir("exp_pro","integration_results.xlsx"), "Sheet1"))
int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"

int_df[!, :r_no3_weight] = int_df[!, :r_no3] ./ int_df[!, :weight]

# group by facies

facies_result = @chain(
    int_df,
    @group_by(facies),
    @summarize(
        mean_r_no3 = mean(r_no3_weight),
        median_r_no3 = median(r_no3_weight),
        std_r_no3 = std(r_no3_weight),
        min_r_no3 = minimum(r_no3_weight),
        max_r_no3 = maximum(r_no3_weight),
    )
)

# Read the velocity data
vel_df = CSV.read(datadir("external","velocity_facies.csv"), DataFrame)

# left join the velocity data with the facies result
facies_result = leftjoin(facies_result, vel_df, on = :facies)

# calculate the length
ρᵦ = 2.65 # g/cm³ or kg/L
# liquid rate
@. facies_result[!, :r_no3_l] = facies_result[!, :mean_r_no3] * ρᵦ * (1 - facies_result[!, :porosity])/ facies_result[!, :porosity] # mol L⁻¹ s⁻¹

# length for denitrifying a given concentration
c_no3 = 50e-3 / 62.0049 # mol L⁻¹
@. facies_result[!, :length] = c_no3 / facies_result[!, :r_no3_l] * facies_result[!, :mean] * 86400 # conversion from m/s to m/d # m


unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result.facies_code = map(facies->facies_code[facies],facies_result.facies)
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
barplot!(ax, facies_result.facies_code, facies_result.length, color = :steelblue)
resize_to_layout!(f)
f
save("plots/length_to_reduce_50mgL-1_NO3-.png", f)