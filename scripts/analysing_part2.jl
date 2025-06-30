using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
using LaTeXStrings

# I want to analyze the rates of nitrate reduction and sulfur oxidation and their correlation

Vw = 80*1e-3 # ml -> L
int_df = DataFrame(XLSX.readtable(datadir("exp_pro","integration_results_part2.xlsx"), "Sheet1"))
int_df[!, :r_no3_weight] = int_df[!, :r_no3] ./ int_df[!, :weight]
int_df[!, :r_so4_no3] = int_df[!, :r_so4] ./ int_df[!, :r_no3]
int_df[!, :r_so4_weight] = int_df[!, :r_so4] ./ int_df[!, :weight]# mol L⁻¹ g⁻¹ d⁻¹
unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))

int_df.facies_code = map(facies->facies_code[facies],int_df.facies)


fontcolor = "#474747"
r_no3 = int_df[!, :r_no3]
r_no3 = convert.(Float64, r_no3)
r_so4_no3 = int_df[!, :r_so4_no3]
r_so4_no3 = convert.(Float64, r_so4_no3)
X = [ones(size(r_no3)) r_no3]
y = r_so4_no3
β = (X'X)\(X'y)
ϵ = y - X*β
s² = ϵ'ϵ/(length(y)-length(β))
σ² = s²*(length(y)-length(β))/length(y)
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSR = sum(ϵ.^2)
R² = 1 - SSR/SST

r_so4_no3
min_s04_part = r_so4_no3/1.6
int_df.min_s04_part = min_s04_part
max_s04_part = r_so4_no3*1.2
int_df.max_s04_part = max_s04_part

f = Figure()#backgroundcolor = :transparent)
ax = Axis(f[1, 1],
    xlabel = "max_s04_part [-]",
    ylabel = "r₀ [mol L⁻¹ g⁻¹ d⁻¹]",
    title = "-",
    # xticks = (1:9, unique_facies),
    #yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    titlesize = 16,
    titlealign = :right,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    #backgroundcolor = :transparent,
    )
hidespines!(ax, :t, :r)
# Plot each facies with a different color
label_values = ["Clay", "Tufa grains", "Calcareous silt", "Tufa & reed", "Silt & moss", "Silt & organic debris",
    "Brown peat", "Black peat"]
labels = Dict(zip(unique_facies, label_values))
for facies in unique_facies
    facies_mask = int_df.facies .== facies
    scatter!(ax, 
             int_df[facies_mask, :max_s04_part], 
             int_df[facies_mask, :r_no3_weight], 
             label = labels[facies],
             markersize = 18)
end
f
# Add a legend
# Legend(f[1,2], ax, "Facies", framevisible = false, position = :lt, orientation = :vertical,
#     titlefont = "Avenir Book", titlesize = 18, titlecolor = fontcolor,
#     labelsize = 16, labelcolor = fontcolor, backgroundcolor = :transparent)
r_no3 = int_df[!, :r_no3]
r_no3 = convert.(Float64, r_no3)
max_s04_part = int_df[!, :max_s04_part]
max_s04_part = convert.(Float64, max_s04_part)
X = [ones(size(max_s04_part)) max_s04_part]
y = r_no3
β = (X'X)\(X'y)
ϵ = y - X*β
s² = ϵ'ϵ/(length(y)-length(β))
σ² = s²*(length(y)-length(β))/length(y)
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSR = sum(ϵ.^2)
R² = 1 - SSR/SST
#lines!(ax, 0:0.01:(maximum(r_no3)*1e6), (β[1] .+ β[2].*(0:0.01:(maximum(r_no3)*1e6))), color = :crimson, linewidth = 2.8)
# text!(ax, 0.5, 0.9,
#     text="R² = $(round(R², digits=2))",
#     color = fontcolor, space = :relative)
resize_to_layout!(f)
f

using StatsBase, GLM
int_df[!,:r_no3] = convert.(Float64, int_df[!,:r_no3])
int_df[!,:r_so4] = convert.(Float64, int_df[!,:r_so4])

lm1 = lm(@formula(r_no3 ~ r_so4), int_df)
display(lm1)

# stoichiometric ratio of the reaction:

f = Figure()#backgroundcolor = :transparent)
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-6}", fontsize = 16)
Label(f[1, 1, BottomRight()], halign = :left, L"\times 10^{-6}", fontsize = 16)
ax = Axis(f[1, 1],
    ylabel = L"$r_{NO_3^-}$ [mol L⁻¹ d⁻¹]",
    xlabel = L"$r_{SO_4^{2-}}$ [mol L⁻¹ d⁻¹]",
    title = "Comparing reaction rates of NO₃⁻ and SO₄²⁻ \n with the stoichiometric ratio",
    # xticks = (1:9, unique_facies),
    #yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    titlesize = 19,
    titlealign = :right,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    #backgroundcolor = :transparent,
    )
hidespines!(ax, :t, :r)
# Plot each facies with a different color
label_values = ["C1: Clay", "T1: Tufa grains", "T2: Calcareous silt", "T4: Tufa & reed", "T6: Silt & moss", "T7: Silt & organic debris",
    "T8: Brown peat", "T9: Black peat"]
labels = Dict(zip(unique_facies, label_values))
for facies in unique_facies
    facies_mask = int_df.facies .== facies
    scatter!(ax, 
             int_df[facies_mask, :r_so4_weight] .* 1e6, 
             int_df[facies_mask, :r_no3_weight] .* 1e6, 
             label = labels[facies],
             markersize = 18)
end
# plot the theoretical lines for the stoichiometric ratio
stoich_ratio = 1.6 # FeS
lines!(ax, 0:0.1:maximum(int_df[!, :r_so4_weight]) .* 1e6, stoich_ratio .* (0:0.1:maximum(int_df[!, :r_so4_weight])* 1e6) , color = :crimson, linewidth = 2.8,
    label = "Stoichiometric Ratio: 8/5")
# label the line
stoich_ratio = 1.4 # FeS2
lines!(ax, 0:0.1:maximum(int_df[!, :r_so4_weight]) .* 1e6, stoich_ratio .* (0:0.1:maximum(int_df[!, :r_so4_weight])* 1e6) , color = :steelblue, linewidth = 2.8,
    label = "Stoichiometric Ratio: 7/5")
stoich_ratio = 1.2 # S0
lines!(ax, 0:0.1:maximum(int_df[!, :r_so4_weight]) .* 1e6, stoich_ratio .* (0:0.1:maximum(int_df[!, :r_so4_weight])* 1e6) , color = :darkgreen, linewidth = 2.8,
    label = "Stoichiometric Ratio: 6/5")
axislegend(ax, position = :lt, merge = true)
resize_to_layout!(f)
f
save(plotsdir("stoichiometric_ratio_no3_so4.png"), f)
save(plotsdir("stoichiometric_ratio_no3_so4.svg"), f)