using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = DataFrame(XLSX.readtable(datadir("exp_pro","integration_results.xlsx"), "Sheet1"))
int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"

int_df[!, :k_no3_weight] = int_df[!, :k_no3] ./ int_df[!, :weight]

# group by facies

facies_result = @chain(
    int_df,
    @group_by(facies),
    @summarize(
        mean_k_no3 = mean(k_no3_weight),
        std_k_no3 = std(k_no3_weight),
        min_k_no3 = minimum(k_no3_weight),
        max_k_no3 = maximum(k_no3_weight),
    )
)

# ggplot(facies_result, @aes(x = facies, y=mean_k_no3, fill = facies)) + geom_col() +
# theme_minimal() + 
# geom_errorbar(@aes(ymin = min_k_no3, ymax = max_k_no3),
# color = fontcolor, linewidth = 0.6) +
# labs(y = "late times NO3⁻ reduction rate (mol L⁻¹ g⁻¹)", x = "Facies")+
# scale_y_continuous(labels = label_scientific.(0:2e-7:1e-6)) +
# theme()

unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result.facies_code = map(facies->facies_code[facies],facies_result.facies)
label_values = ["Clay", "Tufa grains", "Calcareous silt", "Tufa & reed", "Silt & moss", "Silt & organic debris",
    "Brown peat", "Black peat"]
labels = Dict(zip(unique_facies, label_values))

fontcolor = "#474747"
# Labeler = label_scientific()
f = Figure(backgroundcolor = :transparent)
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-6}", fontsize = 16,
    color = fontcolor, font = "Avenir Book")
ax = Axis(f[1, 1],
    #xlabel = "Facies",
    ylabel = "r₀ [mol L⁻¹ g⁻¹ d⁻¹]",
    xticks = (1:8, label_values),
    yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    title = "NO₃⁻ Reduction Rates Across Facies Types",
    titlefont = "Avenir Book",
    titlesize = 20,
    titlealign = :center,
    titlecolor = fontcolor,
    #xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 14,
    xticklabelrotation = deg2rad(35),
    #xlabelrotation = 45,
    yticklabelsize = 16,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    #xlabelcolor = fontcolor,
    ylabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    #xlabelfont = "Avenir Book",
    ylabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    )
hidespines!(ax, :t, :r)
barplot!(ax, facies_result.facies_code, facies_result.mean_k_no3.*1e6, color = :steelblue)
errorbars!(ax, facies_result.facies_code, facies_result.mean_k_no3.*1e6,
    facies_result.mean_k_no3.*1e6-facies_result.min_k_no3.*1e6,
    facies_result.max_k_no3.*1e6 - facies_result.mean_k_no3.*1e6;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
resize_to_layout!(f)
f
save(plotsdir("facies_k_no3.png"), f)
save(plotsdir("facies_k_no3.svg"), f)

f = Figure(backgroundcolor = :transparent)
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-6}", fontsize = 16, 
    color = fontcolor, font = "Avenir Book")
ax = Axis(f[1, 1],
    xlabel = "TOC [%]",
    ylabel = "r₀ [mol L⁻¹ g⁻¹ d⁻¹]",
    # title = "NO₃⁻ Reduction Rates Correlate with TOC",
    # xticks = (1:9, unique_facies),
    yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    # titlesize = 20,
    # titlealign = :right,
    # titlecolor = fontcolor,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 16,
    yticklabelsize = 16,
    xlabelcolor = fontcolor,
    ylabelcolor = fontcolor,
    xticklabelcolor = fontcolor,
    yticklabelcolor = fontcolor,
    xticklabelfont = "Avenir Book",
    yticklabelfont = "Avenir Book",
    xlabelfont = "Avenir Book",
    ylabelfont = "Avenir Book",
    backgroundcolor = :transparent,
    )
hidespines!(ax, :t, :r)
# Plot each facies with a different color

for facies in unique_facies
    facies_mask = int_df.facies .== facies
    scatter!(ax, 
             int_df[facies_mask, "TOC"], 
             int_df[facies_mask, :k_no3_weight] .* 1e6, 
             label = labels[facies],
             markersize = 18)
end

# Add a legend
Legend(f[1,2], ax, "Facies", framevisible = false, position = :lt, orientation = :vertical,
    titlefont = "Avenir Book", titlesize = 18, titlecolor = fontcolor,
    labelsize = 16, labelcolor = fontcolor, backgroundcolor = :transparent)
k_no3 = int_df[!, :k_no3_weight]
TOC = int_df[!, "TOC"]
TOC = convert.(Float64, TOC)
X = [ones(size(TOC)) TOC]
y = k_no3
β = (X'X)\(X'y)
ϵ = y - X*β
s² = ϵ'ϵ/(length(y)-length(β))
σ² = s²*(length(y)-length(β))/length(y)
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSR = sum(ϵ.^2)
R² = 1 - SSR/SST
lines!(ax, 0:0.1:maximum(TOC), (β[1] .+ β[2].*(0:0.1:maximum(TOC))).*1e6, color = :crimson, linewidth = 2.8)
text!(ax, 0.5, 0.9,
    text="R² = $(round(R², digits=2))",
    color = fontcolor, space = :relative)
Label(f[1, 1:2, Top()], halign = :center, "NO₃⁻ Reduction Rates Correlate with TOC",
fontsize = 20,
    color = fontcolor, font = "Avenir Book")
resize_to_layout!(f)
f
save(plotsdir("kno3_toc.svg"), f)
save(plotsdir("kno3_toc.png"), f)
