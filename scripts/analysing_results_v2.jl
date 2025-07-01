using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = CSV.read(datadir("exp_pro","linear_regression_params_v2.csv"), DataFrame)
int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"

# group by facies
facies_result = @chain(
    int_df,
    @group_by(facies),
    @summarize(
        mean_r_no3 = mean(r_no3),
        std_r_no3 = std(r_no3),
        min_r_no3 = minimum(r_no3),
        max_r_no3 = maximum(r_no3),
    )
)
colors_map = :Paired_8

unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result.facies_code = map(facies->facies_code[facies],facies_result.facies)
label_values = ["Clay", "Tufa grains", "Calcareous silt", "Tufa & reed", "Silt & moss", "Silt & organic debris",
    "Brown peat", "Black peat"]
labels = Dict(zip(unique_facies, label_values))

fontcolor = "#474747"
# Labeler = label_scientific()
f = Figure()
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-4}", fontsize = 16,
    color = fontcolor, font = "Avenir Book", padding = (-60, 0, 0, 0))
ax = Axis(f[1, 1],
    #xlabel = "Facies",
    width = 400,
    height = 300,
    ylabel = "r₀ [mol kg⁻¹ d⁻¹]",
    xticks = (1:8, label_values),
    #yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    title = "a.",
    titlefont = "Avenir Book bold",
    titlesize = 20,
    titlealign = :left,
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
barplot!(ax, facies_result.facies_code, facies_result.mean_r_no3.*1e4, color = :steelblue)
errorbars!(ax, facies_result.facies_code, facies_result.mean_r_no3.*1e4,
    facies_result.mean_r_no3.*1e4-facies_result.min_r_no3.*1e4,
    facies_result.max_r_no3.*1e4 - facies_result.mean_r_no3.*1e4;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
Label(f[1, 2, Top()], halign = :left, L"\times 10^{-4}", fontsize = 16, 
    color = fontcolor, font = "Avenir Book", padding = (-70, 0, 0, 0))
ax2 = Axis(f[1, 2],
    width = 400,
    height = 300,
    xlabel = "TOC [%]",
    ylabel = "r₀ [mol kg⁻¹ d⁻¹]",
    title = "b.",
    # xticks = (1:9, unique_facies),
    # yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    titlesize = 20,
    titlealign = :left,
    titlecolor = fontcolor,
    titlefont = "Avenir Book bold",
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
hidespines!(ax2, :t, :r)
# Plot each facies with a different color
c = cgrad(colors_map, (1:length(unique_facies))./length(unique_facies), categorical=true)
k = 1
for facies in unique_facies
    global k
    facies_mask = int_df.facies .== facies
    # filter missing TOC
    println(facies)
    local real_idx = findall(!ismissing, int_df[facies_mask, "TOC"])
    scatter!(ax2, 
             convert.(Float64, int_df[facies_mask, "TOC"][real_idx]), 
             int_df[facies_mask, :r_no3][real_idx] .* 1e4, 
             label = labels[facies],
             color = c[k],
             markersize = 18)
    k += 1
end
linkyaxes!([ax, ax2])

# Add a legend
Legend(f[1,3], ax2, "Facies", framevisible = false, position = :lt, orientation = :vertical,
    titlefont = "Avenir Book", titlesize = 18, titlecolor = fontcolor,
    labelsize = 16, labelcolor = fontcolor, backgroundcolor = :transparent)
r_no3 = int_df[!, :r_no3]


TOC = int_df[!, "TOC"]
real_idx = findall(!ismissing, TOC)
TOC = convert.(Float64, TOC[real_idx])
X = [ones(size(TOC)) TOC]
y = r_no3[real_idx]
β = (X'X)\(X'y)
ϵ = y - X*β
s² = ϵ'ϵ/(length(y)-length(β))
σ² = s²*(length(y)-length(β))/length(y)
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSR = sum(ϵ.^2)
R² = 1 - SSR/SST
lines!(ax2, 0:0.1:maximum(TOC), (β[1] .+ β[2].*(0:0.1:maximum(TOC))).*1e4, color = :crimson, linewidth = 2.8)
text!(ax2, 0.5, 0.9,
    text="R² = $(round(R², digits=2))",
    color = fontcolor, space = :relative)

resize_to_layout!(f)
save(plotsdir("rno3_facies_toc_v2.png"), f)
save(plotsdir("rno3_facies_toc_v2.svg"), f)
save(plotsdir("rno3_facies_toc_v2.pdf"), f)
# Save the dataframe with the results
CSV.write(datadir("exp_pro", "facies_results.csv"), facies_result)