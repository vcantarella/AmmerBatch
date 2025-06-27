using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = DataFrame(XLSX.readtable(datadir("exp_pro","integration_results.xlsx"), "Sheet1"))
int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"
int_df2 = DataFrame(XLSX.readtable(datadir("exp_pro","integration_results_part2.xlsx"), "Sheet1"))
Vw = 80*1e-3 # ml -> L
int_df[!, :r_no3_weight] = int_df[!, :r_no3] ./ (int_df[!, :weight]*1e-3) .* Vw # mol kgs⁻¹ d⁻¹
int_df2[!, :r_no3_weight] = int_df2[!, :r_no3] ./ (int_df2[!, :weight]*1e-3) .* Vw # mol kgs⁻¹ d⁻¹
int_df2 = select(int_df2, Not([:s0, :r_so4, :R²_so4, :r_so4_no3]))
# combine the two dataframes
int_df = vcat(int_df, int_df2)
# group by facies

facies_result = @chain(
    int_df,
    @group_by(facies),
    @summarize(
        mean_r_no3 = mean(r_no3_weight),
        std_r_no3 = std(r_no3_weight),
        min_r_no3 = minimum(r_no3_weight),
        max_r_no3 = maximum(r_no3_weight),
    )
)

# ggplot(facies_result, @aes(x = facies, y=mean_k_no3, fill = facies)) + geom_col() +
# theme_minimal() + 
# geom_errorbar(@aes(ymin = min_k_no3, ymax = max_k_no3),
# color = fontcolor, linewidth = 0.6) +
# labs(y = "late times NO3⁻ reduction rate (mol L⁻¹ g⁻¹)", x = "Facies")+
# scale_y_continuous(labels = label_scientific.(0:2e-7:1e-6)) +
# theme()

colors_map = :Paired_8

unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result.facies_code = map(facies->facies_code[facies],facies_result.facies)
label_values = ["Clay", "Tufa grains", "Calcareous silt", "Tufa & reed", "Silt & moss", "Silt & organic debris",
    "Brown peat", "Black peat"]
labels = Dict(zip(unique_facies, label_values))

fontcolor = "#474747"
# Labeler = label_scientific()
f = Figure(backgroundcolor = :transparent)
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-3}", fontsize = 16,
    color = fontcolor, font = "Avenir Book", padding = (-60, 0, 0, 0))
ax = Axis(f[1, 1],
    #xlabel = "Facies",
    width = 400,
    height = 300,
    ylabel = "r₀ [mol kg⁻¹ d⁻¹]",
    xticks = (1:8, label_values),
    yticks = 1e-1:2e-1:1.2,
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
barplot!(ax, facies_result.facies_code, facies_result.mean_r_no3.*1e3, color = :steelblue)
errorbars!(ax, facies_result.facies_code, facies_result.mean_r_no3.*1e3,
    facies_result.mean_r_no3.*1e3-facies_result.min_r_no3.*1e3,
    facies_result.max_r_no3.*1e3 - facies_result.mean_r_no3.*1e3;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
# resize_to_layout!(f)
# f
# save(plotsdir("facies_k_no3.png"), f)
# save(plotsdir("facies_k_no3.svg"), f)

# f = Figure(backgroundcolor = :transparent)
Label(f[1, 2, Top()], halign = :left, L"\times 10^{-3}", fontsize = 16, 
    color = fontcolor, font = "Avenir Book", padding = (-70, 0, 0, 0))
ax2 = Axis(f[1, 2],
    width = 400,
    height = 300,
    xlabel = "TOC [%]",
    ylabel = "r₀ [mol kg⁻¹ d⁻¹]",
    title = "b.",
    # xticks = (1:9, unique_facies),
    yticks = 1e-1:2e-1:1.2,
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
    real_idx = findall(!ismissing, int_df[facies_mask, "TOC"])
    scatter!(ax2, 
             convert.(Float64, int_df[facies_mask, "TOC"][real_idx]), 
             int_df[facies_mask, :r_no3_weight][real_idx] .* 1e3, 
             label = labels[facies],
             color = c[k],
             markersize = 18)
    k += 1
end

# Add a legend
Legend(f[1,3], ax2, "Facies", framevisible = false, position = :lt, orientation = :vertical,
    titlefont = "Avenir Book", titlesize = 18, titlecolor = fontcolor,
    labelsize = 16, labelcolor = fontcolor, backgroundcolor = :transparent)
r_no3 = int_df[!, :r_no3_weight]


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
lines!(ax2, 0:0.1:maximum(TOC), (β[1] .+ β[2].*(0:0.1:maximum(TOC))).*1e3, color = :crimson, linewidth = 2.8)
text!(ax2, 0.5, 0.9,
    text="R² = $(round(R², digits=2))",
    color = fontcolor, space = :relative)
# Label(f[1, 1:2, Top()], halign = :center, "NO₃⁻ Reduction Rates Correlate with TOC",
# fontsize = 20,
#     color = fontcolor, font = "Avenir Book")
resize_to_layout!(f)
save(plotsdir("rno3_facies_toc.png"), f)
save(plotsdir("rno3_facies_toc.svg"), f)
save(plotsdir("rno3_facies_toc.pdf"), f)
# f
# save(plotsdir("kno3_toc.svg"), f)
# save(plotsdir("kno3_toc.png"), f)

# Save the dataframe with the results
CSV.write(datadir("exp_pro", "facies_results.csv"), facies_result)
# Save the dataframe with the integration results
CSV.write(datadir("exp_pro", "integration_results_final.csv"), int_df)