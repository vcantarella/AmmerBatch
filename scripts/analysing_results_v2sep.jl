using DrWatson
@quickactivate "AmmerBatch"
using Tidier, DataFrames, CSV, Statistics
using CairoMakie
using XLSX
# using CairoMakie
# load integration results
int_df = CSV.read(datadir("exp_pro","linear_regression_params_v2.csv"), DataFrame)
# Remove sample B7 from the dataframe: results have been discarded after intrinsic analysis
int_df = int_df[int_df[!, :sample] .!= "B7", :]

int_df[(int_df[!, :facies].=="C1").||(int_df[!, :facies].=="C2"), :facies] .= "C1"
int_df[!, :exp] = [ifelse(int_df[i, :sample][1] == 'A', 'A', 'B') for i in 1:nrow(int_df)]
# group by facies
facies_result = @chain(
    int_df,
    @group_by(exp, facies),
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
for df in facies_result
    df.facies_code = map(facies->facies_code[facies], df.facies)
end
label_values = ["Clay (FT-C1)", "Tufa grains (FT-01)", "Calcareous silt (FT-02)", "Tufa & reed (FT-04)", "Silt & moss (FT-06)", "Silt & organic debris (FT-07)",
    "Brown peat (FT-08)", "Black peat (FT-09)"]
labels = Dict(zip(unique_facies, label_values))

fontcolor = "#474747"
# Labeler = label_scientific()
f = Figure()
Label(f[2, 1, Top()], halign = :left, L"\times 10^{-4}", fontsize = 18,
    color = fontcolor, font = "Avenir Book", padding = (-60, 0, 0, 0))
axa = Axis(f[2, 1],
    #xlabel = "Facies",
    width = 400,
    height = 300,
    ylabel = L"$r_0$ [mol kg^{-1} d^{-1}]",
    xticks = (1:8, label_values),
    #yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    title = "c.",# Experiment B - facies",
    titlefont = "Avenir Book bold",
    titlesize = 20,
    titlealign = :left,
    titlecolor = fontcolor,
    #xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 14,
    xticklabelrotation = deg2rad(30),
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
# ylims!
hidespines!(axa, :t, :r)
barplot!(axa, facies_result[1].facies_code, facies_result[1].mean_r_no3.*1e4, color = :steelblue)
errorbars!(axa, facies_result[1].facies_code, facies_result[1].mean_r_no3.*1e4,
    facies_result[1].mean_r_no3.*1e4-facies_result[1].min_r_no3.*1e4,
    facies_result[1].max_r_no3.*1e4 - facies_result[1].mean_r_no3.*1e4;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
Label(f[1, 1, Top()], halign = :left, L"\times 10^{-4}", fontsize = 18,
    color = fontcolor, font = "Avenir Book", padding = (-60, 0, 0, 0))
axc = Axis(f[1, 1],
    #xlabel = "Facies",
    width = 400,
    height = 300,
    ylabel = L"$r_0$ [mol kg^{-1} d^{-1}]",
    xticks = (1:8, label_values),
    #yticks = 1e-1:2e-1:1.2,
    xgridvisible = false,
    ygridvisible = false,
    title = "a.",# Experiment A - facies",
    titlefont = "Avenir Book bold",
    titlesize = 20,
    titlealign = :left,
    titlecolor = fontcolor,
    #xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 14,
    xticklabelrotation = deg2rad(30),
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
hidespines!(axc, :t, :r)
barplot!(axc, facies_result[2].facies_code, facies_result[2].mean_r_no3.*1e4, color = :steelblue)
errorbars!(axc, facies_result[2].facies_code, facies_result[2].mean_r_no3.*1e4,
    facies_result[2].mean_r_no3.*1e4-facies_result[2].min_r_no3.*1e4,
    facies_result[2].max_r_no3.*1e4 - facies_result[2].mean_r_no3.*1e4;
    color = fontcolor, linewidth = 0.8, whiskerwidth = 12)
Label(f[2, 2, Top()], halign = :left, L"\times 10^{-4}", fontsize = 18, 
    color = fontcolor, font = "Avenir Book", padding = (-70, 0, 0, 0))
axb = Axis(f[2, 2],
    width = 400,
    height = 300,
    xlabel = "TOC [%]",
    ylabel = L"$r_0$ [mol kg^{-1} d^{-1}]",
    title = "d.",# Experiment B - TOC",
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
hidespines!(axb, :t, :r)
# Plot each facies with a different color
c = cgrad(colors_map, (1:length(unique_facies))./length(unique_facies), categorical=true)
k = 1
for facies in unique_facies
    global k
    facies_mask = int_df.facies .== facies
    # filter missing TOC
    println(facies)
    local real_idx = findall(!ismissing, int_df[facies_mask, "TOC"])
    # filter experiment A
    real_idx = real_idx[int_df[facies_mask, :exp][real_idx] .== 'A']
    scatter!(axb, 
             convert.(Float64, int_df[facies_mask, "TOC"][real_idx]), 
             int_df[facies_mask, :r_no3][real_idx] .* 1e4, 
             label = labels[facies],
             color = c[k],
             markersize = 18)
    k += 1
end


Label(f[1, 2, Top()], halign = :left, L"\times 10^{-4}", fontsize = 18, 
    color = fontcolor, font = "Avenir Book", padding = (-70, 0, 0, 0))
axd = Axis(f[1, 2],
    width = 400,
    height = 300,
    xlabel = "TOC [%]",
    ylabel = L"$r_0$ [mol kg^{-1} d^{-1}]",
    title = "b.",# Experiment A - TOC",
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
hidespines!(axd, :t, :r)
# Plot each facies with a different color
k = 1
for facies in unique_facies
    global k
    facies_mask = int_df.facies .== facies
    # filter missing TOC
    println(facies)
    local real_idx = findall(!ismissing, int_df[facies_mask, "TOC"])
    # filter experiment B
    real_idx = real_idx[int_df[facies_mask, :exp][real_idx] .== 'B']
    scatter!(axd, 
             convert.(Float64, int_df[facies_mask, "TOC"][real_idx]), 
             int_df[facies_mask, :r_no3][real_idx] .* 1e4, 
             label = labels[facies],
             color = c[k],
             markersize = 18)
    k += 1
end
linkyaxes!([axc, axd])
linkyaxes!([axa, axb])
# Add a legend
Legend(f[1:2,3], axb, "Facies Types (FT)", framevisible = false, position = :lt, orientation = :vertical,
    titlefont = "Avenir Book", titlesize = 18, titlecolor = fontcolor,
    labelsize = 16, labelcolor = fontcolor, backgroundcolor = :transparent)



text!(axc, 5, 0.1, text="no sample", color = fontcolor, rotation = deg2rad(90),
    fontsize = 16,
    font = "Avenir Book")
text!(axc, 7, 0.1, text= "no sample", color = fontcolor, rotation = deg2rad(90),
     fontsize = 16,
    font = "Avenir Book")



# Fit a linear regression model to the data for experiment A
r_no3 = int_df[int_df[!,:exp].=='A', :r_no3]
TOC = int_df[int_df[!,:exp].=='A', "TOC"]
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
lines!(axb, 0:0.1:maximum(TOC), (β[1] .+ β[2].*(0:0.1:maximum(TOC))).*1e4, color = :crimson, linewidth = 2.8)
text!(axb, 0.5, 0.8,
    text="R² = $(round(R², digits=2))",
    color = fontcolor, space = :relative)
# Fit a linear regression model to the data for experiment B
r_no3 = int_df[int_df[!,:exp].=='B', :r_no3]
# r_no3 = r_no3[setdiff(1:length(r_no3), 7)]
TOC = int_df[int_df[!,:exp].=='B', "TOC"]
# TOC = TOC[setdiff(1:length(TOC), 7)]
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
lines!(axd, 0:0.1:maximum(TOC), (β[1] .+ β[2].*(0:0.1:maximum(TOC))).*1e4, color = :crimson, linewidth = 2.8)
text!(axd, 0.5, 0.8,
    text="R² = $(round(R², digits=2))",
    color = fontcolor, space = :relative)

resize_to_layout!(f)
save(plotsdir("rno3_facies_toc_v2sep.png"), f)
save(plotsdir("rno3_facies_toc_v2sep.svg"), f)
save(plotsdir("rno3_facies_toc_v2sep.pdf"), f)
# Save the dataframe with the results
CSV.write(datadir("exp_pro", "facies_results.csv"), facies_result)