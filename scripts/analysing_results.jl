using Tidier, DataFrames, CSV, Statistics
# using CairoMakie
# using AlgebraOfGraphics
# load integration results
int_df = CSV.read(datadir("exp_pro","integration_results.csv"), DataFrame)

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
# color = :black, linewidth = 0.6) +
# labs(y = "late times NO3⁻ reduction rate (mol L⁻¹ g⁻¹)", x = "Facies")+
# scale_y_continuous(labels = label_scientific.(0:2e-7:1e-6)) +
# theme()

unique_facies = sort(unique(int_df.facies))
facies_code = Dict(zip(unique_facies, 1:length(unique_facies)))
facies_result.facies_code = map(facies->facies_code[facies],facies_result.facies)

# Labeler = label_scientific()
f = Figure()
ax = Axis(f[1, 1],
    xlabel = "Facies",
    ylabel = "late times NO3⁻ reduction rate (mol L⁻¹ g⁻¹)",
    xticks = (1:9, unique_facies),
    yticks = 1e-7:2e-7:1.2e-6,
    xgridvisible = false,
    ygridvisible = false,
    )
hidespines!(ax, :t, :r)
barplot!(ax, facies_result.facies_code, facies_result.mean_k_no3, color = :steelblue)
errorbars!(ax, facies_result.facies_code, facies_result.mean_k_no3,
    facies_result.mean_k_no3-facies_result.min_k_no3,
    facies_result.max_k_no3 - facies_result.mean_k_no3;
    color = :black, linewidth = 0.8, whiskerwidth = 12)
f
save(plotsdir("facies_k_no3.png"), f)