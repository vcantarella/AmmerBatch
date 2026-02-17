#=
This script is designed to generate a specific, illustrative plot for a single sample ("A308")
to be used as a background for a figure in the manuscript. The plot explains the concept
of zero-order late-time behavior in the batch experiments.

The script performs an analysis similar to `fit_linear_regression_v2.jl`, but is
hardcoded to process only the data for sample "A308".

The main steps are:
1.  **Load Data**:
    - It loads the sample information and the processed data for sample "A308" from the
      `data/exp_pro` directory.
2.  **Data Correction**:
    - It corrects the concentration data for the volume of water removed during sampling,
      just as in the main regression script.
3.  **Identify Zero-Order Phase and Fit Model**:
    - It identifies the time at which the reaction becomes zero-order by analyzing the DOC
      concentration.
    - It performs a linear regression on the nitrate concentration data for the zero-order
      phase to determine the reaction rate and model parameters.
4.  **Generate Plot**:
    - It creates a single, high-quality plot that visualizes:
        - The raw data for nitrate and DOC concentrations.
        - The fitted zero-order model for nitrate reduction.
        - The average DOC concentration during the zero-order phase.
        - A shaded area representing the difference between the data and the model,
          illustrating the initial, non-zero-order phase.
5.  **Save Plot**:
    - The resulting plot is saved in multiple high-quality formats (SVG, PDF, and PNG) in
      the `plots` directory. The file is named `zero_order_explanation`.

This script provides a clear visual explanation of the data analysis method used in the study,
making it a valuable component for the publication.
=#

using DrWatson
@quickactivate "AmmerBatch"
using Roots
using DataFrames, CSV
using Statistics
using CairoMakie
using DataInterpolations
using QuadGK
colors = [:steelblue, :darkseagreen, :crimson]
meas_names = ["NO₃⁻", "DOC", "SO₄²⁻"]

# Part 1:
Vw = 0.08 # L: initial water volume in the batch
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
df_info_p2 = CSV.read(datadir("exp_pro","sample_info_part2.csv"), DataFrame)
# join the two dataframes
df_info = vcat(df_info, df_info_p2)
late_times_params = DataFrame(sample = String[], facies = String[], r_no3 = Float64[], c₀ = Float64[], R² = Float64[], t_zeroorder = Float64[],
                              integral_no3 = Float64[], c_quick = Float64[], t_quick = Float64[],
                              weight = Float64[], TOC = Union{Missing, Float64}[], S_tot = Union{Missing, Float64}[],
                              c₀so4 = Union{Missing, Float64}[], r_so4 = Union{Missing, Float64}[],
                              R²so4 = Union{Missing, Float64}[])
samples = df_info[!, :Sample]

index = findfirst(df_info[!, :Sample] .== "A308")
sample = samples[index]
mass = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1] # g
mass = mass * 1e-3 # convert to kg
df = CSV.read(datadir("exp_pro","$sample.csv"), DataFrame)
# t_change_i = t_change[i]
all_times = df[!, :t]
modelled_no3 = df[!, "NO3-"]
index_no3 = convert.(Bool, 1 .-(ismissing.(modelled_no3)))
index_no3 = findall(index_no3) # indices of non-missing NO₃⁻ values
#modelled_no3 = modelled_no3[index_no3].*1e-3 # convert to mol L⁻¹
modelled_times = all_times[index_no3]
data_fit = LinearInterpolation(modelled_no3[index_no3].*1e-3, modelled_times; extrapolation = ExtrapolationType.Constant)

doc = df[!, "DOC"]
index_doc = convert.(Bool, 1 .-(ismissing.(doc)))
index_doc = findall(index_doc) # indices of non-missing DOC values
#doc = doc[index_doc].*1e-3 # convert to mol L⁻¹
doc_times = all_times[index_doc]
doc_fit = LinearInterpolation(doc[index_doc].*1e-3, doc_times; extrapolation = ExtrapolationType.Constant)
# correct for the reduction in water volume due to sample extraction
vws = Vw * ones(size(all_times))
volume_removed = 0.005 * ones(size(all_times)) # L: volume removed for analysis (4 ml)
mass_removed = data_fit.(all_times) .* volume_removed # mols: mass of NO₃⁻ removed at each time point
doc_removed = doc_fit.(all_times) .* volume_removed # mols: mass of DOC removed at each time point
cum_vol_removed = cumsum(volume_removed)
cum_mass_removed = cumsum(mass_removed) # mols: cumulative mass of NO₃⁻ removed
cum_doc_removed = cumsum(doc_removed) # mols: cumulative mass of DOC removed
vws = Vw .- cum_vol_removed
no3_mol_kg = zeros(length(modelled_times)) # preallocate for NO₃⁻ mol kg⁻¹
global k = 1
for ind in index_no3
    if ind == 1
        no3_mol_kg[1] = modelled_no3[ind]*1e-3 .* 0.08 ./ mass # convert to mol kg⁻¹
    else
        no3_mol_kg[k] = (modelled_no3[ind]*1e-3 .* vws[ind-1] + cum_mass_removed[ind-1]) ./ mass # convert to mol kg⁻¹
    end
    global k += 1
end
doc_mol_kg = zeros(length(doc_times)) # preallocate for DOC mol kg⁻¹
global k = 1
for ind in index_doc
    if ind == 1
        doc_mol_kg[1] = doc[ind]*1e-3 .* 0.08 ./ mass # convert to mol kg⁻¹
    else
        doc_mol_kg[k] = (doc[ind]*1e-3 .* vws[ind-1] + cum_doc_removed[ind-1]) ./ mass # convert to mol kg⁻¹
    end
    global k += 1
end
no3_fit = LinearInterpolation(no3_mol_kg, modelled_times; extrapolation = ExtrapolationType.Constant)
doc_fit = LinearInterpolation(doc_mol_kg, doc_times; extrapolation = ExtrapolationType.Constant)
# find the last two DOC values and calculate the average
doc_last_two = doc_mol_kg[end-2:end]
doc_avg = mean(doc_last_two)
# find the time when the DOC concentration crosses the average
t_change_i = findfirst(doc_fit.(all_times) .<= doc_avg)
if t_change_i === nothing
    t_change_i = length(all_times) # if it never crosses, then we take the last time point
end
t_change = all_times[t_change_i] # convert to time
println("Sample: $(sample), t_change: $(t_change)")
# calculate the slope of the NO₃⁻ reduction
# Linear regression:
X = [ones(size(modelled_times[modelled_times.>t_change])) modelled_times[modelled_times.>t_change]]
y = no3_mol_kg[modelled_times.>t_change]
# linear regression and model fit results:
β = (X'X)\(X'y)
ϵ = y - X*β
s² = ϵ'ϵ/(length(y)-length(β))
σ² = s²*(length(y)-length(β))/length(y)
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSR = sum(ϵ.^2)
R² = 1 - SSR/SST
r_no3 = -β[2] # slope of the NO₃⁻ reduction
c0_model = β[1] # initial concentration of NO₃⁻
model_no3(t) = c0_model - r_no3*t
# calculate when data intercepts the model:
println("Sample: $(sample)")

t_of_zero = t_change > 0 ? find_zeros(t-> no3_fit(t)-model_no3(t), 0, 150)[1] : 0.0
# calculate and plot the integral using Trapz
integral_no3 = quadgk(no3_fit, 0.0, t_of_zero)[1]
c0_data = no3_fit(0.0) # initial concentration
c_quick = t_change > 0 ? c0_data - c0_model : 0.0

integral_function = quadgk(model_no3, 0.0, t_of_zero)[1]
# trapz(times_n, model_no3.(times_n))
resulting_integral = (integral_no3 - integral_function)

t_quick = t_change > 0 ? resulting_integral/c_quick : 0.0
fontcolor = "#474747"
inch = 96
pt = 4/3
cm = inch / 2.54
fig = Figure(resolution = (1500, 900), fontsize = 16)
ax = Axis(fig[1,1],
    xlabel = "Time (days)",
    ylabel = "Concentration (mol kg⁻¹)",
    yticks = 0:0.03:0.1,
    width = 500,
    height = 300,
    xgridvisible = false,
    ygridvisible = false,
    titlesize = 20,
    titlealign = :center,
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
hidespines!(ax, :t, :r)
ylims!(ax, (0, 3.8e-3*Vw/mass)) # mol kg⁻¹
xlims!(ax,(-2, 180))
if t_change > 0
        tfill = 0:0.1:t_of_zero
        fill_between!(ax, tfill, model_no3.(tfill),
            no3_fit.(tfill), color = colors[1], alpha = 0.5)
end
lines!(ax, all_times, model_no3.(all_times), color = colors[1], label = "NO₃⁻",
    linewidth = 2.4)
avg_doc = mean(doc_mol_kg[modelled_times.>t_change])
lines!(ax, all_times, fill(avg_doc, length(all_times)), color = colors[2],
    label = "DOC", linestyle = :dash,
    linewidth = 2.4)
scatter!(ax, modelled_times, no3_mol_kg, color = colors[1], label = meas_names[1],
    marker = :diamond, markersize = 10)
scatter!(ax, doc_times, doc_mol_kg, color = colors[2], label = meas_names[2],
    marker = :diamond, markersize = 10)
axislegend(ax, framevisible = false, position = :rt, orientation = :vertical,
    titlefont = "Avenir Book", titlesize = 18, titlecolor = :black,
    merge = true,
    labelsize = 16, labelcolor = :black, backgroundcolor = :transparent)
resize_to_layout!(fig)
#Legend(grid_plot_figs[grid_plot_num][5,1:2], ax, "Substance", merge = true, framevisible = false, orientation = :horizontal)

# Save in multiple high-quality formats
save(plotsdir("zero_order_explanation.svg"), fig)  # Vector format for scalability
save(plotsdir("zero_order_explanation.pdf"), fig)  # High-quality PDF
save(plotsdir("zero_order_explanation.png"), fig, px_per_unit = 3)  # High DPI PNG

fig