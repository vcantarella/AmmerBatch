#=This script reads the raw experimental data from Excel files (`.xlsx`) located in `data/exp_raw`.It processes the data from two separate experimental parts, converts units, and saves theprocessed data into individual CSV and JLS files for each sample in `data/exp_pro`.It also generates initial concentration plots for each sample and saves them in the `plots` folder.Additionally, it generates aggregated plots by facies for Experiment A and B using the raw data,matching the style of the linear regression analysis plots.=#

using DrWatson
@quickactivate "AmmerBatch"

using DataFrames, XLSX, CairoMakie
using CSV
using Serialization
using DataInterpolations
using Statistics

# -------------------------------------------------------------------------
# SETUP: Load Info and Prepare Aggregated Plots (Experiment A & B)
# -------------------------------------------------------------------------
println("Setting up aggregated plots...")

# Load info for both parts
info = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), "info"))
info_p2 = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed_part2.xlsx"), "info"))
df_info_all = vcat(info, info_p2)

# Constants for conversion and plotting
Vw = 0.08 # L: initial water volume
colors_agg = [:steelblue, :darkseagreen, :crimson] # For aggregated plots
fontcolor = "#474747"
letters_list = 'a':'z'

# Identify unique facies for Experiment A and B
exp_A_samples = df_info_all[startswith.(df_info_all.Sample, "A"), :]
exp_B_samples = df_info_all[startswith.(df_info_all.Sample, "B"), :]

# Merge C2 and C1 for Exp A (Convention)
exp_A_samples[!, "Facies new"] = replace(exp_A_samples[!, "Facies new"], "C2" => "C1") # Handle potential C2 in Facies new if present, though usually it's in Facies
# Ensure we check the correct column. The files seem to have "Facies" and "Facies new".
# fit_linear_regression_v2.jl uses "Facies new".
unique_facies_A = sort(unique(exp_A_samples[!, "Facies new"]))
unique_facies_B = sort(unique(exp_B_samples[!, "Facies new"]))

# Initialize Figures
figure_A = Figure(backgroundcolor = :transparent, size = (1200, 1600)) # Approximate size
figure_B = Figure(backgroundcolor = :transparent, size = (1200, 1600))

# Helper to create axes dict (Copied and adapted from fit_linear_regression_v2.jl)
function create_facies_axes(fig, u_facies, letters, experiment)
    axes = Dict()
    for (i, f) in enumerate(u_facies)
        row = cld(i, 2)
        col = (i-1)%2 + 1
        
        facies_title = "FT-0" * f[2:end] # e.g., T1 -> FT-01
        
        x_label = row == cld(length(u_facies), 2) ? "Time (days)" : ""
        y_label = col == 1 ? "Concentration (mol kg⁻¹)" : ""

        ax = Axis(fig[row, col], 
                  title = "$(letters[i]). $facies_title", 
                  titlealign = :left,
                  titlesize = 14,
                  titlefont = "Avenir Book",
                  titlecolor = fontcolor,
                  xlabel = x_label, 
                  ylabel = y_label,
                  xlabelsize = 14,
                  ylabelsize = 14,
                  xlabelcolor = fontcolor,
                  ylabelcolor = fontcolor,
                  xticklabelsize = 12,
                  yticklabelsize = 12,
                  xticklabelcolor = fontcolor,
                  yticklabelcolor = fontcolor,
                  xticklabelfont = "Avenir Book",
                  yticklabelfont = "Avenir Book",
                  xlabelfont = "Avenir Book",
                  ylabelfont = "Avenir Book",
                  xgridvisible = false,
                  ygridvisible = true,
                  xgridcolor = (:grey64, 0.6),
                  ygridcolor = (:grey64, 0.6),
                  backgroundcolor = :transparent,
                  # Limits adapted from reference
                  limits = ((-2, 126), (0, 3.6)),
                  xticksvisible = row == cld(length(u_facies), 2),
                  xticklabelsvisible = row == cld(length(u_facies), 2),
                  yticksvisible = true,
                  yticklabelsvisible = true,
                  width = 300,
                  height = 150,
                  xticks = 0:25:175,
                  yticks = 0:0.5:3.5,
                  leftspinecolor = fontcolor,
                  rightspinecolor = fontcolor,
                  topspinecolor = fontcolor,
                  bottomspinecolor = fontcolor,
                  xgridwidth = 0.7,
                  ygridwidth = 0.7
                  )
        axes[f] = ax
    end
    return axes
end

axes_A = create_facies_axes(figure_A, unique_facies_A, letters_list, 'A')
axes_B = create_facies_axes(figure_B, unique_facies_B, letters_list, 'B')

# Borehole markers for Experiment A
boreholes = unique([s[1:2] for s in exp_A_samples.Sample])
markers_list = [:circle, :rect, :dtriangle, :diamond, :utriangle, :pentagon, :hexagon]
borehole_markers = Dict(zip(boreholes, markers_list[1:min(length(boreholes), length(markers_list))]))

# Legends
# Legend elements for Experiment A
leg_elements_A = []
push!(leg_elements_A, [LineElement(color = colors_agg[1], linestyle = :solid), MarkerElement(color = colors_agg[1], marker = :circle, markersize = 8)])
push!(leg_elements_A, [LineElement(color = colors_agg[2], linestyle = :solid), MarkerElement(color = colors_agg[2], marker = :circle, markersize = 8)])
leg_labels_A = ["NO₃⁻", "DOC"]
for bh in boreholes
    push!(leg_elements_A, MarkerElement(color = colors_agg[1], marker = get(borehole_markers, bh, :circle), markersize = 8))
    push!(leg_labels_A, bh)
end
Legend(figure_A[cld(length(unique_facies_A), 2) + 1, 1:2], leg_elements_A, leg_labels_A, framevisible = false, orientation = :horizontal, labelsize = 14, labelcolor = fontcolor, labelfont = "Avenir Book", backgroundcolor = :transparent, nbanks = 2)

# Legend elements for Experiment B
leg_elements_B = [
    [LineElement(color = colors_agg[1], linestyle = :solid), MarkerElement(color = colors_agg[1], marker = :circle, markersize = 8)],
    [LineElement(color = colors_agg[2], linestyle = :solid), MarkerElement(color = colors_agg[2], marker = :circle, markersize = 8)]
]
leg_labels_B = ["NO₃⁻", "DOC"]
Legend(figure_B[cld(length(unique_facies_B), 2) + 1, 1:2], leg_elements_B, leg_labels_B, framevisible = false, orientation = :horizontal, labelsize = 14, labelcolor = fontcolor, labelfont = "Avenir Book", backgroundcolor = :transparent)

# Accumulators for stats
facies_counts_A = Dict{String, Int}()
facies_counts_B = Dict{String, Int}()
facies_toc_A = Dict{String, Vector{Float64}}()
facies_toc_B = Dict{String, Vector{Float64}}()


# -------------------------------------------------------------------------
# PROCESSING: Part 1
# -------------------------------------------------------------------------

println("Processing Part 1...")
samples = info[!, 1]
t = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), "t (days)"))
t = t[!, 1]

# Initialize DataFrame for stats (bar plots)
raw_stats = DataFrame(Sample = String[], Facies = String[], 
                      NO3_decrease_pct = Float64[], DOC_decrease_pct = Float64[])

# Measurements setup (Individual plots)
measurements = ["NO3 (mg L)", "DOC (mg L)", "SO4 (mg L)", "NH4 (mg L)", "N2O (ppm)", "NO2- (mg L)"]
meas_names = ["NO3-", "DOC", "SO4-2", "NH4+", "N2O", "NO2-"]
molar_masses = [62.0049, 12.0107, 96.06, 18.038, 44.013, 46.0055]
colors = [:blue, :green, :red, :purple, :orange, :brown]

dfs = []
for (j,measurement) in enumerate(measurements)
    df = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), measurement,"A:O",header=false))
    df = permutedims(df)
    row1 = df[1, :]
    df = df[2:end, :]
    colnames = Symbol.([string(row1[i]) for i in 1:size(df, 2)])
    rename!(df, colnames)
    df = coalesce.(df, missing)
    if measurement != "N2O (ppm)"
        df = df ./ molar_masses[j]
    else
        df = df./10^6 
    end
    push!(dfs, df)
end

# Save info
CSV.write(datadir("exp_pro","sample_info.csv"), info)

for sample in samples
    df_sample = DataFrame(t = t)
    for df in dfs
        df_sample = hcat(df_sample, df[!,sample], makeunique=true)
    end
    rename!(df_sample, [Symbol("t");Symbol.(meas_names)])
    CSV.write(datadir("exp_pro","$(sample).csv"), df_sample)
    serialize(datadir("exp_pro","$(sample).jls"), df_sample)

    # ------------------- Stats Calculation -------------------
    facies_raw = info[info[!, :Sample] .== sample, "Facies new"][1]
    facies = replace(facies_raw, "C2" => "C1") # Normalize for aggregation

    times = df_sample.t
    no3 = df_sample[!, "NO3-"]
    doc = df_sample[!, "DOC"]

    # NO3 Stats
    idx_no3 = .!ismissing.(no3)
    t_no3 = times[idx_no3]
    v_no3 = Float64.(no3[idx_no3])
    no3_pct = NaN
    if length(v_no3) > 1
        itp_no3 = LinearInterpolation(v_no3, t_no3, extrapolation=ExtrapolationType.Constant)
        c0 = itp_no3(0.0)
        c125 = itp_no3(125.0)
        c125 < 0.0 ? c125 = 0.0 : c125
        if c0 != 0
            no3_pct = (c0 - c125) / c0 * 100
        end
    end

    # DOC Stats
    idx_doc = .!ismissing.(doc)
    t_doc = times[idx_doc]
    v_doc = Float64.(doc[idx_doc])
    doc_pct = NaN
    if length(v_doc) > 1
        itp_doc = LinearInterpolation(v_doc, t_doc, extrapolation=ExtrapolationType.Constant)
        c0 = itp_doc(0.0)
        c125 = itp_doc(125.0)
        c125 < 0.0 ? c125 = 0.0 : c125
        if c0 != 0
            doc_pct = (c0 - c125) / c0 * 100
        end
    end
    push!(raw_stats, [sample, facies, no3_pct, doc_pct])

    # ------------------- Individual Plot -------------------
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mmol L⁻¹)",
    title = "$(sample) initial data",
    xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    xgridwidth = 0.4, ygridwidth = 0.4,)
    xlims!(ax,(-2, 180))
    ylims!(ax,(-0.1, 4.1))
    for (i, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df_sample[!,Symbol(meas_name)])
        ts = df_sample[!,:t][setdiff(1:end, missing_idx)]
        values = df_sample[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0; continue; end
        scatter!(ax, ts, values, label = meas_name, color = colors[i], markersize = 10, marker = :diamond)
        lines!(ax, ts, values, label = meas_name, color = colors[i], linestyle = :dash, linewidth = 2)
    end
    axislegend(ax, position = :rt, merge = true)
    save(plotsdir("individual_plots/$(sample)_initial.png"), fig)

    # ------------------- Aggregated Plotting (Experiment A/B) -------------------
    mass = info[info[!, :Sample] .== sample, "dry weight (g)"][1]
    toc_val = info[info[!, :Sample] .== sample, "TOC %"][1]
    
    # Determine Experiment and Axis
    is_exp_A = startswith(sample, "A")
    current_axes = is_exp_A ? axes_A : axes_B
    counts_dict = is_exp_A ? facies_counts_A : facies_counts_B
    toc_dict = is_exp_A ? facies_toc_A : facies_toc_B
    
    if haskey(current_axes, facies)
        ax_f = current_axes[facies]
        
        # Accumulate stats
        counts_dict[facies] = get(counts_dict, facies, 0) + 1
        if !ismissing(toc_val)
            if !haskey(toc_dict, facies); toc_dict[facies] = Float64[]; end
            push!(toc_dict[facies], toc_val)
        end
        
        # Marker
        marker_shape = :circle
        if is_exp_A
             marker_shape = get(borehole_markers, sample[1:2], :circle)
        end

        # Plot NO3 (mmol/L)
        if !isempty(v_no3)
             lines!(ax_f, t_no3, v_no3, color = colors_agg[1], linewidth = 1.5)
             scatter!(ax_f, t_no3, v_no3, color = colors_agg[1], marker = marker_shape, markersize = 10)
        end
        
        # Plot DOC (mmol/L)
        if !isempty(v_doc)
             lines!(ax_f, t_doc, v_doc, color = colors_agg[2], linewidth = 1.5)
             scatter!(ax_f, t_doc, v_doc, color = colors_agg[2], marker = marker_shape, markersize = 10)
        end
    end
end


# -------------------------------------------------------------------------
# PROCESSING: Part 2
# -------------------------------------------------------------------------

println("Processing Part 2...")
samples = info_p2[!, 1]
t = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed_part2.xlsx"), "t (days)"))
t = t[!, 1]

# Measurements for Part 2 (Subset)
measurements = ["NO3 (mg L)", "DOC (mg L)", "SO4 (mg L)"]
meas_names = ["NO3-", "DOC", "SO4-2",]
molar_masses = [62.0049, 12.0107, 96.06,]
colors = [:blue, :green, :red,]

dfs = []
for (j,measurement) in enumerate(measurements)
    df = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed_part2.xlsx"), measurement,"A:M",header=false))
    df = permutedims(df)
    row1 = df[1, :]
    df = df[2:end, :]
    colnames = Symbol.([string(row1[i]) for i in 1:size(df, 2)])
    rename!(df, colnames)
    df = coalesce.(df, missing)
    df = df ./ molar_masses[j] 
    push!(dfs, df)
end

CSV.write(datadir("exp_pro","sample_info_part2.csv"), info_p2)

for sample in samples
    df_sample = DataFrame(t = t)
    for df in dfs
        df_sample = hcat(df_sample, df[!,sample], makeunique=true)
    end
    rename!(df_sample, [Symbol("t");Symbol.(meas_names)])
    CSV.write(datadir("exp_pro","$(sample).csv"), df_sample)
    serialize(datadir("exp_pro","$(sample).jls"), df_sample)

    # ------------------- Stats Calculation -------------------
    facies_raw = info_p2[info_p2[!, :Sample] .== sample, "Facies new"][1]
    facies = facies_raw # No C2 replacement needed for B usually, or handled globally

    times = df_sample.t
    no3 = df_sample[!, "NO3-"]
    doc = df_sample[!, "DOC"]
    so4 = df_sample[!, "SO4-2"]

    idx_no3 = .!ismissing.(no3)
    t_no3 = times[idx_no3]
    v_no3 = Float64.(no3[idx_no3])
    no3_pct = NaN
    if length(v_no3) > 1
        itp_no3 = LinearInterpolation(v_no3, t_no3, extrapolation=ExtrapolationType.Constant)
        c0 = itp_no3(0.0)
        c125 = itp_no3(125.0)
        c125 < 0.0 ? c125 = 0.0 : c125
        if c0 != 0; no3_pct = (c0 - c125) / c0 * 100; end
    end

    idx_doc = .!ismissing.(doc)
    t_doc = times[idx_doc]
    v_doc = Float64.(doc[idx_doc])
    doc_pct = NaN
    if length(v_doc) > 1
        itp_doc = LinearInterpolation(v_doc, t_doc, extrapolation=ExtrapolationType.Constant)
        c0 = itp_doc(0.0)
        c125 = itp_doc(125.0)
        c125 < 0.0 ? c125 = 0.0 : c125
        if c0 != 0; doc_pct = (c0 - c125) / c0 * 100; end
    end
    push!(raw_stats, [sample, facies, no3_pct, doc_pct])

    # ------------------- Individual Plot -------------------
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mmol L⁻¹)",
    title = "$(sample) initial data",
    xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    xgridwidth = 0.4, ygridwidth = 0.4,)
    xlims!(ax,(-2, 180))
    ylims!(ax,(-0.1, 4.1))
    for (i, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df_sample[!,Symbol(meas_name)])
        ts = df_sample[!,:t][setdiff(1:end, missing_idx)]
        values = df_sample[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0; continue; end
        scatter!(ax, ts, values, label = meas_name, color = colors[i], markersize = 10, marker = :diamond)
        lines!(ax, ts, values, label = meas_name, color = colors[i], linestyle = :dash, linewidth = 2)
    end
    axislegend(ax, position = :rt, merge = true)
    save(plotsdir("individual_plots/$(sample)_initial.png"), fig)

    # ------------------- Aggregated Plotting (Experiment A/B) -------------------
    mass = info_p2[info_p2[!, :Sample] .== sample, "dry weight (g)"][1]
    toc_val = info_p2[info_p2[!, :Sample] .== sample, "TOC %"][1]
    
    is_exp_A = startswith(sample, "A")
    current_axes = is_exp_A ? axes_A : axes_B
    counts_dict = is_exp_A ? facies_counts_A : facies_counts_B
    toc_dict = is_exp_A ? facies_toc_A : facies_toc_B
    
    if haskey(current_axes, facies)
        ax_f = current_axes[facies]
        
        counts_dict[facies] = get(counts_dict, facies, 0) + 1
        if !ismissing(toc_val)
            if !haskey(toc_dict, facies); toc_dict[facies] = Float64[]; end
            push!(toc_dict[facies], toc_val)
        end
        
        marker_shape = :circle # No boreholes for B usually, or standard

        if !isempty(v_no3)
             lines!(ax_f, t_no3, v_no3, color = colors_agg[1], linewidth = 1.5)
             scatter!(ax_f, t_no3, v_no3, color = colors_agg[1], marker = marker_shape, markersize = 10)
        end
        
        if !isempty(v_doc)
             lines!(ax_f, t_doc, v_doc, color = colors_agg[2], linewidth = 1.5)
             scatter!(ax_f, t_doc, v_doc, color = colors_agg[2], marker = marker_shape, markersize = 10)
        end
    end
end

# -------------------------------------------------------------------------
# FINALIZE: Add Text and Save Aggregated Plots
# -------------------------------------------------------------------------
println("Finalizing aggregated plots...")

function add_plot_labels(fig, axes_dict, counts_dict, toc_dict, is_exp_A)
    for (facies, count) in counts_dict
        if haskey(axes_dict, facies)
            ax = axes_dict[facies]
            upper_limit = 3.4
            # Count label
            text!(ax, 125/2-5, upper_limit, text = "n=$count", align = (:center, :top), fontsize = 12, font = "Avenir Book", color = fontcolor)

            # Mean SOC label
            if haskey(toc_dict, facies) && !isempty(toc_dict[facies])
                mean_toc = mean(toc_dict[facies])
                text!(ax, 120, upper_limit, text = "mean SOC: $(round(mean_toc, digits=2))%", align = (:right, :top), fontsize = 12, font = "Avenir Book", color = fontcolor)
            end
        end
    end
    resize_to_layout!(fig)
end

add_plot_labels(figure_A, axes_A, facies_counts_A, facies_toc_A, true)
save(plotsdir("experiment_A_facies_raw.png"), figure_A)
save(plotsdir("experiment_A_facies_raw.pdf"), figure_A)

add_plot_labels(figure_B, axes_B, facies_counts_B, facies_toc_B, false)
save(plotsdir("experiment_B_facies_raw.png"), figure_B)
save(plotsdir("experiment_B_facies_raw.pdf"), figure_B)


# # -------------------------------------------------------------------------
# # STATS: Bar Plots (Percentage Decrease)
# # -------------------------------------------------------------------------

# function plot_decrease_stats_raw(stats_df, exp_label, substance, col_mean, col_std)
#     n_facies = nrow(stats_df)
#     facies_names = stats_df.facies
#     means = stats_df[!, col_mean]
#     stds = replace(stats_df[!, col_std], NaN => 0.0)

#     fig = Figure(size = (800, 600))
#     ax = Axis(fig[1, 1], 
#               xlabel = "Facies", 
#               ylabel = "Decrease after 125 days (%)",
#               title = "Experiment $exp_label: $substance Decrease (Raw Data)",
#               xticks = (1:n_facies, facies_names))
    
#     barplot!(ax, 1:n_facies, means, color = :steelblue, strokecolor = :black, strokewidth = 1)
#     errorbars!(ax, 1:n_facies, means, stds, whiskerwidth = 10, color = :black)
    
#     save(plotsdir("experiment_$(exp_label)_$(substance)_decrease_125d_pct_raw.png"), fig)
#     save(plotsdir("experiment_$(exp_label)_$(substance)_decrease_125d_pct_raw.pdf"), fig)
# end

for (exp_label, df_exp) in [("A", filter(row -> startswith(row.Sample, "A"), raw_stats)), 
                            ("B", filter(row -> startswith(row.Sample, "B"), raw_stats))]
    
    if exp_label == "A"
         df_exp.Facies = replace(df_exp.Facies, "C2" => "C1")
    end
    
    stats_exp = combine(groupby(df_exp, :Facies), 
                        :NO3_decrease_pct => mean => :mean_no3_pct, 
                        :NO3_decrease_pct => std => :std_no3_pct,
                        :DOC_decrease_pct => mean => :mean_doc_pct,
                        :DOC_decrease_pct => std => :std_doc_pct)
    
    rename!(stats_exp, :Facies => :facies)
    sort!(stats_exp, :facies)
    
    println("\nExperiment $exp_label Stats (Decrease after 125 days - Raw):")
    display(stats_exp)
    
    CSV.write(datadir("exp_pro", "experiment_$(exp_label)_stats_125d_raw.csv"), stats_exp)
    # plot_decrease_stats_raw(stats_exp, exp_label, "NO3", :mean_no3_pct, :std_no3_pct)
    # plot_decrease_stats_raw(stats_exp, exp_label, "DOC", :mean_doc_pct, :std_doc_pct)
end
