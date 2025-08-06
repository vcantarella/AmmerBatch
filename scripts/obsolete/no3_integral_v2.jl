using DrWatson
using Printf
@quickactivate "AmmerBatch"
using DataFrames, CSV
using Statistics
using CairoMakie
using Trapz
using DataInterpolations
using QuadGK
using Roots
using XLSX
using LaTeXStrings
colors = [:blue, :green]
meas_names = ["NO3-", "DOC"]


R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
linear_regression_params = DataFrame(CSV.File(datadir("exp_pro","linear_regression_params.csv")))
samples = df_info[!, :Sample]
int_df = linear_regression_params
@. int_df[!, :integral_no3] = NaN
@. int_df[!, :c_quick] = NaN
@. int_df[!, :t_quick] = NaN
int_df[!, :weight] = df_info[!, "dry weight (g)"]
int_df[!, :facies] = df_info[!, :Facies]
int_df[!, :TOC] = df_info[!, "TOC %"]

# Besides a grapth for each integration I want to add the axis to larger grid plots (according to size of the number of plots there may be more than on plot)
grid_plots = cld(length(samples), 8)
inch = 96
pt = 4/3
cm = inch / 2.54

grid_plot_figs = [Figure(size=(18cm, 27cm)) for _ in 1:grid_plots]
total_plots = length(samples)
for (i, sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]*1e-3
    r_no3 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "r_no3"][1]
    no30 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "c₀"][1]
    t0 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "t_zeroorder"][1]
    times = df[!, :t]
    no3 = df[!, "NO3-"]
    idxs_n = findall(!ismissing, df[!, "NO3-"])
    idxs_c = findall(!ismissing, df[!, "DOC"])
    times_n = times[idxs_n]
    times_c = times[idxs_c]
    no3 = no3[idxs_n].*1e-3
    doc = df[!, "DOC"][idxs_c].*1e-3
    data_fit = LinearInterpolation(no3, times_n; extrapolation = ExtrapolationType.Constant)
    time_idxc = times_c .> t0
    times_c = times_c[time_idxc]
    doc = doc[time_idxc]
    c0 = no3[1]
    model_no3(t) = no30 - r_no3*t
    # calculate when data intercepts the model:
    println("Sample: $(sample)")
    t_of_zero = t0 > 0 ? find_zeros(t-> data_fit(t)-model_no3(t), 0, 120)[1] : 0.0
    # calculate and plot the integral using Trapz
    integral_no3 = quadgk(data_fit, 0.0, t_of_zero)[1]
    #trapz(times_n, no3)
    c_quick = t0 > 0 ? c0 - no30 : 0.0
    int_df[int_df[!,:sample] .== sample, :c_quick] .= c_quick
    integral_function = quadgk(model_no3, 0.0, t_of_zero)[1]
    # trapz(times_n, model_no3.(times_n))
    resulting_integral = (integral_no3 - integral_function)
    int_df[int_df[!,:sample] .== sample, :integral_no3] .= resulting_integral
    t_quick = t0 > 0 ? resulting_integral/c_quick : 0.0
    int_df[int_df[!,:sample] .== sample, :t_quick] .= t_quick
    doc_ss = mean(doc)
    ### Plotting the results and the data fit.
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    # based on i calculate the grid plot number and position of the axis inside it.
    # 2 cols and 4 rows
    grid_plot_num = ceil(Int, i/8)
    grid_pos = ceil(Int, (i-1) % 8) + 1
    grid_plot_row =   ceil(Int, grid_pos/2)
    grid_plot_col = grid_pos % 2 == 0 ? 2 : 1
    axg = Axis(grid_plot_figs[grid_plot_num][grid_plot_row, grid_plot_col], xlabel = "Time (days)", ylabel = "Concentration (mol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.8e-3))
    xlims!(ax,(-2, 180))
    ylims!(axg, (0, 3.8e-3))
    xlims!(axg,(-2, 180))
    for (j, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df[!,Symbol(meas_name)])
        ts = df[!,:t][setdiff(1:end, missing_idx)]
        values = df[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0
            continue
        end
        values = values.*1e-3
        if j == 1
            lines!(ax, df[!,:t], model_no3.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
            lines!(axg, df[!,:t], model_no3.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
            if t0 > 0
                tfill = 0:0.1:t_of_zero
                fill_between!(ax, tfill, model_no3.(tfill),
                    data_fit.(tfill), color = colors[j], alpha = 0.5)
                fill_between!(axg, tfill, model_no3.(tfill),
                    data_fit.(tfill), color = colors[j], alpha = 0.5)
            end
        elseif j == 2
            lines!(ax, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
            lines!(axg, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
        end
        scatter!(ax, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
        scatter!(axg, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
    end
    if t0 > 0
        text!(ax, 0.02, 0.9,
        text="t_quick = $(round(t_quick,digits=2)) [days]",
        color = :black, space = :relative)
        text!(ax, 0.02, 0.85,
        text="c_quick =  "*@sprintf("%.2E", c_quick)*" [mol L⁻¹]",
        color = :black, space = :relative)
        # text!(axg, 0.5, 0.8,
        # text="t_quick = $(round(t_quick,digits=2)) [days]",
        # color = :black, space = :relative, fontsize = 9)
        # text!(axg, 0.5, 0.7,
        # text="c_quick =  "*@sprintf("%.2E", c_quick)*" [mol L⁻¹]",
        # color = :black, space = :relative, fontsize = 9)
    end
    text!(ax, 0.5, 0.9,
        text="R² = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²"][1], digits=2))",
        color = :black, space = :relative)
    text!(axg, 0.5, 0.9,
        text="R² = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²"][1], digits=2))",
        color = :black, space = :relative)

    axislegend(ax, position = :rt, merge = true, labelsize = 9)
    Legend(grid_plot_figs[grid_plot_num][5,1:2], axg, "Substance", merge = true, framevisible = false, orientation = :horizontal)
    fig
    save(plotsdir("$(sample)_integral.png"), fig)
    
end
[resize_to_layout!(grid_plot_figs[i]) for i in 1:grid_plots]
[save(plotsdir("grid_plot_A_$(i).png"), grid_plot_figs[i], px_per_unit = 400/inch) for i in 1:grid_plots]

XLSX.writetable(datadir("exp_pro","integration_results.xlsx"), int_df; overwrite = true)

## Part 2;
df_info = CSV.read(datadir("exp_pro","sample_info_part2.csv"), DataFrame)
linear_regression_params = DataFrame(CSV.File(datadir("exp_pro","linear_regression_params_part2.csv")))
samples = df_info[!, :Sample]
int_df = linear_regression_params
@. int_df[!, :integral_no3] = NaN
@. int_df[!, :c_quick] = NaN
@. int_df[!, :t_quick] = NaN
@. int_df[!, :r_so4_no3] = NaN
int_df[!, :weight] = df_info[!, "dry weight (g)"]
int_df[!, :facies] = df_info[!, :Facies]
int_df[!, :TOC] = df_info[!, "TOC %"]

colors = [:blue, :green, :red]
meas_names = ["NO3-", "DOC", "SO4-2"]
# Besides a grapth for each integration I want to add the axis to larger grid plots (according to size of the number of plots there may be more than on plot)
grid_plots = cld(length(samples), 8)
grid_plot_figs = [Figure(size=(18cm, 27cm)) for _ in 1:grid_plots]
total_plots = length(samples)


for (i, sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]*1e-3
    r_no3 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "r_no3"][1]
    no30 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "c₀"][1]
    t0 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "t_zeroorder"][1]
    times = df[!, :t]
    no3 = df[!, "NO3-"]
    idxs_n = findall(!ismissing, df[!, "NO3-"])
    idxs_c = findall(!ismissing, df[!, "DOC"])
    times_n = times[idxs_n]
    times_c = times[idxs_c]
    no3 = no3[idxs_n].*1e-3
    doc = df[!, "DOC"][idxs_c].*1e-3
    data_fit = LinearInterpolation(no3, times_n; extrapolation = ExtrapolationType.Constant)
    time_idxc = times_c .> t0
    times_c = times_c[time_idxc]
    doc = doc[time_idxc]
    c0 = no3[1]
    model_no3(t) = no30 - r_no3*t
    # calculate when data intercepts the model:
    println("Sample: $(sample)")
    t_of_zero = t0 > 0 ? find_zeros(t-> data_fit(t)-model_no3(t), 0, 120)[1] : 0.0
    # calculate and plot the integral using Trapz
    integral_no3 = quadgk(data_fit, 0.0, t_of_zero)[1]
    #trapz(times_n, no3)
    c_quick = t0 > 0 ? c0 - no30 : 0.0
    int_df[int_df[!,:sample] .== sample, :c_quick] .= c_quick
    integral_function = quadgk(model_no3, 0.0, t_of_zero)[1]
    # trapz(times_n, model_no3.(times_n))
    resulting_integral = (integral_no3 - integral_function)
    int_df[int_df[!,:sample] .== sample, :integral_no3] .= resulting_integral
    t_quick = t0 > 0 ? resulting_integral/c_quick : 0.0
    int_df[int_df[!,:sample] .== sample, :t_quick] .= t_quick
    doc_ss = mean(doc)
    # Analyse the SO4-2 data:
    # compare the r_so4 with the r_no3
    r_so4 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "r_so4"][1]
    #r_so4 = r_so4 > 0 ? r_so4 : 0.0
    s0 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "s0"][1]
    
    model_so4(t) = s0 + r_so4*t
    @show r_so4/r_no3
    int_df[int_df[!,:sample] .== sample, :r_so4_no3] .= r_so4/r_no3

    ### Plotting the results and the data fit.
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    grid_plot_num = ceil(Int, i/8)
    grid_pos = ceil(Int, (i-1) % 8) + 1
    grid_plot_row =   ceil(Int, grid_pos/2)
    grid_plot_col = grid_pos % 2 == 0 ? 2 : 1
    axg = Axis(grid_plot_figs[grid_plot_num][grid_plot_row, grid_plot_col], xlabel = "Time (days)", ylabel = "Concentration (mol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.8e-3))
    xlims!(ax,(-2, 180))
    ylims!(axg, (0, 3.8e-3))
    xlims!(axg,(-2, 180))
    for (j, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df[!,Symbol(meas_name)])
        ts = df[!,:t][setdiff(1:end, missing_idx)]
        values = df[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0
            continue
        end
        values = values.*1e-3
        if j == 1
            lines!(ax, df[!,:t], model_no3.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
            lines!(axg, df[!,:t], model_no3.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
            if t0 > 0
                tfill = 0:0.1:t_of_zero
                fill_between!(ax, tfill, model_no3.(tfill),
                    data_fit.(tfill), color = colors[j], alpha = 0.5)
                fill_between!(axg, tfill, model_no3.(tfill),
                    data_fit.(tfill), color = colors[j], alpha = 0.5)
            end
        elseif j == 2
            lines!(ax, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
            lines!(axg, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
        elseif j == 3
            lines!(ax, df[!,:t], model_so4.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
            lines!(axg, df[!,:t], model_so4.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                    linewidth = 2)
        end
        scatter!(ax, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
        scatter!(axg, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
        
    end
    if t0 > 0
        text!(ax, 0.02, 0.9,
        text="t_quick = $(round(t_quick,digits=2)) [days]",
        color = :black, space = :relative)
        text!(ax, 0.02, 0.85,
        text="c_quick =  "*@sprintf("%.2E", c_quick)*" [mol L⁻¹]",
        color = :black, space = :relative)
        # text!(axg, 0.2, 0.8,
        # text="t_quick = $(round(t_quick,digits=2)) [days]",
        # color = :black, space = :relative, fontsize = 9)
        # text!(axg, 0.2, 0.7,
        # text="c_quick =  "*@sprintf("%.2E", c_quick)*" [mol L⁻¹]",
        # color = :black, space = :relative, fontsize = 9)
    end
    text!(ax, 0.5, 0.9,
        text="R² = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²"][1], digits=2))",
        color = :black, space = :relative)
    text!(ax, 0.5, 0.85,
        text="R²ₛ = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²_so4"][1], digits=2))",
        color = :black, space = :relative)
    text!(ax, 0.5, 0.80,
        text="r_so4/r_no3 = $(round(r_so4/r_no3, digits=2))",
        color = :black, space = :relative)
    text!(axg, 0.5, 0.9,
        text="R² = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²"][1], digits=2))",
        color = :black, space = :relative)
    # text!(axg, 0.8, 0.8,
    #     text="R²ₛ = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²_so4"][1], digits=2))",
    #     color = :black, space = :relative,  fontsize = 9)
    # text!(axg, 0.8, 0.70,
    #     text=L"$\frac{r_{so4}}{r_{no3}}$ = %$(round(r_so4/r_no3, digits=2))",
    #     color = :black, space = :relative, fontsize = 9)

    axislegend(ax, position = :rt, merge = true)
    Legend(grid_plot_figs[grid_plot_num][5,1:2], axg, "Substance", merge = true, framevisible = false, orientation = :horizontal)
    fig
    save(plotsdir("$(sample)_integral.png"), fig)
    
end
[resize_to_layout!(grid_plot_figs[i]) for i in 1:grid_plots]
[save(plotsdir("grid_plot_B_$(i).png"), grid_plot_figs[i]) for i in 1:grid_plots]
XLSX.writetable(datadir("exp_pro","integration_results_part2.xlsx"), int_df; overwrite = true)