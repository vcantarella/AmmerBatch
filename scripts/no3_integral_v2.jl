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
for (i, sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]*1e-3
    k_no3 = linear_regression_params[linear_regression_params[!, :sample] .== sample, "k_no3"][1]
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
    model_no3(t) = no30 - k_no3*t
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
    ylims!(ax, (0, 3.8e-3))
    xlims!(ax,(-2, 180))
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
            if t0 > 0
                tfill = 0:0.1:t_of_zero
                fill_between!(ax, tfill, model_no3.(tfill),
                    data_fit.(tfill), color = colors[j], alpha = 0.5)
                
            end
        elseif j == 2
            lines!(ax, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
        end
        scatter!(ax, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
        
    end
    if t0 > 0
        text!(ax, 0.02, 0.9,
        text="t_quick = $(round(t_quick,digits=2)) [days]",
        color = :black, space = :relative)
        text!(ax, 0.02, 0.85,
        text="c_quick =  "*@sprintf("%.2E", c_quick)*" [mol L⁻¹]",
        color = :black, space = :relative)
    end
    text!(ax, 0.5, 0.9,
        text="R² = $(round(linear_regression_params[linear_regression_params[!,:sample] .== sample, "R²"][1], digits=2))",
        color = :black, space = :relative)

    axislegend(ax, position = :rt, merge = true)
    fig
    save(plotsdir("$(sample)_integral.png"), fig)
    
end

CSV.write(datadir("exp_pro","integration_results.csv"), int_df)
