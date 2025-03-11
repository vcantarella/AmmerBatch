using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
using DataFrames, CSV
using PRIMA
using Statistics
using OrdinaryDiffEq
using CairoMakie
using nonlinearlstr
using ForwardDiff
using Trapz
colors = [:blue, :green]
meas_names = ["NO3-", "DOC"]

R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
rhs!(du, u, p, t) = doc_model!(du, u, p, t, Kh) # importing the model!!
late_times_params = DataFrame(CSV.File(datadir("exp_pro","late_times_params.csv")))
samples = df_info[!, :Sample]

int_df = DataFrame(sample = String[], integral_no3 = Float64[],
    integral_doc = Float64[], weight = Float64[],
    k_no3 = Float64[], doc_ss = Float64[])

for (i, sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]*1e-3
    k_no3 = late_times_params[late_times_params[!, :sample] .== sample, "k_no3"][1]
    no365 = df[df[!, :t] .== 65, "NO3-"][1]*1e-3
    no30 = no365 + k_no3*65
    model_no3(t) = no30 - k_no3*t #zero_order_model
    times = df[!, :t]
    doc_ss = df[times .== 65, :DOC][1]*1e-3
    no3 = df[!, "NO3-"]
    idxs_n = findall(!ismissing, df[!, "NO3-"])
    idxs_c = findall(!ismissing, df[!, "DOC"])
    times_n = times[idxs_n]
    times_c = times[idxs_c]
    no3 = no3[idxs_n].*1e-3
    doc = df[!, "DOC"][idxs_c].*1e-3
    time_idxn = times_n .≤ 65
    time_idxc = times_c .≤ 65
    times_n = times_n[time_idxn]
    times_c = times_c[time_idxc]
    no3 = no3[time_idxn]
    doc = doc[time_idxc]
    # calculate and plot the integral using Trapz
    integral_no3 = trapz(times_n, no3)
    integral_function = trapz(times_n, model_no3.(times_n))
    resulting_integral = (integral_no3 - integral_function)/weight
    integral_doc = trapz(times_c, doc)
    resulting_doc = (integral_doc - doc_ss*(times_c[end] - times_c[1]))/weight
    println("The integral of NO3- for sample $sample is $resulting_integral  [mol/gsoil]")
    println("The integral of DOC for sample $sample is $resulting_doc  [mol/gsoil]")

    
    ### Plotting the results and the data fit.
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.5e-3))
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
            fill_between!(ax, times_n, model_no3.(times_n),
                no3, color = colors[j], alpha = 0.5)
            lines!(ax, df[!,:t], model_no3.(df[!,:t]), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
        elseif j == 2
            fill_between!(ax, times_c, repeat([doc_ss], length(times_c)),
                doc, color = colors[j], alpha = 0.5)
            lines!(ax, df[!,:t], repeat([doc_ss], length(df[!,:t])), label = meas_name, color = colors[j], linestyle = :dash,
                linewidth = 2)
        end
        scatter!(ax, ts, values, label = meas_name, color = colors[j], markersize = 10,
            marker = :diamond)
        
    end
    text!(ax, 0.02, 0.9,
        text="Integral of NO3⁻ = $(round(resulting_integral,digits=2)) [mol/gsoil]",
        color = :black, space = :relative)
    text!(ax, 0.02, 0.85,
        text="Integral of DOC = $(round(resulting_doc, digits=2)) [mol/gsoil]",
        color = :black, space = :relative)
    axislegend(ax, position = :rt, merge = true)
    fig
    save(plotsdir("$(sample)_integral.png"), fig)
    push!(int_df, (sample, resulting_integral, resulting_doc, weight, k_no3, doc_ss))
    


end

CSV.write(datadir("exp_pro","integration_results.csv"), int_df)
