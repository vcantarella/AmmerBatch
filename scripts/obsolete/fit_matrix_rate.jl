using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("doc_model.jl"))
include(srcdir("const_no3.jl"))
using DataFrames, CSV
using PRIMA
using Statistics
using OrdinaryDiffEq
using CairoMakie
using nonlinearlstr
using ForwardDiff
colors = [:blue, :green, :orange]
meas_names = ["NO3-", "DOC", "N2O"]

R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)

rhs!(du, u, p, t) = const_denitrification!(du, u, p, t, Kh) # importing the model!!
samples = df_info[!, :Sample]
late_times_params = DataFrame(sample = String[], k_no3 = Float64[])
params = DataFrame(sample = String[], doc_f = Float64[], k_no3 = Float64[],
                   k_doc = Float64[], αˡ = Float64[], Kd = Float64[],
                   K_doc = Float64[], K_no3 = Float64[],
                   k_n2o = Float64[],k_n2o_doc = Float64[], 
                   K_n2o = Float64[],cost_f = Float64[])
for sample in samples
    df = CSV.read(datadir("exp_pro","$sample.csv"), DataFrame)
    no30 = df[df[!, :t] .== 65, "NO3-"][1]*1e-3
    doc_f = minimum(skipmissing(df[!, :DOC]))*1e-3
    doc_0 = df[df[!, :t] .== 65, "DOC"][1]*1e-3 .- doc_f
    doc_0 < 0. && (doc_0 = 0.)
    u0 = [no30, doc_0, 0.0, 0.0, 1.0, 1.0]
    modelled_times = df[df[!, :t].≥65, :t]
    modelled_no3 = df[df[!, :t].≥65, "NO3-"]
    index = convert.(Bool, 1 .-(ismissing.(modelled_no3)))
    modelled_no3 = modelled_no3[index].*1e-3
    modelled_times = modelled_times[index]
    tspan = (modelled_times[1], modelled_times[end])
    p0 = [1e-5, 1e-6]
    prob = ODEProblem(rhs!, u0, tspan, p0)
    sol = solve(prob, Tsit5(), saveat = 5)
    
    cost(p) = begin
        k_no3 = p[1]
        k_n2o = 1e-8
        p = [k_no3, k_n2o]
        cost_prob = remake(prob, p = p)
        sol = solve(cost_prob, Tsit5(), saveat = modelled_times)
        solv = vcat(sol.u'...)
        no3_model = solv[:,1] # modelled NO3- concentrations
        cos = sum(abs2, modelled_no3 .- no3_model)
        return cos
    end
    p0 = [p0[1]]
    res = PRIMA.bobyqa(cost, p0, xl = [1e-8], xu = [1e-3])
    p = res[1]
    push!(late_times_params, [sample, p[1]])
    prob = remake(prob, p = [p[1], 1e-8])
    sol = solve(prob, Tsit5(), saveat = modelled_times)
    solv = vcat(sol.u'...)

    # plot the results
    ### Plotting the results and the data fit.
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mmol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.5e-3))
    xlims!(ax,(-2, 180))
    ax2 = Axis(fig[1, 1], yaxisposition = :right, ygridvisible = false,
     xgridvisible = false, ylabelcolor = :orange, ylabel = "N2O [atm]")
    ylims!(ax2, (0, 5.1e-4))
    xlims!(ax2,(-2, 180))
    for (i, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df[!,Symbol(meas_name)])
        ts = df[!,:t][setdiff(1:end, missing_idx)]
        values = df[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0
            continue
        end
        if meas_name == "N2O"
            lines!(ax2, sol.t, solv[:,4], label = meas_name, color = colors[i], linestyle = :dash,)
            scatter!(ax2, ts, values, label = meas_name, color = colors[i], markersize = 10,
            marker = :diamond)
        else
            values = values.*1e-3
            lines!(ax, sol.t, solv[:,i], label = meas_name, color = colors[i], linestyle = :dash,
                linewidth = 2)
            scatter!(ax, ts, values, label = meas_name, color = colors[i], markersize = 10,
                marker = :diamond)
        end
    end
    axislegend(ax, position = :rt, merge = true)
    fig

end
late_times_params
CSV.write(datadir("exp_pro","late_times_params.csv"), late_times_params)
