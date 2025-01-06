using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("doc_model.jl"))
using DataFrames, CSV
using PRIMA
using Statistics
using OrdinaryDiffEq
using CairoMakie
colors = [:blue, :green, :red]
meas_names = ["NO3-", "DOC", "SO4-2"]
samples = ["A308", "A310", "A402","A403", "A405", "A406", "A502", "A503",
           "A504", "A507", "A508", "A509", "A510","A603", "A604",
           "A904", "A905", "A906", "A907", "A908"]

params = DataFrame(sample = String[], r_no3 = Float64[], r_doc = Float64[], αˡ = Float64[],
                   c_eq = Float64[], K_doc = Float64[], K_no3 = Float64[], doc_f = Float64[],
                   doc_0 = Float64[])
for sample in samples
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    doc_min = minimum(skipmissing(df[!, :DOC]))
    doc_max = maximum(skipmissing(df[!, :DOC]))
    # n2o_mean = mean(skipmissing(df[!, :N2O]))   
    so4_mean = mean(skipmissing(df[!, "SO4-2"]))
    doc_0 = 2
    Kd = doc_0/(doc_max-doc_min)
    # Define the initial conditions
    u0 = [df[1, "NO3-"], df[1, :DOC], so4_mean, doc_0]
    model_p = [0.01, 1, 0.1, Kd, 1e-4, 1e-5]
    tspan = (0, maximum(df[!, :t]))
    # define the cost function
    obs_n = convert.(Float64,skipmissing(df[!, "NO3-"]))
    obs_c = convert.(Float64,skipmissing(df[!, "DOC"]))
    idxs_n = findall(!ismissing, df[!, "NO3-"])
    idxs_c = findall(!ismissing, df[!, "DOC"])
    problem = ODEProblem(doc_model, u0, tspan, model_p)
    sol = solve(problem, Tsit5(), saveat=df[!, :t])
    p = [doc_0, doc_min, model_p...]
    xl = [0, doc_min-0.2, 1e-8, 1e-8, 1e-8, 1e-5, 1e-7, 1e-7]
    xu = [4, df[1, :DOC]-0.1, 8, 8, 8, 1e3, 1e-3, 1e-3]
    function cost(p)
        doc_0 = p[1]
        doc_f = p[2]
        model_p = p[3:end]
        u0 = [df[1, "NO3-"], df[1, :DOC]- doc_f, so4_mean, doc_0]
        prob = remake(problem; u0 = u0, p = model_p)
        sol = solve(prob, Rosenbrock23(), saveat=df[!, :t], abstol=1e-7, reltol=1e-7, maxiters=1e5)
        solv = vcat(sol.u'...)
        residuals_n = obs_n .- solv[idxs_n, 1]
        residuals_c = obs_c .-  doc_f .- solv[idxs_c, 2]
        return sum(abs2, residuals_n) + sum(abs2, residuals_c)
    end
    # optimize the parameters
    res = PRIMA.bobyqa(cost, p, xl = repeat([1e-8], length(p)), xu = repeat([10], length(p)))
    p = res[1]
    # solve the ODE
    doc_0 = p[1]
    doc_f = p[2]
    u0 = [df[1, "NO3-"], df[1, :DOC]- doc_f, so4_mean, doc_0]
    model_p = p[3:end]
    prob = ODEProblem(doc_model, u0, tspan, model_p)
    sol = solve(prob, Tsit5(), saveat=df[!, :t])
    solv = vcat(sol.u'...)
    solv[:, 2] .+= doc_f
    # Plot the data for each sample
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mmol L⁻¹)",
    title = "$(sample) Model Results",
    xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    xgridwidth = 0.4, ygridwidth = 0.4,)
    xlims!(ax,(-2, 180))
    ylims!(ax,(-0.1, u0[1]+0.4))
    for (i, meas_name) in enumerate(meas_names)
        missing_idx = findall(ismissing, df[!,Symbol(meas_name)])
        ts = df[!,:t][setdiff(1:end, missing_idx)]
        values = df[setdiff(1:end, missing_idx), Symbol(meas_name)]
        values = convert.(Float64, values)
        if size(values, 1) == 0
            continue
        end
        lines!(ax, sol.t, solv[:,i], label = meas_name, color = colors[i], linestyle = :dash,
        linewidth = 2)
        scatter!(ax, ts, values, label = meas_name, color = colors[i], markersize = 10,
        marker = :diamond)
    end
    axislegend(ax, position = :rt, merge = true)
    fig
    #save the plots
    save(plotsdir("$(sample)_model.png"), fig)
    #save the parameters
    push!(params, (sample, model_p..., doc_f, doc_0))
end

CSV.write(datadir("exp_pro","params_doc_model.csv"), params)