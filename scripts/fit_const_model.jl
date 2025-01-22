using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("const_no3.jl"))
using DataFrames, CSV
using LinearAlgebra
using PRIMA
using Statistics
using OrdinaryDiffEq
using CairoMakie
using ForwardDiff
using SciMLSensitivity
using nonlinearlstr
colors = [:blue, :green, :orange]
meas_names = ["NO3-", "DOC", "N2O"]
samples = ["A301", "A302", "A303", "A306", "A907"]
R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
# reaction rate model
rhs!(du, u, p, t) = const_denitrification(du, u, p, t, Kh)

function find_first_indices(times, target_times)
    indices = []
    for target_time in target_times
        push!(indices, findfirst(t -> t == target_time, times))
    end
    return indices
end

df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
# Defining model callbacks, where the dilution occurs

params = DataFrame(sample = String[], r_no3 = Float64[], r_n2o = Float64[], r_g = Float64[])
for (k,sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    doc_mean = mean(skipmissing(df[!, :DOC]))
    # n2o_mean = mean(skipmissing(df[!, :N2O]))   
    # so4_mean = mean(skipmissing(df[!, "SO4-2"]))
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]

    # Define the initial conditions
    u0 = [df[1, "NO3-"].*1e-3, doc_mean.*1e-3, 0., 0., 0.08, 0.02]
    p = [6e-6, 1e-4, 1e-4]
    if k == 1
        p = [6e-6, 5e-4, 1e-4]
    elseif k == 2
        p = [6e-6, 1e-4, 1e-4]
    elseif k == 3
        p = [6e-6, 2.5e-5, 1e-4]
    elseif k == 4
        p = [6e-6, 1.e-5, 1e-4]
    elseif k == 5
        p = [6e-6, 2e-5, 1e-5]
    end
    tspan = (0, maximum(df[!, :t]))
    # define the cost function
    obs = convert.(Float64,skipmissing(df[!, "NO3-"])).*1e-3
    obs_n2o = convert.(Float64,skipmissing(df[!, "N2O"]))
    idxs = findall(!ismissing, df[!, "NO3-"])
    idxs_n2o = findall(!ismissing, df[!, "N2O"])
    t_liquid = df[idxs, :t]
    t_g = df[idxs_n2o, :t]
    tstops = union(t_liquid, t_g)
    # sort tstops
    sort!(tstops)
    sampling_callback_condition(u, t, integrator) = t ∈ t_liquid
    dilution_callback_condition(u, t, integrator) = t ∈ t_g
    function sampling_affect!(integrator)
        Vg_old = integrator.u[6]
        Vw_old = integrator.u[5]
        Vg_new = Vg_old + 0.004
        integrator.u[4] = integrator.u[4] * (Vg_old / Vg_new)
        integrator.u[6] = Vg_new
        integrator.u[5] = Vw_old - 0.004
    end
    function dilution_affect!(integrator)
        Vg_old = integrator.u[6]
        integrator.u[4] = integrator.u[4] * (1-0.004/Vg_old)
    end
    cb_sampling = DiscreteCallback(sampling_callback_condition, sampling_affect!)
    cb_dilution = DiscreteCallback(dilution_callback_condition, dilution_affect!)
    cb = CallbackSet(cb_sampling, cb_dilution)
    function cost(p)
        prob = ODEProblem(rhs!, u0, tspan, p)
        sol = solve(prob, Tsit5(), saveat=df[:, :t], callback = cb,
         tstops = tstops, abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
        solv = vcat(sol.u'...)
        idx_no3 = find_first_indices(sol.t, df[idxs, :t])
        residuals = obs .- solv[idx_no3, 1]
        idx_n2o = find_first_indices(sol.t, df[idxs_n2o, :t])
        residuals_n2o = obs_n2o .- solv[idx_n2o, 4]
        return sum(abs2, residuals)*10 + sum(abs2, residuals_n2o)
    end

    grad(p) = ForwardDiff.gradient(cost, p)
    hess(p) = ForwardDiff.hessian(cost, p)

    # optimize the parameters
    lb = ones(length(p)).*1e-10
    ub = [1e-4, 1e-3, 1e-3]
    res = nonlinearlstr.bounded_trust_region(cost, grad, hess, p, lb, ub, initial_radius = 1e-5)
    # scale = 1 ./p
    # res = PRIMA.bobyqa(cost, p,
    #     scale = scale,
    #     xl = lb, xu = ub, rhobeg = 1e-8)
    p = res[1]
    # solve the ODE
    prob = ODEProblem(rhs!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=df[:, :t], callback = cb,
         tstops = tstops, abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
    solv = vcat(sol.u'...)
    # Plot the data for each sample
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Concentration (mmol L⁻¹)",
    title = "$(sample) Model Results",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.5e-3))
    ax2 = Axis(fig[1, 1], yaxisposition = :right, ygridvisible = false,
     xgridvisible = false, ylabelcolor = :orange, ylabel = "N2O [atm]")
    ylims!(ax2, (0, 5.1e-4))
    #xlims!(ax,(-2, 180))
    #ylims!(ax,(-0.1, u0[1]+0.4))
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
    #save the plots
    save(plotsdir("$(sample)_model.png"), fig)
    #save the parameters
    push!(params, (sample, p[1], p[2], p[3]))
end

CSV.write(datadir("exp_pro","params_const_model.csv"), params)