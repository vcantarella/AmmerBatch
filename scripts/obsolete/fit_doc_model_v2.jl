using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("doc_model.jl"))
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
rhs!(du, u, p, t) = doc_model!(du, u, p, t, Kh) # importing the model!!
late_times_params = DataFrame(CSV.File(datadir("exp_pro","late_times_params.csv")))

params = DataFrame(sample = String[], doc_f = Float64[], doc_l = Float64[],
                   k_no3 = Float64[], k_doc = Float64[], αˡ = Float64[],
                   Kd = Float64[], rg = Float64[], K_doc = Float64[],
                   K_no3 = Float64[], k_n2o = Float64[], k_n2o_doc = Float64[], 
                   K_n2o = Float64[],cost_f = Float64[])

# I need this function to find the indices of the time points where the samples were taken.
# then I can compare the model results with the data in the cost function.
function find_first_indices(times, target_times)
    indices = []
    for target_time in target_times
        push!(indices, findfirst(t -> t == target_time, times))
    end
    return indices
end




# For the fit to work, we need good initial conditions for each sample.
samples = df_info[!, :Sample]
for (i, sample) in enumerate(samples)
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    weight = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1]*1e-3
    doc_min = minimum(skipmissing(df[!, :DOC]))*1e-3
    doc_max = maximum(skipmissing(df[!, :DOC]))*1e-3
    #so4_mean = mean(skipmissing(df[!, "SO4-2"]))
    Kd = 1e-4
    no30 = df[1, "NO3-"].*1e-3
    doc_f = doc_min
    doc_l = df[1, :DOC].*1e-3 .- doc_f
    doc_l < 0. && (doc_l = 0.)
    if i == 1
        doc_l = 1e-8
    end
    doc_0 = Kd*(doc_l)
    # Define the initial conditions
    u0 = [no30, doc_l, 0.,0., 0.08, 0.02, doc_0, weight]
    # Define the adjustable model parameters
    k_no3 = late_times_params[late_times_params[!, :sample] .== sample, "k_no3"][1]
    if i == 1
        k_doc = 1e-11
        k_n2o_doc = 1e-11
    else
        k_doc = 1e-6
        k_n2o_doc = 1e-6
    end
    αˡ = 1e-2
    Kdoc = 1e-5
    Kno3 = 1e-3
    k_n2o = 1e-7
    Kn2o = 1e-4
    model_p = [k_no3, k_doc, αˡ, Kd, Kdoc, Kno3, k_n2o, k_n2o_doc, Kn2o]
    # adjustable params: k_doc, Kd, Kdoc, Kno3, k_n2o, k_n2o_doc, Kn2o
    adj_p_min = [1e-12, 1e-12, 1e-8, 1e-8, 1e-12, 1e-12, 1e-8]
    adj_p_max = [1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 5e-5, 1e-3]
    tspan = (0, maximum(df[!, :t]))
    # define the cost function
    ## osberved values
    obs_n = convert.(Float64,skipmissing(df[!, "NO3-"])).*1e-3
    obs_c = convert.(Float64,skipmissing(df[!, "DOC"])).*1e-3
    obs_n2o = convert.(Float64,skipmissing(df[!, "N2O"]))
    ## indices of observed values
    idxs_n = findall(!ismissing, df[!, "NO3-"])
    idxs_c = findall(!ismissing, df[!, "DOC"])
    idxs_n2o = findall(!ismissing, df[!, "N2O"])
    # time points where the liquid and gas samples were taken
    t_liquid = df[idxs_n, :t]
    t_g = df[idxs_n2o, :t]
    tstops = union(t_liquid, t_g)
    # sort tstops
    sort!(tstops)
    ## Defining callbacks. These callbacks model when a sample is taken.
    ## When a liquid sample is taken, V of water decreases and V of gas increases.
    ## Does the gas phase concentration is diluted.
    ## When a gas sample is taken. Local gas is diluted by N2 atmosphere.
    ## Callbacks alter the integrator according to DifferentialEquations.jl
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

    ## Defining the ODEProblem and solving it!.
    problem = ODEProblem(rhs!, u0, tspan, model_p)
    sol = solve(problem, Vern7(), saveat=df[!, :t],
     callback = cb, tstops = tstops,
     abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
    # the full parameter set is the model parameters + the bioavailable doc.
    problem_p = [doc_l, k_doc, Kd, Kdoc, Kno3, k_n2o, k_n2o_doc, Kn2o]
    xl = [0., adj_p_min...]
    xu = [doc_max-doc_f, adj_p_max...]
    
    function cost(p)
        doc_l = p[1]
        adjp_p = p[2:end]
        model_p = [k_no3, adjp_p[1], αˡ, adjp_p[2:end]...] # update the model parameters
        Kd = model_p[4]
        u0 = [no30, doc_l, 0.,0., 0.08, 0.02,Kd*doc_l, weight]
        cost_prob = remake(problem; u0 = u0, p = model_p)
        sol = solve(cost_prob, Vern7(), saveat=df[!, :t],
            callback = cb, tstops = tstops,
            abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
        solv = vcat(sol.u'...)
        idx_no3 = find_first_indices(sol.t, df[idxs_n, :t])
        residuals_n = obs_n .- solv[idx_no3, 1]
        idx_doc = find_first_indices(sol.t, df[idxs_c, :t])
        residuals_c = obs_c .-  doc_f .- solv[idx_doc, 2]
        idx_n2o = find_first_indices(sol.t, df[idxs_n2o, :t])
        residuals_n2o = obs_n2o .- solv[idx_n2o, 4]
        return sum(abs2, residuals_n) + sum(abs2, residuals_c)*0.1 + sum(abs2, residuals_n2o)*0.1
    end

    cost(problem_p)

    # optimize the parameters
    println("Fitting $(sample)")
    function residuals(p)
        doc_l = p[1]
        adjp_p = p[2:end]
        model_p = [k_no3, adjp_p[1], αˡ, adjp_p[2:end]...] # update the model parameters
        Kd = model_p[4]
        u0 = [no30, doc_l, 0.,0., 0.08, 0.02,Kd*doc_l, weight]
        cost_prob = remake(problem; u0 = u0, p = model_p)
        sol = solve(cost_prob, Rosenbrock23(), saveat=df[!, :t],
            callback = cb, tstops = tstops,
            abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
        solv = vcat(sol.u'...)
        idx_no3 = find_first_indices(sol.t, df[idxs_n, :t])
        residuals_n = obs_n .- solv[idx_no3, 1]
        idx_doc = find_first_indices(sol.t, df[idxs_c, :t])
        residuals_c = obs_c .-  doc_f .- solv[idx_doc, 2]
        idx_n2o = find_first_indices(sol.t, df[idxs_n2o, :t])
        residuals_n2o = obs_n2o .- solv[idx_n2o, 4]
        return vcat(residuals_n, residuals_c, residuals_n2o)
    end
    residuals(problem_p)
    jaco(p) = ForwardDiff.jacobian(residuals, p)
    res = nonlinearlstr.nlss_bounded_trust_region(residuals, jaco, problem_p, xl, xu, initial_radius = 2e-6)
    p = res[1]
    
    # solve the ODE
    adjp_p = p[2:end]
    model_p = [k_no3, adjp_p[1], αˡ, adjp_p[2:end]...]
    Kd = model_p[4]
    doc_l = p[1]
    u0 = [no30, doc_l, 0.,0., 0.08, 0.02, doc_l*Kd, weight]
    prob = ODEProblem(rhs!, u0, tspan, model_p)
    sol = solve(prob, Vern7(), saveat=0.5,
        callback = cb, tstops = tstops,
        abstol = 1e-9, reltol = 1e-9, maxiters = 10000)
    solv = vcat(sol.u'...)
    solv[:, 2] .+= doc_f

    cost_f = cost(p)


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
    #save the plots
    save(plotsdir("$(sample)_model.png"), fig)
    #save the parameters
   
    push!(params,  (sample, doc_f, doc_l, model_p[1:4]..., 1e-2, model_p[5:end]..., cost_f))
end

CSV.write(datadir("exp_pro","params_doc_model.csv"), params)