using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("const_no3.jl"))
using DataFrames, CSV
using PRIMA
using Statistics
using OrdinaryDiffEq
using CairoMakie
colors = [:blue, :green, :red]
meas_names = ["NO3-", "DOC", "SO4-2"]
samples = ["A301", "A302", "A303", "A306", "A907"]
R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
H = 1/(Kh*R*T) # dimensionless henry's law constant

# Defining model callbacks, where the dilution occurs

params = DataFrame(sample = String[], r_no3 = Float64[])
for sample in samples
    # Load the data
    df = CSV.read(datadir("exp_pro","$(sample).csv"), DataFrame)
    doc_mean = mean(skipmissing(df[!, :DOC]))
    # n2o_mean = mean(skipmissing(df[!, :N2O]))   
    # so4_mean = mean(skipmissing(df[!, "SO4-2"]))

    # Define the initial conditions
    u0 = [df[1, "NO3-"], doc_mean, so4_mean]
    p = [0.01, 0.01, H]
    tspan = (0, maximum(df[!, :t]))
    # define the cost function
    obs = convert.(Float64,skipmissing(df[!, "NO3-"]))
    obs_n2o = convert.(Float64,skipmissing(df[!, "N2O"]))
    idxs = findall(!ismissing, df[!, "NO3-"])
    idxs_n2o = findall(!ismissing, df[!, "N2O"])
    function cost(p)
        prob = ODEProblem(const_denitrification, u0, tspan, p)
        sol = solve(prob, Tsit5(), saveat=df[:, :t])
        solv = vcat(sol.u'...)
        residuals = obs .- solv[idxs, 1]
        residuals_n2o = obs_n2o .- solv[idxs_n2o, 5]

        return sum(abs2, residuals)
    end
    # optimize the parameters
    res = PRIMA.uobyqa(cost, p,)
    p = res[1]
    # solve the ODE
    prob = ODEProblem(const_denitrification, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=df[!, :t])
    solv = vcat(sol.u'...)
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
    push!(params, (sample, p[1]))
end

CSV.write(datadir("exp_pro","params_const_model.csv"), params)