using DrWatson
@quickactivate "AmmerBatch"

# Here you may include files from the source directory
include(srcdir("ode_model.jl"))
using DataFrames, XLSX, CairoMakie
using CSV
using Serialization
println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)
info = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), "info"))
samples = info[!, 1]
t = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), "t (days)"))
t = t[!, 1]
# iterating over the measurements table and extracting the data as Dataframes.
measurements = ["NO3 (mg L)", "DOC (mg L)", "SO4 (mg L)", "NH4 (mg L)", "N2O (ppm)", "NO2- (mg L)"]
meas_names = ["NO3-", "DOC", "SO4-2", "NH4+", "N2O", "NO2-"]
molar_masses = [62.0049, 12.0107, 96.06, 18.038, 44.013, 46.0055]
colors = [:blue, :green, :red, :purple, :orange, :brown]
dfs = []
for (j,measurement) in enumerate(measurements)
    df = DataFrame(XLSX.readtable(datadir("exp_raw","johann_batch_preprocessed.xlsx"), measurement,"A:O",header=false))
    # change the column names
    df = permutedims(df)
    # set the first row as the column names
    row1 = df[1, :]
    df = df[2:end, :]
    colnames = Symbol.([string(row1[i]) for i in 1:size(df, 2)])
    rename!(df, colnames)
    # replace missing values with NaNs
    df = coalesce.(df, missing)
    # divide by the molar masses
    if measurement != "N2O (ppm)"
        df = df ./ molar_masses[j]
    else
        p = df./10^6 # partial pressure in atm (ppm to atm)
        R = 0.0821 # atm L / mol K
        T = 298.15 # K
        V = 0.02 # L
        n = p .* (V / (R * T))
        k_h = 0.024 # mol L⁻¹ atm⁻¹ henry's law constant
        C_w = p .* k_h # mol L⁻¹ water Concentration in equilibrium with the gas phase
        V_w = 0.08 # L water volume
        n_w = C_w .* V_w # mol of N2O in the water
        df = (n .+ n_w)./V_w .*1e3 # mmol L⁻¹ total N2O concentration
    end
    push!(dfs, df)
end
# save info df to csv
CSV.write(datadir("exp_pro","sample_info.csv"), info)
# transpose the dataframes
for sample in samples
    df_sample = DataFrame(t = t)
    for df in dfs
        df_sample = hcat(df_sample, df[!,sample], makeunique=true)
    end
    rename!(df_sample, [Symbol("t");Symbol.(meas_names)])
    CSV.write(datadir("exp_pro","$(sample).csv"), df_sample)
    serialize(datadir("exp_pro","$(sample).jls"), df_sample)

    # Plot the data for each sample
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
        if size(values, 1) == 0
            continue
        end
        scatter!(ax, ts, values, label = meas_name, color = colors[i], markersize = 10,
        marker = :diamond)
        lines!(ax, ts, values, label = meas_name, color = colors[i], linestyle = :dash,
        linewidth = 2)
    end
    axislegend(ax, position = :rt, merge = true)
    save(plotsdir("$(sample)_initial.png"), fig)
end

