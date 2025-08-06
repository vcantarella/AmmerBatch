using DrWatson
@quickactivate "AmmerBatch"
using DataFrames, CSV
using PRIMA
using Statistics
using CairoMakie
using ForwardDiff
colors = [:blue, :green, :orange]
meas_names = ["NO3-", "DOC", "N2O"]

R = 0.0821 # atm L / mol K
T = 298.15 # K
Kh = 0.025 # mol L⁻¹ atm⁻¹ henry's law constant
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
t_change = [0, 0, 0, 0, 65, 65, 10, 20, 50, 30, 0, 0, 30, 65, 0, 50, 65, 20, 105, 10, 30, 80, 0, 100,]
late_times_params = DataFrame(sample = String[], r_no3 = Float64[], c₀ = Float64[], R² = Float64[], t_zeroorder = Float64[])
samples = df_info[!, :Sample]
for (i,sample) in enumerate(samples)
    df = CSV.read(datadir("exp_pro","$sample.csv"), DataFrame)
    t_change_i = t_change[i]
    modelled_times = df[df[!, :t].≥t_change_i, :t]
    modelled_no3 = df[df[!, :t].≥t_change_i, "NO3-"]
    index = convert.(Bool, 1 .-(ismissing.(modelled_no3)))
    modelled_no3 = modelled_no3[index].*1e-3 # convert to mol L⁻¹
    modelled_times = modelled_times[index]
    # preparing linear regression
    X = [ones(size(modelled_times)) modelled_times]
    y = modelled_no3
    # linear regression and model fit results:
    β = (X'X)\(X'y)
    ϵ = y - X*β
    s² = ϵ'ϵ/(length(y)-length(β))
    σ² = s²*(length(y)-length(β))/length(y)
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSR = sum(ϵ.^2)
    R² = 1 - SSR/SST
    push!(late_times_params, [sample, -β[2], β[1], R², t_change_i])

end
late_times_params
CSV.write(datadir("exp_pro","linear_regression_params.csv"), late_times_params)


colors = [:blue, :green, :red]
meas_names = ["NO3-", "DOC", "SO4-2"]
df_info = CSV.read(datadir("exp_pro","sample_info_part2.csv"), DataFrame)
t_change = 1
late_times_params = DataFrame(sample = String[], r_no3 = Float64[], c₀ = Float64[], R² = Float64[], t_zeroorder = Float64[],
    s0 = Float64[], r_so4 = Float64[], R²_so4 = Float64[])
samples = df_info[!, :Sample]
for (i,sample) in enumerate(samples)
    df = CSV.read(datadir("exp_pro","$sample.csv"), DataFrame)
    #t_change_i = t_change[i]
    modelled_times = df[df[!, :t].≥t_change, :t]
    modelled_no3 = df[df[!, :t].≥t_change, "NO3-"]
    index = convert.(Bool, 1 .-(ismissing.(modelled_no3)))
    modelled_no3 = modelled_no3[index].*1e-3 # convert to mol L⁻¹
    modelled_times = modelled_times[index]
    # preparing linear regression
    X = [ones(size(modelled_times)) modelled_times]
    y = modelled_no3
    # linear regression and model fit results:
    β = (X'X)\(X'y)
    ϵ = y - X*β
    s² = ϵ'ϵ/(length(y)-length(β))
    σ² = s²*(length(y)-length(β))/length(y)
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSR = sum(ϵ.^2)
    R² = 1 - SSR/SST
    # sulfur regression:
    modelled_so4 = df[df[!, :t].≥t_change, "SO4-2"]
    s0 = df[1, "SO4-2"]*1e-3 # convert to mol L⁻¹
    index = convert.(Bool, 1 .-(ismissing.(modelled_so4)))
    modelled_times = df[df[!, :t].≥t_change, :t]
    modelled_so4 = modelled_so4[index].*1e-3 # convert to mol L⁻¹
    modelled_times = modelled_times[index]
    y = modelled_so4
    t = modelled_times
    r_so4 = (sum(y.*t)-s0*sum(t))/
        (sum(t.^2))
    ϵ = y .- (s0 .+ t .* r_so4)
    s² = ϵ'ϵ/(length(y)-1) # 1 degree of freedom for the slope
    σ² = s²*(length(y)-1)/length(y)
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSR = sum(ϵ.^2)
    R²_so4 = 1 - SSR/SST
    push!(late_times_params, [sample, -β[2], β[1], R², t_change, s0, r_so4, R²_so4])
end
late_times_params
CSV.write(datadir("exp_pro","linear_regression_params_part2.csv"), late_times_params)
