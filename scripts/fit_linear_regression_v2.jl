using DrWatson
@quickactivate "AmmerBatch"
using Roots
using DataFrames, CSV
using PRIMA
using Statistics
using CairoMakie
using ForwardDiff
using DataInterpolations
using QuadGK
colors = [:blue, :green, :red]
meas_names = ["NO3-", "DOC", "SO4-2"]

# Part 1:
Vw = 0.08 # L: initial water volume in the batch
df_info = CSV.read(datadir("exp_pro","sample_info.csv"), DataFrame)
df_info_p2 = CSV.read(datadir("exp_pro","sample_info_part2.csv"), DataFrame)
# join the two dataframes
df_info = vcat(df_info, df_info_p2)
late_times_params = DataFrame(sample = String[], facies = String[], r_no3 = Float64[], c₀ = Float64[], R² = Float64[], t_zeroorder = Float64[],
                              integral_no3 = Float64[], c_quick = Float64[], t_quick = Float64[],
                              weight = Float64[], TOC = Union{Missing, Float64}[])
samples = df_info[!, :Sample]

# Besides a grapth for each integration I want to add the axis to larger grid plots (according to size of the number of plots there may be more than on plot)
grid_plots = cld(length(samples), 8)
inch = 96
pt = 4/3
cm = inch / 2.54

grid_plot_figs = [Figure(size=(18cm, 27cm)) for _ in 1:grid_plots]
total_plots = length(samples)
for i in 1:total_plots
    sample = samples[i]
    mass = df_info[df_info[!, :Sample] .== sample, "dry weight (g)"][1] # g
    mass = mass * 1e-3 # convert to kg
    df = CSV.read(datadir("exp_pro","$sample.csv"), DataFrame)
    # t_change_i = t_change[i]
    all_times = df[!, :t]
    modelled_no3 = df[!, "NO3-"]
    index_no3 = convert.(Bool, 1 .-(ismissing.(modelled_no3)))
    index_no3 = findall(index_no3) # indices of non-missing NO₃⁻ values
    #modelled_no3 = modelled_no3[index_no3].*1e-3 # convert to mol L⁻¹
    modelled_times = all_times[index_no3]
    data_fit = LinearInterpolation(modelled_no3[index_no3].*1e-3, modelled_times; extrapolation = ExtrapolationType.Constant)
    
    doc = df[!, "DOC"]
    index_doc = convert.(Bool, 1 .-(ismissing.(doc)))
    index_doc = findall(index_doc) # indices of non-missing DOC values
    #doc = doc[index_doc].*1e-3 # convert to mol L⁻¹
    doc_times = all_times[index_doc]
    doc_fit = LinearInterpolation(doc[index_doc].*1e-3, doc_times; extrapolation = ExtrapolationType.Constant)
    # correct for the reduction in water volume due to sample extraction
    vws = Vw * ones(size(all_times))
    volume_removed = 0.005 * ones(size(all_times)) # L: volume removed for analysis (4 ml)
    mass_removed = data_fit.(all_times) .* volume_removed # mols: mass of NO₃⁻ removed at each time point
    doc_removed = doc_fit.(all_times) .* volume_removed # mols: mass of DOC removed at each time point
    cum_vol_removed = cumsum(volume_removed)
    cum_mass_removed = cumsum(mass_removed) # mols: cumulative mass of NO₃⁻ removed
    cum_doc_removed = cumsum(doc_removed) # mols: cumulative mass of DOC removed
    vws = Vw .- cum_vol_removed
    no3_mol_kg = zeros(length(modelled_times)) # preallocate for NO₃⁻ mol kg⁻¹
    k = 1
    for ind in index_no3
        if ind == 1
            no3_mol_kg[1] = modelled_no3[ind]*1e-3 .* 0.08 ./ mass # convert to mol kg⁻¹
        else
            no3_mol_kg[k] = (modelled_no3[ind]*1e-3 .* vws[ind-1] + cum_mass_removed[ind-1]) ./ mass # convert to mol kg⁻¹
        end
        k += 1
    end
    doc_mol_kg = zeros(length(doc_times)) # preallocate for DOC mol kg⁻¹
    k = 1
    for ind in index_doc
        if ind == 1
            doc_mol_kg[1] = doc[ind]*1e-3 .* 0.08 ./ mass # convert to mol kg⁻¹
        else
            doc_mol_kg[k] = (doc[ind]*1e-3 .* vws[ind-1] + cum_doc_removed[ind-1]) ./ mass # convert to mol kg⁻¹
        end
        k += 1
    end
    if sample[1] == 'B' # then we also have sulfate data
        so4 = df[!, "SO4-2"]
        index_so4 = convert.(Bool, 1 .-(ismissing.(so4)))
        index_so4 = findall(index_so4) # indices of non-missing SO₄²⁻ values
        so4_times = all_times[index_so4]
        so4_mol_kg = zeros(length(so4_times)) # preallocate for SO₄²⁻ mol kg⁻¹
        so4_fit = LinearInterpolation(so4[index_so4].*1e-3, so4_times; extrapolation = ExtrapolationType.Constant)
        so4_removed = so4_fit.(all_times) .* volume_removed # mols: mass of SO₄²⁻ removed at each time point
        cum_so4_removed = cumsum(so4_removed) # mols: cumulative mass of SO₄²⁻ removed
        k = 1
        for ind in index_so4
            if ind == 1
                so4_mol_kg[1] = so4[ind]*1e-3 .* 0.08 ./ mass # convert to mol kg⁻¹
            else
                so4_mol_kg[k] = (so4[ind]*1e-3 .* vws[ind-1] + cum_so4_removed[ind-1]) ./ mass # convert to mol kg⁻¹
            end
            k += 1
        end
    end

    # Now I need to define when the experiment became zero-order.
    # I will interpolate the corrected DOC data and then find the time when it crosses the line
    # that represent the constant DOC value (average of the last 2 DOC values).
    # This is when we interpret the experiment as zero-order because the NO₃⁻ reduction is then
    # not limited by the DOC concentration anymore, but to the matrix hydrolysis of e- donors.
    no3_fit = LinearInterpolation(no3_mol_kg, modelled_times; extrapolation = ExtrapolationType.Constant)
    doc_fit = LinearInterpolation(doc_mol_kg, doc_times; extrapolation = ExtrapolationType.Constant)
    # find the last two DOC values and calculate the average
    doc_last_two = doc_mol_kg[end-2:end]
    doc_avg = mean(doc_last_two)
    # find the time when the DOC concentration crosses the average
    t_change_i = findfirst(doc_fit.(all_times) .<= doc_avg)
    if t_change_i === nothing
        t_change_i = length(all_times) # if it never crosses, then we take the last time point
    end
    t_change = all_times[t_change_i] # convert to time
    # handpicking a few samples where my automated criteria did not work
    if sample == "A908" # potential bad_t_change
        global t_change = 100.0
    elseif sample == "A509" # potential bad_t_change
        global t_change = 50.0
    elseif sample == "A508" # potential bad_t_change
        global t_change = 30.0
    elseif sample == "A301" # potential bad_t_change
        global t_change = 50.0  
    end
    println("Sample: $(sample), t_change: $(t_change)")
    # calculate the slope of the NO₃⁻ reduction
    # Linear regression:
    X = [ones(size(modelled_times[modelled_times.>t_change])) modelled_times[modelled_times.>t_change]]
    y = no3_mol_kg[modelled_times.>t_change]
    # linear regression and model fit results:
    β = (X'X)\(X'y)
    ϵ = y - X*β
    s² = ϵ'ϵ/(length(y)-length(β))
    σ² = s²*(length(y)-length(β))/length(y)
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSR = sum(ϵ.^2)
    R² = 1 - SSR/SST
    
    model_no3(t) = β[1] + β[2]*t
    # calculate when data intercepts the model:
    println("Sample: $(sample)")

    t_of_zero = t_change > 0 ? find_zeros(t-> no3_fit(t)-model_no3(t), 0, 150)[1] : 0.0
    # calculate and plot the integral using Trapz
    integral_no3 = quadgk(no3_fit, 0.0, t_of_zero)[1]
    c0 = no3_fit(0.0) # initial concentration
    c_quick = t_change > 0 ? c0 - β[1] : 0.0
    
    integral_function = quadgk(model_no3, 0.0, t_of_zero)[1]
    # trapz(times_n, model_no3.(times_n))
    resulting_integral = (integral_no3 - integral_function)
   
    t_quick = t_change > 0 ? resulting_integral/c_quick : 0.0
    facies = df_info[df_info[!, :Sample] .== sample, "Facies"][1]
    push!(late_times_params, [sample, facies, -β[2], β[1], R², t_change, 
                              resulting_integral, c_quick, t_quick,
                              mass, df_info[df_info[!, :Sample] .== sample, "TOC %"][1]])

    
    # 2 cols and 4 rows
    grid_plot_num = ceil(Int, i/8)
    grid_pos = ceil(Int, (i-1) % 8) + 1
    grid_plot_row =   ceil(Int, grid_pos/2)
    grid_plot_col = grid_pos % 2 == 0 ? 2 : 1
    ax = Axis(grid_plot_figs[grid_plot_num][grid_plot_row, grid_plot_col], xlabel = "Time (days)", ylabel = "Concentration (mol kg⁻¹)",
    title = "$(sample) - facies: $(facies)",
    #xticks = 0:10:170, yticks = 0:0.5:4,
    xgridstyle = :dash, ygridstyle = :dash,
    #xgridwidth = 0.4, ygridwidth = 0.4,
    )
    ylims!(ax, (0, 3.8e-3*Vw/mass)) # mol kg⁻¹
    xlims!(ax,(-2, 180))
    if t_change > 0
            tfill = 0:0.1:t_of_zero
            fill_between!(ax, tfill, model_no3.(tfill),
                no3_fit.(tfill), color = colors[1], alpha = 0.5)
    end
    lines!(ax, all_times, model_no3.(all_times), color = colors[1], label = "Model NO3-")
    scatter!(ax, modelled_times, no3_mol_kg, color = colors[1], label = meas_names[1])
    scatter!(ax, doc_times, doc_mol_kg, color = colors[2], label = meas_names[2])
    if sample[1] == 'B' # then we also have sulfate data
        scatter!(ax, so4_times, so4_mol_kg, color = colors[3], label = meas_names[3])
    end
    text!(ax, 0.5, 0.9,
        text="R² = $(round(R², digits=2))",
        color = :black, space = :relative)
    Legend(grid_plot_figs[grid_plot_num][5,1:2], ax, "Substance", merge = true, framevisible = false, orientation = :horizontal)
    #Label(grid_plot_figs[grid_plot_num][1,1:2, Top()], "Normalized NO3- and DOC plots", fontsize = 18, halign = :center)
end
[resize_to_layout!(grid_plot_figs[i]) for i in 1:grid_plots]
[save(plotsdir("grid_plot_temp_$(i).png"), grid_plot_figs[i], px_per_unit = 400/inch) for i in 1:grid_plots]
CSV.write(datadir("exp_pro","linear_regression_params_v2.csv"), late_times_params)
