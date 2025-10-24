

using CSV
using DataFrames
using Dates
using Statistics
using Distributions
using Plots
using Interpolations
using Printf

"""
SIHRS model for Washington, Mississippi (Mar 2020-Dec 2021) starting from Patient Zero with stochastic simulations
"""
function sihrs_multiple_simulations()
    # Initialize variables at function level
    N = 43000  # Washington County, Mississippi population (2020 Census)
    s0 = 0.0
    i0 = 0.0
    h0 = 0.0
    r0 = 0.0
    d0 = 0.0
    # Load Washington County, MS data for initial conditions
    try
        data_table = CSV.read("washington_mississippi_combined.csv", DataFrame)
        data_table.date = Date.(data_table.date)
        start_date = Date("2020-03-25")
        
        start_idx = findfirst(data_table.date .== start_date)
        if isnothing(start_idx)
            date_diffs = abs.(data_table.date .- start_date)
            start_idx = argmin(date_diffs)
            @warn "Exact start date not found. Using closest date: $(data_table.date[start_idx])"
        end
        
        real_initial_infected = data_table.cases[start_idx]
        real_initial_dead = data_table.deaths[start_idx]
        
        i0 = real_initial_infected / N
        d0 = real_initial_dead / N
        h0 = 0.0
        r0 = 0.0
        s0 = 1.0 - (i0 + h0 + r0 + d0)
        
        @printf("March 25 initial conditions: I=%d, D=%d, H=%d, R=%d, S=%d\n", 
                real_initial_infected, real_initial_dead, 0, 0, round(Int, s0 * N))
        
    catch e
        @warn "Could not load Washington County, MS real data: $(e)"
        real_initial_infected = 5
        real_initial_dead = 0
        
        i0 = real_initial_infected / N
        d0 = real_initial_dead / N
        h0 = 0.0
        r0 = 0.0
        s0 = 1.0 - (i0 + h0 + r0 + d0)
    end

    # Model parameters
    params = Dict(
        "beta" => 0.195,      # infection rate (β > 0) - Updated for Washington, MS
        "gamma" => 0.165,      # I transition rate (γ > 0) - Updated for Washington, MS
        "alpha" => 0.111,       # H transition rate (α > 0)
        "lambda" => 0.0083,   # R transition rate (Λ > 0)
        "pSI" => 1.00,         # probability of S to I (p_{SI} in (0,1])
        "pII" => 0.0,         # probability of I to I (stay infected)
        "pIH" => 0.1614,        # probability of I to H 
        "pIR" => 0.8367,        # probability of I to R 
        "pID" => 0.0019,       # probability of I to D
        "pHH" => 0.00,        # probability of H to H (stay hospitalized)
        "pHR" => 0.846,      # probability of H to R
        "pHD" => 0.154,      # probability of H to D
        "pRR" => 0.02,        # probability of R to R (stay recovered)
        "pRS" => 0.98,        # probability of R to S
        "tmax" => 620,        # simulation end time (extended for Washington, MS data)
        "s0" => s0,           # initial susceptible proportion
        "i0" => i0,           # initial infected proportion
        "h0" => h0,           # initial hospitalized proportion
        "r0" => r0,           # initial recovered proportion
        "d0" => d0            # initial dead proportion
    )
    
    calculated_R0 = (params["beta"] * params["pSI"]) / (params["gamma"] * (1-params["pII"]))
    @printf("Calculated R0 = %.6f \n", calculated_R0)

    validate_parameters(params)
    
    # Validate initial conditions sum to 1
    if abs((params["s0"] + params["i0"] + params["h0"] + params["r0"] + params["d0"]) - 1.0) > 1e-10
        error("Initial conditions must sum to 1")
    end

    num_simulations = 59
    
    if N <= 0
        error("Population size must be positive integer")
    end
    
    # Store results for all simulations
    all_results = Vector{Any}(undef, num_simulations)
    
    # Run multiple simulations
    try
        @printf("Running %d stochastic simulations for N = %d...\n", num_simulations, N)
        
        for sim_idx in 1:num_simulations
            @printf("Running simulation %d/%d...\n", sim_idx, num_simulations)
            all_results[sim_idx] = sihrs_agent_model(N, params)
        end
        
        @printf("All simulations completed!\n")
        
        # Plot all simulations together
        plot_multiple_simulations(all_results, N, params)
        
    catch e
        @printf("Error occurred: %s\n", e)
        rethrow(e)
    end
end

function validate_parameters(params)
    # Validate rates are positive
    if any([params["beta"], params["gamma"], params["alpha"], params["lambda"]] .<= 0)
        error("All rates (beta, gamma, alpha, lambda) must be positive")
    end
    
    # Validate probabilities are in [0,1]
    probs = [params["pSI"], params["pII"], params["pIH"], params["pIR"], params["pID"],
             params["pHH"], params["pHR"], params["pHD"], params["pRR"], params["pRS"]]
    if any(probs .< 0 .|| probs .> 1)
        error("All probabilities must be in [0,1]")
    end
    
    # Validate probability sums
    if abs((params["pII"] + params["pIH"] + params["pIR"] + params["pID"]) - 1.0) > 1e-10
        error("I transition probabilities must sum to 1")
    end
    if abs((params["pHH"] + params["pHR"] + params["pHD"]) - 1.0) > 1e-10
        error("H transition probabilities must sum to 1")
    end
    if abs((params["pRR"] + params["pRS"]) - 1.0) > 1e-10
        error("R transition probabilities must sum to 1")
    end
end

function sihrs_agent_model(N::Int, params::Dict)
    # SIHRS agent-based stochastic model with death
    # Initial conditions
    s0 = round(Int, params["s0"] * N)
    i0 = round(Int, params["i0"] * N)
    h0 = round(Int, params["h0"] * N)
    r0 = round(Int, params["r0"] * N)
    d0 = round(Int, params["d0"] * N)
    
    # Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0 + d0
    if total != N
        # Add or subtract the difference from the largest compartment
        compartments = [s0, i0, h0, r0, d0]
        largest_idx = argmax(compartments)
        compartments[largest_idx] += (N - total)
        s0, i0, h0, r0, d0 = compartments
    end
    
    if (s0 + i0 + h0 + r0 + d0) != N
        error("Initial conditions must sum to N")
    end
    
    max_events = N * 30
    T = zeros(max_events)
    I_prop = zeros(max_events)
    I_count = zeros(Int, max_events)
    H_count = zeros(Int, max_events)
    D_count = zeros(Int, max_events)
    
    # Initialize agent arrays
    S = collect(1:s0)
    I = collect((s0+1):(s0+i0))
    H = collect((s0+i0+1):(s0+i0+h0))
    R = collect((s0+i0+h0+1):(s0+i0+h0+r0))
    D = collect((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))
    

    t = 0.0
    T[1] = 0.0
    event_count = 1
    

    total_pop = s0 + i0 + h0 + r0 + d0
    I_prop[1] = i0 / total_pop
    I_count[1] = i0
    H_count[1] = h0
    D_count[1] = d0
    
    # Main simulation loop
    while !isempty(I) && t < params["tmax"]
        nS = length(S)
        nI = length(I)
        nH = length(H)
        nR = length(R)
        
        infection_rate = params["pSI"] * params["beta"] * nS * nI / N
        to_susceptible_from_R_rate = params["pRS"] * params["lambda"] * nR
        to_hospital_rate = params["gamma"] * nI * params["pIH"]
        to_recovered_from_I_rate = params["gamma"] * nI * params["pIR"]
        to_dead_from_I_rate = params["gamma"] * nI * params["pID"]
        to_recovered_from_H_rate = params["alpha"] * nH * params["pHR"]
        to_dead_from_H_rate = params["alpha"] * nH * params["pHD"]
        
        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate +
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate +
                     to_dead_from_H_rate
        
        if total_rate == 0
            break
        end
        

        dt = rand(Exponential(1 / total_rate))
        t = t + dt
        
        if t > params["tmax"]
            t = params["tmax"]
            event_count = event_count + 1
            T[event_count] = t
            current_total = length(S) + length(I) + length(H) + length(R) + length(D)
            I_prop[event_count] = length(I) / current_total
            I_count[event_count] = length(I)
            D_count[event_count] = length(D)
            break
        end
        
        event_count = event_count + 1
        T[event_count] = t
        

        chance = rand() * total_rate
        if chance < infection_rate
            # S to I transition
            if nS > 0
                num = rand(1:nS)
                infected_agent = S[num]
                deleteat!(S, num)
                push!(I, infected_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate)
            # R to S transition
            if nR > 0
                num = rand(1:nR)
                susceptible_agent = R[num]
                deleteat!(R, num)
                push!(S, susceptible_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate)
            # I to H transition
            if nI > 0
                num = rand(1:nI)
                hospitalized_agent = I[num]
                deleteat!(I, num)
                push!(H, hospitalized_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate)
            # I to R transition
            if nI > 0
                num = rand(1:nI)
                recovered_agent = I[num]
                deleteat!(I, num)
                push!(R, recovered_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate)
            # I to D transition
            if nI > 0
                num = rand(1:nI)
                dead_agent = I[num]
                deleteat!(I, num)
                push!(D, dead_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate)
            # H to R transition
            if nH > 0
                num = rand(1:nH)
                recovered_agent = H[num]
                deleteat!(H, num)
                push!(R, recovered_agent)
            end
        else
            # H to D transition
            if nH > 0
                num = rand(1:nH)
                dead_agent = H[num]
                deleteat!(H, num)
                push!(D, dead_agent)
            end
        end
        
        current_total = length(S) + length(I) + length(H) + length(R) + length(D)
        I_prop[event_count] = length(I) / current_total
        I_count[event_count] = length(I)
        H_count[event_count] = length(H)
        D_count[event_count] = length(D)
    end
    
    # Trim unused preallocated space
    T = T[1:event_count]
    I_prop = I_prop[1:event_count]
    I_count = I_count[1:event_count]
    H_count = H_count[1:event_count]
    D_count = D_count[1:event_count]
    
    # Store results
    result = Dict(
        "N" => N,
        "T" => T,
        "I_count" => I_count,
        "H_count" => H_count,
        "D_count" => D_count,
        "final_time" => t,
        "peak_infected" => maximum(I_count),
        "peak_time" => T[argmax(I_count)],
        "peak_infected_prop" => maximum(I_prop),
        "peak_time_prop" => T[argmax(I_prop)],
        "i_inf" => I_prop[end]
    )
    
    return result
end

function plot_multiple_simulations(all_results, N, params)
    t_grid = collect(0:params["tmax"])
    all_interp_H = zeros(length(all_results), length(t_grid))
    all_interp_D = zeros(length(all_results), length(t_grid))
    
    for i in 1:length(all_results)
        res = all_results[i]
        # Create interpolation function
        if length(res["T"]) > 1
            itp_H = linear_interpolation(res["T"], res["H_count"], extrapolation_bc=0)
            itp_D = linear_interpolation(res["T"], res["D_count"], extrapolation_bc=0)
            
            for j in 1:length(t_grid)
                all_interp_H[i, j] = itp_H(t_grid[j])
                all_interp_D[i, j] = itp_D(t_grid[j])
            end
        end
    end
    
    # Compute active deaths (rolling window)
    window = 14
    all_active_D = zeros(size(all_interp_D))
    for i in 1:size(all_interp_D, 1)
        for t in 1:length(t_grid)
            if t <= window
                all_active_D[i, t] = all_interp_D[i, t]
            else
                all_active_D[i, t] = all_interp_D[i, t] - all_interp_D[i, t-window]
            end
        end
    end
    
    # Calculate statistics for the bandwidth
    valid_sims = 1:size(all_interp_H, 1)
    
    @printf("Using all %d simulations for prediction interval calculation\n", length(valid_sims))
    

    mean_H = mean(all_interp_H[valid_sims, :], dims=1)[1, :]
    lower_H = [quantile(all_interp_H[valid_sims, t], 0.05) for t in 1:size(all_interp_H, 2)]  # 5th percentile
    upper_H = [quantile(all_interp_H[valid_sims, t], 0.95) for t in 1:size(all_interp_H, 2)]  # 95th percentile
    mean_D = mean(all_active_D, dims=1)[1, :]
    lower_D = [quantile(all_active_D[:, t], 0.05) for t in 1:size(all_active_D, 2)]
    upper_D = [quantile(all_active_D[:, t], 0.95) for t in 1:size(all_active_D, 2)]
    
    # Load and process real-world hospitalization data
    population = N
    real_interp_H = zeros(length(t_grid))
    real_interp_D_prop = zeros(length(t_grid))
    simulation_start_date = Date("2020-03-25")  
    
    try
        hosp_data_table = CSV.read("hospitalization_MS_filtered.csv", DataFrame)
        
        # Remove duplicate entries by summing values for each date (different hospitals)
        hosp_data_table = combine(groupby(hosp_data_table, :collection_week), 
                                 :total_adult_and_pediatric_covid_patients => sum => :total_adult_and_pediatric_covid_patients)
        
        # Convert collection_week to Date with proper year handling
        hosp_data_table.collection_week = Date.(hosp_data_table.collection_week, "m/d/yy")
        # Fix the year to be 20xx instead of 00xx
        for i in 1:length(hosp_data_table.collection_week)
            d = hosp_data_table.collection_week[i]
            hosp_data_table.collection_week[i] = Date(2000 + year(d), month(d), day(d))
        end
        
        # Sort by date
        sort!(hosp_data_table, :collection_week)
        
        # Use the total hospitalization data (7-day average of active cases)
        hospitalization_data = hosp_data_table.total_adult_and_pediatric_covid_patients
        hosp_dates = hosp_data_table.collection_week
        
        # Convert to days from March 25 (simulation start)
        simulation_start_date = Date("2020-03-25")
        hosp_days = [Dates.value(Day(d - simulation_start_date)) for d in hosp_dates]
        
        itp_real_H = linear_interpolation(hosp_days, hospitalization_data, extrapolation_bc=NaN)
        real_interp_H = [itp_real_H(t) for t in t_grid]
        
        # For death data, align with March 25 simulation start
        data_table = CSV.read("washington_mississippi_combined.csv", DataFrame)
        data_table.date = Date.(data_table.date)
        start_idx = findfirst(data_table.date .== simulation_start_date)
        if isnothing(start_idx)
            date_diffs = abs.(data_table.date .- simulation_start_date)
            start_idx = argmin(date_diffs)
        end
        dates_from_start = data_table.date[start_idx:end]
        real_interp_D = interp1([Dates.value(Day(d - simulation_start_date)) for d in dates_from_start], 
                               data_table.deaths[start_idx:end], t_grid)
        
        # Compute active deaths for real data
        real_active_D = zeros(length(real_interp_D))
        for t in 1:length(real_interp_D)
            if t <= window
                real_active_D[t] = real_interp_D[t]
            else
                real_active_D[t] = real_interp_D[t] - real_interp_D[t-window]
            end
        end
        real_interp_D_prop = real_active_D / population
        
        @printf("Successfully loaded hospitalization data with %d data points\n", length(hospitalization_data))
        
    catch e
        @printf("Warning: Could not load or process hospitalization data: %s\n", e)
        real_interp_H = zeros(length(t_grid))
        real_interp_D_prop = zeros(length(t_grid))
    end
    
    all_interp_H_prop = all_interp_H / N
    mean_H_prop = mean_H / N
    lower_H_prop = lower_H / N
    upper_H_prop = upper_H / N
    all_active_D_prop = all_active_D / N
    mean_D_prop = mean_D / N
    lower_D_prop = lower_D / N
    upper_D_prop = upper_D / N
    real_interp_H_prop = real_interp_H / population
    
    # Filter out straight lines (extinct simulations)
    plot_sims = Int[]
    for i in 1:size(all_interp_H_prop, 1)
        sim_data = all_interp_H_prop[i, :]
        
        std_val = std(sim_data)
        max_val = maximum(sim_data)
        min_val = minimum(sim_data)
        range_val = max_val - min_val
        
        initial_val = sim_data[1]
        if std_val > 1e-5 || range_val > 1e-4 || max_val > initial_val * 2
            push!(plot_sims, i)
        end
    end
    
    @printf("Plotting %d out of %d simulations (filtered out %d straight lines)\n",
            length(plot_sims), size(all_interp_H_prop, 1), size(all_interp_H_prop, 1) - length(plot_sims))
    
    if !isempty(plot_sims)
        filtered_max = maximum(maximum(all_interp_H_prop[plot_sims, :]))
        all_max = maximum(maximum(all_interp_H_prop))
        @printf("Filtered simulations max: %.6f, All simulations max: %.6f\n", filtered_max, all_max)
    end

    # Create the final plot with the uncertainty envelope
    p1 = plot(t_grid, upper_H_prop, fillrange=lower_H_prop, 
              fillalpha=0.5, fillcolor=:lightblue, linealpha=0,
              label="90% Prediction Interval", linewidth=0)
    # Plot the real hospitalization data as a solid red line (only where it exists)
    valid_real_data = .!isnan.(real_interp_H_prop)
    if any(valid_real_data)
        plot!(t_grid[valid_real_data], real_interp_H_prop[valid_real_data], 
              color=:red, linewidth=2.5, label="Real Hospitalization Data")
    end
    xlabel!("Time (days)")
    ylabel!("Hospitalized Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params["tmax"])
    ylims!(0, maximum([upper_H_prop; real_interp_H_prop]) * 1.1)
    
    tick_interval = 90
    xtick_positions = collect(0:tick_interval:params["tmax"])
    xtick_dates = [simulation_start_date + Day(x) for x in xtick_positions]
    date_labels = [Dates.format(d, "mm/dd") * "/" * string(year(d))[3:4] for d in xtick_dates]
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p1, "SIHRS_Washington_MS_Hospitalization_bandwidth.png")
    
    # --- 6. Create a second figure with all stochastic sims and the real data ---
    p2 = plot()
    # Plot only non-constant stochastic simulations as semi-transparent blue lines
    for i in plot_sims
        plot!(t_grid, all_interp_H_prop[i, :], color=RGBA(0.2, 0.4, 0.8, 0.3), 
              linewidth=1.0, label="")
    end
    valid_real_data = .!isnan.(real_interp_H_prop)
    if any(valid_real_data)
        plot!(t_grid[valid_real_data], real_interp_H_prop[valid_real_data], 
              color=:red, linewidth=2.5, label="Real Hospitalization Data")
    end
    xlabel!("Time (days)")
    ylabel!("Hospitalized Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params["tmax"])
    ylims!(0, maximum([maximum(all_interp_H_prop[plot_sims, :]); maximum(real_interp_H_prop)]) * 1.1)
    
    # Custom legend
    plot!([NaN], [NaN], color=RGBA(0.2, 0.4, 0.8), linewidth=2.5, 
          label="Stochastic Simulations")
    
    date_labels = [Dates.format(d, "mm/dd") * "/" * string(year(d))[3:4] for d in xtick_dates]
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p2, "SIHRS_Washington_MS_Hospitalization_trajectories.png")
end

# Helper function for interpolation
function interp1(x, y, x_new)
    itp = linear_interpolation(x, y, extrapolation_bc=0)
    return [itp(xi) for xi in x_new]
end

# Run the simulation
if abspath(PROGRAM_FILE) == @__FILE__
    sihrs_multiple_simulations()
end 