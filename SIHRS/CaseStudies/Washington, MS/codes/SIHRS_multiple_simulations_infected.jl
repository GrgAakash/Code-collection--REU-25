using CSV
using DataFrames
using Dates
using Statistics
using Distributions
using Plots
using Interpolations
using Printf

# =============================================================================
# TYPE DEFINITIONS
# =============================================================================

# SIHRS model parameters struct for better performance and type safety
struct SIHRSParams
    # Transition rates
    beta::Float64      # infection rate (β > 0)
    gamma::Float64     # I transition rate (γ > 0)
    alpha::Float64     # H transition rate (α > 0)
    lambda::Float64    # R transition rate (Λ > 0)
    
    # Transition probabilities
    pSI::Float64       # probability of S to I
    pII::Float64       # probability of I to I (stay infected)
    pIH::Float64       # probability of I to H
    pIR::Float64       # probability of I to R
    pID::Float64       # probability of I to D
    pHH::Float64       # probability of H to H (stay hospitalized)
    pHR::Float64       # probability of H to R
    pHD::Float64       # probability of H to D
    pRR::Float64       # probability of R to R (stay recovered)
    pRS::Float64       # probability of R to S
    
    # Simulation parameters
    tmax::Int          # simulation end time
    
    # Initial conditions
    s0::Float64        # initial susceptible proportion
    i0::Float64        # initial infected proportion
    h0::Float64        # initial hospitalized proportion
    r0::Float64        # initial recovered proportion
    d0::Float64        # initial dead proportion
end

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Helper function for interpolation
function interp1(x, y, x_new)
    itp = linear_interpolation(x, y, extrapolation_bc=0)
    return [itp(xi) for xi in x_new]
end

# Helper function for date tick setup 
function setup_date_ticks(params::SIHRSParams, start_date_real, tick_interval=90)
    xtick_positions = collect(0:tick_interval:params.tmax)
    xtick_dates = [start_date_real + Day(x) for x in xtick_positions]
    date_labels = [Dates.format(d, "mm/dd") * "/" * string(year(d))[3:4] for d in xtick_dates]
    return xtick_positions, date_labels
end

# Helper function for rolling window calculation 
function compute_rolling_window(data, window_size)
    result = copy(data)
    for t in (window_size + 1):length(data)
        result[t] = data[t] - data[t - window_size]
    end
    return result
end

# Helper function for quantile calculations 
function compute_quantiles(data, q_lower=0.05, q_upper=0.95)
    lower = mapslices(x -> quantile(x, q_lower), data, dims=1)[1, :]
    upper = mapslices(x -> quantile(x, q_upper), data, dims=1)[1, :]
    return lower, upper
end

# =============================================================================
# MAIN SIMULATION FUNCTION
# =============================================================================

"""
SIHRS model for Washington, Mississippi (Mar 2020-Dec 2021) starting from Patient Zero with stochastic simulations
Focuses on infected cases and deaths (like the MATLAB version)
"""
function sihrs_multiple_simulations_infected()
    # Initialize variables at function level
    N = 43000  # Washington County, Mississippi population (2020 Census)
    
    # Initialize default values for initial conditions
    s0 = 0.0
    i0 = 0.0
    h0 = 0.0
    r0 = 0.0
    d0 = 0.0
    
    # --- Load and process real Washington County, MS data for initial conditions ---
    try
        # Load cases and deaths data
        data_table = CSV.read("washington_mississippi_combined.csv", DataFrame)
        data_table.date = Date.(data_table.date)  # Parse dates as Date objects
        start_date = Date("2020-03-25")
        
        # Find start index
        start_idx = findfirst(data_table.date .== start_date)
        if isnothing(start_idx)
            # Find closest date
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
        
    catch e
        @warn "Could not load Washington County, MS real data: $(e)"
        # Fallback to hardcoded values
        real_initial_infected = 1
        real_initial_dead = 0
        
        i0 = real_initial_infected / N
        d0 = real_initial_dead / N
        h0 = 0.0
        r0 = 0.0
        s0 = 1.0 - (i0 + h0 + r0 + d0)
    end

    # Define model parameters using struct for better performance and type safety
    params = SIHRSParams(
        0.255,        # beta: infection rate (β > 0) 
        0.219,        # gamma: I transition rate (γ > 0) 
        0.1,        # alpha: H transition rate (α > 0)
        0.0083,     # lambda: R transition rate (Λ > 0) - Updated for Washington, MS
        1.0,        # pSI: probability of S to I (p_{SI} in (0,1])
        0.00,       # pII: probability of I to I (stay infected)
        0.04,     # pIH: probability of I to H - Updated from P(IH) calculation
        0.959,     # pIR: probability of I to R - Updated to sum to 1
        0.001,     # pID: probability of I to D - Updated from P(ID) calculation
        0.00,       # pHH: probability of H to H (stay hospitalized) - Updated
        0.9982,      # pHR: probability of H to R - Updated
        0.0018,      # pHD: probability of H to D - Updated
        0.02,       # pRR: probability of R to R (stay recovered)
        0.98,       # pRS: probability of R to S
        630,        # tmax: simulation end time (extended for Washington, MS data)
        s0,         # s0: initial susceptible proportion
        i0,         # i0: initial infected proportion
        h0,         # h0: initial hospitalized proportion
        r0,         # r0: initial recovered proportion
        d0          # d0: initial dead proportion
    )
    
    # Verify R0 calculation
    calculated_R0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII)
    @printf("Calculated R0 = %.6f \n", calculated_R0)

    # Validate parameters
    validate_parameters(params)
    
    # Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-10
        error("Initial conditions must sum to 1")
    end

    num_simulations = 59  # Matching MATLAB
    
    # Input validation
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
            all_results[sim_idx] = sihrs_agent_model_infected(N, params)
        end
        
        @printf("All simulations completed!\n")
        
        # Plot all simulations together
        plot_multiple_simulations_infected(all_results, N, params)
        
    catch e
        @printf("Error occurred: %s\n", e)
        rethrow(e)
    end
end

function validate_parameters(params::SIHRSParams)
    # Validate rates are positive
    if any([params.beta, params.gamma, params.alpha, params.lambda] .<= 0)
        error("All rates (beta, gamma, alpha, lambda) must be positive")
    end
    
    # Validate probabilities are in [0,1]
    probs = [params.pSI, params.pII, params.pIH, params.pIR, params.pID,
             params.pHH, params.pHR, params.pHD, params.pRR, params.pRS]
    if any(probs .< 0 .|| probs .> 1)
        error("All probabilities must be in [0,1]")
    end
    
    # Validate probability sums
    if abs((params.pII + params.pIH + params.pIR + params.pID) - 1.0) > 1e-10
        error("I transition probabilities must sum to 1")
    end
    if abs((params.pHH + params.pHR + params.pHD) - 1.0) > 1e-10
        error("H transition probabilities must sum to 1")
    end
    if abs((params.pRR + params.pRS) - 1.0) > 1e-10
        error("R transition probabilities must sum to 1")
    end
end

function sihrs_agent_model_infected(N::Int, params::SIHRSParams)
    # SIHRS agent-based stochastic model with death (focusing on infected and deaths)
    
    # Initial conditions - using params values and ensuring they sum to N
    s0 = round(Int, params.s0 * N)  # susceptible
    i0 = round(Int, params.i0 * N)  # infected
    h0 = round(Int, params.h0 * N)  # hospitalized
    r0 = round(Int, params.r0 * N)  # recovered
    d0 = round(Int, params.d0 * N)  # dead
    
    # Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0 + d0
    if total != N
        # Add or subtract the difference from the largest compartment
        compartments = [s0, i0, h0, r0, d0]
        largest_idx = argmax(compartments)
        compartments[largest_idx] += (N - total)
        s0, i0, h0, r0, d0 = compartments
    end
    
    # Validate initial conditions sum to N
    if (s0 + i0 + h0 + r0 + d0) != N
        error("Initial conditions must sum to N")
    end
    
    # Preallocate arrays for better performance
    # For a 620-day simulation, people can get infected multiple times
    # Each person might have ~10 events: S→I→H→R→S→I→R→S→I→D
    max_events = N * 10  # Estimate maximum number of events
    T = zeros(max_events)
    I_prop = zeros(max_events)
    I_count = zeros(Int, max_events)
    D_count = zeros(Int, max_events)
    
    # Initialize agent arrays
    S = collect(1:s0)
    I = collect((s0+1):(s0+i0))
    H = collect((s0+i0+1):(s0+i0+h0))
    R = collect((s0+i0+h0+1):(s0+i0+h0+r0))
    D = collect((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))
    
    # Initialize time tracking
    t = 0.0
    T[1] = 0.0
    event_count = 1
    
    # Initialize proportion tracking
    total_pop = s0 + i0 + h0 + r0 + d0
    I_prop[1] = i0 / total_pop
    I_count[1] = i0
    D_count[1] = d0
    
    # Main simulation loop
    while !isempty(I) && t < params.tmax
        nS = length(S)
        nI = length(I)
        nH = length(H)
        nR = length(R)
        
        # Calculate event rates according to the mathematical model
        infection_rate = params.pSI * params.beta * nS * nI / N  # S to I rate
        to_susceptible_from_R_rate = params.pRS * params.lambda * nR  # R to S rate
        to_hospital_rate = params.gamma * nI * params.pIH  # I to H rate
        to_recovered_from_I_rate = params.gamma * nI * params.pIR  # I to R rate
        to_dead_from_I_rate = params.gamma * nI * params.pID  # I to D rate
        to_recovered_from_H_rate = params.alpha * nH * params.pHR  # H to R rate
        to_dead_from_H_rate = params.alpha * nH * params.pHD  # H to D rate
        
        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate +
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate +
                     to_dead_from_H_rate
        
        if total_rate == 0
            break
        end
        
        # Time of next event
        dt = rand(Exponential(1 / total_rate))
        t = t + dt
        
        if t > params.tmax
            t = params.tmax
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
        
        # Determine which event occurs
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
        
        # Update tracking arrays
        current_total = length(S) + length(I) + length(H) + length(R) + length(D)
        I_prop[event_count] = length(I) / current_total
        I_count[event_count] = length(I)
        D_count[event_count] = length(D)
    end
    
    # Trim unused preallocated space
    T = T[1:event_count]
    I_prop = I_prop[1:event_count]
    I_count = I_count[1:event_count]
    D_count = D_count[1:event_count]
    
    # Store results
    result = Dict(
        "N" => N,
        "T" => T,
        "I_count" => I_count,
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

function plot_multiple_simulations_infected(all_results, N, params::SIHRSParams)
    # This function creates multiple figures matching the MATLAB version:
    # 1. Uncertainty envelope plot for infected cases
    # 2. All stochastic trajectories for infected cases
    # 3. Active deaths plot
    # 4. Cumulative deaths plot

    # --- 1. Create a common time grid for comparison ---
    t_grid = collect(0:params.tmax)  # daily time steps
    
    # --- 2. Interpolate each simulation's results onto the grid (optimized) ---
    all_interp_I = zeros(length(all_results), length(t_grid))
    all_interp_D = zeros(length(all_results), length(t_grid))
    
    for i in 1:length(all_results)
        res = all_results[i]
        # Create interpolation function and vectorize the evaluation
        if length(res["T"]) > 1
            itp_I = linear_interpolation(res["T"], res["I_count"], extrapolation_bc=0)
            itp_D = linear_interpolation(res["T"], res["D_count"], extrapolation_bc=0)
            
            # Vectorized interpolation - much faster than loop
            all_interp_I[i, :] = itp_I.(t_grid)
            all_interp_D[i, :] = itp_D.(t_grid)
        end
    end
    
    # --- Compute active deaths (rolling window) for each simulation (optimized) ---
    window = 14
    all_active_D = zeros(size(all_interp_D))
    for i in 1:size(all_interp_D, 1)
        all_active_D[i, :] = compute_rolling_window(all_interp_D[i, :], window)
    end
    
    # --- 3. Calculate statistics for the bandwidth (optimized) ---
    mean_I = mean(all_interp_I, dims=1)[1, :]
    lower_I, upper_I = compute_quantiles(all_interp_I)
    mean_D = mean(all_active_D, dims=1)[1, :]
    lower_D, upper_D = compute_quantiles(all_active_D)
    
    # --- 4. Load and process real-world data ---
    population = N  # Use same population as simulation
    real_interp_I = zeros(length(t_grid))
    real_interp_D = zeros(length(t_grid))
    real_interp_D_prop = zeros(length(t_grid))
    real_cumulative_D = zeros(length(t_grid))  # Initialize cumulative deaths
    start_date_real = Date("2020-03-25")  # Default start date
    
    try
        # Load cases and deaths data
        data_table = CSV.read("washington_mississippi_combined.csv", DataFrame)
        data_table.date = Date.(data_table.date)  # Parse dates as Date objects
        start_idx = findfirst(data_table.date .== start_date_real)
        if isnothing(start_idx)
            date_diffs = abs.(data_table.date .- start_date_real)
            start_idx = argmin(date_diffs)
        end
        
        dates_from_start = data_table.date[start_idx:end]
        cumulative_cases_from_start = data_table.cases[start_idx:end]
        recovery_days = 14  # Longer recovery period for rural areas (14 days)
        
        # Calculate active cases (like MATLAB)
        cumulative_shifted = [zeros(Int, recovery_days); cumulative_cases_from_start[1:end-recovery_days]]
        active_cases_count = cumulative_cases_from_start - cumulative_shifted
        # Ensure no negative active cases due to data corrections
        active_cases_count = max.(active_cases_count, 0)
        
        carson_days = [Dates.value(Day(d - dates_from_start[1])) for d in dates_from_start]
        
        # Interpolate active cases to simulation time grid
        itp_real_I = linear_interpolation(carson_days, active_cases_count, extrapolation_bc=0)
        real_interp_I = itp_real_I.(t_grid)
        
        # --- Load daily death data from the dedicated file ---
        # Note: For Washington, MS we'll try to load daily deaths data, but fall back gracefully if not available
        try
            daily_deaths_data = CSV.read("washington_mississippi_daily_deaths.csv", DataFrame)
            daily_deaths_data.date = Date.(daily_deaths_data.date)
            
            # Find start index for daily deaths data
            daily_start_idx = findfirst(daily_deaths_data.date .== start_date_real)
            if isnothing(daily_start_idx)
                date_diffs = abs.(daily_deaths_data.date .== start_date_real)
                daily_start_idx = argmin(date_diffs)
            end
            
            # Get daily deaths from March 25 onwards
            daily_deaths_from_start = daily_deaths_data.daily_deaths[daily_start_idx:end]
            daily_death_dates = daily_deaths_data.date[daily_start_idx:end]
            
            # Calculate days from March 25 for daily deaths
            daily_death_days = [Dates.value(Day(d - start_date_real)) for d in daily_death_dates]
            
            # Interpolate daily deaths to simulation time grid
            itp_real_D = linear_interpolation(daily_death_days, daily_deaths_from_start, extrapolation_bc=0)
            real_interp_D = itp_real_D.(t_grid)
            
            # Use the 7-day moving average from the CSV file for active deaths
            moving_avg_deaths = daily_deaths_data.moving_avg_7day[daily_start_idx:end]
            itp_moving_avg = linear_interpolation(daily_death_days, moving_avg_deaths, extrapolation_bc=0)
            real_active_D = itp_moving_avg.(t_grid)
            real_interp_D_prop = real_active_D / population
            
            # Also get cumulative deaths for the cumulative death plot
            cumulative_deaths_from_start = daily_deaths_data.cumulative_deaths[daily_start_idx:end]
            itp_cumulative = linear_interpolation(daily_death_days, cumulative_deaths_from_start, extrapolation_bc=0)
            real_cumulative_D = itp_cumulative.(t_grid)
        catch e
            @warn "Daily deaths file not available for Washington, MS. Using cumulative data from main file: $(e)"
            # Calculate daily deaths from cumulative data
            cumulative_deaths_from_start = data_table.deaths[start_idx:end]
            daily_deaths_from_start = Float64.([0; diff(cumulative_deaths_from_start)])
            daily_deaths_from_start = max.(daily_deaths_from_start, 0.0)  # Ensure non-negative
            
            # Create 7-day moving average
            window_size = 7
            moving_avg_deaths = Float64.(copy(daily_deaths_from_start))  # Ensure Float64 type
            for i in window_size:length(daily_deaths_from_start)
                moving_avg_deaths[i] = mean(daily_deaths_from_start[max(1, i-window_size+1):i])
            end
            
            # Interpolate to simulation time grid
            carson_days = [Dates.value(Day(d - dates_from_start[1])) for d in dates_from_start]
            itp_moving_avg = linear_interpolation(carson_days, moving_avg_deaths, extrapolation_bc=0)
            real_active_D = itp_moving_avg.(t_grid)
            real_interp_D_prop = real_active_D / population
            
            real_cumulative_D = cumulative_deaths_from_start[end] * ones(length(t_grid))
            if length(cumulative_deaths_from_start) >= length(t_grid)
                itp_cumulative = linear_interpolation(carson_days, cumulative_deaths_from_start, extrapolation_bc=0)
                real_cumulative_D = itp_cumulative.(t_grid)
            end
        end
        
        @printf("Successfully loaded real data with %d data points\n", length(active_cases_count))
        
    catch e
        @printf("Warning: Could not load or process real data: %s\n", e)
        real_interp_I = zeros(length(t_grid))
        real_interp_D_prop = zeros(length(t_grid))
    end
    
    # Convert all counts to proportions for plotting (optimized)
    population_factor = 1.0 / N  # Calculate division factor once
    
    # Convert all arrays to proportions using broadcasting (optimized)
    all_interp_I_prop, mean_I_prop, lower_I_prop, upper_I_prop,
    all_active_D_prop, mean_D_prop, lower_D_prop, upper_D_prop,
    real_interp_I_prop = (all_interp_I, mean_I, lower_I, upper_I,
                         all_active_D, mean_D, lower_D, upper_D,
                         real_interp_I) .* population_factor
    
    # Ensure all arrays are Float64 to avoid type issues
    all_interp_I_prop = Float64.(all_interp_I_prop)
    mean_I_prop = Float64.(mean_I_prop)
    lower_I_prop = Float64.(lower_I_prop)
    upper_I_prop = Float64.(upper_I_prop)
    all_active_D_prop = Float64.(all_active_D_prop)
    mean_D_prop = Float64.(mean_D_prop)
    lower_D_prop = Float64.(lower_D_prop)
    upper_D_prop = Float64.(upper_D_prop)
    real_interp_I_prop = Float64.(real_interp_I_prop)
    real_interp_D_prop = Float64.(real_interp_D_prop)

    # Set up date ticks once for all plots (optimized)
    xtick_positions, date_labels = setup_date_ticks(params, start_date_real, 90)

    # --- 5. Create the final plot with the uncertainty envelope for infected cases ---
    p1 = plot(t_grid, upper_I_prop, fillrange=lower_I_prop, 
              fillalpha=0.5, fillcolor=:lightblue, linealpha=0,
              label="90% Prediction Interval", linewidth=0)
    plot!(t_grid, real_interp_I_prop, color=:red, linewidth=2.5, 
          label="Real Data")
    xlabel!("Time (days)")
    ylabel!("Infected Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params.tmax)
    ylims!(0, maximum([upper_I_prop; real_interp_I_prop]) * 1.1)
    
    # Apply date ticks (reusing calculated values)
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p1, "SIHRS_Washington_MS_Full_Pandemic_bandwidth.png")
    
    # --- 6. Create a second figure with all stochastic sims and the real data (optimized) ---
    p2 = plot()
    
    # Plot all stochastic simulations as a matrix (much faster than loop)
    plot!(t_grid, all_interp_I_prop', color=RGBA(0.2, 0.4, 0.8, 0.3), 
          linewidth=1.0, label="")
    
    # Plot the real data as a solid red line
    plot!(t_grid, real_interp_I_prop, color=:red, linewidth=2.5, 
          label="Real Data")
    xlabel!("Time (days)")
    ylabel!("Infected Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params.tmax)
    ylims!(0, maximum([maximum(all_interp_I_prop); maximum(real_interp_I_prop)]) * 1.1)
    
    # Custom legend
    plot!([NaN], [NaN], color=RGBA(0.2, 0.4, 0.8), linewidth=2.5, 
          label="Stochastic Simulations")
    
    # Apply date ticks (reusing calculated values)
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p2, "SIHRS_Washington_MS_Full_Pandemic_trajectories.png")

    # --- 7. Create a figure for death proportion (active deaths) ---
    p3 = plot(t_grid, upper_D_prop, fillrange=lower_D_prop, 
              fillalpha=0.5, fillcolor=:lightgray, linealpha=0,
              label="90% Prediction Interval", linewidth=0)
    plot!(t_grid, real_interp_D_prop, color=:red, linewidth=2.5, 
          label="Real Data")
    xlabel!("Date")
    ylabel!("Active Death Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params.tmax)
    ylims!(0, maximum([upper_D_prop; real_interp_D_prop]) * 1.1)
    
    # Apply date ticks (reusing calculated values)
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p3, "SIHRS_Washington_MS_Full_Pandemic_active_deaths.png")

    # --- 8. Create a figure for cumulative death proportion ---
    # Use the actual cumulative deaths from the CSV file
    real_interp_D_prop_cumulative = Float64.(real_cumulative_D / population)
    p4 = plot(t_grid, real_interp_D_prop_cumulative, color=:red, linewidth=2.5, 
              label="Real Data")
    xlabel!("Date")
    ylabel!("Cumulative Death Proportion")
    title!("Washington, Mississippi")
    xlims!(0, params.tmax)
    ylims!(0, maximum(real_interp_D_prop_cumulative) * 1.1)
    
    # Apply date ticks (reusing calculated values)
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p4, "SIHRS_Washington_MS_Full_Pandemic_cumulative_deaths.png")
end

# =============================================================================
# EXECUTION
# =============================================================================

# Run the simulation
if abspath(PROGRAM_FILE) == @__FILE__
    sihrs_multiple_simulations_infected()
end 