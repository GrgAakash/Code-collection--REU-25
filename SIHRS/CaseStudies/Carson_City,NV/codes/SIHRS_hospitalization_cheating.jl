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

function interp1(x, y, x_new)
    itp = linear_interpolation(x, y, extrapolation_bc=0)
    return [itp(xi) for xi in x_new]
end

function setup_date_ticks(params::SIHRSParams, start_date_real, tick_interval=30)
    xtick_positions = collect(0:tick_interval:params.tmax)
    xtick_dates = [start_date_real + Day(x) for x in xtick_positions]
    date_labels = [Dates.format(d, "mm/dd") * "/" * string(year(d))[3:4] for d in xtick_dates]
    return xtick_positions, date_labels
end

function compute_quantiles(data, q_lower=0.05, q_upper=0.95)
    lower = mapslices(x -> quantile(x, q_lower), data, dims=1)[1, :]
    upper = mapslices(x -> quantile(x, q_upper), data, dims=1)[1, :]
    return lower, upper
end

# =============================================================================
# MAIN SIMULATION FUNCTION
# =============================================================================

function sihrs_hospitalization_simulation()
    # Initialize variables
    N = 56000  # Carson City, Nevada population
    
    # August 2, 2020 initial conditions (from our calculations)
    s0 = 0.994162
    i0 = 0.001214
    h0 = 0.000141
    r0 = 0.004357
    d0 = 0.000125
    
    # Define model parameters (using existing SIHRS model parameters)
    params = SIHRSParams(
        0.163,      # beta: infection rate (β > 0)
        0.126,      # gamma: I transition rate (γ > 0)
        0.1,        # alpha: H transition rate (α > 0)
        0.0083,     # lambda: R transition rate (Λ > 0)
        1.0,        # pSI: probability of S to I
        0.00,       # pII: probability of I to I (stay infected)
        0.04,       # pIH: probability of I to H
        0.94,       # pIR: probability of I to R
        0.020,      # pID: probability of I to D
        0.01,       # pHH: probability of H to H (stay hospitalized)
        0.9882,     # pHR: probability of H to R
        0.0018,     # pHD: probability of H to D
        0.02,       # pRR: probability of R to R (stay recovered)
        0.98,       # pRS: probability of R to S
        240,        # tmax: simulation end time (8 months from August 2)
        s0,         # s0: initial susceptible proportion
        i0,         # i0: initial infected proportion
        h0,         # h0: initial hospitalized proportion
        r0,         # r0: initial recovered proportion
        d0          # d0: initial dead proportion
    )
    
    # Verify R0 calculation
    calculated_R0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII)
    @printf("Calculated R0 = %.6f \n", calculated_R0)
    
    # Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-5
        error("Initial conditions must sum to 1")
    end
    
    num_simulations = 55  # Matching existing model
    
    # Store results for all simulations
    all_results = Vector{Any}(undef, num_simulations)
    
    # Run multiple simulations
    try
        @printf("Running %d stochastic simulations for N = %d...\n", num_simulations, N)
        @printf("Starting from August 2, 2020 with initial conditions:\n")
        @printf("  S: %.1f, I: %.1f, H: %.1f, R: %.1f, D: %.1f\n", 
                params.s0*N, params.i0*N, params.h0*N, params.r0*N, params.d0*N)
        
        for sim_idx in 1:num_simulations
            @printf("Running simulation %d/%d...\n", sim_idx, num_simulations)
            all_results[sim_idx] = sihrs_agent_model_hospitalization(N, params)
        end
        
        @printf("All simulations completed!\n")
        
        # Plot results and compare with real data
        plot_hospitalization_results(all_results, N, params)
        
    catch e
        @printf("Error occurred: %s\n", e)
        rethrow(e)
    end
end

function sihrs_agent_model_hospitalization(N::Int, params::SIHRSParams)
    # SIHRS agent-based stochastic model focusing on hospitalization
    
    # Initial conditions
    s0 = round(Int, params.s0 * N)
    i0 = round(Int, params.i0 * N)
    h0 = round(Int, params.h0 * N)
    r0 = round(Int, params.r0 * N)
    d0 = round(Int, params.d0 * N)
    
    # Adjust for rounding errors
    total = s0 + i0 + h0 + r0 + d0
    if total != N
        compartments = [s0, i0, h0, r0, d0]
        largest_idx = argmax(compartments)
        compartments[largest_idx] += (N - total)
        s0, i0, h0, r0, d0 = compartments
    end
    
    # Preallocate arrays
    max_events = N * 10
    T = zeros(max_events)
    H_prop = zeros(max_events)
    H_count = zeros(Int, max_events)
    
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
    H_prop[1] = h0 / total_pop
    H_count[1] = h0
    
    # Main simulation loop
    while !isempty(I) && t < params.tmax
        nS = length(S)
        nI = length(I)
        nH = length(H)
        nR = length(R)
        
        # Calculate event rates
        infection_rate = params.pSI * params.beta * nS * nI / N
        to_susceptible_from_R_rate = params.pRS * params.lambda * nR
        to_hospital_rate = params.gamma * nI * params.pIH
        to_recovered_from_I_rate = params.gamma * nI * params.pIR
        to_dead_from_I_rate = params.gamma * nI * params.pID
        to_recovered_from_H_rate = params.alpha * nH * params.pHR
        to_dead_from_H_rate = params.alpha * nH * params.pHD
        
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
            H_prop[event_count] = length(H) / current_total
            H_count[event_count] = length(H)
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
        H_prop[event_count] = length(H) / current_total
        H_count[event_count] = length(H)
    end
    
    # Trim unused preallocated space
    T = T[1:event_count]
    H_prop = H_prop[1:event_count]
    H_count = H_count[1:event_count]
    
    # Store results
    result = Dict(
        "N" => N,
        "T" => T,
        "H_count" => H_count,
        "H_prop" => H_prop,
        "final_time" => t,
        "peak_hospitalized" => maximum(H_count),
        "peak_time" => T[argmax(H_count)],
        "peak_hospitalized_prop" => maximum(H_prop),
        "peak_time_prop" => T[argmax(H_prop)],
        "h_final" => H_prop[end]
    )
    
    return result
end

function plot_hospitalization_results(all_results, N, params::SIHRSParams)
    # Create a common time grid for comparison
    t_grid = collect(0:params.tmax)
    
    # Interpolate each simulation's results onto the grid
    all_interp_H = zeros(length(all_results), length(t_grid))
    
    for i in 1:length(all_results)
        res = all_results[i]
        if length(res["T"]) > 1
            itp_H = linear_interpolation(res["T"], res["H_count"], extrapolation_bc=0)
            all_interp_H[i, :] = itp_H.(t_grid)
        end
    end
    
    # Calculate statistics
    mean_H = mean(all_interp_H, dims=1)[1, :]
    lower_H, upper_H = compute_quantiles(all_interp_H)
    
    # Load real hospitalization data
    real_interp_H = zeros(length(t_grid))
    start_date_real = Date("2020-08-02")
    
    try
        # Load hospitalization data
        hosp_data = CSV.read("hospitalization_Carson_filtered_new.csv", DataFrame)
        
        # Convert collection_week to Date
        hosp_data.collection_week = Date.(hosp_data.collection_week, "m/d/yy")
        
        # Get total hospitalization (adult + pediatric)
        hosp_data.total_hospitalized = hosp_data.total_adult_and_pediatric_covid_patients
        
        # Filter for Carson City (FIPS 32510) and valid data
        carson_data = hosp_data[hosp_data.fips_code .== 32510, :]
        carson_data = carson_data[carson_data.total_hospitalized .> 0, :]
        
        # Calculate days from August 2, 2020
        carson_days = [Dates.value(Day(d - start_date_real)) for d in carson_data.collection_week]
        
        # Interpolate to simulation time grid
        if length(carson_days) > 1
            itp_real_H = linear_interpolation(carson_days, carson_data.total_hospitalized, extrapolation_bc=0)
            real_interp_H = itp_real_H.(t_grid)
        end
        
        @printf("Successfully loaded real hospitalization data with %d data points\n", length(carson_data.total_hospitalized))
        
    catch e
        @printf("Warning: Could not load or process real hospitalization data: %s\n", e)
        real_interp_H = zeros(length(t_grid))
    end
    
    # Convert to proportions
    all_interp_H_prop = all_interp_H / N
    mean_H_prop = mean_H / N
    lower_H_prop = lower_H / N
    upper_H_prop = upper_H / N
    real_interp_H_prop = real_interp_H / N
    
    # Set up date ticks
    xtick_positions, date_labels = setup_date_ticks(params, start_date_real, 30)
    
    # Create the uncertainty envelope plot
    p1 = plot(t_grid, upper_H_prop, fillrange=lower_H_prop, 
              fillalpha=0.5, fillcolor=:lightblue, linealpha=0,
              label="90% Prediction Interval", linewidth=0)
    plot!(t_grid, real_interp_H_prop, color=:red, linewidth=2.5, 
          label="Real Data")
    plot!(t_grid, mean_H_prop, color=:blue, linewidth=2, 
          label="Mean Simulation")
    xlabel!("Time (days)")
    ylabel!("Hospitalized Proportion")
    title!("Carson City Hospitalization: August 2, 2020 Start")
    xlims!(0, params.tmax)
    ylims!(0, maximum([upper_H_prop; real_interp_H_prop]) * 1.1)
    
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p1, "SIHRS_Carson_City_Hospitalization_August2_bandwidth.png")
    
    # Create trajectories plot
    p2 = plot()
    
    # Plot all stochastic simulations
    plot!(t_grid, all_interp_H_prop', color=RGBA(0.2, 0.4, 0.8, 0.3), 
          linewidth=1.0, label="")
    
    # Plot the real data
    plot!(t_grid, real_interp_H_prop, color=:red, linewidth=2.5, 
          label="Real Data")
    plot!(t_grid, mean_H_prop, color=:blue, linewidth=2, 
          label="Mean Simulation")
    
    xlabel!("Time (days)")
    ylabel!("Hospitalized Proportion")
    title!("Carson City Hospitalization: August 2, 2020 Start")
    xlims!(0, params.tmax)
    ylims!(0, maximum([maximum(all_interp_H_prop); maximum(real_interp_H_prop)]) * 1.1)
    
    # Custom legend
    plot!([NaN], [NaN], color=RGBA(0.2, 0.4, 0.8), linewidth=2.5, 
          label="Stochastic Simulations")
    
    xticks!(xtick_positions, date_labels)
    xlabel!("Date (mm/dd/yy)")
    
    savefig(p2, "SIHRS_Carson_City_Hospitalization_August2_trajectories.png")
    
    # Print summary statistics
    println("\n=== Simulation Summary ===")
    println("Peak hospitalized (mean): $(round(maximum(mean_H), digits=1)) people")
    println("Peak hospitalized (real): $(round(maximum(real_interp_H), digits=1)) people")
    println("Final hospitalized (mean): $(round(mean_H[end], digits=1)) people")
    println("Final hospitalized (real): $(round(real_interp_H[end], digits=1)) people")
end

# =============================================================================
# EXECUTION
# =============================================================================

# Run the simulation
if abspath(PROGRAM_FILE) == @__FILE__
    sihrs_hospitalization_simulation()
end 