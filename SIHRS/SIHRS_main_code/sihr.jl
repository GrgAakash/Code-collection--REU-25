using DifferentialEquations
using Statistics
using Plots
using Random
using Distributions
using Printf

"""
SIHR Model - Julia Implementation
Converted from MATLAB sihr.m
This implements a Susceptible-Infected-Hospitalized-Recovered model with agent-based stochastic simulation
"""

# Define model parameters structure
struct SIHRParams
    beta::Float64      # infection rate (β > 0)
    gamma1::Float64    # I to H rate (γ₁ > 0)
    gamma2::Float64    # I to R rate (γ₂ > 0)
    alpha::Float64     # H to R rate (α > 0)
    p1::Float64        # probability of infection (p₁ in (0,1])
    p2::Float64        # probability of leaving I (p₂ in (0,1])
    p3::Float64        # probability of leaving H (p₃ in (0,1])
    ph::Float64        # probability of I to H vs R (p_h in (0,1])
    tmax::Float64      # simulation end time
    s0::Float64        # initial susceptible proportion
    i0::Float64        # initial infected proportion
    h0::Float64        # initial hospitalized proportion
    r0::Float64        # initial recovered proportion
end

function sihr_multiple_populations()
    # Define model parameters
    params = SIHRParams(
        1.30,    # beta - infection rate
        1.0,     # gamma1 - I to H rate
        1.0,     # gamma2 - I to R rate
        1.00,    # alpha - H to R rate
        0.5,     # p1 - probability of infection
        0.5,     # p2 - probability of leaving I
        0.5,     # p3 - probability of leaving H
        0.50,    # ph - probability of I to H vs R
        50.0,    # tmax - simulation end time
        0.96,    # s0 - initial susceptible proportion
        0.04,    # i0 - initial infected proportion
        0.0,     # h0 - initial hospitalized proportion
        0.00     # r0 - initial recovered proportion
    )

    # Validate parameters
    validate_parameters(params)
    
    # Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0) - 1) > 1e-10
        error("Initial conditions must sum to 1")
    end

    # Population sizes to test
    N_values = [316, 3162, 10000]
  
    # Input validation
    if any(N_values .<= 0)
        error("Population sizes must be positive integers")
    end
    
    # Store results for comparison
    results = Vector{Any}(undef, length(N_values))
    
    # Run simulation for each population size
    try
        for idx in 1:length(N_values)
            println("Running simulation for N = $(N_values[idx])...")
            results[idx] = sihr_agent_model(N_values[idx], params)
            println("Completed N = $(N_values[idx])")
        end
        
        # Solve deterministic model
        deterministic_result = solve_deterministic_sihr(params)
        
        # Plot comparison
        plot_comparison(results, N_values, deterministic_result, params)
        
        # Plot ODE only
        plot_ode_only(deterministic_result, params)
        
    catch ME
        println("Error occurred: $(ME)")
        rethrow(ME)
    end
end

function validate_parameters(params::SIHRParams)
    # Validate rates are positive
    if any([params.beta, params.gamma1, params.gamma2, params.alpha] .<= 0)
        error("All rates (beta, gamma1, gamma2, alpha) must be positive")
    end
    
    # Validate probabilities are in (0,1]
    probs = [params.p1, params.p2, params.p3, params.ph]
    if any(probs .<= 0) || any(probs .> 1)
        error("All probabilities must be in (0,1]")
    end
end

function sihr_agent_model(N::Int, params::SIHRParams)
    # SIHR agent-based stochastic model
    if N <= 0
        error("Population size must be positive")
    end
    
    # Initial conditions - using params values and ensuring they sum to N
    s0 = round(Int, params.s0 * N) # susceptible
    i0 = round(Int, params.i0 * N) # infected
    h0 = round(Int, params.h0 * N) # hospitalized
    r0 = round(Int, params.r0 * N) # recovered
    
    # Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0
    if total != N
        # Add or subtract the difference from the largest compartment
        compartments = [s0, i0, h0, r0]
        largest_idx = argmax(compartments)
        if largest_idx == 1
            s0 = s0 + (N - total)
        elseif largest_idx == 2
            i0 = i0 + (N - total)
        elseif largest_idx == 3
            h0 = h0 + (N - total)
        else
            r0 = r0 + (N - total)
        end
    end
    
    # Validate initial conditions sum to N
    if (s0 + i0 + h0 + r0) != N
        error("Initial conditions must sum to N")
    end
    
    # Preallocate arrays for better performance
    max_events = N * 10 # Estimate maximum number of events
    T = zeros(max_events)
    S_prop = zeros(max_events)
    I_prop = zeros(max_events)
    H_prop = zeros(max_events)
    R_prop = zeros(max_events)
    I_count = zeros(Int, max_events)
    
    # Initialize agent arrays
    S = collect(1:s0)
    I = collect((s0+1):(s0+i0))
    H = collect((s0+i0+1):(s0+i0+h0))
    R = if r0 > 0
        collect((s0+i0+h0+1):(s0+i0+h0+r0))
    else
        Int[]
    end
    
    # Initialize time tracking
    t = 0.0
    T[1] = 0.0
    event_count = 1
    
    # Initialize proportion tracking
    total_pop = s0 + i0 + h0 + r0
    S_prop[1] = s0 / total_pop
    I_prop[1] = i0 / total_pop
    H_prop[1] = h0 / total_pop
    R_prop[1] = r0 / total_pop
    I_count[1] = i0
    
    # Main simulation loop
    while !isempty(I) && t < params.tmax
        nI = length(I)
        nS = length(S)
        nH = length(H)
        
        # Calculate event rates according to the mathematical model
        infection_rate = params.p1 * params.beta * nS * nI / N  # S to I rate
        to_hospital_rate = params.p2 * params.ph * params.gamma1 * nI  # I to H rate
        to_recovered_from_I_rate = params.p2 * (1-params.ph) * params.gamma2 * nI  # I to R rate
        to_recovered_from_H_rate = params.p3 * params.alpha * nH  # H to R rate
        
        total_rate = infection_rate + to_hospital_rate + to_recovered_from_I_rate + to_recovered_from_H_rate
        
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
            current_total = length(S) + length(I) + length(H) + length(R)
            S_prop[event_count] = length(S) / current_total
            I_prop[event_count] = length(I) / current_total
            H_prop[event_count] = length(H) / current_total
            R_prop[event_count] = length(R) / current_total
            I_count[event_count] = length(I)
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
        elseif chance < (infection_rate + to_hospital_rate)
            # I to H transition
            if nI > 0
                num = rand(1:nI)
                hospitalized_agent = I[num]
                deleteat!(I, num)
                push!(H, hospitalized_agent)
            end
        elseif chance < (infection_rate + to_hospital_rate + to_recovered_from_I_rate)
            # I to R transition
            if nI > 0
                num = rand(1:nI)
                recovered_agent = I[num]
                deleteat!(I, num)
                push!(R, recovered_agent)
            end
        else
            # H to R transition
            if nH > 0
                num = rand(1:nH)
                recovered_agent = H[num]
                deleteat!(H, num)
                push!(R, recovered_agent)
            end
        end
        
        # Update tracking arrays
        current_total = length(S) + length(I) + length(H) + length(R)
        S_prop[event_count] = length(S) / current_total
        I_prop[event_count] = length(I) / current_total
        H_prop[event_count] = length(H) / current_total
        R_prop[event_count] = length(R) / current_total
        I_count[event_count] = length(I)
    end
    
    # Trim unused preallocated space
    T = T[1:event_count]
    S_prop = S_prop[1:event_count]
    I_prop = I_prop[1:event_count]
    H_prop = H_prop[1:event_count]
    R_prop = R_prop[1:event_count]
    I_count = I_count[1:event_count]
    
    # Calculate peak infected and peak time
    peak_infected = maximum(I_count)
    peak_time = T[findfirst(x -> x == peak_infected, I_count)]
    
    # Store results
    return (
        N = N,
        T = T,
        S_prop = S_prop,
        I_prop = I_prop,
        H_prop = H_prop,
        R_prop = R_prop,
        I_count = I_count,
        final_time = t,
        peak_infected = peak_infected,
        peak_time = peak_time,
        s_inf = S_prop[end],
        i_inf = I_prop[end],
        h_inf = H_prop[end],
        r_inf = R_prop[end]
    )
end

function solve_deterministic_sihr(params::SIHRParams)
    # Solve the deterministic SIHR model using ODE solver
    
    # Time span
    tspan = (0.0, params.tmax)
    
    # Initial conditions vector
    y0 = [params.s0, params.i0, params.h0, params.r0]
    
    # Define the ODE system exactly as in the mathematical model
    function ode_system!(du, u, p, t)
        s, i, h, r = u
        beta, gamma1, gamma2, alpha, p1, p2, p3, ph = p
        
        du[1] = -p1 * beta * s * i                                                          # ds/dt
        du[2] = p1 * beta * s * i - p2 * (ph * gamma1 + (1-ph) * gamma2) * i              # di/dt
        du[3] = p2 * ph * gamma1 * i - p3 * alpha * h                                      # dh/dt
        du[4] = p2 * (1-ph) * gamma2 * i + p3 * alpha * h                                  # dr/dt
    end
    
    # Parameters for ODE
    p = [params.beta, params.gamma1, params.gamma2, params.alpha, 
         params.p1, params.p2, params.p3, params.ph]
    
    # Solve the ODE system
    prob = ODEProblem(ode_system!, y0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)
    
    # Extract results
    T = sol.t
    S_prop = sol[1, :]
    I_prop = sol[2, :]
    H_prop = sol[3, :]
    R_prop = sol[4, :]
    
    # Find peak infected and peak time
    peak_infected_prop, peak_idx = findmax(I_prop)
    peak_time = T[peak_idx]
    final_time = T[end]
    
    # Calculate R0
    R0 = params.p1 * params.beta / (params.p2 * (params.ph * params.gamma1 + (1-params.ph) * params.gamma2))
    
    # Store asymptotic values
    s_inf = S_prop[end]
    i_inf = I_prop[end]
    h_inf = H_prop[end]
    r_inf = R_prop[end]
    
    # Verify asymptotic behavior
    if i_inf > 1e-6
        @warn "i(∞) may not be approaching 0 as expected"
    end
    if h_inf > 1e-6
        @warn "h(∞) may not be approaching 0 as expected"
    end
    
    # Verify the asymptotic relation from the theorem
    s_inf_theoretical = params.s0 * exp(-R0 * (params.s0 + params.i0 - s_inf))
    if abs(s_inf - s_inf_theoretical) > 1e-6
        @warn "Asymptotic s(∞) relation may not be satisfied"
    end
    
    return (
        T = T,
        S_prop = S_prop,
        I_prop = I_prop,
        H_prop = H_prop,
        R_prop = R_prop,
        peak_infected_prop = peak_infected_prop,
        peak_time = peak_time,
        final_time = final_time,
        R0 = R0,
        s_inf = s_inf,
        i_inf = i_inf,
        h_inf = h_inf,
        r_inf = r_inf
    )
end

function plot_comparison(results, N_values, det_result, params::SIHRParams)
    # Create comparison plots including deterministic solution
    
    # Colors for different population sizes and deterministic
    colors = [:blue, :green, :red] # More distinctive colors
    det_color = :purple  # Purple for deterministic
    
    # Create main comparison plot
    p = plot(layout=(2, 2), size=(1600, 1200))
    
    # Plot 1: Susceptible Proportion Over Time
    plot!(p, subplot=1, xlabel="Time", ylabel="Proportion Susceptible", 
          title="Susceptible Proportion Over Time", grid=true, xlims=(0, params.tmax))
    plot_handles = []
    for i in 1:length(results)
        h = plot!(p, results[i].T, results[i].S_prop, subplot=1, 
                  color=colors[i], linewidth=1.5, label="N=$(N_values[i])")
        push!(plot_handles, h)
    end
    h_det = plot!(p, det_result.T, det_result.S_prop, subplot=1, 
                  linestyle=:dash, color=det_color, linewidth=2, label="Deterministic")
    push!(plot_handles, h_det)
    
    # Plot 2: Infected Proportion Over Time
    plot!(p, subplot=2, xlabel="Time", ylabel="Proportion Infected", 
          title="Infected Proportion Over Time", grid=true, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p, results[i].T, results[i].I_prop, subplot=2, 
              color=colors[i], linewidth=1.5, label="")
    end
    plot!(p, det_result.T, det_result.I_prop, subplot=2, 
          linestyle=:dash, color=det_color, linewidth=2, label="")
    
    # Plot 3: Hospitalized Proportion Over Time
    plot!(p, subplot=3, xlabel="Time", ylabel="Proportion Hospitalized", 
          title="Hospitalized Proportion Over Time", grid=true, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p, results[i].T, results[i].H_prop, subplot=3, 
              color=colors[i], linewidth=1.5, label="")
    end
    plot!(p, det_result.T, det_result.H_prop, subplot=3, 
          linestyle=:dash, color=det_color, linewidth=2, label="")
    
    # Plot 4: Recovered Proportion Over Time
    plot!(p, subplot=4, xlabel="Time", ylabel="Proportion Recovered", 
          title="Recovered Proportion Over Time", grid=true, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p, results[i].T, results[i].R_prop, subplot=4, 
              color=colors[i], linewidth=1.5, label="")
    end
    plot!(p, det_result.T, det_result.R_prop, subplot=4, 
          linestyle=:dash, color=det_color, linewidth=2, label="")
    
    # Add legend to first subplot
    plot!(p, subplot=1, legend=:best)
    
    # Save the figure
    savefig(p, "SIHR_simulation_results.png")
    display(p)
    
    # Print summary statistics
    println("\n=== SIMULATION SUMMARY ===")
    println("Population Size | Peak Infected | Peak Time | s(∞) | i(∞) | h(∞) | r(∞)")
    println("----------------|---------------|-----------|-------|-------|-------|-------")
    for i in 1:length(results)
        @printf("%15d | %13d | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f\n",
            results[i].N, results[i].peak_infected, results[i].peak_time,
            results[i].s_inf, results[i].i_inf, results[i].h_inf, results[i].r_inf)
    end
    @printf("%15s | %13.4f | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f\n",
        "Deterministic", det_result.peak_infected_prop, det_result.peak_time,
        det_result.s_inf, det_result.i_inf, det_result.h_inf, det_result.r_inf)
    
    # Print asymptotic analysis
    println("\n=== ASYMPTOTIC ANALYSIS ===")
    @printf("R₀ = %.4f\n", det_result.R0)
    s_inf_theoretical = params.s0 * exp(-det_result.R0 * (params.s0 + params.i0 - det_result.s_inf))
    @printf("Theoretical s(∞) = %.6f\n", s_inf_theoretical)
    @printf("Numerical s(∞) = %.6f\n", det_result.s_inf)
    @printf("i(∞) = %.6f (should be ≈ 0)\n", det_result.i_inf)
    @printf("h(∞) = %.6f (should be = 0)\n", det_result.h_inf)
    @printf("r(∞) = %.6f\n", det_result.r_inf)
end

function plot_ode_only(det_result, params::SIHRParams)
    # Create a separate plot showing only the ODE solution for SIHR
    
    # Color scheme for ODE plots
    ode_colors = [:blue, :orange, :yellow, :purple] # Blue, Orange, Yellow, Purple
    
    # Plot all four compartments on a single chart
    p_ode = plot(size=(1000, 700))
    plot!(p_ode, det_result.T, det_result.S_prop, color=ode_colors[1], linewidth=2.5, label="Susceptible (S)")
    plot!(p_ode, det_result.T, det_result.I_prop, color=ode_colors[2], linewidth=2.5, label="Infected (I)")
    plot!(p_ode, det_result.T, det_result.H_prop, color=ode_colors[3], linewidth=2.5, label="Hospitalized (H)")
    plot!(p_ode, det_result.T, det_result.R_prop, color=ode_colors[4], linewidth=2.5, label="Recovered (R)")
    
    # Customize the plot
    plot!(p_ode, xlabel="Time", ylabel="Proportion of Population", 
          title="SIHR Model - Deterministic ODE Solution", 
          grid=true, xlims=(0, params.tmax), ylims=(0, 1), legend=:right)
    
    # Save the figure
    savefig(p_ode, "SIHR_ODE_only.png")
    display(p_ode)
    
    # Print ODE-specific summary
    println("\n=== ODE SOLUTION SUMMARY ===")
    @printf("Peak Infected Proportion: %.4f at time %.2f\n", det_result.peak_infected_prop, det_result.peak_time)
    @printf("R₀: %.4f\n", det_result.R0)
    println("Asymptotic Values:")
    @printf("  s(∞) = %.6f\n", det_result.s_inf)
    @printf("  i(∞) = %.6f\n", det_result.i_inf)
    @printf("  h(∞) = %.6f\n", det_result.h_inf)
    @printf("  r(∞) = %.6f\n", det_result.r_inf)
end

# Run the simulation if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    sihr_multiple_populations()
end
