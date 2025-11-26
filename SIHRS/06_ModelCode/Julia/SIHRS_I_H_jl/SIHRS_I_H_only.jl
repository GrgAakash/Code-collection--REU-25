using DifferentialEquations
using Plots
using Random
using Statistics
using LinearAlgebra
using Printf
using Distributions

"""
SIHRS Epidemic Model - I and H Compartments Only

This script focuses specifically on plotting the Infected (I) and Hospitalized (H) 
compartments from the SIHRS epidemic model with death.
"""

# Define the parameter structure
struct SIHRSParams
    β::Float64      # infection rate (β > 0)
    γ::Float64      # I transition rate (γ > 0)
    α::Float64      # H transition rate (α > 0)
    λ::Float64      # R transition rate (λ > 0) - immunity period of 4 months
    pSI::Float64    # probability of S to I (p_{SI} in (0,1])
    pII::Float64    # probability of I to I (stay infected)
    pIH::Float64    # probability of I to H
    pIR::Float64    # probability of I to R
    pID::Float64    # probability of I to D
    pHH::Float64    # probability of H to H (stay hospitalized)
    pHR::Float64    # probability of H to R
    pHD::Float64    # probability of H to D
    pRR::Float64    # probability of R to R (stay recovered)
    pRS::Float64    # probability of R to S
    tmax::Float64   # simulation end time
    s0::Float64     # initial susceptible proportion
    i0::Float64     # initial infected proportion
    h0::Float64     # initial hospitalized proportion
    r0::Float64     # initial recovered proportion
    d0::Float64     # initial dead proportion
end

# Default parameters matching the original SIHRS model
function default_params()
    return SIHRSParams(
        0.212,   # β
        0.10,    # γ
        0.09,    # α
        0.0083,  # λ
        1.0,     # pSI
        0.0,     # pII
        0.035,   # pIH
        0.945,   # pIR
        0.02,    # pID
        0.0,     # pHH
        0.92,    # pHR
        0.08,    # pHD
        0.02,    # pRR
        0.98,    # pRS
        1000.0,  # tmax
        0.96,    # s0
        0.04,    # i0
        0.0,     # h0
        0.0,     # r0
        0.0      # d0
    )
end

"""
Validate SIHRS parameters
"""
function validate_parameters(params::SIHRSParams)
    # Validate rates are positive
    if any([params.β, params.γ, params.α, params.λ] .<= 0)
        error("All rates (β, γ, α, λ) must be positive")
    end
    
    # Validate probabilities are in [0,1]
    probs = [params.pSI, params.pII, params.pIH, params.pIR, params.pID,
             params.pHH, params.pHR, params.pHD, params.pRR, params.pRS]
    if any(probs .< 0 .|| probs .> 1)
        error("All probabilities must be in [0,1]")
    end
    
    # Validate probability sums
    if abs((params.pII + params.pIH + params.pIR + params.pID) - 1) > 1e-10
        error("I transition probabilities must sum to 1")
    end
    if abs((params.pHH + params.pHR + params.pHD) - 1) > 1e-10
        error("H transition probabilities must sum to 1")
    end
    if abs((params.pRR + params.pRS) - 1) > 1e-10
        error("R transition probabilities must sum to 1")
    end
end

"""
Validate initial conditions sum to 1
"""
function validate_initial_conditions(params::SIHRSParams)
    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1) > 1e-10
        error("Initial conditions must sum to 1")
    end
end

"""
SIHRS agent-based stochastic model with death
"""
function sihrs_agent_model(N::Int, params::SIHRSParams)
    # Input validation
    if N <= 0
        error("Population size must be positive")
    end
    
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
    
    # Initialize agent arrays
    S = collect(1:s0)           # susceptible agents
    I = collect((s0+1):(s0+i0)) # infected agents
    H = collect((s0+i0+1):(s0+i0+h0)) # hospitalized agents
    R = collect((s0+i0+h0+1):(s0+i0+h0+r0)) # recovered agents
    D = collect((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0)) # dead agents
    
    # Initialize tracking arrays
    T = Float64[0.0]           # time points
    S_prop = Float64[params.s0] # susceptible proportions
    I_prop = Float64[params.i0] # infected proportions
    H_prop = Float64[params.h0] # hospitalized proportions
    R_prop = Float64[params.r0] # recovered proportions
    D_prop = Float64[params.d0] # dead proportions
    I_count = Int[i0]           # infected counts
    H_count = Int[h0]           # hospitalized counts
    
    # Initialize time tracking
    t = 0.0
    event_count = 1
    
    # Main simulation loop
    while !isempty(I) && t < params.tmax
        nS = length(S)
        nI = length(I)
        nH = length(H)
        nR = length(R)
        
        # Calculate event rates according to the mathematical model
        infection_rate = params.pSI * params.β * nS * nI / N  # S to I rate
        to_susceptible_from_R_rate = params.pRS * params.λ * nR  # R to S rate
        to_hospital_rate = params.γ * nI * params.pIH  # I to H rate
        to_recovered_from_I_rate = params.γ * nI * params.pIR  # I to R rate
        to_dead_from_I_rate = params.γ * nI * params.pID  # I to D rate
        to_recovered_from_H_rate = params.α * nH * params.pHR  # H to R rate
        to_dead_from_H_rate = params.α * nH * params.pHD  # H to D rate
        
        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate +
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate +
                     to_dead_from_H_rate
        
        if total_rate == 0
            break
        end
        
        # Time of next event
        dt = rand(Exponential(1 / total_rate))
        t += dt
        
        if t > params.tmax
            t = params.tmax
            push!(T, t)
            current_total = length(S) + length(I) + length(H) + length(R) + length(D)
            push!(S_prop, length(S) / current_total)
            push!(I_prop, length(I) / current_total)
            push!(H_prop, length(H) / current_total)
            push!(R_prop, length(R) / current_total)
            push!(D_prop, length(D) / current_total)
            push!(I_count, length(I))
            push!(H_count, length(H))
            break
        end
        
        push!(T, t)
        
        # Determine which event occurs
        chance = rand() * total_rate
        if chance < infection_rate
            # S to I transition
            if !isempty(S)
                idx = rand(1:length(S))
                infected_agent = S[idx]
                deleteat!(S, idx)
                push!(I, infected_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate)
            # R to S transition
            if !isempty(R)
                idx = rand(1:length(R))
                susceptible_agent = R[idx]
                deleteat!(R, idx)
                push!(S, susceptible_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate)
            # I to H transition
            if !isempty(I)
                idx = rand(1:length(I))
                hospitalized_agent = I[idx]
                deleteat!(I, idx)
                push!(H, hospitalized_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate)
            # I to R transition
            if !isempty(I)
                idx = rand(1:length(I))
                recovered_agent = I[idx]
                deleteat!(I, idx)
                push!(R, recovered_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate)
            # I to D transition
            if !isempty(I)
                idx = rand(1:length(I))
                dead_agent = I[idx]
                deleteat!(I, idx)
                push!(D, dead_agent)
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate)
            # H to R transition
            if !isempty(H)
                idx = rand(1:length(H))
                recovered_agent = H[idx]
                deleteat!(H, idx)
                push!(R, recovered_agent)
            end
        else
            # H to D transition
            if !isempty(H)
                idx = rand(1:length(H))
                dead_agent = H[idx]
                deleteat!(H, idx)
                push!(D, dead_agent)
            end
        end
        
        # Update tracking arrays
        current_total = length(S) + length(I) + length(H) + length(R) + length(D)
        push!(S_prop, length(S) / current_total)
        push!(I_prop, length(I) / current_total)
        push!(H_prop, length(H) / current_total)
        push!(R_prop, length(R) / current_total)
        push!(D_prop, length(D) / current_total)
        push!(I_count, length(I))
        push!(H_count, length(H))
    end
    
    # Store results
    result = Dict{String, Any}()
    result["N"] = N
    result["T"] = T
    result["S_prop"] = S_prop
    result["I_prop"] = I_prop
    result["H_prop"] = H_prop
    result["R_prop"] = R_prop
    result["D_prop"] = D_prop
    result["I_count"] = I_count
    result["H_count"] = H_count
    result["final_time"] = t
    result["peak_infected"] = maximum(I_count)
    result["peak_hospitalized"] = maximum(H_count)
    result["peak_time"] = T[findfirst(x -> x == maximum(I_count), I_count)]
    result["peak_h_time"] = T[findfirst(x -> x == maximum(H_count), H_count)]
    
    # Calculate and store asymptotic values
    result["s_inf"] = S_prop[end]
    result["i_inf"] = I_prop[end]
    result["h_inf"] = H_prop[end]
    result["r_inf"] = R_prop[end]
    result["d_inf"] = D_prop[end]
    
    return result
end

"""
Solve the deterministic SIHRS model using ODE solver
"""
function solve_deterministic_sihrs(params::SIHRSParams)
    # Time span
    tspan = (0.0, params.tmax)
    
    # Initial conditions vector
    y0 = [params.s0, params.i0, params.h0, params.r0, params.d0]
    
    # Define the ODE system exactly as in the mathematical model
    function ode_system!(dy, y, p, t)
        dy[1] = -params.β * y[1] * y[2] * params.pSI + params.pRS * params.λ * y[4]  # ds/dt
        dy[2] = params.β * y[1] * y[2] * params.pSI - params.γ * (1 - params.pII) * y[2]  # di/dt
        dy[3] = params.pIH * params.γ * y[2] - params.α * (1 - params.pHH) * y[3]  # dh/dt
        dy[4] = params.pIR * params.γ * y[2] + params.pHR * params.α * y[3] - params.pRS * params.λ * y[4]  # dr/dt
        dy[5] = params.pID * params.γ * y[2] + params.pHD * params.α * y[3]  # dd/dt
    end
    
    # Solve the ODE system
    prob = ODEProblem(ode_system!, y0, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)
    
    # Verify conservation
    sum_y = [sum(u) for u in sol.u]
    if any(abs.(sum_y .- 1) .> 1e-6)
        @warn "Conservation of population not satisfied"
    end
    
    # Store results
    det_result = Dict{String, Any}()
    det_result["T"] = sol.t
    det_result["S_prop"] = [u[1] for u in sol.u]
    det_result["I_prop"] = [u[2] for u in sol.u]
    det_result["H_prop"] = [u[3] for u in sol.u]
    det_result["R_prop"] = [u[4] for u in sol.u]
    det_result["D_prop"] = [u[5] for u in sol.u]
    
    # Find peak infected and peak time
    peak_infected_prop, peak_idx = findmax(det_result["I_prop"])
    det_result["peak_infected_prop"] = peak_infected_prop
    det_result["peak_time"] = det_result["T"][peak_idx]
    
    # Find peak hospitalized and peak time
    peak_hospitalized_prop, peak_h_idx = findmax(det_result["H_prop"])
    det_result["peak_hospitalized_prop"] = peak_hospitalized_prop
    det_result["peak_h_time"] = det_result["T"][peak_h_idx]
    
    det_result["final_time"] = det_result["T"][end]
    
    # Calculate R0
    det_result["R0"] = params.pSI * params.β / (params.γ * (1 - params.pII))
    
    # Store asymptotic values
    det_result["s_inf"] = det_result["S_prop"][end]
    det_result["i_inf"] = det_result["I_prop"][end]
    det_result["h_inf"] = det_result["H_prop"][end]
    det_result["r_inf"] = det_result["R_prop"][end]
    det_result["d_inf"] = det_result["D_prop"][end]
    
    return det_result
end

"""
Create focused plots for I and H compartments only
"""
function plot_I_H_compartments(results, N_values, det_result, params)
    # Colors for different population sizes and deterministic
    colors = [:red, :orange, :darkred]  # Red tones for I and H
    det_color = :purple  # Purple for deterministic
    
    # 1. Create combined I and H ODE plot (deterministic only)
    p_ode_combined = plot(title="Deterministic I and H Compartments", 
                          xlabel="Time", ylabel="Proportion", 
                          grid=true, xlims=(0, params.tmax), size=(1000,700))
    
    # Plot deterministic I and H curves together
    plot!(p_ode_combined, det_result["T"], det_result["I_prop"], 
          color=:red, linewidth=3, label="Infected (I)")
    plot!(p_ode_combined, det_result["T"], det_result["H_prop"], 
          color=:lime, linewidth=3, label="Hospitalized (H)")
    
    # Add legend and save
    plot!(p_ode_combined, legend=:topright)
    savefig(p_ode_combined, "SIHRS_ODE_I_H_combined.png")
    
    # 2. Create individual plots for each population size N
    for i in 1:length(results)
        # Create figure for this N value with I and H separately
        p_n = plot(title="SIHRS Stochastic Model - N = $(N_values[i])", 
                   xlabel="Time (days)", ylabel="Proportion", 
                   grid=true, xlims=(0, params.tmax), size=(1000,700))
        
        # Plot I compartment for this N
        plot!(p_n, results[i]["T"], results[i]["I_prop"], 
              color=:red, linewidth=2, label="Infected (I)")
        
        # Plot H compartment for this N  
        plot!(p_n, results[i]["T"], results[i]["H_prop"], 
              color=:lime, linewidth=2, label="Hospitalized (H)")
        
        # Customize the plot
        plot!(p_n, legend=:topright)
        
        # Add parameter annotations
        param_text = "R₀=$(round(det_result["R0"], digits=2)), β=$(params.β), γ=$(params.γ), α=$(params.α), λ=$(params.λ)"
        annotate!(p_n, 0.02, 0.98, text(param_text, 10, :top, :left))
        
        # Save the figure
        savefig(p_n, "SIHRS_stochastic_N$(N_values[i])_I_H.png")
    end
    
    return p_ode_combined
end

"""
Main function to run SIHRS I and H focused analysis
"""
function sihrs_I_H_analysis()
    # Define model parameters
    params = default_params()
    
    # Validate parameters
    validate_parameters(params)
    
    # Validate initial conditions sum to 1
    validate_initial_conditions(params)
    
    # Population sizes to test
    N_values = [316, 3162, 10000]
    
    # Input validation
    if any(N_values .<= 0)
        error("Population sizes must be positive integers")
    end
    
    # Store results for comparison
    results = Vector{Dict{String, Any}}(undef, length(N_values))
    deterministic_result = nothing
    
    # Run simulation for each population size
    try
        for idx in 1:length(N_values)
            println("Running simulation for N = $(N_values[idx])...")
            results[idx] = sihrs_agent_model(N_values[idx], params)
            println("Completed N = $(N_values[idx])")
        end
        
        # Solve deterministic model
        println("Solving deterministic model...")
        deterministic_result = solve_deterministic_sihrs(params)
        println("Deterministic model completed")
        
        # Plot I and H focused analysis
        println("Generating I and H focused plots...")
        plot_I_H_compartments(results, N_values, deterministic_result, params)
        println("Plots generated and saved")
        
    catch e
        println("Error occurred: $(e)")
        rethrow(e)
    end
    
    # Print focused summary statistics
    println("\n=== I AND H COMPARTMENT ANALYSIS ===")
    println("Population Size | Peak I | Peak H | Peak I Time | Peak H Time | I(∞) | H(∞)")
    println("----------------|---------|---------|-------------|-------------|-------|-------")
    for i in 1:length(results)
        println(@sprintf("%15d | %7d | %7d | %11.2f | %11.2f | %5.3f | %5.3f",
            results[i]["N"], results[i]["peak_infected"], results[i]["peak_hospitalized"],
            results[i]["peak_time"], results[i]["peak_h_time"],
            results[i]["i_inf"], results[i]["h_inf"]))
    end
    
    println(@sprintf("%15s | %7.4f | %7.4f | %11.2f | %11.2f | %5.3f | %5.3f",
        "Deterministic", deterministic_result["peak_infected_prop"], 
        deterministic_result["peak_hospitalized_prop"],
        deterministic_result["peak_time"], deterministic_result["peak_h_time"],
        deterministic_result["i_inf"], deterministic_result["h_inf"]))
    
    # Print key metrics
    println("\n=== KEY METRICS ===")
    println("R₀ = $(round(deterministic_result["R0"], digits=4))")
    println("Peak I proportion = $(round(deterministic_result["peak_infected_prop"], digits=4))")
    println("Peak H proportion = $(round(deterministic_result["peak_hospitalized_prop"], digits=4))")
    println("I peak time = $(round(deterministic_result["peak_time"], digits=2)) days")
    println("H peak time = $(round(deterministic_result["peak_h_time"], digits=2)) days")
    
    return results, deterministic_result
end

# Run the analysis if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    sihrs_I_H_analysis()
end
