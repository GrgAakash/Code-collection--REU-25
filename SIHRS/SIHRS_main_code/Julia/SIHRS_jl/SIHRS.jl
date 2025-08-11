using DifferentialEquations
using Plots
using Random
using Statistics
using LinearAlgebra
using Printf
using Distributions

"""
SIHRS Epidemic Model with Death - Julia Implementation

This module provides both stochastic (agent-based) and deterministic (ODE) 
implementations of the SIHRS epidemic model with death.
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

# Default parameters matching the MATLAB version
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
    result["final_time"] = t
    result["peak_infected"] = maximum(I_count)
    result["peak_time"] = T[findfirst(x -> x == maximum(I_count), I_count)]
    
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
Create comparison plots including deterministic solution
"""
function plot_comparison(results, N_values, det_result, params)
    # Colors for different population sizes and deterministic
    colors = [:blue, :green, :red]  # More distinctive colors
    det_color = :purple  # Purple for deterministic
    
    # Create comparison plots
    p1 = plot(title="Susceptible Proportion Over Time", xlabel="Time", ylabel="Proportion Susceptible", 
              grid=true, legend=false, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p1, results[i]["T"], results[i]["S_prop"], color=colors[i], linewidth=1.5, 
              label="N=$(N_values[i])")
    end
    plot!(p1, det_result["T"], det_result["S_prop"], linestyle=:dash, color=det_color, 
          linewidth=2, label="Deterministic")
    
    p2 = plot(title="Infected Proportion Over Time", xlabel="Time", ylabel="Proportion Infected", 
              grid=true, legend=false, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p2, results[i]["T"], results[i]["I_prop"], color=colors[i], linewidth=1.5)
    end
    plot!(p2, det_result["T"], det_result["I_prop"], linestyle=:dash, color=det_color, linewidth=2)
    
    p3 = plot(title="Hospitalized Proportion Over Time", xlabel="Time", ylabel="Proportion Hospitalized", 
              grid=true, legend=false, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p3, results[i]["T"], results[i]["H_prop"], color=colors[i], linewidth=1.5)
    end
    plot!(p3, det_result["T"], det_result["H_prop"], linestyle=:dash, color=det_color, linewidth=2)
    
    p4 = plot(title="Recovered Proportion Over Time", xlabel="Time", ylabel="Proportion Recovered", 
              grid=true, legend=false, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p4, results[i]["T"], results[i]["R_prop"], color=colors[i], linewidth=1.5)
    end
    plot!(p4, det_result["T"], det_result["R_prop"], linestyle=:dash, color=det_color, linewidth=2)
    
    p5 = plot(title="Dead Proportion Over Time", xlabel="Time", ylabel="Proportion Dead", 
              grid=true, legend=false, xlims=(0, params.tmax))
    for i in 1:length(results)
        plot!(p5, results[i]["T"], results[i]["D_prop"], color=colors[i], linewidth=1.5)
    end
    plot!(p5, det_result["T"], det_result["D_prop"], linestyle=:dash, color=det_color, linewidth=2)
    
    # Create deterministic curves plot
    p6 = plot(title="Deterministic SIHRS with Death Model Dynamics", xlabel="Time", ylabel="Population Proportion", 
              grid=true, xlims=(0, params.tmax))
    plot!(p6, det_result["T"], det_result["S_prop"], linewidth=2, color=:blue, label="Susceptible")
    plot!(p6, det_result["T"], det_result["I_prop"], linewidth=2, color=:red, label="Infected")
    plot!(p6, det_result["T"], det_result["H_prop"], linewidth=2, color=:lime, label="Hospitalized")
    plot!(p6, det_result["T"], det_result["R_prop"], linewidth=2, color=:magenta, label="Recovered")
    plot!(p6, det_result["T"], det_result["D_prop"], linewidth=2, color=:black, label="Dead")
    
    # Combine all plots
    comparison_plot = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1200,900))
    
    # Save plots
    savefig(comparison_plot, "SIHRS_simulation_results.png")
    savefig(p6, "SIHRS_deterministic_only.png")
    
    # Create individual ODE plots for each compartment
    # Susceptible (S) only
    p_S_only = plot(title="Deterministic Susceptible (S) Compartment", 
                    xlabel="Time", ylabel="Proportion Susceptible", 
                    grid=true, xlims=(0, params.tmax), size=(800,600))
    plot!(p_S_only, det_result["T"], det_result["S_prop"], 
          color=:blue, linewidth=3, label="Susceptible")
    savefig(p_S_only, "SIHRS_ODE_S_only.png")
    
    # Infected (I) only
    p_I_only = plot(title="Deterministic Infected (I) Compartment", 
                    xlabel="Time", ylabel="Proportion Infected", 
                    grid=true, xlims=(0, params.tmax), size=(800,600))
    plot!(p_I_only, det_result["T"], det_result["I_prop"], 
          color=:red, linewidth=3, label="Infected")
    savefig(p_I_only, "SIHRS_ODE_I_only.png")
    
    # Hospitalized (H) only
    p_H_only = plot(title="Deterministic Hospitalized (H) Compartment", 
                    xlabel="Time", ylabel="Proportion Hospitalized", 
                    grid=true, xlims=(0, params.tmax), size=(800,600))
    plot!(p_H_only, det_result["T"], det_result["H_prop"], 
          color=:lime, linewidth=3, label="Hospitalized")
    savefig(p_H_only, "SIHRS_ODE_H_only.png")
    
    # Recovered (R) only
    p_R_only = plot(title="Deterministic Recovered (R) Compartment", 
                    xlabel="Time", ylabel="Proportion Recovered", 
                    grid=true, xlims=(0, params.tmax), size=(800,600))
    plot!(p_R_only, det_result["T"], det_result["R_prop"], 
          color=:magenta, linewidth=3, label="Recovered")
    savefig(p_R_only, "SIHRS_ODE_R_only.png")
    
    # Dead (D) only
    p_D_only = plot(title="Deterministic Dead (D) Compartment", 
                    xlabel="Time", ylabel="Proportion Dead", 
                    grid=true, xlims=(0, params.tmax), size=(800,600))
    plot!(p_D_only, det_result["T"], det_result["D_prop"], 
          color=:black, linewidth=3, label="Dead")
    savefig(p_D_only, "SIHRS_ODE_D_only.png")
    
    # Create separate figures for each population size N
    for i in 1:length(results)
        p_n = plot(title="SIHRS Stochastic Model - N = $(N_values[i])", xlabel="Time (days)", ylabel="Proportion", 
                   grid=true, xlims=(0, params.tmax), ylims=(0,1), size=(800,600))
        plot!(p_n, results[i]["T"], results[i]["S_prop"], color=:blue, linewidth=2, label="Susceptible (S)")
        plot!(p_n, results[i]["T"], results[i]["I_prop"], color=:red, linewidth=2, label="Infected (I)")
        plot!(p_n, results[i]["T"], results[i]["H_prop"], color=:lime, linewidth=2, label="Hospitalized (H)")
        plot!(p_n, results[i]["T"], results[i]["R_prop"], color=:magenta, linewidth=2, label="Recovered (R)")
        plot!(p_n, results[i]["T"], results[i]["D_prop"], color=:black, linewidth=2, label="Dead (D)")
        
        # Add parameter annotations
        param_text = "R₀=$(round(det_result["R0"], digits=2)), β=$(params.β), γ=$(params.γ), α=$(params.α), λ=$(params.λ)"
        annotate!(p_n, 0.02, 0.98, text(param_text, 10, :top, :left))
        
        savefig(p_n, "SIHRS_stochastic_N$(N_values[i]).png")
    end
    
    return comparison_plot
end

"""
Main function to run SIHRS multiple population simulations
"""
function sihrs_multiple_populations()
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
        
        # Plot comparison
        println("Generating plots...")
        plot_comparison(results, N_values, deterministic_result, params)
        println("Plots generated and saved")
        
    catch e
        println("Error occurred: $(e)")
        rethrow(e)
    end
    
    # Print summary statistics
    println("\n=== SIMULATION SUMMARY ===")
    println("Population Size | Peak Infected | Peak Time | s(∞) | i(∞) | h(∞) | r(∞) | d(∞)")
    println("----------------|---------------|-----------|-------|-------|-------|-------|-------")
    for i in 1:length(results)
        println(@sprintf("%15d | %13d | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f",
            results[i]["N"], results[i]["peak_infected"], results[i]["peak_time"],
            results[i]["s_inf"], results[i]["i_inf"], results[i]["h_inf"], 
            results[i]["r_inf"], results[i]["d_inf"]))
    end
    println(@sprintf("%15s | %13.4f | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f",
        "Deterministic", deterministic_result["peak_infected_prop"], deterministic_result["peak_time"],
        deterministic_result["s_inf"], deterministic_result["i_inf"], deterministic_result["h_inf"],
        deterministic_result["r_inf"], deterministic_result["d_inf"]))
    
    # Print asymptotic analysis
    println("\n=== ASYMPTOTIC ANALYSIS ===")
    println("R₀ = $(round(deterministic_result["R0"], digits=4))")
    println("s(∞) = $(round(deterministic_result["s_inf"], digits=6))")
    println("i(∞) = $(round(deterministic_result["i_inf"], digits=6))")
    println("h(∞) = $(round(deterministic_result["h_inf"], digits=6))")
    println("r(∞) = $(round(deterministic_result["r_inf"], digits=6))")
    println("d(∞) = $(round(deterministic_result["d_inf"], digits=6))")
    
    return results, deterministic_result
end

# Run the simulation if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    sihrs_multiple_populations()
end
