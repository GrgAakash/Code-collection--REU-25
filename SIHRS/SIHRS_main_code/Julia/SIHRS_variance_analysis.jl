using DifferentialEquations
using Statistics
using Plots
using Random
using Distributions

"""
SIHRS Model Variance Analysis in Julia
This function analyzes the stochastic scaling properties of the SIHRS model
and validates the variance scaling conjecture: V_N^(l)(t) ∝ (1/N) × population_proportions
"""
function SIHRS_variance_analysis()
    # Model Parameters (matching SIHRS.m)
    params = Dict(
        "beta" => 0.212,           # infection rate (β > 0)
        "gamma" => 0.10,           # I transition rate (γ > 0)
        "alpha" => 0.1,            # H transition rate (α > 0)
        "lambda" => 0.0083,        # R transition rate (Λ > 0) immunity period of 4 months
        "pSI" => 1.0,              # probability of S to I (p_{SI} in (0,1])
        "pII" => 0.0,              # probability of I to I (stay infected)
        "pIH" => 0.04,             # probability of I to H
        "pIR" => 0.959,            # probability of I to R
        "pID" => 0.001,            # probability of I to D
        "pHH" => 0.01,             # probability of H to H (stay hospitalized)
        "pHR" => 0.9882,           # probability of H to R
        "pHD" => 0.0018,           # probability of H to D
        "pRR" => 0.02,             # probability of R to R (stay recovered)
        "pRS" => 0.98,             # probability of R to S
        "tmax" => 1000,            # simulation end time
        "s0" => 0.96,              # initial susceptible proportion
        "i0" => 0.04,              # initial infected proportion
        "h0" => 0.0,               # initial hospitalized proportion
        "r0" => 0.0,               # initial recovered proportion
        "d0" => 0.0                # initial dead proportion
    )
    
    # Validate parameters
    validate_parameters(params)
    
    # Test different population sizes to see stochastic effects
    N_values = [316, 3162]
    num_trials = 15  # Number of trials per batch
    num_batches = 5  # Number of batches for each N
    
    # Time grid for analysis
    t_grid = range(0, params["tmax"], length=1000)
    
    # Pre-allocate storage for variance analysis (5 batches × 15 trials each)
    S_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches)
    I_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches)
    H_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches)
    R_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches)
    D_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches)
    
    # Pre-allocate storage for batch-averaged expected variance
    S_expected_var = zeros(length(N_values), length(t_grid), num_batches)
    I_expected_var = zeros(length(N_values), length(t_grid), num_batches)
    H_expected_var = zeros(length(N_values), length(t_grid), num_batches)
    R_expected_var = zeros(length(N_values), length(t_grid), num_batches)
    D_expected_var = zeros(length(N_values), length(t_grid), num_batches)
    
    # Cell arrays to store theoretical variance data
    v_s_theory = Vector{Vector{Float64}}(undef, length(N_values))
    v_i_theory = Vector{Vector{Float64}}(undef, length(N_values))
    v_h_theory = Vector{Vector{Float64}}(undef, length(N_values))
    v_r_theory = Vector{Vector{Float64}}(undef, length(N_values))
    v_d_theory = Vector{Vector{Float64}}(undef, length(N_values))
    
    # Store ODE solutions for time-averaged analysis
    s_ode_all = Vector{Vector{Float64}}(undef, length(N_values))
    i_ode_all = Vector{Vector{Float64}}(undef, length(N_values))
    h_ode_all = Vector{Vector{Float64}}(undef, length(N_values))
    r_ode_all = Vector{Vector{Float64}}(undef, length(N_values))
    d_ode_all = Vector{Vector{Float64}}(undef, length(N_values))
    
    # Main Simulation Loop
    for idx in 1:length(N_values)
        N = N_values[idx]
        println("Running variance analysis for N = $N...")
        
        # Run 5 batches of 15 trials each
        for batch in 1:num_batches
            println("  Batch $batch/$num_batches:")
            
            # Run 15 trials for this batch
            for trial in 1:num_trials
                println("    Trial $trial/$num_trials")
                
                # Run single stochastic simulation
                result = sihrs_agent_model(N, params)
                
                # Interpolate variance data onto common time grid
                S_var_all[idx, :, trial, batch] = interpolate_variance(result.T, result.vs, t_grid)
                I_var_all[idx, :, trial, batch] = interpolate_variance(result.T, result.vi, t_grid)
                H_var_all[idx, :, trial, batch] = interpolate_variance(result.T, result.vh, t_grid)
                R_var_all[idx, :, trial, batch] = interpolate_variance(result.T, result.vr, t_grid)
                D_var_all[idx, :, trial, batch] = interpolate_variance(result.T, result.vd, t_grid)
            end
            
            # Calculate expected variance for this batch (average across 15 trials)
            S_expected_var[idx, :, batch] = mean(S_var_all[idx, :, :, batch], dims=3)[:, 1, 1]
            I_expected_var[idx, :, batch] = mean(I_var_all[idx, :, :, batch], dims=3)[:, 1, 1]
            H_expected_var[idx, :, batch] = mean(H_var_all[idx, :, :, batch], dims=3)[:, 1, 1]
            R_expected_var[idx, :, batch] = mean(R_var_all[idx, :, :, batch], dims=3)[:, 1, 1]
            D_expected_var[idx, :, batch] = mean(D_var_all[idx, :, :, batch], dims=3)[:, 1, 1]
            
            println("    Batch $batch expected variance calculated")
        end
        
        # Theoretical Variance Calculation using ODE solution
        # Solve deterministic ODE for theoretical comparison
        det_result = solve_deterministic_sihrs(params)
        
        # Interpolate ODE solution to our time grid
        s_ode = interpolate_variance(det_result.T, det_result.S_prop, t_grid)
        i_ode = interpolate_variance(det_result.T, det_result.I_prop, t_grid)
        h_ode = interpolate_variance(det_result.T, det_result.H_prop, t_grid)
        r_ode = interpolate_variance(det_result.T, det_result.R_prop, t_grid)
        d_ode = interpolate_variance(det_result.T, det_result.D_prop, t_grid)
        
        # Store ODE solutions for time-averaged analysis
        s_ode_all[idx] = s_ode
        i_ode_all[idx] = i_ode
        h_ode_all[idx] = i_ode
        r_ode_all[idx] = r_ode
        d_ode_all[idx] = d_ode
        
        # Calculate theoretical variance using the mathematical formulas
        # V_N^(1)(t) = (1/N)(β s(t) i(t) p_SI + Λ r(t) p_RS)
        v_s_theory[idx] = (1/N) * (params["beta"] * s_ode .* i_ode * params["pSI"] + params["lambda"] * r_ode * params["pRS"])
        
        # V_N^(2)(t) = (1/N)(γ i(t) (p_IH + p_IR + p_ID) + β s(t) i(t) p_SI)
        v_i_theory[idx] = (1/N) * (params["gamma"] * i_ode .* (params["pIH"] + params["pIR"] + params["pID"]) + params["beta"] * s_ode .* i_ode * params["pSI"])
        
        # V_N^(3)(t) = (1/N)(γ i(t) p_IH + α h(t) (p_HR + p_HD))
        v_h_theory[idx] = (1/N) * (params["gamma"] * i_ode .* params["pIH"] + params["alpha"] * h_ode .* (params["pHR"] + params["pHD"]))
        
        # V_N^(4)(t) = (1/N)(γ i(t) p_IR + α h(t) p_HR + Λ r(t) p_RS)
        v_r_theory[idx] = (1/N) * (params["gamma"] * i_ode .* params["pIR"] + params["alpha"] * h_ode .* params["pHR"] + params["lambda"] * r_ode .* params["pRS"])
        
        # V_N^(5)(t) = (1/N)(γ i(t) p_ID + α h(t) p_HD)
        v_d_theory[idx] = (1/N) * (params["gamma"] * i_ode .* params["pID"] + params["alpha"] * h_ode .* params["pHD"])
    end
    
    # Time-Averaged Variance Analysis
    # Define time intervals for analysis [T1, T2]
    time_intervals = [[100, 200], [300, 400], [500, 600], [700, 800], [800, 900], [900, 1000]]
    
    println("\n=== TIME-AVERAGED VARIANCE ANALYSIS ===")
    for (interval_idx, interval) in enumerate(time_intervals)
        T1, T2 = interval[1], interval[2]
        
        println("\nTime interval [$T1, $T2]:")
        
        # Find indices for this time interval
        t_indices = findall(t -> t >= T1 && t <= T2, t_grid)
        
        for n_idx in 1:length(N_values)
            N = N_values[n_idx]
            println("  N = $N:")
            
            # Calculate time-averaged simulated variance
            V_bar_S = mean(mean(S_expected_var[n_idx, t_indices, :], dims=2))
            V_bar_I = mean(mean(I_expected_var[n_idx, t_indices, :], dims=2))
            V_bar_H = mean(mean(H_expected_var[n_idx, t_indices, :], dims=2))
            V_bar_R = mean(mean(R_expected_var[n_idx, t_indices, :], dims=2))
            V_bar_D = mean(mean(D_expected_var[n_idx, t_indices, :], dims=2))
            
            # Calculate time-averaged theoretical variance
            v_bar_S = mean(v_s_theory[n_idx][t_indices])
            v_bar_I = mean(v_i_theory[n_idx][t_indices])
            v_bar_H = mean(v_h_theory[n_idx][t_indices])
            v_bar_R = mean(v_r_theory[n_idx][t_indices])
            v_bar_D = mean(v_d_theory[n_idx][t_indices])
            
            # Calculate time-averaged population proportions
            s_avg = mean(s_ode_all[n_idx][t_indices])
            i_avg = mean(i_ode_all[n_idx][t_indices])
            h_avg = mean(h_ode_all[n_idx][t_indices])
            r_avg = mean(r_ode_all[n_idx][t_indices])
            d_avg = mean(d_ode_all[n_idx][t_indices])
            
            println("    S: V̄_sim=$V_bar_S, V̄_theory=$v_bar_S, s̄=$s_avg")
            println("    I: V̄_sim=$V_bar_I, V̄_theory=$v_bar_I, ī=$i_avg")
            println("    H: V̄_sim=$V_bar_H, V̄_theory=$v_bar_H, h̄=$h_avg")
            println("    R: V̄_sim=$V_bar_R, V̄_theory=$v_bar_R, r̄=$r_avg")
            println("    D: V̄_sim=$V_bar_D, V̄_theory=$v_bar_D, d̄=$d_avg")
        end
    end
    
    # Plotting Expected Variance Analysis
    plot_variance_error_bars(N_values, t_grid, S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var,
        v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory)
    
    println("\nExpected variance analysis completed successfully!")
    println("Following the new conjecture: E[V_N^(l)(t)] proportional to (1/N) x population_proportions")
end

"""
Helper function to validate parameters
"""
function validate_parameters(params)
    # Check that all rates are positive
    if any([params["beta"], params["gamma"], params["alpha"], params["lambda"]] .<= 0)
        error("All rates (beta, gamma, alpha, lambda) must be positive")
    end
    
    # Probabilities should be between 0 and 1
    probs = [params["pSI"], params["pII"], params["pIH"], params["pIR"], params["pID"],
             params["pHH"], params["pHR"], params["pHD"], params["pRR"], params["pRS"]]
    if any(probs .< 0) || any(probs .> 1)
        error("All probabilities must be in [0,1]")
    end
    
    # Make sure probability sums add up correctly
    if abs((params["pII"] + params["pIH"] + params["pIR"] + params["pID"]) - 1) > 1e-10
        error("I transition probabilities must sum to 1")
    end
    if abs((params["pHH"] + params["pHR"] + params["pHD"]) - 1) > 1e-10
        error("H transition probabilities must sum to 1")
    end
    if abs((params["pRR"] + params["pRS"]) - 1) > 1e-10
        error("R transition probabilities must sum to 1")
    end
end

"""
SIHRS agent-based stochastic model with variance tracking
"""
function sihrs_agent_model(N, params)
    # Initial conditions
    s0 = round(Int, params["s0"] * N)
    i0 = round(Int, params["i0"] * N)
    h0 = round(Int, params["h0"] * N)
    r0 = round(Int, params["r0"] * N)
    d0 = round(Int, params["d0"] * N)
    
    # Initialize arrays
    max_events = N * 30
    T = zeros(max_events)
    S_prop = zeros(max_events)
    I_prop = zeros(max_events)
    H_prop = zeros(max_events)
    R_prop = zeros(max_events)
    D_prop = zeros(max_events)
    
    # Variance tracking arrays
    vs = zeros(max_events)
    vi = zeros(max_events)
    vh = zeros(max_events)
    vr = zeros(max_events)
    vd = zeros(max_events)
    
    # Initialize
    S = collect(1:s0)
    I = collect((s0+1):(s0+i0))
    H = Int[]
    R = Int[]
    D = Int[]
    
    t = 0.0
    event_count = 1
    
    # Record initial state
    T[1] = 0
    S_prop[1] = s0 / N
    I_prop[1] = i0 / N
    H_prop[1] = 0
    R_prop[1] = 0
    D_prop[1] = 0
    
    # Calculate initial variance
    vs[1] = (params["beta"] * (s0/N) * (i0/N) * params["pSI"] + params["lambda"] * (r0/N) * params["pRS"]) / N
    vi[1] = (params["gamma"] * (i0/N) * (params["pIH"] + params["pIR"] + params["pID"]) + params["beta"] * (s0/N) * (i0/N) * params["pSI"]) / N
    vh[1] = (params["gamma"] * (i0/N) * params["pIH"] + params["alpha"] * 0 * (params["pHR"] + params["pHD"])) / N
    vr[1] = (params["gamma"] * (i0/N) * params["pIR"] + params["alpha"] * 0 * params["pHR"] + params["lambda"] * (r0/N) * params["pRS"]) / N
    vd[1] = (params["gamma"] * (i0/N) * params["pID"] + params["alpha"] * 0 * params["pHD"]) / N
    
    # Main simulation loop
    while t < params["tmax"] && event_count < max_events
        ns = length(S)
        ni = length(I)
        nh = length(H)
        nr = length(R)
        nd = length(D)
        
        # Calculate event rates
        infection_rate = params["beta"] * ns * ni / N * params["pSI"]
        i_to_h_rate = params["gamma"] * ni * params["pIH"]
        i_to_r_rate = params["gamma"] * ni * params["pIR"]
        i_to_d_rate = params["gamma"] * ni * params["pID"]
        h_to_r_rate = params["alpha"] * nh * params["pHR"]
        h_to_d_rate = params["alpha"] * nh * params["pHD"]
        r_to_s_rate = params["lambda"] * nr * params["pRS"]
        
        total_rate = infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate + h_to_d_rate + r_to_s_rate
        
        if total_rate == 0
            break
        end
        
        # Time to next event
        dt = rand(Exponential(1 / total_rate))
        t += dt
        
        if t > params["tmax"]
            t = params["tmax"]
        end
        
        event_count += 1
        T[event_count] = t
        
        # Determine which event occurs
        r = rand() * total_rate
        if r < infection_rate
            # S -> I
            if ns > 0
                push!(I, pop!(S))
            end
        elseif r < infection_rate + i_to_h_rate
            # I -> H
            if ni > 0
                push!(H, pop!(I))
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate
            # I -> R
            if ni > 0
                push!(R, pop!(I))
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate
            # I -> D
            if ni > 0
                push!(D, pop!(I))
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate
            # H -> R
            if nh > 0
                push!(R, pop!(H))
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate + h_to_d_rate
            # H -> D
            if nh > 0
                push!(D, pop!(H))
            end
        else
            # R -> S
            if nr > 0
                push!(S, pop!(R))
            end
        end
        
        # Update proportions
        ns = length(S)
        ni = length(I)
        nh = length(H)
        nr = length(R)
        nd = length(D)
        
        S_prop[event_count] = ns / N
        I_prop[event_count] = ni / N
        H_prop[event_count] = nh / N
        R_prop[event_count] = nr / N
        D_prop[event_count] = nd / N
        
        # Calculate variance at this time step
        vs[event_count] = (params["beta"] * (ns/N) * (ni/N) * params["pSI"] + params["lambda"] * (nr/N) * params["pRS"]) / N
        vi[event_count] = (params["gamma"] * (ni/N) * (params["pIH"] + params["pIR"] + params["pID"]) + params["beta"] * (ns/N) * (ni/N) * params["pSI"]) / N
        vh[event_count] = (params["gamma"] * (ni/N) * params["pIH"] + params["alpha"] * (nh/N) * (params["pHR"] + params["pHD"])) / N
        vr[event_count] = (params["gamma"] * (ni/N) * params["pIR"] + params["alpha"] * (nh/N) * params["pHR"] + params["lambda"] * (nr/N) * params["pRS"]) / N
        vd[event_count] = (params["gamma"] * (ni/N) * params["pID"] + params["alpha"] * (nh/N) * params["pHD"]) / N
    end
    
    # Trim arrays
    T = T[1:event_count]
    S_prop = S_prop[1:event_count]
    I_prop = I_prop[1:event_count]
    H_prop = H_prop[1:event_count]
    R_prop = R_prop[1:event_count]
    D_prop = D_prop[1:event_count]
    vs = vs[1:event_count]
    vi = vi[1:event_count]
    vh = vh[1:event_count]
    vr = vr[1:event_count]
    vd = vd[1:event_count]
    
    # Store results
    return (N=N, T=T, S_prop=S_prop, I_prop=I_prop, H_prop=H_prop, R_prop=R_prop, D_prop=D_prop,
            vs=vs, vi=vi, vh=vh, vr=vr, vd=vd, final_time=t)
end

"""
Solve the deterministic SIHRS model using ODE solver
"""
function solve_deterministic_sihrs(params)
    # Define the ODE system
    function ode_system!(du, u, p, t)
        s, i, h, r, d = u
        beta, gamma, alpha, lambda, pSI, pII, pIH, pIR, pID, pHH, pHR, pHD, pRR, pRS = p
        
        du[1] = -beta * s * i * pSI + pRS * lambda * r  # ds/dt
        du[2] = beta * s * i * pSI - gamma * (1 - pII) * i  # di/dt
        du[3] = pIH * gamma * i - alpha * (1 - pHH) * h  # dh/dt
        du[4] = pIR * gamma * i + pHR * alpha * h - pRS * lambda * r  # dr/dt
        du[5] = pID * gamma * i + pHD * alpha * h  # dd/dt
    end
    
    # Parameters for ODE
    p = [params["beta"], params["gamma"], params["alpha"], params["lambda"],
         params["pSI"], params["pII"], params["pIH"], params["pIR"], params["pID"],
         params["pHH"], params["pHR"], params["pHD"], params["pRR"], params["pRS"]]
    
    # Initial conditions
    u0 = [params["s0"], params["i0"], params["h0"], params["r0"], params["d0"]]
    
    # Time span
    tspan = (0.0, params["tmax"])
    
    # Solve ODE
    prob = ODEProblem(ode_system!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)
    
    return (T=sol.t, S_prop=sol[1,:], I_prop=sol[2,:], H_prop=sol[3,:], R_prop=sol[4,:], D_prop=sol[5,:])
end

"""
Interpolate variance data onto common time grid
"""
function interpolate_variance(T, V, t_grid)
    # Simple linear interpolation
    result = zeros(length(t_grid))
    for (i, t) in enumerate(t_grid)
        if t <= T[1]
            result[i] = V[1]
        elseif t >= T[end]
            result[i] = V[end]
        else
            # Find the two closest time points
            idx = findfirst(x -> x >= t, T)
            if idx === nothing
                result[i] = V[end]
            elseif idx == 1
                result[i] = V[1]
            else
                t1, t2 = T[idx-1], T[idx]
                v1, v2 = V[idx-1], V[idx]
                # Linear interpolation
                result[i] = v1 + (v2 - v1) * (t - t1) / (t2 - t1)
            end
        end
    end
    return result
end

"""
Create error bar plots showing final expected variance (average of 5 batches) vs theoretical
"""
function plot_variance_error_bars(N_values, t_grid, S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var,
                                v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory)
    # Create error bar plots showing final expected variance (average of 5 batches) vs theoretical
    # Following the new conjecture: E[V_N^(l)(t)] proportional to (1/N) x population_proportions
    dt = 2  # Half width for time intervals around each midpoint
    midpoints = 5:dt*2:1000  # Time points for analysis (covering most of 1000-day simulation)
    compartments = ["Susceptible", "Infected", "Hospitalized", "Recovered", "Dead"]
    all_sim_vars = [S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var]
    all_theory_vars = [v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory]
    colors = [:blue, :orange, :green, :red, :purple]
    
    # Create separate figures for each compartment
    for comp_idx in 1:length(compartments)
        p = plot(layout=(1, length(N_values)), size=(1600, 800))
        
        for n_idx in 1:length(N_values)
            N = N_values[n_idx]
            
            # Get batch data for this compartment and N
            batch_data = all_sim_vars[comp_idx][n_idx, :, :]  # [time_points x 5_batches]
            theory_data = all_theory_vars[comp_idx][n_idx]
            
            # Plot theoretical variance curve (smooth line)
            # Only add label to first subplot, use empty string for others
            if n_idx == 1
                plot!(p, t_grid, theory_data, subplot=n_idx, 
                      linestyle=:dash, color=:black, linewidth=2, 
                      label="Theoretical Curve")
            else
                plot!(p, t_grid, theory_data, subplot=n_idx, 
                      linestyle=:dash, color=:black, linewidth=2, 
                      label="")
            end
            
            # Create arrays to store all error bar data for this N
            error_bar_times = Float64[]
            error_bar_mins = Float64[]
            error_bar_maxs = Float64[]
            
            # Calculate and collect error bar data
            for t_idx in 1:length(midpoints)
                t_mid = midpoints[t_idx]
                t_min = t_mid - dt
                t_max = t_mid + dt
                
                # Find time points in window
                in_window_indices = findall(t -> t >= t_min && t < t_max, t_grid)
                if isempty(in_window_indices)
                    continue
                end
                
                # Get batch values in this time window
                batch_vals_in_window = batch_data[in_window_indices, :]
                
                # Calculate min and max across batches in this time window
                batch_mins = minimum(batch_vals_in_window, dims=1)[1, :]  # Min of each batch
                batch_maxs = maximum(batch_vals_in_window, dims=1)[1, :]  # Max of each batch
                
                # Calculate final expected variance bounds (average across 5 batches)
                y_min_avg = mean(batch_mins)
                y_max_avg = mean(batch_maxs)
                
                if !isnan(y_min_avg) && !isnan(y_max_avg)
                    push!(error_bar_times, t_mid)
                    push!(error_bar_mins, y_min_avg)
                    push!(error_bar_maxs, y_max_avg)
                end
            end
            
            # Plot all error bars at once with proper legend handling
            if !isempty(error_bar_times)
                # Only add label to first subplot, use empty string for others
                if n_idx == 1
                    plot!(p, error_bar_times, error_bar_mins, subplot=n_idx, 
                          color=colors[n_idx], linewidth=1.5, 
                          label="Expected Variance Range")
                else
                    plot!(p, error_bar_times, error_bar_mins, subplot=n_idx, 
                          color=colors[n_idx], linewidth=1.5, 
                          label="", legend=false)
                end
                
                # Add the rest without labels - explicitly set legend=false
                plot!(p, error_bar_times, error_bar_maxs, subplot=n_idx, 
                      color=colors[n_idx], linewidth=1.5, 
                      label="", legend=false)
                
                # Add vertical lines connecting min to max without labels - explicitly set legend=false
                for i in 1:length(error_bar_times)
                    plot!(p, [error_bar_times[i], error_bar_times[i]], 
                          [error_bar_mins[i], error_bar_maxs[i]], 
                          subplot=n_idx, color=colors[n_idx], linewidth=1.5, 
                          label="", legend=false)
                end
            end
            
            # Set subplot properties
            plot!(p, subplot=n_idx, xlabel="Time", ylabel="Variance", 
                  title="$(compartments[comp_idx]) - N = $N (5 batches x 15 runs)",
                  grid=true, legend=:best)
        end
        
        # Add overall title
        plot!(p, title="$(compartments[comp_idx]) - Variance Analysis", titlefontsize=16)
        
        # Display plot
        display(p)
        
        # Save plot
        savefig(p, "SIHRS_$(compartments[comp_idx])_final_expected_variance.png")
    end
    
    println("Created separate figures for each compartment:")
    for comp_idx in 1:length(compartments)
        println("  - $(compartments[comp_idx]): SIHRS_$(compartments[comp_idx])_final_expected_variance.png")
    end
end

# Run the analysis
if abspath(PROGRAM_FILE) == @__FILE__
    SIHRS_variance_analysis()
end
