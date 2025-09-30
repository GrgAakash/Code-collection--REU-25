% Generate CSV file with SIHRS population data and D_N values
% This script runs the SIHRS simulation and outputs the data to CSV

clear all;
close all;

function generate_sihrs_csv()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching disp_SIHRS_renormalized.m
    params.beta = 0.12;           % Infection rate (β > 0)
    params.pSI = 1.0;              % Infection probability (S to I)
    params.pII = 0.0;              % probability of I to I (stay infected)
    params.pIH = 0.04;             % probability of I to H
    params.pIR = 0.959;            % probability of I to R
    params.pID = 0.001;            % probability of I to D
    params.pHH = 0.01;             % probability of H to H (stay hospitalized)
    params.pHR = 0.9882;           % probability of H to R
    params.pHD = 0.0018;           % probability of H to D
    params.pRR = 0.02;             % probability of R to R (stay recovered)
    params.pRS = 0.98;             % probability of R to S
    params.gamma = 0.1;    % Adjusted for exact critical hit at S = 142/300
    params.alpha = 0.1;            % Hospitalized transition rate (α > 0)
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0) immunity period of 4 months
    params.T = 300;               % Total simulation time
    params.dt = 0.01;             % Finer time step for better blow-up capture
    params.N_values = [1600]; 
    params.initial_s = 0.96;       % Initial susceptible fraction
    params.initial_i = 0.04;       % Initial infected fraction
    params.initial_h = 0;          % Initial hospitalized fraction
    params.initial_r = 0;          % Initial recovered fraction
    params.initial_d = 0;          % Initial dead fraction
    params.n_runs = 1;            % Number of stochastic runs
    
    % Calculate R0 from parameters
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('Calculated R₀ = %.4f from parameters\n', R0);
    
    % Precompute time vector
    t = 0:params.dt:params.T;
    
    % Run simulation for each population size
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d, R0 = %.2f...\n', params.n_runs, N, R0);
        
        % Run multiple stochastic simulations
        result = run_multiple_gillespie_renormalized(N, params, t);
        
        % Calculate D_N values using pure mathematical formulas
        result = calculate_dn_values(result, N, params);
        
        % Save population data and D_N values to CSV file
        population_data = table(t', result.S_mean', result.I_mean', result.H_mean', result.R_mean', result.D_mean', ...
                               result.D1', result.D2', result.D3', result.D4', result.D5', ...
            'VariableNames', {'Time', 'S_mean', 'I_mean', 'H_mean', 'R_mean', 'D_mean', ...
                              'D1_Susceptible', 'D2_Infected', 'D3_Hospitalized', 'D4_Recovered', 'D5_Dead'});
        csv_filename = sprintf('SIHRS_population_and_DN_data_N%d_R0_%.2f.csv', N, R0);
        writetable(population_data, csv_filename);
        fprintf('Population data and D_N values saved to: %s\n', csv_filename);
        
        % Print summary statistics
        fprintf('\nSummary for N = %d:\n', N);
        fprintf('  Susceptible: Mean D_N = %.6f, Max D_N = %.6f\n', mean(result.D1), max(result.D1));
        fprintf('  Infected: Mean D_N = %.6f, Max D_N = %.6f\n', mean(result.D2), max(result.D2));
        fprintf('  Hospitalized: Mean D_N = %.6f, Max D_N = %.6f\n', mean(result.D3), max(result.D3));
        fprintf('  Recovered: Mean D_N = %.6f, Max D_N = %.6f\n', mean(result.D4), max(result.D4));
        fprintf('  Dead: Mean D_N = %.6f, Max D_N = %.6f\n', mean(result.D5), max(result.D5));
    end
end

function result = run_multiple_gillespie_renormalized(N, params, t)
    % Run multiple Gillespie simulations and aggregate results
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    H_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    D_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        S_all(run, :) = S_interp;
        I_all(run, :) = I_interp;
        H_all(run, :) = H_interp;
        R_all(run, :) = R_interp;
        D_all(run, :) = D_interp;
    end
    
    % Compute mean proportions
    result.S_mean = mean(S_all, 1);
    result.I_mean = mean(I_all, 1);
    result.H_mean = mean(H_all, 1);
    result.R_mean = mean(R_all, 1);
    result.D_mean = mean(D_all, 1);
end

function result = calculate_dn_values(result, N, params)
    % Calculate D_N values using pure mathematical formulas - no safety checks!
    
    % D_N^(1)(t) for Susceptible (S)
    % (D_N^(1))^2 = (1/N) * (β_SI * i_N(t)/s_N(t) + λ_RS * r_N(t)/(s_N(t))^2)
    result.D1_squared = (1/N) * (params.beta * params.pSI * result.I_mean ./ result.S_mean + ...
                                params.lambda * params.pRS * result.R_mean ./ (result.S_mean.^2));
    result.D1 = sqrt(result.D1_squared);
    
    % D_N^(2)(t) for Infected (I)
    % (D_N^(2))^2 = (1/N) * ((γ_IH + γ_IR + γ_ID)/i_N(t) + β_SI * s_N(t)/i_N(t))
    result.D2_squared = (1/N) * ((params.gamma * (params.pIH + params.pIR + params.pID)) ./ result.I_mean + ...
                                params.beta * params.pSI * result.S_mean ./ result.I_mean);
    result.D2 = sqrt(result.D2_squared);
    
    % D_N^(3)(t) for Hospitalized (H)
    % (D_N^(3))^2 = (1/N) * ((α_HR + α_HD)/h_N(t) + γ_IH * i_N(t)/(h_N(t))^2)
    result.D3_squared = (1/N) * ((params.alpha * (params.pHR + params.pHD)) ./ result.H_mean + ...
                                params.gamma * params.pIH * result.I_mean ./ (result.H_mean.^2));
    result.D3 = sqrt(result.D3_squared);
    
    % D_N^(4)(t) for Recovered (R)
    % (D_N^(4))^2 = (1/N) * (λ_RS/r_N(t) + γ_IR * i_N(t)/(r_N(t))^2 + α_HR * h_N(t)/(r_N(t))^2)
    result.D4_squared = (1/N) * (params.lambda * params.pRS ./ result.R_mean + ...
                                params.gamma * params.pIR * result.I_mean ./ (result.R_mean.^2) + ...
                                params.alpha * params.pHR * result.H_mean ./ (result.R_mean.^2));
    result.D4 = sqrt(result.D4_squared);
    
    % D_N^(5)(t) for Dead (D)
    % V5 = (1/N) * (pID * γ * i_N(t) + pHD * α * h_N(t))
    % D_N^(5) = √(V5)/d_N(t)
    V5 = (1/N) * (params.pID * params.gamma * result.I_mean + params.pHD * params.alpha * result.H_mean);
    result.D5 = sqrt(V5) ./ result.D_mean;
end

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, params)
    % Initialize populations
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
    % Verify population conservation
    if abs(S + I + H + R + D - N) > 1
        error('Initial populations do not sum to N');
    end
    
    % Preallocate arrays
    max_events = round(10 * params.T * (params.beta + params.gamma + params.alpha + params.lambda) * N);
    S_hist = zeros(1, max_events);
    I_hist = zeros(1, max_events);
    H_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events);
    D_hist = zeros(1, max_events);
    time_pts = zeros(1, max_events);
    
    S_hist(1) = S;
    I_hist(1) = I;
    H_hist(1) = H;
    R_hist(1) = R;
    D_hist(1) = D;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    % Gillespie algorithm for SIHRS with death
    while current_time < params.T
        
        % Calculate event rates
        si_rate = (params.beta / N) * S * I * params.pSI; % rate of S->I
        ir_rate = params.gamma * I * params.pIR; % rate of I->R
        ih_rate = params.gamma * I * params.pIH; % rate of I->H
        id_rate = params.gamma * I * params.pID; % rate of I->D
        hr_rate = params.alpha * H * params.pHR; % rate of H->R
        hd_rate = params.alpha * H * params.pHD; % rate of H->D
        rs_rate = params.lambda * R * params.pRS; % rate of R->S
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        
        if total_rate == 0
            break;
        end
        
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        
        if current_time > params.T
            break;
        end
        
        % Determine which event occurs
        chance = rand * total_rate;
        if chance < si_rate
            % S->I transition
            if S > 0
                S = S - 1;
                I = I + 1;
            end
        elseif chance < (si_rate + ir_rate)
            % I->R transition
            if I > 0
                I = I - 1;
                R = R + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate)
            % I->H transition
            if I > 0
                I = I - 1;
                H = H + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            % I->D transition
            if I > 0
                I = I - 1;
                D = D + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            % H->R transition
            if H > 0
                H = H - 1;
                R = R + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            % H->D transition
            if H > 0
                H = H - 1;
                D = D + 1;
            end
        else
            % R->S transition
            if R > 0
                R = R - 1;
                S = S + 1;
            end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S;
        I_hist(event_count) = I;
        H_hist(event_count) = H;
        R_hist(event_count) = R;
        D_hist(event_count) = D;
        time_pts(event_count) = current_time;
    end
    
    % Trim arrays
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    D_hist = D_hist(1:event_count);
    time_pts = time_pts(1:event_count);
    
    % Verify population conservation
    assert(all(S_hist + I_hist + H_hist + R_hist + D_hist == N), 'Population not conserved');
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    % Interpolate to fixed time grid using 'previous'
    try
        S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
        I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
        H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
        R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
        D_interp = interp1(time_pts, D_hist, t, 'previous') / N;
        
        % Handle values beyond last event
        S_interp(t > max(time_pts)) = S_hist(end) / N;
        I_interp(t > max(time_pts)) = I_hist(end) / N;
        H_interp(t > max(time_pts)) = H_hist(end) / N;
        R_interp(t > max(time_pts)) = R_hist(end) / N;
        D_interp(t > max(time_pts)) = D_hist(end) / N;
    catch e
        error('Interpolation failed: %s', e.message);
    end
end

% Run the simulation and generate CSV
generate_sihrs_csv();
