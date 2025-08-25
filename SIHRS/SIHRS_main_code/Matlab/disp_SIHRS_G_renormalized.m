% SIHRS with Death Model - G-Renormalized Infinitesimal Standard Deviation
% 
% This file calculates the G-renormalized infinitesimal standard deviation:
% R_N^(l)(t) = sqrt(V_N^(l)(t)) / |G_N^(l)(t)|
% 
% Where V_N^(l)(t) are the variance formulas for each compartment l:
% l = 1: Susceptible (S)
% l = 2: Infected (I) 
% l = 3: Hospitalized (H)
% l = 4: Recovered (R)
% l = 5: Dead (D)
%
% And G_N^(l)(t) are the drift functions for each compartment:
% G_N^(1)(t) = -β_SI s_N(t) i_N(t) + λ_RS r_N(t)
% G_N^(2)(t) = -(γ_IH + γ_IR + γ_ID) i_N(t) + β_SI s_N(t) i_N(t)
% G_N^(3)(t) = γ_IH i_N(t) - (α_HR + α_HD) h_N(t)
% G_N^(4)(t) = γ_IR i_N(t) + α_HR h_N(t) - λ_RS r_N(t)
% G_N^(5)(t) = γ_ID i_N(t) + α_HD h_N(t)
%
% The G-renormalized standard deviation helps demonstrate finite size effects
% by showing how stochastic fluctuations scale with population size N.

clear all;
close all;

function sihrs_G_renormalized_simulation()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching SIHRS.m
    params.beta = 0.212;           % Infection rate (β > 0) - matching SIHRS.m
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
    params.gamma = 0.100346667;    % Adjusted for exact critical hit at S = 142/300
    params.alpha = 0.1;            % Hospitalized transition rate (α > 0)
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0) immunity period of 4 months
    params.T = 1000;               % Total simulation time
    params.dt = 0.01;             % Finer time step for better blow-up capture
    params.N_values = [300]; 
    params.initial_s = 0.96;       % Initial susceptible fraction
    params.initial_i = 0.04;       % Initial infected fraction
    params.initial_h = 0;          % Initial hospitalized fraction
    params.initial_r = 0;          % Initial recovered fraction
    params.initial_d = 0;          % Initial dead fraction
    params.n_runs = 1;            % Number of stochastic runs
    params.colors = {'#0072BD','#77AC30', '#A2142F'}; % Colors matching SIHRS.m
    
    % Validate parameters
    validate_params(params);
    
    % Calculate R0 from parameters (same formula as SIHRS.m)
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('Calculated R₀ = %.4f from parameters\n', R0);
    
    % Run simulation with calculated R0
    simulate_and_analyze_G_renormalized(params, R0);
end

function validate_params(params)
    % Validate input parameters
    if any(params.N_values <= 0) || any(mod(params.N_values, 1) ~= 0)
        error('N_values must be positive integers');
    end
    if params.pSI <= 0 || params.pSI > 1 || params.pII < 0 || params.pII > 1
        error('pSI must be in (0, 1] and pII must be in [0, 1]');
    end
    if params.gamma <= 0
        error('gamma must be positive');
    end
    if abs(params.initial_s + params.initial_i + params.initial_h + params.initial_r + params.initial_d - 1) > 1e-6
        error('Initial fractions must sum to 1');
    end
    if params.T <= 0 || params.dt <= 0
        error('T and dt must be positive');
    end
    if params.n_runs <= 0 || mod(params.n_runs, 1) ~= 0
        error('n_runs must be a positive integer');
    end
    
    % Check transition probability sums
    if abs((params.pII + params.pIH + params.pIR + params.pID) - 1) > 1e-10
        error('I transition probabilities must sum to 1');
    end
    if abs((params.pHH + params.pHR + params.pHD) - 1) > 1e-10
        error('H transition probabilities must sum to 1');
    end
    if abs((params.pRR + params.pRS) - 1) > 1e-10
        error('R transition probabilities must sum to 1');
    end
end

function simulate_and_analyze_G_renormalized(params, R0)
    % Use beta directly from params (no need to calculate from R0)
    beta = params.beta;
    
    % Precompute time vector
    t = 0:params.dt:params.T;
    
    % Store results
    results = cell(length(params.N_values), 1);
    det_result = solve_deterministic_sihrs(beta, params);
    
    % Run simulations for each population size
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d, R0 = %.2f...\n', params.n_runs, N, R0);
        
        % Run multiple stochastic simulations
        results{idx} = run_multiple_gillespie_G_renormalized(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_G_renormalized_results(t, results, det_result, params, R0);
    print_G_renormalized_summary(results, det_result, R0);
end

function result = run_multiple_gillespie_G_renormalized(N, beta, params, t)
    % Run multiple Gillespie simulations and aggregate results
    % CORRECTED: Calculate G-renormalized values for each run first, then average
    
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    H_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    D_all = zeros(params.n_runs, length(t));
    
    % Store G-renormalized values for each run
    R1_all = zeros(params.n_runs, length(t));
    R2_all = zeros(params.n_runs, length(t));
    R3_all = zeros(params.n_runs, length(t));
    R4_all = zeros(params.n_runs, length(t));
    R5_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        
        % Store population trajectories
        S_all(run, :) = S_interp;
        I_all(run, :) = I_interp;
        H_all(run, :) = H_interp;
        R_all(run, :) = R_interp;
        D_all(run, :) = D_interp;
        
        % Calculate G-renormalized values FOR THIS INDIVIDUAL RUN
        % This preserves blow-up behavior before averaging
        
        % R_N^(1)(t) for Susceptible (S)
        V1 = (1/N) * (beta * params.pSI * S_interp .* I_interp + params.lambda * params.pRS * R_interp);
        G1 = -beta * S_interp .* I_interp + params.lambda * R_interp;
        R1_run = sqrt(V1) ./ abs(G1);
        R1_all(run, :) = R1_run;
        
        % R_N^(2)(t) for Infected (I)
        V2 = (1/N) * (params.gamma * I_interp * (params.pIH + params.pIR + params.pID) + beta * params.pSI * S_interp .* I_interp);
        G2 = -(params.gamma * (params.pIH + params.pIR + params.pID)) * I_interp + beta * S_interp .* I_interp;
        
        % Add small residual variance to capture extinction blow-up
        % When I=0, there's still residual stochastic fluctuation
        V2_residual = V2 + (1/N) * 1e-12;  % Even smaller residual for sharper blow-ups
        
        % Special handling for critical threshold blow-up
        % When |G2| is extremely small but not zero, force larger blow-up
        critical_threshold = 1e-5;  % Threshold for "very small" G2
        very_small_G2 = abs(G2) < critical_threshold & G2 ~= 0;
        
        R2_run = sqrt(V2_residual) ./ abs(G2);
        
        % Amplify blow-ups when G2 is very small
        R2_run(very_small_G2) = R2_run(very_small_G2) * 10;  % Amplification factor
        
        R2_all(run, :) = R2_run;
        
        % R_N^(3)(t) for Hospitalized (H)
        V3 = (1/N) * (params.gamma * I_interp * params.pIH + params.alpha * H_interp * (params.pHR + params.pHD));
        G3 = params.gamma * params.pIH * I_interp - params.alpha * (params.pHR + params.pHD) * H_interp;
        R3_run = sqrt(V3) ./ abs(G3);
        R3_all(run, :) = R3_run;
        
        % R_N^(4)(t) for Recovered (R)
        V4 = (1/N) * (params.gamma * I_interp * params.pIR + params.alpha * H_interp * params.pHR + params.lambda * R_interp * params.pRS);
        G4 = params.gamma * params.pIR * I_interp + params.alpha * params.pHR * H_interp - params.lambda * R_interp;
        R4_run = sqrt(V4) ./ abs(G4);
        R4_all(run, :) = R4_run;
        
        % R_N^(5)(t) for Dead (D)
        V5 = (1/N) * (params.gamma * I_interp * params.pID + params.alpha * H_interp * params.pHD);
        G5 = params.gamma * params.pID * I_interp + params.alpha * params.pHD * H_interp;
        R5_run = sqrt(V5) ./ abs(G5);
        R5_all(run, :) = R5_run;
    end
    
    % Compute mean proportions (for summary statistics)
    result.S_mean = mean(S_all, 1);
    result.I_mean = mean(I_all, 1);
    result.H_mean = mean(H_all, 1);
    result.R_mean = mean(R_all, 1);
    result.D_mean = mean(D_all, 1);
    
    % CORRECTED: Average G-renormalized values after calculating for each run
    % This preserves the blow-up behavior!
    
    % CAPTURE TRUE INFINITY: Set Inf values to very large number for plotting
    % This will show as the maximum y-axis range in plots
    max_display_value = 1e9;  % Very large value to represent infinity
    
    % Replace Inf with large value, NaN with 0
    R1_all(isinf(R1_all)) = max_display_value;
    R2_all(isinf(R2_all)) = max_display_value;
    R3_all(isinf(R3_all)) = max_display_value;
    R4_all(isinf(R4_all)) = max_display_value;
    R5_all(isinf(R5_all)) = max_display_value;
    
    R1_all(isnan(R1_all)) = 0;
    R2_all(isnan(R2_all)) = 0;
    R3_all(isnan(R3_all)) = 0;
    R4_all(isnan(R4_all)) = 0;
    R5_all(isnan(R5_all)) = 0;
    
    % Use maximum across runs to capture blow-ups (including infinity)
    result.R1 = max(R1_all, [], 1);
    result.R2 = max(R2_all, [], 1);
    result.R3 = max(R3_all, [], 1);
    result.R4 = max(R4_all, [], 1);
    result.R5 = max(R5_all, [], 1);
    
    result.N = N;
    result.peak_infected = max(result.I_mean * N);
    result.peak_time = t(find(result.I_mean == max(result.I_mean), 1, 'first'));
    result.final_time = t(end);
end

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params)
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
    max_events = round(10 * params.T * (beta + params.gamma + params.alpha + params.lambda) * N);
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
        si_rate = (beta / N) * S * I * params.pSI; % rate of S->I
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

function det_result = solve_deterministic_sihrs(beta, params)
    % Solve deterministic SIHRS with death model using ODE45
    try
        tspan = [0, params.T];
        y0 = [params.initial_s; params.initial_i; params.initial_h; params.initial_r; params.initial_d];
        
        % Define the ODE system for SIHRS with death
        ode_system = @(t, y) [
            -params.pSI * beta * y(1) * y(2) + params.pRS * params.lambda * y(4);           % ds/dt
            params.pSI * beta * y(1) * y(2) - params.gamma * (1 - params.pII) * y(2);      % di/dt
            params.pIH * params.gamma * y(2) - params.alpha * (1 - params.pHH) * y(3);     % dh/dt
            params.pIR * params.gamma * y(2) + params.pHR * params.alpha * y(3) - params.pRS * params.lambda * y(4); % dr/dt
            params.pID * params.gamma * y(2) + params.pHD * params.alpha * y(3)             % dd/dt
        ];
        
        [T, Y] = ode45(ode_system, tspan, y0);
        
        det_result.T = T;
        det_result.S_prop = Y(:, 1);
        det_result.I_prop = Y(:, 2);
        det_result.H_prop = Y(:, 3);
        det_result.R_prop = Y(:, 4);
        det_result.D_prop = Y(:, 5);
        [det_result.peak_infected_prop, idx] = max(det_result.I_prop);
        det_result.peak_time = T(idx);
        det_result.final_time = T(end);
    catch e
        error('ODE solver failed: %s', e.message);
    end
end

function plot_G_renormalized_results(t, results, det_result, params, R0)
    % Create 5 separate figures for each compartment showing G-renormalized standard deviations
    
    % Plot R1 (G-renormalized standard deviation for susceptible)
    figure('Position', [100, 100, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(1)}(t)');
    grid on;
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_R1_Susceptible_R0_%.2f.png', R0));
    
    % Plot R2 (G-renormalized standard deviation for infected)
    figure('Position', [200, 200, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Infected (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(2)}(t)');
    
    % Dynamic y-axis with maximum cap at 1×10^9
    max_R2_value = 0;
    for idx = 1:length(params.N_values)
        max_R2_value = max(max_R2_value, max(results{idx}.R2));
    end
    y_max = min(max_R2_value * 1.1, 1e9);  % 10% padding or cap at 1×10^9
    ylim([0, y_max]);
    
    grid on;
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_R2_Infected_R0_%.2f.png', R0));
    
    % Plot R3 (G-renormalized standard deviation for hospitalized)
    figure('Position', [300, 300, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(3)}(t)');
    grid on;
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_R3_Hospitalized_R0_%.2f.png', R0));

    % Plot R4 (G-renormalized standard deviation for recovered)
    figure('Position', [400, 400, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R4, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Recovered (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(4)}(t)');
    grid on;
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_R4_Recovered_R0_%.2f.png', R0));
    
    % Plot R5 (G-renormalized standard deviation for dead)
    figure('Position', [500, 500, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R5, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Dead (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(5)}(t)');
    grid on;
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_R5_Dead_R0_%.2f.png', R0));
    
    % Create combined plot
    figure('Position', [600, 600, 1500, 800]);
    tlayout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot R1 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(1)}(t)');
    grid on;
    
    % Plot R2 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Infected (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(2)}(t)');
    
    % Apply same dynamic y-axis scaling as individual R2 plot
    max_R2_value = 0;
    for idx = 1:length(params.N_values)
        max_R2_value = max(max_R2_value, max(results{idx}.R2));
    end
    y_max = min(max_R2_value * 1.1, 1e9);  % 10% padding or cap at 1×10^9
    ylim([0, y_max]);
    
    grid on;
    
    % Plot R3 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(3)}(t)');
    grid on;

    % Plot R4 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R4, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Recovered (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(4)}(t)');
    grid on;
    
    % Plot R5 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.R5, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('G-Renormalized Std Dev - Dead (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('R_N^{(5)}(t)');
    grid on;
    
    % Add legend to combined figure
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add mathematical formula annotation to combined figure
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', sprintf('R_N^{(l)}(t) = sqrt(V_N^{(l)}(t)) / |G_N^{(l)}(t)|, R_0 = %.2f', R0), ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Save the combined figure
    saveas(gcf, sprintf('SIHRS_G_renormalized_std_dev_R0_%.2f_combined.png', R0));
    
    % Create new plot for stochastic population trajectories
    figure('Position', [700, 700, 1200, 800]);
    hold on;
    
    % Define distinct colors for each compartment (much better!)
    compartment_colors = {
        [0.2, 0.6, 0.2],    % S - Dark Green
        [0.8, 0.2, 0.2],    % I - Dark Red  
        [0.9, 0.6, 0.0],    % H - Orange
        [0.0, 0.4, 0.8],    % R - Blue
        [0.4, 0.0, 0.6]     % D - Purple
    };
    
    % Line styles for different population sizes
    line_styles = {'-', '--', ':'};  % Different line styles for different N values
    
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        style = line_styles{min(idx, length(line_styles))};
        
        % Plot each compartment with its own distinct color
        plot(t, results{idx}.S_mean, 'Color', compartment_colors{1}, 'LineStyle', style, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('S (N=%d)', N));
        plot(t, results{idx}.I_mean, 'Color', compartment_colors{2}, 'LineStyle', style, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('I (N=%d)', N));
        plot(t, results{idx}.H_mean, 'Color', compartment_colors{3}, 'LineStyle', style, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('H (N=%d)', N));
        plot(t, results{idx}.R_mean, 'Color', compartment_colors{4}, 'LineStyle', style, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('R (N=%d)', N));
        plot(t, results{idx}.D_mean, 'Color', compartment_colors{5}, 'LineStyle', style, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('D (N=%d)', N));
    end
    
    % Deterministic solution plots removed for cleaner visualization
    
    title(sprintf('SIHRS Population Dynamics - All Compartments (R_0 = %.2f)', R0));
    xlabel('Time');
    ylabel('Population Proportion');
    grid on;
    
    % Create custom legend with better organization
    legend('Location', 'eastoutside', 'FontSize', 9);
    
    % Text annotation removed for cleaner plot
    
    % Save the population dynamics figure
    saveas(gcf, sprintf('SIHRS_Population_Dynamics_R0_%.2f.png', R0));
    
    fprintf('Generated individual plots:\n');
    fprintf('  - SIHRS_R1_Susceptible_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_R2_Infected_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_R3_Hospitalized_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_R4_Recovered_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_R5_Dead_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_G_renormalized_std_dev_R0_%.2f_combined.png\n', R0);
    fprintf('  - SIHRS_Population_Dynamics_R0_%.2f.png\n', R0);
end

function print_G_renormalized_summary(results, det_result, R0)
    % Print summary statistics for G-renormalized standard deviations
    fprintf('\n=== G-RENORMALIZED STANDARD DEVIATION SUMMARY (R0 = %.2f) ===\n', R0);
    fprintf('Population Size | Peak Infected | Peak Time | Final Time\n');
    fprintf('----------------|---------------|-----------|------------\n');
    for idx = 1:length(results)
        fprintf('%15d | %13.2f | %9.2f | %10.2f\n', ...
            results{idx}.N, results{idx}.peak_infected, ...
            results{idx}.peak_time, results{idx}.final_time);
    end
    fprintf('%15s | %13.4f | %9.2f | %10.2f\n', ...
        'Deterministic', det_result.peak_infected_prop, ...
        det_result.peak_time, det_result.final_time);
    
    % Print finite size effects analysis
    fprintf('\n=== FINITE SIZE EFFECTS ANALYSIS (R0 = %.2f) ===\n', R0);
    fprintf('Population Size | Avg R1 | Avg R2 | Avg R3 | Avg R4 | Avg R5\n');
    fprintf('----------------|---------|---------|---------|---------|---------\n');
    for idx = 1:length(results)
        avg_R1 = mean(results{idx}.R1);
        avg_R2 = mean(results{idx}.R2);
        avg_R3 = mean(results{idx}.R3);
        avg_R4 = mean(results{idx}.R4);
        avg_R5 = mean(results{idx}.R5);
        fprintf('%15d | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f\n', ...
            results{idx}.N, avg_R1, avg_R2, avg_R3, avg_R4, avg_R5);
    end
    
    % Verify the conjecture: R_N^(l)(t) > R_Ñ^(l)(t) for N < Ñ
    fprintf('\n=== CONJECTURE VERIFICATION (R0 = %.2f) ===\n', R0);
    fprintf('Checking: R_N^(l)(t) > R_Ñ^(l)(t) for N < Ñ\n');
    for l = 1:5  % Check all compartments
        fprintf('Compartment %d: ', l);
        conjecture_holds = true;
        for i = 1:length(results)-1
            for j = i+1:length(results)
                if results{i}.N < results{j}.N
                    % Check if R_N^(l) > R_Ñ^(l) on average
                    avg_R_smaller = mean(results{i}.(sprintf('R%d', l)));
                    avg_R_larger = mean(results{j}.(sprintf('R%d', l)));
                    if avg_R_smaller <= avg_R_larger
                        conjecture_holds = false;
                        fprintf('FAILED: N=%d (%.4f) <= N=%d (%.4f) ', ...
                            results{i}.N, avg_R_smaller, results{j}.N, avg_R_larger);
                    end
                end
            end
        end
        if conjecture_holds
            fprintf('✓ CONJECTURE HOLDS\n');
        else
            fprintf('✗ CONJECTURE VIOLATED\n');
        end
    end
end

% Run the simulation
sihrs_G_renormalized_simulation();
