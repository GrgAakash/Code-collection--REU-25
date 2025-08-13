% SIHRS with Death Model - Renormalized Infinitesimal Standard Deviation
% 
% This file calculates the renormalized infinitesimal standard deviation:
% D_N^(l)(t) = sqrt(V_N^(l)(t)) / |m_N^(l)(t)|
% 
% Where V_N^(l)(t) are the variance formulas for each compartment l:
% l = 1: Susceptible (S)
% l = 2: Infected (I) 
% l = 3: Hospitalized (H)
% l = 4: Recovered (R)
% l = 5: Dead (D)
%
% The renormalized standard deviation helps demonstrate finite size effects
% by showing how stochastic fluctuations scale with population size N.

clear all;
close all;

function sihrs_renormalized_simulation()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching SIHRS.m
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
    params.gamma = 0.1;            % Infection transition rate (γ > 0)
    params.alpha = 0.1;            % Hospitalized transition rate (α > 0)
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0) immunity period of 4 months
    params.T = 1000;               % Total simulation time
    params.dt = 0.01;              % Time step for interpolation
    params.N_values = [316, 3162, 10000]; % Population sizes - matching SIHRS.m
    params.initial_s = 0.96;       % Initial susceptible fraction
    params.initial_i = 0.04;       % Initial infected fraction
    params.initial_h = 0;          % Initial hospitalized fraction
    params.initial_r = 0;          % Initial recovered fraction
    params.initial_d = 0;          % Initial dead fraction
    params.n_runs = 40;            % Number of stochastic runs
    params.colors = {'#0072BD', '#77AC30', '#A2142F'}; % Colors matching SIHRS.m
    params.R0_values = [2.12]; % R0 values
    
    % Validate parameters
    validate_params(params);
    
    % Run simulations for each R0
    for r_idx = 1:length(params.R0_values)
        R0 = params.R0_values(r_idx);
        simulate_and_analyze_renormalized(params, R0);
    end
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

function simulate_and_analyze_renormalized(params, R0)
    % Calculate beta from R0
    beta = R0 * params.gamma / params.pSI;
    
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
        results{idx} = run_multiple_gillespie_renormalized(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_renormalized_results(t, results, det_result, params, R0);
    print_renormalized_summary(results, det_result, R0);
end

function result = run_multiple_gillespie_renormalized(N, beta, params, t)
    % Run multiple Gillespie simulations and aggregate results
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    H_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    D_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
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
    
    % Compute renormalized infinitesimal standard deviations D_N^(l)(t)
    % Using the mathematical formulas provided
    % Note: We ensure denominators are non-zero to avoid division by zero
    
    % D_N^(1)(t) for Susceptible (S)
    % (D_N^(1))^2 = (1/N) * (β_SI * i_N(t)/s_N(t) + λ_RS * r_N(t)/(s_N(t))^2)
    % Check if s_N(t) is non-zero
    valid_S = result.S_mean > 1e-10;  % Threshold for non-zero values
    result.D1_squared = zeros(size(result.S_mean));
    result.D1_squared(valid_S) = (1/N) * (beta * params.pSI * result.I_mean(valid_S) ./ result.S_mean(valid_S) + ...
                                    params.lambda * params.pRS * result.R_mean(valid_S) ./ (result.S_mean(valid_S).^2));
    result.D1 = sqrt(result.D1_squared);
    
    % D_N^(2)(t) for Infected (I)
    % (D_N^(2))^2 = (1/N) * ((γ_IH + γ_IR + γ_ID)/i_N(t) + β_SI * s_N(t)/i_N(t))
    % Check if i_N(t) is non-zero
    valid_I = result.I_mean > 1e-10;  % Threshold for non-zero values
    result.D2_squared = zeros(size(result.I_mean));
    result.D2_squared(valid_I) = (1/N) * ((params.gamma * (params.pIH + params.pIR + params.pID)) ./ result.I_mean(valid_I) + ...
                                    beta * params.pSI * result.S_mean(valid_I) ./ result.I_mean(valid_I));
    result.D2 = sqrt(result.D2_squared);
    
    % D_N^(3)(t) for Hospitalized (H)
    % (D_N^(3))^2 = (1/N) * ((α_HR + α_HD)/h_N(t) + γ_IH * i_N(t)/(h_N(t))^2)
    % Check if h_N(t) is non-zero
    valid_H = result.H_mean > 1e-10;  % Threshold for non-zero values
    result.D3_squared = zeros(size(result.H_mean));
    result.D3_squared(valid_H) = (1/N) * ((params.alpha * (params.pHR + params.pHD)) ./ result.H_mean(valid_H) + ...
                                    params.gamma * params.pIH * result.I_mean(valid_H) ./ (result.H_mean(valid_H).^2));
    result.D3 = sqrt(result.D3_squared);
    
    % D_N^(4)(t) for Recovered (R)
    % (D_N^(4))^2 = (1/N) * (λ_RS/r_N(t) + γ_IR * i_N(t)/(r_N(t))^2 + α_HR * h_N(t)/(r_N(t))^2)
    % Check if r_N(t) is non-zero
    valid_R = result.R_mean > 1e-10;  % Threshold for non-zero values
    result.D4_squared = zeros(size(result.R_mean));
    result.D4_squared(valid_R) = (1/N) * (params.lambda * params.pRS ./ result.R_mean(valid_R) + ...
                                    params.gamma * params.pIR * result.I_mean(valid_R) ./ (result.R_mean(valid_R).^2) + ...
                                    params.alpha * params.pHR * result.H_mean(valid_R) ./ (result.R_mean(valid_R).^2));
    result.D4 = sqrt(result.D4_squared);
    
    % D_N^(5)(t) for Dead (D) - not in the mathematical formulas but included for completeness
    % Using the variance formula from the original code
    % Check if d_N(t) is non-zero
    valid_D = result.D_mean > 1e-10;  % Threshold for non-zero values
    V5 = (1/N) * (params.pID * params.gamma * result.I_mean + params.pHD * params.alpha * result.H_mean);
    result.D5 = zeros(size(result.D_mean));
    result.D5(valid_D) = sqrt(V5(valid_D)) ./ result.D_mean(valid_D);
    
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

function plot_renormalized_results(t, results, det_result, params, R0)
    % Create figure with tiled layout for 5 compartments showing renormalized standard deviations
    figure('Position', [100, 100, 1500, 600]);
    tlayout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot D1 (renormalized standard deviation for susceptible)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Renormalized Std Dev - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('D_N^{(1)}(t)');
    grid on;
    
    % Plot D2 (renormalized standard deviation for infected)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Renormalized Std Dev - Infected (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('D_N^{(2)}(t)');
    grid on;
    
    % Plot D3 (renormalized standard deviation for hospitalized)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Renormalized Std Dev - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('D_N^{(3)}(t)');
    grid on;

    % Plot D4 (renormalized standard deviation for recovered)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D4, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Renormalized Std Dev - Recovered (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('D_N^{(4)}(t)');
    grid on;
    ylim([0, 0.05]); % Set y-axis maximum to 0.05
    
    % Plot D5 (renormalized standard deviation for dead)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D5, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Renormalized Std Dev - Dead (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('D_N^{(5)}(t)');
    grid on;
    
    % Add legend
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add mathematical formula annotation
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', sprintf('D_N^{(l)}(t) = sqrt(V_N^{(l)}(t)) / |m_N^{(l)}(t)|, R_0 = %.2f', R0), ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Save the figure
    saveas(gcf, sprintf('SIHRS_renormalized_std_dev_R0_%.2f.png', R0));
end

function print_renormalized_summary(results, det_result, R0)
    % Print summary statistics for renormalized standard deviations
    fprintf('\n=== RENORMALIZED STANDARD DEVIATION SUMMARY (R0 = %.2f) ===\n', R0);
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
    fprintf('Population Size | Avg D1 | Avg D2 | Avg D3 | Avg D4 | Avg D5\n');
    fprintf('----------------|---------|---------|---------|---------|---------\n');
    for idx = 1:length(results)
        avg_D1 = mean(results{idx}.D1);
        avg_D2 = mean(results{idx}.D2);
        avg_D3 = mean(results{idx}.D3);
        avg_D4 = mean(results{idx}.D4);
        avg_D5 = mean(results{idx}.D5);
        fprintf('%15d | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f\n', ...
            results{idx}.N, avg_D1, avg_D2, avg_D3, avg_D4, avg_D5);
    end
    
    % Verify the conjecture: D_N^(l)(t) > D_Ñ^(l)(t) for N < Ñ
    fprintf('\n=== CONJECTURE VERIFICATION (R0 = %.2f) ===\n', R0);
    fprintf('Checking: D_N^(l)(t) > D_Ñ^(l)(t) for N < Ñ\n');
    for l = 1:4  % Check compartments 1-4 as specified in the conjecture
        fprintf('Compartment %d: ', l);
        conjecture_holds = true;
        for i = 1:length(results)-1
            for j = i+1:length(results)
                if results{i}.N < results{j}.N
                    % Check if D_N^(l) > D_Ñ^(l) on average
                    avg_D_smaller = mean(results{i}.(sprintf('D%d', l)));
                    avg_D_larger = mean(results{j}.(sprintf('D%d', l)));
                    if avg_D_smaller <= avg_D_larger
                        conjecture_holds = false;
                        fprintf('FAILED: N=%d (%.4f) <= N=%d (%.4f) ', ...
                            results{i}.N, avg_D_smaller, results{j}.N, avg_D_larger);
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
sihrs_renormalized_simulation();
