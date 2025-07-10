
% Model Parameters:
% - N: Population size
% - beta: Infection rate (R0 = beta when p1/p2 = 1, gamma = 1)
% - p1: Infection probability
% - p2: Recovery probability
% - gamma: Recovery rate
% - initial_s, initial_i, initial_r: Initial fractions (s + i + r = 1)
% - T: Maximum simulation time
% - dt: Time step for interpolation
%
% Variance Formulas:
% - V1 = (1/N) * p1 * beta * S_prop * I_prop (susceptible)
% - V2 = (1/N) * (p1 * beta * S_prop * I_prop + p2 * gamma * I_prop) (infected)
% - V3 = (1/N) * p2 * gamma * I_prop (recovered)
clear all;
close all;
function sir_simulation_improved()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters
    params.p1 = 0.5;              % Infection probability
    params.p2 = 0.5;              % Recovery probability
    params.gamma = 1;             % Recovery rate
    params.T = 23;                % Total simulation time (from Code 1)
    params.dt = 0.01;             % Time step for interpolation
    params.N_values = round([10^2.5, 10^3, 10^3.5, 10^4]); % Population sizes
    params.initial_s = 0.96;      % Initial susceptible fraction
    params.initial_i = 0.04;      % Initial infected fraction
    params.initial_r = 0;         % Initial recovered fraction
    params.n_runs = 40;          % Number of stochastic runs
    params.colors = {'#0000FF', '#FF0000', '#008000', '#FF00FF'}; % Colors for plotting
    params.R0_values = [0.95, 1.3]; % R0 values
    
    % Validate parameters
    validate_params(params);
    
    % Run simulations for each R0
    for r_idx = 1:length(params.R0_values)
        R0 = params.R0_values(r_idx);
        simulate_and_analyze(params, R0);
    end
end

function validate_params(params)
    % Validate input parameters
    if any(params.N_values <= 0) || any(mod(params.N_values, 1) ~= 0)
        error('N_values must be positive integers');
    end
    if params.p1 <= 0 || params.p1 > 1 || params.p2 <= 0 || params.p2 > 1
        error('p1 and p2 must be in (0, 1]');
    end
    if params.gamma <= 0
        error('gamma must be positive');
    end
    if abs(params.initial_s + params.initial_i + params.initial_r - 1) > 1e-6
        error('Initial fractions must sum to 1');
    end
    if params.T <= 0 || params.dt <= 0
        error('T and dt must be positive');
    end
    if params.n_runs <= 0 || mod(params.n_runs, 1) ~= 0
        error('n_runs must be a positive integer');
    end
end

function simulate_and_analyze(params, R0)
    % Calculate beta from R0
    beta = R0 * params.gamma;
    
    % Precompute time vector
    t = 0:params.dt:params.T;
    
    % Store results
    results = cell(length(params.N_values), 1);
    det_result = solve_deterministic_sir(beta, params);
    
    % Run simulations for each population size
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d, R0 = %.2f...\n', params.n_runs, N, R0);
        
        % Run multiple stochastic simulations
        results{idx} = run_multiple_gillespie(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_results(t, results, det_result, params, R0);
    print_summary(results, det_result, R0);
end

function result = run_multiple_gillespie(N, beta, params, t)
    % Run multiple Gillespie simulations and aggregate results
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, R_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, R_interp] = interpolate_results(S_hist, I_hist, R_hist, time_pts, t, N);
        S_all(run, :) = S_interp;
        I_all(run, :) = I_interp;
        R_all(run, :) = R_interp;
    end
    
    % Compute mean and standard deviation
    result.S_mean = mean(S_all, 1);
    result.I_mean = mean(I_all, 1);
    result.R_mean = mean(R_all, 1);
    
    % Compute infinitesimal standard deviations
    V1 = (1/N) * params.p1 * beta * result.S_mean .* result.I_mean;
    V2 = (1/N) * (params.p1 * beta * result.S_mean .* result.I_mean + params.p2 * params.gamma * result.I_mean);
    V3 = (1/N) * params.p2 * params.gamma * result.I_mean;
    
    result.sigma1 = sqrt(V1);
    result.sigma2 = sqrt(V2);
    result.sigma3 = sqrt(V3);
    result.N = N;
    result.peak_infected = max(result.I_mean * N);
    result.peak_time = t(find(result.I_mean == max(result.I_mean), 1, 'first'));
    result.final_time = t(end);
end

function [S_hist, I_hist, R_hist, time_pts] = gillespie_sim(N, beta, params)
    % Initialize populations
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    R = round(N * params.initial_r);
    
    % Verify population conservation
    if abs(S + I + R - N) > 1
        error('Initial populations do not sum to N');
    end
    
    % Preallocate arrays
    max_events = round(10 * params.T * (beta + params.gamma) * N);
    S_hist = zeros(1, max_events);
    I_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events);
    time_pts = zeros(1, max_events);
    
    S_hist(1) = S;
    I_hist(1) = I;
    R_hist(1) = R;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    % Gillespie algorithm
    while current_time < params.T && I > 0
        a1 = (beta / N) * S * I * params.p1;
        a2 = params.gamma * I * params.p2;
        a0 = a1 + a2;
        
        if a0 <= 0
            break;
        end
        
        tau = -log(rand) / a0;
        current_time = current_time + tau;
        
        if current_time > params.T
            break;
        end
        
        if rand < a1 / a0
            if S > 0
                S = S - 1;
                I = I + 1;
            end
        else
            if I > 0
                I = I - 1;
                R = R + 1;
            end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S;
        I_hist(event_count) = I;
        R_hist(event_count) = R;
        time_pts(event_count) = current_time;
    end
    
    % Trim arrays
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    time_pts = time_pts(1:event_count);
    
    % Verify population conservation
    assert(all(S_hist + I_hist + R_hist == N), 'Population not conserved');
end

function [S_interp, I_interp, R_interp] = interpolate_results(S_hist, I_hist, R_hist, time_pts, t, N)
    % Interpolate to fixed time grid using 'previous'
    try
        S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
        I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
        R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
        
        % Handle values beyond last event
        S_interp(t > max(time_pts)) = S_hist(end) / N;
        I_interp(t > max(time_pts)) = I_hist(end) / N;
        R_interp(t > max(time_pts)) = R_hist(end) / N;
    catch e
        error('Interpolation failed: %s', e.message);
    end
end

function det_result = solve_deterministic_sir(beta, params)
    % Solve deterministic SIR model using ODE45
    try
        tspan = [0, params.T];
        y0 = [params.initial_s; params.initial_i; params.initial_r];
        ode_system = @(t, y) [
            -params.p1 * beta * y(2) * y(1);
            params.p1 * beta * y(2) * y(1) - params.p2 * params.gamma * y(2);
            params.p2 * params.gamma * y(2)
        ];
        
        [T, Y] = ode45(ode_system, tspan, y0);
        
        det_result.T = T;
        det_result.S_prop = Y(:, 1);
        det_result.I_prop = Y(:, 2);
        det_result.R_prop = Y(:, 3);
        [det_result.peak_infected_prop, idx] = max(det_result.I_prop);
        det_result.peak_time = T(idx);
        det_result.final_time = T(end);
    catch e
        error('ODE solver failed: %s', e.message);
    end
end

function plot_results(t, results, det_result, params, R0)
    % Create figure with tiled layout
    figure('Position', [100, 100, 1200, 400]);
    tlayout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot sigma1 (susceptible)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.sigma1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Susceptible Std Dev (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('\sigma_N^{(1)}(t)');
    grid on;
    
    % Plot sigma2 (infected)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.sigma2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Infected Std Dev (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('\sigma_N^{(2)}(t)');
    grid on;
    
    % Plot sigma3 (recovered)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.sigma3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Recovered Std Dev (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('\sigma_N^{(3)}(t)');
    grid on;
    
    % Add legend
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add beta and gamma annotation
    annotation('textbox', [0.05, 0.02, 0.1, 0.05], ...
        'String', sprintf('β = %.2f, γ = %.1f', R0, params.gamma), ...
        'FontSize', 10, 'EdgeColor', 'none');
end

function print_summary(results, det_result, R0)
    % Print summary statistics
    fprintf('\n=== SIMULATION SUMMARY (R0 = %.2f) ===\n', R0);
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
    
    % Convergence analysis
    fprintf('\n=== CONVERGENCE ANALYSIS (R0 = %.2f) ===\n', R0);
    fprintf('Population Size | Peak Infected (Normalized) | Difference from Deterministic\n');
    fprintf('----------------|----------------------------|------------------------------\n');
    for idx = 1:length(results)
        normalized_peak = results{idx}.peak_infected / results{idx}.N;
        diff_from_det = abs(normalized_peak - det_result.peak_infected_prop);
        fprintf('%15d | %26.4f | %28.4f\n', ...
            results{idx}.N, normalized_peak, diff_from_det);
    end
end

% Run the simulation
sir_simulation_improved();