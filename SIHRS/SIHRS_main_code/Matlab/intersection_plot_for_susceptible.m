% Intersection Plot for Susceptible Compartment (R1)
% This script creates individual plots for each N showing S(t) vs Theoretical Threshold
% Clear visualization of when blow-up conditions are met for the susceptible compartment
% SAVED VERSION - Do not delete!

clear all;
close all;

% Set random seed for reproducibility
rng(1);

% Parameters for SIHRS with death
beta = 0.212;           % Infection rate
pSI = 1.0;              % Infection probability (S to I)
pII = 0.0;              % probability of I to I (stay infected)
pIH = 0.1060;           % probability of I to H
pIR = 0.8921;           % probability of I to R
pID = 0.0019;           % probability of I to D
pHH = 0.00;             % probability of H to H (stay hospitalized)
pHR = 0.846;            % probability of H to R
pHD = 0.154;            % probability of H to D
pRR = 0.02;             % probability of R to R (stay recovered)
pRS = 0.98;             % probability of R to S
gamma = 0.100346667;    % Adjusted for exact critical hit at S = 142/300
alpha = 0.1;            % Hospitalized transition rate
lambda = 0.0083;        % Recovered to susceptible rate (lambda)
T = 200;               % Total simulation time
dt = 0.01;              % Time step for integration
N_values = [1600, 3000]; % Multiple population sizes for comparison

% Calculate R0 from parameters
R0 = pSI * beta / (gamma * (1 - pII));
fprintf('Calculated R0 = %.4f from parameters\n', R0);

% Precompute time vector
t = 0:dt:T;

% Store results
results = cell(length(N_values), 1);

% Run simulations for each population size
for idx = 1:length(N_values)
    N = N_values(idx);
    fprintf('Running simulation for N = %d, R0 = %.2f...\n', N, R0);
    
    % Run Gillespie simulation
    [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta);
    
    % Interpolate to fixed time grid
    [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
    
    % Calculate the theoretical blow-up condition: (lambda * p_RS * r_N(t)) / (beta * p_SI * i_N(t))
    % This is equivalent to: (lambda * p_RS * R(t)) / (beta * p_SI * I(t))
    theoretical_threshold = (lambda * pRS * R_interp) ./ (beta * pSI * I_interp);
    
    % Find intersection points (where S(t) ≈ theoretical threshold)
    tolerance = 0.001;  % 0.001% tolerance for intersection detection
    intersection_points = abs(S_interp - theoretical_threshold) <= tolerance;
    intersection_indices = find(intersection_points);
    
    fprintf('  Found %d intersection points\n', length(intersection_indices));
    
    % Store results
    results{idx}.N = N;
    results{idx}.S_interp = S_interp;
    results{idx}.I_interp = I_interp;
    results{idx}.R_interp = R_interp;
    results{idx}.theoretical_threshold = theoretical_threshold;
    results{idx}.intersection_points = intersection_points;
    results{idx}.intersection_indices = intersection_indices;
    
    fprintf('Completed N = %d\n', N);
end

% Create separate plots for each population size
for idx = 1:length(N_values)
    N = results{idx}.N;
    
    % Create individual figure for this N - MAIN PLOT ONLY
    figure('Position', [100 + 100*idx, 100 + 100*idx, 800, 600]);
    
    % Main intersection plot - ONLY THIS ONE
    hold on;
    % Plot S(t) as solid line
    plot(t, results{idx}.S_interp, 'Color', [0.2, 0.5, 0.8], 'LineWidth', 2, 'DisplayName', 'S(t)');
    % Plot theoretical threshold as dashed line
    plot(t, results{idx}.theoretical_threshold, '--', 'Color', [0.8, 0.2, 0.5], 'LineWidth', 2, 'DisplayName', 'Theoretical Threshold');
    
    % Highlight intersection points
    if ~isempty(results{idx}.intersection_indices)
        plot(t(results{idx}.intersection_indices), results{idx}.S_interp(results{idx}.intersection_indices), 'o', ...
            'Color', [0.9, 0.6, 0.1], 'MarkerSize', 8, 'MarkerFaceColor', [0.9, 0.6, 0.1], ...
            'DisplayName', 'Intersections');
    end
    
    title(sprintf('N = %d: S(t) vs Theoretical Threshold', N));
    subtitle(sprintf('(lambda * p_{RS} * r_N(t)) / (beta * p_SI * i_N(t)), R0 = %.2f', R0));
    xlabel('Time'); ylabel('Population Fraction');
    grid on;
    legend('Location', 'best');
    
    % Save individual figure with script name prefix
    saveas(gcf, sprintf('intersection_plot_for_susceptible_N%d_R0_%.2f.png', N, R0));
    
    fprintf('Created separate plot for N = %d\n', N);
end





title('Combined: S(t) vs Theoretical Threshold for All N Values');
subtitle(sprintf('(lambda * p_RS * r_N(t)) / (beta * p_SI * i_N(t)), R0 = %.2f', R0));
xlabel('Time'); ylabel('Population Fraction');
grid on;
legend('Location', 'best');



% Print summary
fprintf('\n=== INDIVIDUAL INTERSECTION PLOTS SUMMARY ===\n');
fprintf('Created individual plots for each population size N\n');
fprintf('Each plot shows:\n');
fprintf('  1. S(t) vs Theoretical Threshold (lambda * p_RS * r_N(t)) / (beta * p_SI * i_N(t))\n');
fprintf('  2. Difference |S(t) - Theoretical Threshold|\n');
fprintf('  3. Zoomed intersection analysis\n');
fprintf('  4. Threshold components: lambda * p_RS * R(t) and beta * p_SI * I(t)\n\n');

for idx = 1:length(N_values)
    N = results{idx}.N;
    fprintf('N = %d:\n', N);
    fprintf('  Total intersection points: %d\n', length(results{idx}.intersection_indices));
    if ~isempty(results{idx}.intersection_indices)
        fprintf('  First intersection at t = %.2f\n', t(results{idx}.intersection_indices(1)));
        fprintf('  Last intersection at t = %.2f\n', t(results{idx}.intersection_indices(end)));
    end
    fprintf('  Plot saved as: intersection_plot_for_susceptible_N%d_R0_%.2f.png\n', N, R0);
    fprintf('\n');
end


fprintf('\n=== SCRIPT SAVED PERMANENTLY ===\n');
fprintf('This script is now saved as intersection_plot_for_susceptible.m\n');
fprintf('You can run it anytime to regenerate the susceptible intersection plots!\n');

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta)
    % Initialize populations
    S = round(N * 0.96);
    I = round(N * 0.04);
    H = 0;
    R = 0;
    D = 0;
    
    % Preallocate arrays
    max_events = round(10 * 1000 * (beta + 0.100346667 + 0.1 + 0.0083) * N);
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
    
    % Gillespie algorithm
    while current_time < 1000
        % Calculate event rates
        si_rate = (beta / N) * S * I * 1.0;
        ir_rate = 0.100346667 * I * 0.959;
        ih_rate = 0.100346667 * I * 0.04;
        id_rate = 0.100346667 * I * 0.001;
        hr_rate = 0.1 * H * 0.9882;
        hd_rate = 0.1 * H * 0.0018;
        rs_rate = 0.0083 * R * 0.98;
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        
        if total_rate == 0
            break;
        end
        
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        
        if current_time > 1000
            break;
        end
        
        % Determine which event occurs
        chance = rand * total_rate;
        if chance < si_rate
            if S > 0; S = S - 1; I = I + 1; end
        elseif chance < (si_rate + ir_rate)
            if I > 0; I = I - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate)
            if I > 0; I = I - 1; H = H + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            if I > 0; I = I - 1; D = D + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            if H > 0; H = H - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            if H > 0; H = H - 1; D = D + 1; end
        else
            if R > 0; R = R - 1; S = S + 1; end
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
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    % Interpolate to fixed time grid
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
end
