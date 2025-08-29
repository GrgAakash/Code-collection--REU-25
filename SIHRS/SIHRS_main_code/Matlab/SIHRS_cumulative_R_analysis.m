% SIHRS Cumulative G-Renormalized Standard Deviation Analysis
%
% This file calculates the cumulative G-renormalized infinitesimal standard deviation:
% Cumulative_R_N^(l)(t) = ∫₀ᵗ R_N^(l)(τ) dτ
%
% Where R_N^(l)(t) = sqrt(V_N^(l)(t)) / |G_N^(l)(t)|
%
% The cumulative R_N shows the total accumulated relative stochastic noise
% from the beginning of the simulation up to time t. This provides insight
% into how finite size effects build up over the course of the epidemic.

clear all;
close all;

function sihrs_cumulative_R_analysis()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching main analysis
    params.beta = 0.212;           % Infection rate (β > 0)
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
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0)
    params.T = 1000;               % Total simulation time
    params.dt = 0.01;             % Time step for integration
    params.N_values = [1600,3000]; % Multiple population sizes for comparison
    params.initial_s = 0.96;       % Initial susceptible fraction
    params.initial_i = 0.04;       % Initial infected fraction
    params.initial_h = 0;          % Initial hospitalized fraction
    params.initial_r = 0;          % Initial recovered fraction
    params.initial_d = 0;          % Initial dead fraction
    params.n_runs = 1;           % Multiple runs for better statistics (like main simulation)
    params.colors = {'#0072BD','#77AC30', '#A2142F'}; % Colors for different N
    
    % Validate parameters
    validate_params(params);
    
    % Calculate R0 from parameters
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('Calculated R₀ = %.4f from parameters\n', R0);
    
    % Run cumulative analysis
    simulate_and_analyze_cumulative_R(params, R0);
end

function validate_params(params)
    % Validate input parameters (same as main analysis)
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
end

function simulate_and_analyze_cumulative_R(params, R0)
    % Use beta directly from params
    beta = params.beta;
    
    % Precompute time vector
    t = 0:params.dt:params.T;
    
    % Store results
    results = cell(length(params.N_values), 1);
    
    % Run simulations for each population size
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d, R0 = %.2f...\n', params.n_runs, N, R0);
        
        % Run multiple stochastic simulations and calculate cumulative R
        results{idx} = run_cumulative_R_analysis(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_cumulative_R_results(t, results, params, R0);
    print_cumulative_R_summary(results, R0);
end

function result = run_cumulative_R_analysis(N, beta, params, t)
    % Run multiple Gillespie simulations and calculate cumulative R values
    
    % Store cumulative R values for each run
    CumR1_all = zeros(params.n_runs, length(t));
    CumR2_all = zeros(params.n_runs, length(t));
    CumR3_all = zeros(params.n_runs, length(t));
    CumR4_all = zeros(params.n_runs, length(t));
    CumR5_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        
        % Calculate instantaneous G-renormalized values FOR THIS RUN
        
        % R_N^(1)(t) for Susceptible (S)
        V1 = (1/N) * (beta * params.pSI * S_interp .* I_interp + params.lambda * params.pRS * R_interp);
        G1 = -beta * S_interp .* I_interp + params.lambda * R_interp;
        R1_run = sqrt(V1) ./ abs(G1);
        
        % R_N^(2)(t) for Infected (I)
        V2 = (1/N) * (params.gamma * I_interp * (params.pIH + params.pIR + params.pID) + beta * params.pSI * S_interp .* I_interp);
        G2 = -(params.gamma * (params.pIH + params.pIR + params.pID)) * I_interp + beta * S_interp .* I_interp;
        V2_residual = V2 + (1/N) * 1e-12;  % Small residual for extinction handling
        R2_run = sqrt(V2_residual) ./ abs(G2);
        
        % R_N^(3)(t) for Hospitalized (H)
        V3 = (1/N) * (params.gamma * I_interp * params.pIH + params.alpha * H_interp * (params.pHR + params.pHD));
        G3 = params.gamma * params.pIH * I_interp - params.alpha * (params.pHR + params.pHD) * H_interp;
        R3_run = sqrt(V3) ./ abs(G3);
        
        % R_N^(4)(t) for Recovered (R)
        V4 = (1/N) * (params.gamma * I_interp * params.pIR + params.alpha * H_interp * params.pHR + params.lambda * R_interp * params.pRS);
        G4 = params.gamma * params.pIR * I_interp + params.alpha * params.pHR * H_interp - params.lambda * R_interp;
        R4_run = sqrt(V4) ./ abs(G4);
        
        % R_N^(5)(t) for Dead (D)
        V5 = (1/N) * (params.gamma * I_interp * params.pID + params.alpha * H_interp * params.pHD);
        G5 = params.gamma * params.pID * I_interp + params.alpha * params.pHD * H_interp;
        R5_run = sqrt(V5) ./ abs(G5);
        
        % Handle NaN values (set to 0)
        R1_run(isnan(R1_run)) = 0;
        R2_run(isnan(R2_run)) = 0;
        R3_run(isnan(R3_run)) = 0;
        R4_run(isnan(R4_run)) = 0;
        R5_run(isnan(R5_run)) = 0;
        
        % For cumulative calculation: handle both infinities AND large finite values as blow-ups
        % Once any R2 becomes infinite OR very large, cumulative should be infinite thereafter
        
        % Store original infinity locations
        inf_R1 = isinf(R1_run);
        inf_R2 = isinf(R2_run);
        inf_R3 = isinf(R3_run);
        inf_R4 = isinf(R4_run);
        inf_R5 = isinf(R5_run);
        
        % Also detect large finite values as blow-ups (critical threshold blow-ups)
        blowup_threshold = 1e6;  % Any R2 > 1e6 is considered a blow-up
        large_R1 = R1_run > blowup_threshold;
        large_R2 = R2_run > blowup_threshold;
        large_R3 = R3_run > blowup_threshold;
        large_R4 = R4_run > blowup_threshold;
        large_R5 = R5_run > blowup_threshold;
        
        % Combine infinity and large value blow-ups
        blowup_R1 = inf_R1 | large_R1;
        blowup_R2 = inf_R2 | large_R2;
        blowup_R3 = inf_R3 | large_R3;
        blowup_R4 = inf_R4 | large_R4;
        blowup_R5 = inf_R5 | large_R5;
        
        % Temporarily replace blow-ups with large finite values for integration
        temp_max_value = 1e6;
        R1_temp = R1_run; R1_temp(blowup_R1) = temp_max_value;
        R2_temp = R2_run; R2_temp(blowup_R2) = temp_max_value;
        R3_temp = R3_run; R3_temp(blowup_R3) = temp_max_value;
        R4_temp = R4_run; R4_temp(blowup_R4) = temp_max_value;
        R5_temp = R5_run; R5_temp(blowup_R5) = temp_max_value;
        
        % Calculate cumulative integrals
        CumR1_run = cumtrapz(t, R1_temp);
        CumR2_run = cumtrapz(t, R2_temp);
        CumR3_run = cumtrapz(t, R3_temp);
        CumR4_run = cumtrapz(t, R4_temp);
        CumR5_run = cumtrapz(t, R5_temp);
        
        % Handle blow-ups: once any R_N hits blow-up threshold, cumulative should be infinity thereafter
        if any(blowup_R1)
            first_blowup_idx = find(blowup_R1, 1);
            CumR1_run(first_blowup_idx:end) = inf;
        end
        if any(blowup_R2)
            first_blowup_idx = find(blowup_R2, 1);
            CumR2_run(first_blowup_idx:end) = inf;
        end
        if any(blowup_R3)
            first_blowup_idx = find(blowup_R3, 1);
            CumR3_run(first_blowup_idx:end) = inf;
        end
        if any(blowup_R4)
            first_blowup_idx = find(blowup_R4, 1);
            CumR4_run(first_blowup_idx:end) = inf;
        end
        if any(blowup_R5)
            first_blowup_idx = find(blowup_R5, 1);
            CumR5_run(first_blowup_idx:end) = inf;
        end
        
        CumR1_all(run, :) = CumR1_run;
        CumR2_all(run, :) = CumR2_run;
        CumR3_all(run, :) = CumR3_run;
        CumR4_all(run, :) = CumR4_run;
        CumR5_all(run, :) = CumR5_run;
    end
    
    % Use maximum to preserve blow-ups (same logic as main simulation)
    result.CumR1 = max(CumR1_all, [], 1);
    result.CumR2 = max(CumR2_all, [], 1);
    result.CumR3 = max(CumR3_all, [], 1);
    result.CumR4 = max(CumR4_all, [], 1);
    result.CumR5 = max(CumR5_all, [], 1);
    
    % Store final cumulative values for summary
    result.final_CumR1 = result.CumR1(end);
    result.final_CumR2 = result.CumR2(end);
    result.final_CumR3 = result.CumR3(end);
    result.final_CumR4 = result.CumR4(end);
    result.final_CumR5 = result.CumR5(end);
    
    result.N = N;
end

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params)
    % Gillespie simulation (same as main analysis)
    % Initialize populations
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
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
    
    % Gillespie algorithm
    while current_time < params.T
        % Calculate event rates
        si_rate = (beta / N) * S * I * params.pSI;
        ir_rate = params.gamma * I * params.pIR;
        ih_rate = params.gamma * I * params.pIH;
        id_rate = params.gamma * I * params.pID;
        hr_rate = params.alpha * H * params.pHR;
        hd_rate = params.alpha * H * params.pHD;
        rs_rate = params.lambda * R * params.pRS;
        
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
function plot_cumulative_R_results(t, results, params, R0)
    % Create plots for cumulative G-renormalized standard deviations
    
    % Create combined plot for all cumulative R values
    figure('Position', [100, 100, 1500, 800]);
    tlayout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot Cumulative R1 (Susceptible)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(1)} - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(1)}(τ) dτ');
    grid on;
    
    % Plot Cumulative R2 (Infected) - Linear Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        CumR2_plot = results{idx}.CumR2;
        % Handle infinities for plotting - cap at y-axis max to create "fake" blow-up effect
        if any(isinf(CumR2_plot))
            % Find the first blow-up point
            first_inf_idx = find(isinf(CumR2_plot), 1);
            
            if first_inf_idx > 1
                % Set y_max to 1×10^7 to create "fake" blow-up ceiling
                y_max = 1e7;
                
                % Replace all infinite values with y_max (visual ceiling)
                CumR2_plot(isinf(CumR2_plot)) = y_max;
                ylim([0, y_max]);
                
                % Add visual indication of blow-up
                blow_up_time = t(first_inf_idx);
                fprintf('Infected blow-up detected at t = %.2f for N = %d\n', blow_up_time, results{idx}.N);
            else
                % If blow-up happens immediately, use the same ceiling
                y_max = 1e7;
                CumR2_plot(isinf(CumR2_plot)) = y_max;
                ylim([0, y_max]);
            end
        end
        plot(t, CumR2_plot, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(2)} - Infected (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(2)}(τ) dτ');
    grid on;
    
    % Plot Cumulative R3 (Hospitalized)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(3)} - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(3)}(τ) dτ');
    grid on;
    
    % Plot Cumulative R4 (Recovered)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR4, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(4)} - Recovered (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(4)}(τ) dτ');
    grid on;
    
    % Plot Cumulative R5 (Dead)
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR5, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(5)} - Dead (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(5)}(τ) dτ');
    grid on;
    
    % Add legend
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add mathematical formula annotation
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', sprintf('Cumulative R_N^{(l)}(t) = ∫₀ᵗ R_N^{(l)}(τ) dτ, R_0 = %.2f', R0), ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Save the combined figure
    saveas(gcf, sprintf('SIHRS_Cumulative_R_R0_%.2f_combined.png', R0));
    
    % Create LOG SCALE version of cumulative R analysis for dramatic blow-up visualization
    figure('Position', [200, 200, 1500, 800]);
    tlayout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot Cumulative R1 (Susceptible) - Log Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(1)} - Susceptible (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(1)}(τ) dτ');
    set(gca, 'YScale', 'log');
    ylim([0.1, 1000]);  % Adjust range for log scale
    grid on;
    
    % Plot Cumulative R2 (Infected) - Log Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        CumR2_plot = results{idx}.CumR2;
        % Handle infinities for plotting - cap at y-axis max for log scale
        if any(isinf(CumR2_plot))
            % Find the first blow-up point
            first_inf_idx = find(isinf(CumR2_plot), 1);
            
            if first_inf_idx > 1
                % Set y_max to 1×10^7 for log scale visualization
                y_max = 1e7;
                
                % Replace all infinite values with y_max (visual ceiling)
                CumR2_plot(isinf(CumR2_plot)) = y_max;
                
                % Add visual indication of blow-up
                blow_up_time = t(first_inf_idx);
                fprintf('Infected cumulative blow-up detected at t = %.2f for N = %d (Log Scale)\n', blow_up_time, results{idx}.N);
            else
                % If blow-up happens immediately, use the same ceiling
                y_max = 1e7;
                CumR2_plot(isinf(CumR2_plot)) = y_max;
            end
        end
        plot(t, CumR2_plot, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(2)} - Infected (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(2)}(τ) dτ');
    
    % Set log scale for dramatic blow-up visualization
    set(gca, 'YScale', 'log');
    ylim([1, 1e7]);  % Log scale can't start at 0
    
    grid on;
    
    % Plot Cumulative R3 (Hospitalized) - Log Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(3)} - Hospitalized (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(3)}(τ) dτ');
    set(gca, 'YScale', 'log');
    ylim([0.1, 1000]);  % Adjust range for log scale
    grid on;
    
    % Plot Cumulative R4 (Recovered) - Log Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR4, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(4)} - Recovered (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(4)}(τ) dτ');
    set(gca, 'YScale', 'log');
    ylim([0.1, 1000]);  % Adjust range for log scale
    grid on;
    
    % Plot Cumulative R5 (Dead) - Log Scale
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.CumR5, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(5)} - Dead (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(5)}(τ) dτ');
    set(gca, 'YScale', 'log');
    ylim([0.1, 1000]);  % Adjust range for log scale
    grid on;
    
    % Add legend
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add mathematical formula annotation for log scale
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', sprintf('Cumulative R_N^{(l)}(t) = ∫₀ᵗ R_N^{(l)}(τ) dτ, R_0 = %.2f (Log Scale)', R0), ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Save the combined log scale figure
    saveas(gcf, sprintf('SIHRS_Cumulative_R_R0_%.2f_combined_LOG.png', R0));
    
    % Create individual log-scale plots for key compartments
    % Individual Cumulative R2 (Infected) - Log Scale for dramatic blow-up visualization
    figure('Position', [300, 300, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        CumR2_plot = results{idx}.CumR2;
        % Handle infinities for plotting
        if any(isinf(CumR2_plot))
            CumR2_plot(isinf(CumR2_plot)) = 1e7;  % Cap at y_max for visualization
        end
        plot(t, CumR2_plot, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Cumulative R_N^{(2)} - Infected (R_0 = %.2f) - Log Scale', R0));
    xlabel('Time'); ylabel('∫₀ᵗ R_N^{(2)}(τ) dτ');
    
    % Set log scale for dramatic blow-up visualization
    set(gca, 'YScale', 'log');
    ylim([1, 1e7]);  % Log scale can't start at 0
    
    grid on;
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_Cumulative_R2_Infected_Log_R0_%.2f.png', R0));
    
    fprintf('Generated cumulative R analysis plots:\n');
    fprintf('  - SIHRS_Cumulative_R_R0_%.2f_combined.png (Linear Scale)\n', R0);
    fprintf('  - SIHRS_Cumulative_R_R0_%.2f_combined_LOG.png (Log Scale)\n', R0);
    fprintf('  - SIHRS_Cumulative_R2_Infected_Log_R0_%.2f.png (Individual Log Scale)\n', R0);
end

function print_cumulative_R_summary(results, R0)
    % Print summary statistics for cumulative G-renormalized standard deviations
    fprintf('\n=== CUMULATIVE G-RENORMALIZED STANDARD DEVIATION SUMMARY (R0 = %.2f) ===\n', R0);
    fprintf('Population Size | Final CumR1 | Final CumR2 | Final CumR3 | Final CumR4 | Final CumR5\n');
    fprintf('----------------|-------------|-------------|-------------|-------------|-------------\n');
    for idx = 1:length(results)
        fprintf('%15d | %11.2f | %11.2f | %11.2f | %11.2f | %11.2f\n', ...
            results{idx}.N, results{idx}.final_CumR1, results{idx}.final_CumR2, ...
            results{idx}.final_CumR3, results{idx}.final_CumR4, results{idx}.final_CumR5);
    end
    
    % Analyze finite size scaling of cumulative effects
    fprintf('\n=== FINITE SIZE SCALING ANALYSIS ===\n');
    fprintf('Cumulative R values should scale as 1/sqrt(N) for finite size effects\n');
    if length(results) >= 2
        for l = 1:5
            fprintf('Compartment %d scaling:\n', l);
            for i = 1:length(results)-1
                for j = i+1:length(results)
                    N1 = results{i}.N;
                    N2 = results{j}.N;
                    CumR1 = results{i}.(sprintf('final_CumR%d', l));
                    CumR2 = results{j}.(sprintf('final_CumR%d', l));
                    
                    theoretical_ratio = sqrt(N2/N1);
                    actual_ratio = CumR1/CumR2;
                    
                    fprintf('  N=%d vs N=%d: Theoretical ratio=%.3f, Actual ratio=%.3f\n', ...
                        N1, N2, theoretical_ratio, actual_ratio);
                end
            end
        end
    end
end
% Run the cumulative analysis
sihrs_cumulative_R_analysis();


