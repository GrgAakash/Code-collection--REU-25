% SIHRS with Death Model - D_N Renormalized Infinitesimal Standard Deviation
% 
% This file calculates the D_N renormalized infinitesimal standard deviation:
% D_N^(l)(t) = sqrt(V_N^(l)(t)) / |m_N^(l)(t)|
% 
% Where V_N^(l)(t) are the variance formulas for each compartment l:
% l = 1: Susceptible (S) - m_N^(1)(t) = s_N(t)
% l = 2: Infected (I) - m_N^(2)(t) = i_N(t) 
% l = 3: Hospitalized (H) - m_N^(3)(t) = h_N(t)
% l = 4: Recovered (R) - m_N^(4)(t) = r_N(t)
% l = 5: Dead (D) - m_N^(5)(t) = d_N(t)
%
% Variance formulas:
% V_N^(1)(t) = (1/N) * (β_SI * s_N(t) * i_N(t) + λ_RS * r_N(t))
% V_N^(2)(t) = (1/N) * ((γ_IH + γ_IR + γ_ID) * i_N(t) + β_SI * s_N(t) * i_N(t))
% V_N^(3)(t) = (1/N) * (γ_IH * i_N(t) + h_N(t) * (α_HR + α_HD))
% V_N^(4)(t) = (1/N) * (γ_IR * i_N(t) + α_HR * h_N(t) + λ_RS * r_N(t))
% V_N^(5)(t) = (1/N) * (γ_ID * i_N(t) + α_HD * h_N(t))

clear all;
close all;

function sihrs_DN_renormalized_simulation()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching SIHRS.m
    params.beta = 0.212;           % Infection rate (β > 0)
    params.pSI = 1.0;              % Infection probability (S to I)
    params.pII = 0.0;              % probability of I to I (stay infected)
    params.pIH = 0.1060;           % probability of I to H
    params.pIR = 0.8921;           % probability of I to R
    params.pID = 0.0019;           % probability of I to D
    params.pHH = 0.00;             % probability of H to H (stay hospitalized)
    params.pHR = 0.846;            % probability of H to R
    params.pHD = 0.154;            % probability of H to D
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
    
    % Calculate R0 from parameters (same formula as SIHRS.m)
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('Calculated R₀ = %.4f from parameters\n', R0);
    
    % Run simulation with calculated R0
    simulate_and_analyze_DN_renormalized(params, R0);
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

function simulate_and_analyze_DN_renormalized(params, R0)
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
        results{idx} = run_multiple_gillespie_DN_renormalized(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_DN_renormalized_results(t, results, det_result, params, R0);
    print_DN_renormalized_summary(results, det_result, R0);
end

function result = run_multiple_gillespie_DN_renormalized(N, beta, params, t)
    % Run multiple Gillespie simulations and aggregate results
    % Calculate D_N renormalized values for each run first, then average
    
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    H_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    D_all = zeros(params.n_runs, length(t));
    
    % Store D_N renormalized values for each run
    D1_all = zeros(params.n_runs, length(t));
    D2_all = zeros(params.n_runs, length(t));
    D3_all = zeros(params.n_runs, length(t));
    D4_all = zeros(params.n_runs, length(t));
    D5_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        
        % Store population trajectories
        S_all(run, :) = S_interp;
        I_all(run, :) = I_interp;
        H_all(run, :) = H_interp;
        R_all(run, :) = R_interp;
        D_all(run, :) = D_interp;
        
        % Calculate D_N renormalized values FOR THIS INDIVIDUAL RUN
        % This preserves blow-up behavior before averaging
        
        % D_N^(1)(t) for Susceptible (S) - m_N^(1)(t) = s_N(t)
        % V_N^(1)(t) = (1/N) * (β_SI * s_N(t) * i_N(t) + λ_RS * r_N(t))
        V1 = (1/N) * (beta * params.pSI * S_interp .* I_interp + params.lambda * params.pRS * R_interp);
        D1_run = sqrt(V1) ./ S_interp;  % D_N^(1) = sqrt(V1) / s_N(t)
        D1_all(run, :) = D1_run;
        
        % D_N^(2)(t) for Infected (I) - m_N^(2)(t) = i_N(t)
        % V_N^(2)(t) = (1/N) * ((γ_IH + γ_IR + γ_ID) * i_N(t) + β_SI * s_N(t) * i_N(t))
        V2 = (1/N) * ((params.gamma * (params.pIH + params.pIR + params.pID)) * I_interp + beta * params.pSI * S_interp .* I_interp);
        D2_run = sqrt(V2) ./ I_interp;  % D_N^(2) = sqrt(V2) / i_N(t)
        D2_all(run, :) = D2_run;
        
        % D_N^(3)(t) for Hospitalized (H) - m_N^(3)(t) = h_N(t)
        % V_N^(3)(t) = (1/N) * (γ_IH * i_N(t) + h_N(t) * (α_HR + α_HD))
        V3 = (1/N) * (params.gamma * params.pIH * I_interp + H_interp * (params.alpha * (params.pHR + params.pHD)));
        D3_run = sqrt(V3) ./ H_interp;  % D_N^(3) = sqrt(V3) / h_N(t)
        D3_all(run, :) = D3_run;
        
        % D_N^(4)(t) for Recovered (R) - m_N^(4)(t) = r_N(t)
        % V_N^(4)(t) = (1/N) * (γ_IR * i_N(t) + α_HR * h_N(t) + λ_RS * r_N(t))
        V4 = (1/N) * (params.gamma * params.pIR * I_interp + params.alpha * params.pHR * H_interp + params.lambda * params.pRS * R_interp);
        % Piecewise: D_N^(4)(t) = 0 if r_N(t) = 0, otherwise sqrt(V4)/r_N(t)
        D4_run = zeros(size(R_interp));
        valid_R = R_interp > 0;
        D4_run(valid_R) = sqrt(V4(valid_R)) ./ R_interp(valid_R);
        D4_all(run, :) = D4_run;
        
        % D_N^(5)(t) for Dead (D) - m_N^(5)(t) = d_N(t)
        % V_N^(5)(t) = (1/N) * (γ_ID * i_N(t) + α_HD * h_N(t))
        V5 = (1/N) * (params.gamma * params.pID * I_interp + params.alpha * params.pHD * H_interp);
        % Piecewise: D_N^(5)(t) = 0 if d_N(t) = 0, otherwise sqrt(V5)/d_N(t)
        D5_run = zeros(size(D_interp));
        valid_D = D_interp > 0;
        D5_run(valid_D) = sqrt(V5(valid_D)) ./ D_interp(valid_D);
        D5_all(run, :) = D5_run;
    end
    
    % Compute mean proportions (for summary statistics)
    result.S_mean = mean(S_all, 1);
    result.I_mean = mean(I_all, 1);
    result.H_mean = mean(H_all, 1);
    result.R_mean = mean(R_all, 1);
    result.D_mean = mean(D_all, 1);
    
    % CAPTURE TRUE INFINITY: Set Inf values to very large number for plotting
    % This will show as the maximum y-axis range in plots
    max_display_value = 1e9;  % Very large value to represent infinity
    
    % Store original values for blow-up detection
    result.D1_original = D1_all;
    result.D2_original = D2_all;
    result.D3_original = D3_all;
    result.D4_original = D4_all;
    result.D5_original = D5_all;
    
    % Replace Inf with large value, NaN with 0 for display
    D1_all(isinf(D1_all)) = max_display_value;
    D2_all(isinf(D2_all)) = max_display_value;
    D3_all(isinf(D3_all)) = max_display_value;
    D4_all(isinf(D4_all)) = max_display_value;
    D5_all(isinf(D5_all)) = max_display_value;
    
    D1_all(isnan(D1_all)) = 0;
    D2_all(isnan(D2_all)) = 0;
    D3_all(isnan(D3_all)) = 0;
    D4_all(isnan(D4_all)) = 0;
    D5_all(isnan(D5_all)) = 0;
    
    % Use maximum across runs to capture blow-ups (including infinity)
    result.D1 = max(D1_all, [], 1);
    result.D2 = max(D2_all, [], 1);
    result.D3 = max(D3_all, [], 1);
    result.D4 = max(D4_all, [], 1);
    result.D5 = max(D5_all, [], 1);
    
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
            % If no more events possible, add final state and continue to final time
            event_count = event_count + 1;
            S_hist(event_count) = S;
            I_hist(event_count) = I;
            H_hist(event_count) = H;
            R_hist(event_count) = R;
            D_hist(event_count) = D;
            time_pts(event_count) = params.T;
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
    
    % Debug: Print simulation end time
    fprintf('  Simulation ended at t = %.2f (target: %.2f)\n', time_pts(end), params.T);
    
    % Verify population conservation
    assert(all(S_hist + I_hist + H_hist + R_hist + D_hist == N), 'Population not conserved');
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    % Interpolate to fixed time grid using 'previous'
    try
        % Ensure we have data points at the full time range
        if max(time_pts) < max(t)
            % Add final state at the end of time range
            time_pts = [time_pts, max(t)];
            S_hist = [S_hist, S_hist(end)];
            I_hist = [I_hist, I_hist(end)];
            H_hist = [H_hist, H_hist(end)];
            R_hist = [R_hist, R_hist(end)];
            D_hist = [D_hist, D_hist(end)];
        end
        
        S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
        I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
        H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
        R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
        D_interp = interp1(time_pts, D_hist, t, 'previous') / N;
        
        % Handle values beyond last event (should not be needed now)
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

function plot_DN_renormalized_results(t, results, det_result, params, R0)
    % Create 5 separate figures for each compartment showing D_N renormalized standard deviations
    
    % Plot D1 (D_N renormalized standard deviation for susceptible) - Linear Scale
    figure('Position', [100, 100, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('D_N - Susceptible (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(1)}(t)');
    grid on;
    % Calculate dynamic y-axis maximum: 1.2 times the peak value across all N values
    max_D1 = max(cellfun(@(r) max(r.D1), results));
    ylim([0, max_D1 * 1.2]);
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_DN1_Susceptible_Linear_R0_%.2f.png', R0));
    
    % Plot D2 (D_N renormalized standard deviation for infected) - Linear Scale
    figure('Position', [200, 200, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        D2_infected = results{idx}.D2;
        D2_original = results{idx}.D2_original;
        
        % Find finite values and plot them
        finite_mask = isfinite(D2_infected);
        if any(finite_mask)
            plot(t(finite_mask), D2_infected(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values in original data and draw vertical lines for blow-ups
        nan_inf_mask = isnan(D2_original) | isinf(D2_original);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t < t(start_idx) & finite_mask, 1, 'last');
                    
                    if ~isempty(before_idx)
                        % Draw vertical line from last finite value to infinity
                        last_finite_value = D2_infected(before_idx);
                        plot([t(before_idx), t(before_idx)], [last_finite_value, 1e9], ...
                             'Color', params.colors{idx}, 'LineWidth', 2.0);
                    end
                    
                    % Draw vertical line at the start of the blow-up region
                    plot([t(start_idx), t(start_idx)], [0, 1e9], ...
                         'Color', params.colors{idx}, 'LineWidth', 2.0);
                end
            end
        end
    end
    title(sprintf('D_N - Infected (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(2)}(t)');
    grid on;
    % Set y-axis range: 1.1 times peak for finite values, but allow blow-ups to go to top
    max_D2 = max(cellfun(@(r) max(r.D2(isfinite(r.D2))), results));
    % Use 1.1 times peak value for infected compartment
    ylim([0, max_D2 * 1.1]); % Set to 1.1 times peak value
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_DN2_Infected_Linear_R0_%.2f.png', R0));
    
    % Plot D3 (D_N renormalized standard deviation for hospitalized) - Linear Scale
    figure('Position', [300, 300, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for hospitalized
        t_hosp = t(2:end);
        D3_hosp = results{idx}.D3(2:end);
        
        % Find finite values and plot them
        finite_mask = isfinite(D3_hosp);
        if any(finite_mask)
            plot(t_hosp(finite_mask), D3_hosp(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values and draw horizontal lines connecting vertical spikes
        nan_inf_mask = isnan(D3_hosp) | isinf(D3_hosp);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t_hosp < t_hosp(start_idx) & finite_mask, 1, 'last');
                    % Find the first finite value after this gap
                    after_idx = find(t_hosp > t_hosp(end_idx) & finite_mask, 1, 'first');
                    
                    if ~isempty(before_idx) && ~isempty(after_idx)
                        % Draw horizontal line at the top connecting the spikes
                        plot([t_hosp(before_idx), t_hosp(after_idx)], [1e9, 1e9], ...
                             'Color', params.colors{idx}, 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
    title(sprintf('D_N - Hospitalized (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(3)}(t)');
    grid on;
    % Set y-axis range: 1.2 times peak for finite values, but allow blow-ups to go to top
    max_D3 = max(D3_hosp(finite_mask));
    % Always use 1e9 for hospitalized since it has blow-ups
    ylim([0, 1e9]); % Set to exactly 1e9 so spikes reach the top
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_DN3_Hospitalized_Linear_R0_%.2f.png', R0));
    
    % Plot D4 (D_N renormalized standard deviation for recovered) - Linear Scale
    figure('Position', [400, 400, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for recovered (will blow up since r_N(0) = 0)
        t_rec = t(2:end);
        D4_rec = results{idx}.D4(2:end);
        
        % Find finite values and plot them
        finite_mask = isfinite(D4_rec);
        if any(finite_mask)
            plot(t_rec(finite_mask), D4_rec(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values and draw horizontal lines connecting vertical spikes
        nan_inf_mask = isnan(D4_rec) | isinf(D4_rec);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t_rec < t_rec(start_idx) & finite_mask, 1, 'last');
                    % Find the first finite value after this gap
                    after_idx = find(t_rec > t_rec(end_idx) & finite_mask, 1, 'first');
                    
                    if ~isempty(before_idx) && ~isempty(after_idx)
                        % Draw horizontal line at the top connecting the spikes
                        plot([t_rec(before_idx), t_rec(after_idx)], [1e6, 1e6], ...
                             'Color', params.colors{idx}, 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
    title(sprintf('D_N - Recovered (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(4)}(t)');
    grid on;
    % Set y-axis range to 1.2 times the peak (like Susceptible and Infected)
    max_D4 = max(D4_rec(finite_mask));
    if isempty(max_D4) || isnan(max_D4) || isinf(max_D4)
        ylim([0, 1]); % Default range if no valid data
    else
        ylim([0, max_D4 * 1.2]);
    end
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_DN4_Recovered_Linear_R0_%.2f.png', R0));
    
    % Plot D5 (D_N renormalized standard deviation for dead) - Linear Scale
    figure('Position', [500, 500, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for dead
        t_dead = t(2:end);
        D5_dead = results{idx}.D5(2:end);
        
        % Find finite values and plot them
        finite_mask = isfinite(D5_dead);
        if any(finite_mask)
            plot(t_dead(finite_mask), D5_dead(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values and draw horizontal lines connecting vertical spikes
        nan_inf_mask = isnan(D5_dead) | isinf(D5_dead);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t_dead < t_dead(start_idx) & finite_mask, 1, 'last');
                    % Find the first finite value after this gap
                    after_idx = find(t_dead > t_dead(end_idx) & finite_mask, 1, 'first');
                    
                    if ~isempty(before_idx) && ~isempty(after_idx)
                        % Draw horizontal line at the top connecting the spikes
                        plot([t_dead(before_idx), t_dead(after_idx)], [1e6, 1e6], ...
                             'Color', params.colors{idx}, 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
    title(sprintf('D_N - Dead (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(5)}(t)');
    grid on;
    % Set y-axis range to 1.2 times the peak (like Susceptible and Infected)
    max_D5 = max(D5_dead(finite_mask));
    if isempty(max_D5) || isnan(max_D5) || isinf(max_D5)
        ylim([0, 1]); % Default range if no valid data
    else
        ylim([0, max_D5 * 1.2]);
    end
    legend(legend_labels, 'Location', 'best');
    saveas(gcf, sprintf('SIHRS_DN5_Dead_Linear_R0_%.2f.png', R0));
    
    % Create combined plot
    figure('Position', [600, 600, 1500, 800]);
    tlayout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot D1 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('D_N - Susceptible (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(1)}(t)');
    grid on;
    ylim([0, max_D1 * 1.2]);
    
    % Plot D2 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        D2_infected = results{idx}.D2;
        D2_original = results{idx}.D2_original;
        
        % Find finite values and plot them
        finite_mask = isfinite(D2_infected);
        if any(finite_mask)
            plot(t(finite_mask), D2_infected(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values in original data and draw vertical lines for blow-ups
        nan_inf_mask = isnan(D2_original) | isinf(D2_original);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t < t(start_idx) & finite_mask, 1, 'last');
                    
                    if ~isempty(before_idx)
                        % Draw vertical line from last finite value to infinity
                        last_finite_value = D2_infected(before_idx);
                        plot([t(before_idx), t(before_idx)], [last_finite_value, 1e9], ...
                             'Color', params.colors{idx}, 'LineWidth', 2.0);
                    end
                    
                    % Draw vertical line at the start of the blow-up region
                    plot([t(start_idx), t(start_idx)], [0, 1e9], ...
                         'Color', params.colors{idx}, 'LineWidth', 2.0);
                end
            end
        end
    end
    title(sprintf('D_N - Infected (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(2)}(t)');
    grid on;
    % Set y-axis range: 1.1 times peak for finite values, but allow blow-ups to go to top
    max_D2_combined = max(cellfun(@(r) max(r.D2(isfinite(r.D2))), results));
    % Use 1.1 times peak value for infected compartment
    ylim([0, max_D2_combined * 1.1]); % Set to 1.1 times peak value
    
    % Plot D3 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for hospitalized
        t_hosp = t(2:end);
        D3_hosp = results{idx}.D3(2:end);
        
        % Handle NaN/Inf values by replacing with large number for horizontal line effect
        D3_plot = D3_hosp;
        D3_plot(isnan(D3_plot) | isinf(D3_plot)) = 1e9;
        
        plot(t_hosp, D3_plot, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('D_N - Hospitalized (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(3)}(t)');
    grid on;
    % Set y-axis range: 1.2 times peak for finite values, but allow blow-ups to go to top
    max_D3_combined = max(D3_plot(D3_plot < 1e9)); % Exclude the 1e9 blow-up values
    if any(D3_plot == 1e9)
        % If there are blow-ups, set y-axis so spikes go all the way to the top
        ylim([0, 1e9]); % Set to exactly 1e9 so spikes reach the top
    else
        % Normal scaling if no blow-ups
        ylim([0, max_D3_combined * 1.2]);
    end

    % Plot D4 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for recovered (will blow up since r_N(0) = 0)
        t_rec = t(2:end);
        D4_rec = results{idx}.D4(2:end);
        
        % Find finite values and plot them
        finite_mask = isfinite(D4_rec);
        if any(finite_mask)
            plot(t_rec(finite_mask), D4_rec(finite_mask), 'Color', params.colors{idx}, 'LineWidth', 1.5);
        end
        
        % Find NaN/Inf values and draw horizontal lines connecting vertical spikes
        nan_inf_mask = isnan(D4_rec) | isinf(D4_rec);
        if any(nan_inf_mask)
            % Find consecutive NaN/Inf regions
            nan_inf_indices = find(nan_inf_mask);
            if ~isempty(nan_inf_indices)
                % Group consecutive NaN/Inf indices
                diff_indices = diff(nan_inf_indices);
                gap_starts = [nan_inf_indices(1); nan_inf_indices(find(diff_indices > 1) + 1)];
                gap_ends = [nan_inf_indices(find(diff_indices > 1)); nan_inf_indices(end)];
                
                for gap_idx = 1:length(gap_starts)
                    start_idx = gap_starts(gap_idx);
                    end_idx = gap_ends(gap_idx);
                    
                    % Find the last finite value before this gap
                    before_idx = find(t_rec < t_rec(start_idx) & finite_mask, 1, 'last');
                    % Find the first finite value after this gap
                    after_idx = find(t_rec > t_rec(end_idx) & finite_mask, 1, 'first');
                    
                    if ~isempty(before_idx) && ~isempty(after_idx)
                        % Draw horizontal line at the top connecting the spikes
                        plot([t_rec(before_idx), t_rec(after_idx)], [1e6, 1e6], ...
                             'Color', params.colors{idx}, 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
    title(sprintf('D_N - Recovered (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(4)}(t)');
    grid on;
    % Set y-axis range to 1.2 times the peak (like Susceptible and Infected)
    max_D4_combined = max(D4_rec(finite_mask));
    if isempty(max_D4_combined) || isnan(max_D4_combined) || isinf(max_D4_combined)
        ylim([0, 1]); % Default range if no valid data
    else
        ylim([0, max_D4_combined * 1.2]);
    end
    
    % Plot D5 in combined figure
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        % Exclude t=0 data for dead
        t_dead = t(2:end);
        D5_dead = results{idx}.D5(2:end);
        
        % Handle NaN/Inf values by replacing with large number for horizontal line effect
        D5_plot = D5_dead;
        D5_plot(isnan(D5_plot) | isinf(D5_plot)) = 1e6;
        
        plot(t_dead, D5_plot, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('D_N - Dead (R_0 = %.2f) - Linear Scale', R0));
    xlabel('Time'); ylabel('D_N^{(5)}(t)');
    grid on;
    % Set y-axis range to 1.2 times the peak (like Susceptible and Infected)
    max_D5_combined = max(D5_plot);
    if isempty(max_D5_combined) || isnan(max_D5_combined) || isinf(max_D5_combined)
        ylim([0, 1]); % Default range if no valid data
    else
        ylim([0, max_D5_combined * 1.2]);
    end
    
    % Add legend to combined figure
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    % Add mathematical formula annotation to combined figure
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', sprintf('D_N^{(l)}(t) = sqrt(V_N^{(l)}(t)) / m_N^{(l)}(t), R_0 = %.2f', R0), ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Save the combined figure
    saveas(gcf, sprintf('SIHRS_DN_renormalized_std_dev_R0_%.2f_combined.png', R0));
    
    % Create population dynamics plot for each N
    figure('Position', [700, 700, 1200, 800]);
    hold on;
    
    % Define distinct colors for each compartment
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
    
    title(sprintf('SIHRS Population Dynamics - All Compartments (R_0 = %.2f)', R0));
    xlabel('Time');
    ylabel('Population Proportion');
    grid on;
    
    % Create custom legend with better organization
    legend('Location', 'eastoutside', 'FontSize', 9);
    
    % Save the population dynamics figure
    saveas(gcf, sprintf('SIHRS_DN_Population_Dynamics_R0_%.2f.png', R0));
    
    fprintf('Generated individual plots:\n');
    fprintf('  - SIHRS_DN1_Susceptible_Linear_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_DN2_Infected_Linear_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_DN3_Hospitalized_Linear_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_DN4_Recovered_Linear_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_DN5_Dead_Linear_R0_%.2f.png\n', R0);
    fprintf('  - SIHRS_DN_renormalized_std_dev_R0_%.2f_combined.png\n', R0);
    fprintf('  - SIHRS_DN_Population_Dynamics_R0_%.2f.png\n', R0);
end

function print_DN_renormalized_summary(results, det_result, R0)
    % Print summary statistics for D_N renormalized standard deviations
    fprintf('\n=== D_N RENORMALIZED STANDARD DEVIATION SUMMARY (R0 = %.2f) ===\n', R0);
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
end

% Run the simulation
sihrs_DN_renormalized_simulation();
