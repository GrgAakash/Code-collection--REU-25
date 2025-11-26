function SIHRS_variance_analysis()
    %% SIHRS Model Square Root of Variance Analysis
    % This function analyzes the stochastic scaling properties of the SIHRS model
    % and validates the sqrt(variance) scaling conjecture: sqrt(V_N^(l)(t)) ∝ sqrt((1/N) × population_proportions)
    
    %% Model Parameters (matching SIHRS.m)
    params = struct(...
        'beta', 0.212, ...           % infection rate (β > 0)
        'gamma', 0.10, ...           % I transition rate (γ > 0)
        'alpha', 0.1, ...            % H transition rate (α > 0)
        'lambda', 0.0083, ...        % R transition rate (Λ > 0) immunity period of 4 months
        'pSI', 1.0, ...              % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0, ...              % probability of I to I (stay infected)
        'pIH', 0.1060, ...           % probability of I to H
        'pIR', 0.8921, ...           % probability of I to R
        'pID', 0.0019, ...           % probability of I to D
        'pHH', 0.00, ...             % probability of H to H (stay hospitalized)
        'pHR', 0.846, ...            % probability of H to R
        'pHD', 0.154, ...            % probability of H to D
        'pRR', 0.02, ...             % probability of R to R (stay recovered)
        'pRS', 0.98, ...             % probability of R to S
        'tmax', 1000, ...            % simulation end time
        's0', 0.96, ...              % initial susceptible proportion
        'i0', 0.04, ...              % initial infected proportion
        'h0', 0.0, ...               % initial hospitalized proportion
        'r0', 0.0, ...               % initial recovered proportion
        'd0', 0.0 ...                % initial dead proportion
    );
    
    % Validate parameters
    validateParameters(params);
    
    % Test different population sizes to see stochastic effects
    N_values = [1600, 3162];
    num_trials = 15; % Number of trials per batch
    num_batches = 5; % Number of batches for each N
    
    % Time grid for analysis
    t_grid = linspace(0, params.tmax, 1000);
    
    % Pre-allocate storage for variance analysis (5 batches x 15 trials each)
    S_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
    I_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
    H_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
    R_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
    D_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
    
    % Pre-allocate storage for batch-averaged expected variance
    S_expected_var = zeros(length(N_values), length(t_grid), num_batches);
    I_expected_var = zeros(length(N_values), length(t_grid), num_batches);
    H_expected_var = zeros(length(N_values), length(t_grid), num_batches);
    R_expected_var = zeros(length(N_values), length(t_grid), num_batches);
    D_expected_var = zeros(length(N_values), length(t_grid), num_batches);
    
    % Cell arrays to store theoretical variance data
    v_s_theory = cell(length(N_values), 1);
    v_i_theory = cell(length(N_values), 1);
    v_h_theory = cell(length(N_values), 1);
    v_r_theory = cell(length(N_values), 1);
    v_d_theory = cell(length(N_values), 1);
    
    %% Main Simulation Loop
    for idx = 1:length(N_values)
        N = N_values(idx);
        fprintf('Running variance analysis for N = %d...\n', N);
        
        % Run 5 batches of 15 trials each
        for batch = 1:num_batches
            fprintf('  Batch %d/%d:\n', batch, num_batches);
            
            % Run 15 trials for this batch
            for trial = 1:num_trials
                fprintf('    Trial %d/%d\n', trial, num_trials);
                
                % Run single stochastic simulation
                result = sihrs_agent_model(N, params);
                
                % Interpolate sqrt(variance) data onto common time grid
                S_var_all(idx, :, trial, batch) = sqrt(interp1(result.T, result.vs, t_grid, 'linear', 0));
                I_var_all(idx, :, trial, batch) = sqrt(interp1(result.T, result.vi, t_grid, 'linear', 0));
                H_var_all(idx, :, trial, batch) = sqrt(interp1(result.T, result.vh, t_grid, 'linear', 0));
                R_var_all(idx, :, trial, batch) = sqrt(interp1(result.T, result.vr, t_grid, 'linear', 0));
                D_var_all(idx, :, trial, batch) = sqrt(interp1(result.T, result.vd, t_grid, 'linear', 0));
            end
            
            % Calculate expected variance for this batch (average across 15 trials)
            S_expected_var(idx, :, batch) = mean(S_var_all(idx, :, :, batch), 3);
            I_expected_var(idx, :, batch) = mean(I_var_all(idx, :, :, batch), 3);
            H_expected_var(idx, :, batch) = mean(H_var_all(idx, :, :, batch), 3);
            R_expected_var(idx, :, batch) = mean(R_var_all(idx, :, :, batch), 3);
            D_expected_var(idx, :, batch) = mean(D_var_all(idx, :, :, batch), 3);
            
            fprintf('    Batch %d expected variance calculated\n', batch);
        end
        
        %% Theoretical Variance Calculation using ODE solution
        % Solve deterministic ODE for theoretical comparison
        det_result = solve_deterministic_sihrs(params);
        
        % Interpolate ODE solution to our time grid
        s_ode = interp1(det_result.T, det_result.S_prop, t_grid, 'linear', 0);
        i_ode = interp1(det_result.T, det_result.I_prop, t_grid, 'linear', 0);
        h_ode = interp1(det_result.T, det_result.H_prop, t_grid, 'linear', 0);
        r_ode = interp1(det_result.T, det_result.R_prop, t_grid, 'linear', 0);
        d_ode = interp1(det_result.T, det_result.D_prop, t_grid, 'linear', 0);
        
        % Calculate theoretical sqrt(variance) using the mathematical formulas
        % sqrt(V_N^(1)(t)) = sqrt((1/N)(β s(t) i(t) p_SI + Λ r(t) p_RS))
        v_s_theory{idx} = sqrt((1/N) * (params.beta * s_ode .* i_ode * params.pSI + params.lambda * r_ode * params.pRS));
        
        % sqrt(V_N^(2)(t)) = sqrt((1/N)(γ i(t) (p_IH + p_IR + p_ID) + β s(t) i(t) p_SI))
        v_i_theory{idx} = sqrt((1/N) * (params.gamma * i_ode * (params.pIH + params.pIR + params.pID) + params.beta * s_ode .* i_ode * params.pSI));
        
        % sqrt(V_N^(3)(t)) = sqrt((1/N)(γ i(t) p_IH + α h(t) (p_HR + p_HD)))
        v_h_theory{idx} = sqrt((1/N) * (params.gamma * i_ode * params.pIH + params.alpha * h_ode * (params.pHR + params.pHD)));
        
        % sqrt(V_N^(4)(t)) = sqrt((1/N)(γ i(t) p_IR + α h(t) p_HR + Λ r(t) p_RS))
        v_r_theory{idx} = sqrt((1/N) * (params.gamma * i_ode * params.pIR + params.alpha * h_ode * params.pHR + params.lambda * r_ode * params.pRS));
        
        % sqrt(V_N^(5)(t)) = sqrt((1/N)(γ i(t) p_ID + α h(t) p_HD))
        v_d_theory{idx} = sqrt((1/N) * (params.gamma * i_ode * params.pID + params.alpha * h_ode * params.pHD));
    end
    
    %% Time-Averaged Square Root of Variance Analysis
    % Define time intervals for analysis [T1, T2]
    time_intervals = {[100, 200], [300, 400], [500, 600], [700, 800],[800, 900],[900, 1000]};
    
    fprintf('\n=== TIME-AVERAGED SQRT(VARIANCE) ANALYSIS ===\n');
    for interval_idx = 1:length(time_intervals)
        T1 = time_intervals{interval_idx}(1);
        T2 = time_intervals{interval_idx}(2);
        
        fprintf('\nTime interval [%.0f, %.0f]:\n', T1, T2);
        
        % Find indices for this time interval
        t_indices = t_grid >= T1 & t_grid <= T2;
        
        for n_idx = 1:length(N_values)
            N = N_values(n_idx);
            fprintf('  N = %d:\n', N);
            
            % Calculate time-averaged simulated sqrt(variance)
            V_bar_S = mean(mean(S_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_I = mean(mean(I_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_H = mean(mean(H_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_R = mean(mean(R_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_D = mean(mean(D_var_all(n_idx, t_indices, :), 3), 2);
            
            % Calculate time-averaged theoretical sqrt(variance)
            v_bar_S = mean(v_s_theory{n_idx}(t_indices));
            v_bar_I = mean(v_i_theory{n_idx}(t_indices));
            v_bar_H = mean(v_h_theory{n_idx}(t_indices));
            v_bar_R = mean(v_r_theory{n_idx}(t_indices));
            v_bar_D = mean(v_d_theory{n_idx}(t_indices));
            
            % Calculate time-averaged population proportions
            s_avg = mean(s_ode(t_indices));
            i_avg = mean(i_ode(t_indices));
            h_avg = mean(h_ode(t_indices));
            r_avg = mean(r_ode(t_indices));
            d_avg = mean(d_ode(t_indices));
            
            fprintf('    S: V̄_sim=%.6f, V̄_theory=%.6f, s̄=%.4f\n', V_bar_S, v_bar_S, s_avg);
            fprintf('    I: V̄_sim=%.6f, V̄_theory=%.6f, ī=%.4f\n', V_bar_I, v_bar_I, i_avg);
            fprintf('    H: V̄_sim=%.6f, V̄_theory=%.6f, h̄=%.4f\n', V_bar_H, v_bar_H, h_avg);
            fprintf('    R: V̄_sim=%.6f, V̄_theory=%.6f, r̄=%.4f\n', V_bar_R, v_bar_R, r_avg);
            fprintf('    D: V̄_sim=%.6f, V̄_theory=%.6f, d̄=%.4f\n', V_bar_D, v_bar_D, d_avg);
        end
    end
    
    %% Plotting Expected Square Root of Variance Analysis
    plot_variance_error_bars(N_values, t_grid, S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var, ...
        v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory);
    
    fprintf('\nExpected sqrt(variance) analysis completed successfully!\n');
    fprintf('Following the new conjecture: E[sqrt(V_N^(l)(t))] proportional to sqrt((1/N) x population_proportions)\n');
end

%% Helper Functions
function validateParameters(params)
    % Check that all rates are positive
    if any([params.beta, params.gamma, params.alpha, params.lambda] <= 0)
        error('All rates (beta, gamma, alpha, lambda) must be positive');
    end
    
    % Probabilities should be between 0 and 1
    probs = [params.pSI, params.pII, params.pIH, params.pIR, params.pID, ...
             params.pHH, params.pHR, params.pHD, params.pRR, params.pRS];
    if any(probs < 0 | probs > 1)
        error('All probabilities must be in [0,1]');
    end
    
    % Make sure probability sums add up correctly
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

function result = sihrs_agent_model(N, params)
    % SIHRS agent-based stochastic model with variance tracking
    % This is a modified version that tracks variance at each step
    
    % Initial conditions
    s0 = round(params.s0 * N);
    i0 = round(params.i0 * N);
    h0 = round(params.h0 * N);
    r0 = round(params.r0 * N);
    d0 = round(params.d0 * N);
    
    % Initialize arrays
    max_events = N * 30;
    T = zeros(1, max_events);
    S_prop = zeros(1, max_events);
    I_prop = zeros(1, max_events);
    H_prop = zeros(1, max_events);
    R_prop = zeros(1, max_events);
    D_prop = zeros(1, max_events);
    
    % Variance tracking arrays
    vs = zeros(1, max_events);
    vi = zeros(1, max_events);
    vh = zeros(1, max_events);
    vr = zeros(1, max_events);
    vd = zeros(1, max_events);
    
    % Initialize
    S = 1:s0;
    I = (s0+1):(s0+i0);
    H = [];
    R = [];
    D = [];
    
    t = 0;
    event_count = 1;
    
    % Record initial state
    T(1) = 0;
    S_prop(1) = s0 / N;
    I_prop(1) = i0 / N;
    H_prop(1) = 0;
    R_prop(1) = 0;
    D_prop(1) = 0;
    
            % Calculate initial variance (note: we still calculate variance here, sqrt is applied later)
        vs(1) = (params.beta * (s0/N) * (i0/N) * params.pSI + params.lambda * (r0/N) * params.pRS) / N;
        vi(1) = (params.gamma * (i0/N) * (params.pIH + params.pIR + params.pID) + params.beta * (s0/N) * (i0/N) * params.pSI) / N;
        vh(1) = (params.gamma * (i0/N) * params.pIH + params.alpha * 0 * (params.pHR + params.pHD)) / N;
        vr(1) = (params.gamma * (i0/N) * params.pIR + params.alpha * 0 * params.pHR + params.lambda * (r0/N) * params.pRS) / N;
        vd(1) = (params.gamma * (i0/N) * params.pID + params.alpha * 0 * params.pHD) / N;
    
    % Main simulation loop
    while t < params.tmax && event_count < max_events
        ns = numel(S);
        ni = numel(I);
        nh = numel(H);
        nr = numel(R);
        nd = numel(D);
        
        % Calculate event rates
        infection_rate = params.beta * ns * ni / N * params.pSI;
        i_to_h_rate = params.gamma * ni * params.pIH;
        i_to_r_rate = params.gamma * ni * params.pIR;
        i_to_d_rate = params.gamma * ni * params.pID;
        h_to_r_rate = params.alpha * nh * params.pHR;
        h_to_d_rate = params.alpha * nh * params.pHD;
        r_to_s_rate = params.lambda * nr * params.pRS;
        
        total_rate = infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate + h_to_d_rate + r_to_s_rate;
        
        if total_rate == 0
            break;
        end
        
        % Time to next event
        dt = exprnd(1 / total_rate);
        t = t + dt;
        
        if t > params.tmax
            t = params.tmax;
        end
        
        event_count = event_count + 1;
        T(event_count) = t;
        
        % Determine which event occurs
        r = rand * total_rate;
        if r < infection_rate
            % S -> I
            if ns > 0
                I(end+1) = S(end);
                S(end) = [];
            end
        elseif r < infection_rate + i_to_h_rate
            % I -> H
            if ni > 0
                H(end+1) = I(end);
                I(end) = [];
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate
            % I -> R
            if ni > 0
                R(end+1) = I(end);
                I(end) = [];
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate
            % I -> D
            if ni > 0
                D(end+1) = I(end);
                I(end) = [];
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate
            % H -> R
            if nh > 0
                R(end+1) = H(end);
                H(end) = [];
            end
        elseif r < infection_rate + i_to_h_rate + i_to_r_rate + i_to_d_rate + h_to_r_rate + h_to_d_rate
            % H -> D
            if nh > 0
                D(end+1) = H(end);
                H(end) = [];
            end
        else
            % R -> S
            if nr > 0
                S(end+1) = R(end);
                R(end) = [];
            end
        end
        
        % Update proportions
        ns = numel(S);
        ni = numel(I);
        nh = numel(H);
        nr = numel(R);
        nd = numel(D);
        
        S_prop(event_count) = ns / N;
        I_prop(event_count) = ni / N;
        H_prop(event_count) = nh / N;
        R_prop(event_count) = nr / N;
        D_prop(event_count) = nd / N;
        
        % Calculate variance at this time step (note: we still calculate variance here, sqrt is applied later)
        vs(event_count) = (params.beta * (ns/N) * (ni/N) * params.pSI + params.lambda * (nr/N) * params.pRS) / N;
        vi(event_count) = (params.gamma * (ni/N) * (params.pIH + params.pIR + params.pID) + params.beta * (ns/N) * (ni/N) * params.pSI) / N;
        vh(event_count) = (params.gamma * (ni/N) * params.pIH + params.alpha * (nh/N) * (params.pHR + params.pHD)) / N;
        vr(event_count) = (params.gamma * (ni/N) * params.pIR + params.alpha * (nh/N) * params.pHR + params.lambda * (nr/N) * params.pRS) / N;
        vd(event_count) = (params.gamma * (ni/N) * params.pID + params.alpha * (nh/N) * params.pHD) / N;
    end
    
    % Trim arrays
    T = T(1:event_count);
    S_prop = S_prop(1:event_count);
    I_prop = I_prop(1:event_count);
    H_prop = H_prop(1:event_count);
    R_prop = R_prop(1:event_count);
    D_prop = D_prop(1:event_count);
    vs = vs(1:event_count);
    vi = vi(1:event_count);
    vh = vh(1:event_count);
    vr = vr(1:event_count);
    vd = vd(1:event_count);
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.H_prop = H_prop;
    result.R_prop = R_prop;
    result.D_prop = D_prop;
    result.vs = vs;
    result.vi = vi;
    result.vh = vh;
    result.vr = vr;
    result.vd = vd;
    result.final_time = t;
end

function det_result = solve_deterministic_sihrs(params)
    % Solve the deterministic SIHRS model using ODE45
    tspan = [0, params.tmax];
    y0 = [params.s0; params.i0; params.h0; params.r0; params.d0];
    
    % Define the ODE system
    ode_system = @(t, y) [
        -params.beta * y(1) * y(2) * params.pSI + params.pRS * params.lambda * y(4); % ds/dt
        params.beta * y(1) * y(2) * params.pSI - params.gamma * (1 - params.pII) * y(2); % di/dt
        params.pIH * params.gamma * y(2) - params.alpha * (1 - params.pHH) * y(3); % dh/dt
        params.pIR * params.gamma * y(2) + params.pHR * params.alpha * y(3) - params.pRS * params.lambda * y(4); % dr/dt
        params.pID * params.gamma * y(2) + params.pHD * params.alpha * y(3) % dd/dt
    ];
    
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [T, Y] = ode45(ode_system, tspan, y0, options);
    
    det_result.T = T;
    det_result.S_prop = Y(:, 1);
    det_result.I_prop = Y(:, 2);
    det_result.H_prop = Y(:, 3);
    det_result.R_prop = Y(:, 4);
    det_result.D_prop = Y(:, 5);
end



function plot_variance_error_bars(N_values, t_grid, S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var, ...
                                v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory)
    % Create error bar plots showing final expected sqrt(variance) (average of 5 batches) vs theoretical
    % Following the new conjecture: E[sqrt(V_N^(l)(t))] proportional to sqrt((1/N) x population_proportions)
    dt = 2; % Half width for time intervals around each midpoint (adjusted for 700-day simulation)
    all_midpoints = 5:dt*2:1000; % All possible time points for analysis
    error_bar_interval = 4; % Show error bar every nth window
    midpoints = all_midpoints(1:error_bar_interval:end); % Reduced set of midpoints for error bars
    compartments = {'Susceptible', 'Infected', 'Hospitalized', 'Recovered', 'Dead'};
    all_sim_vars = {S_expected_var, I_expected_var, H_expected_var, R_expected_var, D_expected_var};
    all_theory_vars = {v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory};
    colors = lines(length(N_values));
    
    % Create separate figures for each compartment
    for comp_idx = 1:numel(compartments)
        figure('Position', [100 + (comp_idx-1)*50, 100 + (comp_idx-1)*50, 1600, 800]);
        
        for n_idx = 1:numel(N_values)
            subplot(1, length(N_values), n_idx);
            hold on;
            grid on;
            N = N_values(n_idx);
            
            % Get batch data for this compartment and N
            batch_data = squeeze(all_sim_vars{comp_idx}(n_idx, :, :)); % [time_points x 5_batches]
            theory_data = all_theory_vars{comp_idx}{n_idx};

            % Plot theoretical sqrt(variance) curve (smooth line)
            plot(t_grid, theory_data, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Theoretical Curve');
            
            for t_idx = 1:length(midpoints)
                t_mid = midpoints(t_idx);
                t_min = t_mid - dt;
                t_max = t_mid + dt; 

                in_window_indices = t_grid >= t_min & t_grid < t_max;
                if ~any(in_window_indices)
                    continue;
                end
             
                % Get batch values in this time window
                batch_vals_in_window = batch_data(in_window_indices, :);
                
                % Calculate min and max across batches in this time window
                batch_mins = min(batch_vals_in_window, [], 1);  % Min of each batch
                batch_maxs = max(batch_vals_in_window, [], 1);  % Max of each batch

                % Calculate final expected sqrt(variance) bounds (average across 5 batches)
                y_min_avg = mean(batch_mins(~isnan(batch_mins)));  % Final E[min sqrt(V_N^(l)(t))]
                y_max_avg = mean(batch_maxs(~isnan(batch_maxs)));  % Final E[max sqrt(V_N^(l)(t))]

                if ~isnan(y_min_avg) && ~isnan(y_max_avg)
                    % Plot error bar showing final expected sqrt(variance) range
                    plot([t_mid, t_mid], [y_min_avg, y_max_avg], 'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', 'Expected Sqrt(Variance) Range');
                end
            end
        
            xlabel('Time');
            ylabel('Sqrt(Variance)');
            title([compartments{comp_idx}, ' - N = ', num2str(N), ' (5 batches x 15 runs)']);
            legend('Theoretical Curve', 'Expected Sqrt(Variance) Range', 'Location', 'best');
            
            % Add inset zoom plot for interesting region
            add_inset_zoom(gca, t_grid, theory_data, midpoints, batch_data, colors(n_idx,:), dt, N);
            
            hold off;
        end
        
        % Add title for each compartment figure
        sgtitle([compartments{comp_idx}, ' - Sqrt(Variance) Analysis'], 'FontSize', 16);
        
        % Save each compartment figure separately
        saveas(gcf, ['SIHRS_', compartments{comp_idx}, '_final_expected_sqrtvariance.png']);
    end
    
    fprintf('Created separate figures for each compartment:\n');
    for comp_idx = 1:numel(compartments)
        fprintf('  - %s: SIHRS_%s_final_expected_sqrtvariance.png\n', compartments{comp_idx}, compartments{comp_idx});
    end
end

function add_inset_zoom(main_ax, t_grid, theory_data, midpoints, batch_data, error_color, dt, N)
    % Add an inset zoom plot to show detailed view of an interesting region
    
    % Define zoom region (adjust these values based on your data)
    zoom_start = 100;  % Start time for zoom region
    zoom_end = 300;    % End time for zoom region
    
    % Find indices for zoom region
    zoom_indices = t_grid >= zoom_start & t_grid <= zoom_end;
    if sum(zoom_indices) == 0
        return; % Skip if no data in zoom region
    end
    
    % Get current figure and main axes position
    main_pos = get(main_ax, 'Position');
    
    % Create inset axes (positioned in bottom-right corner of main plot)
    inset_width = 0.25;  % Width as fraction of main plot
    inset_height = 0.25; % Height as fraction of main plot
    inset_x = main_pos(1) + main_pos(3) - inset_width - 0.05;  % Right side with margin
    inset_y = main_pos(2) + 0.05;  % Bottom with margin
    
    inset_ax = axes('Position', [inset_x, inset_y, inset_width, inset_height]);
    
    % Plot theoretical curve in zoom region
    plot(inset_ax, t_grid(zoom_indices), theory_data(zoom_indices), '--', ...
         'Color', 'k', 'LineWidth', 1.5);
    hold(inset_ax, 'on');
    
    % Add error bars in zoom region (only midpoints within zoom range)
    zoom_midpoints = midpoints(midpoints >= zoom_start & midpoints <= zoom_end);
    
    for t_idx = 1:length(zoom_midpoints)
        t_mid = zoom_midpoints(t_idx);
        t_min = t_mid - dt;
        t_max = t_mid + dt;
        
        % Find data points in this time window
        in_window_indices = t_grid >= t_min & t_grid < t_max;
        if ~any(in_window_indices)
            continue;
        end
        
        % Get batch values in this time window
        batch_vals_in_window = batch_data(in_window_indices, :);
        
        % Calculate min and max across batches
        batch_mins = min(batch_vals_in_window, [], 1);
        batch_maxs = max(batch_vals_in_window, [], 1);
        
        % Calculate final expected bounds
        y_min_avg = mean(batch_mins(~isnan(batch_mins)));
        y_max_avg = mean(batch_maxs(~isnan(batch_maxs)));
        
        if ~isnan(y_min_avg) && ~isnan(y_max_avg)
            % Plot error bar
            plot(inset_ax, [t_mid, t_mid], [y_min_avg, y_max_avg], ...
                 'Color', error_color, 'LineWidth', 1.2);
        end
    end
    
    % Customize inset plot
    set(inset_ax, 'FontSize', 8);
    xlabel(inset_ax, 'Time', 'FontSize', 8);
    ylabel(inset_ax, 'Sqrt(Var)', 'FontSize', 8);
    title(inset_ax, sprintf('Zoom: t=[%d,%d]', zoom_start, zoom_end), 'FontSize', 9);
    grid(inset_ax, 'on');
    
    % Set limits for zoom region
    xlim(inset_ax, [zoom_start, zoom_end]);
    
    % Add a box around the inset to make it stand out
    set(inset_ax, 'Box', 'on', 'LineWidth', 1.2);
    
    hold(inset_ax, 'off');
    
    % Optional: Add connection lines from main plot to inset
    % This creates the "magnifying glass" effect
    main_xlim = get(main_ax, 'XLim');
    main_ylim = get(main_ax, 'YLim');
    
    % Draw connection lines (commented out as they can be cluttering)
    % You can uncomment these if you want the connection lines
    % annotation('line', [zoom_start/main_xlim(2), inset_x], ...
    %           [0.1, inset_y], 'Color', 'red', 'LineStyle', '--');
    % annotation('line', [zoom_end/main_xlim(2), inset_x+inset_width], ...
    %           [0.1, inset_y], 'Color', 'red', 'LineStyle', '--');
end
