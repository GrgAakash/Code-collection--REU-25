function SIHRS_variance_analysis()
    %% SIHRS Model Variance Analysis
    % This function analyzes the stochastic scaling properties of the SIHRS model
    % and validates the variance scaling conjecture: V_N^(l)(t) ∝ (1/N) × population_proportions
    
    %% Model Parameters (matching SIHRS.m)
    params = struct(...
        'beta', 0.212, ...           % infection rate (β > 0)
        'gamma', 0.10, ...           % I transition rate (γ > 0)
        'alpha', 0.1, ...            % H transition rate (α > 0)
        'lambda', 0.0083, ...        % R transition rate (Λ > 0) immunity period of 4 months
        'pSI', 1.0, ...              % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0, ...              % probability of I to I (stay infected)
        'pIH', 0.04, ...             % probability of I to H
        'pIR', 0.959, ...            % probability of I to R
        'pID', 0.001, ...            % probability of I to D
        'pHH', 0.01, ...             % probability of H to H (stay hospitalized)
        'pHR', 0.9882, ...           % probability of H to R
        'pHD', 0.0018, ...           % probability of H to D
        'pRR', 0.02, ...             % probability of R to R (stay recovered)
        'pRS', 0.98, ...             % probability of R to S
        'tmax', 100, ...            % simulation end time
        's0', 0.96, ...              % initial susceptible proportion
        'i0', 0.04, ...              % initial infected proportion
        'h0', 0.0, ...               % initial hospitalized proportion
        'r0', 0.0, ...               % initial recovered proportion
        'd0', 0.0 ...                % initial dead proportion
    );
    
    % Validate parameters
    validateParameters(params);
    
    % Test different population sizes to see stochastic effects
    N_values = [316, 3162, 10000];
    num_trials = 15; % Number of trials for each N (increased for better statistical power)
    
    % Time grid for analysis
    t_grid = linspace(0, params.tmax, 1000);
    
    % Pre-allocate storage for variance analysis
    S_var_all = zeros(length(N_values), length(t_grid), num_trials);
    I_var_all = zeros(length(N_values), length(t_grid), num_trials);
    H_var_all = zeros(length(N_values), length(t_grid), num_trials);
    R_var_all = zeros(length(N_values), length(t_grid), num_trials);
    D_var_all = zeros(length(N_values), length(t_grid), num_trials);
    
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
        
        % Run multiple stochastic simulations
        for trial = 1:num_trials
            fprintf('  Trial %d/%d\n', trial, num_trials);
            
            % Run single stochastic simulation
            result = sihrs_agent_model(N, params);
            
            % Interpolate variance data onto common time grid
            S_var_all(idx, :, trial) = interp1(result.T, result.vs, t_grid, 'linear', 0);
            I_var_all(idx, :, trial) = interp1(result.T, result.vi, t_grid, 'linear', 0);
            H_var_all(idx, :, trial) = interp1(result.T, result.vh, t_grid, 'linear', 0);
            R_var_all(idx, :, trial) = interp1(result.T, result.vr, t_grid, 'linear', 0);
            D_var_all(idx, :, trial) = interp1(result.T, result.vd, t_grid, 'linear', 0);
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
        
        % Calculate theoretical variance using the mathematical formulas
        % V_N^(1)(t) = (1/N)(β s(t) i(t) p_SI + Λ r(t) p_RS)
        v_s_theory{idx} = (1/N) * (params.beta * s_ode .* i_ode * params.pSI + params.lambda * r_ode * params.pRS);
        
        % V_N^(2)(t) = (1/N)(γ i(t) (p_IH + p_IR + p_ID) + β s(t) i(t) p_SI)
        v_i_theory{idx} = (1/N) * (params.gamma * i_ode * (params.pIH + params.pIR + params.pID) + params.beta * s_ode .* i_ode * params.pSI);
        
        % V_N^(3)(t) = (1/N)(γ i(t) p_IH + α h(t) (p_HR + p_HD))
        v_h_theory{idx} = (1/N) * (params.gamma * i_ode * params.pIH + params.alpha * h_ode * (params.pHR + params.pHD));
        
        % V_N^(4)(t) = (1/N)(γ i(t) p_IR + α h(t) p_HR + Λ r(t) p_RS)
        v_r_theory{idx} = (1/N) * (params.gamma * i_ode * params.pIR + params.alpha * h_ode * params.pHR + params.lambda * r_ode * params.pRS);
        
        % V_N^(5)(t) = (1/N)(γ i(t) p_ID + α h(t) p_HD)
        v_d_theory{idx} = (1/N) * (params.gamma * i_ode * params.pID + params.alpha * h_ode * params.pHD);
    end
    
    %% Time-Averaged Variance Analysis
    % Define time intervals for analysis [T1, T2]
    time_intervals = {[100, 200], [300, 400], [500, 600], [700, 800]};
    
    fprintf('\n=== TIME-AVERAGED VARIANCE ANALYSIS ===\n');
    for interval_idx = 1:length(time_intervals)
        T1 = time_intervals{interval_idx}(1);
        T2 = time_intervals{interval_idx}(2);
        
        fprintf('\nTime interval [%.0f, %.0f]:\n', T1, T2);
        
        % Find indices for this time interval
        t_indices = t_grid >= T1 & t_grid <= T2;
        
        for n_idx = 1:length(N_values)
            N = N_values(n_idx);
            fprintf('  N = %d:\n', N);
            
            % Calculate time-averaged simulated variance
            V_bar_S = mean(mean(S_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_I = mean(mean(I_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_H = mean(mean(H_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_R = mean(mean(R_var_all(n_idx, t_indices, :), 3), 2);
            V_bar_D = mean(mean(D_var_all(n_idx, t_indices, :), 3), 2);
            
            % Calculate time-averaged theoretical variance
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
    
    %% Plotting Expected Variance Analysis
    plot_variance_error_bars(N_values, t_grid, S_var_all, I_var_all, H_var_all, R_var_all, D_var_all, ...
        v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory);
    
    fprintf('\nExpected variance analysis completed successfully!\n');
    fprintf('Following the new conjecture: E[V_N^(l)(t)] proportional to (1/N) x population_proportions\n');
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
    max_events = N * 10;
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
    
    % Calculate initial variance
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
        
        % Calculate variance at this time step
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



function plot_variance_error_bars(N_values, t_grid, S_var_all, I_var_all, H_var_all, R_var_all, D_var_all, ...
                                v_s_theory, v_i_theory, v_h_theory, v_r_theory, v_d_theory)
    % Create error bar plots showing expected variance ranges with theoretical reference
    % Following the new conjecture: E[V_N^(l)(t)] proportional to (1/N) x population_proportions
    dt = 2; % Half width for time intervals around each midpoint
    midpoints = 5:dt*2:95; % Time points for analysis (covering most of 100-day simulation)
    compartments = {'Susceptible', 'Infected', 'Hospitalized', 'Recovered', 'Dead'};
    all_sim_vars = {S_var_all, I_var_all, H_var_all, R_var_all, D_var_all};
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
            
            all_trials_data = squeeze(all_sim_vars{comp_idx}(n_idx, :, :))';
            theory_data = all_theory_vars{comp_idx}{n_idx};

            for t_idx = 1:length(midpoints)
                t_mid = midpoints(t_idx);
                t_min = t_mid - dt;
                t_max = t_mid + dt; 

                in_window_indices = t_grid >= t_min & t_grid < t_max;
                if ~any(in_window_indices)
                    continue;
                end
             
                vals_in_window = all_trials_data(:, in_window_indices);
                
                % Calculate expected variance following the new conjecture:
                % For each trial, find min and max in the time window
                trial_mins = min(vals_in_window, [], 2);  % Min of each trial
                trial_maxs = max(vals_in_window, [], 2);  % Max of each trial

                % Average across trials to get expected variance bounds
                y_min_avg = mean(trial_mins(~isnan(trial_mins)));  % E[min V_N^(l)(t)]
                y_max_avg = mean(trial_maxs(~isnan(trial_maxs)));  % E[max V_N^(l)(t)]

                if ~isnan(y_min_avg) && ~isnan(y_max_avg)
                    % Plot error bar showing expected variance range
                    plot([t_mid, t_mid], [y_min_avg, y_max_avg], 'Color', colors(n_idx,:), 'LineWidth', 1.5);
                    
                    % Plot theoretical expected variance value
                    [~, time_index] = min(abs(t_grid - t_mid));
                    y_theory = theory_data(time_index); 
                    plot(t_mid, y_theory, '_', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
                end
            end
        
            xlabel('Time');
            ylabel('Expected Variance E[V_N^{(l)}(t)]');
            title([compartments{comp_idx}, ' - N = ', num2str(N)]);
            legend('Simulated Expected Range', 'Theoretical Expected Value', 'Location', 'best');
            hold off;
        end
        
        % Add title for each compartment figure
        sgtitle([compartments{comp_idx}, ' - Expected Variance Analysis'], 'FontSize', 16);
        
        % Save each compartment figure separately
        saveas(gcf, ['SIHRS_', compartments{comp_idx}, '_expected_variance.png']);
    end
    
    fprintf('Created separate figures for each compartment:\n');
    for comp_idx = 1:numel(compartments)
        fprintf('  - %s: SIHRS_%s_expected_variance.png\n', compartments{comp_idx}, compartments{comp_idx});
    end
end
