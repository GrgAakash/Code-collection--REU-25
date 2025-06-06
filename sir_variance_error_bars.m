% SIR Variance Error Bars
% Caden McCann / Amy Carr
% This function generates a plot for each compartment + N, and compares the
% theoretical variance from ODE system with the simulated variance.

function sir_variance_error_bars()
    %% Model Parameters
    beta = 0.95;
    gamma = 1;
    p1 = 0.5;
    p2 = 0.5;
    max_T = 23;
    s_prop = 0.96;
    i_prop = 0.04;
    r_prop = 1 - s_prop - i_prop;
    Ns = [316, 1000, 3162, 10000];
    num_trials = 20;
    
    %% Agent-Based Simulations
    t_grid = linspace(0, max_T, 1000);
    S_var_all = zeros(length(Ns), length(t_grid), num_trials);
    I_var_all = zeros(length(Ns), length(t_grid), num_trials);
    R_var_all = zeros(length(Ns), length(t_grid), num_trials);
    v_s_dot = cell(length(Ns), 1);
    v_i_dot = cell(length(Ns), 1);
    v_r_dot = cell(length(Ns), 1);
    
    %% Repeat simulation for each N
    for idx = 1:length(Ns)
        for trial = 1:num_trials
        % initializations for simulation
        N = Ns(idx);
        s0_abs = round(s_prop * N);
        i0_abs = round(i_prop * N);
        r0_abs = round(r_prop * N);
        S = 1:s0_abs;
        I = (s0_abs+1):(s0_abs+i0_abs);
        R = [];
        t = 0;
        T = 0;
        S_prop = s0_abs / N;
        I_prop = i0_abs / N;
        R_prop = r0_abs / N;
        S_var = 0;
        I_var = 0;
        R_var = 0;
        I_count = numel(I);  % initial infected count
        T = 0;
        I_count_list = I_count;  % store history of I counts

        %% Begin Simulation
        while ~isempty(I) && t < max_T
            nI = numel(I);
            if nI < 1
                nI = 0;
            end    
            nS = numel(S);

            % Clock configuration
            infection_rate = beta * nS * nI / (N^1.0);
            recovery_rate = gamma * nI;
            event_rate = infection_rate + recovery_rate;
            dt = exprnd(1 / event_rate);
            t = t + dt;

            % Instantaneous variance formulas
            S_var(end+1) = p1 * beta * (nI/N) * (nS/N) / N;
            I_var(end+1) = (p1 * beta * (nI/N) * (nS/N) + p2 * gamma * (nI/N)) / N;
            R_var(end+1) = p2 * gamma * (nI/N) / N;

            % Split clock
            if rand < (infection_rate / event_rate)
                if nS > 0 && rand < p1
                    num = randi(nS);
                    infected = S(num);
                    S(num) = [];
                    I(end+1) = infected;
                end
            else
                if rand < p2
                    num = randi(nI);
                    recovered = I(num);
                    I(num) = [];
                    R(end+1) = recovered;
                end
            end

            % Update counts
            N_current = numel(S) + numel(I) + numel(R);
            T(end+1) = t;
            S_prop(end+1) = numel(S) / N_current;
            I_prop(end+1) = numel(I) / N_current;
            R_prop(end+1) = numel(R) / N_current;
            I_count_list(end+1) = numel(I);
            
        end

        % Fix first value since it was initially 0
        S_var(1) = S_var(2);
        I_var(1) = I_var(2);
        R_var(1) = R_var(2);
        
        % Append final point at (23,0) if simulation ended early
        if T(end) < 23
            T(end+1) = 23;
            S_var(end+1) = 0;
            I_var(end+1) = 0;
            R_var(end+1) = 0;
        end

        % Interpolate points into lines, and store them in large arrays
        S_interp_var = interp1(T, S_var, t_grid, 'linear', 0);
        I_interp_var = interp1(T, I_var, t_grid, 'linear', 0);
        R_interp_var = interp1(T, R_var, t_grid, 'linear', 0);
        S_var_all(idx, :, trial) = S_interp_var;
        I_var_all(idx, :, trial) = I_interp_var;
        R_var_all(idx, :, trial) = R_interp_var;
        end

        %% Deterministic ODE Model    
        % Solve ODE system
        [~, y] = ode45(@(t, y) sir_system(t, y, beta, gamma, p1, p2), t_grid, [s_prop, i_prop, r_prop]);
        s = y(:, 1);
        i = y(:, 2);
    
        % Compute instantaneous variance for ODE system
        v_s_dot{idx} = p1 * beta * s .* i / N;
        v_i_dot{idx} = (p1 * beta * s .* i + p2 * gamma * i) / N;
        v_r_dot{idx} = p2 * gamma * i / N;
    end

    %% Graphing
    plot_graphs(Ns, t_grid, S_var_all, I_var_all, R_var_all, v_s_dot, v_i_dot, v_r_dot);
end

%% Code for ODE system
function dydt = sir_system(~, y, beta, gamma, p1, p2)
    s = y(1);
    i = y(2);

    ds = -p1 * beta * s * i;
    di = p1 * beta * s * i - p2 * gamma * i;
    dr = p2 * gamma * i;

    dydt = [ds; di; dr];
end

%% Code for plotting graphs
function plot_graphs(Ns, t_grid, S_var_all, I_var_all, R_var_all, v_s_dot, v_i_dot, v_r_dot)
    dt = 1.4; % Half width for interval around each midpoint
    midpoints = 0.5:dt*2:22.5;
    compartments = {'Susceptible', 'Infected', 'Recovered'};
    colors = lines(length(Ns));
    
    for comp = 1:3 % change to comp = 2:2 for just infected etc
        for idx = 1:numel(Ns)
            figure;
            hold on;
            N = Ns(idx);

            % Retrieve correct data points for this graph
            switch comp
                case 1
                    all_trials = squeeze(S_var_all(idx, :, :))';
                    theory_data = v_s_dot{idx};
                case 2
                    all_trials = squeeze(I_var_all(idx, :, :))';
                    theory_data = v_i_dot{idx};
                case 3
                    all_trials = squeeze(R_var_all(idx, :, :))';
                    theory_data = v_r_dot{idx};
            end

            % For each midpoint, calculate min/max over the interval
            for t = 1:length(midpoints)
                % Window
                t_mid = midpoints(t);
                t_min = t_mid - dt;
                t_max = t_mid + dt; 

                % Find the data points in [window)
                in_window = t_grid >= t_min & t_grid < t_max;
                if ~any(in_window)
                    continue;
                end
             
                vals = all_trials(:, in_window);  % [num_trials x num_timepoints_in_window]
                window_vals = vals(:);            % Flatten to a single vector

                % For each trial, get min and max of this window
                trial_mins = zeros(size(vals, 1), 1);
                trial_maxs = zeros(size(vals, 1), 1);

                for t_idx = 1:size(vals,1)
                    v = vals(t_idx, :);
                    if isempty(v)
                        trial_mins(t_idx) = Nan;
                        trial_maxs(t_idx) = NaN;
                    else
                        trial_mins(t_idx) = min(v);
                        trial_maxs(t_idx) = max(v);
                    end
                end

                % Remove NaNs in case any trials failed
                trial_mins = trial_mins(~isnan(trial_mins));
                trial_maxs = trial_maxs(~isnan(trial_maxs));

                if isempty(trial_mins)
                    fprintf("No valid points at midpoint %.2f\n", t_mid);
                end

                y_min_avg = mean(trial_mins);
                y_max_avg = mean(trial_maxs);

                % Plot error bar
                x = t_mid;
                plot([x,x], [y_min_avg, y_max_avg], 'Color', colors(idx,:), 'LineWidth', 1.5);
                y_theory = theory_data(t_grid == min(t_grid(t_grid >= t_mid)));
                plot(x, y_theory, '_', 'Color', colors(idx, :), 'MarkerSize', 12, 'LineWidth', 2);
            end
        
            % configure graphs
            xlabel('Time (days)');
            ylabel(['Variance']);
            title(['Simulated and Theoretical Variance for ', compartments{comp}, ', N = ', num2str(N)]);
            grid on;
            hold off;
        end
    end
end