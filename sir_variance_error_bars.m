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
    
    %% Agent-Based Simulations
    t_grid = linspace(0, max_T, 1000);
    S_var_all = zeros(length(Ns), length(t_grid));
    I_var_all = zeros(length(Ns), length(t_grid));
    R_var_all = zeros(length(Ns), length(t_grid));
    v_s_dot = cell(length(Ns), 1);
    v_i_dot = cell(length(Ns), 1);
    v_r_dot = cell(length(Ns), 1);
    
    %% Repeat simulation for each N
    for idx = 1:length(Ns)
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
        S_interp_var = interp1(T, S_var, t_grid, 'linear');
        I_interp_var = interp1(T, I_var, t_grid, 'linear');
        R_interp_var = interp1(T, R_var, t_grid, 'linear');
        S_var_all(idx, :) = S_interp_var;
        I_var_all(idx, :) = I_interp_var;
        R_var_all(idx, :) = R_interp_var;

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
                    sim_data = S_var_all(idx, :);
                    theory_data = v_s_dot{idx};
                case 2
                    sim_data = I_var_all(idx, :);
                    theory_data = v_i_dot{idx};
                case 3
                    sim_data = R_var_all(idx, :);
                    theory_data = v_r_dot{idx};
            end

            % For each midpoint, calculate min/max over the interval
            for k = 1:length(midpoints)
                % Window
                t_mid = midpoints(k);
                t_min = t_mid - dt;
                t_max = t_mid + dt; 

                % Find the data points in [window)
                in_window = t_grid >= t_min & t_grid < t_max;
                window_vals = sim_data(in_window);
                if isempty(window_vals)
                    continue; % No bars for empty interval
                end

                % Draw error bar between min and max centered at midpoint
                x = t_mid;
                y_min = min(window_vals);
                y_max = max(window_vals);
                %disp([x, y_min*10e4, y_max*10e4]);
                plot([x, x], [y_min, y_max], 'Color', colors(idx,:), 'LineWidth', 1.5);

                % Draw dash at theoretical line intercept
                y_theory = theory_data(t_grid == min(t_grid(t_grid >= t_mid)));
                plot(x, y_theory, '_', 'Color', colors(idx,:), 'MarkerSize', 12, 'LineWidth', 2);
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
