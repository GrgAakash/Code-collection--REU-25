function sihmd_variance_error_bars() 
    %% Model Parameters
    constants.max_t = 30;
    constants.susceptible_event_rate = 1.5; 
    constants.s_i_prob = 1;

    constants.infected_event_rate = 0.5;
    constants.i_m_prob = 0.095; 
    constants.i_d_prob = 0.005;
    constants.i_h_prob = 0.005;
    constants.i_s_prob = 1 - constants.i_m_prob - constants.i_h_prob - constants.i_d_prob;
    
    constants.hospitalized_event_rate = 1; 
    constants.h_m_prob = 0.6;
    constants.h_d_prob = 0.15;
    constants.h_s_prob = 1 - constants.h_m_prob - constants.h_d_prob;
    
    s0_prop = 0.96;
    i0_prop = 0.04;
    m0_prop = 1 - s0_prop - i0_prop;
    
    Ns = [316, 1000, 3162, 10000];
    num_trials = 3; % for each N

    %% Data Storage Initialization
    t_grid = linspace(0, constants.max_t, 1000); % Common time grid for interpolation
    
    % Pre-allocate matrices to store interpolated variance data for all trials
    S_var_all = zeros(length(Ns), length(t_grid), num_trials);
    I_var_all = zeros(length(Ns), length(t_grid), num_trials);
    H_var_all = zeros(length(Ns), length(t_grid), num_trials);
    M_var_all = zeros(length(Ns), length(t_grid), num_trials);
    D_var_all = zeros(length(Ns), length(t_grid), num_trials);
    
    % Cell arrays to store theoretical variance data
    v_s_dot = cell(length(Ns), 1);
    v_i_dot = cell(length(Ns), 1);
    v_h_dot = cell(length(Ns), 1);
    v_m_dot = cell(length(Ns), 1);
    v_d_dot = cell(length(Ns), 1);

    %% Main Simulation Loop
    for idx = 1:length(Ns)
        constants.N = Ns(idx);
        fprintf('Running simulations for N = %d...\n', constants.N);
        
        for trial = 1:num_trials
            fprintf('  Trial %d/%d\n', trial, num_trials);
            
            % Initial populations for this trial
            I0 = round(i0_prop * constants.N);
            M0 = round(m0_prop * constants.N);
            S0 = constants.N - I0 - M0;

            S = 1:S0;
            I = (S0+1):(S0+I0);
            H = [];
            M = (S0+I0+1):constants.N;
            D = [];

            % Initialize data structure for this single trial
            data.ns = S0; data.ni = I0; data.nh = 0; data.nm = M0; data.nd = 0;
            data.vs = 0; data.vi = 0; data.vh = 0; data.vm = 0; data.vd = 0;
            data.T = 0;

            % Run the stochastic simulation for one trial
            data = simulate_ipc(S, I, H, M, D, constants, data);
            
            % Fix first variance value which is initialized to 0
            if numel(data.vs) > 1
                data.vs(1) = data.vs(2); data.vi(1) = data.vi(2); data.vh(1) = data.vh(2); 
                data.vm(1) = data.vm(2); data.vd(1) = data.vd(2);
            end

            % Interpolate results onto the common time grid
            S_var_all(idx, :, trial) = interp1(data.T, data.vs, t_grid, 'linear', 0);
            I_var_all(idx, :, trial) = interp1(data.T, data.vi, t_grid, 'linear', 0);
            H_var_all(idx, :, trial) = interp1(data.T, data.vh, t_grid, 'linear', 0);
            M_var_all(idx, :, trial) = interp1(data.T, data.vm, t_grid, 'linear', 0);
            D_var_all(idx, :, trial) = interp1(data.T, data.vd, t_grid, 'linear', 0);
        end

        %% Deterministic ODE Model for Theoretical Variance
        N = Ns(idx);
        y0 = [s0_prop; i0_prop; 0; m0_prop; 0]; 
        [~, y] = ode45(@(t, y) sihmd_system(t, y, constants), t_grid, y0);
        s = y(:, 1);
        i = y(:, 2);
        h = y(:, 3);
    
        % Compute instantaneous variance for ODE system based on provided formulas
        v_s_dot{idx} = (constants.susceptible_event_rate*s.*i*constants.s_i_prob + constants.infected_event_rate*i*constants.i_s_prob + constants.hospitalized_event_rate*h*constants.h_s_prob) / N;
        v_i_dot{idx} = (constants.susceptible_event_rate*s.*i*constants.s_i_prob + constants.infected_event_rate*i) / N;
        v_h_dot{idx} = (constants.infected_event_rate*i*constants.i_h_prob + constants.hospitalized_event_rate*h) / N;
        v_m_dot{idx} = (constants.infected_event_rate*i*constants.i_m_prob + constants.hospitalized_event_rate*h*constants.h_m_prob) / N;
        v_d_dot{idx} = (constants.infected_event_rate*i*constants.i_d_prob + constants.hospitalized_event_rate*h*constants.h_d_prob) / N;
    end
    
    %% Graphing
    plot_variance_error_bars(Ns, t_grid, S_var_all, I_var_all, H_var_all, M_var_all, D_var_all, ...
        v_s_dot, v_i_dot, v_h_dot, v_m_dot, v_d_dot);
end

%% Simulation Core
function data = simulate_ipc(S, I, H, M, D, constants, data) 
    ns = numel(S);
    ni = numel(I);
    nh = numel(H);

    % Calculate total event rates for the master clock
    s_clock = constants.susceptible_event_rate * ns * ni / constants.N;
    i_clock = constants.infected_event_rate * ni;
    h_clock = constants.hospitalized_event_rate * nh;
    master_clock = s_clock + i_clock + h_clock;

    % Base case for recursion: end if time limit is reached, no events can occur, or all compartments are empty
    if data.T(end) >= constants.max_t || master_clock == 0 || (ns == 0 && ni == 0 && nh == 0 && numel(M) == 0 && numel(D) == 0)
        if data.T(end) < constants.max_t
            data.T(end+1) = constants.max_t;
            % Append last known value to the end to make interpolation flat
            data.ns(end+1)=data.ns(end); data.ni(end+1)=data.ni(end); data.nh(end+1)=data.nh(end); data.nm(end+1)=data.nm(end); data.nd(end+1)=data.nd(end);
            data.vs(end+1)=data.vs(end); data.vi(end+1)=data.vi(end); data.vh(end+1)=data.vh(end); data.vm(end+1)=data.vm(end); data.vd(end+1)=data.vd(end);
        end
        return;
    end

    % Determine time until next event from an exponential distribution
    dt = exprnd(1 / master_clock);
    data.T(end+1) = data.T(end) + dt;

    % Record current state before the event
    data.ns(end+1) = ns; 
    data.ni(end+1) = ni; 
    data.nh(end+1) = nh; 
    data.nm(end+1) = numel(M); 
    data.nd(end+1) = numel(D);
    
    % Calculate and record instantaneous variance for each compartment
    v_s = (constants.susceptible_event_rate*ns/constants.N*ni/constants.N*constants.s_i_prob + constants.infected_event_rate*ni/constants.N*constants.i_s_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_s_prob) / constants.N;
    v_i = (constants.susceptible_event_rate*ns/constants.N*ni/constants.N*constants.s_i_prob + constants.infected_event_rate*ni/constants.N) / constants.N;
    v_h = (constants.infected_event_rate*ni/constants.N*constants.i_h_prob + constants.hospitalized_event_rate*nh/constants.N) / constants.N;
    v_m = (constants.infected_event_rate*ni/constants.N*constants.i_m_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_m_prob) / constants.N;
    v_d = (constants.infected_event_rate*ni/constants.N*constants.i_d_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_d_prob) / constants.N;
    data.vs(end+1) = v_s; data.vi(end+1) = v_i; data.vh(end+1) = v_h; data.vm(end+1) = v_m; data.vd(end+1) = v_d;

    % Determine which event occurs based on relative rates
    r = rand;
    if (r < s_clock / master_clock)
        if rand < constants.s_i_prob
            I(end+1) = S(end);
            S(end) = [];
        end
    elseif (r < (s_clock + i_clock) / master_clock)
        r2 = rand;
        if r2 < constants.i_h_prob
            H(end+1) = I(end); I(end) = [];
        elseif r2 < (constants.i_h_prob + constants.i_m_prob)
            M(end+1) = I(end); I(end) = [];
        elseif r2 < (constants.i_h_prob + constants.i_m_prob + constants.i_d_prob)
            D(end+1) = I(end); I(end) = []; 
        else
            S(end+1) = I(end); I(end) = [];
        end
    else
        r2 = rand;
        if r2 < constants.h_s_prob
            S(end+1) = H(end); H(end) = [];
        elseif r2 < (constants.h_s_prob + constants.h_m_prob)
            M(end+1) = H(end); H(end) = [];
        else
            D(end+1) = H(end); H(end) = [];
        end
    end

    % Recursive call to continue the simulation
    data = simulate_ipc(S, I, H, M, D, constants, data);
end

%% ODE System
function dydt = sihmd_system(~, y, constants)
    s = y(1);
    i = y(2);
    h = y(3);

    beta = constants.susceptible_event_rate;
    gamma = constants.infected_event_rate;
    alpha = constants.hospitalized_event_rate;
    
    ds = -beta * s * i + gamma * constants.i_s_prob * i + alpha * constants.h_s_prob * h;
    di = beta * s * i - gamma * i;
    dh = gamma * constants.i_h_prob * i - alpha * h;
    dm = gamma * constants.i_m_prob * i + alpha * constants.h_m_prob * h;
    dd = gamma * constants.i_d_prob * i + alpha * constants.h_d_prob * h;

    dydt = [ds; di; dh; dm; dd];
end

%% Plotting 
function plot_variance_error_bars(Ns, t_grid, S_var_all, I_var_all, H_var_all, M_var_all, D_var_all, ...
                                v_s_dot, v_i_dot, v_h_dot, v_m_dot, v_d_dot)
    dt = 0.4; % Half width for interval around each midpoint
    midpoints = 0.5:dt*2:max(t_grid);
    compartments = {'Susceptible', 'Infected', 'Hospitalized', 'Immune', 'Dead'};
    all_sim_vars = {S_var_all, I_var_all, H_var_all, M_var_all, D_var_all};
    all_theory_vars = {v_s_dot, v_i_dot, v_h_dot, v_m_dot, v_d_dot};
    colors = lines(length(Ns));
    
    for comp_idx = 1:numel(compartments)
        for n_idx = 1:numel(Ns)
            figure;
            hold on;
            grid on;
            N = Ns(n_idx);
            
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
                
                trial_mins = min(vals_in_window, [], 2);
                trial_maxs = max(vals_in_window, [], 2);

                y_min_avg = mean(trial_mins(~isnan(trial_mins)));
                y_max_avg = mean(trial_maxs(~isnan(trial_maxs)));

                if ~isnan(y_min_avg) && ~isnan(y_max_avg)
                    % Plot error bar for simulated variance range
                    plot([t_mid, t_mid], [y_min_avg, y_max_avg], 'Color', colors(n_idx,:), 'LineWidth', 1.5);
                    
                    % Find the theoretical value at the closest time point and plot it
                    [~, time_index] = min(abs(t_grid - t_mid));
                    y_theory = theory_data(time_index); 
                    plot(t_mid, y_theory, '_', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
                end
            end
        
            xlabel('Time (days)');
            ylabel('Variance');
            title(['Simulated and Theoretical Variance for ', compartments{comp_idx}, ', N = ', num2str(N)]);
            legend('Simulated Range', 'Theoretical Value', 'Location', 'best');
            hold off;
        end
    end
end
