function SIHRS_hospitalized_washingtonparamsBUT_smallN()
% Simple SIHRS model with Washington MS parameters

    % Population and initial conditions
    N = 5591;
    
    % Initial proportions based on real washington MS data 0.000894884466408
    i0 = 0.000111304038111;
    h0 = 0.000894884466408;
    r0 = 0.0;
    d0 = 0.0;
    s0 = 1.0 - (i0 + h0 + r0 + d0);  
    
    % Model parameters - Washington MS parameters
    params = struct(...
        'beta', 0.194, 'gamma', 0.165, 'alpha', 0.111, 'lambda', 0.0083, ...
        'pSI', 1.00, 'pII', 0.0, 'pIH', 0.1614, 'pIR', 0.8367, 'pID', 0.0019, ...
        'pHH', 0.00, 'pHR', 0.846, 'pHD', 0.154, 'pRR', 0.02, 'pRS', 0.98, ...
        'tmax', 505, ...
        's0', s0, 'i0', i0, 'h0', h0, 'r0', r0, 'd0', d0 ...
    );

    
    % Run simulations
    num_simulations = 9;
    fprintf('Running %d simulations for N = %d (Washington MS params)...\n', num_simulations, N);
    
    all_results = cell(num_simulations, 1);
    for sim_idx = 1:num_simulations
        all_results{sim_idx} = sihrs_agent_model(N, params);
    end
    
    % Plot results
    plot_multiple_simulations_washington(all_results, N, params);
end


function result = sihrs_agent_model(N, params)
    % SIHRS agent-based stochastic model with death
    % Initial conditions
    s0 = round(params.s0 * N);
    i0 = round(params.i0 * N);
    h0 = round(params.h0 * N);
    r0 = round(params.r0 * N);
    d0 = round(params.d0 * N);

    % Adjust for rounding errors
    total = s0 + i0 + h0 + r0 + d0;
    if total ~= N
        compartments = [s0, i0, h0, r0, d0];
        [~, largest_idx] = max(compartments);
        compartments(largest_idx) = compartments(largest_idx) + (N - total);
        s0 = compartments(1);
        i0 = compartments(2);
        h0 = compartments(3);
        r0 = compartments(4);
        d0 = compartments(5);
    end

    if (s0 + i0 + h0 + r0 + d0) ~= N
        error('Initial conditions must sum to N');
    end

    max_events = N * 30;
    T = zeros(max_events, 1);
    I_prop = zeros(max_events, 1);
    I_count = zeros(max_events, 1);
    H_count = zeros(max_events, 1);
    D_count = zeros(max_events, 1);

    % Initialize agent lists
    S = (1:s0)';
    I = ((s0+1):(s0+i0))';
    H = ((s0+i0+1):(s0+i0+h0))';
    R = ((s0+i0+h0+1):(s0+i0+h0+r0))';
    D = ((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))';

    % Initialize time and event tracking
    t = 0.0;
    T(1) = 0.0;
    event_count = 1;

    % Record initial state
    total_pop = s0 + i0 + h0 + r0 + d0;
    I_prop(1) = i0 / total_pop;
    I_count(1) = i0;
    H_count(1) = h0;
    D_count(1) = d0;

    % Main simulation loop
    while ~isempty(I) && t < params.tmax
        nS = length(S);
        nI = length(I);
        nH = length(H);
        nR = length(R);

        infection_rate = params.pSI * params.beta * nS * nI / N;
        to_susceptible_from_R_rate = params.pRS * params.lambda * nR;
        to_hospital_rate = params.gamma * nI * params.pIH;
        to_recovered_from_I_rate = params.gamma * nI * params.pIR;
        to_dead_from_I_rate = params.gamma * nI * params.pID;
        to_recovered_from_H_rate = params.alpha * nH * params.pHR;
        to_dead_from_H_rate = params.alpha * nH * params.pHD;

        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate + ...
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate + ...
                     to_dead_from_H_rate;

        if total_rate == 0
            break;
        end

        % Sample time to next event
        dt = exprnd(1 / total_rate);
        t = t + dt;

        if t > params.tmax
            t = params.tmax;
            event_count = event_count + 1;
            T(event_count) = t;
            current_total = length(S) + length(I) + length(H) + length(R) + length(D);
            I_prop(event_count) = length(I) / current_total;
            I_count(event_count) = length(I);
            H_count(event_count) = length(H);
            D_count(event_count) = length(D);
            break;
        end

        event_count = event_count + 1;
        T(event_count) = t;

        % Determine which event occurs
        chance = rand() * total_rate;
        if chance < infection_rate
            % S to I transition
            if nS > 0
                num = randi(nS);
                infected_agent = S(num);
                S(num) = [];
                I = [I; infected_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate)
            % R to S transition
            if nR > 0
                num = randi(nR);
                susceptible_agent = R(num);
                R(num) = [];
                S = [S; susceptible_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate)
            % I to H transition
            if nI > 0
                num = randi(nI);
                hospitalized_agent = I(num);
                I(num) = [];
                H = [H; hospitalized_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate)
            % I to R transition
            if nI > 0
                num = randi(nI);
                recovered_agent = I(num);
                I(num) = [];
                R = [R; recovered_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate)
            % I to D transition
            if nI > 0
                num = randi(nI);
                dead_agent = I(num);
                I(num) = [];
                D = [D; dead_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate)
            % H to R transition
            if nH > 0
                num = randi(nH);
                recovered_agent = H(num);
                H(num) = [];
                R = [R; recovered_agent];
            end
        else
            % H to D transition
            if nH > 0
                num = randi(nH);
                dead_agent = H(num);
                H(num) = [];
                D = [D; dead_agent];
            end
        end

        % Record current state
        current_total = length(S) + length(I) + length(H) + length(R) + length(D);
        I_prop(event_count) = length(I) / current_total;
        I_count(event_count) = length(I);
        H_count(event_count) = length(H);
        D_count(event_count) = length(D);
    end

    % Trim arrays to actual size
    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    I_count = I_count(1:event_count);
    H_count = H_count(1:event_count);
    D_count = D_count(1:event_count);

    % Return result structure
    result = struct(...
        'N', N, ...
        'T', T, ...
        'I_count', I_count, ...
        'H_count', H_count, ...
        'D_count', D_count, ...
        'final_time', t, ...
        'peak_infected', max(I_count), ...
        'peak_time', T(argmax(I_count)), ...
        'peak_infected_prop', max(I_prop), ...
        'peak_time_prop', T(argmax(I_prop)), ...
        'i_inf', I_prop(end) ...
    );
end

function plot_multiple_simulations_washington(all_results, N, params)
    t_grid = (0:params.tmax)';
    all_interp_H = zeros(length(all_results), length(t_grid));
    all_interp_D = zeros(length(all_results), length(t_grid));

    for i = 1:length(all_results)
        res = all_results{i};

        if length(res.T) > 1
            itp_H = griddedInterpolant(res.T, res.H_count, 'linear', 'none');
            itp_D = griddedInterpolant(res.T, res.D_count, 'linear', 'none');

            for j = 1:length(t_grid)
                if t_grid(j) >= min(res.T) && t_grid(j) <= max(res.T)
                    all_interp_H(i, j) = itp_H(t_grid(j));
                    all_interp_D(i, j) = itp_D(t_grid(j));
                end
            end
        end
    end

    % Calculate statistics
    window = 14;
    all_active_D = zeros(size(all_interp_D));
    for i = 1:size(all_interp_D, 1)
        for t = 1:length(t_grid)
            if t <= window
                all_active_D(i, t) = all_interp_D(i, t);
            else
                all_active_D(i, t) = all_interp_D(i, t) - all_interp_D(i, t-window);
            end
        end
    end

    valid_sims = 1:size(all_interp_H, 1);
    fprintf('Using all %d simulations for prediction interval calculation\n', length(valid_sims));

    mean_H = mean(all_interp_H(valid_sims, :), 1)';
    lower_H = quantile(all_interp_H(valid_sims, :), 0.05, 1)';
    upper_H = quantile(all_interp_H(valid_sims, :), 0.95, 1)';

    population = N;
    simulation_start_date = datetime('2020-08-30');

    all_interp_H_prop = all_interp_H / N;
    mean_H_prop = mean_H / N;
    lower_H_prop = lower_H / N;
    upper_H_prop = upper_H / N;

    % First plot: Bandwidth
    figure;
    fill([t_grid; flipud(t_grid)], [upper_H_prop; flipud(lower_H_prop)], ...
         [0.7 0.9 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington MS Parameters - Small N Simulation');
    xlim([0, params.tmax]);
    ylim([0, max(upper_H_prop) * 1.1]);

    tick_interval = 90;
    xtick_positions = 0:tick_interval:params.tmax;
    xtick_dates = simulation_start_date + days(xtick_positions);
    date_labels = cellstr(datestr(xtick_dates, 'mm/dd/yy'));
    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_SmallN_Hospitalization_bandwidth.png');

    % Second plot: Trajectories
    figure;

    for i = 1:size(all_interp_H_prop, 1)
        if max(all_interp_H_prop(i, :)) > min(all_interp_H_prop(i, :)) * 1.1
            plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
            hold on;
        end
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington MS Parameters - Small N Simulation');
    xlim([0, params.tmax]);
    ylim([0, max(all_interp_H_prop(:)) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    % Create legend handles in correct order to fix colors
    h1 = plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2.5);
    legend(h1, {'Stochastic Simulations'}, 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_SmallN_Hospitalization_trajectories.png');

    % Third plot: ODE vs Stochastic Comparison
    figure;
    
    % Plot stochastic simulations (lighter blue)
    for i = 1:size(all_interp_H_prop, 1)
        if max(all_interp_H_prop(i, :)) > min(all_interp_H_prop(i, :)) * 1.1
            plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.2], 'LineWidth', 0.8);
            hold on;
        end
    end
    
    % Compute and plot ODE solution
    ode_H = solve_ode_hospitalization(params, N);
    ode_H_prop = ode_H / N;
    plot(t_grid, ode_H_prop, 'g-', 'LineWidth', 3.0);

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington MS Parameters - Small N: ODE vs Stochastic');
    xlim([0, params.tmax]);
    ylim([0, max([max(all_interp_H_prop(:)); max(ode_H_prop)]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    % Create legend handles in correct order
    h1 = plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2.5);
    h2 = plot(NaN, NaN, 'g-', 'LineWidth', 3.0);
    legend([h1, h2], {'Stochastic Simulations', 'ODE Solution'}, 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_SmallN_Hospitalization_combined.png');
end

function ode_H = solve_ode_hospitalization(params, N)
    % Solve the deterministic ODE system for hospitalization
    tspan = [0, params.tmax];
    
    % Initial conditions
    S0 = params.s0 * N;
    I0 = params.i0 * N;
    H0 = params.h0 * N;
    R0 = params.r0 * N;
    D0 = params.d0 * N;
    
    y0 = [S0; I0; H0; R0; D0];
    
    % Solve ODE
    [t, y] = ode45(@(t, y) sihrs_ode(t, y, params, N), tspan, y0);
    
    % Interpolate H compartment to match time grid
    t_grid = (0:params.tmax)';
    ode_H = interp1(t, y(:, 3), t_grid, 'linear', 'extrap');
    
    % Ensure non-negative values
    ode_H = max(ode_H, 0);
end

function dydt = sihrs_ode(t, y, params, N)
    % SIHRS ODE system
    S = y(1);
    I = y(2);
    H = y(3);
    R = y(4);
    D = y(5);
    
    % Rates
    infection_rate = params.beta * params.pSI * S * I / N;
    to_susceptible_from_R_rate = params.lambda * params.pRS * R;
    to_hospital_rate = params.gamma * params.pIH * I;
    to_recovered_from_I_rate = params.gamma * params.pIR * I;
    to_dead_from_I_rate = params.gamma * params.pID * I;
    to_recovered_from_H_rate = params.alpha * params.pHR * H;
    to_dead_from_H_rate = params.alpha * params.pHD * H;
    
    % Differential equations
    dSdt = -infection_rate + to_susceptible_from_R_rate;
    dIdt = infection_rate - to_hospital_rate - to_recovered_from_I_rate - to_dead_from_I_rate;
    dHdt = to_hospital_rate - to_recovered_from_H_rate - to_dead_from_H_rate;
    dRdt = to_recovered_from_I_rate + to_recovered_from_H_rate - to_susceptible_from_R_rate;
    dDdt = to_dead_from_I_rate + to_dead_from_H_rate;
    
    dydt = [dSdt; dIdt; dHdt; dRdt; dDdt];
end

function idx = argmax(x)
    [~, idx] = max(x);
end
