function SIHRS_hospitalized_aug30()
% SIHRS model for Washington, Mississippi (Aug 2020-Dec 2021) starting from August 30, 2020 with stochastic simulations

    % Initialize variables at function level
    N = 44922;  % Washington County, Mississippi population (2020 Census)
    s0 = 0.0;
    i0 = 0.0;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;

    % Load Washington County, Mississippi data for initial conditions
    try

        data_table = readtable('washington_mississippi_active_cases_aug30.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_date = datetime('2020-08-30');


        start_idx = find(data_table.date == start_date, 1);
        if isempty(start_idx)

            date_diffs = abs(datenum(data_table.date) - datenum(start_date));
            [~, start_idx] = min(date_diffs);
            warning('Exact start date not found. Using closest date: %s', ...
                    datestr(data_table.date(start_idx)));
        end


        real_initial_infected = data_table.active_cases(start_idx);
        real_initial_dead = data_table.cumulative_deaths(start_idx);
        

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = real_initial_hospitalized / N;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);

        fprintf('August 30 initial conditions: I=%d, D=%d, H=%.1f, R=%d, S=%d\n', ...
                real_initial_infected, real_initial_dead, real_initial_hospitalized, 0, round(s0 * N));

    catch ME
        warning('Could not load Washington County, MS real data: %s', ME.message);
        real_initial_infected = 166.0;
        real_initial_dead = 60;
        real_initial_hospitalized = 40.2;  

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = real_initial_hospitalized / N;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);
    end

    % Model parameters
    params = struct(...
        'beta', 0.194,      ... % infection rate (β > 0) - Updated for Washington, MS
        'gamma', 0.165,     ... % I transition rate (γ > 0) - Updated for Washington, MS
        'alpha', 0.111,       ... % H transition rate (α > 0)
        'lambda', 0.0083,   ... % R transition rate (Λ > 0) - Updated for Washington, MS
        'pSI', 1.00,        ... % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0,         ... % probability of I to I (stay infected)
        'pIH', 0.1614,      ... % probability of I to H - Updated from P(IH) calculation
        'pIR', 0.8367,      ... % probability of I to R - Updated to sum to 1
        'pID', 0.0019,      ... % probability of I to D - Updated from P(ID) calculation
        'pHH', 0.00,        ... % probability of H to H (stay hospitalized) - Updated
        'pHR', 0.846,       ... % probability of H to R - Updated
        'pHD', 0.154,       ... % probability of H to D - Updated
        'pRR', 0.02,        ... % probability of R to R (stay recovered)
        'pRS', 0.98,        ... % probability of R to S
        'tmax', 620,        ... % simulation end time (extended for Washington, MS data)
        's0', s0,           ... % initial susceptible proportion
        'i0', i0,           ... % initial infected proportion
        'h0', h0,           ... % initial hospitalized proportion
        'r0', r0,           ... % initial recovered proportion
        'd0', d0            ... % initial dead proportion
    );


    calculated_R0 = (params.beta * params.pSI) / (params.gamma * (1 - params.pII));
    fprintf('Calculated R0 = %.6f \n', calculated_R0);


    validate_parameters(params);


    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-10
        error('Initial conditions must sum to 1');
    end

    num_simulations = 9;


    if N <= 0
        error('Population size must be positive integer');
    end


    all_results = cell(num_simulations, 1);


    try
        fprintf('Running %d stochastic simulations for N = %d...\n', num_simulations, N);

        for sim_idx = 1:num_simulations
            fprintf('Running simulation %d/%d...\n', sim_idx, num_simulations);
            all_results{sim_idx} = sihrs_agent_model(N, params);
        end

        fprintf('All simulations completed!\n');


        plot_multiple_simulations_aug30(all_results, N, params);

    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
end

function validate_parameters(params)
    % Validate rates are positive
    if any([params.beta, params.gamma, params.alpha, params.lambda] <= 0)
        error('All rates (beta, gamma, alpha, lambda) must be positive');
    end

    % Validate probabilities are in [0,1]
    probs = [params.pSI, params.pII, params.pIH, params.pIR, params.pID, ...
             params.pHH, params.pHR, params.pHD, params.pRR, params.pRS];
    if any(probs < 0 | probs > 1)
        error('All probabilities must be in [0,1]');
    end

    % Validate probability sums
    if abs((params.pII + params.pIH + params.pIR + params.pID) - 1.0) > 1e-10
        error('I transition probabilities must sum to 1');
    end
    if abs((params.pHH + params.pHR + params.pHD) - 1.0) > 1e-10
        error('H transition probabilities must sum to 1');
    end
    if abs((params.pRR + params.pRS) - 1.0) > 1e-10
        error('R transition probabilities must sum to 1');
    end
end

function result = sihrs_agent_model(N, params)
    % SIHRS agent-based stochastic model with death
    % Initial conditions
    s0 = round(params.s0 * N);
    i0 = round(params.i0 * N);
    h0 = round(params.h0 * N);
    r0 = round(params.r0 * N);
    d0 = round(params.d0 * N);


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


    S = (1:s0)';
    I = ((s0+1):(s0+i0))';
    H = ((s0+i0+1):(s0+i0+h0))';
    R = ((s0+i0+h0+1):(s0+i0+h0+r0))';
    D = ((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))';


    t = 0.0;
    T(1) = 0.0;
    event_count = 1;


    total_pop = s0 + i0 + h0 + r0 + d0;
    I_prop(1) = i0 / total_pop;
    I_count(1) = i0;
    H_count(1) = h0;
    D_count(1) = d0;


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


        current_total = length(S) + length(I) + length(H) + length(R) + length(D);
        I_prop(event_count) = length(I) / current_total;
        I_count(event_count) = length(I);
        H_count(event_count) = length(H);
        D_count(event_count) = length(D);
    end


    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    I_count = I_count(1:event_count);
    H_count = H_count(1:event_count);
    D_count = D_count(1:event_count);


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

function plot_multiple_simulations_aug30(all_results, N, params)
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
    mean_D = mean(all_active_D, 1)';
    lower_D = quantile(all_active_D, 0.05, 1)';
    upper_D = quantile(all_active_D, 0.95, 1)';

    population = N;
    real_interp_H = zeros(length(t_grid), 1);
    real_interp_D_prop = zeros(length(t_grid), 1);
    simulation_start_date = datetime('2020-08-30');

    try

        hosp_data_table = readtable('hospitalization_MS_filtered.csv');


        hosp_data_table.collection_week = datetime(hosp_data_table.collection_week, 'InputFormat', 'M/d/yy');
        

        for i = 1:height(hosp_data_table)
            if year(hosp_data_table.collection_week(i)) < 2000
                hosp_data_table.collection_week(i) = hosp_data_table.collection_week(i) + years(2000);
            end
        end


        [unique_dates, ~, idx] = unique(hosp_data_table.collection_week);
        aggregated_data = accumarray(idx, hosp_data_table.total_adult_and_pediatric_covid_patients, [], @sum);
        

        [hosp_dates, sort_idx] = sort(unique_dates);
        hospitalization_data = aggregated_data(sort_idx);


        hosp_days = days(hosp_dates - simulation_start_date);


        real_interp_H = interp1(hosp_days, hospitalization_data, t_grid, 'linear', NaN);
        

        real_interp_H_for_plot = real_interp_H;
        real_interp_H(isnan(real_interp_H)) = 0;


        data_table = readtable('washington_mississippi_combined_aug30.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_idx = find(data_table.date == simulation_start_date, 1);
        if isempty(start_idx)
            date_diffs = abs(datenum(data_table.date) - datenum(simulation_start_date));
            [~, start_idx] = min(date_diffs);
        end
        dates_from_start = data_table.date(start_idx:end);
        real_interp_D = interp1(days(dates_from_start - simulation_start_date), ...
                               data_table.deaths(start_idx:end), t_grid, 'linear', 0);


        real_active_D = zeros(length(real_interp_D), 1);
        for t = 1:length(real_interp_D)
            if t <= window
                real_active_D(t) = real_interp_D(t);
            else
                real_active_D(t) = real_interp_D(t) - real_interp_D(t-window);
            end
        end
        real_interp_D_prop = real_active_D / population;

        fprintf('Successfully loaded hospitalization data with %d unique dates and %d total data points\n', length(hospitalization_data), height(hosp_data_table));
        fprintf('Hospitalization data range: %.1f to %.1f patients\n', min(hospitalization_data), max(hospitalization_data));
        fprintf('Date range: %s to %s\n', datestr(min(hosp_dates)), datestr(max(hosp_dates)));

    catch ME
        fprintf('Warning: Could not load or process hospitalization data: %s\n', ME.message);
        real_interp_H = zeros(length(t_grid), 1);
        real_interp_D_prop = zeros(length(t_grid), 1);
    end


    all_interp_H_prop = all_interp_H / N;
    mean_H_prop = mean_H / N;
    lower_H_prop = lower_H / N;
    upper_H_prop = upper_H / N;
    real_interp_H_prop = real_interp_H / population;


    figure;
    fill([t_grid; flipud(t_grid)], [upper_H_prop; flipud(lower_H_prop)], ...
         [0.7 0.9 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;


    valid_real_data = ~isnan(real_interp_H_for_plot / population);
    if any(valid_real_data)
        valid_H_prop = real_interp_H_for_plot(valid_real_data) / population;
        plot(t_grid(valid_real_data), valid_H_prop, ...
             'r-', 'LineWidth', 2.5);
        fprintf('Plotting real hospitalization data: %d valid points, range %.6f to %.6f\n', ...
                sum(valid_real_data), min(valid_H_prop), max(valid_H_prop));
    else
        fprintf('Warning: No valid real hospitalization data to plot\n');
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington, Mississippi - Aug 30 Start');
    xlim([0, params.tmax]);

    if any(valid_real_data)
        max_real_H = max(real_interp_H_for_plot(valid_real_data) / population);
        ylim([0, max([upper_H_prop; max_real_H]) * 1.1]);
    else
        ylim([0, max(upper_H_prop) * 1.1]);
    end


    tick_interval = 90;
    xtick_positions = 0:tick_interval:params.tmax;
    xtick_dates = simulation_start_date + days(xtick_positions);
    date_labels = cellstr(datestr(xtick_dates, 'mm/dd/yy'));
    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Hospitalization Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Aug30_Hospitalization_bandwidth.png');


    figure;

    for i = 1:size(all_interp_H_prop, 1)
        if max(all_interp_H_prop(i, :)) > min(all_interp_H_prop(i, :)) * 1.1
            plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
            hold on;
        end
    end


    valid_real_data = ~isnan(real_interp_H_for_plot / population);
    if any(valid_real_data)
        valid_H_prop = real_interp_H_for_plot(valid_real_data) / population;
        plot(t_grid(valid_real_data), valid_H_prop, ...
             'r-', 'LineWidth', 2.5);
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington, Mississippi - Aug 30 Start');
    xlim([0, params.tmax]);

    if any(valid_real_data)
        max_real_H = max(real_interp_H_for_plot(valid_real_data) / population);
        ylim([0, max([max(all_interp_H_prop(:)); max_real_H]) * 1.1]);
    else
        ylim([0, max(all_interp_H_prop(:)) * 1.1]);
    end


    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('Stochastic Simulations', 'Real Hospitalization Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Aug30_Hospitalization_trajectories.png');

    % Third figure: ODE + Real World + Stochastic Simulations Combined
    figure;
    
    % Plot stochastic simulations (lighter blue)
    for i = 1:size(all_interp_H_prop, 1)
        if max(all_interp_H_prop(i, :)) > min(all_interp_H_prop(i, :)) * 1.1
            plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.2], 'LineWidth', 0.8);
            hold on;
        end
    end
    
    % Compute and plot ODE solution
    ode_H = solve_ode_hospitalization_aug30(params, N);
    ode_H_prop = ode_H / N;
    plot(t_grid, ode_H_prop, 'g-', 'LineWidth', 3.0);
    
    % Plot real world data
    if any(valid_real_data)
        valid_H_prop = real_interp_H_for_plot(valid_real_data) / population;
        plot(t_grid(valid_real_data), valid_H_prop, ...
             'r-', 'LineWidth', 2.5);
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Washington, MS - Aug 30 Start: ODE vs Stochastic vs Real Data');
    xlim([0, params.tmax]);

    if any(valid_real_data)
        max_real_H = max(real_interp_H_for_plot(valid_real_data) / population);
        ylim([0, max([max(all_interp_H_prop(:)); max(ode_H_prop); max_real_H]) * 1.1]);
    else
        ylim([0, max([max(all_interp_H_prop(:)); max(ode_H_prop)]) * 1.1]);
    end

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    % Create legend handles in correct order
    h1 = plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2.5);
    h2 = plot(NaN, NaN, 'g-', 'LineWidth', 3.0);
    h3 = plot(NaN, NaN, 'r-', 'LineWidth', 2.5);
    legend([h1, h2, h3], {'Stochastic Simulations', 'ODE Solution', 'Real Hospitalization Data'}, 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Aug30_Hospitalization_combined.png');
end

function ode_H = solve_ode_hospitalization_aug30(params, N)
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
    [t, y] = ode45(@(t, y) sihrs_ode_aug30(t, y, params, N), tspan, y0);
    
    % Interpolate H compartment to match time grid
    t_grid = (0:params.tmax)';
    ode_H = interp1(t, y(:, 3), t_grid, 'linear', 'extrap');
    
    % Ensure non-negative values
    ode_H = max(ode_H, 0);
end

function dydt = sihrs_ode_aug30(t, y, params, N)
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
