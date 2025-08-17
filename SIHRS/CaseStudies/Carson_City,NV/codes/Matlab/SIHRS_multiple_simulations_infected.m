function SIHRS_multiple_simulations_infected()
% SIHRS model for Carson City (Mar 2020-Dec 2021) starting from Patient Zero with stochastic simulations

    % Initialize variables at function level
    N = 56000;  % Carson City, Nevada population
    s0 = 0.0;
    i0 = 0.0;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;

    % Load Carson City, Nevada data for initial conditions
    try

        data_table = readtable('carson_city_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_date = datetime('2020-03-25');


        start_idx = find(data_table.date == start_date, 1);
        if isempty(start_idx)

            date_diffs = abs(datenum(data_table.date) - datenum(start_date));
            [~, start_idx] = min(date_diffs);
            warning('Exact start date not found. Using closest date: %s', ...
                    datestr(data_table.date(start_idx)));
        end

        real_initial_infected = data_table.cases(start_idx);
        real_initial_dead = data_table.deaths(start_idx);

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = 0.0;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);

        fprintf('March 25 initial conditions: I=%d, D=%d, H=%d, R=%d, S=%d\n', ...
                real_initial_infected, real_initial_dead, 0, 0, round(s0 * N));

    catch ME
        warning('Could not load Carson City real data: %s', ME.message);

        real_initial_infected = 1;
        real_initial_dead = 0;

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = 0.0;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);
    end

    % Model parameters
    params = struct(...
        'beta', 0.157,      ... % infection rate (β > 0)
        'gamma', 0.127,     ... % I transition rate (γ > 0) 
        'alpha', 0.111,       ... % H transition rate (α > 0)
        'lambda', 0.0083,    ... % R transition rate (Λ > 0)  
        'pSI', 1.0,         ... % probability of S to I (p_{SI} in (0,1])
        'pII', 0.00,        ... % probability of I to I (stay infected)
        'pIH', 0.1060,        ... % probability of I to H 
        'pIR', 0.8921,       ... % probability of I to R 
        'pID', 0.0019,       ... % probability of I to D
        'pHH', 0.01,        ... % probability of H to H (stay hospitalized)
        'pHR', 0.836,      ... % probability of H to R
        'pHD', 0.154,      ... % probability of H to D
        'pRR', 0.02,        ... % probability of R to R (stay recovered)
        'pRS', 0.98,        ... % probability of R to S
        'tmax', 620,        ... % simulation end time (extended for Carson City, NV data)
        's0', s0,           ... % initial susceptible proportion
        'i0', i0,           ... % initial infected proportion
        'h0', h0,           ... % initial hospitalized proportion
        'r0', r0,           ... % initial recovered proportion
        'd0', d0            ... % initial dead proportion
    );




    calculated_R0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated R0 = %.6f \n', calculated_R0);


    validate_parameters(params);


    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-10
        error('Initial conditions must sum to 1');
    end

    num_simulations = 59;


    if N <= 0
        error('Population size must be positive integer');
    end


    all_results = cell(num_simulations, 1);


    try
        fprintf('Running %d stochastic simulations for N = %d...\n', num_simulations, N);

        for sim_idx = 1:num_simulations
            fprintf('Running simulation %d/%d...\n', sim_idx, num_simulations);
            all_results{sim_idx} = sihrs_agent_model_infected(N, params);
        end

        fprintf('All simulations completed!\n');


        plot_multiple_simulations_infected(all_results, N, params);

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

function result = sihrs_agent_model_infected(N, params)
    % SIHRS agent-based stochastic model with death (focusing on infected and deaths)
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
        D_count(event_count) = length(D);
    end


    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    I_count = I_count(1:event_count);
    D_count = D_count(1:event_count);


    result = struct(...
        'N', N, ...
        'T', T, ...
        'I_count', I_count, ...
        'D_count', D_count, ...
        'final_time', t, ...
        'peak_infected', max(I_count), ...
        'peak_time', T(argmax(I_count)), ...
        'peak_infected_prop', max(I_prop), ...
        'peak_time_prop', T(argmax(I_prop)), ...
        'i_inf', I_prop(end) ...
    );
end

function plot_multiple_simulations_infected(all_results, N, params)
    t_grid = (0:params.tmax)';
    all_interp_I = zeros(length(all_results), length(t_grid));
    all_interp_D = zeros(length(all_results), length(t_grid));

    for i = 1:length(all_results)
        res = all_results{i};

        if length(res.T) > 1
            itp_I = griddedInterpolant(res.T, res.I_count, 'linear', 'none');
            itp_D = griddedInterpolant(res.T, res.D_count, 'linear', 'none');

            for j = 1:length(t_grid)
                if t_grid(j) >= min(res.T) && t_grid(j) <= max(res.T)
                    all_interp_I(i, j) = itp_I(t_grid(j));
                    all_interp_D(i, j) = itp_D(t_grid(j));
                end
            end
        end
    end


    window = 14;
    all_active_D = zeros(size(all_interp_D));
    for i = 1:size(all_interp_D, 1)
        all_active_D(i, :) = compute_rolling_window(all_interp_D(i, :), window);
    end


    mean_I = mean(all_interp_I, 1)';
    lower_I = quantile(all_interp_I, 0.05, 1)';
    upper_I = quantile(all_interp_I, 0.95, 1)';
    mean_D = mean(all_active_D, 1)';
    lower_D = quantile(all_active_D, 0.05, 1)';
    upper_D = quantile(all_active_D, 0.95, 1)';


    population = N;  % Use same population as simulation
    real_interp_I = zeros(length(t_grid), 1);
    real_interp_D = zeros(length(t_grid), 1);
    real_interp_D_prop = zeros(length(t_grid), 1);
    real_cumulative_D = zeros(length(t_grid), 1);
    start_date_real = datetime('2020-03-25');  % Default start date

    try

        data_table = readtable('carson_city_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_idx = find(data_table.date == start_date_real, 1);
        if isempty(start_idx)
            date_diffs = abs(datenum(data_table.date) - datenum(start_date_real));
            [~, start_idx] = min(date_diffs);
        end

        dates_from_start = data_table.date(start_idx:end);
        cumulative_cases_from_start = data_table.cases(start_idx:end);
        recovery_days = 14;  % Longer recovery period for rural areas (14 days)


        cumulative_shifted = [zeros(recovery_days, 1); cumulative_cases_from_start(1:end-recovery_days)];
        active_cases_count = cumulative_cases_from_start - cumulative_shifted;

        active_cases_count = max(active_cases_count, 0);

        carson_days = days(dates_from_start - dates_from_start(1));


        real_interp_I = interp1(carson_days, active_cases_count, t_grid, 'linear', 0);


        daily_deaths_data = readtable('carson_city_daily_deaths.csv');
        daily_deaths_data.date = datetime(daily_deaths_data.date, 'InputFormat', 'yyyy-MM-dd');


        daily_start_idx = find(daily_deaths_data.date == start_date_real, 1);
        if isempty(daily_start_idx)
            date_diffs = abs(datenum(daily_deaths_data.date) - datenum(start_date_real));
            [~, daily_start_idx] = min(date_diffs);
        end


        daily_deaths_from_start = daily_deaths_data.daily_deaths(daily_start_idx:end);
        daily_death_dates = daily_deaths_data.date(daily_start_idx:end);


        daily_death_days = days(daily_death_dates - start_date_real);


        real_interp_D = interp1(daily_death_days, daily_deaths_from_start, t_grid, 'linear', 0);


        moving_avg_deaths = daily_deaths_data.moving_avg_7day(daily_start_idx:end);
        real_active_D = interp1(daily_death_days, moving_avg_deaths, t_grid, 'linear', 0);
        real_interp_D_prop = real_active_D / population;


        cumulative_deaths_from_start = daily_deaths_data.cumulative_deaths(daily_start_idx:end);
        real_cumulative_D = interp1(daily_death_days, cumulative_deaths_from_start, t_grid, 'linear', 0);

        fprintf('Successfully loaded real data with %d data points\n', length(active_cases_count));

    catch ME
        fprintf('Warning: Could not load or process real data: %s\n', ME.message);
        real_interp_I = zeros(length(t_grid), 1);
        real_interp_D_prop = zeros(length(t_grid), 1);
    end


    population_factor = 1.0 / N;  % Calculate division factor once

    all_interp_I_prop = all_interp_I * population_factor;
    mean_I_prop = mean_I * population_factor;
    lower_I_prop = lower_I * population_factor;
    upper_I_prop = upper_I * population_factor;
    all_active_D_prop = all_active_D * population_factor;
    mean_D_prop = mean_D * population_factor;
    lower_D_prop = lower_D * population_factor;
    upper_D_prop = upper_D * population_factor;
    real_interp_I_prop = real_interp_I * population_factor;


    tick_interval = 90;
    xtick_positions = 0:tick_interval:params.tmax;
    xtick_dates = start_date_real + days(xtick_positions);
    date_labels = cellstr(datestr(xtick_dates, 'mm/dd/yy'));


    figure;
    fill([t_grid; flipud(t_grid)], [upper_I_prop; flipud(lower_I_prop)], ...
         [0.7 0.9 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    plot(t_grid, real_interp_I_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Time (days)');
    ylabel('Infected Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max([upper_I_prop; real_interp_I_prop]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Full_Pandemic_bandwidth.png');


    figure;

    for i = 1:size(all_interp_I_prop, 1)
        plot(t_grid, all_interp_I_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
        hold on;
    end


    plot(t_grid, real_interp_I_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Time (days)');
    ylabel('Infected Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max([max(max(all_interp_I_prop)), max(real_interp_I_prop)]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    % Custom legend
    plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2.5);
    legend('Stochastic Simulations', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Full_Pandemic_trajectories.png');

    % --- 7. Create a figure for death proportion (active deaths) ---
    figure;
    fill([t_grid; flipud(t_grid)], [upper_D_prop; flipud(lower_D_prop)], ...
         [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    plot(t_grid, real_interp_D_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Date');
    ylabel('Active Death Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max([upper_D_prop; real_interp_D_prop]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Full_Pandemic_active_deaths.png');

    % --- 8. Create a figure for cumulative death proportion ---
    real_interp_D_prop_cumulative = real_cumulative_D / population;
    figure;
    plot(t_grid, real_interp_D_prop_cumulative, 'r-', 'LineWidth', 2.5);
    xlabel('Date');
    ylabel('Cumulative Death Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max(real_interp_D_prop_cumulative) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Full_Pandemic_cumulative_deaths.png');
end

% Helper function for rolling window calculation
function result = compute_rolling_window(data, window_size)
    result = data;  % MATLAB automatically copies arrays
    for t = (window_size + 1):length(data)
        result(t) = data(t) - data(t - window_size);
    end
end

function idx = argmax(x)
    [~, idx] = max(x);
end
