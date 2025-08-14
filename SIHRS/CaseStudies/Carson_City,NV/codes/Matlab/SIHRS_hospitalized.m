function SIHRS_hospitalized()
% SIHRS model for Carson City (Mar 2020-Nov 2021) starting from Patient Zero with stochastic simulations
% MATLAB version of the Julia SIHRS_hospitalized.jl

    % Initialize variables
    N = 56000;  % Carson City, Nevada population
    s0 = 0.0;
    i0 = 0.0;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;

    % --- Load and process real Carson City data for initial conditions ---
    try
        % Load cases and deaths data
        data_table = readtable('carson_city_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_date = datetime('2020-03-25');

        % Find start index for March 25
        start_idx = find(data_table.date == start_date, 1);
        if isempty(start_idx)
            % Find closest date
            date_diffs = abs(datenum(data_table.date) - datenum(start_date));
            [~, start_idx] = min(date_diffs);
            warning('Exact start date not found. Using closest date: %s', ...
                    datestr(data_table.date(start_idx)));
        end

        % Use March 25 real data for initial conditions
        real_initial_infected = data_table.cases(start_idx);
        real_initial_dead = data_table.deaths(start_idx);

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = 0.0;  % No hospitalization data in March
        r0 = 0.0;  % No recovered data in March
        s0 = 1.0 - (i0 + h0 + r0 + d0);

        fprintf('March 25 initial conditions: I=%d, D=%d, H=%d, R=%d, S=%d\n', ...
                real_initial_infected, real_initial_dead, 0, 0, round(s0 * N));

    catch ME
        warning('Could not load Carson City real data: %s', ME.message);
        % Fallback to hardcoded values for March 25
        real_initial_infected = 3;  % March 25 had 3 cases
        real_initial_dead = 0;      % March 25 had 0 deaths

        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = 0.0;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);
    end

    % Define model parameters structure for SIHRS with death model
    params = struct(...
        'beta', 0.183,      ... % infection rate (β > 0) - DECREASED for later peak
        'gamma', 0.150,     ... % I transition rate (γ > 0) - DECREASED for later peak 
        'alpha', 0.1,       ... % H transition rate (α > 0)
        'lambda', 0.0083,   ... % R transition rate (Λ > 0)
        'pSI', 1.00,        ... % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0,         ... % probability of I to I (stay infected)
        'pIH', 0.04,        ... % probability of I to H 
        'pIR', 0.959,       ... % probability of I to R 
        'pID', 0.001,       ... % probability of I to D
        'pHH', 0.01,        ... % probability of H to H (stay hospitalized)
        'pHR', 0.9882,      ... % probability of H to R
        'pHD', 0.0018,      ... % probability of H to D
        'pRR', 0.02,        ... % probability of R to R (stay recovered)
        'pRS', 0.98,        ... % probability of R to S
        'tmax', 620,        ... % simulation end time (long enough for all data)
        's0', s0,           ... % initial susceptible proportion
        'i0', i0,           ... % initial infected proportion
        'h0', h0,           ... % initial hospitalized proportion
        'r0', r0,           ... % initial recovered proportion
        'd0', d0            ... % initial dead proportion
    );

    % Verify R0 calculation
    calculated_R0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated R0 = %.6f \n', calculated_R0);

    % Validate parameters
    validate_parameters(params);

    % Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-10
        error('Initial conditions must sum to 1');
    end

    num_simulations = 5;

    % Input validation
    if N <= 0
        error('Population size must be positive integer');
    end

    % Store results for all simulations
    all_results = cell(num_simulations, 1);

    % Run multiple simulations
    try
        fprintf('Running %d stochastic simulations for N = %d...\n', num_simulations, N);

        for sim_idx = 1:num_simulations
            fprintf('Running simulation %d/%d...\n', sim_idx, num_simulations);
            all_results{sim_idx} = sihrs_agent_model(N, params);
        end

        fprintf('All simulations completed!\n');

        % Plot all simulations together
        plot_multiple_simulations(all_results, N, params);

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

    % Initial conditions - using params values and ensuring they sum to N
    s0 = round(params.s0 * N);  % susceptible
    i0 = round(params.i0 * N);  % infected
    h0 = round(params.h0 * N);  % hospitalized
    r0 = round(params.r0 * N);  % recovered
    d0 = round(params.d0 * N);  % dead

    % Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0 + d0;
    if total ~= N
        % Add or subtract the difference from the largest compartment
        compartments = [s0, i0, h0, r0, d0];
        [~, largest_idx] = max(compartments);
        compartments(largest_idx) = compartments(largest_idx) + (N - total);
        s0 = compartments(1);
        i0 = compartments(2);
        h0 = compartments(3);
        r0 = compartments(4);
        d0 = compartments(5);
    end

    % Validate initial conditions sum to N
    if (s0 + i0 + h0 + r0 + d0) ~= N
        error('Initial conditions must sum to N');
    end

    % Preallocate arrays for better performance
    max_events = N * 10;  % Estimate maximum number of events
    T = zeros(max_events, 1);
    I_prop = zeros(max_events, 1);
    I_count = zeros(max_events, 1);
    H_count = zeros(max_events, 1);
    D_count = zeros(max_events, 1);

    % Initialize agent arrays
    S = (1:s0)';
    I = ((s0+1):(s0+i0))';
    H = ((s0+i0+1):(s0+i0+h0))';
    R = ((s0+i0+h0+1):(s0+i0+h0+r0))';
    D = ((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))';

    % Initialize time tracking
    t = 0.0;
    T(1) = 0.0;
    event_count = 1;

    % Initialize proportion tracking
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

        % Calculate event rates according to the mathematical model
        infection_rate = params.pSI * params.beta * nS * nI / N;  % S to I rate
        to_susceptible_from_R_rate = params.pRS * params.lambda * nR;  % R to S rate
        to_hospital_rate = params.gamma * nI * params.pIH;  % I to H rate
        to_recovered_from_I_rate = params.gamma * nI * params.pIR;  % I to R rate
        to_dead_from_I_rate = params.gamma * nI * params.pID;  % I to D rate
        to_recovered_from_H_rate = params.alpha * nH * params.pHR;  % H to R rate
        to_dead_from_H_rate = params.alpha * nH * params.pHD;  % H to D rate

        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate + ...
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate + ...
                     to_dead_from_H_rate;

        if total_rate == 0
            break;
        end

        % Time of next event
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

        % Update tracking arrays
        current_total = length(S) + length(I) + length(H) + length(R) + length(D);
        I_prop(event_count) = length(I) / current_total;
        I_count(event_count) = length(I);
        H_count(event_count) = length(H);
        D_count(event_count) = length(D);
    end

    % Trim unused preallocated space
    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    I_count = I_count(1:event_count);
    H_count = H_count(1:event_count);
    D_count = D_count(1:event_count);

    % Store results
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

function plot_multiple_simulations(all_results, N, params)
    % This function creates two figures:
    % 1. The uncertainty envelope plot (and shaded prediction interval).
    % 2. The histogram of peak infection values (like Figure 4B).

    % --- 1. Create a common time grid for comparison ---
    t_grid = (0:params.tmax)';  % daily time steps

    % --- 2. Interpolate each simulation's results onto the grid ---
    all_interp_H = zeros(length(all_results), length(t_grid));
    all_interp_D = zeros(length(all_results), length(t_grid));

    for i = 1:length(all_results)
        res = all_results{i};
        % Create interpolation function
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

    % --- Compute active deaths (rolling window) for each simulation ---
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

    % --- 3. Calculate statistics for the bandwidth (using all simulations) ---
    % Use all simulations without filtering
    valid_sims = 1:size(all_interp_H, 1);

    fprintf('Using all %d simulations for prediction interval calculation\n', length(valid_sims));

    % Calculate statistics using all simulations
    mean_H = mean(all_interp_H(valid_sims, :), 1)';
    lower_H = quantile(all_interp_H(valid_sims, :), 0.05, 1)';  % 5th percentile
    upper_H = quantile(all_interp_H(valid_sims, :), 0.95, 1)';  % 95th percentile
    mean_D = mean(all_active_D, 1)';
    lower_D = quantile(all_active_D, 0.05, 1)';
    upper_D = quantile(all_active_D, 0.95, 1)';

    % --- 4. Load and process real-world hospitalization data ---
    population = N;  % Use same population as simulation
    real_interp_H = zeros(length(t_grid), 1);
    real_interp_D_prop = zeros(length(t_grid), 1);
    simulation_start_date = datetime('2020-03-25');  % Simulation start date

    try
        % Load hospitalization data
        hosp_data_table = readtable('hospitalization_Carson_filtered_new.csv');

        % Convert collection_week to datetime with proper year handling
        hosp_data_table.collection_week = datetime(hosp_data_table.collection_week, 'InputFormat', 'M/d/yy');

        % Sort by date
        hosp_data_table = sortrows(hosp_data_table, 'collection_week');

        % Use the total hospitalization data
        hospitalization_data = hosp_data_table.total_adult_and_pediatric_covid_patients;
        hosp_dates = hosp_data_table.collection_week;

        % Convert to days from March 25 (simulation start)
        hosp_days = days(hosp_dates - simulation_start_date);

        % Interpolate hospitalization data to the simulation time grid
        real_interp_H = interp1(hosp_days, hospitalization_data, t_grid, 'linear', 0);

        % For death data, align with March 25 simulation start
        data_table = readtable('carson_city_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_idx = find(data_table.date == simulation_start_date, 1);
        if isempty(start_idx)
            date_diffs = abs(datenum(data_table.date) - datenum(simulation_start_date));
            [~, start_idx] = min(date_diffs);
        end
        dates_from_start = data_table.date(start_idx:end);
        real_interp_D = interp1(days(dates_from_start - simulation_start_date), ...
                               data_table.deaths(start_idx:end), t_grid, 'linear', 0);

        % Compute active deaths for real data
        real_active_D = zeros(length(real_interp_D), 1);
        for t = 1:length(real_interp_D)
            if t <= window
                real_active_D(t) = real_interp_D(t);
            else
                real_active_D(t) = real_interp_D(t) - real_interp_D(t-window);
            end
        end
        real_interp_D_prop = real_active_D / population;

        fprintf('Successfully loaded hospitalization data with %d data points\n', length(hospitalization_data));

    catch ME
        fprintf('Warning: Could not load or process hospitalization data: %s\n', ME.message);
        real_interp_H = zeros(length(t_grid), 1);
        real_interp_D_prop = zeros(length(t_grid), 1);
    end

    % Convert all counts to proportions for plotting
    all_interp_H_prop = all_interp_H / N;
    mean_H_prop = mean_H / N;
    lower_H_prop = lower_H / N;
    upper_H_prop = upper_H / N;
    real_interp_H_prop = real_interp_H / population;

    % --- 5. Create the final plot with the uncertainty envelope ---
    figure;
    fill([t_grid; flipud(t_grid)], [upper_H_prop; flipud(lower_H_prop)], ...
         [0.7 0.9 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;

    % Only plot real hospitalization data where it exists (not zero)
    valid_real_data = real_interp_H_prop > 0;
    if any(valid_real_data)
        plot(t_grid(valid_real_data), real_interp_H_prop(valid_real_data), ...
             'r-', 'LineWidth', 2.5);
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max([upper_H_prop; real_interp_H_prop]) * 1.1]);

    % Set x-ticks every 90 days to avoid overlapping
    tick_interval = 90;
    xtick_positions = 0:tick_interval:params.tmax;
    xtick_dates = simulation_start_date + days(xtick_positions);
    date_labels = cellstr(datestr(xtick_dates, 'mm/dd/yy'));
    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Hospitalization Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Hospitalization_bandwidth.png');

    % --- 6. Create a second figure with all stochastic sims and the real data ---
    figure;
    % Plot only non-constant stochastic simulations as semi-transparent blue lines
    for i = 1:size(all_interp_H_prop, 1)
        if max(all_interp_H_prop(i, :)) > min(all_interp_H_prop(i, :)) * 1.1  % Filter constant lines
            plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
            hold on;
        end
    end

    % Plot the real hospitalization data as a solid red line
    valid_real_data = real_interp_H_prop > 0;
    if any(valid_real_data)
        plot(t_grid(valid_real_data), real_interp_H_prop(valid_real_data), ...
             'r-', 'LineWidth', 2.5);
    end

    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('Carson City, NV');
    xlim([0, params.tmax]);
    ylim([0, max([max(all_interp_H_prop(:)); max(real_interp_H_prop)]) * 1.1]);

    % Set x-ticks
    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('Stochastic Simulations', 'Real Hospitalization Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Carson_City_Hospitalization_trajectories.png');
end

function idx = argmax(x)
    [~, idx] = max(x);
end
