function SIHRS_multiple_simulations_infected()
% SIHRS model for Washington, Mississippi (Mar 2020-Dec 2021) 
% starting from Patient Zero with stochastic simulations
% Focuses on infected cases and deaths (like the Julia version)
% MATLAB version of SIHRS_multiple_simulations_infected.jl

    % Initialize variables at function level
    N = 43000;  % Washington County, Mississippi population (2020 Census)
    
    % Initialize default values for initial conditions
    s0 = 0.0;
    i0 = 0.0;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;
    
    % --- Load and process real Washington County, MS data for initial conditions ---
    try
        % Load cases and deaths data
        data_table = readtable('washington_mississippi_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_date = datetime('2020-03-25');
        
        % Find start index
        start_idx = find(data_table.date == start_date, 1);
        if isempty(start_idx)
            % Find closest date
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
        
    catch ME
        warning('Could not load Washington County, MS real data: %s', ME.message);
        % Fallback to hardcoded values
        real_initial_infected = 1;
        real_initial_dead = 0;
        
        i0 = real_initial_infected / N;
        d0 = real_initial_dead / N;
        h0 = 0.0;
        r0 = 0.0;
        s0 = 1.0 - (i0 + h0 + r0 + d0);
    end

    % Define model parameters structure for SIHRS with death model
    params = struct(...
        'beta', 0.172,        ... % beta: infection rate (β > 0) 
        'gamma', 0.141,       ... % gamma: I transition rate (γ > 0) 
        'alpha', 0.1,         ... % alpha: H transition rate (α > 0)
        'lambda', 0.0083,     ... % lambda: R transition rate (Λ > 0) - Updated for Washington, MS
        'pSI', 1.0,           ... % pSI: probability of S to I (p_{SI} in (0,1])
        'pII', 0.00,          ... % pII: probability of I to I (stay infected)
        'pIH', 0.1614,        ... % pIH: probability of I to H - Updated from P(IH) calculation
        'pIR', 0.8367,        ... % pIR: probability of I to R - Updated to sum to 1
        'pID', 0.0019,        ... % pID: probability of I to D - Updated from P(ID) calculation
        'pHH', 0.00,          ... % pHH: probability of H to H (stay hospitalized) - Updated
        'pHR', 0.846,         ... % pHR: probability of H to R - Updated
        'pHD', 0.154,         ... % pHD: probability of H to D - Updated
        'pRR', 0.02,          ... % pRR: probability of R to R (stay recovered)
        'pRS', 0.98,          ... % pRS: probability of R to S
        'tmax', 620,          ... % tmax: simulation end time (matching Carson City)
        's0', s0,             ... % s0: initial susceptible proportion
        'i0', i0,             ... % i0: initial infected proportion
        'h0', h0,             ... % h0: initial hospitalized proportion
        'r0', r0,             ... % r0: initial recovered proportion
        'd0', d0              ... % d0: initial dead proportion
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

    num_simulations = 9;  % Matching Carson City
    
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
            all_results{sim_idx} = sihrs_agent_model_infected(N, params);
        end
        
        fprintf('All simulations completed!\n');
        
        % Plot all simulations together
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
    % For a 620-day simulation, people can get infected multiple times
    % Each person might have ~10 events: S→I→H→R→S→I→R→S→I→D
    max_events = N * 10;  % Estimate maximum number of events
    T = zeros(max_events, 1);
    I_prop = zeros(max_events, 1);
    I_count = zeros(max_events, 1);
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
        D_count(event_count) = length(D);
    end
    
    % Trim unused preallocated space
    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    I_count = I_count(1:event_count);
    D_count = D_count(1:event_count);
    
    % Store results
    result = struct(...
        'N', N, ...
        'T', T, ...
        'I_count', I_count, ...
        'D_count', D_count, ...
        'final_time', t, ...
        'peak_infected', max(I_count), ...
        'peak_time', T(find(I_count == max(I_count), 1)), ...
        'peak_infected_prop', max(I_prop), ...
        'peak_time_prop', T(find(I_prop == max(I_prop), 1)), ...
        'i_inf', I_prop(end) ...
    );
end

function plot_multiple_simulations_infected(all_results, N, params)
    % This function creates multiple figures matching the Julia version:
    % 1. Uncertainty envelope plot for infected cases
    % 2. All stochastic trajectories for infected cases
    % 3. Active deaths plot
    % 4. Cumulative deaths plot

    % --- 1. Create a common time grid for comparison ---
    t_grid = (0:params.tmax)';  % daily time steps
    
    % --- 2. Interpolate each simulation's results onto the grid ---
    all_interp_I = zeros(length(all_results), length(t_grid));
    all_interp_D = zeros(length(all_results), length(t_grid));
    
    for i = 1:length(all_results)
        res = all_results{i};
        % Create interpolation function and vectorize the evaluation
        if length(res.T) > 1
            itp_I = griddedInterpolant(res.T, res.I_count, 'linear', 'none');
            itp_D = griddedInterpolant(res.T, res.D_count, 'linear', 'none');
            
            % Vectorized interpolation - much faster than loop
            for j = 1:length(t_grid)
                if t_grid(j) >= min(res.T) && t_grid(j) <= max(res.T)
                    all_interp_I(i, j) = itp_I(t_grid(j));
                    all_interp_D(i, j) = itp_D(t_grid(j));
                end
            end
        end
    end
    
    % --- Compute active deaths (rolling window) for each simulation ---
    window = 14;
    all_active_D = zeros(size(all_interp_D));
    for i = 1:size(all_interp_D, 1)
        all_active_D(i, :) = compute_rolling_window(all_interp_D(i, :), window);
    end
    
    % --- 3. Calculate statistics for the bandwidth ---
    mean_I = mean(all_interp_I, 1)';
    lower_I = quantile(all_interp_I, 0.05, 1)';  % 5th percentile
    upper_I = quantile(all_interp_I, 0.95, 1)';  % 95th percentile
    mean_D = mean(all_active_D, 1)';
    lower_D = quantile(all_active_D, 0.05, 1)';
    upper_D = quantile(all_active_D, 0.95, 1)';
    
        % --- 4. Load and process real-world data ---
    population = N;  % Use same population as simulation
    real_interp_I = zeros(length(t_grid), 1);
    real_interp_D = zeros(length(t_grid), 1);
    real_interp_D_prop = zeros(length(t_grid), 1);
    real_cumulative_D = zeros(length(t_grid), 1);
    start_date_real = datetime('2020-03-25');  % Default start date

    try
        % Load cases and deaths data
        data_table = readtable('washington_mississippi_combined.csv');
        data_table.date = datetime(data_table.date, 'InputFormat', 'yyyy-MM-dd');
        start_idx = find(data_table.date == start_date_real, 1);
        if isempty(start_idx)
            date_diffs = abs(datenum(data_table.date) - datenum(start_date_real));
            [~, start_idx] = min(date_diffs);
        end

        dates_from_start = data_table.date(start_idx:end);
        cumulative_cases_from_start = data_table.cases(start_idx:end);
        recovery_days = 14;  % Longer recovery period for rural areas (14 days)

        % Calculate active cases (like Julia)
        cumulative_shifted = [zeros(recovery_days, 1); cumulative_cases_from_start(1:end-recovery_days)];
        active_cases_count = cumulative_cases_from_start - cumulative_shifted;
        % Ensure no negative active cases due to data corrections
        active_cases_count = max(active_cases_count, 0);

        carson_days = days(dates_from_start - dates_from_start(1));

        % Interpolate active cases to simulation time grid
        real_interp_I = interp1(carson_days, active_cases_count, t_grid, 'linear', 0);

        % Load daily death data from the dedicated file
        daily_deaths_data = readtable('washington_mississippi_daily_deaths.csv');
        daily_deaths_data.date = datetime(daily_deaths_data.date, 'InputFormat', 'yyyy-MM-dd');

        % Find start index for daily deaths data
        daily_start_idx = find(daily_deaths_data.date == start_date_real, 1);
        if isempty(daily_start_idx)
            date_diffs = abs(datenum(daily_deaths_data.date) - datenum(start_date_real));
            [~, daily_start_idx] = min(date_diffs);
        end

        % Get daily deaths from March 25 onwards
        daily_deaths_from_start = daily_deaths_data.daily_deaths(daily_start_idx:end);
        daily_death_dates = daily_deaths_data.date(daily_start_idx:end);

        % Calculate days from March 25 for daily deaths
        daily_death_days = days(daily_death_dates - start_date_real);

        % Interpolate daily deaths to simulation time grid
        real_interp_D = interp1(daily_death_days, daily_deaths_from_start, t_grid, 'linear', 0);

        % Use the 7-day moving average from the CSV file for active deaths
        moving_avg_deaths = daily_deaths_data.moving_avg_7day(daily_start_idx:end);
        real_active_D = interp1(daily_death_days, moving_avg_deaths, t_grid, 'linear', 0);
        real_interp_D_prop = real_active_D / population;

        % Also get cumulative deaths for the cumulative death plot
        cumulative_deaths_from_start = daily_deaths_data.cumulative_deaths(daily_start_idx:end);
        real_cumulative_D = interp1(daily_death_days, cumulative_deaths_from_start, t_grid, 'linear', 0);

        % Debug: Print data ranges and sizes
        fprintf('Successfully loaded real data with %d data points\n', length(active_cases_count));
        fprintf('Active cases range: %.1f to %.1f\n', min(active_cases_count), max(active_cases_count));
        fprintf('Daily deaths range: %.1f to %.1f\n', min(daily_deaths_from_start), max(daily_deaths_from_start));
        fprintf('Date range: %s to %s\n', datestr(min(dates_from_start)), datestr(max(dates_from_start)));
        fprintf('Daily death days range: %.1f to %.1f\n', min(daily_death_days), max(daily_death_days));
        fprintf('Simulation time grid: %.1f to %.1f\n', min(t_grid), max(t_grid));
        fprintf('Real interp I range: %.6f to %.6f\n', min(real_interp_I), max(real_interp_I));
        fprintf('Real interp D range: %.6f to %.6f\n', min(real_interp_D), max(real_interp_D));
        fprintf('Real active D range: %.6f to %.6f\n', min(real_active_D), max(real_active_D));

    catch ME
        fprintf('Warning: Could not load or process real data: %s\n', ME.message);
        real_interp_I = zeros(length(t_grid), 1);
        real_interp_D_prop = zeros(length(t_grid), 1);
    end

    % Convert all counts to proportions for plotting
    population_factor = 1.0 / N;  % Calculate division factor once
    
    % Convert all arrays to proportions
    all_interp_I_prop = all_interp_I * population_factor;
    mean_I_prop = mean_I * population_factor;
    lower_I_prop = lower_I * population_factor;
    upper_I_prop = upper_I * population_factor;
    all_active_D_prop = all_active_D * population_factor;
    mean_D_prop = mean_D * population_factor;
    lower_D_prop = lower_D * population_factor;
    upper_D_prop = upper_D * population_factor;
    real_interp_I_prop = real_interp_I * population_factor;
    
    % Debug: Check the data before plotting
    fprintf('After conversion to proportions:\n');
    fprintf('Real interp I prop range: %.6f to %.6f\n', min(real_interp_I_prop), max(real_interp_I_prop));
    fprintf('Real interp D prop range: %.6f to %.6f\n', min(real_interp_D_prop), max(real_interp_D_prop));
    fprintf('Real cumulative D prop range: %.6f to %.6f\n', min(real_cumulative_D/population), max(real_cumulative_D/population));
    
    % Check for NaN or Inf values
    fprintf('NaN count in real_interp_I_prop: %d\n', sum(isnan(real_interp_I_prop)));
    fprintf('NaN count in real_interp_D_prop: %d\n', sum(isnan(real_interp_D_prop)));
    fprintf('Zero count in real_interp_I_prop: %d\n', sum(real_interp_I_prop == 0));
    fprintf('Zero count in real_interp_D_prop: %d\n', sum(real_interp_D_prop == 0));

    % Set up date ticks
    tick_interval = 90;
    xtick_positions = 0:tick_interval:params.tmax;
    xtick_dates = start_date_real + days(xtick_positions);
    date_labels = cellstr(datestr(xtick_dates, 'mm/dd/yy'));

    % --- 5. Create the final plot with the uncertainty envelope for infected cases ---
    figure;
    fill([t_grid; flipud(t_grid)], [upper_I_prop; flipud(lower_I_prop)], ...
         [0.7 0.9 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    
    % Debug: Check what we're plotting
    fprintf('Plotting infected cases:\n');
    fprintf('t_grid size: %d, real_interp_I_prop size: %d\n', length(t_grid), length(real_interp_I_prop));
    fprintf('t_grid range: %.1f to %.1f\n', min(t_grid), max(t_grid));
    fprintf('real_interp_I_prop range: %.6f to %.6f\n', min(real_interp_I_prop), max(real_interp_I_prop));
    
    plot(t_grid, real_interp_I_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Time (days)');
    ylabel('Infected Proportion');
    title('Washington, Mississippi');
    xlim([0, params.tmax]);
    ylim([0, max([upper_I_prop; real_interp_I_prop]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Full_Pandemic_bandwidth.png');
    
    % --- 6. Create a second figure with all stochastic sims and the real data ---
    figure;
    % Plot all stochastic simulations as semi-transparent blue lines
    for i = 1:size(all_interp_I_prop, 1)
        plot(t_grid, all_interp_I_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
        hold on;
    end

    % Plot the real data as a solid red line
    plot(t_grid, real_interp_I_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Time (days)');
    ylabel('Infected Proportion');
    title('Washington, Mississippi');
    xlim([0, params.tmax]);
    ylim([0, max([max(max(all_interp_I_prop)), max(real_interp_I_prop)]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    % Custom legend
    plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2.5);
    legend('Stochastic Simulations', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Full_Pandemic_trajectories.png');

    % --- 7. Create a figure for death proportion (active deaths) ---
    figure;
    fill([t_grid; flipud(t_grid)], [upper_D_prop; flipud(lower_D_prop)], ...
         [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    
    % Debug: Check what we're plotting for active deaths
    fprintf('Plotting active deaths:\n');
    fprintf('t_grid size: %d, real_interp_D_prop size: %d\n', length(t_grid), length(real_interp_D_prop));
    fprintf('real_interp_D_prop range: %.6f to %.6f\n', min(real_interp_D_prop), max(real_interp_D_prop));
    fprintf('Non-zero values in real_interp_D_prop: %d\n', sum(real_interp_D_prop > 0));
    
    plot(t_grid, real_interp_D_prop, 'r-', 'LineWidth', 2.5);
    xlabel('Date');
    ylabel('Active Death Proportion');
    title('Washington, Mississippi');
    xlim([0, params.tmax]);
    ylim([0, max([upper_D_prop; real_interp_D_prop]) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('90% Prediction Interval', 'Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Full_Pandemic_active_deaths.png');

    % --- 8. Create a figure for cumulative death proportion ---
    real_interp_D_prop_cumulative = real_cumulative_D / population;
    figure;
    plot(t_grid, real_interp_D_prop_cumulative, 'r-', 'LineWidth', 2.5);
    xlabel('Date');
    ylabel('Cumulative Death Proportion');
    title('Washington, Mississippi');
    xlim([0, params.tmax]);
    ylim([0, max(real_interp_D_prop_cumulative) * 1.1]);

    xticks(xtick_positions);
    xticklabels(date_labels);
    xlabel('Date (mm/dd/yy)');

    legend('Real Data', 'Location', 'best');
    saveas(gcf, 'SIHRS_Washington_MS_Full_Pandemic_cumulative_deaths.png');
end

function result = compute_rolling_window(data, window_size)
    % Helper function for rolling window calculation
    result = data;  % MATLAB automatically copies arrays
    for t = (window_size + 1):length(data)
        result(t) = data(t) - data(t - window_size);
    end
end
