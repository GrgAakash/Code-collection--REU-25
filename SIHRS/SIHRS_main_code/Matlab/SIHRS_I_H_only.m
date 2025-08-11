function sihrs_I_H_analysis()
    % SIHRS Epidemic Model - I and H Compartments Only (MATLAB Version)
    %
    % This script focuses specifically on plotting the Infected (I) and 
    % Hospitalized (H) compartments from the SIHRS epidemic model with death.
    
    % Define model parameters structure for SIHRS with death model
    params = struct(...
        'beta', 0.212, ...    % infection rate (β > 0)
        'gamma', 0.10, ...   % I transition rate (γ > 0)
        'alpha', 0.09, ...   % H transition rate (α > 0)
        'lambda', 0.01, ...   % R transition rate (λ > 0) immunity period
        'pSI', 1.0, ...     % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0, ...     % probability of I to I (stay infected)
        'pIH', 0.1, ...      % probability of I to H
        'pIR', 0.88, ...     % probability of I to R
        'pID', 0.02, ...    % probability of I to D
        'pHH', 0.0, ...     % probability of H to H (stay hospitalized)
        'pHR', 0.92, ...     % probability of H to R
        'pHD', 0.08, ...     % probability of H to D
        'pRR', 0.02, ...     % probability of R to R (stay recovered)
        'pRS', 0.98, ...     % probability of R to S
        'tmax', 1000, ...    % simulation end time
        's0', 0.96, ...      % initial susceptible proportion
        'i0', 0.04, ...      % initial infected proportion
        'h0', 0.0, ...      % initial hospitalized proportion
        'r0', 0.0, ...      % initial recovered proportion
        'd0', 0.0 ...       % initial dead proportion
    ); 
    
    % Validate parameters
    validateParameters(params);
    
    % Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1) > 1e-10
        error('Initial conditions must sum to 1');
    end
    
    % Population sizes to test
    N_values = [316, 3162, 10000];
    
    % Input validation
    if any(N_values <= 0)
        error('Population sizes must be positive integers');
    end
    
    % Store results for comparison
    results = cell(length(N_values), 1);
    
    % Run simulation for each population size
    try
        for idx = 1:length(N_values)
            fprintf('Running simulation for N = %d...\n', N_values(idx));
            results{idx} = sihrs_agent_model(N_values(idx), params);
            fprintf('Completed N = %d\n', N_values(idx));
        end
        
        % Solve deterministic model
        fprintf('Solving deterministic model...\n');
        deterministic_result = solve_deterministic_sihrs(params);
        fprintf('Deterministic model completed\n');
        
        % Plot I and H focused analysis
        fprintf('Generating I and H focused plots...\n');
        plot_I_H_compartments(results, N_values, deterministic_result, params);
        fprintf('Plots generated and saved\n');
        
    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
    
    % Print focused summary statistics
    fprintf('\n=== I AND H COMPARTMENT ANALYSIS ===\n');
    fprintf('Population Size | Peak I | Peak H | Peak I Time | Peak H Time | I(∞) | H(∞)\n');
    fprintf('----------------|---------|---------|-------------|-------------|-------|-------\n');
    for i = 1:length(results)
        fprintf('%15d | %7d | %7d | %11.2f | %11.2f | %5.3f | %5.3f\n', ...
            results{i}.N, results{i}.peak_infected, results{i}.peak_hospitalized, ...
            results{i}.peak_time, results{i}.peak_h_time, ...
            results{i}.i_inf, results{i}.h_inf);
    end
    
    fprintf('%15s | %7.4f | %7.4f | %11.2f | %11.2f | %5.3f | %5.3f\n', ...
        'Deterministic', deterministic_result.peak_infected_prop, ...
        deterministic_result.peak_hospitalized_prop, ...
        deterministic_result.peak_time, deterministic_result.peak_h_time, ...
        deterministic_result.i_inf, deterministic_result.h_inf);
    
    % Print key metrics
    fprintf('\n=== KEY METRICS ===\n');
    fprintf('R₀ = %.4f\n', deterministic_result.R0);
    fprintf('Peak I proportion = %.4f\n', deterministic_result.peak_infected_prop);
    fprintf('Peak H proportion = %.4f\n', deterministic_result.peak_hospitalized_prop);
    fprintf('I peak time = %.2f days\n', deterministic_result.peak_time);
    fprintf('H peak time = %.2f days\n', deterministic_result.peak_h_time);
    
    % Return results for further analysis
    varargout{1} = results;
    varargout{2} = deterministic_result;
end

function validateParameters(params)
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
    % SIHRS agent-based stochastic model with death
    validateattributes(N, {'numeric'}, {'positive', 'integer', 'scalar'});
    
    % Initial conditions - using params values and ensuring they sum to N
    s0 = round(params.s0 * N); % susceptible
    i0 = round(params.i0 * N); % infected
    h0 = round(params.h0 * N); % hospitalized
    r0 = round(params.r0 * N); % recovered
    d0 = round(params.d0 * N); % dead
    
    % Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0 + d0;
    if total ~= N
        % Add or subtract the difference from the largest compartment
        [~, largest_idx] = max([s0, i0, h0, r0, d0]);
        switch largest_idx
            case 1
                s0 = s0 + (N - total);
            case 2
                i0 = i0 + (N - total);
            case 3
                h0 = h0 + (N - total);
            case 4
                r0 = r0 + (N - total);
            case 5
                d0 = d0 + (N - total);
        end
    end
    
    % Validate initial conditions sum to N
    if (s0 + i0 + h0 + r0 + d0) ~= N
        error('Initial conditions must sum to N');
    end
    
    % Initialize agent arrays
    S = 1:s0;           % susceptible agents
    I = (s0+1):(s0+i0); % infected agents
    H = (s0+i0+1):(s0+i0+h0); % hospitalized agents
    R = (s0+i0+h0+1):(s0+i0+h0+r0); % recovered agents
    D = (s0+i0+h0+r0+1):(s0+i0+h0+r0+d0); % dead agents
    
    % Initialize tracking arrays
    T = 0;           % time points
    S_prop = params.s0; % susceptible proportions
    I_prop = params.i0; % infected proportions
    H_prop = params.h0; % hospitalized proportions
    R_prop = params.r0; % recovered proportions
    D_prop = params.d0; % dead proportions
    I_count = i0;           % infected counts
    H_count = h0;           % hospitalized counts
    
    % Initialize time tracking
    t = 0;
    event_count = 1;
    
    % Main simulation loop
    while ~isempty(I) && t < params.tmax
        nS = numel(S);
        nI = numel(I);
        nH = numel(H);
        nR = numel(R);
        
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
            current_total = numel(S) + numel(I) + numel(H) + numel(R) + numel(D);
            S_prop(event_count) = numel(S) / current_total;
            I_prop(event_count) = numel(I) / current_total;
            H_prop(event_count) = numel(H) / current_total;
            R_prop(event_count) = numel(R) / current_total;
            D_prop(event_count) = numel(D) / current_total;
            I_count(event_count) = numel(I);
            H_count(event_count) = numel(H);
            break;
        end
        
        event_count = event_count + 1;
        T(event_count) = t;
        
        % Determine which event occurs
        chance = rand * total_rate;
        if chance < infection_rate
            % S to I transition
            if nS > 0
                num = randi([1, nS]);
                infected_agent = S(num);
                S(num) = [];
                I(end+1) = infected_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate)
            % R to S transition
            if nR > 0
                num = randi([1, nR]);
                susceptible_agent = R(num);
                R(num) = [];
                S(end+1) = susceptible_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate)
            % I to H transition
            if nI > 0
                num = randi([1, nI]);
                hospitalized_agent = I(num);
                I(num) = [];
                H(end+1) = hospitalized_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate)
            % I to R transition
            if nI > 0
                num = randi([1, nI]);
                recovered_agent = I(num);
                I(num) = [];
                R(end+1) = recovered_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate)
            % I to D transition
            if nI > 0
                num = randi([1, nI]);
                dead_agent = I(num);
                I(num) = [];
                D(end+1) = dead_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate)
            % H to R transition
            if nH > 0
                num = randi([1, nH]);
                recovered_agent = H(num);
                H(num) = [];
                R(end+1) = recovered_agent;
            end
        else
            % H to D transition
            if nH > 0
                num = randi([1, nH]);
                dead_agent = H(num);
                H(num) = [];
                D(end+1) = dead_agent;
            end
        end
        
        % Update tracking arrays
        current_total = numel(S) + numel(I) + numel(H) + numel(R) + numel(D);
        S_prop(event_count) = numel(S) / current_total;
        I_prop(event_count) = numel(I) / current_total;
        H_prop(event_count) = numel(H) / current_total;
        R_prop(event_count) = numel(R) / current_total;
        D_prop(event_count) = numel(D) / current_total;
        I_count(event_count) = numel(I);
        H_count(event_count) = numel(H);
    end
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.H_prop = H_prop;
    result.R_prop = R_prop;
    result.D_prop = D_prop;
    result.I_count = I_count;
    result.H_count = H_count;
    result.final_time = t;
    result.peak_infected = max(I_count);
    result.peak_hospitalized = max(H_count);
    result.peak_time = T(find(I_count == max(I_count), 1, 'first'));
    result.peak_h_time = T(find(H_count == max(H_count), 1, 'first'));
    
    % Calculate and store asymptotic values
    result.s_inf = S_prop(end);
    result.i_inf = I_prop(end);
    result.h_inf = H_prop(end);
    result.r_inf = R_prop(end);
    result.d_inf = D_prop(end);
end

function det_result = solve_deterministic_sihrs(params)
    % Solve the deterministic SIHRS model using ODE45
    
    % Time span
    tspan = [0, params.tmax];
    
    % Initial conditions vector
    y0 = [params.s0; params.i0; params.h0; params.r0; params.d0];
    
    % Define the ODE system exactly as in the mathematical model
    ode_system = @(t, y) [
        -params.beta * y(1) * y(2) * params.pSI + params.pRS * params.lambda * y(4); % ds/dt
        params.beta * y(1) * y(2) * params.pSI - params.gamma * (1 - params.pII) * y(2); % di/dt
        params.pIH * params.gamma * y(2) - params.alpha * (1 - params.pHH) * y(3); % dh/dt
        params.pIR * params.gamma * y(2) + params.pHR * params.alpha * y(3) - params.pRS * params.lambda * y(4); % dr/dt
        params.pID * params.gamma * y(2) + params.pHD * params.alpha * y(3) % dd/dt
    ];
    
    % Set ODE options for better accuracy
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    % Solve the ODE system
    [T, Y] = ode45(ode_system, tspan, y0, options);
    
    % Verify conservation
    sum_y = sum(Y, 2);
    if any(abs(sum_y - 1) > 1e-6)
        warning('Conservation of population not satisfied');
    end
    
    % Store results
    det_result.T = T;
    det_result.S_prop = Y(:, 1);
    det_result.I_prop = Y(:, 2);
    det_result.H_prop = Y(:, 3);
    det_result.R_prop = Y(:, 4);
    det_result.D_prop = Y(:, 5);
    
    % Find peak infected and peak time
    [peak_infected_prop, peak_idx] = max(det_result.I_prop);
    det_result.peak_infected_prop = peak_infected_prop;
    det_result.peak_time = T(peak_idx);
    
    % Find peak hospitalized and peak time
    [peak_hospitalized_prop, peak_h_idx] = max(det_result.H_prop);
    det_result.peak_hospitalized_prop = peak_hospitalized_prop;
    det_result.peak_h_time = T(peak_h_idx);
    
    det_result.final_time = T(end);
    
    % Calculate R0
    det_result.R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    
    % Store asymptotic values
    det_result.s_inf = det_result.S_prop(end);
    det_result.i_inf = det_result.I_prop(end);
    det_result.h_inf = det_result.H_prop(end);
    det_result.r_inf = det_result.R_prop(end);
    det_result.d_inf = det_result.D_prop(end);
end

function plot_I_H_compartments(results, N_values, det_result, params)
    % Create focused plots for I and H compartments only
    
    % 1. Create combined I and H ODE plot (deterministic only)
    figure('Position', [100, 100, 1000, 700]);
    plot(det_result.T, det_result.I_prop, 'r-', 'LineWidth', 3, 'DisplayName', 'Infected (I)');
    hold on;
    plot(det_result.T, det_result.H_prop, 'Color', [0, 1, 0], 'LineWidth', 3, 'DisplayName', 'Hospitalized (H)');
    xlabel('Time');
    ylabel('Proportion');
    title('Deterministic I and H Compartments');
    grid on;
    xlim([0, params.tmax]);
    legend('Location', 'northeast');
    saveas(gcf, 'SIHRS_ODE_I_H_combined.png');
    
    % 2. Create individual plots for each population size N
    for i = 1:length(results)
        % Create figure for this N value with I and H separately
        figure('Position', [100, 100, 1000, 700]);
        
        % Plot I compartment for this N
        plot(results{i}.T, results{i}.I_prop, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected (I)');
        hold on;
        
        % Plot H compartment for this N  
        plot(results{i}.T, results{i}.H_prop, 'Color', [0, 1, 0], 'LineWidth', 2, 'DisplayName', 'Hospitalized (H)');
        
        % Customize the plot
        xlabel('Time (days)');
        ylabel('Proportion');
        title(sprintf('SIHRS Stochastic Model - N = %d', N_values(i)));
        legend('Location', 'northeast');
        grid on;
        xlim([0, params.tmax]);
        
        % Add parameter annotations
        param_text = sprintf('R₀=%.2f, β=%.4f, γ=%.4f, α=%.4f, λ=%.4f', ...
            det_result.R0, params.beta, params.gamma, params.alpha, params.lambda);
        text(0.02, 0.98, param_text, 'Units', 'normalized', 'FontSize', 10, ...
            'VerticalAlignment', 'top', 'BackgroundColor', 'white');
        
        % Save the figure
        saveas(gcf, sprintf('SIHRS_stochastic_N%d_I_H.png', N_values(i)));
    end
end

% Run the analysis if this file is executed directly
if ~exist('OCTAVE_VERSION', 'builtin')
    % This is MATLAB, not Octave
    sihrs_I_H_analysis();
end
