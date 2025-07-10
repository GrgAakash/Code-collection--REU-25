function sihrs_multiple_populations()
    % Define model parameters structure for SIHRS with death model
  params = struct(...
        'beta',1.11, ...    % infection rate (β > 0)
        'gamma', 0.143, ...   % I transition rate (γ > 0)
        'alpha', 0.09, ...   % H transition rate (α > 0)
        'lambda', 0.0083, ...  % R transition rate (Λ > 0) immunity period of 4 months is assumed
        'pSI', 1.0, ...     % probability of S to I (p_{SI} in (0,1])
        'pII', 0.0, ...     % probability of I to I (stay infected)
        'pIH', 0.033, ...    % probability of I to H
        'pIR', 0.966, ...     % probability of I to R
        'pID', 0.001, ...    % probability of I to D
        'pHH', 0.0, ...     % probability of H to H (stay hospitalized)
        'pHR', 0.97, ...     % probability of H to R
        'pHD', 0.03, ...     % probability of H to D
        'pRR', 0, ...     % probability of R to R (stay recovered)
        'pRS', 1, ...     % probability of R to S
        'tmax', 200, ...    % simulation end time
        's0', 0.7, ...      % initial susceptible proportion
        'i0', 0.3, ...      % initial infected proportion
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
        deterministic_result = solve_deterministic_sihrs(params);
        
        % Plot comparison
        plot_comparison(results, N_values, deterministic_result, params);
        
    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
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
    
    % Preallocate arrays for better performance
    max_events = N * 10; % Estimate maximum number of events
    T = zeros(1, max_events);
    S_prop = zeros(1, max_events);
    I_prop = zeros(1, max_events);
    H_prop = zeros(1, max_events);
    R_prop = zeros(1, max_events);
    D_prop = zeros(1, max_events);
    I_count = zeros(1, max_events);
    
    % Initialize agent arrays with preallocation
    S = zeros(1, N);
    S(1:s0) = 1:s0;
    S = S(1:s0);
    
    I = zeros(1, N);
    I(1:i0) = (s0+1):(s0+i0);
    I = I(1:i0);
    
    H = zeros(1, N);
    H(1:h0) = (s0+i0+1):(s0+i0+h0);
    H = H(1:h0);
    
    R = zeros(1, N);
    R(1:r0) = (s0+i0+h0+1):(s0+i0+h0+r0);
    R = R(1:r0);
    
    D = zeros(1, N);
    D(1:d0) = (s0+i0+h0+r0+1):(s0+i0+h0+r0+d0);
    D = D(1:d0);
    
    % Initialize time tracking
    t = 0;
    T(1) = 0;
    event_count = 1;
    
    % Initialize proportion tracking
    total_pop = s0 + i0 + h0 + r0 + d0;
    S_prop(1) = s0 / total_pop;
    I_prop(1) = i0 / total_pop;
    H_prop(1) = h0 / total_pop;
    R_prop(1) = r0 / total_pop;
    D_prop(1) = d0 / total_pop;
    I_count(1) = i0;
    
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
    end
    
    % Trim unused preallocated space
    T = T(1:event_count);
    S_prop = S_prop(1:event_count);
    I_prop = I_prop(1:event_count);
    H_prop = H_prop(1:event_count);
    R_prop = R_prop(1:event_count);
    D_prop = D_prop(1:event_count);
    I_count = I_count(1:event_count);
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.H_prop = H_prop;
    result.R_prop = R_prop;
    result.D_prop = D_prop;
    result.I_count = I_count;
    result.final_time = t;
    result.peak_infected = max(I_count);
    result.peak_time = T(find(I_count == max(I_count), 1, 'first'));
    
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

function plot_comparison(results, N_values, det_result, params)
    % Create comparison plots including deterministic solution
    
    % Create two figures
    figure('Position', [100, 100, 1920, 1440]);  % Original comparison plots
    
    % Use tiledlayout for better spacing control
    t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Colors for different population sizes and deterministic
    colors = {'#0072BD', '#77AC30', '#A2142F'}; % More distinctive colors
    det_color = '#7E2F8E';  % Purple for deterministic
    
    % Plot 1: Susceptible Proportion Over Time
    nexttile;
    hold on;
    plot_handles = [];
    for i = 1:length(results)
        h = plot(results{i}.T, results{i}.S_prop, 'Color', colors{i}, 'LineWidth', 1.5);
        plot_handles = [plot_handles, h];
    end
    h_det = plot(det_result.T, det_result.S_prop, '--', 'Color', det_color, 'LineWidth', 2);
    plot_handles = [plot_handles, h_det];
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Susceptible', 'FontSize', 14);
    title('Susceptible Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    
    % Plot 2: Infected Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.I_prop, 'Color', colors{i}, 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.I_prop, '--', 'Color', det_color, 'LineWidth', 2);
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Infected', 'FontSize', 14);
    title('Infected Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    
    % Plot 3: Hospitalized Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.H_prop, 'Color', colors{i}, 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.H_prop, '--', 'Color', det_color, 'LineWidth', 2);
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Hospitalized', 'FontSize', 14);
    title('Hospitalized Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    
    % Plot 4: Recovered Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.R_prop, 'Color', colors{i}, 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.R_prop, '--', 'Color', det_color, 'LineWidth', 2);
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Recovered', 'FontSize', 14);
    title('Recovered Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    
    % Plot 5: Dead Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.D_prop, 'Color', colors{i}, 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.D_prop, '--', 'Color', det_color, 'LineWidth', 2);
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Dead', 'FontSize', 14);
    title('Dead Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    
    % Add a single legend below all plots
    lgd = legend(plot_handles, [arrayfun(@(x) sprintf('N=%d', x), N_values, 'UniformOutput', false), {'Deterministic'}], ...
        'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12);
    lgd.Layout.Tile = 'south';
    
    % Add parameter display including R0
    param_text = sprintf('R₀=%.2f, β=%.2f, γ=%.2f, α=%.2f, Λ=%.2f\np_{SI}=%.2f, p_{II}=%.2f, p_{IH}=%.2f, p_{IR}=%.2f, p_{ID}=%.2f\np_{HH}=%.2f, p_{HR}=%.2f, p_{HD}=%.2f, p_{RR}=%.2f, p_{RS}=%.2f', ...
        det_result.R0, params.beta, params.gamma, params.alpha, params.lambda, ...
        params.pSI, params.pII, params.pIH, params.pIR, params.pID, ...
        params.pHH, params.pHR, params.pHD, params.pRR, params.pRS);
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', param_text, ...
        'FontSize', 14, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    
    % Save the figure
    saveas(gcf, 'SIHRS_simulation_results.png');
    
    % Create new figure for deterministic curves only
    figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Plot all deterministic curves together
    plot(det_result.T, det_result.S_prop, 'LineWidth', 2, 'DisplayName', 'Susceptible');
    plot(det_result.T, det_result.I_prop, 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(det_result.T, det_result.H_prop, 'LineWidth', 2, 'DisplayName', 'Hospitalized');
    plot(det_result.T, det_result.R_prop, 'LineWidth', 2, 'DisplayName', 'Recovered');
    plot(det_result.T, det_result.D_prop, 'LineWidth', 2, 'DisplayName', 'Dead');
    
    % Customize the plot
    xlabel('Time', 'FontSize', 14);
    ylabel('Population Proportion', 'FontSize', 14);
    title('Deterministic SIHRS with Death Model Dynamics', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    legend('Location', 'east', 'FontSize', 12);
    
    % Add R0 and parameters text
    text(0.02, -0.15, sprintf('R₀=%.2f, β=%.2f, γ=%.2f, α=%.2f, Λ=%.2f', det_result.R0, params.beta, params.gamma, params.alpha, params.lambda), ...
        'Units', 'normalized', 'FontSize', 12);
    
    % Save the deterministic plot
    saveas(gcf, 'SIHRS_deterministic_only.png');
    
    % Create separate figures for each population size N
    for i = 1:length(results)
        % Create figure for this N value
        figure('Position', [100, 100, 1000, 700]);
        
        % Plot all compartments for this N value
        plot(results{i}.T, results{i}.S_prop, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible (S)');
        hold on;
        plot(results{i}.T, results{i}.I_prop, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected (I)');
        plot(results{i}.T, results{i}.H_prop, 'm-', 'LineWidth', 2, 'DisplayName', 'Hospitalized (H)');
        plot(results{i}.T, results{i}.R_prop, 'g-', 'LineWidth', 2, 'DisplayName', 'Recovered (R)');
        plot(results{i}.T, results{i}.D_prop, 'k-', 'LineWidth', 2, 'DisplayName', 'Dead (D)');
        
        % Customize the plot
        xlabel('Time (days)', 'FontSize', 12);
        ylabel('Proportion', 'FontSize', 12);
        title(sprintf('SIHRS Stochastic Model - N = %d', N_values(i)), 'FontSize', 14);
        legend('Location', 'eastoutside', 'FontSize', 10);
        grid on;
        xlim([0, params.tmax]);
        ylim([0, 1]);
        
        % Add parameter annotations
        param_text = sprintf('R₀=%.2f, β=%.4f, γ=%.4f, α=%.4f, λ=%.4f', ...
            det_result.R0, params.beta, params.gamma, params.alpha, params.lambda);
        text(0.02, 0.98, param_text, 'Units', 'normalized', 'FontSize', 10, ...
            'VerticalAlignment', 'top', 'BackgroundColor', 'white');
        
        % Save the figure
        saveas(gcf, sprintf('SIHRS_stochastic_N%d.png', N_values(i)));
    end
    
    % Print summary statistics
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('Population Size | Peak Infected | Peak Time | s(∞) | i(∞) | h(∞) | r(∞) | d(∞)\n');
    fprintf('----------------|---------------|-----------|-------|-------|-------|-------|-------\n');
    for i = 1:length(results)
        fprintf('%15d | %13d | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
            results{i}.N, results{i}.peak_infected, results{i}.peak_time, ...
            results{i}.s_inf, results{i}.i_inf, results{i}.h_inf, results{i}.r_inf, results{i}.d_inf);
    end
    fprintf('%15s | %13.4f | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
        'Deterministic', det_result.peak_infected_prop, det_result.peak_time, ...
        det_result.s_inf, det_result.i_inf, det_result.h_inf, det_result.r_inf, det_result.d_inf);
    
    % Print asymptotic analysis
    fprintf('\n=== ASYMPTOTIC ANALYSIS ===\n');
    fprintf('R₀ = %.4f\n', det_result.R0);
    fprintf('s(∞) = %.6f\n', det_result.s_inf);
    fprintf('i(∞) = %.6f\n', det_result.i_inf);
    fprintf('h(∞) = %.6f\n', det_result.h_inf);
    fprintf('r(∞) = %.6f\n', det_result.r_inf);
    fprintf('d(∞) = %.6f\n', det_result.d_inf);
end

% Run the simulation 