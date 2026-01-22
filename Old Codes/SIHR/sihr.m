clear all;
close all;
function sir_multiple_populations()
    % Define model parameters structure
    params = struct(...
        'beta', 1.30, ...    % infection rate (β > 0)
        'gamma1', 1.0, ...  % I to H rate (γ₁ > 0)
        'gamma2', 1.0, ...  % I to R rate (γ₂ > 0)
        'alpha', 1.00, ...   % H to R rate (α > 0)
        'p1', 0.5, ...      % probability of infection (p₁ in (0,1])
        'p2', 0.5, ...      % probability of leaving I (p₂ in (0,1])
        'p3', 0.5, ...      % probability of leaving H (p₃ in (0,1])
        'ph', 0.50, ...      % probability of I to H vs R (p_h in (0,1])
        'tmax', 50, ...     % simulation end time
        's0', 0.96, ...     % initial susceptible proportion
        'i0', 0.04, ...     % initial infected proportion
        'h0', 0.0, ...     % initial hospitalized proportion
        'r0', 0.00 ...      % initial recovered proportion
    );

    % Validate parameters
    validateParameters(params);
    
    % Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.r0) - 1) > 1e-10
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
            results{idx} = sir_agent_model(N_values(idx), params);
            fprintf('Completed N = %d\n', N_values(idx));
        end
        
        % Solve deterministic model
        deterministic_result = solve_deterministic_sir(params);
        
        % Plot comparison
        plot_comparison(results, N_values, deterministic_result, params);
        
        % Plot ODE only
        plot_ode_only(deterministic_result, params);
        
    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
end

function validateParameters(params)
    % Validate rates are positive
    if any([params.beta, params.gamma1, params.gamma2, params.alpha] <= 0)
        error('All rates (beta, gamma1, gamma2, alpha) must be positive');
    end
    
    % Validate probabilities are in (0,1]
    probs = [params.p1, params.p2, params.p3, params.ph];
    if any(probs <= 0 | probs > 1)
        error('All probabilities must be in (0,1]');
    end
end

function result = sir_agent_model(N, params)
    % SIHR agent-based stochastic model
    validateattributes(N, {'numeric'}, {'positive', 'integer', 'scalar'});
    
    % Initial conditions - using params values and ensuring they sum to N
    s0 = round(params.s0 * N); % susceptible
    i0 = round(params.i0 * N); % infected
    h0 = round(params.h0 * N); % hospitalized
    r0 = round(params.r0 * N); % recovered
    
    % Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + r0;
    if total ~= N
        % Add or subtract the difference from the largest compartment
        [~, largest_idx] = max([s0, i0, h0, r0]);
        switch largest_idx
            case 1
                s0 = s0 + (N - total);
            case 2
                i0 = i0 + (N - total);
            case 3
                h0 = h0 + (N - total);
            case 4
                r0 = r0 + (N - total);
        end
    end
    
    % Validate initial conditions sum to N
    if (s0 + i0 + h0 + r0) ~= N
        error('Initial conditions must sum to N');
    end
    
    % Preallocate arrays for better performance
    max_events = N * 10; % Estimate maximum number of events
    T = zeros(1, max_events);
    S_prop = zeros(1, max_events);
    I_prop = zeros(1, max_events);
    H_prop = zeros(1, max_events);
    R_prop = zeros(1, max_events);
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
    if r0 > 0
        R(1:r0) = (s0+i0+h0+1):(s0+i0+h0+r0);
        R = R(1:r0);
    else
        R = zeros(1, 0);
    end
    
    % Initialize time tracking
    t = 0;
    T(1) = 0;
    event_count = 1;
    
    % Initialize proportion tracking
    total_pop = s0 + i0 + h0 + r0;
    S_prop(1) = s0 / total_pop;
    I_prop(1) = i0 / total_pop;
    H_prop(1) = h0 / total_pop;
    R_prop(1) = r0 / total_pop;
    I_count(1) = i0;
    
    % Main simulation loop
    while ~isempty(I) && t < params.tmax
        nI = numel(I);
        nS = numel(S);
        nH = numel(H);
        
        % Calculate event rates according to the mathematical model
        infection_rate = params.p1 * params.beta * nS * nI / N;  % S to I rate
        to_hospital_rate = params.p2 * params.ph * params.gamma1 * nI;  % I to H rate
        to_recovered_from_I_rate = params.p2 * (1-params.ph) * params.gamma2 * nI;  % I to R rate
        to_recovered_from_H_rate = params.p3 * params.alpha * nH;  % H to R rate
        
        total_rate = infection_rate + to_hospital_rate + to_recovered_from_I_rate + to_recovered_from_H_rate;
        
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
            current_total = numel(S) + numel(I) + numel(H) + numel(R);
            S_prop(event_count) = numel(S) / current_total;
            I_prop(event_count) = numel(I) / current_total;
            H_prop(event_count) = numel(H) / current_total;
            R_prop(event_count) = numel(R) / current_total;
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
        elseif chance < (infection_rate + to_hospital_rate)
            % I to H transition
            if nI > 0
                num = randi([1, nI]);
                hospitalized_agent = I(num);
                I(num) = [];
                H(end+1) = hospitalized_agent;
            end
        elseif chance < (infection_rate + to_hospital_rate + to_recovered_from_I_rate)
            % I to R transition
            if nI > 0
                num = randi([1, nI]);
                recovered_agent = I(num);
                I(num) = [];
                R(end+1) = recovered_agent;
            end
        else
            % H to R transition
            if nH > 0
                num = randi([1, nH]);
                recovered_agent = H(num);
                H(num) = [];
                R(end+1) = recovered_agent;
            end
        end
        
        % Update tracking arrays
        current_total = numel(S) + numel(I) + numel(H) + numel(R);
        S_prop(event_count) = numel(S) / current_total;
        I_prop(event_count) = numel(I) / current_total;
        H_prop(event_count) = numel(H) / current_total;
        R_prop(event_count) = numel(R) / current_total;
        I_count(event_count) = numel(I);
    end
    
    % Trim unused preallocated space
    T = T(1:event_count);
    S_prop = S_prop(1:event_count);
    I_prop = I_prop(1:event_count);
    H_prop = H_prop(1:event_count);
    R_prop = R_prop(1:event_count);
    I_count = I_count(1:event_count);
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.H_prop = H_prop;
    result.R_prop = R_prop;
    result.I_count = I_count;
    result.final_time = t;
    result.peak_infected = max(I_count);
    result.peak_time = T(find(I_count == max(I_count), 1, 'first'));
    
    % Calculate and store asymptotic values
    result.s_inf = S_prop(end);
    result.i_inf = I_prop(end);
    result.h_inf = H_prop(end);
    result.r_inf = R_prop(end);
end

function det_result = solve_deterministic_sir(params)
    % Solve the deterministic SIHR model using ODE45
    
    % Time span
    tspan = [0, params.tmax];
    
    % Initial conditions vector
    y0 = [params.s0; params.i0; params.h0; params.r0];
    
    % Define the ODE system exactly as in the mathematical model
    ode_system = @(t, y) [
        -params.p1 * params.beta * y(1) * y(2);                                                          % ds/dt
        params.p1 * params.beta * y(1) * y(2) - params.p2 * (params.ph * params.gamma1 + (1-params.ph) * params.gamma2) * y(2);  % di/dt
        params.p2 * params.ph * params.gamma1 * y(2) - params.p3 * params.alpha * y(3);                  % dh/dt
        params.p2 * (1-params.ph) * params.gamma2 * y(2) + params.p3 * params.alpha * y(3)              % dr/dt
    ];
    
    % Set ODE options for better accuracy
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    % Solve the ODE system
    [T, Y] = ode45(ode_system, tspan, y0, options);
    
    % Store results
    det_result.T = T;
    det_result.S_prop = Y(:, 1);
    det_result.I_prop = Y(:, 2);
    det_result.H_prop = Y(:, 3);
    det_result.R_prop = Y(:, 4);
    
    % Find peak infected and peak time
    [peak_infected_prop, peak_idx] = max(det_result.I_prop);
    det_result.peak_infected_prop = peak_infected_prop;
    det_result.peak_time = T(peak_idx);
    det_result.final_time = T(end);
    
    % Calculate R0
    det_result.R0 = params.p1 * params.beta / (params.p2 * (params.ph * params.gamma1 + (1-params.ph) * params.gamma2));
    
    % Store asymptotic values
    det_result.s_inf = det_result.S_prop(end);
    det_result.i_inf = det_result.I_prop(end);
    det_result.h_inf = det_result.H_prop(end);
    det_result.r_inf = det_result.R_prop(end);
    
    % Verify asymptotic behavior
    if det_result.i_inf > 1e-6
        warning('i(∞) may not be approaching 0 as expected');
    end
    if det_result.h_inf > 1e-6
        warning('h(∞) may not be approaching 0 as expected');
    end
    
    % Verify the asymptotic relation from the theorem
    s_inf_theoretical = params.s0 * exp(-det_result.R0 * (params.s0 + params.i0 - det_result.s_inf));
    if abs(det_result.s_inf - s_inf_theoretical) > 1e-6
        warning('Asymptotic s(∞) relation may not be satisfied');
    end
end

function plot_comparison(results, N_values, det_result, params)
    % Create comparison plots including deterministic solution
    
    % Main figure - increased size
    fig = figure('Position', [100, 100, 1920, 1440]);
    
    % Use tiledlayout for better spacing control
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
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
    
    % Add a single legend below all plots
    lgd = legend(plot_handles, [arrayfun(@(x) sprintf('N=%d', x), N_values, 'UniformOutput', false), {'Deterministic'}], ...
        'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12);
    lgd.Layout.Tile = 'south';
    
    % Add parameter display including R0
    param_text = sprintf('R₀=%.2f, β=%.2f, γ₁=%.2f, γ₂=%.2f, α=%.2f\np₁=%.2f, p₂=%.2f, p₃=%.2f, p_h=%.2f', ...
        det_result.R0, params.beta, params.gamma1, params.gamma2, params.alpha, ...
        params.p1, params.p2, params.p3, params.ph);
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', param_text, ...
        'FontSize', 14, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    
    % Save the figure
    saveas(fig, 'SIHR_simulation_results.png');
    
    % Print summary statistics
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('Population Size | Peak Infected | Peak Time | s(∞) | i(∞) | h(∞) | r(∞)\n');
    fprintf('----------------|---------------|-----------|-------|-------|-------|-------\n');
    for i = 1:length(results)
        fprintf('%15d | %13d | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
            results{i}.N, results{i}.peak_infected, results{i}.peak_time, ...
            results{i}.s_inf, results{i}.i_inf, results{i}.h_inf, results{i}.r_inf);
    end
    fprintf('%15s | %13.4f | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
        'Deterministic', det_result.peak_infected_prop, det_result.peak_time, ...
        det_result.s_inf, det_result.i_inf, det_result.h_inf, det_result.r_inf);
    
    % Print asymptotic analysis
    fprintf('\n=== ASYMPTOTIC ANALYSIS ===\n');
    fprintf('R₀ = %.4f\n', det_result.R0);
    fprintf('Theoretical s(∞) = %.6f\n', params.s0 * exp(-det_result.R0 * (params.s0 + params.i0 - det_result.s_inf)));
    fprintf('Numerical s(∞) = %.6f\n', det_result.s_inf);
    fprintf('i(∞) = %.6f (should be ≈ 0)\n', det_result.i_inf);
    fprintf('h(∞) = %.6f (should be = 0)\n', det_result.h_inf);
    fprintf('r(∞) = %.6f\n', det_result.r_inf);
end

function plot_ode_only(det_result, params)
    % Create a separate plot showing only the ODE solution for SIHR
    
    % Create figure for ODE only
    fig_ode = figure('Position', [100, 100, 1000, 700]);
    
    % Color scheme for ODE plots
    ode_colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; % Blue, Orange, Yellow, Purple
    
    % Plot all four compartments on a single chart
    hold on;
    plot(det_result.T, det_result.S_prop, 'Color', ode_colors{1}, 'LineWidth', 2.5, 'DisplayName', 'Susceptible (S)');
    plot(det_result.T, det_result.I_prop, 'Color', ode_colors{2}, 'LineWidth', 2.5, 'DisplayName', 'Infected (I)');
    plot(det_result.T, det_result.H_prop, 'Color', ode_colors{3}, 'LineWidth', 2.5, 'DisplayName', 'Hospitalized (H)');
    plot(det_result.T, det_result.R_prop, 'Color', ode_colors{4}, 'LineWidth', 2.5, 'DisplayName', 'Recovered (R)');
    hold off;
    
    % Customize the plot
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion of Population', 'FontSize', 14);
    title('SIHR Model - Deterministic ODE Solution', 'FontSize', 16, 'FontWeight', 'bold');
    grid on;
    xlim([0, params.tmax]);
    ylim([0, 1]);
    
    % Add legend
    legend('Location', 'eastoutside', 'FontSize', 12);
    
    % Save the figure
    saveas(fig_ode, 'SIHR_ODE_only.png');
    
    % Print ODE-specific summary
    fprintf('\n=== ODE SOLUTION SUMMARY ===\n');
    fprintf('Peak Infected Proportion: %.4f at time %.2f\n', det_result.peak_infected_prop, det_result.peak_time);
    fprintf('R₀: %.4f\n', det_result.R0);
    fprintf('Asymptotic Values:\n');
    fprintf('  s(∞) = %.6f\n', det_result.s_inf);
    fprintf('  i(∞) = %.6f\n', det_result.i_inf);
    fprintf('  h(∞) = %.6f\n', det_result.h_inf);
    fprintf('  r(∞) = %.6f\n', det_result.r_inf);
end

% Run the simulation
sir_multiple_populations();