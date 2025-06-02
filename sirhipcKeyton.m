clear all;
close all;

function sir_multiple_populations()
    % Population sizes to test
    N_values = [316, 3162, 10000];
    
    % Store results for comparison
    results = cell(length(N_values), 1);
    
    % Run simulation for each population size
    for idx = 1:length(N_values)
        fprintf('Running simulation for N = %d...\n', N_values(idx));
        results{idx} = sir_agent_model(N_values(idx));
        fprintf('Completed N = %d\n', N_values(idx));
    end
    
    % Solve deterministic model
    deterministic_result = solve_deterministic_sir();
    
    % Plot comparison
    plot_comparison(results, N_values, deterministic_result);
end

function result = sir_agent_model(N)
    % SIR agent-based stochastic model
    
    % Initial conditions (same proportions as original)
    s0 = round(0.96 * N); % susceptible agents at time 0
    i0 = round(0.04 * N); % infected at time 0
    h0 = 0; % hospitalized at time 0
    r0 = 0;               % recovered at time 0
    
    % Model parameters
    beta = 0.95;    % infection constant
    gamma = 1;      % infectious to hospitalized or recovered constant
    alpha = 1; %hospitalized to recovered constant
    p1 = 0.5;       % odds of infection
    p2 = 0.5;       % odds of moving from infectious to hospitalied or recovered
    p3 = 0.5; % odds of moving from hospitalized to recoverd
    ph = 0.5; % odds of moving from infectious to hospitalized

    
    % Initialize agent arrays
    S = 1:s0;       % array of susceptible agent indices
    I = (s0+1):(s0+i0); % array of infected agent indices
    H = []; %array of hospitalized agent indices (initially empty)
    R = [];         % array of recovered agent indices (initially empty)
    
    % Initialize time tracking
    t = 0;
    T = [0];
    
    % Initialize proportion tracking
    total_pop = s0 + i0 + r0 + h0;
    S_prop = [s0 / total_pop];
    I_prop = [i0 / total_pop];
    R_prop = [r0 / total_pop];
    H_prop = [h0 / total_pop];
    I_count = [i0];
    
    % Main simulation loop
    while ~isempty(I) && t < 25
        nI = numel(I);
        nS = numel(S);
        nH = numel(H);
        
        % Calculate event rates
        infection_rate = beta * nS * nI / N;
        recovery_rate = gamma * nI;
        hospital_rate = alpha * nH;
        event_rate = infection_rate + recovery_rate + hospital_rate;
        
        if event_rate == 0
            break;
        end
        
        % Time of next event (exponential distribution)
        dt = exprnd(1 / event_rate);
        t = t + dt;
        if t > 25
            t = 25; % Cap time at 25
            T(end+1) = t;
            current_total = numel(S) + numel(I) + numel(H) + numel(R);
            S_prop(end+1) = numel(S) / current_total;
            I_prop(end+1) = numel(I) / current_total;
            H_prop(end+1) = numel(H) / current_total;
            R_prop(end+1) = numel(R) / current_total;
            I_count(end+1) = numel(I);
            break;
        end
        T(end+1) = t;
        % Determine which event occurs
        chance = rand;
        if chance < (infection_rate / event_rate)
            % Infection event
            if nS > 0 && rand < p1
                num = randi([1, nS]);
                infected_agent = S(num);
                S(num) = [];
                I(end+1) = infected_agent;
            end
        elseif chance < (infection_rate / event_rate + recovery_rate / event_rate)
            % Recovery event
            if nI > 0 && rand < p2
                if rand < ph
                    num = randi([1,nI]);
                    hospitalized_agent = I(num);
                    I(num) = [];
                    H(end+1) = hospitalized_agent;
                else
                    num = randi([1, nI]);
                    recovered_agent = I(num);
                    I(num) = [];
                    R(end+1) = recovered_agent;
                end
            end
          
        else
            % Hospital event
            if nH > 0
                if rand < p3
                    num = randi([1,nH]);
                    recovered_agent = H(num);
                    H(num) = [];
                    R(end+1) = recovered_agent;
                end
            end
        end
        
        % Update tracking arrays
        current_total = numel(S) + numel(I) + numel(R) +numel(H);
        S_prop(end+1) = numel(S) / current_total;
        I_prop(end+1) = numel(I) / current_total;
        R_prop(end+1) = numel(R) / current_total;
        H_prop(end+1) = numel(H) / current_total;
        I_count(end+1) = numel(I);
    end
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.R_prop = R_prop;
    result.H_prop = H_prop;
    result.I_count = I_count;
    result.final_time = t;
    result.peak_infected = max(I_count);
    result.peak_time = T(find(I_count == max(I_count), 1, 'first')); % Handle ties by taking first peak
    end


function det_result = solve_deterministic_sir()
    % Solve the deterministic SIR model using ODE45
    
    % Model parameters (same as stochastic model)
    beta = 0.95;    % infection constant
    gamma = 1;
    alpha = 1;
    p1 = 0.5;
    p2 = 0.5;
    p3 = 0.5;
    ph = 0.5;
    
    % Initial conditions (proportions)
    s0 = 0.96;
    i0 = 0.04;
    h0 = 0.0;
    r0 = 0.0;
    
    % Time span (adjusted to 0 to 25)
    tspan = [0, 25];
    
    % Initial conditions vector
    y0 = [s0; i0; h0; r0];
    
    % Define the ODE system
    
    
    ode_system = @(t, y) [
        -p1*beta*y(1)*y(2);           % ds/dt
        p1*beta*y(1)*y(2) - p2*gamma*y(2);  % di/dt
        ph*p2*gamma*y(2) - alpha*p3*y(3); % dh/dt
        (1-ph)*p2*gamma*y(2) + alpha*p3*y(3)           % dr/dt
    ];
    
    % Solve the ODE system
    [T, Y] = ode45(ode_system, tspan, y0);
    
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
end

function plot_comparison(results, N_values, det_result)
    % Create comparison plots including deterministic solution
    
    % Main figure - increased size
    figure('Position', [100, 100, 1920, 1440]);
    
    % Use tiledlayout for better spacing control
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Colors for different population sizes and deterministic
    colors = ['b', 'g', 'r'];
    det_color = 'm';  % Magenta for deterministic
    
    % Plot 1: Susceptible Proportion Over Time
    nexttile;
    hold on;
    plot_handles = [];
    for i = 1:length(results)
        h = plot(results{i}.T, results{i}.S_prop, [colors(i) '-'], 'LineWidth', 1.5);
        plot_handles = [plot_handles, h];
    end
    h_det = plot(det_result.T, det_result.S_prop, [det_color '--'], 'LineWidth', 2);
    plot_handles = [plot_handles, h_det];
    xlabel('Time', 'Color', 'k', 'FontSize', 14);
    ylabel('Proportion Susceptible', 'Color', 'k', 'FontSize', 14);
    title('Susceptible Proportion Over Time', 'Color', 'k', 'FontSize', 16);
    grid on;
    xlim([0, 25]);
    
    % Plot 2: Infected Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.I_prop, [colors(i) '-'], 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.I_prop, [det_color '--'], 'LineWidth', 2);
    xlabel('Time', 'Color', 'k', 'FontSize', 14);
    ylabel('Proportion Infected', 'Color', 'k', 'FontSize', 14);
    title('Infected Proportion Over Time', 'Color', 'k', 'FontSize', 16);
    grid on;
    xlim([0, 25]);
    
    % Plot 3: Recovered Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.R_prop, [colors(i) '-'], 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.R_prop, [det_color '--'], 'LineWidth', 2);
    xlabel('Time', 'Color', 'k', 'FontSize', 14);
    ylabel('Proportion Recovered', 'Color', 'k', 'FontSize', 14);
    title('Recovered Proportion Over Time', 'Color', 'k', 'FontSize', 16);
    grid on;
    xlim([0, 25]);

    % Plot 4: Hospitalized Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.H_prop, [colors(i) '-'], 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.H_prop, [det_color '--'], 'LineWidth', 2);
    xlabel('Time', 'Color', 'k', 'FontSize', 14);
    ylabel('Proportion Hospialized', 'Color', 'k', 'FontSize', 14);
    title('Recovered Proportion Over Time', 'Color', 'k', 'FontSize', 16);
    grid on;
    xlim([0, 25]);
    
    % Add a single legend below all plots
    lgd = legend(plot_handles, [arrayfun(@(x) sprintf('N=%d', x), N_values, 'UniformOutput', false), {'Deterministic'}], ...
        'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12);
    lgd.Layout.Tile = 'south';
    
    % Add beta and gamma display at the bottom-left
    annotation('textbox', [0.05, 0.02, 0.1, 0.05], ...
        'String', sprintf('β = %.1f, γ = %.1f', 0.95, 1.0), ...
        'FontSize', 14, ...
        'Color', 'k', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    
    % Print summary statistics
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('Population Size | Peak Infected | Peak Time | Final Time\n');
    fprintf('----------------|---------------|-----------|------------\n');
    for i = 1:length(results)
        fprintf('%15d | %13d | %9.2f | %10.2f\n', ...
            results{i}.N, results{i}.peak_infected, ...
            results{i}.peak_time, results{i}.final_time);
    end
    fprintf('%15s | %13.4f | %9.2f | %10.2f\n', ...
        'Deterministic', det_result.peak_infected_prop, ...
        det_result.peak_time, det_result.final_time);
    
    % Convergence analysis
    fprintf('\n=== CONVERGENCE ANALYSIS ===\n');
    fprintf('Population Size | Peak Infected (Normalized) | Difference from Deterministic\n');
    fprintf('----------------|----------------------------|------------------------------\n');
    for i = 1:length(results)
        normalized_peak = results{i}.peak_infected / results{i}.N;
        diff_from_det = abs(normalized_peak - det_result.peak_infected_prop);
        fprintf('%15d | %26.4f | %28.4f\n', ...
            results{i}.N, normalized_peak, diff_from_det);
    end
    fprintf('%15s | %25.4f | %28s\n', ...
        'Deterministic', det_result.peak_infected_prop, 'Reference');
    
    % Efficiency check
    fprintf('\n=== EFFICIENCY CHECK ===\n');
    for i = 1:length(results)
        fprintf('Population Size N = %d: Peak Infected = %d, Final Time = %.2f\n', ...
            results{i}.N, results{i}.peak_infected, results{i}.final_time);
    end
    fprintf('Deterministic: Peak Infected (proportion) = %.2f, Final Time = %.2f\n', ...
        det_result.peak_infected_prop, det_result.final_time);
end

% Run the simulation
sir_multiple_populations();