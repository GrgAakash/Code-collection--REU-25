function sihmd_dynamic_parameter()
    % Define model parameters structure
    params = struct(...
        'beta', 1.0, ...    % infection rate (β > 0)
        'hosp_threshold', 0.05, ...     % point at which pIH changes
        'gamma', 0.5, ...   % I transition rate (γ > 0)
        'alpha', 1.0, ...   % H transition rate (α > 0)
        'pI', 1.0, ...      % probability of S to I (p_I in (0,1])
        'pIS', 0.6, ... % probability of I to S after threshold
        'pIH1', 0.08, ...    % probability of I to H before threshold
        'pIH2', 0.00000001, ...    % probability of I to H after threshold
        'pIM', 0.2, ...     % probability of I to M
        'pID', 0.010175, ... % probability of I to D
        'pII1', 1 - 0.6 - 0.08 - 0.2 - 0.010175, ... % probability of I to I
        'pII2', 1 - 0.6 - 0.00000001 - 0.2 - 0.010175, ... % probability of I to I after threshold
        'pHS', 0.25, ...    % probability of H to S
        'pHM', 0.6, ...     % probability of H to M
        'pHD', 0.15, ...    % probability of H to D
        'tmax', 30, ...     % simulation end time
        's0', 0.98, ...      % initial susceptible proportion
        'i0', 0.02, ...      % initial infected proportion
        'h0', 0.0, ...      % initial hospitalized proportion
        'm0', 0.0, ...      % initial immune proportion
        'd0', 0.0);     % initial dead proportion
    % Validate parameters
    validateParameters(params);
    
    % Validate initial conditions sum to 1
    if abs((params.s0 + params.i0 + params.h0 + params.m0 + params.d0) - 1) > 1e-10
        error('Initial conditions must sum to 1');
    end

    % Population sizes to test
    N_values = [1749];
  
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
            results{idx} = sihmd_agent_model(N_values(idx), params);
            fprintf('Completed N = %d\n', N_values(idx));
        end
        
        % Solve deterministic model
        deterministic_result = solve_deterministic_sihmd(params);
        
        % Plot comparison
        plot_comparison(results, N_values, deterministic_result, params);
        
    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
end

function validateParameters(params)
    % Validate rates are positive
    if any([params.beta, params.gamma, params.alpha] <= 0)
        error('All rates (beta, gamma, alpha) must be positive');
    end
    
    % Validate probabilities are in (0,1]
    probs = [params.pI, params.pIS, params.pIH1, params.pIH2, params.pIM, params.pID, params.pII1, params.pII2, params.pHS, params.pHM, params.pHD];
    if any(probs <= 0 | probs > 1)
        error('All probabilities must be in (0,1]');
    end
    
    % Validate probability sums
    if abs((params.pIS + params.pIH1 + params.pIM + params.pID + params.pII1) - 1) > 1e-10
        error('I transition probabilities must sum to 1');
    end
    if abs((params.pIS + params.pIH2 + params.pIM + params.pID + params.pII2) - 1) > 1e-10
        error('I transition probabilities must sum to 1');
    end
    if abs((params.pHS + params.pHM + params.pHD) - 1) > 1e-10
        error('H transition probabilities must sum to 1');
    end
end

function result = sihmd_agent_model(N, params)
    % SIHMD agent-based stochastic model
    validateattributes(N, {'numeric'}, {'positive', 'integer', 'scalar'});
    
    % Initial conditions - using params values and ensuring they sum to N
    s0 = round(params.s0 * N); % susceptible
    i0 = round(params.i0 * N); % infected
    h0 = round(params.h0 * N); % hospitalized
    m0 = round(params.m0 * N); % immune
    d0 = round(params.d0 * N); % dead
    
    % Adjust for rounding errors to ensure sum is exactly N
    total = s0 + i0 + h0 + m0 + d0;
    if total ~= N
        % Add or subtract the difference from the largest compartment
        [~, largest_idx] = max([s0, i0, h0, m0, d0]);
        switch largest_idx
            case 1
                s0 = s0 + (N - total);
            case 2
                i0 = i0 + (N - total);
            case 3
                h0 = h0 + (N - total);
            case 4
                m0 = m0 + (N - total);
            case 5
                d0 = d0 + (N - total);
        end
    end
    
    % Validate initial conditions sum to N
    if (s0 + i0 + h0 + m0 + d0) ~= N
        error('Initial conditions must sum to N');
    end
    
    % Preallocate arrays for better performance
    max_events = N * 10; % Estimate maximum number of events
    T = zeros(1, max_events);
    S_prop = zeros(1, max_events);
    I_prop = zeros(1, max_events);
    H_prop = zeros(1, max_events);
    M_prop = zeros(1, max_events);
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
    
    M = zeros(1, N);
    M(1:m0) = (s0+i0+h0+1):(s0+i0+h0+m0);
    M = M(1:m0);
    
    D = zeros(1, N);
    D(1:d0) = (s0+i0+h0+m0+1):(s0+i0+h0+m0+d0);
    D = D(1:d0);
    
    % Initialize time tracking
    t = 0;
    T(1) = 0;
    event_count = 1;
    
    % Initialize proportion tracking
    total_pop = s0 + i0 + h0 + m0 + d0;
    S_prop(1) = s0 / total_pop;
    I_prop(1) = i0 / total_pop;
    H_prop(1) = h0 / total_pop;
    M_prop(1) = m0 / total_pop;
    D_prop(1) = d0 / total_pop;
    I_count(1) = i0;
    
    % Main simulation loop
    while ~isempty(I) && t < params.tmax
        % Check if we need to change p_IH
        if H_prop(event_count) >= params.hosp_threshold
            pIH = params.pIH2;
            pII = params.pII2;
        else
            pIH = params.pIH1;
            pII = params.pII1;
        end
        nS = numel(S);
        nI = numel(I);
        nH = numel(H);
        
        % Calculate event rates according to the mathematical model
        infection_rate = params.pI * params.beta * nS * nI / N;  % S to I rate
        to_susceptible_from_I_rate = params.gamma * nI * params.pIS;  % I to S rate
        to_hospital_rate = params.gamma * nI * pIH;  % I to H rate
        to_immune_from_I_rate = params.gamma * nI * params.pIM;  % I to M rate
        to_dead_from_I_rate = params.gamma * nI * params.pID;  % I to D rate
        to_susceptible_from_H_rate = params.alpha * nH * params.pHS;  % H to S rate
        to_immune_from_H_rate = params.alpha * nH * params.pHM;  % H to M rate
        to_dead_from_H_rate = params.alpha * nH * params.pHD;  % H to D rate
        
        total_rate = infection_rate + to_susceptible_from_I_rate + to_hospital_rate + ...
                     to_immune_from_I_rate + to_dead_from_I_rate + to_susceptible_from_H_rate + ...
                     to_immune_from_H_rate + to_dead_from_H_rate;
        
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
            current_total = numel(S) + numel(I) + numel(H) + numel(M) + numel(D);
            S_prop(event_count) = numel(S) / current_total;
            I_prop(event_count) = numel(I) / current_total;
            H_prop(event_count) = numel(H) / current_total;
            M_prop(event_count) = numel(M) / current_total;
            D_prop(event_count) = numel(D) / current_total;
            I_count(event_count) = numel(I);
            break;
        end
        
        event_count = event_count + 1;
        T(event_count) = t;
        
        % Determine which event occurs
        % NOTE: To allow I agents to stay I after clock advancement with probability pII,
        % we must include an explicit "I to I" event in the event selection.
        % First, define pIS and pII for this time step:
        if H_prop(event_count) >= params.hosp_threshold
            pIH = params.pIH2;
            pII = params.pII2;
        else
            pIH = params.pIH1;
            pII = params.pII1;
        end

        % Recalculate rates for I transitions (including I->I)
        to_susceptible_from_I_rate = params.gamma * nI * params.pIS;  % I to S
        to_hospital_rate = params.gamma * nI * pIH;             % I to H
        to_immune_from_I_rate = params.gamma * nI * params.pIM; % I to M
        to_dead_from_I_rate = params.gamma * nI * params.pID;   % I to D
        to_stay_infected_rate = params.gamma * nI * pII;        % I to I (no change)

        % Update total_rate to include I->I
        total_rate = infection_rate + to_susceptible_from_I_rate + to_hospital_rate + ...
                 to_immune_from_I_rate + to_dead_from_I_rate + to_stay_infected_rate + ...
                 to_susceptible_from_H_rate + to_immune_from_H_rate + to_dead_from_H_rate;

        chance = rand * total_rate;
        if chance < infection_rate
            % S to I transition
            if nS > 0
            num = randi([1, nS]);
            infected_agent = S(num);
            S(num) = [];
            I(end+1) = infected_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate)
            % I to S transition
            if nI > 0
            num = randi([1, nI]);
            susceptible_agent = I(num);
            I(num) = [];
            S(end+1) = susceptible_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate)
            % I to H transition
            if nI > 0
            num = randi([1, nI]);
            hospitalized_agent = I(num);
            I(num) = [];
            H(end+1) = hospitalized_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate + to_immune_from_I_rate)
            % I to M transition
            if nI > 0
            num = randi([1, nI]);
            immune_agent = I(num);
            I(num) = [];
            M(end+1) = immune_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate + to_immune_from_I_rate + to_dead_from_I_rate)
            % I to D transition
            if nI > 0
            num = randi([1, nI]);
            dead_agent = I(num);
            I(num) = [];
            D(end+1) = dead_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate + to_immune_from_I_rate + to_dead_from_I_rate + to_stay_infected_rate)
            % I to I transition (agent stays infected, no change)
            % No action needed, just advance time
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate + to_immune_from_I_rate + to_dead_from_I_rate + to_stay_infected_rate + to_susceptible_from_H_rate)
            % H to S transition
            if nH > 0
            num = randi([1, nH]);
            susceptible_agent = H(num);
            H(num) = [];
            S(end+1) = susceptible_agent;
            end
        elseif chance < (infection_rate + to_susceptible_from_I_rate + to_hospital_rate + to_immune_from_I_rate + to_dead_from_I_rate + to_stay_infected_rate + to_susceptible_from_H_rate + to_immune_from_H_rate)
            % H to M transition
            if nH > 0
            num = randi([1, nH]);
            immune_agent = H(num);
            H(num) = [];
            M(end+1) = immune_agent;
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
        current_total = numel(S) + numel(I) + numel(H) + numel(M) + numel(D);
        S_prop(event_count) = numel(S) / current_total;
        I_prop(event_count) = numel(I) / current_total;
        H_prop(event_count) = numel(H) / current_total;
        M_prop(event_count) = numel(M) / current_total;
        D_prop(event_count) = numel(D) / current_total;
        I_count(event_count) = numel(I);
    end
    
    % Trim unused preallocated space
    T = T(1:event_count);
    S_prop = S_prop(1:event_count);
    I_prop = I_prop(1:event_count);
    H_prop = H_prop(1:event_count);
    M_prop = M_prop(1:event_count);
    D_prop = D_prop(1:event_count);
    I_count = I_count(1:event_count);
    
    % Store results
    result.N = N;
    result.T = T;
    result.S_prop = S_prop;
    result.I_prop = I_prop;
    result.H_prop = H_prop;
    result.M_prop = M_prop;
    result.D_prop = D_prop;
    result.I_count = I_count;
    result.final_time = t;
    result.peak_infected = max(I_count);
    result.peak_time = T(find(I_count == max(I_count), 1, 'first'));
    
    % Calculate and store asymptotic values
    result.s_inf = S_prop(end);
    result.i_inf = I_prop(end);
    result.h_inf = H_prop(end);
    result.m_inf = M_prop(end);
    result.d_inf = D_prop(end);
end

function det_result = solve_deterministic_sihmd(params)
    % Solve the deterministic SIHMD model using ODE45

    % Time span
    tspan = [0, params.tmax];

    % Define time-dependent pIH function
    % Use a persistent variable to store the first time H_prop crosses the threshold
    persistent t_threshold
    t_threshold = [];

    % ODE system with dynamic pIH and pIS based on H proportion
    function dydt = ode_system_with_dynamic_pIH(t, y)
        H_prop = y(3);
        if isempty(t_threshold) && (H_prop >= params.hosp_threshold)
            t_threshold = t;
        end
        if ~isempty(t_threshold) && (t >= t_threshold)
            pIH = params.pIH2;
            pII = params.pII2;
        else
            pIH = params.pIH1;
            pII = params.pII1;
        end
        % In the ODE model, the probability of staying infected (pII) reduces the outflow from I.
        % All transitions out of I are multiplied by (1 - pII).
        gamma_eff = params.gamma * (1 - pII);
        dydt = [
            -params.beta * y(1) * y(2) * params.pI + gamma_eff * y(2) * params.pIS + params.alpha * y(3) * params.pHS; % ds/dt
            params.beta * y(1) * y(2) * params.pI - gamma_eff * y(2); % di/dt
            gamma_eff * y(2) * pIH - params.alpha * y(3);             % dh/dt
            gamma_eff * y(2) * params.pIM + params.alpha * y(3) * params.pHM; % dm/dt
            gamma_eff * y(2) * params.pID + params.alpha * y(3) * params.pHD  % dd/dt
        ];
    end

    % Initial conditions vector
    y0 = [params.s0; params.i0; params.h0; params.m0; params.d0];
    ode_system = @ode_system_with_dynamic_pIH;

    % Define the ODE system using time-dependent pIH
    % ode_system = @(t, y) [
    %     -params.beta * y(1) * y(2) * params.pI + params.gamma * y(2) * pIS_func(t) + params.alpha * y(3) * params.pHS; % ds/dt
    %     params.beta * y(1) * y(2) * params.pI - params.gamma * y(2); % di/dt
    %     params.gamma * y(2) * pIH_func(t) - params.alpha * y(3);                                                        % dh/dt
    %     params.gamma * y(2) * params.pIM + params.alpha * y(3) * params.pHM;                                           % dm/dt
    %     params.gamma * y(2) * params.pID + params.alpha * y(3) * params.pHD                                            % dd/dt
    % ];

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
    det_result.M_prop = Y(:, 4);
    det_result.D_prop = Y(:, 5);

    % Find peak infected and peak time
    [peak_infected_prop, peak_idx] = max(det_result.I_prop);
    det_result.peak_infected_prop = peak_infected_prop;
    det_result.peak_time = T(peak_idx);
    det_result.final_time = T(end);

    % Calculate R0 using beta
    det_result.R0 = params.pI * params.beta / params.gamma;

    % Store asymptotic values
    det_result.s_inf = det_result.S_prop(end);
    det_result.i_inf = det_result.I_prop(end);
    det_result.h_inf = det_result.H_prop(end);
    det_result.m_inf = det_result.M_prop(end);
    det_result.d_inf = det_result.D_prop(end);

    % Verify asymptotic behavior
    if det_result.i_inf > 1e-6
        warning('i(∞) may not be approaching 0 as expected');
    end
    if det_result.h_inf > 1e-6
        warning('h(∞) may not be approaching 0 as expected');
    end
end

function plot_comparison(results, N_values, det_result, params)
    % Create comparison plots including deterministic solution
    
    % Create two figures
    figure('Position', [100, 100, 1920, 1440]);  % Original comparison plots
    
    % Use tiledlayout for better spacing control
    t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Colors for different population sizes and deterministic
    colors = {'#0072BD', '#77AC30', '#A2142F', '#FF00FF'}; % More distinctive colors
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
    
    % Plot 4: Immune Proportion Over Time
    nexttile;
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.M_prop, 'Color', colors{i}, 'LineWidth', 1.5);
    end
    plot(det_result.T, det_result.M_prop, '--', 'Color', det_color, 'LineWidth', 2);
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Immune', 'FontSize', 14);
    title('Immune Proportion Over Time', 'FontSize', 16);
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
    param_text = sprintf('R₀=%.2f, β=%.2f, γ=%.2f, α=%.2f\np_I=%.2f, p_II1=%.2f, p_II2=%.2f, p_IH1=%.2f, p_IH2 = %.2f, p_IM=%.2f, p_ID=%.2f\np_HS=%.2f, p_HM=%.2f, p_HD=%.2f', ...
        det_result.R0, params.beta, params.gamma, params.alpha, ...
        params.pI, params.pII1, params.pII2, params.pIH1, params.pIH2, params.pIM, params.pID, ...
        params.pHS, params.pHM, params.pHD);
    annotation('textbox', [0.05, 0.02, 0.9, 0.05], ...
        'String', param_text, ...
        'FontSize', 14, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
    
    
    % Create new figure for deterministic curves only
    figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Plot all deterministic curves together
    plot(det_result.T, det_result.S_prop, 'LineWidth', 2, 'DisplayName', 'Susceptible');
    plot(det_result.T, det_result.I_prop, 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(det_result.T, det_result.H_prop, 'LineWidth', 2, 'DisplayName', 'Hospitalized');
    plot(det_result.T, det_result.M_prop, 'LineWidth', 2, 'DisplayName', 'Immune');
    plot(det_result.T, det_result.D_prop, 'LineWidth', 2, 'DisplayName', 'Dead');
    
    % Customize the plot
    xlabel('Time', 'FontSize', 14);
    ylabel('Population Proportion', 'FontSize', 14);
    title('Deterministic SIHMD Model Dynamics', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    legend('Location', 'east', 'FontSize', 12);
    
    % Add R0 and parameters text
    text(0.02, -0.15, sprintf('R₀=%.2f, β=%.2f, γ=%.2f, α=%.2f', det_result.R0, params.beta, params.gamma, params.alpha), ...
        'Units', 'normalized', 'FontSize', 12);
    
    % Plot just infected population over time (all N and deterministic)
    figure('Position', [200, 200, 900, 500]);
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.I_prop, 'Color', colors{i}, 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', results{i}.N));
    end
    plot(det_result.T, det_result.I_prop, '--', 'Color', det_color, 'LineWidth', 2.5, 'DisplayName', 'Deterministic');
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion Infected', 'FontSize', 14);
    title('Infected Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    legend('Location', 'northeast', 'FontSize', 12);
    saveas(gcf, 'SIHMD_infected_only.png');

    % Plot just hospitalized population over time (all N and deterministic)
    figure('Position', [200, 200, 900, 500]);
    hold on;
    for i = 1:length(results)
        plot(results{i}.T, results{i}.H_prop, 'Color', colors{i}, 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', results{i}.N));
    end
    plot(det_result.T, det_result.H_prop, '--', 'Color', det_color, 'LineWidth', 2.5, 'DisplayName', 'Deterministic');
    xlabel('Time', 'FontSize', 14);
    ylabel('Proportion of Population', 'FontSize', 14);
    title('Hospitalized Proportion Over Time', 'FontSize', 16);
    grid on;
    xlim([0, params.tmax]);
    legend('Location', 'northeast', 'FontSize', 12);
    saveas(gcf, 'SIHMD_hospitalized_only.png');
    
    % Print summary statistics
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('Population Size | Peak Infected | Peak Time | s(∞) | i(∞) | h(∞) | m(∞) | d(∞)\n');
    fprintf('----------------|---------------|-----------|-------|-------|-------|-------|-------\n');
    for i = 1:length(results)
        fprintf('%15d | %13d | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
            results{i}.N, results{i}.peak_infected, results{i}.peak_time, ...
            results{i}.s_inf, results{i}.i_inf, results{i}.h_inf, results{i}.m_inf, results{i}.d_inf);
    end
    fprintf('%15s | %13.4f | %9.2f | %5.3f | %5.3f | %5.3f | %5.3f | %5.3f\n', ...
        'Deterministic', det_result.peak_infected_prop, det_result.peak_time, ...
        det_result.s_inf, det_result.i_inf, det_result.h_inf, det_result.m_inf, det_result.d_inf);
    
    % Print asymptotic analysis
    fprintf('\n=== ASYMPTOTIC ANALYSIS ===\n');
    fprintf('R₀ = %.4f\n', det_result.R0);
    fprintf('i(∞) = %.6f (should be ≈ 0)\n', det_result.i_inf);
    fprintf('h(∞) = %.6f (should be ≈ 0)\n', det_result.h_inf);
    fprintf('m(∞) = %.6f\n', det_result.m_inf);
    fprintf('d(∞) = %.6f\n', det_result.d_inf);

    % Export simulation data to CSV for Python plotting
    % Export agent-based results (first N only, for simplicity)
    if ~isempty(results)
        sim_table = table(results{1}.T', results{1}.S_prop', results{1}.I_prop', results{1}.H_prop', ...
                          results{1}.M_prop', results{1}.D_prop', results{1}.I_count', ...
                          'VariableNames', {'Time', 'S', 'I', 'H', 'M', 'D', 'I_count'});
        writetable(sim_table, sprintf('SIHMD_agent_N%d.csv', results{1}.N));
    end

    % Export deterministic results
    det_table = table(det_result.T, det_result.S_prop, det_result.I_prop, det_result.H_prop, ...
                      det_result.M_prop, det_result.D_prop, ...
                      'VariableNames', {'Time', 'S_prop', 'I_prop', 'H_prop', 'M_prop', 'D_prop'});
    writetable(det_table, 'SIHMD_deterministic.csv');

end