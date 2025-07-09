clear all;
close all;
params = struct(...
        'beta', 30.00, ...    % infection rate (β > 0)
        'gamma1', 1.5, ...  % I to H rate (γ₁ > 0)
        'gamma2', 1.5, ...  % I to R rate (γ₂ > 0)
        'alpha', 0.40, ...   % H to R rate (α > 0)
        'p1', 0.5, ...      % probability of infection (p₁ in (0,1])
        'p2', 0.5, ...      % probability of leaving I (p₂ in (0,1])
        'p3', 0.5, ...      % probability of leaving H (p₃ in (0,1])
        'ph', 0.80, ...      % probability of I to H vs R (p_h in (0,1])
        'tmax', 30, ...     % simulation end time
        's0', 0.96, ...     % initial susceptible proportion
        'i0', 0.04, ...     % initial infected proportion
        'h0', 0.00, ...     % initial hospitalized proportion
        'r0', 0.00 ...      % initial recovered proportion
    );

% Validate parameters
validateParams(params);

% Calculate R0
R0 = params.p1 * params.beta / (params.p2 * (params.ph * params.gamma1 + (1-params.ph) * params.gamma2));
fprintf('R0 = %.4f\n', R0);

% Solve and plot
deterministic_result = solve_deterministic_sir(params);
plot_sihr_analysis(params, deterministic_result);

function validateParams(params)
    % Validate rates are positive
    if any([params.beta, params.gamma1, params.gamma2, params.alpha] <= 0)
        error('All rates (beta, gamma1, gamma2, alpha) must be positive');
    end
    
    % Validate probabilities are in (0,1]
    probs = [params.p1, params.p2, params.p3, params.ph];
    if any(probs <= 0 | probs > 1)
        error('All probabilities must be in (0,1]');
    end
    
    % Validate initial conditions sum to 1
    total = params.s0 + params.i0 + params.h0 + params.r0;
    if abs(total - 1) > 1e-10
        error('Initial conditions must sum to 1');
    end
end

function det_result = solve_deterministic_sir(params)
    % Time span
    tspan = [0, params.tmax];
    
    % Initial conditions vector
    y0 = [params.s0; params.i0; params.h0; params.r0];
    
    % Define the ODE system
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
    
    % Find peaks and their times
    [det_result.peak_I, det_result.peak_I_idx] = max(det_result.I_prop);
    [det_result.peak_H, det_result.peak_H_idx] = max(det_result.H_prop);
    det_result.peak_I_time = T(det_result.peak_I_idx);
    det_result.peak_H_time = T(det_result.peak_H_idx);
end

function plot_sihr_analysis(params, det_result)
    % Create figure with two subplots
    figure('Position', [100, 100, 1200, 800]);
    
    % Extract data
    T = det_result.T;
    S = det_result.S_prop;
    I = det_result.I_prop;
    H = det_result.H_prop;
    R = det_result.R_prop;
    
    % Calculate derivatives
    dS_dt = zeros(size(T));
    dI_dt = zeros(size(T));
    dH_dt = zeros(size(T));
    dR_dt = zeros(size(T));
    
    for idx = 1:length(T)
        dS_dt(idx) = -params.p1 * params.beta * S(idx) * I(idx);
        dI_dt(idx) = params.p1 * params.beta * S(idx) * I(idx) - ...
                     params.p2 * (params.ph * params.gamma1 + (1-params.ph) * params.gamma2) * I(idx);
        dH_dt(idx) = params.p2 * params.ph * params.gamma1 * I(idx) - params.p3 * params.alpha * H(idx);
        dR_dt(idx) = params.p2 * (1-params.ph) * params.gamma2 * I(idx) + params.p3 * params.alpha * H(idx);
    end
    
    % Plot derivatives
    subplot(2,1,1);
    plot(T, dS_dt, 'b-', 'LineWidth', 2); hold on;
    plot(T, dI_dt, 'r-', 'LineWidth', 2);
    plot(T, dH_dt, 'g-', 'LineWidth', 2);
    plot(T, dR_dt, 'm-', 'LineWidth', 2);
    
    % Find and mark critical points (where derivatives cross zero)
    for i = 1:length(T)-1
        % For dS/dt
        if dS_dt(i)*dS_dt(i+1) <= 0
            plot(T(i), 0, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
        end
        % For dI/dt
        if dI_dt(i)*dI_dt(i+1) <= 0
            plot(T(i), 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        end
        % For dH/dt
        if dH_dt(i)*dH_dt(i+1) <= 0
            plot(T(i), 0, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        end
        % For dR/dt
        if dR_dt(i)*dR_dt(i+1) <= 0
            plot(T(i), 0, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
    
    % Find and mark peaks and troughs of derivatives
    for i = 2:length(T)-1
        % For dH/dt peaks
        if dH_dt(i) > dH_dt(i-1) && dH_dt(i) > dH_dt(i+1)
            plot(T(i), dH_dt(i), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
        elseif dH_dt(i) < dH_dt(i-1) && dH_dt(i) < dH_dt(i+1)
            plot(T(i), dH_dt(i), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
        end
        
        % For dI/dt peaks
        if dI_dt(i) > dI_dt(i-1) && dI_dt(i) > dI_dt(i+1)
            plot(T(i), dI_dt(i), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
        elseif dI_dt(i) < dI_dt(i-1) && dI_dt(i) < dI_dt(i+1)
            plot(T(i), dI_dt(i), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
        end
    end
    
    xlabel('Time', 'FontSize', 12);
    ylabel('Rate of Change', 'FontSize', 12);
    title('SIHR Model Derivatives', 'FontSize', 14);
    legend('dS/dt', 'dI/dt', 'dH/dt', 'dR/dt', 'Location', 'east');
    grid on;
    
    % Plot SIHR dynamics
    subplot(2,1,2);
    plot(T, S, 'b-', 'LineWidth', 2); hold on;
    plot(T, I, 'r-', 'LineWidth', 2);
    plot(T, H, 'g-', 'LineWidth', 2);
    plot(T, R, 'm-', 'LineWidth', 2);
    
    % Add dots at peak points
    plot(det_result.peak_I_time, det_result.peak_I, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot(det_result.peak_H_time, det_result.peak_H, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    
    xlabel('Time', 'FontSize', 12);
    ylabel('Population Proportion', 'FontSize', 12);
    title('SIHR Model Dynamics', 'FontSize', 14);
    legend('S(t)', 'I(t)', 'H(t)', 'R(t)', 'Location', 'east');
    grid on;
end