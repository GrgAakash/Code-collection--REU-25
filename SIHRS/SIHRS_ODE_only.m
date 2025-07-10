function sihrs_ode_only()
    % SIHRS Model - ODE Solution Only
    % Based on SIHRS.m parameters
    
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

    % Solve deterministic model
    deterministic_result = solve_deterministic_sihrs(params);
    
    % Plot ODE results
    plot_ode_results(deterministic_result, params);
    
    % Print analysis
    print_analysis(deterministic_result, params);
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

function plot_ode_results(det_result, params)
    % Plot ODE results - main plot only
    
    % Create single figure with all compartments
    figure('Position', [100, 100, 1000, 600]);
    
    % Plot all compartments together
    plot(det_result.T, det_result.S_prop, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible');
    hold on;
    plot(det_result.T, det_result.I_prop, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(det_result.T, det_result.H_prop, 'm-', 'LineWidth', 2, 'DisplayName', 'Hospitalized');
    plot(det_result.T, det_result.R_prop, 'g-', 'LineWidth', 2, 'DisplayName', 'Recovered');
    plot(det_result.T, det_result.D_prop, 'k-', 'LineWidth', 2, 'DisplayName', 'Dead');
    
    xlabel('Time (days)', 'FontSize', 14);
    ylabel('Proportion', 'FontSize', 14);
    title('SIHRS Model - All Compartments', 'FontSize', 16);
    legend('Location', 'eastoutside', 'FontSize', 12);
    grid on;
    xlim([0, params.tmax]);
    ylim([0, 1]);
    
    % Add parameter text
    param_text = sprintf('R₀=%.2f, β=%.4f, γ=%.4f, α=%.3f, Λ=%.4f', ...
        det_result.R0, params.beta, params.gamma, params.alpha, params.lambda);
    text(0.02, 0.98, param_text, 'Units', 'normalized', 'FontSize', 12, ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'white');
    
    % Save the figure
    saveas(gcf, 'SIHRS_ODE_main_plot.png');
end

function print_analysis(det_result, params)
    % Print analysis results
    
    fprintf('\n=== SIHRS ODE ANALYSIS ===\n');
    fprintf('R₀ = %.4f\n', det_result.R0);
    fprintf('Peak infected proportion = %.4f at day %.1f\n', det_result.peak_infected_prop, det_result.peak_time);
    
    fprintf('\nAsymptotic Values (t = %.0f days):\n', params.tmax);
    fprintf('s(∞) = %.6f\n', det_result.s_inf);
    fprintf('i(∞) = %.6f\n', det_result.i_inf);
    fprintf('h(∞) = %.6f\n', det_result.h_inf);
    fprintf('r(∞) = %.6f\n', det_result.r_inf);
    fprintf('d(∞) = %.6f\n', det_result.d_inf);
    
    fprintf('\nParameter Summary:\n');
    fprintf('β = %.4f (infection rate)\n', params.beta);
    fprintf('γ = %.4f (recovery rate)\n', params.gamma);
    fprintf('α = %.3f (hospital recovery rate)\n', params.alpha);
    fprintf('λ = %.4f (immunity loss rate)\n', params.lambda);
    fprintf('Immunity duration = %.1f days\n', 1/params.lambda);
    
    fprintf('\nDeath Rate Analysis:\n');
    death_rate_infection = params.gamma * params.pID;
    death_rate_hospital = params.alpha * params.pHD;
    total_death_rate = death_rate_infection + death_rate_hospital;
    fprintf('Death rate from infection = %.6f\n', death_rate_infection);
    fprintf('Death rate from hospitalization = %.6f\n', death_rate_hospital);
    fprintf('Total death rate = %.6f\n', total_death_rate);
end 