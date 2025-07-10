function modsihrs_ode_only()
    % SIHRS Model - ODE Solution Only (Simplified)
    % Based on modSIHRS.m parameters
    
    % Initial conditions
    s0 = 0.96;
    i0 = 0.04;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;
    
    % Model parameters (SIHRS)
    max_t = 500;
    beta = 0.5855;        % infection rate
    gamma = 0.0714;     % recovery rate
    alpha = 0.17;      % hospital recovery rate
    lambda = 0.0083;      % immunity loss rate
    
    % Probabilities (must sum to 1 for each transition)
    pSI = 1.0;          % S to I
    pII = 0.00;         % I to I (stay infected)
    pIH = 0.015;         % I to H
    pIR = 0.984;         % I to R
    pID = 0.001;         % I to D
    pHH = 0.0;         % H to H (stay hospitalized)
    pHR = 0.97;         % H to R
    pHD = 0.03;         % H to D
    pRR = 0.2;          % R to R (stay recovered)
    pRS = 0.8;          % R to S
    
    % Validate probability constraints
    pII_sum = pII + pIH + pIR + pID;
    pHH_sum = pHH + pHR + pHD;
    pRR_sum = pRR + pRS;
    
    if abs(pII_sum - 1) > 1e-10
        error('I transition probabilities must sum to 1, got %.6f', pII_sum);
    end
    if abs(pHH_sum - 1) > 1e-10
        error('H transition probabilities must sum to 1, got %.6f', pHH_sum);
    end
    if abs(pRR_sum - 1) > 1e-10
        error('R transition probabilities must sum to 1, got %.6f', pRR_sum);
    end

    % Solve deterministic model
    [T, Y] = solve_ode_system(max_t, beta, gamma, alpha, lambda, pSI, pII, pIH, pIR, pID, pHH, pHR, pHD, pRR, pRS, s0, i0, h0, r0, d0);
    
    % Plot results
    plot_results(T, Y, max_t, beta, gamma, alpha, lambda);
    
    % Print key results
    R0 = pSI * beta / (gamma * (1 - pII));
    [peak_infected, peak_idx] = max(Y(:, 2));
    
    fprintf('\n=== SIHRS ODE ANALYSIS (modSIHRS) ===\n');
    fprintf('R₀ = %.4f\n', R0);
    fprintf('Peak infected proportion = %.4f at day %.1f\n', peak_infected, T(peak_idx));
    fprintf('Asymptotic values (t = %.0f days):\n', max_t);
    fprintf('s(∞) = %.6f, i(∞) = %.6f, h(∞) = %.6f, r(∞) = %.6f, d(∞) = %.6f\n', ...
        Y(end, 1), Y(end, 2), Y(end, 3), Y(end, 4), Y(end, 5));
end

function [T, Y] = solve_ode_system(max_t, beta, gamma, alpha, lambda, pSI, pII, pIH, pIR, pID, pHH, pHR, pHD, pRR, pRS, s0, i0, h0, r0, d0)
    % Solve the deterministic SIHRS model using ODE45
    
    % Time span
    tspan = [0, max_t];
    
    % Initial conditions vector
    y0 = [s0; i0; h0; r0; d0];
    
    % Define the ODE system
    ode_system = @(t, y) [
        -beta * y(1) * y(2) * pSI + pRS * lambda * y(4); % ds/dt
        beta * y(1) * y(2) * pSI - gamma * (1 - pII) * y(2); % di/dt
        pIH * gamma * y(2) - alpha * (1 - pHH) * y(3); % dh/dt
        pIR * gamma * y(2) + pHR * alpha * y(3) - pRS * lambda * y(4); % dr/dt
        pID * gamma * y(2) + pHD * alpha * y(3) % dd/dt
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
end

function plot_results(T, Y, max_t, beta, gamma, alpha, lambda)
    % Plot ODE results in a single comprehensive figure
    
    % Create figure
    figure('Position', [100, 100, 1000, 700]);
    
    % Plot all compartments together
    plot(T, Y(:, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible (S)');
    hold on;
    plot(T, Y(:, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'Infected (I)');
    plot(T, Y(:, 3), 'm-', 'LineWidth', 2, 'DisplayName', 'Hospitalized (H)');
    plot(T, Y(:, 4), 'g-', 'LineWidth', 2, 'DisplayName', 'Recovered (R)');
    plot(T, Y(:, 5), 'k-', 'LineWidth', 2, 'DisplayName', 'Dead (D)');
    
    xlabel('Time (days)', 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 12);
    title('SIHRS Model - All Compartments (modSIHRS)', 'FontSize', 14);
    legend('Location', 'eastoutside', 'FontSize', 10);
    grid on;
    xlim([0, max_t]);
    ylim([0, 1]);
    
    % Add parameter annotations
    param_text = sprintf('R₀=%.2f, β=%.4f, γ=%.4f, α=%.4f, λ=%.4f', ...
        beta * 1.0 / (gamma * (1 - 0.01)), beta, gamma, alpha, lambda);
    text(0.02, 0.98, param_text, 'Units', 'normalized', 'FontSize', 10, ...
        'VerticalAlignment', 'top', 'BackgroundColor', 'white');
    
    % Save the figure
    saveas(gcf, 'modSIHRS_ODE_simple.png');
end 