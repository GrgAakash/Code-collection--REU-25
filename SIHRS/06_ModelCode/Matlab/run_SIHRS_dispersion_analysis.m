% SIHRS Dispersion Time Derivative Analysis (Equation 5.17)
%
% The time derivative of the dispersion quantity:
% d/dt [ sigma(M_N^{(l)}(t)) / E[m_N^{(l)}(t)] ]
%
% Where sigma(M_N^{(l)}(t)) = sqrt( ∫_0^t E[V_N^{(l)}(w)] dw )
% 
% Using the formula from Equation 5.17 in the paper:
% d/dt = [ E[m] * E[V] - 2 * (sigma(M))^2 * E[G] ] / [ 2 * sigma(M) * (E[m])^2 ]

clear all;
close all;

function sihrs_dispersion_derivative_analysis()
    rng(1);
    
    % Parameters
    params.beta = 0.212;
    params.pSI = 1.0;
    params.pII = 0.0;
    params.pIH = 0.1060;
    params.pIR = 0.8921;
    params.pID = 0.0019;
    params.pHH = 0.00;
    params.pHR = 0.846;
    params.pHD = 0.154;
    params.pRR = 0.02;
    params.pRS = 0.98;
    params.gamma = 0.100346667;
    params.alpha = 0.1;
    params.lambda = 0.0083;
    params.T = 1000;
    params.dt = 0.01;
    params.N_values = [1600, 3000, 6000];
    params.initial_s = 0.96;
    params.initial_i = 0.04;
    params.initial_h = 0;
    params.initial_r = 0;
    params.initial_d = 0;
    params.n_runs = 100;
    params.colors = {'#0072BD','#77AC30', '#A2142F'};
    
    validate_params(params);
    
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('R₀ = %.4f\n', R0);
    
    simulate_and_analyze_dispersion_derivative(params, R0);
end

function validate_params(params)
    if any(params.N_values <= 0) || any(mod(params.N_values, 1) ~= 0)
        error('N_values must be positive integers');
    end
    if params.pSI <= 0 || params.pSI > 1 || params.pII < 0 || params.pII > 1
        error('pSI must be in (0, 1] and pII must be in [0, 1]');
    end
    if params.gamma <= 0
        error('gamma must be positive');
    end
    if abs(params.initial_s + params.initial_i + params.initial_h + params.initial_r + params.initial_d - 1) > 1e-6
        error('Initial fractions must sum to 1');
    end
    if params.T <= 0 || params.dt <= 0
        error('T and dt must be positive');
    end
    if params.n_runs <= 0 || mod(params.n_runs, 1) ~= 0
        error('n_runs must be a positive integer');
    end
end

function simulate_and_analyze_dispersion_derivative(params, R0)
    beta = params.beta;
    t = 0:params.dt:params.T;
    results = cell(length(params.N_values), 1);
    
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d...\n', params.n_runs, N);
        results{idx} = run_dispersion_derivative_analysis(N, beta, params, t);
        fprintf('Done N = %d\n', N);
    end
    
    plot_dispersion_derivative_results(t, results, params, R0);
end

function result = run_dispersion_derivative_analysis(N, beta, params, t)
    % Calculate m, G, V for each path, then average
    m_S_all = zeros(params.n_runs, length(t));
    m_I_all = zeros(params.n_runs, length(t));
    m_H_all = zeros(params.n_runs, length(t));
    G_S_all = zeros(params.n_runs, length(t));
    G_I_all = zeros(params.n_runs, length(t));
    G_H_all = zeros(params.n_runs, length(t));
    V_S_all = zeros(params.n_runs, length(t));
    V_I_all = zeros(params.n_runs, length(t));
    V_H_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        
        m_S_all(run, :) = S_interp;
        m_I_all(run, :) = I_interp;
        m_H_all(run, :) = H_interp;
        
    
        G_S_all(run, :) = -beta * S_interp .* I_interp + params.lambda * R_interp;
        G_I_all(run, :) = beta * S_interp .* I_interp - params.gamma * I_interp;
        G_H_all(run, :) = params.gamma * params.pIH * I_interp - params.alpha * H_interp;
        
    
        V_S_all(run, :) = (1/N) * (beta * params.pSI * S_interp .* I_interp + params.lambda * params.pRS * R_interp);
        V_I_all(run, :) = (1/N) * (params.gamma * I_interp * (params.pIH + params.pIR + params.pID) + beta * params.pSI * S_interp .* I_interp);
        V_H_all(run, :) = (1/N) * (params.gamma * I_interp * params.pIH + params.alpha * H_interp * (params.pHR + params.pHD));
    end
    
    % Average to get E[m], E[G], E[V]
    Em_S = mean(m_S_all, 1);
    Em_I = mean(m_I_all, 1);
    Em_H = mean(m_H_all, 1);
    EG_S = mean(G_S_all, 1);
    EG_I = mean(G_I_all, 1);
    EG_H = mean(G_H_all, 1);
    EV_S = mean(V_S_all, 1);
    EV_I = mean(V_I_all, 1);
    EV_H = mean(V_H_all, 1);
    
    % Compute dispersion
    intEV_S = cumtrapz(t, EV_S);
    intEV_I = cumtrapz(t, EV_I);
    intEV_H = cumtrapz(t, EV_H);
    
    sigmaM_S = sqrt(max(intEV_S, 0));
    sigmaM_I = sqrt(max(intEV_I, 0));
    sigmaM_H = sqrt(max(intEV_H, 0));
    
    Disp_S = sigmaM_S ./ Em_S;
    Disp_I = sigmaM_I ./ Em_I;
    Disp_H = sigmaM_H ./ Em_H;
    
    % Derivative via Eq. 5.17
    dDisp_S = compute_dispersion_derivative(t, Em_S, EG_S, EV_S);
    dDisp_I = compute_dispersion_derivative(t, Em_I, EG_I, EV_I);
    dDisp_H = compute_dispersion_derivative(t, Em_H, EG_H, EV_H);
    
    % Store
    result.Disp1 = Disp_S;
    result.Disp2 = Disp_I;
    result.Disp3 = Disp_H;
    result.dDisp1 = dDisp_S;
    result.dDisp2 = dDisp_I;
    result.dDisp3 = dDisp_H;
    result.N = N;
end

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params)
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
    max_events = round(10 * params.T * (beta + params.gamma + params.alpha + params.lambda) * N);
    S_hist = zeros(1, max_events);
    I_hist = zeros(1, max_events);
    H_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events);
    D_hist = zeros(1, max_events);
    time_pts = zeros(1, max_events);
    
    S_hist(1) = S;
    I_hist(1) = I;
    H_hist(1) = H;
    R_hist(1) = R;
    D_hist(1) = D;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    while current_time < params.T
        si_rate = (beta / N) * S * I * params.pSI;
        ir_rate = params.gamma * I * params.pIR;
        ih_rate = params.gamma * I * params.pIH;
        id_rate = params.gamma * I * params.pID;
        hr_rate = params.alpha * H * params.pHR;
        hd_rate = params.alpha * H * params.pHD;
        rs_rate = params.lambda * R * params.pRS;
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        
        if total_rate == 0
            break;
        end
        
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        
        if current_time > params.T
            break;
        end
        
        chance = rand * total_rate;
        if chance < si_rate
            if S > 0; S = S - 1; I = I + 1; end
        elseif chance < (si_rate + ir_rate)
            if I > 0; I = I - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate)
            if I > 0; I = I - 1; H = H + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            if I > 0; I = I - 1; D = D + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            if H > 0; H = H - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            if H > 0; H = H - 1; D = D + 1; end
        else
            if R > 0; R = R - 1; S = S + 1; end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S;
        I_hist(event_count) = I;
        H_hist(event_count) = H;
        R_hist(event_count) = R;
        D_hist(event_count) = D;
        time_pts(event_count) = current_time;
    end
    
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    D_hist = D_hist(1:event_count);
    time_pts = time_pts(1:event_count);
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
    I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
    H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
    R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
    D_interp = interp1(time_pts, D_hist, t, 'previous') / N;
    
    S_interp(t > max(time_pts)) = S_hist(end) / N;
    I_interp(t > max(time_pts)) = I_hist(end) / N;
    H_interp(t > max(time_pts)) = H_hist(end) / N;
    R_interp(t > max(time_pts)) = R_hist(end) / N;
    D_interp(t > max(time_pts)) = D_hist(end) / N;
end

function dDisp_dt = compute_dispersion_derivative(t, Em, EG, EV)
    t  = t(:);
    Em = Em(:);
    EG = EG(:);
    EV = EV(:);
    
    intEV = cumtrapz(t, EV);
    sigmaM = sqrt(max(intEV, 0));
    
    numer = Em .* EV - 2 .* (sigmaM.^2) .* EG;
    denom = 2 .* sigmaM .* (Em.^2);
    
    dDisp_dt = numer ./ denom;
    invalid = (denom == 0) | ~isfinite(denom) | (sigmaM == 0);
    dDisp_dt(invalid) = NaN;
end

function plot_dispersion_derivative_results(t, results, params, R0)
    figure('Position', [100, 100, 1500, 900]);
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Susceptible - Dispersion
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.Disp1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Dispersion - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('σ(M)/E[m]');
    grid on;
    
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.dDisp1, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Derivative - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('d/dt [σ(M)/E[m]]');
    grid on;
    
    % Infected
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.Disp2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Dispersion - Infected (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('σ(M)/E[m]');
    grid on;
    
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.dDisp2, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Derivative - Infected (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('d/dt [σ(M)/E[m]]');
    grid on;
    
    % Hospitalized
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.Disp3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Dispersion - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('σ(M)/E[m]');
    grid on;
    
    nexttile;
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.dDisp3, 'Color', params.colors{idx}, 'LineWidth', 1.5);
    end
    title(sprintf('Derivative - Hospitalized (R_0 = %.2f)', R0));
    xlabel('Time'); ylabel('d/dt [σ(M)/E[m]]');
    grid on;
    
    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    lgd = legend(legend_labels, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';
    
    saveas(gcf, sprintf('SIHRS_Dispersion_Derivative_R0_%.2f_combined.png', R0));
    fprintf('Saved: SIHRS_Dispersion_Derivative_R0_%.2f_combined.png\n', R0);
end

% Run the analysis
sihrs_dispersion_derivative_analysis();

