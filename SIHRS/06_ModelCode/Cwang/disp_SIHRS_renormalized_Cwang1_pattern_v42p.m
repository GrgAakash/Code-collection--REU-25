%% 
% SIHRS with Death Model - Renormalized Infinitesimal Standard Deviation
% 
% This file calculates the renormalized infinitesimal standard deviation:
% D_N^(l)(t) = sqrt(V_N^(l)(t)) / |m_N^(l)(t)|
% 
% Where V_N^(l)(t) are the variance formulas for each compartment l:
% l = 1: Susceptible (S)
% l = 2: Infected (I) 
% l = 3: Hospitalized (H)
% l = 4: Recovered (R)
% l = 5: Dead (D)
%
% The renormalized standard deviation helps demonstrate finite size effects
% by showing how stochastic fluctuations scale with population size N.

clear all;
close all;

function sihrs_renormalized_simulation()
    % Set random seed for reproducibility
    rng(1);
    
    % Centralized parameters for SIHRS with death - matching SIHRS.m
    params.pSI = 1.0;              % Infection probability (S to I)
    params.pII = 0.0;              % probability of I to I (stay infected)
    params.pIH = 0.04;             % probability of I to H
    params.pIR = 0.959;            % probability of I to R
    params.pID = 0.001;            % probability of I to D
    params.pHH = 0.01;             % probability of H to H (stay hospitalized)
    params.pHR = 0.9882;           % probability of H to R
    params.pHD = 0.0018;           % probability of H to D
    params.pRR = 0.02;             % probability of R to R (stay recovered)
    params.pRS = 0.98;             % probability of R to S
    params.gamma = 0.1;            % Infection transition rate (γ > 0)
    params.alpha = 0.1;            % Hospitalized transition rate (α > 0)
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0) immunity period of 4 months
    params.T = 1000;               % Total simulation time
    params.dt = 0.01;              % Time step for interpolation

    % Population sizes - matching SIHRS.m
    params.N_values = [650, 1300, 2900, 6000, 15000, 20000, 40000];

    % Colors for each N
    params.colors = num2cell(lines(numel(params.N_values)), 2);   % each cell is [r g b]

    % Initial fractions
    params.initial_s = 0.96;       
    params.initial_i = 0.035;      
    params.initial_h = 0.005;      
    params.initial_r = 0.000;          
    params.initial_d = 0.000;          

    params.n_runs   = 30;          % Number of stochastic runs (used for means below)
    params.R0_values = [2.12];     % R0 values
    
    % Validate parameters
    validate_params(params);
    
    % Run simulations for each R0
    for r_idx = 1:length(params.R0_values)
        R0 = params.R0_values(r_idx);
        simulate_and_analyze_renormalized(params, R0);
    end
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

function simulate_and_analyze_renormalized(params, R0)
    % Calculate beta from R0
    beta = R0 * params.gamma / params.pSI;
    
    % Precompute time vector
    t = 0:params.dt:params.T;
    
    % Store results
    results = cell(length(params.N_values), 1);
    det_result = solve_deterministic_sihrs(beta, params);
    
    % Run simulations for each population size
    for idx = 1:length(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d simulations for N = %d, R0 = %.2f...\n', params.n_runs, N, R0);
        results{idx} = run_multiple_gillespie_renormalized(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results
    plot_renormalized_results(t, results, det_result, params, R0);
end

function result = run_multiple_gillespie_renormalized(N, beta, params, t)
    % Run multiple Gillespie simulations and aggregate results
    S_all = zeros(params.n_runs, length(t));
    I_all = zeros(params.n_runs, length(t));
    H_all = zeros(params.n_runs, length(t));
    R_all = zeros(params.n_runs, length(t));
    D_all = zeros(params.n_runs, length(t));
    
    for run = 1:params.n_runs
        [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params);
        [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);
        S_all(run, :) = S_interp;
        I_all(run, :) = I_interp;
        H_all(run, :) = H_interp;
        R_all(run, :) = R_interp;
        D_all(run, :) = D_interp;
    end
    
    % Compute mean proportions
    result.S_mean = mean(S_all, 1);
    result.I_mean = mean(I_all, 1);
    result.H_mean = mean(H_all, 1);
    result.R_mean = mean(R_all, 1);
    result.D_mean = mean(D_all, 1);
    
    % Compute renormalized infinitesimal standard deviations D_N^(l)(t)
    % D1 (S)
    valid_S = result.S_mean >= 0;
    result.D1_squared(valid_S) = (1/N) * ( ...
        beta * params.pSI * result.I_mean(valid_S) .* result.S_mean(valid_S) + ...
        params.lambda * params.pRS * result.R_mean(valid_S) );
    result.D1 = sqrt(result.D1_squared);
    result.C1 = [0, cumsum(result.D1(1:end-1) .* diff(t))];

    % D2 (I)
    valid_I = result.I_mean >= 0;
    result.D2_squared(valid_I) = (1/N) * ( ...
        (params.gamma * (params.pIH + params.pIR + params.pID)) .* result.I_mean(valid_I) + ...
        beta * params.pSI * result.S_mean(valid_I) .* result.I_mean(valid_I));
    result.D2 = sqrt(result.D2_squared);
    result.C2 = [0, cumsum(result.D2(1:end-1) .* diff(t))];

    % D3 (H)
    valid_H = result.H_mean >= 0;
    result.D3_squared(valid_H) = (1/N) * ( ...
        (params.alpha * (params.pHR + params.pHD)) .* result.H_mean(valid_H) + ...
        params.gamma * params.pIH * result.I_mean(valid_H) );
    result.D3 = sqrt(result.D3_squared);
    result.C3 = [0, cumsum(result.D3(1:end-1) .* diff(t))];

    % D4 (R)
    valid_R = result.R_mean >= 0;
    result.D4_squared(valid_R) = (1/N) * ( ...
        params.lambda * params.pRS * result.R_mean(valid_R) + ...
        params.gamma * params.pIR * result.I_mean(valid_R) + ...
        params.alpha * params.pHR * result.H_mean(valid_R) );
    result.D4 = sqrt(result.D4_squared);
    result.C4 = [0, cumsum(result.D4(1:end-1) .* diff(t))];

    % D5 (D)
    valid_D = result.D_mean >= 0;
    V5 = (1/N) * (params.pID * params.gamma * result.I_mean + params.pHD * params.alpha * result.H_mean);
    result.D5(valid_D) = sqrt(V5(valid_D));
    result.C5 = [0, cumsum(result.D5(1:end-1) .* diff(t))];

    result.N = N;
    result.peak_infected = max(result.I_mean * N);
    result.peak_time = t(find(result.I_mean == max(result.I_mean), 1, 'first'));
    result.final_time = t(end);
end

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, beta, params)
    % Initialize populations
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
    if abs(S + I + H + R + D - N) > 1
        error('Initial populations do not sum to N');
    end
    
    % Preallocate arrays
    max_events = round(10 * params.T * (beta + params.gamma + params.alpha + params.lambda) * N);
    S_hist = zeros(1, max_events);
    I_hist = zeros(1, max_events);
    H_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events);
    D_hist = zeros(1, max_events);
    time_pts = zeros(1, max_events);
    
    S_hist(1) = S; I_hist(1) = I; H_hist(1) = H; R_hist(1) = R; D_hist(1) = D;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    % Gillespie algorithm for SIHRS with death
    while current_time < params.T
        si_rate = (beta / N) * S * I * params.pSI;  % S->I
        ir_rate = params.gamma * I * params.pIR;    % I->R
        ih_rate = params.gamma * I * params.pIH;    % I->H
        id_rate = params.gamma * I * params.pID;    % I->D
        hr_rate = params.alpha * H * params.pHR;    % H->R
        hd_rate = params.alpha * H * params.pHD;    % H->D
        rs_rate = params.lambda * R * params.pRS;   % R->S
        
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
            if S > 0, S = S - 1; I = I + 1; end
        elseif chance < (si_rate + ir_rate)
            if I > 0, I = I - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate)
            if I > 0, I = I - 1; H = H + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            if I > 0, I = I - 1; D = D + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            if H > 0, H = H - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            if H > 0, H = H - 1; D = D + 1; end
        else
            if R > 0, R = R - 1; S = S + 1; end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S; I_hist(event_count) = I; H_hist(event_count) = H;
        R_hist(event_count) = R; D_hist(event_count) = D; time_pts(event_count) = current_time;
    end
    
    % Trim arrays
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    D_hist = D_hist(1:event_count);
    time_pts = time_pts(1:event_count);
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    try
        S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
        I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
        H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
        R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
        D_interp = interp1(time_pts, D_hist, t, 'previous') / N;
        
        % Handle values beyond last event
        S_interp(t > max(time_pts)) = S_hist(end) / N;
        I_interp(t > max(time_pts)) = I_hist(end) / N;
        H_interp(t > max(time_pts)) = H_hist(end) / N;
        R_interp(t > max(time_pts)) = R_hist(end) / N;
        D_interp(t > max(time_pts)) = D_hist(end) / N;
    catch e
        error('Interpolation failed: %s', e.message);
    end
end

function det_result = solve_deterministic_sihrs(beta, params)
    try
        tspan = [0, params.T];
        y0 = [params.initial_s; params.initial_i; params.initial_h; params.initial_r; params.initial_d];
        
        % ODE system
        ode_system = @(t, y) [
            -params.pSI * beta * y(1) * y(2) + params.pRS * params.lambda * y(4);           % ds/dt
             params.pSI * beta * y(1) * y(2) - params.gamma * (1 - params.pII) * y(2);      % di/dt
             params.pIH * params.gamma * y(2) - params.alpha * (1 - params.pHH) * y(3);     % dh/dt
             params.pIR * params.gamma * y(2) + params.pHR * params.alpha * y(3) - params.pRS * params.lambda * y(4); % dr/dt
             params.pID * params.gamma * y(2) + params.pHD * params.alpha * y(3)            % dd/dt
        ];
        
        [T, Y] = ode45(ode_system, tspan, y0);
        
        det_result.T = T;
        det_result.S_prop = Y(:, 1);
        det_result.I_prop = Y(:, 2);
        det_result.H_prop = Y(:, 3);
        det_result.R_prop = Y(:, 4);
        det_result.D_prop = Y(:, 5);
        [det_result.peak_infected_prop, idx] = max(det_result.I_prop);
        det_result.peak_time = T(idx);
        det_result.final_time = T(end);
    catch e
        error('ODE solver failed: %s', e.message);
    end
end

function plot_renormalized_results(t, results, det_result, params, R0)
    % ------------------- D1 (log10) -------------------
 
     

    % =====================================================================
    % NEW: Five independent blocks for ratios Dℓ(t) / (ODE compartment)(t)
    % =====================================================================

    % -------------------- Block 1: D1 / S_ODE --------------------
    FIG_POS = [140, 140, 800, 500];
    LEGEND_FONTSIZE = 18;
    AXIS_FONTSIZE   = 26;
    XLIM_USE = [0, 1000];   % [] to skip
    YLIM_USE = [];          % [] to skip
    LOG_Y    = false;       % true -> log y-axis

    S_det = interp1(det_result.T, det_result.S_prop, t, 'linear', 'extrap');
    tmask = (t >= 1) & (S_det > 0);

    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.C1;
        ratio = nan(size(t));
        ratio(tmask) = Dvals(tmask) ./ S_det(tmask);
        plot(t, ratio, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('D_1(t) / S_{ODE}(t)', 'Interpreter','tex', 'FontSize', 20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize', AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end
    if ~isempty(YLIM_USE), ylim(YLIM_USE); end
    if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false), ...
        %   'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_RCumatio_D1_over_S_R0_%.2f.png', R0));

    % -------------------- Block 2: D2 / I_ODE --------------------
    FIG_POS = [180, 180, 800, 500];
    LEGEND_FONTSIZE = 20;
    AXIS_FONTSIZE   = 28;
    XLIM_USE = [0, 1000];
    YLIM_USE = [];
    LOG_Y    = false;

    I_det = interp1(det_result.T, det_result.I_prop, t, 'linear', 'extrap');
    tmask = (t >= 1) & (I_det > 0);

    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.C2;
        ratio = nan(size(t));
        ratio(tmask) = Dvals(tmask) ./ I_det(tmask);
        plot(t, ratio, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('D_2(t) / I_{ODE}(t)','Interpreter','tex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize', AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end
    if ~isempty(YLIM_USE), ylim(YLIM_USE); end
    if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false), ...
         %  'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D2_over_I_R0_%.2f.png', R0));

    % -------------------- Block 3: D3 / H_ODE --------------------
    FIG_POS = [220, 220, 800, 500];
    LEGEND_FONTSIZE = 18;
    AXIS_FONTSIZE   = 26;
    XLIM_USE = [0, 1000];
    YLIM_USE = [];     % set [] to skip
    LOG_Y    = false;

    H_det = interp1(det_result.T, det_result.H_prop, t, 'linear', 'extrap');
    tmask = (t >= 1) & (H_det > 0);

    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.C3;
        ratio = nan(size(t));
        ratio(tmask) = Dvals(tmask) ./ H_det(tmask);
        plot(t, ratio, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('D_3(t) / H_{ODE}(t)','Interpreter','tex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize', AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end
    if ~isempty(YLIM_USE), ylim(YLIM_USE); end
    if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false), ...
         %  'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D3_over_H_R0_%.2f.png', R0));

    % -------------------- Block 4: D4 / R_ODE --------------------
    FIG_POS = [260, 260, 800, 500];
    LEGEND_FONTSIZE = 16;
    AXIS_FONTSIZE   = 24;
    XLIM_USE = [10, 1000];
    YLIM_USE = [];
    LOG_Y    = false;

    R_det = interp1(det_result.T, det_result.R_prop, t, 'linear', 'extrap');
    tmask = (t >= 1) & (R_det > 0);

    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.C4;
        ratio = nan(size(t));
        ratio(tmask) = Dvals(tmask) ./ R_det(tmask);
        plot(t, ratio, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('D_4(t) / R_{ODE}(t)','Interpreter','tex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize', AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end
    if ~isempty(YLIM_USE), ylim(YLIM_USE); end
    if LOG_Y, set(gca,'YScale','log'); end
    %legend(arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false), ...
        %   'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D4_over_R_R0_%.2f.png', R0));

    % -------------------- Block 5: D5 / D_ODE --------------------
    FIG_POS = [300, 300, 800, 500];
    LEGEND_FONTSIZE = 22;
    AXIS_FONTSIZE   = 28;
    XLIM_USE = [10, 1000];
    YLIM_USE = [];
    LOG_Y    = false;

    D_det = interp1(det_result.T, det_result.D_prop, t, 'linear', 'extrap');
    tmask = (t >= 1) & (D_det > 0);

    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.C5;
        ratio = nan(size(t));
        ratio(tmask) = Dvals(tmask) ./ D_det(tmask);
        plot(t, ratio, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('D_5(t) / D_{ODE}(t)','Interpreter','tex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize', AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end
    if ~isempty(YLIM_USE), ylim(YLIM_USE); end
    if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false), ...
       %    'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D5_over_D_R0_%.2f.png', R0));

    % =====================================================================
    % NEW: Five independent blocks for CUMULATIVE integrals of the ratios
    %      ∫ Dℓ(s) / ODEℓ(s) ds, starting at t>=1 and ODEℓ>0
    % =====================================================================

    % ---- Cum ∫ D1 / S_ODE ----
    FIG_POS = [140, 140, 800, 500];
    LEGEND_FONTSIZE = 18; AXIS_FONTSIZE = 26; XLIM_USE = [0, 1000]; YLIM_USE = []; LOG_Y = false;
    S_det = interp1(det_result.T, det_result.S_prop, t, 'linear', 'extrap');
    valid_den = (t >= 1) & (S_det > 0);
    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.D1; C = nan(size(t));
        if any(valid_den)
            tv = t(valid_den); rv = Dvals(valid_den) ./ S_det(valid_den);
            Cv = cumtrapz(tv, rv); C(valid_den) = Cv;
        end
        plot(t, C, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize',20);
    ylabel('$\int D_1/S_{\mathrm{ODE}}\,ds$','Interpreter','latex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize',AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end; if ~isempty(YLIM_USE), ylim(YLIM_USE); end; if LOG_Y, set(gca,'YScale','log'); end
%    legend(arrayfun(@(n)sprintf('N = %d',n), params.N_values,'UniformOutput',false), 'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D1_over_Sp_R0_%.2f.png', R0));

    % ---- Cum ∫ D2 / I_ODE ----
    FIG_POS = [180, 180, 800, 500];
    LEGEND_FONTSIZE = 20; AXIS_FONTSIZE = 28; XLIM_USE = [0, 1000]; YLIM_USE = []; LOG_Y = false;
    I_det = interp1(det_result.T, det_result.I_prop, t, 'linear', 'extrap');
    valid_den = (t >= 1) & (I_det > 0);
    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.D2; C = nan(size(t));
        if any(valid_den)
            tv = t(valid_den); rv = Dvals(valid_den) ./ I_det(valid_den);
            Cv = cumtrapz(tv, rv); C(valid_den) = Cv;
        end
        plot(t, C, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize',20);
    ylabel('$\int D_2/I_{\mathrm{ODE}}\,ds$','Interpreter','latex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize',AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end; if ~isempty(YLIM_USE), ylim(YLIM_USE); end; if LOG_Y, set(gca,'YScale','log'); end
 %   legend(arrayfun(@(n)sprintf('N = %d',n), params.N_values,'UniformOutput',false), 'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D2_over_I_R0_%.2f.png', R0));

    % ---- Cum ∫ D3 / H_ODE ----
    FIG_POS = [220, 220, 800, 500];
    LEGEND_FONTSIZE = 18; AXIS_FONTSIZE = 26; XLIM_USE = [0, 1000]; YLIM_USE = []; LOG_Y = false;
    H_det = interp1(det_result.T, det_result.H_prop, t, 'linear', 'extrap');
    valid_den = (t >= 1) & (H_det > 0);
    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.D3; C = nan(size(t));
        if any(valid_den)
            tv = t(valid_den); rv = Dvals(valid_den) ./ H_det(valid_den);
            Cv = cumtrapz(tv, rv); C(valid_den) = Cv;
        end
        plot(t, C, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize',20);
    ylabel('$\int D_3/H_{\mathrm{ODE}}\,ds$','Interpreter','latex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize',AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end; if ~isempty(YLIM_USE), ylim(YLIM_USE); end; if LOG_Y, set(gca,'YScale','log'); end
 %   legend(arrayfun(@(n)sprintf('N = %d',n), params.N_values,'UniformOutput',false), 'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D3_over_Hp_R0_%.2f.png', R0));

    % ---- Cum ∫ D4 / R_ODE ----
    FIG_POS = [260, 260, 800, 500];
    LEGEND_FONTSIZE = 16; AXIS_FONTSIZE = 24; XLIM_USE = [10, 1000]; YLIM_USE = []; LOG_Y = false;
    R_det = interp1(det_result.T, det_result.R_prop, t, 'linear', 'extrap');
    valid_den = (t >= 1) & (R_det > 0);
    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.D4; C = nan(size(t));
        if any(valid_den)
            tv = t(valid_den); rv = Dvals(valid_den) ./ R_det(valid_den);
            Cv = cumtrapz(tv, rv); C(valid_den) = Cv;
        end
        plot(t, C, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize',20);
    ylabel('$\int D_4/R_{\mathrm{ODE}}\,ds$','Interpreter','latex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize',AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end; if ~isempty(YLIM_USE), ylim(YLIM_USE); end; if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n)sprintf('N = %d',n), params.N_values,'UniformOutput',false), 'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D4_over_Rp_R0_%.2f.png', R0));

    % ---- Cum ∫ D5 / D_ODE ----
    FIG_POS = [300, 300, 800, 500];
    LEGEND_FONTSIZE = 22; AXIS_FONTSIZE = 28; XLIM_USE = [10, 1000]; YLIM_USE = []; LOG_Y = false;
    D_det = interp1(det_result.T, det_result.D_prop, t, 'linear', 'extrap');
    valid_den = (t >= 1) & (D_det > 0);
    figure('Position', FIG_POS); hold on;
    for idx = 1:length(params.N_values)
        Dvals = results{idx}.D5; C = nan(size(t));
        if any(valid_den)
            tv = t(valid_den); rv = Dvals(valid_den) ./ D_det(valid_den);
            Cv = cumtrapz(tv, rv); C(valid_den) = Cv;
        end
        plot(t, C, 'LineWidth', 3, 'Color', params.colors{idx});
    end
    xlabel('Time (days)','FontSize',20);
    ylabel('$\int D_5/D_{\mathrm{ODE}}\,ds$','Interpreter','latex','FontSize',20);
    grid on; set(gca,'LooseInset',get(gca,'TightInset')); set(gca,'FontSize',AXIS_FONTSIZE);
    if ~isempty(XLIM_USE), xlim(XLIM_USE); end; if ~isempty(YLIM_USE), ylim(YLIM_USE); end; if LOG_Y, set(gca,'YScale','log'); end
  %  legend(arrayfun(@(n)sprintf('N = %d',n), params.N_values,'UniformOutput',false), 'Location','northeast','FontSize',LEGEND_FONTSIZE);
    saveas(gcf, sprintf('SIHRS_CumRatio_D5_over_Dp_R0_%.2f.png', R0));

    % -------------------- Also create a combined plot for reference (optional) --------------------
    % (left as-is or add new combined views below if needed)

end

% Run the simulation
sihrs_renormalized_simulation();
