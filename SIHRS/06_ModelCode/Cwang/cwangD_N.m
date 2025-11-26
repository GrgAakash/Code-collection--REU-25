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
    params.pIH = 0.1060;           % probability of I to H
    params.pIR = 0.8921;           % probability of I to R
    params.pID = 0.0019;           % probability of I to D
    params.pHH = 0.00;             % probability of H to H (stay hospitalized)
    params.pHR = 0.846;            % probability of H to R
    params.pHD = 0.154;            % probability of H to D
    params.pRR = 0.02;             % probability of R to R (stay recovered)
    params.pRS = 0.98;             % probability of R to S
    params.gamma = 0.1;            % Infection transition rate (γ > 0)
    params.alpha = 0.1;            % Hospitalized transition rate (α > 0)
    params.lambda = 0.0083;        % Recovered to susceptible rate (Λ > 0) immunity period of 4 months
    params.T = 1000;               % Total simulation time
    params.dt = 0.01;              % Time step for interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pattern formation plot. only choose the I and H, and cummulative, all
%parameters, small and large, sweep
   %%%%%%%%%%%%%%%%
   
    params.N_values = [350, 600, 1600, 3000, 6000, 10000]; % Population sizes - matching SIHRS.m
    
  % params.N_values = [300, 1600, 10000];

 % params.N_values = [300, 1600, 3000];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%params.N_values = [250, 350, 800, 1600, 3000];

%params.N_values = [250, 350, 450];
% params.N_values = [350, 800, 1600, 3000, 10000];
 %params.N_values = [800, 1600, 3000];

%params.N_values = [1600, 3000, 6000, 10000];
%Cwang:  you only need to delete all the instances of params.colors, to
%allow infinitely many Ns you want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right after you set params.N_values:
params.colors = num2cell(lines(numel(params.N_values)), 2);   % each cell is [r g b]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   params.initial_s = 0.96;       % Initial susceptible fraction
    params.initial_i = 0.035;       % Initial infected fraction
    params.initial_h = 0.005;          % Initial hospitalized fraction
    params.initial_r = 0;          % Initial recovered fraction
    params.initial_d = 0;          % Initial dead fraction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    params.n_runs = 20;            % Number of stochastic runs
    %params.colors = {'#0072BD', '#77AC30', '#A2142F'}; % Colors matching SIHRS.m
    params.R0_values = [2.12]; % R0 values
    
    % Validate parameters
    validate_params(params);
    
    % Run simulations for each R0
    for r_idx = 1:length(params.R0_values)
        R0 = params.R0_values(r_idx);
        simulate_and_analyze_renormalized(params, R0);
    end
end

function validate_params(params)
    % Validate input parameters
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
    
    % Check transition probability sums
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
        
        % Run multiple stochastic simulations
        results{idx} = run_multiple_gillespie_renormalized(N, beta, params, t);
        fprintf('Completed N = %d\n', N);
    end
    
    % Plot results and print summary
    plot_renormalized_results(t, results, det_result, params, R0);
%    print_renormalized_summary(results, det_result, R0);
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
    % Using the mathematical formulas provided
    % Note: We ensure denominators are non-zero to avoid division by zero
    
    % D_N^(1)(t) for Susceptible (S)
    % (D_N^(1))^2 = (1/N) * (β_SI * i_N(t)/s_N(t) + λ_RS * r_N(t)/(s_N(t))^2)
    % Check if s_N(t) is non-zero
    valid_S = result.S_mean > 1e-10;  % Threshold for non-zero values
    result.D1_squared = zeros(size(result.S_mean));
    result.D1_squared(valid_S) = (1/N) * (beta * params.pSI * result.I_mean(valid_S) ./ result.S_mean(valid_S) + ...
                                    params.lambda * params.pRS * result.R_mean(valid_S) ./ (result.S_mean(valid_S).^2));
    result.D1 = sqrt(result.D1_squared);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>> ADD: cumulative D1 over time (integral from 0 to t)
  %  result.C1 = cumtrapz(t, result.D1);         % cumulative integral
  result.C1 = [0, cumsum(result.D1(1:end-1) .* diff(t))];
 

%chg = [true, diff(result.D1)~=0];  result.C1 = cumsum(result.D1 .* chg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % D_N^(2)(t) for Infected (I)
    % (D_N^(2))^2 = (1/N) * ((γ_IH + γ_IR + γ_ID)/i_N(t) + β_SI * s_N(t)/i_N(t))
    % Check if i_N(t) is non-zero
    valid_I = result.I_mean > 1e-10;  % Threshold for non-zero values
    result.D2_squared = zeros(size(result.I_mean));
    result.D2_squared(valid_I) = (1/N) * ((params.gamma * (params.pIH + params.pIR + params.pID)) ./ result.I_mean(valid_I) + ...
                                    beta * params.pSI * result.S_mean(valid_I) ./ result.I_mean(valid_I));
    result.D2 = sqrt(result.D2_squared);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 result.C2 = [0, cumsum(result.D2(1:end-1) .* diff(t))];
 
% chg = [true, diff(result.D2)~=0];  result.C2 = cumsum(result.D2 .* chg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % D_N^(3)(t) for Hospitalized (H)
    % (D_N^(3))^2 = (1/N) * ((α_HR + α_HD)/h_N(t) + γ_IH * i_N(t)/(h_N(t))^2)
    % Check if h_N(t) is non-zero
    valid_H = result.H_mean > 1e-10;  % Threshold for non-zero values
    result.D3_squared = zeros(size(result.H_mean));
    result.D3_squared(valid_H) = (1/N) * ((params.alpha * (params.pHR + params.pHD)) ./ result.H_mean(valid_H) + ...
                                    params.gamma * params.pIH * result.I_mean(valid_H) ./ (result.H_mean(valid_H).^2));
    result.D3 = sqrt(result.D3_squared);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 result.C3 = [0, cumsum(result.D3(1:end-1) .* diff(t))];
% result.C3 = cumulative_step_sum(result.D3);      % <-- jump-sum, no dt
 
%chg = [true, diff(result.D3)~=0];  result.C3 = cumsum(result.D3 .* chg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % D_N^(4)(t) for Recovered (R)
    % (D_N^(4))^2 = (1/N) * (λ_RS/r_N(t) + γ_IR * i_N(t)/(r_N(t))^2 + α_HR * h_N(t)/(r_N(t))^2)
    % Check if r_N(t) is non-zero
    valid_R = result.R_mean > 1e-10;  % Threshold for non-zero values
    result.D4_squared = zeros(size(result.R_mean));
    result.D4_squared(valid_R) = (1/N) * (params.lambda * params.pRS ./ result.R_mean(valid_R) + ...
                                    params.gamma * params.pIR * result.I_mean(valid_R) ./ (result.R_mean(valid_R).^2) + ...
                                    params.alpha * params.pHR * result.H_mean(valid_R) ./ (result.R_mean(valid_R).^2));
    result.D4 = sqrt(result.D4_squared);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 result.C4 = [0, cumsum(result.D4(1:end-1) .* diff(t))];
% result.C4 = cumulative_step_sum(result.D4);      % <-- jump-sum, no dt
%chg = [true, diff(result.D4)~=0];  result.C4 = cumsum(result.D4 .* chg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % D_N^(5)(t) for Dead (D) - not in the mathematical formulas but included for completeness
    % Using the variance formula from the original code
    % Check if d_N(t) is non-zero
    valid_D = result.D_mean > 1e-10;  % Threshold for non-zero values
    V5 = (1/N) * (params.pID * params.gamma * result.I_mean + params.pHD * params.alpha * result.H_mean);
    result.D5 = zeros(size(result.D_mean));
    result.D5(valid_D) = sqrt(V5(valid_D)) ./ result.D_mean(valid_D);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 result.C5 = [0, cumsum(result.D5(1:end-1) .* diff(t))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % result.C5 = cumulative_step_sum(result.D5); 

 %chg = [true, diff(result.D5)~=0];  result.C5 = cumsum(result.D5 .* chg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % Verify population conservation
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
    
    S_hist(1) = S;
    I_hist(1) = I;
    H_hist(1) = H;
    R_hist(1) = R;
    D_hist(1) = D;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    % Gillespie algorithm for SIHRS with death
    while current_time < params.T
        % Calculate event rates
        si_rate = (beta / N) * S * I * params.pSI; % rate of S->I
        ir_rate = params.gamma * I * params.pIR; % rate of I->R
        ih_rate = params.gamma * I * params.pIH; % rate of I->H
        id_rate = params.gamma * I * params.pID; % rate of I->D
        hr_rate = params.alpha * H * params.pHR; % rate of H->R
        hd_rate = params.alpha * H * params.pHD; % rate of H->D
        rs_rate = params.lambda * R * params.pRS; % rate of R->S
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        
        if total_rate == 0
            break;
        end
        
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        
        if current_time > params.T
            break;
        end
        
        % Determine which event occurs
        chance = rand * total_rate;
        if chance < si_rate
            % S->I transition
            if S > 0
                S = S - 1;
                I = I + 1;
            end
        elseif chance < (si_rate + ir_rate)
            % I->R transition
            if I > 0
                I = I - 1;
                R = R + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate)
            % I->H transition
            if I > 0
                I = I - 1;
                H = H + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            % I->D transition
            if I > 0
                I = I - 1;
                D = D + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            % H->R transition
            if H > 0
                H = H - 1;
                R = R + 1;
            end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            % H->D transition
            if H > 0
                H = H - 1;
                D = D + 1;
            end
        else
            % R->S transition
            if R > 0
                R = R - 1;
                S = S + 1;
            end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S;
        I_hist(event_count) = I;
        H_hist(event_count) = H;
        R_hist(event_count) = R;
        D_hist(event_count) = D;
        time_pts(event_count) = current_time;
    end
    
    % Trim arrays
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    D_hist = D_hist(1:event_count);
    time_pts = time_pts(1:event_count);
    
    % Verify population conservation
    assert(all(S_hist + I_hist + H_hist + R_hist + D_hist == N), 'Population not conserved');
end

function [S_interp, I_interp, H_interp, R_interp, D_interp] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
    % Interpolate to fixed time grid using 'previous'
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
    % Solve deterministic SIHRS with death model using ODE45
    try
        tspan = [0, params.T];
        y0 = [params.initial_s; params.initial_i; params.initial_h; params.initial_r; params.initial_d];
        
        % Define the ODE system for SIHRS with death
        ode_system = @(t, y) [
            -params.pSI * beta * y(1) * y(2) + params.pRS * params.lambda * y(4);           % ds/dt
            params.pSI * beta * y(1) * y(2) - params.gamma * (1 - params.pII) * y(2);      % di/dt
            params.pIH * params.gamma * y(2) - params.alpha * (1 - params.pHH) * y(3);     % dh/dt
            params.pIR * params.gamma * y(2) + params.pHR * params.alpha * y(3) - params.pRS * params.lambda * y(4); % dr/dt
            params.pID * params.gamma * y(2) + params.pHD * params.alpha * y(3)             % dd/dt
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
   
% Create 5 separate figures for each compartment showing renormalized standard deviations
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot D1 (renormalized standard deviation for susceptible)
    figure('Position', [100, 100, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D1,  'LineWidth', 3);
    end
   % title(sprintf('Renormalized Std Dev - Susceptible (R_0 = %.2f)', R0));
    xlabel('Time (days)','FontSize', 20);
    ylabel('$\mathcal{D}_N^{(1)}(t)$', 'Interpreter','latex', 'FontSize', 20);
    grid on;
    
 %   xlim([0.1, 1000]);


set(gca,'LooseInset',get(gca,'TightInset'));

set(gca,'FontSize',30);
     %'FontName', 'Courier');

    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);

     

    % Save the figure
    saveas(gcf, sprintf('SIHRS_D1_Susceptible_R0_%.2f.png', R0));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NEW: Plot cumulative D1 (integral of D1 over time)
    figure('Position', [150, 150, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.C1,   'LineWidth', 3);
    end
    xlabel('Time (days)','FontSize', 20);
    ylabel('$\int_0^t \mathcal{D}_N^{(1)}(s)\,ds$', 'Interpreter','latex', 'FontSize', 20);
    grid on;

    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gca,'FontSize',30);

    legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);

    % Save the figure
    saveas(gcf, sprintf('SIHRS_CumD1_Susceptible_R0_%.2f.png', R0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: Cumulative D2 (integral of D2 over time)
figure('Position', [220, 220, 800, 500]);
hold on;
for idx = 1:length(params.N_values)
    plot(t, results{idx}.C2,   'LineWidth', 3);
end
xlabel('Time (days)','FontSize', 30);
ylabel('$\int_0^t \mathcal{D}_N^{(2)}(s)\,ds$', 'Interpreter','latex', 'FontSize', 20);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'FontSize',30);
legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
saveas(gcf, sprintf('SIHRS_CumD2_Infected_R0_%.2f.png', R0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot D2 (renormalized standard deviation for infected)
    figure('Position', [200, 200, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D2,   'LineWidth', 3);
    end
   % title(sprintf('Renormalized Std Dev - Infected (R_0 = %.2f)', R0));
xlabel('Time (days)','FontSize', 30);
    ylabel('$\mathcal{D}_N^{(2)}(t)$', 'Interpreter','latex', 'FontSize', 20);
    grid on;

%xlim([0.1, 1000]);


    set(gca,'LooseInset',get(gca,'TightInset'));
     
set(gca,'FontSize',30);
     %'FontName', 'Courier');



       legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
  



    % Save the figure
    saveas(gcf, sprintf('SIHRS_D2_Infected_R0_%.2f.png', R0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot D3 (renormalized standard deviation for hospitalized)
    figure('Position', [300, 300, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, results{idx}.D3,   'LineWidth', 4);
    end
  %  title(sprintf('Renormalized Std Dev - Hospitalized (R_0 = %.2f)', R0));
     
    
    
    
    xlabel('Time (days)','FontSize', 30);
    ylabel('$\mathcal{D}_N^{(3)}(t)$', 'Interpreter','latex', 'FontSize', 20);
    grid on;


 % Set y-axis maximum to 1
  %  ylim([0, 0.6]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 xlim([10, 1000]);
%xlim([0, 1000]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontSize',30);
     %'FontName', 'Courier');

    set(gca,'LooseInset',get(gca,'TightInset'));
    % axis square;

       legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
  
    
    
    % Save the figure
   
   saveas(gcf, sprintf('SIHRS_D3_Hospitalized_R0_%.2f.png', R0));


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: Cumulative D3 (integral of D3 over time)
figure('Position', [320, 320, 800, 500]);
hold on;
for idx = 1:length(params.N_values)
    plot(t, results{idx}.C3,   'LineWidth', 3);
end
xlabel('Time (days)','FontSize', 20);
ylabel('$\int_0^t \mathcal{D}_N^{(3)}(s)\,ds$', 'Interpreter','latex', 'FontSize', 20);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 xlim([10, 1000]);  % match your D3 x-limits if desired

%xlim([0, 1000]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'FontSize',30);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
saveas(gcf, sprintf('SIHRS_CumD3_Hospitalized_R0_%.2f.png', R0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Plot D4 (renormalized standard deviation for recovered)
    figure('Position', [400, 400, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, log10(results{idx}.D4),   'LineWidth', 1.5);
    end
    %title(sprintf('Renormalized Std Dev - Recovered (R_0 = %.2f)', R0));
  
   
      xlabel('Time (days)','FontSize', 20);
    ylabel('$\log{\left(\mathcal{D}_N^{(4)}(t)\right)}$', 'Interpreter','latex', 'FontSize', 20);
    grid on;


 % Set y-axis maximum to 0.07
  %  ylim([0, 0.07]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([10, 1000]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gca,'LooseInset',get(gca,'TightInset'));


    % axis square;

 set(gca,'FontSize',20);
     %'FontName', 'Courier');

       legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
  
    
    
    
    
    
    % Save the figure
    saveas(gcf, sprintf('SIHRS_D4_Recoveredlog_R0_%.2f.png', R0));
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: Cumulative D4 (integral of D4 over time)
figure('Position', [420, 420, 800, 500]);
hold on;
for idx = 1:length(params.N_values)
    plot(t, results{idx}.C4,   'LineWidth', 3);
end
xlabel('Time (days)','FontSize', 10);
ylabel('$\int_0^t \mathcal{D}_N^{(4)}(s)\,ds$', 'Interpreter','latex', 'FontSize', 10);
grid on;
xlim([10, 1000]);  % match your D4 x-limits if desired
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'FontSize',20);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 15);
saveas(gcf, sprintf('SIHRS_CumD4_Recovered_R0_%.2f.png', R0));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot D5 (renormalized standard deviation for dead)
    figure('Position', [500, 500, 800, 500]);
    hold on;
    for idx = 1:length(params.N_values)
        plot(t, log10(results{idx}.D5),   'LineWidth', 1.5);
    end
   % title(sprintf('Renormalized Std Dev - Dead (R_0 = %.2f)', R0));
    
     
    
    
    
    xlabel('Time (days)','FontSize', 30);
    ylabel('$\log{\left(\mathcal{D}_N^{(5)}(t)\right)}$', 'Interpreter','latex', 'FontSize', 10);
    grid on;


 % Set y-axis maximum to 1
 %   ylim([0, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([10, 1000]);

    set(gca,'LooseInset',get(gca,'TightInset'))
    % axis square;

       legend_labels = arrayfun(@(n) sprintf('N = %d', n), params.N_values, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'northwest', 'FontSize', 20);
  
    
    
    
    % Save the figure
    saveas(gcf, sprintf('SIHRS_D5_Deadlog_R0_%.2f.png', R0));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW: Cumulative D5 (integral of D5 over time)
figure('Position', [520, 520, 800, 500]);
hold on;
for idx = 1:length(params.N_values)
    plot(t, results{idx}.C5, 'Color', params.colors{idx}, 'LineWidth', 3);
end
xlabel('Time (days)','FontSize', 10);
ylabel('$\int_0^t \mathcal{D}_N^{(5)}(s)\,ds$', 'Interpreter','latex', 'FontSize', 10);
grid on;
xlim([10, 1000]);  % match your D5 x-limits if desired
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'FontSize',20);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 15);
saveas(gcf, sprintf('SIHRS_CumD5_Dead_R0_%.2f.png', R0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Also create a combined plot for reference
 
end

 

% Run the simulation
sihrs_renormalized_simulation();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
