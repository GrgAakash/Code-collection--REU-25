% SIHRS Model Implementation
% Based on the mathematical model:
% s' = -β*pSI*s*i + pRS*Λ*r
% i' = β*pSI*s*i - γ*(1-pII)*i  
% h' = pIH*γ*i - α*(1-pHH)*h
% r' = pIR*γ*i + pHR*α*h - pRS*Λ*r
% d' = pID*γ*i + pHD*α*h

function sihrs_model()
    % Initial conditions
    s0 = 0.96;
    i0 = 0.04;
    h0 = 0.0;
    r0 = 0.0;
    d0 = 0.0;
    
    % Population sizes to test
    Ns = [316, 1000, 3162, 10000];
    
    % Store data for all simulations
    all_data.ns = {}; all_data.ni = {}; all_data.nh = {}; 
    all_data.nr = {}; all_data.nd = {};
    all_data.vs = {}; all_data.vi = {}; all_data.vh = {}; 
    all_data.vr = {}; all_data.vd = {};
    all_data.T = {};

    % Model parameters (SIHRS)
    constants.max_t = 500;
    constants.beta = 0.5855;        % infection rate
    constants.gamma = 0.0714;     % recovery rate
    constants.alpha = 0.17;      % hospital recovery rate
    constants.lambda = 0.0012;      % immunity loss rate
    
    % Probabilities (must sum to 1 for each transition)
    constants.pSI = 1.0;          % S to I
    constants.pII = 0.00;         % I to I (stay infected)
    constants.pIH = 0.015;         % I to H
    constants.pIR = 0.984;         % I to R
    constants.pID = 0.001;         % I to D
    constants.pHH = 0.00;         % H to H (stay hospitalized)
    constants.pHR = 0.97;         % H to R
    constants.pHD = 0.03;         % H to D
    constants.pRR = 0.05;          % R to R (stay recovered)
    constants.pRS = 0.95;          % R to S
    
    % Validate probability constraints
    validate_probabilities(constants);

    %% Loop for each N
    for i = 1:numel(Ns)
        constants.N = Ns(i);
        I0 = round(i0 * constants.N);
        H0 = round(h0 * constants.N);
        R0 = round(r0 * constants.N);
        D0 = round(d0 * constants.N);
        S0 = constants.N - I0 - H0 - R0 - D0; % ensure sum equals N

        % Arrays of agents
        S = 1:S0;
        I = (S0+1):(S0+I0);
        H = (S0+I0+1):(S0+I0+H0);
        R = (S0+I0+H0+1):(S0+I0+H0+R0);
        D = (S0+I0+H0+R0+1):constants.N;

        % Initialize data structure
        data.ns = S0;
        data.ni = I0;
        data.nh = H0;
        data.nr = R0;
        data.nd = D0;
        data.vs = 0;
        data.vi = 0;
        data.vh = 0;
        data.vr = 0;
        data.vd = 0;
        data.T = 0;

        %% Perform Simulations
        data = simulate_sihrs(S, I, H, R, D, constants, data);
        
        % Store results
        all_data.ns{end+1} = data.ns;
        all_data.ni{end+1} = data.ni;
        all_data.nh{end+1} = data.nh;
        all_data.nr{end+1} = data.nr;
        all_data.nd{end+1} = data.nd;
        all_data.vs{end+1} = data.vs;
        all_data.vi{end+1} = data.vi;
        all_data.vh{end+1} = data.vh;
        all_data.vr{end+1} = data.vr;
        all_data.vd{end+1} = data.vd;
        all_data.T{end+1} = data.T;

        %% Graphing
        plot_proportions_vs_time(data, constants.N);
    end

    %% Deterministic Model (ODE solution)
    y0 = [s0; i0; h0; r0; d0]; 
    tspan = [0 constants.max_t];
    [t_det, y_det] = ode45(@(t, y) sihrs_system(t, y, constants), tspan, y0);
    
    % Plot comparisons
    plot_compartment_vs_time(all_data.T, all_data.ns, t_det, y_det(:,1), Ns, "Susceptible Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.ni, t_det, y_det(:,2), Ns, "Infected Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nh, t_det, y_det(:,3), Ns, "Hospitalized Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nr, t_det, y_det(:,4), Ns, "Recovered Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nd, t_det, y_det(:,5), Ns, "Dead Proportion", constants.max_t);
    
    % Calculate and display R0
    R0 = (constants.pSI * constants.beta) / (constants.gamma * (1 - constants.pII));
    fprintf('R0 = %.2f\n', R0);
    
    % Display asymptotic values
    fprintf('Asymptotic Values (t=%.0f):\n', constants.max_t);
    fprintf('s(∞) = %.4f\n', y_det(end,1));
    fprintf('i(∞) = %.4f\n', y_det(end,2));
    fprintf('h(∞) = %.4f\n', y_det(end,3));
    fprintf('r(∞) = %.4f\n', y_det(end,4));
    fprintf('d(∞) = %.4f\n', y_det(end,5));
end

function validate_probabilities(constants)
    % Validate probability constraints
    pII_sum = constants.pII + constants.pIH + constants.pIR + constants.pID;
    pHH_sum = constants.pHH + constants.pHR + constants.pHD;
    pRR_sum = constants.pRR + constants.pRS;
    
    if abs(pII_sum - 1) > 1e-10
        error('I transition probabilities must sum to 1, got %.6f', pII_sum);
    end
    if abs(pHH_sum - 1) > 1e-10
        error('H transition probabilities must sum to 1, got %.6f', pHH_sum);
    end
    if abs(pRR_sum - 1) > 1e-10
        error('R transition probabilities must sum to 1, got %.6f', pRR_sum);
    end
end

function data = simulate_sihrs(S, I, H, R, D, constants, data) 
    % Add safety check to prevent infinite recursion
    if length(data.T) > 100000  % Safety limit
        fprintf('Simulation stopped: too many events\n');
        return;
    end
    
    ns = numel(S);
    ni = numel(I);
    nh = numel(H);
    nr = numel(R);
    nd = numel(D);

    % Calculate event rates according to SIHRS model
    infection_rate = constants.pSI * constants.beta * ns * ni / constants.N;
    recovery_rate = constants.gamma * ni * (1 - constants.pII);
    hospital_rate = constants.alpha * nh * (1 - constants.pHH);
    immunity_loss_rate = constants.pRS * constants.lambda * nr;
    
    master_clock = infection_rate + recovery_rate + hospital_rate + immunity_loss_rate;

    % Check termination conditions
    if data.T(end) >= constants.max_t || master_clock == 0 || (ni == 0 && nh == 0)
        if data.T(end) < constants.max_t
            data.ns(end+1) = data.ns(end);
            data.ni(end+1) = data.ni(end);
            data.nh(end+1) = data.nh(end);
            data.nr(end+1) = data.nr(end);
            data.nd(end+1) = data.nd(end);
            data.vs(end+1) = data.vs(end);
            data.vi(end+1) = data.vi(end);
            data.vh(end+1) = data.vh(end);
            data.vr(end+1) = data.vr(end);
            data.vd(end+1) = data.vd(end);
            data.T(end+1) = constants.max_t;
        end
        return;
    end

    % Record current state
    data.ns(end+1) = ns; 
    data.ni(end+1) = ni; 
    data.nh(end+1) = nh; 
    data.nr(end+1) = nr; 
    data.nd(end+1) = nd;
    
    % Calculate variances (simplified)
    data.vs(end+1) = infection_rate / constants.N;
    data.vi(end+1) = recovery_rate / constants.N;
    data.vh(end+1) = hospital_rate / constants.N;
    data.vr(end+1) = immunity_loss_rate / constants.N;
    data.vd(end+1) = (constants.pID * constants.gamma * ni + constants.pHD * constants.alpha * nh) / constants.N;

    % Time to next event
    dt = exprnd(1 / master_clock);
    data.T(end+1) = data.T(end) + dt;

    % Determine which event occurs
    r = rand;
    if (r < infection_rate / master_clock)                    % S to I
        if ns > 0
            I(end+1) = S(end);
            S(end) = [];
        end
    elseif (r < (infection_rate + recovery_rate) / master_clock)  % I transitions
        if ni > 0
            r2 = rand;
            if r2 < constants.pIH                              % I to H
                H(end+1) = I(end); 
                I(end) = [];
            elseif r2 < (constants.pIH + constants.pIR)        % I to R
                R(end+1) = I(end);
                I(end) = [];
            elseif r2 < (constants.pIH + constants.pIR + constants.pID) % I to D
                D(end+1) = I(end);
                I(end) = [];
            else                                               % I to I (stay)
                % No change
            end
        end
    elseif (r < (infection_rate + recovery_rate + hospital_rate) / master_clock) % H transitions
        if nh > 0
            r2 = rand;
            if r2 < constants.pHR                              % H to R
                R(end+1) = H(end);
                H(end) = [];
            elseif r2 < (constants.pHR + constants.pHD)        % H to D
                D(end+1) = H(end);
                H(end) = [];
            else                                               % H to H (stay)
                % No change
            end
        end
    else                                                       % R to S
        if nr > 0
            S(end+1) = R(end);
            R(end) = [];
        end
    end

    % Recursive call with updated state
    data = simulate_sihrs(S, I, H, R, D, constants, data);
end

function dydt = sihrs_system(~, y, constants)
    % SIHRS ODE system
    s = y(1);
    i = y(2);
    h = y(3);
    r = y(4);
    d = y(5);

    % ODEs according to the mathematical model
    ds = -constants.beta * constants.pSI * s * i + constants.pRS * constants.lambda * r;
    di = constants.beta * constants.pSI * s * i - constants.gamma * (1 - constants.pII) * i;
    dh = constants.pIH * constants.gamma * i - constants.alpha * (1 - constants.pHH) * h;
    dr = constants.pIR * constants.gamma * i + constants.pHR * constants.alpha * h - constants.pRS * constants.lambda * r;
    dd = constants.pID * constants.gamma * i + constants.pHD * constants.alpha * h;

    dydt = [ds; di; dh; dr; dd];
end

function plot_proportions_vs_time(data, N) 
    figure;
    hold on;
    plot(data.T, data.ns / N, 'b', 'DisplayName', 'Susceptible', 'LineWidth', 1.5);
    plot(data.T, data.ni / N, 'r', 'DisplayName', 'Infected', 'LineWidth', 1.5);
    plot(data.T, data.nh / N, 'm', 'DisplayName', 'Hospitalized', 'LineWidth', 1.5);
    plot(data.T, data.nr / N, 'g', 'DisplayName', 'Recovered', 'LineWidth', 1.5);
    plot(data.T, data.nd / N, 'k', 'DisplayName', 'Dead', 'LineWidth', 1.5);
    ylabel('Proportion');
    legend('Location', 'eastoutside');
    title(['SIHRS Model - All Proportions vs Time (N = ', num2str(N), ')']);
    xlabel('Time');
    grid on;
    hold off;
end

function plot_compartment_vs_time(sim_t, sim_data, theo_t, theo_data, Ns, compartment, max_t)
    figure;
    hold on;
    colors = lines(length(Ns));
    legend_labels = cell(1, length(Ns)+1);

    for i = 1:length(Ns)
        plot(sim_t{i}, sim_data{i} / Ns(i), 'Color', colors(i,:), 'LineWidth', 1.5);
        legend_labels{i} = ['N = ', num2str(Ns(i))];
    end
    
    if ~isempty(theo_data)
        plot(theo_t, theo_data, '--', 'Color', 'k', 'LineWidth', 2);
        legend_labels{length(Ns)+1} = 'ODE';
    end

    ylabel(compartment);
    title([compartment, ' vs Time']);
    xlabel('Time');
    legend(legend_labels);
    xlim([0 max_t]);
    grid on;
    hold off;
end 