function SIHRS_multiple_simulations_IPC_vs_ODE()

    % Initialize variables at function level
    N = 300;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s0 = 0.60444078;
    i0 = 0.20422149;
   % h0 = 0.19337719;
   h0=1-s0-i0;
    r0 = 0.0;
    d0 = 0.0;

    % Model parameters
    params = struct(...
        'beta', 1.90875822,      ... % infection rate (β > 0)
        'gamma', 1,     ... % I transition rate (γ > 0) 
        'alpha', 1,       ... % H transition rate (α > 0)
        'lambda', 0.0083,    ... % R transition rate (Λ > 0)  
        'pSI', 0.5,         ... % probability of S to I (p_{SI} in (0,1]) %%%%%%%%%%%%%%%%%%%%%%%% cwang: becareful! 
        'pII', 0.5,        ... % probability of I to I (stay infected)
        'pIH', 0.4,        ... % probability of I to H 
        'pIR', 0.1,       ... % probability of I to R 
        'pID', 0.0,       ... % probability of I to D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cwagn: always zero
        'pHH', 0.5,        ... % probability of H to H (stay hospitalized)
        'pHR', 0.5,      ... % probability of H to R
        'pHD', 0,      ... % probability of H to D    %%%%%%%%%%%%%%%%%%%%%% cWANG: ALWASY MAKE IT ZERO
         'pRR',1,        ... % probability of R to R (stay recovered) Cwang: alwasy make this 1
        'pRS', 0,        ... % probability of R to S    %%%%%%%%%%%%%%%%%%%Cwagng: alway make this 0
        'tmax', 10,        ... % simulation end time )  %%%%%%%%%%%%%%%%%%%%%%%%Cwagn: we need to choose
        's0', s0,           ... % initial susceptible proportion
        'i0', i0,           ... % initial infected proportion
        'h0', h0,           ... % initial hospitalized proportion
        'r0', r0,           ... % initial recovered proportion
        'd0', d0            ... % initial dead proportion
    );




    calculated_R0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated R0 = %.6f \n', calculated_R0);


 calculated_sigma0 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated sigma0 = %.6f \n', calculated_sigma0);

 calculated_sigma1 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated sigma0 = %.6f \n', calculated_sigma1);


  calculated_sigma2 = (params.beta * params.pSI) / params.gamma * (1 - params.pII);
    fprintf('Calculated sigma0 = %.6f \n', calculated_sigma2);



    validate_parameters(params);


    if abs((params.s0 + params.i0 + params.h0 + params.r0 + params.d0) - 1.0) > 1e-10
        error('Initial conditions must sum to 1');
    end

    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_simulations = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if N <= 0
        error('Population size must be positive integer');
    end


    all_results = cell(num_simulations, 1);


    try
        fprintf('Running %d stochastic simulations for N = %d...\n', num_simulations, N);

        for sim_idx = 1:num_simulations
            fprintf('Running simulation %d/%d...\n', sim_idx, num_simulations);
            all_results{sim_idx} = sihrs_agent_model_infected(N, params);
        end

        fprintf('All simulations completed!\n');


        plot_multiple_simulations_infected(all_results, N, params);

    catch ME
        fprintf('Error occurred: %s\n', ME.message);
        rethrow(ME);
    end
end

function validate_parameters(params)
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
    if abs((params.pII + params.pIH + params.pIR + params.pID) - 1.0) > 1e-10
        error('I transition probabilities must sum to 1');
    end
    if abs((params.pHH + params.pHR + params.pHD) - 1.0) > 1e-10
        error('H transition probabilities must sum to 1');
    end
    if abs((params.pRR + params.pRS) - 1.0) > 1e-10
        error('R transition probabilities must sum to 1');
    end
end

function result = sihrs_agent_model_infected(N, params)
    % SIHRS agent-based stochastic model with death (focusing on infected and deaths)
    % Initial conditions
    s0 = round(params.s0 * N);
    i0 = round(params.i0 * N);
    h0 = round(params.h0 * N);
    r0 = round(params.r0 * N);
    d0 = round(params.d0 * N);


    total = s0 + i0 + h0 + r0 + d0;
    if total ~= N

        compartments = [s0, i0, h0, r0, d0];
        [~, largest_idx] = max(compartments);
        compartments(largest_idx) = compartments(largest_idx) + (N - total);
        s0 = compartments(1);
        i0 = compartments(2);
        h0 = compartments(3);
        r0 = compartments(4);
        d0 = compartments(5);
    end


    if (s0 + i0 + h0 + r0 + d0) ~= N
        error('Initial conditions must sum to N');
    end

    max_events = N * 30;
    T = zeros(max_events, 1);
    I_prop = zeros(max_events, 1);
    H_prop = zeros(max_events, 1);
    I_count = zeros(max_events, 1);
    H_count = zeros(max_events, 1);
   


    S = (1:s0)';
    I = ((s0+1):(s0+i0))';
    H = ((s0+i0+1):(s0+i0+h0))';
    R = ((s0+i0+h0+1):(s0+i0+h0+r0))';
    D = ((s0+i0+h0+r0+1):(s0+i0+h0+r0+d0))';


    t = 0.0;
    T(1) = 0.0;
    event_count = 1;


    total_pop = s0 + i0 + h0 + r0 + d0;
    I_prop(1) = i0 / total_pop;
    H_prop(1) = h0 / total_pop;
    I_count(1) = i0;
    H_count(1) = h0;


    while ~isempty(I) && t < params.tmax
        nS = length(S);
        nI = length(I);
        nH = length(H);
        nR = length(R);

        infection_rate = params.pSI * params.beta * nS * nI / N;
        to_susceptible_from_R_rate = params.pRS * params.lambda * nR;
        to_hospital_rate = params.gamma * nI * params.pIH;
        to_recovered_from_I_rate = params.gamma * nI * params.pIR;
        to_dead_from_I_rate = params.gamma * nI * params.pID;
        to_recovered_from_H_rate = params.alpha * nH * params.pHR;
        to_dead_from_H_rate = params.alpha * nH * params.pHD;

        total_rate = infection_rate + to_susceptible_from_R_rate + to_hospital_rate + ...
                     to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate + ...
                     to_dead_from_H_rate;

        if total_rate == 0
            break;
        end


        dt = exprnd(1 / total_rate);
        t = t + dt;

        if t > params.tmax
            t = params.tmax;
            event_count = event_count + 1;
            T(event_count) = t;
            current_total = length(S) + length(I) + length(H) + length(R) + length(D);
            I_prop(event_count) = length(I) / current_total;
            H_prop(event_count) = length(H) / current_total;
            I_count(event_count) = length(I);
            H_count(event_count) = length(H);
            break;
        end

        event_count = event_count + 1;
        T(event_count) = t;


        chance = rand() * total_rate;
        if chance < infection_rate
            % S to I transition
            if nS > 0
                num = randi(nS);
                infected_agent = S(num);
                S(num) = [];
                I = [I; infected_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate)
            % R to S transition
            if nR > 0
                num = randi(nR);
                susceptible_agent = R(num);
                R(num) = [];
                S = [S; susceptible_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate)
            % I to H transition
            if nI > 0
                num = randi(nI);
                hospitalized_agent = I(num);
                I(num) = [];
                H = [H; hospitalized_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate)
            % I to R transition
            if nI > 0
                num = randi(nI);
                recovered_agent = I(num);
                I(num) = [];
                R = [R; recovered_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate)
            % I to D transition
            if nI > 0
                num = randi(nI);
                dead_agent = I(num);
                I(num) = [];
                D = [D; dead_agent];
            end
        elseif chance < (infection_rate + to_susceptible_from_R_rate + to_hospital_rate + to_recovered_from_I_rate + to_dead_from_I_rate + to_recovered_from_H_rate)
            % H to R transition
            if nH > 0
                num = randi(nH);
                recovered_agent = H(num);
                H(num) = [];
                R = [R; recovered_agent];
            end
        else
            % H to D transition
            if nH > 0
                num = randi(nH);
                dead_agent = H(num);
                H(num) = [];
                D = [D; dead_agent];
            end
        end


        current_total = length(S) + length(I) + length(H) + length(R) + length(D);
        I_prop(event_count) = length(I) / current_total;
        H_prop(event_count) = length(H) / current_total;
        I_count(event_count) = length(I);
        H_count(event_count) = length(H);
    end


    T = T(1:event_count);
    I_prop = I_prop(1:event_count);
    H_prop = H_prop(1:event_count);
    I_count = I_count(1:event_count);
    H_count = H_count(1:event_count);


    result = struct(...
        'N', N, ...
        'T', T, ...
        'I_prop', I_prop, ...
        'H_prop', H_prop, ...
        'I_count', I_count, ...
        'H_count', H_count, ...
        'final_time', t, ...
        'peak_infected', max(I_count), ...
        'peak_time', T(argmax(I_count)), ...
        'peak_infected_prop', max(I_prop), ...
        'peak_time_prop', T(argmax(I_prop)), ...
        'i_inf', I_prop(end) ...
    );
end

function plot_multiple_simulations_infected(all_results, N, params)
    t_grid = (0:params.tmax)';
    all_interp_I_prop = zeros(length(all_results), length(t_grid));
    all_interp_H_prop = zeros(length(all_results), length(t_grid));

    for i = 1:length(all_results)
        res = all_results{i};

        if length(res.T) > 1
            itp_I = griddedInterpolant(res.T, res.I_prop, 'linear', 'none');
            itp_H = griddedInterpolant(res.T, res.H_prop, 'linear', 'none');

            for j = 1:length(t_grid)
                if t_grid(j) >= min(res.T) && t_grid(j) <= max(res.T)
                    all_interp_I_prop(i, j) = itp_I(t_grid(j));
                    all_interp_H_prop(i, j) = itp_H(t_grid(j));
                end
            end
        end
    end


    %%%%% Figure of Infected proportions %%%%%%%
    figure;

    for i = 1:size(all_interp_I_prop, 1)
        plot(t_grid, all_interp_I_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
        hold on;
    end

    det_result = solve_deterministic_sihrs(params);
    det_I_prop = det_result.I_prop;

    i_det = plot(det_result.T, det_I_prop, 'r-', 'LineWidth', 2.5); %%%% use det here instead of real %%%
    xlabel('Time (days)');
    ylabel('Infected Proportion');
    title('IPC simulations vs ODE');
    xlim([0, params.tmax]);
    ylim([0, max([max(max(all_interp_I_prop)), max(det_I_prop)]) * 1.1]);


    % Custom legend
    stoch=plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
    legend([stoch, i_det],{'Stochastic Simulations', 'ODE Data'});
    saveas(gcf, 'SIHRS_infected_ODE_vs_IPC.png');


    %%%%%%% Figure of Hospitalized proportion%%%%%%%%%%
    figure;

    for i = 1:size(all_interp_H_prop, 1)
        plot(t_grid, all_interp_H_prop(i, :), 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
        hold on;
    end

    det_result = solve_deterministic_sihrs(params);
    det_H_prop = det_result.H_prop;

    h_det=plot(det_result.T, det_H_prop, 'g-', 'LineWidth', 2.5); %%%% use det here instead of real %%%
    xlabel('Time (days)');
    ylabel('Hospitalized Proportion');
    title('IPC simulations vs ODE');
    xlim([0, params.tmax]);
    ylim([0, max([max(max(all_interp_H_prop)), max(det_H_prop)]) * 1.1]);


    % Custom legend
    h_stoch=plot(NaN, NaN, 'Color', [0.2, 0.4, 0.8, 0.3], 'LineWidth', 1.0);
    legend([h_stoch, h_det],{'Stochastic Simulations', 'ODE Data'});
    saveas(gcf, 'SIHRS_hospitalized_ODE_vs_IPC.png');

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

function idx = argmax(x)
    [~, idx] = max(x);
end
