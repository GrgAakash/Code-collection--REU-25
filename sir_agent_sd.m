function sir_agent_sd()
    % Model Parameters
    beta = 0.95;
    gamma = 1;
    p1 = 0.5;
    p2 = 0.5;
    max_T = 23;
    t_grid = linspace(0, max_T, 1000);

    %% Deterministic (ODE) Model
    s0 = 0.96;
    i0 = 0.04;
    r0 = 0;
    y0 = [s0; i0; r0];

    [t_ode, y_ode] = ode45(@(t, y) sir_system(t, y, beta, gamma, p1, p2), [0 max_T], y0);
    s_ode = y_ode(:,1);
    i_ode = y_ode(:,2);
    r_ode = y_ode(:,3);
    i_ode_interp = interp1(t_ode, i_ode, t_grid, 'linear', 'extrap');

    %% Agent-Based Simulations
    Ns = [316, 1000, 3162, 10000];
    I_all = zeros(length(Ns), length(t_grid));
    S_var_all = zeros(length(Ns), length(t_grid));
    I_var_all = zeros(length(Ns), length(t_grid));
    R_var_all = zeros(length(Ns), length(t_grid));
    colors = lines(length(Ns));
    
    I_count_all_316 = zeros(1, length(t_grid));  % Only track for N = 316

    for idx = 1:length(Ns)
        N = Ns(idx);
        s0_abs = round(0.96 * N);
        i0_abs = round(0.04 * N);
        r0_abs = 0;

        S = 1:s0_abs;
        I = (s0_abs+1):(s0_abs+i0_abs);
        R = [];

        t = 0;
        T = 0;
        S_prop = s0_abs / N;
        I_prop = i0_abs / N;
        R_prop = r0_abs / N;
        S_var = 0;
        I_var = 0;
        R_var = 0;

        I_count = numel(I);  % initial infected count
        T = 0;
        I_count_list = I_count;  % store history of I counts


        while ~isempty(I) && t < max_T
            nI = numel(I);
            nS = numel(S);

            infection_rate = beta * nS * nI / (N^1.0); % density-dependent
            recovery_rate = gamma * nI;
            event_rate = infection_rate + recovery_rate;

            if event_rate == 0
                break;
            end

            dt = exprnd(1 / event_rate);
            t = t + dt;

            S_var(end+1) = p1 * beta * (nI/N) * (nS/N) / N;
            I_var(end+1) = (p1 * beta * (nI/N) * (nS/N) + p2 * gamma * (nI/N)) / N;
            R_var(end+1) = p2 * gamma * (nI/N) / N;

            if rand < (infection_rate / event_rate)
                if nS > 0 && rand < p1
                    num = randi(nS);
                    infected = S(num);
                    S(num) = [];
                    I(end+1) = infected;
                end
            else
                if rand < p2
                    num = randi(nI);
                    recovered = I(num);
                    I(num) = [];
                    R(end+1) = recovered;
                end
            end

            N_current = numel(S) + numel(I) + numel(R);
            T(end+1) = t;
            S_prop(end+1) = numel(S) / N_current;
            I_prop(end+1) = numel(I) / N_current;
            R_prop(end+1) = numel(R) / N_current;

            I_count_list(end+1) = numel(I);
            
        end

        I_interp = interp1(T, I_prop, t_grid, 'linear', 'extrap');
        I_interp = max(0, min(1, I_interp));
        I_all(idx, :) = I_interp;
        S_interp = interp1(T, S_var, t_grid, 'linear', 'extrap');
        I_interp_var = interp1(T, I_var, t_grid, 'linear', 'extrap');
        R_interp = interp1(T, R_var, t_grid, 'linear', 'extrap');
        S_var_all(idx, :) = S_interp;
        I_var_all(idx, :) = I_interp_var;
        R_var_all(idx, :) = R_interp;

        if N == 316
            I_interp_count = interp1(T, I_count_list, t_grid, 'linear');
            I_count_all_316(:) = I_interp_count;
        end
    end

    %% Plotting
    % Plot S variance over time
    figure;
    hold on;
    for idx = 1:length(Ns)
        plot(t_grid, sqrt(S_var_all(idx, :)), 'Color', colors(idx,:), 'LineWidth', 1.2);
    end
    xlabel('Time (days)');
    ylabel('SD of Susceptible Proportion');
    title('SD of Susceptible Over Time');
    legend(arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false));
    grid on;
    hold off;
    
    % Plot I variance over time
    figure;
    hold on;
    for idx = 1:length(Ns)
        plot(t_grid, sqrt(I_var_all(idx, :)), 'Color', colors(idx,:), 'LineWidth', 1.2);
    end
    xlabel('Time (days)');
    ylabel('SD of Infected Proportion');
    title('SD of Infected Over Time');
    legend(arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false));
    grid on;
    hold off;
    
    % Plot R variance over time
    figure;
    hold on;
    for idx = 1:length(Ns)
        plot(t_grid, sqrt(R_var_all(idx, :)), 'Color', colors(idx,:), 'LineWidth', 1.2);
    end
    xlabel('Time (days)');
    ylabel('SD of Recovered Proportion');
    title('SD of Recovered Over Time');
    legend(arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false));
    grid on;
    hold off;

    % Plot Infected count for N = 316
    figure;
    plot(t_grid, I_count_all_316, 'r', 'LineWidth', 1.5);
    xlabel('Time (days)');
    ylabel('Infected Agent Count');
    title('Infected Agent Count Over Time (N = 316)');
    grid on;

end

function dydt = sir_system(t, y, beta, gamma, p1, p2)
    s = y(1);
    i = y(2);
    r = y(3);

    ds = -p1 * beta * s * i;
    di = p1 * beta * s * i - p2 * gamma * i;
    dr = p2 * gamma * i;

    dydt = [ds; di; dr];
end
