function sir_combined_model()
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
    colors = lines(length(Ns));

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
        end

        I_interp = interp1(T, I_prop, t_grid, 'linear', 'extrap');
        I_interp = max(0, min(1, I_interp));
        I_all(idx, :) = I_interp;
    end

    %% Plotting
    figure;

    % Top subplot: Deterministic SIR model
    subplot(2,1,1);
    hold on;
    plot(t_ode, s_ode, '-k', 'LineWidth', 1.5);
    plot(t_ode, i_ode, 'Color', [0.6 0 0.8], 'LineWidth', 1.5);
    plot(t_ode, r_ode, 'Color', [0.4 0.8 1], 'LineWidth', 1.5);
    xlabel('Time (days)');
    ylabel('Proportion of population');
    title('Deterministic SIR Model (β = 0.95)');
    legend({'Susceptible', 'Infected', 'Recovered'});
    grid on;
    hold off;

    % Bottom subplot: Agent-based vs Deterministic Infected
    subplot(2,1,2);
    hold on;
    for idx = 1:length(Ns)
        plot(t_grid, I_all(idx, :), 'Color', colors(idx,:), 'LineWidth', 1);
    end
    plot(t_grid, i_ode_interp, '--k', 'LineWidth', 2); % ODE reference
    legend([arrayfun(@(n) sprintf('N = %d', n), Ns, 'UniformOutput', false), {'ODE Model'}]);
    xlabel('Time (days)');
    ylabel('Proportion Infected');
    title('Infected Proportion: Agent-Based vs ODE Model');
    grid on;
    hold off;
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
