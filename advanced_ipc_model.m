%        -> d
%             -> m
% s -> i -> h -> d
%             -> s
%        -> m
%        -> s
% Agents start susceptible, infected, or immune
%
% Infections spread at the rate of the possible pairs between s and i still
% Infected agents can move to one of four compartments: immune, susceptible 
% (currently healthy but able to get sick again), hospitalized, and dead
%
% Immune and dead populations never change
% 
% Hospitalized agents can move to one of three compartments: immune, 
% susceptible (currently healthy but able to get sick again), and dead

% rates provided by https://data.who.int/dashboards/covid19/cases
function advanced_ipc_model() 
    s0 = 0.96;
    i0 = 0.04;
    m0 = 1 - s0 - i0; % assume that at time 0, no one is immune
    Ns = [316,1000,3162,10000];

    constants.susceptible_event_rate = 1.3; 
    constants.s_i_prob = 1; % susceptible agents can only become infected

    constants.infected_event_rate = 1;
    constants.i_m_prob = 0.2; % 20% of people who catch covid have long-term immunity
    constants.i_d_prob = 0.010175; % 1052449 deaths with covid as primary factor, out of 103436829 cases
    constants.i_h_prob = 0.04; % 4% of people who catch covid become hospitalized
    constants.i_s_prob = 1 - constants.i_m_prob - constants.i_h_prob - ...
        constants.i_d_prob; % the rest of infected become healthy, but not forever
    
    constants.hospitalized_event_rate = 1; 
    constants.h_m_prob = 0.6; % since the infection was more severe, agents develop stronger immunity
    constants.h_d_prob = 0.15; % on average, 15% of hospitalized patients die
    constants.h_s_prob = 1 - constants.h_m_prob - ...
        constants.h_d_prob; % the rest become healthy, but not forever
    
    constants.max_t = 30;
    
    % Deterministic setup
    beta = constants.susceptible_event_rate;
    gamma = constants.infected_event_rate;   % infected event rate
    alpha = constants.hospitalized_event_rate;   % hospitalized event rate

    y0 = [s0; i0; 0; m0; 0]; % Initial conditions: s, i, h, m, d
    tspan = [0 constants.max_t];

    % Solve ODE
    [t_det, y_det] = ode45(@(t, y) sihmd_system(t, y, beta, gamma, alpha, ...
        constants.i_s_prob, constants.i_h_prob, constants.i_m_prob, constants.i_d_prob, ...
        constants.h_s_prob, constants.h_m_prob, constants.h_d_prob), tspan, y0);

% Plot
plot_deterministic(t_det, y_det);


    for i = 1:numel(Ns)
        constants.N = Ns(i);
        I0 = round(i0 * constants.N);
        M0 = round(m0 * constants.N);
        S0 = constants.N - I0 - M0; % defined like this to ensure the agent sum is correct

        % Arrays of agents
        S = 1:S0;
        I = (S0+1):(S0+I0);
        H = [];
        M = (S0+I0+1):constants.N;
        D = [];

        data.ns = S0;
        data.ni = I0;
        data.nh = 0;
        data.nm = M0;
        data.nd = 0;
        data.T = 0;

        %% Perform Simulations
        data = simulate_ipc(S, I, H, M, D, constants, data);

        %% Feed the function 'data' to plot all compartments, or just name a single one
        %plot_compartment_vs_time(data.T, data, constants.N, "All");
        plot_compartment_vs_time(data.T, data.ns, constants.N, "Susceptible")
    end

end

function data = simulate_ipc(S, I, H, M, D, constants, data) 
    ns = numel(S);
    ni = numel(I);
    nh = numel(H);
    nm = numel(M);
    nd = numel(D);

    s_clock = constants.susceptible_event_rate * ns * ni / constants.N;
    i_clock = constants.infected_event_rate * ni;
    h_clock = constants.hospitalized_event_rate * nh;
    master_clock = s_clock + i_clock + h_clock;

    if data.T(end) >= constants.max_t || master_clock == 0
        disp("Done!");
        return;
    end

    data.ns(end+1) = ns;
    data.ni(end+1) = ni;
    data.nh(end+1) = nh;
    data.nm(end+1) = nm;
    data.nd(end+1) = nd;

    dt = exprnd(1 / master_clock);
    data.T(end+1) = data.T(end) + dt;

    r = rand;
    if (r < s_clock / master_clock)                 % transition from s to i
        I(end+1) = S(end);
        S(end) = [];
    elseif (r < (s_clock + i_clock) / master_clock) % transition from i to [s h m d]
        r2 = rand;
        if r2 < constants.i_h_prob                              % i to h
            H(end+1) = I(end); 
            I(end) = [];
        elseif r2 < (constants.i_h_prob + constants.i_m_prob)   % i to m
            M(end+1) = I(end);
            I(end) = [];
        elseif r2 < (constants.i_h_prob + constants.i_m_prob + constants.i_d_prob) % i to d
            D(end+1) = I(end);
            I(end) = []; 
        else                                                    % i to s
            S(end+1) = I(end);
            I(end) = [];
        end
    else                                            % transition from h to [s m d]
        r2 = rand;
        if r2 < constants.h_s_prob                              % h to s
            S(end+1) = H(end);
            H(end) = [];
        elseif r2 < (constants.h_s_prob + constants.h_m_prob)   % h to m
            M(end+1) = H(end);
            H(end) = [];
        else                                                    % h to d
            D(end+1) = H(end);
            H(end) = [];
        end
    end

    data = simulate_ipc(S, I, H, M, D, constants, data);
end

function plot_compartment_vs_time(time, data, N, compartment)
    figure;
    hold on;

    if isstruct(data)
        % Full data struct provided
        plot(time, data.ns / N, 'b', 'DisplayName', 'Susceptible');
        plot(time, data.ni / N, 'r', 'DisplayName', 'Infected');
        plot(time, data.nh / N, 'k', 'DisplayName', 'Hospitalized');
        plot(time, data.nm / N, 'g', 'DisplayName', 'Immune');
        plot(time, data.nd / N, 'm', 'DisplayName', 'Dead');
        ylabel('Proportion');
        legend('Location', 'eastoutside');
        title(['Compartment Proportions vs Time (N = ', num2str(N), ')']);
    elseif isnumeric(data)
        % Single compartment provided
        plot(time, data / N, 'b', 'DisplayName', 'Single Compartment');
        ylabel('Proportion');
        legend(compartment);
        title([compartment, ' Proportions vs Time (N = ', num2str(N), ')']);
    else
        error('Input must be either a struct or a numeric array.');
    end

    xlabel('Time');
    grid on;
    hold off;
end

function dydt = sihmd_system(t, y, beta, gamma, alpha, p_IS, p_IH, p_IM, p_ID, p_HS, p_HM, p_HD)
    s = y(1);
    i = y(2);
    h = y(3);
    m = y(4);
    d = y(5);

    ds = -beta*s*i + gamma*p_IS*i + alpha*p_HS*h;
    di = beta*i*s - gamma*i;
    dh = gamma*p_IH*i - p_HM*alpha*h;
    dm = gamma*p_IM*i + alpha*p_HM*h;
    dd = gamma*p_ID*i + alpha*p_HD*h;

    dydt = [ds; di; dh; dm; dd];
end

function plot_deterministic(t, y)
    s = y(:,1);
    i = y(:,2);
    h = y(:,3);
    m = y(:,4);
    d = y(:,5);

    figure;
    plot(t, s, '-k'); hold on;
    plot(t, i, 'Color', [0.6 0 0.8]);
    plot(t, h, 'Color', [0.2 0.8 0]);
    plot(t, m, 'Color', [0 0 1]);
    plot(t, d, 'Color', [0.4 0.8 1]);

    title('Deterministic SIHMD Model');
    xlabel('Time (days)');
    ylabel('Proportion of Population');
    legend({'Susceptible', 'Infected', 'Hospitalized', 'Immune', 'Dead'});
    grid on;
end
