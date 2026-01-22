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
    all_data.ns = {}; all_data.ni = {}; all_data.nh = {}; all_data.nm = {}; all_data.nd = {};
    all_data.vs = {}; all_data.vi = {}; all_data.vh = {}; all_data.vm = {}; all_data.vd = {};
    all_data.T = {};

    constants.max_t = 30;
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
   

    %% Loop for each N
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
        data.vs = 0;
        data.vi = 0;
        data.vh = 0;
        data.vm = 0;
        data.vd = 0;
        data.T = 0;

        %% Perform Simulations
        data = simulate_ipc(S, I, H, M, D, constants, data);
        all_data.ns{end+1}=data.ns;all_data.ni{end+1}=data.ni;all_data.nh{end+1}=data.nh;all_data.nm{end+1}=data.nm;all_data.nd{end+1}=data.nd;
        data.vs(1) = data.vs(2); data.vi(1) = data.vi(2); data.vh(1) = data.vh(2); data.vm(1) = data.vm(2); data.vd(1) = data.vd(2);
        all_data.vs{end+1}=data.vs;all_data.vi{end+1}=data.vi;all_data.vh{end+1}=data.vh;all_data.vm{end+1}=data.vm;all_data.vd{end+1}=data.vd;
        all_data.T{end+1} = data.T;

        %% Graphing
        plot_proportions_vs_time(data, constants.N);
    end


    %% Deterministic Model
    beta = constants.susceptible_event_rate;
    gamma = constants.infected_event_rate;
    alpha = constants.hospitalized_event_rate;
    y0 = [s0; i0; 0; m0; 0]; 
    tspan = [0 constants.max_t];
    [t_det, y_det] = ode45(@(t, y) sihmd_system(t, y, constants, beta, gamma, alpha), tspan, y0);
    
    plot_compartment_vs_time(all_data.T, all_data.ns, t_det, y_det(:,1), Ns, "Susceptible Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.ni, t_det, y_det(:,2), Ns, "Infected Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nh, t_det, y_det(:,3), Ns, "Hospitalized Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nm, t_det, y_det(:,4), Ns, "Immune Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.nd, t_det, y_det(:,5), Ns, "Dead Proportion", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.vs, [], [], Ns, "Susceptible Variance", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.vi, [], [], Ns, "Infected Variance", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.vh, [], [], Ns, "Hospitalized Variance", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.vm, [], [], Ns, "Immune Variance", constants.max_t);
    plot_compartment_vs_time(all_data.T, all_data.vd, [], [], Ns, "Dead Variance", constants.max_t);
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
        if data.T(end) < constants.max_t
            data.ns(end+1)=data.ns(end);data.ni(end+1)=data.ni(end);data.nh(end+1)=data.nh(end);data.nm(end+1)=data.nm(end);data.nd(end+1)=data.nd(end);
            data.vs(end+1)=data.vs(end);data.vi(end+1)=data.vi(end);data.vh(end+1)=data.vh(end);data.vm(end+1)=data.vm(end);data.vd(end+1)=data.vd(end);
            data.T(end+1) = constants.max_t;
        end
        disp("Done!");
        return;
    end

    data.ns(end+1) = ns; 
    data.ni(end+1) = ni; 
    data.nh(end+1) = nh; 
    data.nm(end+1) = nm; 
    data.nd(end+1) = nd;
    data.vs(end+1) = (constants.susceptible_event_rate*ns/constants.N*ni/constants.N*constants.s_i_prob + constants.infected_event_rate*ni/constants.N*constants.i_s_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_s_prob) / constants.N;
    data.vi(end+1) = (constants.susceptible_event_rate*ns/constants.N*ni/constants.N*constants.s_i_prob + constants.infected_event_rate) / constants.N;
    data.vh(end+1) = (constants.infected_event_rate*ni/constants.N*constants.i_h_prob + constants.hospitalized_event_rate*nh/constants.N) / constants.N;
    data.vm(end+1) = (constants.infected_event_rate*ni/constants.N*constants.i_m_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_m_prob) / constants.N;
    data.vd(end+1) = (constants.infected_event_rate*ni/constants.N*constants.i_d_prob + constants.hospitalized_event_rate*nh/constants.N*constants.h_d_prob) / constants.N;

    dt = exprnd(1 / master_clock);
    data.T(end+1) = data.T(end) + dt;

    r = rand;
    if (r < s_clock / master_clock)                 % transition from s to i
        r2 = rand;
        if r2 < constants.s_i_prob                              % always s to i
            I(end+1) = S(end);
            S(end) = [];
        end
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

function dydt = sihmd_system(~, y, constants, beta, gamma, alpha)
    s = y(1);
    i = y(2);
    h = y(3);

    % Using constants from struct
    ds = -beta * s * i + gamma * constants.i_s_prob * i + alpha * constants.h_s_prob * h;
    di = beta * s * i - gamma * i;
    dh = gamma * constants.i_h_prob * i - alpha * h;
    dm = gamma * constants.i_m_prob * i + alpha * constants.h_m_prob * h;
    dd = gamma * constants.i_d_prob * i + alpha * constants.h_d_prob * h;

    dydt = [ds; di; dh; dm; dd];
end

function plot_proportions_vs_time(data, N) 
    figure;
    hold on;
    plot(data.T, data.ns / N, 'b', 'DisplayName', 'Susceptible', 'LineWidth', 1.5);
    plot(data.T, data.ni / N, 'r', 'DisplayName', 'Infected', 'LineWidth', 1.5);
    plot(data.T, data.nh / N, 'm', 'DisplayName', 'Hospitalized', 'LineWidth', 1.5);
    plot(data.T, data.nm / N, 'g', 'DisplayName', 'Immune', 'LineWidth', 1.5);
    plot(data.T, data.nd / N, 'k', 'DisplayName', 'Dead', 'LineWidth', 1.5);
    ylabel('Proportion');
    legend('Location', 'eastoutside');
    title(['All Proportions vs Time (N = ', num2str(N), ')']);
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
    plot(theo_t, theo_data, 'Color', 'k', 'LineWidth', 1.5);
    legend_labels{length(Ns)+1} = ['ODE'];

    ylabel(compartment);
    title([compartment, ' vs Time']);
    xlabel('Time');
    legend(legend_labels);
    xlim([0 max_t]);
    grid on;
    hold off;
end
