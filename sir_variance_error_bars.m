% SIR Variance Error Bars
% Caden McCann / Amy Carr
% This function generates a plot for each compartment, and compares the
% theoretical variance (determined by equations with starting parameters)
% with the simulated variance (determined by numerical integration over
% subintervals).

function sir_variance_error_bars()
    %% Model Parameters
    beta = 0.95; % beta*p1 / gamma*p2 is R0 reproduction number
    gamma = 1;
    p1 = 0.5;
    p2 = 0.5;
    max_T = 23;
    t_interval = 10; % for variance sampling
    t_grid = linspace(0, max_T, 1000);
    s0 = 0.96; % start with 96% susceptible
    i0 = 0.04;
    Ns = [316, 1000, 3162, 10000]; % varying population sizes

    %% Agent-Based Simulations
    S_var_all = zeros(length(Ns), length(t_grid));
    I_var_all = zeros(length(Ns), length(t_grid));
    R_var_all = zeros(length(Ns), length(t_grid));
    V_avg_all = zeros(length(Ns), 3, max_T-t_interval+1); % used to graph variance bars
    
    for idx = 1:length(Ns) % loop over each sample size N
        
        % actual starting populations
        N = Ns(idx);
        s0_abs = round(0.96 * N);
        i0_abs = round(0.04 * N);
        r0_abs = 0;

        % arrays containing the actual agents
        S = 1:s0_abs;
        I = (s0_abs+1):(s0_abs+i0_abs);
        R = [];

        % other initializations
        t = 0; % current time
        T = 0; % array of significant times
        S_prop = s0_abs / N;
        I_prop = i0_abs / N;
        R_prop = r0_abs / N;
        S_var = 0;
        I_var = 0;
        R_var = 0;

        %% Perform simulation for one N population
        while ~isempty(I) && t < max_T
            nI = numel(I);
            if nI < 1
                nI = 0;
            end    
            nS = numel(S);

            infection_rate = beta * nS * nI / (N^1.0); % density-dependent
            recovery_rate = gamma * nI;
            event_rate = infection_rate + recovery_rate; % master t clock
            
            dt = exprnd(1 / event_rate); % generate random waiting time
            t = t + dt; % exact time of next event

            % calculate infinitesimal variances at this t
            S_var(end+1) = p1 * beta * (nI/N) * (nS/N) / N;
            I_var(end+1) = (p1 * beta * (nI/N) * (nS/N) + p2 * gamma * (nI/N)) / N;
            R_var(end+1) = p2 * gamma * (nI/N) / N;

            if rand < (infection_rate / event_rate) % split master t clock into i clock
                if nS > 0 && rand < p1 % infection event actually happens
                    num = randi(nS); % get a random agent number
                    infected = S(num); % retrieve that agent
                    S(num) = []; % remove from susceptible population
                    I(end+1) = infected; % add to infected population
                end
            else % split master t clock into r clock
                if rand < p2 % recovery event actually happens
                    num = randi(nI); % get a random agent number
                    recovered = I(num); % retrieve that agent
                    I(num) = []; % remove from infected population
                    R(end+1) = recovered; % add to recovered population
                end
            end

            % update the counts and record relevant data for graphing
            N_current = numel(S) + numel(I) + numel(R);
            T(end+1) = t;
            S_prop(end+1) = numel(S) / N_current;
            I_prop(end+1) = numel(I) / N_current;
            R_prop(end+1) = numel(R) / N_current;            
        end

        %% If the simulation ended early (none infected at T=23), pad the rest of the data with zeros
        if t < max_T
            T(end+1) = max_T;
            S_var(end+1) = 0;
            I_var(end+1) = 0;
            R_var(end+1) = 0;
        end

        %% Create interpolated lines to smooth out the data points
        S_interp_var = interp1(T, S_var, t_grid, 'linear', 'extrap');
        I_interp_var = interp1(T, I_var, t_grid, 'linear', 'extrap');
        R_interp_var = interp1(T, R_var, t_grid, 'linear', 'extrap');
        S_var_all(idx, :) = S_interp_var;
        I_var_all(idx, :) = I_interp_var;
        R_var_all(idx, :) = R_interp_var;

        %% Find the integral of the variance over multiple intervals
        for interval_num = 1:max_T-t_interval+1 % 14 integrals when T=23
            t_start = interval_num-1;
            t_end = t_start + t_interval; % first integral over the interval [0,10]

            % variances [0.014 0.009 0.007 0.004] (for susceptible)
            % times     [  1.1  3.2   5.9   7.4 ]
            % start 3 end 8
            interval_mask = (t_grid >= t_start) & (t_grid <= t_end); % [0 1 1 0]
            t_sub = t_grid(interval_mask); % [3.2 5.9 7.4]

            % change in variance over the given interval for each compartment
            S_var_sub = S_interp_var(interval_mask); % [ 0.009  0.007 0.004 ]
            I_var_sub = I_interp_var(interval_mask);
            R_var_sub = R_interp_var(interval_mask);

            % numerical integration for each compartment, saving the average variances over each interval
            %                                (5.9-3.2) * (0.009+0.007)/2  /  (6-3)
            V_avg_all(idx, 1, interval_num) = 1*trapz(t_sub, S_var_sub) / (t_end - t_start);
            V_avg_all(idx, 2, interval_num) = 1*trapz(t_sub, I_var_sub) / (t_end - t_start);
            V_avg_all(idx, 3, interval_num) = 1*trapz(t_sub, R_var_sub) / (t_end - t_start);
        end 
    end

    %% Find the highest and lowest average variance over each interval for each compartment
    V_range = zeros(numel(Ns), 3, 2);
    for n = 1:numel(Ns) 
        for compartment = 1:3
            variance_column = squeeze(V_avg_all(n, compartment, :));
            V_range(n, compartment, 1) = min(variance_column);
            V_range(n, compartment, 2) = max(variance_column);
        end
    end
    
    compartments = ["Susceptible", "Infected", "Recovered"];
    intercepts = [log10(p1*beta*i0*s0); log10(p1*beta*i0*s0 + p2*gamma*i0); log10(p2*gamma*i0)];

    %% Graph the plots
    for compartment = 1:3
        figure;
        hold on;
        
        Ns_fine = logspace(log10(min(Ns)), log10(max(Ns)), 100);  % log scale
        slope = -1;
        intercept = intercepts(compartment);  % already log10 of coefficient
        
        y_theoretical = 10^intercept * Ns_fine.^slope;
        
        plot(Ns_fine, y_theoretical, 'LineWidth', 2);  % Red dashed line

        % create the vertical bar for each N
        for n = 1:numel(Ns) 
            x = Ns(n);
            y_min = V_range(n, compartment, 1);
            y_max = V_range(n, compartment, 2);
            line([x x], [y_min y_max], 'color', 'black', 'LineWidth', 1, 'Marker', '_');
        end

        
        set(gca, 'XScale', 'log'); % convert graph to log-log scale
        set(gca, 'YScale', 'log');
        xticks(Ns);
        xticklabels(string(round(log10(Ns),2)));
        xlabel('Log_{10}(N)');
        yticks = get(gca, 'YTick');
        yticklabels(string(log10(yticks)));
        ylabel('Log_{10}(Variance)');
        title(sprintf('Min-Max Variance Range: %s', compartments(compartment)));
        grid off;
        xlim([100,31622]); % set the window size
        ylim([1e-7,5e-4]);
    end
end
