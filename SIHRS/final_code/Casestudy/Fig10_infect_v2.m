%% Figure 10: Infection Case Study (Carson City, NV)

clear all; close all;

N = 58639;  % Carson City population

% Load initial conditions from data
try
    data = readtable('carson_city_combined.csv');
    data.date = datetime(data.date, 'InputFormat', 'yyyy-MM-dd');
    start_date = datetime('2020-03-25');
    
    idx = find(data.date == start_date, 1);
    if isempty(idx)
        [~, idx] = min(abs(datenum(data.date) - datenum(start_date)));
    end
    
    i0 = data.cases(idx) / N;
    d0 = data.deaths(idx) / N;
    h0 = 0; r0 = 0;
    s0 = 1 - (i0 + h0 + r0 + d0);
    
    fprintf('Initial: I=%d, D=%d, S=%d\n', data.cases(idx), data.deaths(idx), round(s0*N));
catch
    i0 = 3/N; d0 = 0; h0 = 0; r0 = 0;
    s0 = 1 - (i0 + h0 + r0 + d0);
end

p.beta = 0.15711428571; p.gamma = 0.1429; p.alpha = 0.111; p.lambda = 0.0083;
p.pSI = 1.0; p.pII = 0.1;
p.pIH = 0.092; p.pIR = 0.8080; p.pID = 0;
p.pHH = 0.0; p.pHR = 0.846; p.pHD = 0.154;
p.pRR = 0.02; p.pRS = 0.98;
p.tmax = 620;
p.s0 = s0; p.i0 = i0; p.h0 = h0; p.r0 = r0; p.d0 = d0;

R0 = p.beta * p.pSI / p.gamma * (1 - p.pII);
fprintf('R0 = %.4f\n', R0);

nsims = 20;
res = cell(nsims, 1);
for k = 1:nsims
    fprintf('Sim %d/%d\n', k, nsims);
    res{k} = run_stoch(N, p);
end

plot_results(res, N, p);

%% Gillespie
function r = run_stoch(N, p)
    S = round(p.s0*N); I = round(p.i0*N); H = round(p.h0*N);
    R = round(p.r0*N); D = round(p.d0*N);
    
    tot = S+I+H+R+D;
    if tot ~= N, S = S + (N-tot); end
    
    maxev = N * 30;
    Th = zeros(maxev,1); Ih = zeros(maxev,1); Dh = zeros(maxev,1);
    
    Th(1)=0; Ih(1)=I; Dh(1)=D;
    ev = 1; tc = 0;
    
    while I > 0 && tc < p.tmax
        r_si = p.pSI * p.beta * S * I / N;
        r_rs = p.pRS * p.lambda * R;
        r_ih = p.gamma * I * p.pIH;
        r_ir = p.gamma * I * p.pIR;
        r_id = p.gamma * I * p.pID;
        r_hr = p.alpha * H * p.pHR;
        r_hd = p.alpha * H * p.pHD;
        
        tot = r_si + r_rs + r_ih + r_ir + r_id + r_hr + r_hd;
        if tot == 0, break; end
        
        tau = -log(rand)/tot;
        tc = tc + tau;
        if tc > p.tmax
            tc = p.tmax;
            ev = ev+1;
            Th(ev)=tc; Ih(ev)=I; Dh(ev)=D;
            break;
        end
        
        ev = ev+1;
        Th(ev) = tc;
        
        ch = rand * tot;
        if ch < r_si
            S=S-1; I=I+1;
        elseif ch < r_si+r_rs
            R=R-1; S=S+1;
        elseif ch < r_si+r_rs+r_ih
            I=I-1; H=H+1;
        elseif ch < r_si+r_rs+r_ih+r_ir
            I=I-1; R=R+1;
        elseif ch < r_si+r_rs+r_ih+r_ir+r_id
            I=I-1; D=D+1;
        elseif ch < r_si+r_rs+r_ih+r_ir+r_id+r_hr
            H=H-1; R=R+1;
        else
            H=H-1; D=D+1;
        end
        
        Ih(ev)=I; Dh(ev)=D;
    end
    
    r.T = Th(1:ev); r.I = Ih(1:ev); r.D = Dh(1:ev);
end

%% Plotting
function plot_results(res, N, p)
    t_grid = (0:p.tmax)';
    nsims = length(res);
    all_I = zeros(nsims, length(t_grid));
    all_D = zeros(nsims, length(t_grid));
    
    for k = 1:nsims
        if length(res{k}.T) > 1
            itp_I = griddedInterpolant(res{k}.T, res{k}.I, 'linear', 'none');
            itp_D = griddedInterpolant(res{k}.T, res{k}.D, 'linear', 'none');
            for j = 1:length(t_grid)
                if t_grid(j) >= min(res{k}.T) && t_grid(j) <= max(res{k}.T)
                    all_I(k,j) = itp_I(t_grid(j));
                    all_D(k,j) = itp_D(t_grid(j));
                end
            end
        end
    end
    
    % Active deaths (14-day window)
    window = 14;
    all_active_D = zeros(size(all_D));
    for k = 1:nsims
        for t = 1:length(t_grid)
            if t <= window
                all_active_D(k,t) = all_D(k,t);
            else
                all_active_D(k,t) = all_D(k,t) - all_D(k,t-window);
            end
        end
    end
    
    mean_I = mean(all_I, 1)';
    lower_I = quantile(all_I, 0.05, 1)';
    upper_I = quantile(all_I, 0.95, 1)';
    
    % Load real data
    start_date = datetime('2020-03-25');
    real_I = zeros(length(t_grid), 1);
    
    try
        data = readtable('carson_city_combined.csv');
        data.date = datetime(data.date, 'InputFormat', 'yyyy-MM-dd');
        idx = find(data.date == start_date, 1);
        if isempty(idx)
            [~, idx] = min(abs(datenum(data.date) - datenum(start_date)));
        end
        
        dates = data.date(idx:end);
        cumcases = data.cases(idx:end);
        recovery = 14;
        shifted = [zeros(recovery,1); cumcases(1:end-recovery)];
        active = max(cumcases - shifted, 0);
        
        days_from_start = days(dates - dates(1));
        real_I = interp1(days_from_start, active, t_grid, 'linear', 0);
    catch ME
        fprintf('Could not load data: %s\n', ME.message);
    end
    
    % ODE solution
    ode_I = solve_ode(p, N, t_grid);
    
    % Convert to proportions
    mean_I_prop = mean_I / N;
    lower_I_prop = lower_I / N;
    upper_I_prop = upper_I / N;
    real_I_prop = real_I / N;
    ode_I_prop = ode_I / N;
    all_I_prop = all_I / N;
    
    % Date labels
    tick_pos = 0:90:p.tmax;
    tick_dates = start_date + days(tick_pos);
    date_labels = cellstr(datestr(tick_dates, 'mm/dd/yy'));
    
    %% Figure 1: 90% envelope + ODE + real
    figure;
    fill([t_grid; flipud(t_grid)], [upper_I_prop; flipud(lower_I_prop)], ...
         [1 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;
    plot(t_grid, ode_I_prop, '--', 'Color', '#7E2F8E', 'LineWidth', 4);
    plot(t_grid, real_I_prop, 'g-', 'LineWidth', 3);
    xlabel('Date (mm/dd/yy)'); ylabel('the infected compartment fraction');
    xlim([0 p.tmax]);
    ylim([0, max([upper_I_prop; real_I_prop]) * 1.1]);
    xticks(tick_pos); xticklabels(date_labels);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 15);
    legend('$90\%$ envelope', '$i(t)$', 'real data', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 20);
    saveas(gcf, 'Fig10_infect_envelope.png');
    
    %% Figure 2: Individual trajectories + real
    figure;
    for k = 1:nsims
        plot(t_grid, all_I_prop(k,:), 'Color', [1 0.5 0.5], 'LineWidth', 1);
        hold on;
    end
    plot(t_grid, real_I_prop, 'g-', 'LineWidth', 2.5);
    xlabel('Date (mm/dd/yy)'); ylabel('the infected compartment fraction');
    xlim([0 p.tmax]);
    ylim([0, max([max(all_I_prop(:)); max(real_I_prop)]) * 1.1]);
    xticks(tick_pos); xticklabels(date_labels);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 15);
    h1 = plot(NaN, NaN, 'Color', [1 0.5 0.5], 'LineWidth', 2);
    h2 = plot(NaN, NaN, 'g-', 'LineWidth', 3);
    legend([h1, h2], {'$\textbf{i}_N(t)$', 'real data'}, 'Interpreter', 'latex', 'FontSize', 20);
    saveas(gcf, 'Fig10_infect_trajectories.png');
end

%% ODE solver
function I = solve_ode(p, N, t_grid)
    y0 = [p.s0*N; p.i0*N; p.h0*N; p.r0*N; p.d0*N];
    [t, y] = ode45(@(t,y) sihrs_ode(t,y,p,N), [0 p.tmax], y0, ...
        odeset('RelTol',1e-8,'AbsTol',1e-10));
    I = max(interp1(t, y(:,2), t_grid, 'linear', 'extrap'), 0);
end

function dy = sihrs_ode(~, y, p, N)
    S=y(1); I=y(2); H=y(3); R=y(4);
    dy = zeros(5,1);
    dy(1) = -p.beta*p.pSI*S*I/N + p.lambda*p.pRS*R;
    dy(2) = p.beta*p.pSI*S*I/N - p.gamma*(p.pIH+p.pIR+p.pID)*I;
    dy(3) = p.gamma*p.pIH*I - p.alpha*(p.pHR+p.pHD)*H;
    dy(4) = p.gamma*p.pIR*I + p.alpha*p.pHR*H - p.lambda*p.pRS*R;
    dy(5) = p.gamma*p.pID*I + p.alpha*p.pHD*H;
end

