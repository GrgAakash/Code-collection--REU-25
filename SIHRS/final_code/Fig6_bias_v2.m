%% 
% SIHRS - Bias and Variance Accumulation (E_N and F_N)
% log(sqrt(E_N)) and sqrt(F_N) for all compartments

clear all; close all;
rng(1);

% params
p.beta = 0.213882; p.pSI = 1.0; p.pII = 0.1;
p.pIH = 0.07143; p.pIR = 0.82857; p.pID = 0;
p.pHH = 0.00; p.pHR = 0.846; p.pHD = 0.154;
p.pRR = 0.02; p.pRS = 0.98;
p.gamma = 0.1429; p.alpha = 0.111; p.lambda = 0.0083;
p.T = 1000; p.dt = 0.01;

p.Nvals = [2900, 1300, 550];
p.cols = num2cell(lines(numel(p.Nvals)), 2);

p.s0 = 0.96; p.i0 = 0.035; p.h0 = 0.005; p.r0 = 0.0; p.d0 = 0.0;
p.nruns = 70;

% check
if abs(p.s0+p.i0+p.h0+p.r0+p.d0 - 1) > 1e-6
    error('initial conditions dont sum to 1');
end

R0 = p.pSI * p.beta / (p.gamma * (1 - p.pII));
beta = p.beta;
t = 0:p.dt:p.T;
comp = {'S','I','H','R','D'};

% get ODE solution
[~, Yode] = solve_ode(beta, p);
sode = interp1(linspace(0,p.T,size(Yode,1)), Yode(:,1), t);
iode = interp1(linspace(0,p.T,size(Yode,1)), Yode(:,2), t);
hode = interp1(linspace(0,p.T,size(Yode,1)), Yode(:,3), t);
rode = interp1(linspace(0,p.T,size(Yode,1)), Yode(:,4), t);
dode = interp1(linspace(0,p.T,size(Yode,1)), Yode(:,5), t);
ode_sol = {sode, iode, hode, rode, dode};

res = cell(length(p.Nvals), 1);
for k = 1:length(p.Nvals)
    N = p.Nvals(k);
    fprintf('N=%d ... ', N);
    res{k} = runsim(N, beta, p, t, ode_sol);
    fprintf('done\n');
end

% log(sqrt(E_N)) plots
for c = 1:5
    figure('Position',[100 100 800 600]); hold on;
    for k = 1:length(p.Nvals)
        logSqrtE = log(sqrt(max(res{k}.E{c}, 1e-20)));
        plot(t, logSqrtE, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('N=%d',p.Nvals(k)));
    end
    xlabel('Time (days)','Interpreter','latex');
    ylabel(sprintf('$\\log\\sqrt{\\mathcal{E}_N^{(%d)}}$',c),'Interpreter','latex');
    title(sprintf('Bias: $\\mathbf{%s}_N(t)$',lower(comp{c})),'Interpreter','latex');
    legend('Location','southeast'); grid on;
    set(gca,'FontSize',20); xlim([0 1000]);
    saveas(gcf, sprintf('logSqrtE_%s_R0_%.2f.png', comp{c}, R0));
end

% sqrt(F_N) plots
for c = 1:5
    figure('Position',[100 100 800 600]); hold on;
    for k = 1:length(p.Nvals)
        sqrtF = sqrt(max(res{k}.F{c}, 0));
        plot(t, sqrtF, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('N=%d',p.Nvals(k)));
    end
    xlabel('Time (days)','Interpreter','latex');
    ylabel(sprintf('$\\sqrt{\\mathcal{F}_N^{(%d)}}$',c),'Interpreter','latex');
    title(sprintf('Var Accumulation: $\\mathbf{%s}_N(t)$',lower(comp{c})),'Interpreter','latex');
    legend('Location','northwest'); grid on;
    set(gca,'FontSize',20); xlim([0 1000]);
    saveas(gcf, sprintf('sqrtF_%s_R0_%.2f.png', comp{c}, R0));
end

fprintf('\nFig6 done! 10 plots\n');

%% sim function
function result = runsim(N, beta, p, t, ode_sol)
    n = p.nruns;
    % accumulators for mean and MSD
    sum_s = zeros(1,length(t)); sum_i = zeros(1,length(t));
    sum_h = zeros(1,length(t)); sum_r = zeros(1,length(t)); sum_d = zeros(1,length(t));
    sum_s2 = zeros(1,length(t)); sum_i2 = zeros(1,length(t));
    sum_h2 = zeros(1,length(t)); sum_r2 = zeros(1,length(t)); sum_d2 = zeros(1,length(t));
    
    for run = 1:n
        [sh,ih,hh,rh,dh,th] = gillespie(N, beta, p);
        si = interp1(th, sh, t, 'previous', 'extrap')/N;
        ii = interp1(th, ih, t, 'previous', 'extrap')/N;
        hi = interp1(th, hh, t, 'previous', 'extrap')/N;
        ri = interp1(th, rh, t, 'previous', 'extrap')/N;
        di = interp1(th, dh, t, 'previous', 'extrap')/N;
        
        sum_s = sum_s + si; sum_i = sum_i + ii; sum_h = sum_h + hi;
        sum_r = sum_r + ri; sum_d = sum_d + di;
        sum_s2 = sum_s2 + si.^2; sum_i2 = sum_i2 + ii.^2; sum_h2 = sum_h2 + hi.^2;
        sum_r2 = sum_r2 + ri.^2; sum_d2 = sum_d2 + di.^2;
    end
    
    % means
    mean_s = sum_s/n; mean_i = sum_i/n; mean_h = sum_h/n;
    mean_r = sum_r/n; mean_d = sum_d/n;
    
    % variance
    var_s = sum_s2/n - mean_s.^2; var_i = sum_i2/n - mean_i.^2;
    var_h = sum_h2/n - mean_h.^2; var_r = sum_r2/n - mean_r.^2;
    var_d = sum_d2/n - mean_d.^2;
    
    % E_N = integral (E[m_N] - m_ODE)^2 dt (bias)
    E_s = [0, cumsum((mean_s(1:end-1) - ode_sol{1}(1:end-1)).^2 .* diff(t))];
    E_i = [0, cumsum((mean_i(1:end-1) - ode_sol{2}(1:end-1)).^2 .* diff(t))];
    E_h = [0, cumsum((mean_h(1:end-1) - ode_sol{3}(1:end-1)).^2 .* diff(t))];
    E_r = [0, cumsum((mean_r(1:end-1) - ode_sol{4}(1:end-1)).^2 .* diff(t))];
    E_d = [0, cumsum((mean_d(1:end-1) - ode_sol{5}(1:end-1)).^2 .* diff(t))];
    
    % F_N = integral Var(m_N) dt (variance accumulation)
    F_s = [0, cumsum(max(var_s(1:end-1),0) .* diff(t))];
    F_i = [0, cumsum(max(var_i(1:end-1),0) .* diff(t))];
    F_h = [0, cumsum(max(var_h(1:end-1),0) .* diff(t))];
    F_r = [0, cumsum(max(var_r(1:end-1),0) .* diff(t))];
    F_d = [0, cumsum(max(var_d(1:end-1),0) .* diff(t))];
    
    result.E = {E_s, E_i, E_h, E_r, E_d};
    result.F = {F_s, F_i, F_h, F_r, F_d};
    result.N = N;
end

%% gillespie
function [sh,ih,hh,rh,dh,th] = gillespie(N, beta, p)
    S=round(N*p.s0); I=round(N*p.i0); H=round(N*p.h0); R=round(N*p.r0); D=round(N*p.d0);
    maxev = round(10*p.T*(beta+p.gamma+p.alpha+p.lambda)*N);
    sh=zeros(1,maxev); ih=zeros(1,maxev); hh=zeros(1,maxev);
    rh=zeros(1,maxev); dh=zeros(1,maxev); th=zeros(1,maxev);
    sh(1)=S; ih(1)=I; hh(1)=H; rh(1)=R; dh(1)=D; th(1)=0;
    ev=1; tc=0;
    
    while tc < p.T
        r_si = (beta/N)*S*I*p.pSI;
        r_ir = p.gamma*I*p.pIR; r_ih = p.gamma*I*p.pIH; r_id = p.gamma*I*p.pID;
        r_hr = p.alpha*H*p.pHR; r_hd = p.alpha*H*p.pHD;
        r_rs = p.lambda*R*p.pRS;
        tot = r_si+r_ir+r_ih+r_id+r_hr+r_hd+r_rs;
        if tot==0, break; end
        
        tau = -log(rand)/tot;
        tc = tc + tau;
        if tc > p.T, break; end
        
        ch = rand*tot;
        if ch < r_si
            S=S-1; I=I+1;
        elseif ch < r_si+r_ir
            I=I-1; R=R+1;
        elseif ch < r_si+r_ir+r_ih
            I=I-1; H=H+1;
        elseif ch < r_si+r_ir+r_ih+r_id
            I=I-1; D=D+1;
        elseif ch < r_si+r_ir+r_ih+r_id+r_hr
            H=H-1; R=R+1;
        elseif ch < r_si+r_ir+r_ih+r_id+r_hr+r_hd
            H=H-1; D=D+1;
        else
            R=R-1; S=S+1;
        end
        ev = ev+1;
        sh(ev)=S; ih(ev)=I; hh(ev)=H; rh(ev)=R; dh(ev)=D; th(ev)=tc;
    end
    sh=sh(1:ev); ih=ih(1:ev); hh=hh(1:ev); rh=rh(1:ev); dh=dh(1:ev); th=th(1:ev);
end

%% ODE solver
function [T, Y] = solve_ode(beta, p)
    y0 = [p.s0, p.i0, p.h0, p.r0, p.d0];
    [T, Y] = ode45(@(t,y) sihrs_ode(t,y,beta,p), [0 p.T], y0);
end

function dydt = sihrs_ode(~, y, beta, p)
    S=y(1); I=y(2); H=y(3); R=y(4);
    dS = -p.pSI*beta*S*I + p.pRS*p.lambda*R;
    dI = p.pSI*beta*S*I - p.gamma*(1-p.pII)*I;
    dH = p.pIH*p.gamma*I - p.alpha*(1-p.pHH)*H;
    dR = p.pIR*p.gamma*I + p.pHR*p.alpha*H - p.pRS*p.lambda*R;
    dD = p.pID*p.gamma*I + p.pHD*p.alpha*H;
    dydt = [dS; dI; dH; dR; dD];
end

