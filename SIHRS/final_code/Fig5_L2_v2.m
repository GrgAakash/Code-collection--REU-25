%% 
% SIHRS - L2 Discrepancy Analysis (D_N and R_N)
% Plots log(D_N) and log(R_N) for all compartments
% Based on Figure 5 methodology

clear all; close all;

function run_L2_analysis()
    rng(1);
    
    % ---- Parameters ----
    p.beta = 0.213882;
    p.pSI = 1.0;   p.pII = 0.1;
    p.pIH = 0.07143; p.pIR = 0.82857; p.pID = 0;
    p.pHH = 0.0;    p.pHR = 0.846;  p.pHD = 0.154;
    p.pRR = 0.02;   p.pRS = 0.98;
    p.gamma = 0.1429;
    p.alpha = 0.111;
    p.lambda = 0.0083;   % immunity waning
    p.T = 1000;
    p.dt = 0.01;
    
    p.Nvals = [2900, 1300, 550];
    p.cols = num2cell(lines(numel(p.Nvals)), 2);
    
    % initial conditions
    p.s0 = 0.96;  p.i0 = 0.035;  p.h0 = 0.005;
    p.r0 = 0.0;   p.d0 = 0.0;
    
    p.nruns = 70;
    
    % check params
    if abs(p.s0+p.i0+p.h0+p.r0+p.d0 - 1) > 1e-6
        error('initial conditions dont sum to 1');
    end
    
    R0 = p.pSI * p.beta / (p.gamma * (1 - p.pII));
    beta = p.beta;
    t = 0:p.dt:p.T;
    
    % get ODE solution
    [Tode, Yode] = solve_ode(beta, p);
    s_ode = interp1(Tode, Yode(:,1), t, 'linear', 'extrap');
    i_ode = interp1(Tode, Yode(:,2), t, 'linear', 'extrap');
    h_ode = interp1(Tode, Yode(:,3), t, 'linear', 'extrap');
    r_ode = interp1(Tode, Yode(:,4), t, 'linear', 'extrap');
    d_ode = interp1(Tode, Yode(:,5), t, 'linear', 'extrap');
    
    ode_all = {s_ode, i_ode, h_ode, r_ode, d_ode};
    
    % run for each N
    res = cell(length(p.Nvals), 1);
    for k = 1:length(p.Nvals)
        N = p.Nvals(k);
        fprintf('N=%d ... ', N);
        
        msd_s = zeros(1,length(t)); msd_i = zeros(1,length(t));
        msd_h = zeros(1,length(t)); msd_r = zeros(1,length(t));
        msd_d = zeros(1,length(t));
        
        for run = 1:p.nruns
            [sh,ih,hh,rh,dh,th] = gillespie(N, beta, p);
            si = interp1(th, sh, t, 'previous', 'extrap')/N;
            ii = interp1(th, ih, t, 'previous', 'extrap')/N;
            hi = interp1(th, hh, t, 'previous', 'extrap')/N;
            ri = interp1(th, rh, t, 'previous', 'extrap')/N;
            di = interp1(th, dh, t, 'previous', 'extrap')/N;
            
            msd_s = msd_s + (si - s_ode).^2;
            msd_i = msd_i + (ii - i_ode).^2;
            msd_h = msd_h + (hi - h_ode).^2;
            msd_r = msd_r + (ri - r_ode).^2;
            msd_d = msd_d + (di - d_ode).^2;
        end
        
        res{k}.msd = {msd_s/p.nruns, msd_i/p.nruns, msd_h/p.nruns, msd_r/p.nruns, msd_d/p.nruns};
        res{k}.N = N;
        fprintf('done\n');
    end
    
    % plotting
    comp = {'S','I','H','R','D'};
    
    % log(D_N) plots
    for c = 1:5
        figure('Position',[100 100 800 600]); hold on;
        for k = 1:length(p.Nvals)
            logD = log(sqrt(res{k}.msd{c}));
            plot(t, logD, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('N=%d',p.Nvals(k)));
        end
        xlabel('Time (days)','Interpreter','latex');
        ylabel(sprintf('$\\log(\\mathcal{D}_N^{(%d)})$',c),'Interpreter','latex');
        title(sprintf('Statistics related to $\\mathbf{%s}_N(t)$',lower(comp{c})),'Interpreter','latex');
        legend('Location','southeast'); grid on;
        set(gca,'FontSize',20); xlim([0 1000]);
        saveas(gcf, sprintf('logD_%s_R0_%.2f.png', comp{c}, R0));
    end
    
    % log(R_N) plots  
    for c = 1:5
        figure('Position',[100 100 800 600]); hold on;
        for k = 1:length(p.Nvals)
            sqrtmsd = sqrt(res{k}.msd{c});
            logR = log(sqrtmsd ./ abs(ode_all{c}));
            plot(t, logR, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('N=%d',p.Nvals(k)));
        end
        xlabel('Time (days)','Interpreter','latex');
        ylabel(sprintf('$\\log(\\mathcal{R}_N^{(%d)})$',c),'Interpreter','latex');
        title(sprintf('Statistics related to $\\mathbf{%s}_N(t)$',lower(comp{c})),'Interpreter','latex');
        legend('Location','southeast'); grid on;
        set(gca,'FontSize',20); xlim([0 1000]);
        saveas(gcf, sprintf('logR_%s_R0_%.2f.png', comp{c}, R0));
    end
    
    disp('Done - 10 plots saved');
end

%% Gillespie
function [S,I,H,R,D,T] = gillespie(N, beta, p)
    s = round(N*p.s0); i = round(N*p.i0); h = round(N*p.h0);
    r = round(N*p.r0); d = round(N*p.d0);
    
    maxev = round(10*p.T*(beta+p.gamma+p.alpha+p.lambda)*N);
    S = zeros(1,maxev); I = S; H = S; R = S; D = S; T = S;
    S(1)=s; I(1)=i; H(1)=h; R(1)=r; D(1)=d; T(1)=0;
    cnt = 1; t = 0;
    
    while t < p.T
        r1 = (beta/N)*s*i*p.pSI;
        r2 = p.gamma*i*p.pIR;
        r3 = p.gamma*i*p.pIH;
        r4 = p.gamma*i*p.pID;
        r5 = p.alpha*h*p.pHR;
        r6 = p.alpha*h*p.pHD;
        r7 = p.lambda*r*p.pRS;
        rtot = r1+r2+r3+r4+r5+r6+r7;
        
        if rtot == 0, break; end
        
        tau = -log(rand)/rtot;
        t = t + tau;
        if t > p.T, break; end
        
        u = rand*rtot;
        if u < r1
            s=s-1; i=i+1;
        elseif u < r1+r2
            i=i-1; r=r+1;
        elseif u < r1+r2+r3
            i=i-1; h=h+1;
        elseif u < r1+r2+r3+r4
            i=i-1; d=d+1;
        elseif u < r1+r2+r3+r4+r5
            h=h-1; r=r+1;
        elseif u < r1+r2+r3+r4+r5+r6
            h=h-1; d=d+1;
        else
            r=r-1; s=s+1;
        end
        
        cnt = cnt+1;
        S(cnt)=s; I(cnt)=i; H(cnt)=h; R(cnt)=r; D(cnt)=d; T(cnt)=t;
    end
    S=S(1:cnt); I=I(1:cnt); H=H(1:cnt); R=R(1:cnt); D=D(1:cnt); T=T(1:cnt);
end

%% ODE solver
function [T,Y] = solve_ode(beta, p)
    y0 = [p.s0; p.i0; p.h0; p.r0; p.d0];
    odefun = @(t,y) [...
        -beta*y(1)*y(2)*p.pSI + p.pRS*p.lambda*y(4);
        beta*y(1)*y(2)*p.pSI - p.gamma*(1-p.pII)*y(2);
        p.pIH*p.gamma*y(2) - p.alpha*(1-p.pHH)*y(3);
        p.pIR*p.gamma*y(2) + p.pHR*p.alpha*y(3) - p.pRS*p.lambda*y(4);
        p.pID*p.gamma*y(2) + p.pHD*p.alpha*y(3)];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [T,Y] = ode45(odefun, [0 p.T], y0, opts);
end

run_L2_analysis();

