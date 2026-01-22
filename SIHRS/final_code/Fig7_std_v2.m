%% 
% SIHRS - Standard Deviation sigma(m_N)
% Instantaneous std of compartment fractions

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

res = cell(length(p.Nvals), 1);
for k = 1:length(p.Nvals)
    N = p.Nvals(k);
    fprintf('N=%d ... ', N);
    res{k} = runsim(N, beta, p, t);
    fprintf('done\n');
end

% sigma(m_N) plots
for c = 1:5
    figure('Position',[100 100 800 600]); hold on;
    for k = 1:length(p.Nvals)
        sigma = res{k}.std{c};
        plot(t, sigma, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('N=%d',p.Nvals(k)));
    end
    xlabel('Time (days)','Interpreter','latex');
    ylabel(sprintf('$\\sigma(\\mathbf{m}_N^{(%d)}(t))$',c),'Interpreter','latex');
    title(sprintf('Statistics related to $\\mathbf{%s}_N(t)$.',lower(comp{c})),'Interpreter','latex');
    legend('Location','northeast'); grid on;
    set(gca,'FontSize',20); xlim([5 1000]);
    saveas(gcf, sprintf('STD_%s_R0_%.2f.png', comp{c}, R0));
end

fprintf('\nFig7 done! 5 plots\n');

%% sim function
function result = runsim(N, beta, p, t)
    n = p.nruns;
    S_all = zeros(n, length(t)); I_all = zeros(n, length(t));
    H_all = zeros(n, length(t)); R_all = zeros(n, length(t)); D_all = zeros(n, length(t));
    
    for run = 1:n
        [sh,ih,hh,rh,dh,th] = gillespie(N, beta, p);
        S_all(run,:) = interp1(th, sh, t, 'previous', 'extrap')/N;
        I_all(run,:) = interp1(th, ih, t, 'previous', 'extrap')/N;
        H_all(run,:) = interp1(th, hh, t, 'previous', 'extrap')/N;
        R_all(run,:) = interp1(th, rh, t, 'previous', 'extrap')/N;
        D_all(run,:) = interp1(th, dh, t, 'previous', 'extrap')/N;
    end
    
    result.std = {std(S_all,0,1), std(I_all,0,1), std(H_all,0,1), std(R_all,0,1), std(D_all,0,1)};
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

