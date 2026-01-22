%% 
% SIHRS - Martingale Variance: sigma(G^cum) and sigma(M)
% Based on decomposition: m_N(t) = m_0 + G^cum + M

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

% sigma(G^cum) plots
for c = 1:5
    figure('Position',[100 100 800 600]); hold on;
    for k = 1:length(p.Nvals)
        sg = res{k}.Gcum_std{c};
        plot(t, sg, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('N=%d',p.Nvals(k)));
    end
    xlabel('Time (days)','Interpreter','latex');
    ylabel(sprintf('$\\sigma(\\mathcal{G}_N^{(%d)}(t))$',c),'Interpreter','latex');
    title(sprintf('Statistics related to $\\mathbf{%s}_N(t)$.',lower(comp{c})),'Interpreter','latex');
    legend('Location','northeast'); grid on;
    set(gca,'FontSize',20); xlim([5 1000]);
    saveas(gcf, sprintf('VarGcum%s_R0_%.2f.png', comp{c}, R0));
end

% sigma(M) plots
for c = 1:5
    figure('Position',[100 100 800 600]); hold on;
    for k = 1:length(p.Nvals)
        sm = res{k}.M_std{c};
        plot(t, sm, 'Color', p.cols{k}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('N=%d',p.Nvals(k)));
    end
    xlabel('Time (days)','Interpreter','latex');
    ylabel(sprintf('$\\sigma(\\mathcal{M}_N^{(%d)}(t))$',c),'Interpreter','latex');
    title(sprintf('Statistics related to $\\mathbf{%s}_N(t)$.',lower(comp{c})),'Interpreter','latex');
    legend('Location','northeast'); grid on;
    set(gca,'FontSize',20); xlim([5 1000]);
    saveas(gcf, sprintf('VarMart%s_R0_%.2f.png', comp{c}, R0));
end

fprintf('\nFig8 done! 10 plots\n');

%% sim function
function result = runsim(N, beta, p, t)
    n = p.nruns;
    
    % accumulators for G^cum
    sGc = zeros(1,length(t)); sGc2 = zeros(1,length(t));
    sGi = zeros(1,length(t)); sGi2 = zeros(1,length(t));
    sGh = zeros(1,length(t)); sGh2 = zeros(1,length(t));
    sGr = zeros(1,length(t)); sGr2 = zeros(1,length(t));
    sGd = zeros(1,length(t)); sGd2 = zeros(1,length(t));
    
    % accumulators for M
    sMs = zeros(1,length(t)); sMs2 = zeros(1,length(t));
    sMi = zeros(1,length(t)); sMi2 = zeros(1,length(t));
    sMh = zeros(1,length(t)); sMh2 = zeros(1,length(t));
    sMr = zeros(1,length(t)); sMr2 = zeros(1,length(t));
    sMd = zeros(1,length(t)); sMd2 = zeros(1,length(t));
    
    for run = 1:n
        [sh,ih,hh,rh,dh,th] = gillespie(N, beta, p);
        si = interp1(th, sh, t, 'previous', 'extrap')/N;
        ii = interp1(th, ih, t, 'previous', 'extrap')/N;
        hi = interp1(th, hh, t, 'previous', 'extrap')/N;
        ri = interp1(th, rh, t, 'previous', 'extrap')/N;
        di = interp1(th, dh, t, 'previous', 'extrap')/N;
        
        % drift G_l(t)
        Gs = -p.pSI*beta.*si.*ii + p.pRS*p.lambda.*ri;
        Gi = p.pSI*beta.*si.*ii - p.gamma*(1-p.pII).*ii;
        Gh = p.pIH*p.gamma.*ii - p.alpha*(1-p.pHH).*hi;
        Gr = p.pIR*p.gamma.*ii + p.pHR*p.alpha.*hi - p.pRS*p.lambda.*ri;
        Gd = p.pID*p.gamma.*ii + p.pHD*p.alpha.*hi;
        
        % G^cum = integral G dt
        Gcs = [0, cumsum(Gs(1:end-1).*diff(t))];
        Gci = [0, cumsum(Gi(1:end-1).*diff(t))];
        Gch = [0, cumsum(Gh(1:end-1).*diff(t))];
        Gcr = [0, cumsum(Gr(1:end-1).*diff(t))];
        Gcd = [0, cumsum(Gd(1:end-1).*diff(t))];
        
        % M = m - m_0 - G^cum
        Ms = si - p.s0 - Gcs;
        Mi = ii - p.i0 - Gci;
        Mh = hi - p.h0 - Gch;
        Mr = ri - p.r0 - Gcr;
        Md = di - p.d0 - Gcd;
        
        % accumulate G^cum
        sGc = sGc + Gcs; sGc2 = sGc2 + Gcs.^2;
        sGi = sGi + Gci; sGi2 = sGi2 + Gci.^2;
        sGh = sGh + Gch; sGh2 = sGh2 + Gch.^2;
        sGr = sGr + Gcr; sGr2 = sGr2 + Gcr.^2;
        sGd = sGd + Gcd; sGd2 = sGd2 + Gcd.^2;
        
        % accumulate M
        sMs = sMs + Ms; sMs2 = sMs2 + Ms.^2;
        sMi = sMi + Mi; sMi2 = sMi2 + Mi.^2;
        sMh = sMh + Mh; sMh2 = sMh2 + Mh.^2;
        sMr = sMr + Mr; sMr2 = sMr2 + Mr.^2;
        sMd = sMd + Md; sMd2 = sMd2 + Md.^2;
    end
    
    % Var = E[X^2] - E[X]^2
    result.Gcum_std = {sqrt(max(sGc2/n-(sGc/n).^2,0)), sqrt(max(sGi2/n-(sGi/n).^2,0)), ...
                       sqrt(max(sGh2/n-(sGh/n).^2,0)), sqrt(max(sGr2/n-(sGr/n).^2,0)), ...
                       sqrt(max(sGd2/n-(sGd/n).^2,0))};
    result.M_std = {sqrt(max(sMs2/n-(sMs/n).^2,0)), sqrt(max(sMi2/n-(sMi/n).^2,0)), ...
                    sqrt(max(sMh2/n-(sMh/n).^2,0)), sqrt(max(sMr2/n-(sMr/n).^2,0)), ...
                    sqrt(max(sMd2/n-(sMd/n).^2,0))};
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

