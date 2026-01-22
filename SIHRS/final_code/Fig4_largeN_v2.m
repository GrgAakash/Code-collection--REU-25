%% Figure 4: Large N Comparison

clear all; close all;

p.beta = 0.213882; p.gamma = 0.1429; p.alpha = 0.111; p.lambda = 0.0083;
p.pSI = 1.0; p.pII = 0.1;
p.pIH = 0.07143; p.pIR = 0.82857; p.pID = 0;
p.pHH = 0.0; p.pHR = 0.846; p.pHD = 0.154;
p.pRR = 0.02; p.pRS = 0.98;
p.tmax = 1000;
p.s0 = 0.96; p.i0 = 0.035; p.h0 = 0.005; p.r0 = 0.0; p.d0 = 0.0;

if abs(p.s0+p.i0+p.h0+p.r0+p.d0 - 1) > 1e-10
    error('initial conditions dont sum to 1');
end

Nvals = [15000];

res = cell(length(Nvals), 1);
for k = 1:length(Nvals)
    fprintf('N=%d ... ', Nvals(k));
    res{k} = run_stoch(Nvals(k), p);
    fprintf('done\n');
end

ode = solve_ode(p);
R0 = p.pSI * p.beta / (p.gamma * (1 - p.pII));
fprintf('R0 = %.4f\n', R0);

%% Stochastic S,I,R for each N

for k = 1:length(res)
    figure('Position', [100 100 1000 700]); hold on;
    plot(res{k}.T, res{k}.S, 'b', 'LineWidth', 4, 'DisplayName', '$\textbf{s}_N(t)$');
    plot(res{k}.T, res{k}.I, 'r', 'LineWidth', 4, 'DisplayName', '$\textbf{i}_N(t)$');
    plot(res{k}.T, res{k}.R, 'g', 'LineWidth', 4, 'DisplayName', '$\textbf{r}_N(t)$');
    xlabel('Time (days)', 'FontSize', 30);
    ylabel('$\textbf{s}_N(t)$, $\textbf{i}_N(t)$, $\textbf{r}_N(t)$', 'Interpreter', 'latex', 'FontSize', 30);
    grid on; xlim([0 p.tmax]); ylim([0 1]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
    saveas(gcf, sprintf('largeN_stoch_N%d_sir.png', Nvals(k)));
end

%% ODE plots

figure('Position', [100 100 800 600]); hold on;
plot(ode.T, ode.H, 'Color', '#A2142F', 'LineWidth', 4, 'DisplayName', '$h(t)$');
xlabel('Time (days)', 'FontSize', 30);
ylabel('$h(t)$', 'Interpreter', 'latex', 'FontSize', 30);
grid on; xlim([0 p.tmax]);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
saveas(gcf, 'largeN_det_h.png');

figure('Position', [100 100 800 600]); hold on;
plot(ode.T, ode.S, 'r', 'LineWidth', 4, 'DisplayName', '$s(t)$');
plot(ode.T, ode.R, 'b', 'LineWidth', 4, 'DisplayName', '$r(t)$');
xlabel('Time (days)', 'FontSize', 30);
ylabel('$s(t)$, $r(t)$', 'Interpreter', 'latex', 'FontSize', 30);
grid on; xlim([0 p.tmax]);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
saveas(gcf, 'largeN_det_sr.png');

figure('Position', [100 100 800 600]); hold on;
plot(ode.T, ode.I, 'Color', [1 0.5 0.5], 'LineWidth', 4, 'DisplayName', '$i(t)$');
plot(ode.T, ode.D, 'k', 'LineWidth', 4, 'DisplayName', '$d(t)$');
xlabel('Time (days)', 'FontSize', 30);
ylabel('$i(t)$, $d(t)$', 'Interpreter', 'latex', 'FontSize', 30);
grid on; xlim([0 p.tmax]);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
saveas(gcf, 'largeN_det_id.png');

%% Stochastic S,R for each N

for k = 1:length(res)
    figure('Position', [100 100 1000 700]); hold on;
    plot(res{k}.T, res{k}.S, 'Color', '#EDB120', 'LineWidth', 4, 'DisplayName', '$\textbf{s}_N(t)$');
    plot(res{k}.T, res{k}.R, 'g', 'LineWidth', 4, 'DisplayName', '$\textbf{r}_N(t)$');
    xlabel('Time (days)', 'FontSize', 30);
    ylabel('$\textbf{s}_N(t)$, $\textbf{r}_N(t)$', 'Interpreter', 'latex', 'FontSize', 30);
    grid on; xlim([0 p.tmax]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
    saveas(gcf, sprintf('largeN_stoch_N%d_sr.png', Nvals(k)));
end

%% Stochastic I,D for each N

for k = 1:length(res)
    figure('Position', [100 100 1000 700]); hold on;
    plot(res{k}.T, res{k}.I, 'Color', '#7E2F8E', 'LineWidth', 4, 'DisplayName', '$\textbf{i}_N(t)$');
    plot(res{k}.T, res{k}.D, 'c', 'LineWidth', 4, 'DisplayName', '$\textbf{d}_N(t)$');
    xlabel('Time (days)', 'FontSize', 30);
    ylabel('$\textbf{i}_N(t)$, $\textbf{d}_N(t)$', 'Interpreter', 'latex', 'FontSize', 30);
    grid on; xlim([0 p.tmax]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
    saveas(gcf, sprintf('largeN_stoch_N%d_id.png', Nvals(k)));
end

%% Stochastic H for each N

for k = 1:length(res)
    figure('Position', [100 100 1000 700]); hold on;
    plot(res{k}.T, res{k}.H, 'Color', '#77AC30', 'LineWidth', 4, 'DisplayName', '$\mathbf{h}_N(t)$');
    xlabel('Time (days)', 'FontSize', 30);
    ylabel('$\mathbf{h}_N(t)$', 'Interpreter', 'latex', 'FontSize', 30);
    grid on; xlim([0 p.tmax]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 40);
    set(gca, 'LooseInset', get(gca,'TightInset'), 'FontSize', 30);
    saveas(gcf, sprintf('largeN_H_N%d.png', Nvals(k)));
end

fprintf('\n=== SUMMARY ===\n');
fprintf('N\t\tpeak_I\t\tpeak_t\t\ts_inf\t\td_inf\n');
for k = 1:length(res)
    fprintf('%d\t\t%d\t\t%.1f\t\t%.3f\t\t%.3f\n', ...
        res{k}.N, res{k}.peak_I, res{k}.peak_t, res{k}.S(end), res{k}.D(end));
end

%% Gillespie
function r = run_stoch(N, p)
    S = round(p.s0*N); I = round(p.i0*N); H = round(p.h0*N);
    R = round(p.r0*N); D = round(p.d0*N);
    
    tot = S+I+H+R+D;
    if tot ~= N, S = S + (N-tot); end
    
    maxev = N * 100;
    Th = zeros(1,maxev); Sh = zeros(1,maxev); Ih = zeros(1,maxev);
    Hh = zeros(1,maxev); Rh = zeros(1,maxev); Dh = zeros(1,maxev);
    
    Th(1)=0; Sh(1)=S/N; Ih(1)=I/N; Hh(1)=H/N; Rh(1)=R/N; Dh(1)=D/N;
    ev = 1; tc = 0;
    
    while tc < p.tmax
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
            Th(ev)=tc; Sh(ev)=S/N; Ih(ev)=I/N; Hh(ev)=H/N; Rh(ev)=R/N; Dh(ev)=D/N;
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
        
        Sh(ev)=S/N; Ih(ev)=I/N; Hh(ev)=H/N; Rh(ev)=R/N; Dh(ev)=D/N;
    end
    
    r.T = Th(1:ev); r.S = Sh(1:ev); r.I = Ih(1:ev);
    r.H = Hh(1:ev); r.R = Rh(1:ev); r.D = Dh(1:ev);
    r.N = N;
    [r.peak_I, idx] = max(Ih(1:ev)*N);
    r.peak_t = Th(idx);
end

%% ODE
function r = solve_ode(p)
    y0 = [p.s0; p.i0; p.h0; p.r0; p.d0];
    [T, Y] = ode45(@(t,y) sihrs_ode(t,y,p), [0 p.tmax], y0, ...
        odeset('RelTol',1e-8,'AbsTol',1e-10));
    r.T = T; r.S = Y(:,1); r.I = Y(:,2); r.H = Y(:,3); r.R = Y(:,4); r.D = Y(:,5);
end

function dy = sihrs_ode(~, y, p)
    S=y(1); I=y(2); H=y(3); R=y(4);
    dy = zeros(5,1);
    dy(1) = -p.beta*S*I*p.pSI + p.pRS*p.lambda*R;
    dy(2) = p.beta*S*I*p.pSI - p.gamma*(1-p.pII)*I;
    dy(3) = p.pIH*p.gamma*I - p.alpha*(1-p.pHH)*H;
    dy(4) = p.pIR*p.gamma*I + p.pHR*p.alpha*H - p.pRS*p.lambda*R;
    dy(5) = p.pID*p.gamma*I + p.pHD*p.alpha*H;
end

