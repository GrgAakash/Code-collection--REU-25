% Calculate exactly when G2 and V2 become zero for infected compartment
clear all; close all;

% Same parameters as main file
rng(1);  % Same seed for reproducibility

params.beta = 0.212;
params.pSI = 1.0;
params.pII = 0.0;
params.pIH = 0.04;
params.pIR = 0.959;
params.pID = 0.001;
params.gamma = 0.1;
params.alpha = 0.1;
params.lambda = 0.0083;
params.T = 1000;
params.dt = 0.01;
params.initial_s = 0.96;
params.initial_i = 0.04;
params.initial_h = 0;
params.initial_r = 0;
params.initial_d = 0;

N = 300;
beta = params.beta;
t = 0:params.dt:params.T;

fprintf('=== WHEN DO G2 AND V2 BECOME ZERO? ===\n');
fprintf('For Infected Compartment (2nd compartment)\n\n');

% Mathematical formulas from main code:
fprintf('Mathematical formulas:\n');
fprintf('V2(t) = (1/N) * (gamma * I(t) * (pIH + pIR + pID) + beta * pSI * S(t) * I(t))\n');
fprintf('G2(t) = -(gamma * (pIH + pIR + pID)) * I(t) + beta * S(t) * I(t)\n');
fprintf('G2(t) = I(t) * [beta * S(t) - gamma * (pIH + pIR + pID)]\n\n');

% Calculate critical values
gamma_eff = params.gamma * (params.pIH + params.pIR + params.pID);
S_critical = gamma_eff / beta;

fprintf('Key parameters:\n');
fprintf('gamma_eff = gamma * (pIH + pIR + pID) = %.3f\n', gamma_eff);
fprintf('beta = %.3f\n', beta);
fprintf('S_critical = gamma_eff/beta = %.3f\n', S_critical);
fprintf('pSI = %.3f\n', params.pSI);
fprintf('N = %d\n\n', N);

% Run simulation (same as main code)
[S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim_zeros(N, beta, params);

% Interpolate (same as main code)
S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
D_interp = interp1(time_pts, D_hist, t, 'previous') / N;

% Handle values beyond last event
S_interp(t > max(time_pts)) = S_hist(end) / N;
I_interp(t > max(time_pts)) = I_hist(end) / N;
H_interp(t > max(time_pts)) = H_hist(end) / N;
R_interp(t > max(time_pts)) = R_hist(end) / N;
D_interp(t > max(time_pts)) = D_hist(end) / N;

% Calculate V2 and G2 (same as main code)
V2 = (1/N) * (params.gamma * I_interp * (params.pIH + params.pIR + params.pID) + beta * params.pSI * S_interp .* I_interp);
G2 = -(params.gamma * (params.pIH + params.pIR + params.pID)) * I_interp + beta * S_interp .* I_interp;

fprintf('=== ANALYSIS OF ZEROS ===\n');

% Find when G2 = 0
G2_zero_idx = find(abs(G2) < 1e-12);
if ~isempty(G2_zero_idx)
    G2_zero_times = t(G2_zero_idx);
    fprintf('G2 = 0 at %d time points\n', length(G2_zero_idx));
    fprintf('First G2 = 0 at t = %.2f\n', G2_zero_times(1));
    fprintf('Last G2 = 0 at t = %.2f\n', G2_zero_times(end));
    
    % Check conditions when G2 = 0
    first_zero_idx = G2_zero_idx(1);
    fprintf('At first G2 = 0 (t = %.2f):\n', G2_zero_times(1));
    fprintf('  I = %.6f\n', I_interp(first_zero_idx));
    fprintf('  S = %.6f\n', S_interp(first_zero_idx));
    fprintf('  beta * S = %.6f\n', beta * S_interp(first_zero_idx));
    fprintf('  gamma_eff = %.6f\n', gamma_eff);
    fprintf('  V2 = %.2e\n', V2(first_zero_idx));
else
    fprintf('G2 never becomes exactly zero\n');
end

% Find when V2 = 0
V2_zero_idx = find(abs(V2) < 1e-12);
if ~isempty(V2_zero_idx)
    V2_zero_times = t(V2_zero_idx);
    fprintf('\nV2 = 0 at %d time points\n', length(V2_zero_idx));
    fprintf('First V2 = 0 at t = %.2f\n', V2_zero_times(1));
    fprintf('Last V2 = 0 at t = %.2f\n', V2_zero_times(end));
    
    % Check conditions when V2 = 0
    first_V2_zero_idx = V2_zero_idx(1);
    fprintf('At first V2 = 0 (t = %.2f):\n', V2_zero_times(1));
    fprintf('  I = %.6f\n', I_interp(first_V2_zero_idx));
    fprintf('  S = %.6f\n', S_interp(first_V2_zero_idx));
    fprintf('  G2 = %.2e\n', G2(first_V2_zero_idx));
else
    fprintf('\nV2 never becomes exactly zero\n');
end

% Find when I = 0
I_zero_idx = find(I_interp == 0);
if ~isempty(I_zero_idx)
    I_zero_time = t(I_zero_idx(1));
    fprintf('\nI = 0 starting at t = %.2f\n', I_zero_time);
    fprintf('At I = 0:\n');
    fprintf('  S = %.6f\n', S_interp(I_zero_idx(1)));
    fprintf('  G2 = %.2e\n', G2(I_zero_idx(1)));
    fprintf('  V2 = %.2e\n', V2(I_zero_idx(1)));
else
    fprintf('\nI never becomes exactly zero\n');
end

% Mathematical analysis
fprintf('\n=== MATHEMATICAL CONDITIONS ===\n');
fprintf('G2 = 0 when:\n');
fprintf('  Case 1: I = 0 (infection extinct)\n');
fprintf('  Case 2: I > 0 AND beta * S = gamma_eff (%.3f)\n', gamma_eff);
fprintf('         This means S = %.3f (critical threshold)\n', S_critical);

fprintf('\nV2 = 0 when:\n');
fprintf('  V2 = (1/N) * I * [gamma * (pIH + pIR + pID) + beta * pSI * S]\n');
fprintf('  Since all terms are positive when I > 0:\n');
fprintf('  V2 = 0 ONLY when I = 0\n');

% Check if S ever crosses critical threshold while I > 0
active_idx = I_interp > 0;
S_active = S_interp(active_idx);
t_active = t(active_idx);

S_crosses_critical = S_active <= S_critical;
if any(S_crosses_critical)
    cross_idx = find(S_crosses_critical, 1, 'first');
    cross_time = t_active(cross_idx);
    cross_S = S_active(cross_idx);
    cross_I = I_interp(active_idx);
    cross_I_val = cross_I(cross_idx);
    
    fprintf('\nS crosses critical threshold while I > 0:\n');
    fprintf('  Time: t = %.2f\n', cross_time);
    fprintf('  S = %.6f (critical = %.6f)\n', cross_S, S_critical);
    fprintf('  I = %.6f\n', cross_I_val);
    
    % Find corresponding G2 value
    full_idx = find(t == cross_time);
    if ~isempty(full_idx)
        fprintf('  G2 = %.2e\n', G2(full_idx));
        fprintf('  V2 = %.2e\n', V2(full_idx));
    end
else
    fprintf('\nS never crosses critical threshold while I > 0\n');
end

% Plot the results
figure('Position', [100, 100, 1400, 1000]);

subplot(2, 3, 1);
plot(t, S_interp, 'b-', 'LineWidth', 2); hold on;
plot(t, I_interp, 'r-', 'LineWidth', 2);
yline(S_critical, 'k--', 'LineWidth', 1, 'DisplayName', sprintf('S_{critical} = %.3f', S_critical));
xlabel('Time'); ylabel('Fraction'); title('Population Dynamics');
legend('S(t)', 'I(t)', 'S_{critical}', 'Location', 'best');
grid on;

subplot(2, 3, 2);
plot(t, G2, 'g-', 'LineWidth', 2);
xlabel('Time'); ylabel('G2(t)'); title('Drift Function G2(t)');
yline(0, 'k--', 'Alpha', 0.5);
if ~isempty(G2_zero_idx)
    plot(t(G2_zero_idx), G2(G2_zero_idx), 'ro', 'MarkerSize', 4);
end
grid on;

subplot(2, 3, 3);
plot(t, V2, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('V2(t)'); title('Variance Function V2(t)');
if ~isempty(V2_zero_idx)
    plot(t(V2_zero_idx), V2(V2_zero_idx), 'ro', 'MarkerSize', 4);
end
grid on;

subplot(2, 3, 4);
semilogy(t, abs(G2) + 1e-15, 'g-', 'LineWidth', 2);
xlabel('Time'); ylabel('|G2(t)|'); title('|G2(t)| (Log Scale)');
grid on;

subplot(2, 3, 5);
semilogy(t, V2 + 1e-15, 'b-', 'LineWidth', 2);
xlabel('Time'); ylabel('V2(t)'); title('V2(t) (Log Scale)');
grid on;

subplot(2, 3, 6);
active_mask = I_interp > 0;
if any(active_mask)
    R2_active = sqrt(V2(active_mask)) ./ abs(G2(active_mask));
    semilogy(t(active_mask), R2_active, 'm-', 'LineWidth', 2);
end
xlabel('Time'); ylabel('R2(t)'); title('R2(t) when I > 0 (Log Scale)');
grid on;

sgtitle('When do G2 and V2 become Zero?');

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim_zeros(N, beta, params)
    % Same Gillespie as main code
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
    max_events = round(10 * params.T * (beta + params.gamma + params.alpha + params.lambda) * N);
    S_hist = zeros(1, max_events); I_hist = zeros(1, max_events); H_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events); D_hist = zeros(1, max_events); time_pts = zeros(1, max_events);
    
    S_hist(1) = S; I_hist(1) = I; H_hist(1) = H; R_hist(1) = R; D_hist(1) = D; time_pts(1) = 0;
    event_count = 1; current_time = 0;
    
    while current_time < params.T
        si_rate = (beta / N) * S * I * params.pSI;
        ir_rate = params.gamma * I * params.pIR;
        ih_rate = params.gamma * I * params.pIH;
        id_rate = params.gamma * I * params.pID;
        hr_rate = params.alpha * H * params.pHR;
        hd_rate = params.alpha * H * params.pHD;
        rs_rate = params.lambda * R * params.pRS;
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        if total_rate == 0; break; end
        
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        if current_time > params.T; break; end
        
        chance = rand * total_rate;
        if chance < si_rate
            if S > 0; S = S - 1; I = I + 1; end
        elseif chance < (si_rate + ir_rate)
            if I > 0; I = I - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate)
            if I > 0; I = I - 1; H = H + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            if I > 0; I = I - 1; D = D + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            if H > 0; H = H - 1; R = R + 1; end
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            if H > 0; H = H - 1; D = D + 1; end
        else
            if R > 0; R = R - 1; S = S + 1; end
        end
        
        event_count = event_count + 1;
        S_hist(event_count) = S; I_hist(event_count) = I; H_hist(event_count) = H;
        R_hist(event_count) = R; D_hist(event_count) = D; time_pts(event_count) = current_time;
    end
    
    S_hist = S_hist(1:event_count); I_hist = I_hist(1:event_count); H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count); D_hist = D_hist(1:event_count); time_pts = time_pts(1:event_count);
end
