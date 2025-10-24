% Export G values for all compartments to CSV file
% This is the exact code used to create SIHRS_G_values_N300_T300.csv
clear all; close all;

% Same parameters as main file
rng(1);
params.beta = 0.212;
params.gamma = 0.1;
params.pIH = 0.1060;
params.pIR = 0.8921;
params.pID = 0.0019;
params.pSI = 1.0;
params.T = 300;  % Using T = 300 as per your change
params.dt = 0.01;
params.initial_s = 0.96;
params.initial_i = 0.04;
params.initial_h = 0;
params.initial_r = 0;
params.initial_d = 0;
params.lambda = 0.0083;
params.pRS = 0.98;
params.alpha = 0.1;
params.pHR = 0.9882;
params.pHD = 0.0018;
params.pRR = 0.02;

N = 300;
t = 0:params.dt:params.T;

fprintf('Running Gillespie simulation...\n');
% Run one Gillespie simulation
[S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, params);

fprintf('Interpolating results...\n');
% Interpolate to fixed grid
S_interp = interp1(time_pts, S_hist, t, 'previous') / N;
I_interp = interp1(time_pts, I_hist, t, 'previous') / N;
H_interp = interp1(time_pts, H_hist, t, 'previous') / N;
R_interp = interp1(time_pts, R_hist, t, 'previous') / N;
D_interp = interp1(time_pts, D_hist, t, 'previous') / N;

% Handle NaN values
S_interp(isnan(S_interp)) = S_hist(end) / N;
I_interp(isnan(I_interp)) = I_hist(end) / N;
H_interp(isnan(H_interp)) = H_hist(end) / N;
R_interp(isnan(R_interp)) = R_hist(end) / N;
D_interp(isnan(D_interp)) = D_hist(end) / N;

fprintf('Calculating G values...\n');
% Calculate all G values (drift functions)
beta = params.beta;

% G1 (Susceptible)
G1 = -beta * S_interp .* I_interp + params.lambda * R_interp;

% G2 (Infected) 
G2 = -(params.gamma * (params.pIH + params.pIR + params.pID)) * I_interp + beta * S_interp .* I_interp;

% G3 (Hospitalized)
G3 = -(params.alpha * (params.pHR + params.pHD)) * H_interp + params.gamma * params.pIH * I_interp;

% G4 (Recovered)
G4 = -(params.lambda * params.pRS) * R_interp + params.gamma * params.pIR * I_interp + params.alpha * params.pHR * H_interp;

% G5 (Dead)
G5 = params.gamma * params.pID * I_interp + params.alpha * params.pHD * H_interp;

fprintf('Creating CSV data...\n');
% Create data matrix
csv_data = [t', G1', G2', G3', G4', G5'];

% Create header
header = {'Time', 'G1_Susceptible', 'G2_Infected', 'G3_Hospitalized', 'G4_Recovered', 'G5_Dead'};

% Write to CSV file
filename = 'SIHRS_G_values_N300_T300.csv';
fprintf('Writing to %s...\n', filename);

% Open file for writing
fid = fopen(filename, 'w');

% Write header
fprintf(fid, '%s', header{1});
for i = 2:length(header)
    fprintf(fid, ',%s', header{i});
end
fprintf(fid, '\n');

% Write data
for i = 1:length(t)
    fprintf(fid, '%.2f,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
            csv_data(i,1), csv_data(i,2), csv_data(i,3), csv_data(i,4), csv_data(i,5), csv_data(i,6));
end

fclose(fid);

fprintf('\n=== CSV FILE CREATED ===\n');
fprintf('Filename: %s\n', filename);
fprintf('Rows: %d (including header)\n', length(t) + 1);
fprintf('Columns: %d\n', size(csv_data, 2));
fprintf('Time range: %.2f to %.2f\n', t(1), t(end));
fprintf('Time step: %.2f\n', params.dt);

% Show some sample data
fprintf('\n=== SAMPLE DATA ===\n');
fprintf('Time    | G1        | G2        | G3        | G4        | G5\n');
fprintf('--------|-----------|-----------|-----------|-----------|----------\n');
sample_indices = [1, 1001, 2001, 5001, 10001, 14207, length(t)]; % Various time points including extinction
for idx = sample_indices
    if idx <= length(t)
        fprintf('%7.2f | %9.2e | %9.2e | %9.2e | %9.2e | %9.2e\n', ...
                t(idx), G1(idx), G2(idx), G3(idx), G4(idx), G5(idx));
    end
end

% Find key events
extinction_idx = find(I_interp == 0, 1, 'first');
if ~isempty(extinction_idx)
    fprintf('\n=== KEY EVENT ===\n');
    fprintf('Extinction at t = %.2f (row %d):\n', t(extinction_idx), extinction_idx + 1); % +1 for header
    fprintf('G1 = %.2e, G2 = %.2e, G3 = %.2e, G4 = %.2e, G5 = %.2e\n', ...
            G1(extinction_idx), G2(extinction_idx), G3(extinction_idx), G4(extinction_idx), G5(extinction_idx));
end

fprintf('\nCSV file ready for analysis!\n');

function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, params)
    % Initialize populations
    S = round(N * params.initial_s);
    I = round(N * params.initial_i);
    H = round(N * params.initial_h);
    R = round(N * params.initial_r);
    D = round(N * params.initial_d);
    
    % Pre-allocate arrays
    max_events = 100000;
    S_hist = zeros(1, max_events);
    I_hist = zeros(1, max_events);
    H_hist = zeros(1, max_events);
    R_hist = zeros(1, max_events);
    D_hist = zeros(1, max_events);
    time_pts = zeros(1, max_events);
    
    % Initial state
    S_hist(1) = S; I_hist(1) = I; H_hist(1) = H; R_hist(1) = R; D_hist(1) = D;
    time_pts(1) = 0;
    event_count = 1;
    current_time = 0;
    
    while current_time < params.T && I > 0
        % Calculate rates
        si_rate = (params.beta / N) * S * I * params.pSI;
        ir_rate = params.gamma * I * params.pIR;
        ih_rate = params.gamma * I * params.pIH;
        id_rate = params.gamma * I * params.pID;
        hr_rate = params.alpha * H * params.pHR;
        hd_rate = params.alpha * H * params.pHD;
        rs_rate = params.lambda * R * params.pRS;
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        
        if total_rate == 0
            break;
        end
        
        % Time to next event
        tau = -log(rand) / total_rate;
        current_time = current_time + tau;
        
        if current_time > params.T
            break;
        end
        
        % Choose which event occurs
        chance = rand * total_rate;
        
        if chance < si_rate
            S = S - 1; I = I + 1;
        elseif chance < (si_rate + ir_rate)
            I = I - 1; R = R + 1;
        elseif chance < (si_rate + ir_rate + ih_rate)
            I = I - 1; H = H + 1;
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate)
            I = I - 1; D = D + 1;
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate)
            H = H - 1; R = R + 1;
        elseif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate)
            H = H - 1; D = D + 1;
        else
            R = R - 1; S = S + 1;
        end
        
        % Record state
        event_count = event_count + 1;
        S_hist(event_count) = S;
        I_hist(event_count) = I;
        H_hist(event_count) = H;
        R_hist(event_count) = R;
        D_hist(event_count) = D;
        time_pts(event_count) = current_time;
    end
    
    % Trim arrays
    S_hist = S_hist(1:event_count);
    I_hist = I_hist(1:event_count);
    H_hist = H_hist(1:event_count);
    R_hist = R_hist(1:event_count);
    D_hist = D_hist(1:event_count);
    time_pts = time_pts(1:event_count);
end
