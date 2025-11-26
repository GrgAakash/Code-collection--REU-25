% Export sqrt(Variance) values for all compartments to CSV file
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

fprintf('Calculating sqrt(Variance) values...\n');
% Calculate all sqrt(Variance) values
beta = params.beta;

% sqrt(V1) for Susceptible
V1 = (1/N) * (beta * params.pSI * S_interp .* I_interp + params.lambda * R_interp * params.pRS);
sqrt_V1 = sqrt(V1);

% sqrt(V2) for Infected
V2 = (1/N) * (params.gamma * I_interp * (params.pIH + params.pIR + params.pID) + beta * params.pSI * S_interp .* I_interp);
sqrt_V2 = sqrt(V2);

% sqrt(V2_residual) - same as used in main code for R2 calculation
V2_residual = V2 + (1/N) * 1e-10;  % Small residual variance
sqrt_V2_residual = sqrt(V2_residual);

% sqrt(V3) for Hospitalized
V3 = (1/N) * (params.gamma * I_interp * params.pIH + params.alpha * H_interp * (params.pHR + params.pHD));
sqrt_V3 = sqrt(V3);

% sqrt(V4) for Recovered
V4 = (1/N) * (params.gamma * I_interp * params.pIR + params.alpha * H_interp * params.pHR + params.lambda * R_interp * params.pRS);
sqrt_V4 = sqrt(V4);

% sqrt(V5) for Dead
V5 = (1/N) * (params.gamma * I_interp * params.pID + params.alpha * H_interp * params.pHD);
sqrt_V5 = sqrt(V5);

fprintf('Creating CSV data...\n');
% Create data matrix
csv_data = [t', sqrt_V1', sqrt_V2', sqrt_V2_residual', sqrt_V3', sqrt_V4', sqrt_V5'];

% Create header
header = {'Time', 'sqrt_V1_Susceptible', 'sqrt_V2_Infected', 'sqrt_V2_Residual', 'sqrt_V3_Hospitalized', 'sqrt_V4_Recovered', 'sqrt_V5_Dead'};

% Write to CSV file
filename = 'SIHRS_sqrt_Variance_N300_T300.csv';
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
    fprintf(fid, '%.2f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
            csv_data(i,1), csv_data(i,2), csv_data(i,3), csv_data(i,4), csv_data(i,5), csv_data(i,6), csv_data(i,7));
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
fprintf('Time    | sqrt_V1   | sqrt_V2   | sqrt_V2_R | sqrt_V3   | sqrt_V4   | sqrt_V5\n');
fprintf('--------|-----------|-----------|-----------|-----------|-----------|----------\n');
sample_indices = [1, 1001, 2001, 5001, 10001, 14207, length(t)]; % Various time points including extinction
for idx = sample_indices
    if idx <= length(t)
        fprintf('%7.2f | %9.2e | %9.2e | %9.2e | %9.2e | %9.2e | %9.2e\n', ...
                t(idx), sqrt_V1(idx), sqrt_V2(idx), sqrt_V2_residual(idx), sqrt_V3(idx), sqrt_V4(idx), sqrt_V5(idx));
    end
end

% Find key events
extinction_idx = find(I_interp == 0, 1, 'first');
if ~isempty(extinction_idx)
    fprintf('\n=== KEY EVENT ===\n');
    fprintf('Extinction at t = %.2f (row %d):\n', t(extinction_idx), extinction_idx + 1); % +1 for header
    fprintf('sqrt_V1 = %.2e, sqrt_V2 = %.2e, sqrt_V2_R = %.2e\n', ...
            sqrt_V1(extinction_idx), sqrt_V2(extinction_idx), sqrt_V2_residual(extinction_idx));
    fprintf('sqrt_V3 = %.2e, sqrt_V4 = %.2e, sqrt_V5 = %.2e\n', ...
            sqrt_V3(extinction_idx), sqrt_V4(extinction_idx), sqrt_V5(extinction_idx));
end

% Show variance behavior
fprintf('\n=== VARIANCE ANALYSIS ===\n');
fprintf('At extinction (t = %.2f):\n', t(extinction_idx));
fprintf('- sqrt_V2 (natural) = %.2e (goes to 0 when I=0)\n', sqrt_V2(extinction_idx));
fprintf('- sqrt_V2_residual  = %.2e (stays > 0 due to residual)\n', sqrt_V2_residual(extinction_idx));
fprintf('- Difference shows residual effect for R2 calculation\n');

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
