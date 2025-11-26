%% 
% SIHRS (all 5 compartments): Qℓ and RQℓ with path-average-first construction
% -------------------------------------------------------------------------
% For ℓ = 1..5 (S,I,H,R,D), per path p:
%   Vℓ^(p)(t), Gℓ^(p)(t)  (instantaneous variance & drift in fraction space)
% Path-average-first:
%   meanVℓ(t) = ⟨Vℓ^(p)(t)⟩_p,   meanGℓ(t) = ⟨Gℓ^(p)(t)⟩_p
%   meanCumVℓ(t) = ∫_0^t meanVℓ(τ)dτ
%   meanCumGℓ(t) = ℓ(0) + ∫_0^t meanGℓ(τ)dτ     % includes matching initial fraction
%
% Qℓ(t)  = meanCumGℓ(t) .* meanVℓ(t) - 2 .* meanCumVℓ(t) .* meanGℓ(t)
% RQℓ(t) = Qℓ(t) ./ ( 2 .* sqrt(meanCumVℓ(t)) .* (meanCumGℓ(t).^2) )
%
% Per-path (no averaging) analogues use CumVℓ^(p) and CumGℓ^(p) (with same initial add-in)
% and the SAME normalization form with sqrt(CumVℓ^(p)).
%
% RQ plots use a symmetric log transform (symlog) only; Q plots stay in linear scale.
%
% Output (per N) generated for Qℓ_mean, Qℓ_paths, RQℓ_mean_symlog, RQℓ_paths_symlog.

clear; close all; clc;
sihrs_all_Q_RQ_main();

% ========================= MAIN =========================
function sihrs_all_Q_RQ_main()
    rng(1); % reproducible

    % ---- Parameters ----
    params.beta   = 0.212;
    params.pSI    = 1.0;
    params.pII    = 0.0;
    params.pIH    = 0.04;
    params.pIR    = 0.959;
    params.pID    = 0.001;
    params.pHH    = 0.01;
    params.pHR    = 0.9882;
    params.pHD    = 0.0018;
    params.pRR    = 0.02;
    params.pRS    = 0.98;
    params.gamma  = 0.100346667;
    params.alpha  = 0.1;
    params.lambda = 0.0083;

    params.T      = 1000;
    params.dt     = 0.01;

    params.N_values     = [600, 2900, 6000];
    params.num_paths    = 10;     % paths simulated per N
    params.n_plot_paths = 10;     % how many random I/H paths to draw (<= num_paths)

    % Initial fractions (sum to 1)
    params.initial_s = 0.96;
    params.initial_i = 0.035;
    params.initial_h = 0.005;
    params.initial_r = 0.0;
    params.initial_d = 0.0;

    % symlog setting (RQ plots only)
    params.symlog_linthresh = 1e-5;

    validate_params(params);

    % Info: R0
    R0 = params.pSI * params.beta / (params.gamma * (1 - params.pII));
    fprintf('Calculated R0 = %.4f\n', R0);

    % Time grid
    t = 0:params.dt:params.T;

    % Deterministic reference (fractions) for I/H overlays
    det = solve_deterministic_sihrs(params);

    % Loop over N
    base_seed = 3636;
    results = cell(numel(params.N_values), 1);

    for idx = 1:numel(params.N_values)
        N = params.N_values(idx);
        fprintf('Running %d stochastic path(s) for N = %d ...\n', params.num_paths, N);

        % Accumulators for mean-first construction (ℓ=1..5)
        sumV = zeros(5, numel(t));
        sumG = zeros(5, numel(t));

        % Per-path storage for Q and RQ (no averaging)
        Q_paths  = zeros(5, params.num_paths, numel(t));   % (ℓ, p, t)
        RQ_paths = zeros(5, params.num_paths, numel(t));   % (ℓ, p, t)

        % For overlay & path plots
        first_S = []; first_I = []; first_H = []; first_R = []; first_D = [];
        nplot = min(params.n_plot_paths, params.num_paths);
        pick  = randperm(params.num_paths, nplot);
        I_plot = zeros(nplot, numel(t));
        H_plot = zeros(nplot, numel(t));
        ipos = 0;

        for p = 1:params.num_paths
            stream = RandStream('Threefry','Seed', base_seed + 1000*idx + p);

            % Gillespie (counts)
            [S_hist,I_hist,H_hist,R_hist,D_hist,time_pts] = gillespie_sim(N, params, stream);

            % Interpolate to uniform grid (fractions)
            [S,I,H,R,D] = interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N);

            % First path store for ODE overlays
            if p == 1
                first_S = S; first_I = I; first_H = H; first_R = R; first_D = D;
            end

            % Subset for I/H random-path figures
            if any(pick == p)
                ipos = ipos + 1;
                I_plot(ipos,:) = I;
                H_plot(ipos,:) = H;
            end

            % Instantaneous (per-path) for all ℓ
            [V1,G1,V2,G2,V3,G3,V4,G4,V5,G5] = inst_all(S,I,H,R,D,N,params);

            % Pack
            V = [V1; V2; V3; V4; V5];
            G = [G1; G2; G3; G4; G5];

            % Accumulate for means
            sumV = sumV + V;
            sumG = sumG + G;

            % Per-path cumulatives (CumG adds initial fractions)
            init = [params.initial_s, params.initial_i, params.initial_h, params.initial_r, params.initial_d];
            for ell = 1:5
                CumV_p = cumtrapz(t, V(ell,:));
                CumG_p = init(ell) + cumtrapz(t, G(ell,:));

                % Per-path Qℓ (no averaging)
                Qp = CumG_p .* V(ell,:) - 2 .* CumV_p .* G(ell,:);
                Q_paths(ell,p,:) = Qp;

                % Per-path RQℓ (no averaging) with sqrt(CumV)
                eps_den = 1e-15;
                CumV_clamped = max(CumV_p, 0);
                denp = 2 .* sqrt(CumV_clamped) .* (CumG_p.^2);
                RQ_paths(ell,p,:) = Qp ./ max(denp, eps_den);
            end
        end

        % MEAN across paths, THEN integrate (path-average-first)
        m = params.num_paths;
        meanV = sumV / m;     % 5 x T
        meanG = sumG / m;     % 5 x T

        meanCumV = zeros(5, numel(t));
        meanCumG = zeros(5, numel(t));
        init = [params.initial_s, params.initial_i, params.initial_h, params.initial_r, params.initial_d];
        for ell = 1:5
            meanCumV(ell,:) = cumtrapz(t, meanV(ell,:));
            meanCumG(ell,:) = init(ell) + cumtrapz(t, meanG(ell,:));
        end

        % Mean-based Q and RQ for all ℓ
        Q_mean  = zeros(5, numel(t));
        RQ_mean = zeros(5, numel(t));
        eps_den = 1e-15;
        for ell = 1:5
            Q_mean(ell,:)  = meanCumG(ell,:) .* meanV(ell,:) - 2 .* meanCumV(ell,:) .* meanG(ell,:);
            den            = 2 .* sqrt(max(meanCumV(ell,:),0)) .* (meanCumG(ell,:).^2);
            RQ_mean(ell,:) = Q_mean(ell,:) ./ max(den, eps_den);
        end

        % Bundle results for this N
        R = struct();
        R.N = N; R.t = t;

        % first path for overlay
        R.S1 = first_S; R.I1 = first_I; R.H1 = first_H; R.R1 = first_R; R.D1 = first_D;

        % random-path sets for I/H figures
        R.I_paths_plot = I_plot(1:ipos, :);
        R.H_paths_plot = H_plot(1:ipos, :);

        % mean-first series and mean Q/RQ for all ℓ
        R.meanV = meanV; R.meanG = meanG; R.meanCumV = meanCumV; R.meanCumG = meanCumG;
        R.Q_mean = Q_mean; R.RQ_mean = RQ_mean;

        % per-path Q/RQ (no averaging) for all ℓ
        R.Q_paths  = Q_paths;
        R.RQ_paths = RQ_paths;

        results{idx} = R;
        fprintf('Completed N = %d with %d path(s)\n', N, m);
    end

    % Plots:
    plot_overlay_IH_first_path(results, det);     % first-path vs ODE (I/H)
    plot_IH_random_paths(results);                % random I/H paths with mean overlay
    plot_Q_RQ_all(results, params);               % Qℓ (linear) & RQℓ (symlog) for ℓ=1..5

    fprintf('All plots generated.\n');
end

% ========================= HELPERS =========================
function validate_params(params)
    if any(params.N_values <= 0) || any(mod(params.N_values,1) ~= 0)
        error('N_values must be positive integers');
    end
    if params.pSI <= 0 || params.pSI > 1 || params.pII < 0 || params.pII > 1
        error('pSI in (0,1], pII in [0,1]');
    end
    if any([params.gamma, params.alpha, params.lambda] <= 0)
        error('gamma, alpha, lambda must be positive');
    end
    if abs(params.initial_s + params.initial_i + params.initial_h + params.initial_r + params.initial_d - 1) > 1e-12
        error('Initial fractions must sum to 1');
    end
    if params.T <= 0 || params.dt <= 0
        error('T and dt must be positive');
    end
    if params.num_paths < 1 || mod(params.num_paths,1) ~= 0
        error('num_paths must be a positive integer');
    end
end

% Instantaneous drift & variance for ALL compartments (per path, fractions)
function [V1,G1,V2,G2,V3,G3,V4,G4,V5,G5] = inst_all(S, I, H, R, D, N, P)
    % S (ℓ=1)
    V1 = (1/N) * ( P.beta*P.pSI .* S .* I + P.lambda*P.pRS .* R );
    G1 = -P.beta .* S .* I + P.lambda .* R;

    % I (ℓ=2)
    V2 = (1/N) * ( P.gamma*(P.pIH + P.pIR + P.pID) .* I + P.beta*P.pSI .* S .* I );
    G2 = -P.gamma*(P.pIH + P.pIR + P.pID) .* I + P.beta .* S .* I;

    % H (ℓ=3)
    V3 = (1/N) * ( P.gamma*P.pIH .* I + P.alpha*(P.pHR + P.pHD) .* H );
    G3 =  P.gamma*P.pIH .* I - P.alpha*(P.pHR + P.pHD) .* H;

    % R (ℓ=4)
    V4 = (1/N) * ( P.gamma*P.pIR .* I + P.alpha*P.pHR .* H + P.lambda*P.pRS .* R );
    G4 =  P.gamma*P.pIR .* I + P.alpha*P.pHR .* H - P.lambda .* R;

    % D (ℓ=5)
    V5 = (1/N) * ( P.gamma*P.pID .* I + P.alpha*P.pHD .* H );
    G5 =  P.gamma*P.pID .* I + P.alpha*P.pHD .* H;
end

% Gillespie (counts)
function [S_hist, I_hist, H_hist, R_hist, D_hist, time_pts] = gillespie_sim(N, P, stream)
    S = round(N * P.initial_s);
    I = round(N * P.initial_i);
    H = round(N * P.initial_h);
    R = round(N * P.initial_r);
    D = round(N * P.initial_d);

    max_events = max(1, round(10 * P.T * (P.beta + P.gamma + P.alpha + P.lambda) * N));
    S_hist = zeros(1,max_events); I_hist = zeros(1,max_events);
    H_hist = zeros(1,max_events); R_hist = zeros(1,max_events);
    D_hist = zeros(1,max_events); time_pts = zeros(1,max_events);

    S_hist(1)=S; I_hist(1)=I; H_hist(1)=H; R_hist(1)=R; D_hist(1)=D; time_pts(1)=0;
    event_count = 1; current_time = 0;

    while current_time < P.T
        si_rate = (P.beta / N) * S * I * P.pSI;
        ir_rate = P.gamma * I * P.pIR;
        ih_rate = P.gamma * I * P.pIH;
        id_rate = P.gamma * I * P.pID;
        hr_rate = P.alpha * H * P.pHR;
        hd_rate = P.alpha * H * P.pHD;
        rs_rate = P.lambda * R * P.pRS;

        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate;
        if total_rate <= 0, break; end

        tau = -log(rand(stream)) / total_rate;
        current_time = current_time + tau;
        if current_time > P.T, break; end

        u = rand(stream) * total_rate;
        if u < si_rate
            if S>0, S=S-1; I=I+1; end
        elseif u < si_rate + ir_rate
            if I>0, I=I-1; R=R+1; end
        elseif u < si_rate + ir_rate + ih_rate
            if I>0, I=I-1; H=H+1; end
        elseif u < si_rate + ir_rate + ih_rate + id_rate
            if I>0, I=I-1; D=D+1; end
        elseif u < si_rate + ir_rate + ih_rate + id_rate + hr_rate
            if H>0, H=H-1; R=R+1; end
        elseif u < si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate
            if H>0, H=H-1; D=D+1; end
        else
            if R>0, R=R-1; S=S+1; end
        end

        event_count = event_count + 1;
        if event_count > numel(S_hist)
            grow = round(0.5*numel(S_hist))+1;
            S_hist=[S_hist,zeros(1,grow)]; I_hist=[I_hist,zeros(1,grow)];
            H_hist=[H_hist,zeros(1,grow)]; R_hist=[R_hist,zeros(1,grow)];
            D_hist=[D_hist,zeros(1,grow)]; time_pts=[time_pts,zeros(1,grow)];
        end
        S_hist(event_count)=S; I_hist(event_count)=I; H_hist(event_count)=H;
        R_hist(event_count)=R; D_hist(event_count)=D; time_pts(event_count)=current_time;
    end

    S_hist=S_hist(1:event_count); I_hist=I_hist(1:event_count);
    H_hist=H_hist(1:event_count); R_hist=R_hist(1:event_count);
    D_hist=D_hist(1:event_count); time_pts=time_pts(1:event_count);
end

% Interpolate counts -> fractions on uniform grid (left-continuous)
function [S,I,H,R,D] = interpolate_results(S_hist,I_hist,H_hist,R_hist,D_hist,time_pts,t,N)
    S = interp1(time_pts, S_hist, t, 'previous', 'extrap')/N;
    I = interp1(time_pts, I_hist, t, 'previous', 'extrap')/N;
    H = interp1(time_pts, H_hist, t, 'previous', 'extrap')/N;
    R = interp1(time_pts, R_hist, t, 'previous', 'extrap')/N;
    D = interp1(time_pts, D_hist, t, 'previous', 'extrap')/N;

    % clamp to [0,1]
    S=max(S,0); I=max(I,0); H=max(H,0); R=max(R,0); D=max(D,0);
end

% Deterministic ODE (fractions) for overlay
function det = solve_deterministic_sihrs(P)
    tspan = [0, P.T];
    y0 = [P.initial_s; P.initial_i; P.initial_h; P.initial_r; P.initial_d];

    function dydt = f(~, y)
        S = y(1); I = y(2); H = y(3); R = y(4); D = y(5);
        infection = P.beta * P.pSI * S * I;
        back_to_S = P.lambda * P.pRS * R;

        dS = -infection + back_to_S;
        dI =  infection - P.gamma * (P.pIH + P.pIR + P.pID) * I;
        dH =  P.gamma * P.pIH * I - P.alpha * (1 - P.pHH) * H;
        dR =  P.gamma * P.pIR * I + P.alpha * P.pHR * H - back_to_S;
        dD =  P.gamma * P.pID * I + P.alpha * P.pHD * H;
        dydt = [dS; dI; dH; dR; dD];
    end

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [T,Y] = ode45(@f, tspan, y0, opts);

    det.T=T; det.S=Y(:,1); det.I=Y(:,2); det.H=Y(:,3); det.R=Y(:,4); det.D=Y(:,5);
end

% ---------- PLOTS (I/H overlays & random paths) ----------
function plot_overlay_IH_first_path(results, det)
    for k = 1:numel(results)
        N = results{k}.N; t = results{k}.t;

        % I(t) overlay
        figure('Position',[100,100,1000,700]); hold on;
        plot(t, results{k}.I1, '-', 'LineWidth', 2.5, 'Color', [1 0.5 0.5], 'DisplayName', sprintf('i_N(t), N=%d', N));
        plot(det.T, det.I, '--', 'Color', [0.49 0.18 0.56], 'LineWidth', 2.8, 'DisplayName', 'i(t) ODE');
        xlabel('Time (days)'); ylabel('I fraction'); grid on; legend('Location','northeast');
        set(gca,'FontSize',26);
        saveas(gcf, sprintf('I_vs_ODE_N%d.png', N));

        % H(t) overlay
        figure('Position',[100,100,1000,700]); hold on;
        plot(t, results{k}.H1, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('h_N(t), N=%d', N));
        plot(det.T, det.H, '--', 'LineWidth', 2.8, 'DisplayName', 'h(t) ODE');
        xlabel('Time (days)'); ylabel('H fraction'); grid on; legend('Location','northeast');
        set(gca,'FontSize',26);
        saveas(gcf, sprintf('H_vs_ODE_N%d.png', N));
    end
end

function plot_IH_random_paths(results)
    for k = 1:numel(results)
        N = results{k}.N; t = results{k}.t;

        % I(t) random paths
        figure('Position',[200,200,1000,650]); hold on;
        Iset = results{k}.I_paths_plot; P = size(Iset,1); cmap = lines(max(P,7));
        for p = 1:P
            plot(t, Iset(p,:), '-', 'LineWidth', 1.1, 'Color', cmap(mod(p-1,size(cmap,1))+1,:));
        end
        if P > 0, plot(t, mean(Iset,1), 'k-', 'LineWidth', 2.4, 'DisplayName', 'shown-path mean'); end
        title(sprintf('Infected fraction I(t) — %d random paths (N=%d)', P, N));
        xlabel('Time (days)'); ylabel('I(t)'); grid on; legend('Location','northeastoutside');
        set(gca,'FontSize',22);
        saveas(gcf, sprintf('I_paths_random_with_mean_N%d.png', N));

        % H(t) random paths
        figure('Position',[220,220,1000,650]); hold on;
        Hset = results{k}.H_paths_plot; P = size(Hset,1); cmap = lines(max(P,7));
        for p = 1:P
            plot(t, Hset(p,:), '-', 'LineWidth', 1.1, 'Color', cmap(mod(p-1,size(cmap,1))+1,:));
        end
        if P > 0, plot(t, mean(Hset,1), 'k-', 'LineWidth', 2.4, 'DisplayName', 'shown-path mean'); end
        title(sprintf('Hospitalized fraction H(t) — %d random paths (N=%d)', P, N));
        xlabel('Time (days)'); ylabel('H(t)'); grid on; legend('Location','northeastoutside');
        set(gca,'FontSize',22);
        saveas(gcf, sprintf('H_paths_random_with_mean_N%d.png', N));
    end
end

% ---------- PLOTS (Qℓ linear; RQℓ symlog) ----------
function plot_Q_RQ_all(results, params)
    linthresh = params.symlog_linthresh;
    names = {'S','I','H','R','D'};
    for k = 1:numel(results)
        N = results{k}.N; t = results{k}.t;

        % ----- Qℓ mean (linear) -----
        for ell = 1:5
            figure('Position',[160+ell*10,160+ell*10,900,600]);
            plot(t, squeeze(results{k}.Q_mean(ell,:)), 'LineWidth', 2.8);
            xlabel('Time (days)'); ylabel(sprintf('Q_{%d}^{mean}(t)', ell)); grid on;
            legend(sprintf('Mean Q_{%d}, N=%d', ell, N),'Location','northwest');
            title(sprintf('Q_{%d} (%s)', ell, names{ell}));
            set(gca,'FontSize',22);
            saveas(gcf, sprintf('Q%d_mean_N%d.png', ell, N));
        end

        % ----- Qℓ per-path (linear) -----
        for ell = 1:5
            figure('Position',[170+ell*10,170+ell*10,900,600]); hold on;
            Qp = squeeze(results{k}.Q_paths(ell,:,:));  % (P x T)
            Pn = size(Qp,1);
            cmap = lines(max(Pn,7));
            for p = 1:Pn
                plot(t, Qp(p,:), '-', 'LineWidth', 1.3, 'Color', cmap(mod(p-1,size(cmap,1))+1,:));
            end
            xlabel('Time (days)'); ylabel(sprintf('Q_{%d}^{(p)}(t)', ell)); grid on;
            legend(arrayfun(@(p)sprintf('path %d',p),1:Pn,'UniformOutput',false),'Location','northeastoutside');
            title(sprintf('Q_{%d} Paths (%s), N=%d', ell, names{ell}, N));
            set(gca,'FontSize',20);
            saveas(gcf, sprintf('Q%d_paths_N%d.png', ell, N));
        end

        % ----- RQℓ mean (symlog) -----
        for ell = 1:5
            figure('Position',[180+ell*10,180+ell*10,900,600]);
            y = squeeze(results{k}.RQ_mean(ell,:)); y(~isfinite(y)) = NaN;
            plot(t, symlog_transform(y, linthresh), 'LineWidth', 2.8);
            xlabel('Time (days)'); ylabel(sprintf('symlog_{10}(RQ_{%d}^{mean}), a=%g', ell, linthresh)); grid on;
            legend(sprintf('Mean RQ_{%d}, N=%d', ell, N),'Location','northwest');
            title(sprintf('RQ_{%d} (%s) symlog', ell, names{ell}));
            set(gca,'FontSize',22);
            saveas(gcf, sprintf('RQ%d_mean_symlog_N%d.png', ell, N));
        end

        % ----- RQℓ per-path (symlog) -----
        for ell = 1:5
            figure('Position',[190+ell*10,190+ell*10,900,600]); hold on;
            RQp = squeeze(results{k}.RQ_paths(ell,:,:)); % (P x T)
            Pn = size(RQp,1);
            cmap = lines(max(Pn,7));
            for p = 1:Pn
                yy = RQp(p,:); yy(~isfinite(yy)) = NaN;
                plot(t, symlog_transform(yy, linthresh), '-', 'LineWidth', 1.3, 'Color', cmap(mod(p-1,size(cmap,1))+1,:));
            end
            xlabel('Time (days)'); ylabel(sprintf('symlog_{10}(RQ_{%d}^{(p)}), a=%g', ell, linthresh)); grid on;
            legend(arrayfun(@(p)sprintf('path %d',p),1:Pn,'UniformOutput',false),'Location','northeastoutside');
            title(sprintf('RQ_{%d} Paths (%s) symlog, N=%d', ell, names{ell}, N));
            set(gca,'FontSize',20);
            saveas(gcf, sprintf('RQ%d_paths_symlog_N%d.png', ell, N));
        end
    end
end

% ---------- symlog helper (Option 1) ----------
function yplot = symlog_transform(y, linthresh)
    % Symmetric log base-10 transform with linear region around zero.
    % yplot = sign(y) * log10(1 + |y| / linthresh)
    yplot = sign(y) .* log10(1 + abs(y) / linthresh);
end
