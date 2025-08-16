# SIHRS Variance Analysis - Complete Technical Documentation

## Overview
This document explains the complete methodology behind the SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible) model variance analysis, including how plots are generated and what they represent.

## Mathematical Background

### The New Conjecture
The analysis validates the conjecture:
```
E[V_N^(l)(t)] ∝ (1/N) × population_proportions
```

Where:
- **E[V_N^(l)(t)]**: Expected variance of compartment l at time t
- **N**: Population size
- **population_proportions**: s(t), i(t), h(t), r(t), d(t) from ODE solution

### Theoretical Variance Formulas
```
V_N^(1)(t) = (1/N)(β s(t) i(t) pSI + λ r(t) pRS)           # Susceptible
V_N^(2)(t) = (1/N)(γ i(t) (pIH + pIR + pID) + β s(t) i(t) pSI)  # Infected
V_N^(3)(t) = (1/N)(γ i(t) pIH + α h(t) (pHR + pHD))        # Hospitalized
V_N^(4)(t) = (1/N)(γ i(t) pIR + α h(t) pHR + λ r(t) pRS)   # Recovered
V_N^(5)(t) = (1/N)(γ i(t) pID + α h(t) pHD)                # Dead
```

## Data Generation Process

### 1. Simulation Structure
```
For each N value (316, 3162, 10000):
├── Batch 1: 15 trials
├── Batch 2: 15 trials
├── Batch 3: 15 trials
├── Batch 4: 15 trials
└── Batch 5: 15 trials

Total: 5 batches × 15 trials = 75 simulations per N
```

### 2. Data Storage Structure
```matlab
% 4D arrays: [N_values × time_points × trials × batches]
S_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
I_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
H_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
R_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);
D_var_all = zeros(length(N_values), length(t_grid), num_trials, num_batches);

% 3D arrays: [N_values × time_points × batches] - averaged across trials
S_expected_var = zeros(length(N_values), length(t_grid), num_batches);
I_expected_var = zeros(length(N_values), length(t_grid), num_batches);
H_expected_var = zeros(length(N_values), length(t_grid), num_batches);
R_expected_var = zeros(length(N_values), length(t_grid), num_batches);
D_expected_var = zeros(length(N_values), length(t_grid), num_batches);
```

### 3. Time Grid
```matlab
t_grid = linspace(0, params.tmax, 1000);  % 1000 time points from 0 to 700 days
```

## Data Processing Steps

### Step 1: Individual Trial Simulation
```matlab
for trial = 1:num_trials  % 15 trials per batch
    result = sihrs_agent_model(N, params);
    
    % Interpolate variance data onto common time grid
    S_var_all(idx, :, trial, batch) = interp1(result.T, result.vs, t_grid, 'linear', 0);
    I_var_all(idx, :, trial, batch) = interp1(result.T, result.vi, t_grid, 'linear', 0);
    H_var_all(idx, :, trial, batch) = interp1(result.T, result.vh, t_grid, 'linear', 0);
    R_var_all(idx, :, trial, batch) = interp1(result.T, result.vr, t_grid, 'linear', 0);
    D_var_all(idx, :, trial, batch) = interp1(result.T, result.vd, t_grid, 'linear', 0);
end
```

### Step 2: Batch Averaging (Expected Variance)
```matlab
% Calculate expected variance for this batch (average across 15 trials)
S_expected_var(idx, :, batch) = mean(S_var_all(idx, :, :, batch), 3);
I_expected_var(idx, :, batch) = mean(I_var_all(idx, :, :, batch), 3);
H_expected_var(idx, :, batch) = mean(H_var_all(idx, :, :, batch), 3);
R_expected_var(idx, :, batch) = mean(R_var_all(idx, :, :, batch), 3);
D_expected_var(idx, :, batch) = mean(D_var_all(idx, :, :, batch), 3);
```

**What this means**: For each time point, take the average of variance across 15 trials to get one "expected variance curve" per batch.

### Step 3: Theoretical Variance Calculation
```matlab
% Solve deterministic ODE for theoretical comparison
det_result = solve_deterministic_sihrs(params);

% Interpolate ODE solution to our time grid
s_ode = interp1(det_result.T, det_result.S_prop, t_grid, 'linear', 0);
i_ode = interp1(det_result.T, det_result.I_prop, t_grid, 'linear', 0);
h_ode = interp1(det_result.T, det_result.H_prop, t_grid, 'linear', 0);
r_ode = interp1(det_result.T, det_result.R_prop, t_grid, 'linear', 0);
d_ode = interp1(det_result.T, det_result.D_prop, t_grid, 'linear', 0);

% Calculate theoretical variance using mathematical formulas
v_s_theory{idx} = (1/N) * (params.beta * s_ode .* i_ode * params.pSI + params.lambda * r_ode * params.pRS);
v_i_theory{idx} = (1/N) * (params.gamma * i_ode * (params.pIH + params.pIR + params.pID) + params.beta * s_ode .* i_ode * params.pSI);
v_h_theory{idx} = (1/N) * (params.gamma * i_ode * params.pIH + params.alpha * h_ode * (params.pHR + params.pHD));
v_r_theory{idx} = (1/N) * (params.gamma * i_ode * params.pIR + params.alpha * h_ode * params.pHR + params.lambda * r_ode * params.pRS);
v_d_theory{idx} = (1/N) * (params.gamma * i_ode * params.pID + params.alpha * h_ode * params.pHD);
```

## Plotting Methodology

### 1. Time Window Configuration
```matlab
dt = 2; % Half width for time intervals around each midpoint
midpoints = 5:dt*2:650; % Time points: 5, 9, 13, 17, ..., 649
```

**What this means**:
- **dt = 2**: Creates 4-day time windows (t-2 to t+2)
- **midpoints**: Analysis points every 4 days from t=5 to t=649
- **Total analysis points**: 162 time windows

### 2. Time Window Processing
For each midpoint (e.g., t=25):
```matlab
t_mid = midpoints(t_idx);    % Current midpoint (e.g., t=25)
t_min = t_mid - dt;          % Window start (e.g., t=23)
t_max = t_mid + dt;          % Window end (e.g., t=27)

in_window_indices = t_grid >= t_min & t_grid < t_max;
```

**Time window [23, 27)**: Collects data from time points 23, 24, 25, 26

### 3. Error Bar Calculation
```matlab
% Get batch values in this time window
batch_vals_in_window = batch_data(in_window_indices, :);
% Shape: [~4 time points × 5 batches]

% Calculate min and max across batches in this time window
batch_mins = min(batch_vals_in_window, [], 1);  % Min of each batch
batch_maxs = max(batch_vals_in_window, [], 1);  % Max of each batch

% Calculate final expected variance bounds (average across 5 batches)
y_min_avg = mean(batch_mins(~isnan(batch_mins)));  % Final E[min V_N^(l)(t)]
y_max_avg = mean(batch_maxs(~isnan(batch_maxs)));  % Final E[max V_N^(l)(t)]
```

**Step-by-step breakdown**:
1. **Extract data** from time window [23, 27) for all 5 batches
2. **For each batch**: Find min and max variance in the 4-day window
3. **Across batches**: Average the 5 min values and 5 max values
4. **Result**: One error bar showing expected variance range

### 4. Plotting Elements
```matlab
% Plot error bar showing final expected variance range
plot([t_mid, t_mid], [y_min_avg, y_max_avg], 'Color', colors(n_idx,:), 'LineWidth', 1.5);

% Plot theoretical expected variance value
[~, time_index] = min(abs(t_grid - t_mid));
y_theory = theory_data(time_index); 
plot(t_mid, y_theory, '_', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
```

**What gets plotted**:
- **Vertical line**: Error bar from y_min_avg to y_max_avg
- **Horizontal marker**: Theoretical variance value at that time point

## Output Files

### 1. Individual Compartment Figures
- **Susceptible**: `SIHRS_Susceptible_final_expected_variance.png`
- **Infected**: `SIHRS_Infected_final_expected_variance.png`
- **Hospitalized**: `SIHRS_Hospitalized_final_expected_variance.png`
- **Recovered**: `SIHRS_Recovered_final_expected_variance.png`
- **Dead**: `SIHRS_Dead_final_expected_variance.png`

### 2. Figure Layout
Each figure contains 3 subplots arranged horizontally:
- **Left**: N = 316 (small population, high stochasticity)
- **Middle**: N = 3162 (medium population, moderate stochasticity)
- **Right**: N = 10000 (large population, low stochasticity)

## Data Flow Summary

```
Stochastic Simulation (75 runs per N)
           ↓
    Variance Calculation
           ↓
    Trial Averaging (15 trials → 1 batch result)
           ↓
    5 Batch Results
           ↓
    Time Window Analysis (4-day windows)
           ↓
    Error Bar Calculation (min/max across batches)
           ↓
    Final Plot (Error bars + Theoretical markers)
```

## Key Parameters

### Simulation Parameters
- **Population sizes**: N = [316, 3162, 10000]
- **Trials per batch**: 15
- **Number of batches**: 5
- **Total runs per N**: 75
- **Simulation time**: 700 days
- **Time grid points**: 1000

### Plotting Parameters
- **Time window width**: 4 days (dt = 2)
- **Analysis points**: 162 (every 4 days from t=5 to t=649)
- **Figure dimensions**: 1600 × 800 pixels
- **Marker size**: 12 (horizontal lines)
- **Line width**: 1.5 (error bars), 2 (theoretical markers)

## Expected Results

### 1. Variance Scaling
- **N = 316**: High variance, high stochasticity, potential extinction
- **N = 3162**: Medium variance, moderate stochasticity, good agreement
- **N = 10000**: Low variance, low stochasticity, excellent agreement

### 2. Error Bar Behavior
- **Height decreases** as N increases (1/N scaling)
- **Width decreases** as N increases (less batch-to-batch variability)
- **Agreement with theory** improves as N increases

### 3. Theoretical Validation
- **Small N**: Simulated variance may deviate from theory due to stochastic effects
- **Large N**: Simulated variance should closely match theoretical predictions
- **Overall trend**: Convergence to theoretical 1/N scaling relationship

## Troubleshooting

### Common Issues
1. **Error bars too small**: Check if batches are too similar (increase batch diversity)
2. **Theoretical markers too large**: Reduce MarkerSize parameter
3. **Plots too cluttered**: Increase dt for wider time windows
4. **Memory issues**: Reduce number of trials or batches

### Performance Optimization
- **Reduce trials**: 10 instead of 15 (faster, less statistical power)
- **Reduce batches**: 3 instead of 5 (faster, less robust)
- **Increase dt**: 5 instead of 2 (fewer analysis points, faster)

## Mathematical Validation

### Expected Behavior
1. **1/N scaling**: Variance should decrease proportionally with population size
2. **Time evolution**: Variance should follow disease dynamics
3. **Stochastic convergence**: Large N should approach deterministic behavior
4. **Batch consistency**: Error bars should shrink with increasing N

### Quality Metrics
- **Theoretical agreement**: How well simulated results match theory
- **Stochastic scaling**: How well 1/N relationship holds
- **Batch consistency**: How similar results are across batches
- **Temporal stability**: How stable variance is over time

This analysis provides a comprehensive validation of the SIHRS model's stochastic properties and their relationship to deterministic predictions.
