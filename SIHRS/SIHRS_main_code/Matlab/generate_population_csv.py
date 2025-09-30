#!/usr/bin/env python3
"""
Generate CSV file with SIHRS population data and D_N values
Python equivalent of generate_population_csv.m
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random

def generate_sihrs_csv():
    """Main function to run SIHRS simulation and generate CSV"""
    # Set random seed for reproducibility
    np.random.seed(1)
    random.seed(1)
    
    # Parameters matching generate_population_csv.m
    params = {
        'beta': 0.12,           # Infection rate (β > 0)
        'pSI': 1.0,             # Infection probability (S to I)
        'pII': 0.0,             # probability of I to I (stay infected)
        'pIH': 0.04,            # probability of I to H
        'pIR': 0.959,           # probability of I to R
        'pID': 0.001,           # probability of I to D
        'pHH': 0.01,            # probability of H to H (stay hospitalized)
        'pHR': 0.9882,          # probability of H to R
        'pHD': 0.0018,          # probability of H to D
        'pRR': 0.02,            # probability of R to R (stay recovered)
        'pRS': 0.98,            # probability of R to S
        'gamma': 0.1,           # Adjusted for exact critical hit at S = 142/300
        'alpha': 0.1,           # Hospitalized transition rate (α > 0)
        'lambda': 0.0083,       # Recovered to susceptible rate (Λ > 0) immunity period of 4 months
        'T': 300,               # Total simulation time
        'dt': 0.01,             # Finer time step for better blow-up capture
        'N_values': [1600], 
        'initial_s': 0.96,      # Initial susceptible fraction
        'initial_i': 0.04,      # Initial infected fraction
        'initial_h': 0,         # Initial hospitalized fraction
        'initial_r': 0,         # Initial recovered fraction
        'initial_d': 0,         # Initial dead fraction
        'n_runs': 1             # Number of stochastic runs
    }
    
    # Calculate R0 from parameters
    R0 = params['pSI'] * params['beta'] / (params['gamma'] * (1 - params['pII']))
    print(f'Calculated R₀ = {R0:.4f} from parameters')
    
    # Precompute time vector
    t = np.arange(0, params['T'] + params['dt'], params['dt'])
    
    # Run simulation for each population size
    for idx, N in enumerate(params['N_values']):
        print(f'Running {params["n_runs"]} simulations for N = {N}, R0 = {R0:.2f}...')
        
        # Run multiple stochastic simulations
        result = run_multiple_gillespie_renormalized(N, params, t)
        
        # Calculate D_N values using pure mathematical formulas
        result = calculate_dn_values(result, N, params)
        
        # Save population data and D_N values to CSV file
        population_data = pd.DataFrame({
            'Time': t,
            'S_mean': result['S_mean'],
            'I_mean': result['I_mean'],
            'H_mean': result['H_mean'],
            'R_mean': result['R_mean'],
            'D_mean': result['D_mean'],
            'D1_Susceptible': result['D1'],
            'D2_Infected': result['D2'],
            'D3_Hospitalized': result['D3'],
            'D4_Recovered': result['D4'],
            'D5_Dead': result['D5']
        })
        
        csv_filename = f'SIHRS_population_and_DN_data_N{N}_R0_{R0:.2f}.csv'
        population_data.to_csv(csv_filename, index=False)
        print(f'Population data and D_N values saved to: {csv_filename}')
        
        # Print summary statistics
        print(f'\nSummary for N = {N}:')
        print(f'  Susceptible: Mean D_N = {np.mean(result["D1"]):.6f}, Max D_N = {np.max(result["D1"]):.6f}')
        print(f'  Infected: Mean D_N = {np.mean(result["D2"]):.6f}, Max D_N = {np.max(result["D2"]):.6f}')
        print(f'  Hospitalized: Mean D_N = {np.mean(result["D3"]):.6f}, Max D_N = {np.max(result["D3"]):.6f}')
        print(f'  Recovered: Mean D_N = {np.mean(result["D4"]):.6f}, Max D_N = {np.max(result["D4"]):.6f}')
        print(f'  Dead: Mean D_N = {np.mean(result["D5"]):.6f}, Max D_N = {np.max(result["D5"]):.6f}')
        
        # Check for zero infected population
        min_infected = np.min(result['I_mean'])
        zero_infected_indices = np.where(result['I_mean'] == 0)[0]
        print(f'\nInfected population analysis:')
        print(f'  Minimum infected proportion: {min_infected:.6f}')
        print(f'  Number of time points with zero infected: {len(zero_infected_indices)}')
        if len(zero_infected_indices) > 0:
            print(f'  First time with zero infected: t = {t[zero_infected_indices[0]]:.2f}')
            print(f'  Last time with zero infected: t = {t[zero_infected_indices[-1]]:.2f}')

def run_multiple_gillespie_renormalized(N, params, t):
    """Run multiple Gillespie simulations and aggregate results"""
    S_all = np.zeros((params['n_runs'], len(t)))
    I_all = np.zeros((params['n_runs'], len(t)))
    H_all = np.zeros((params['n_runs'], len(t)))
    R_all = np.zeros((params['n_runs'], len(t)))
    D_all = np.zeros((params['n_runs'], len(t)))
    
    for run in range(params['n_runs']):
        S_hist, I_hist, H_hist, R_hist, D_hist, time_pts = gillespie_sim(N, params)
        S_interp, I_interp, H_interp, R_interp, D_interp = interpolate_results(
            S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N)
        S_all[run, :] = S_interp
        I_all[run, :] = I_interp
        H_all[run, :] = H_interp
        R_all[run, :] = R_interp
        D_all[run, :] = D_interp
    
    # Compute mean proportions
    result = {
        'S_mean': np.mean(S_all, axis=0),
        'I_mean': np.mean(I_all, axis=0),
        'H_mean': np.mean(H_all, axis=0),
        'R_mean': np.mean(R_all, axis=0),
        'D_mean': np.mean(D_all, axis=0)
    }
    return result

def calculate_dn_values(result, N, params):
    """Calculate D_N values using pure mathematical formulas - no safety checks!"""
    
    # D_N^(1)(t) for Susceptible (S)
    # (D_N^(1))^2 = (1/N) * (β_SI * i_N(t)/s_N(t) + λ_RS * r_N(t)/(s_N(t))^2)
    result['D1_squared'] = (1/N) * (params['beta'] * params['pSI'] * result['I_mean'] / result['S_mean'] + 
                                   params['lambda'] * params['pRS'] * result['R_mean'] / (result['S_mean']**2))
    result['D1'] = np.sqrt(result['D1_squared'])
    
    # D_N^(2)(t) for Infected (I)
    # (D_N^(2))^2 = (1/N) * ((γ_IH + γ_IR + γ_ID)/i_N(t) + β_SI * s_N(t)/i_N(t))
    result['D2_squared'] = (1/N) * ((params['gamma'] * (params['pIH'] + params['pIR'] + params['pID'])) / result['I_mean'] + 
                                   params['beta'] * params['pSI'] * result['S_mean'] / result['I_mean'])
    result['D2'] = np.sqrt(result['D2_squared'])
    
    # D_N^(3)(t) for Hospitalized (H)
    # (D_N^(3))^2 = (1/N) * ((α_HR + α_HD)/h_N(t) + γ_IH * i_N(t)/(h_N(t))^2)
    result['D3_squared'] = (1/N) * ((params['alpha'] * (params['pHR'] + params['pHD'])) / result['H_mean'] + 
                                   params['gamma'] * params['pIH'] * result['I_mean'] / (result['H_mean']**2))
    result['D3'] = np.sqrt(result['D3_squared'])
    
    # D_N^(4)(t) for Recovered (R)
    # (D_N^(4))^2 = (1/N) * (λ_RS/r_N(t) + γ_IR * i_N(t)/(r_N(t))^2 + α_HR * h_N(t)/(r_N(t))^2)
    result['D4_squared'] = (1/N) * (params['lambda'] * params['pRS'] / result['R_mean'] + 
                                   params['gamma'] * params['pIR'] * result['I_mean'] / (result['R_mean']**2) + 
                                   params['alpha'] * params['pHR'] * result['H_mean'] / (result['R_mean']**2))
    result['D4'] = np.sqrt(result['D4_squared'])
    
    # D_N^(5)(t) for Dead (D)
    # V5 = (1/N) * (pID * γ * i_N(t) + pHD * α * h_N(t))
    # D_N^(5) = √(V5)/d_N(t)
    V5 = (1/N) * (params['pID'] * params['gamma'] * result['I_mean'] + params['pHD'] * params['alpha'] * result['H_mean'])
    result['D5'] = np.sqrt(V5) / result['D_mean']
    
    return result

def gillespie_sim(N, params):
    """Gillespie simulation for SIHRS with death"""
    # Initialize populations
    S = int(N * params['initial_s'])
    I = int(N * params['initial_i'])
    H = int(N * params['initial_h'])
    R = int(N * params['initial_r'])
    D = int(N * params['initial_d'])
    
    # Verify population conservation
    if abs(S + I + H + R + D - N) > 1:
        raise ValueError('Initial populations do not sum to N')
    
    # Preallocate arrays
    max_events = int(10 * params['T'] * (params['beta'] + params['gamma'] + params['alpha'] + params['lambda']) * N)
    S_hist = np.zeros(max_events)
    I_hist = np.zeros(max_events)
    H_hist = np.zeros(max_events)
    R_hist = np.zeros(max_events)
    D_hist = np.zeros(max_events)
    time_pts = np.zeros(max_events)
    
    S_hist[0] = S
    I_hist[0] = I
    H_hist[0] = H
    R_hist[0] = R
    D_hist[0] = D
    time_pts[0] = 0
    event_count = 1
    current_time = 0
    
    # Gillespie algorithm for SIHRS with death
    while current_time < params['T']:
        # Calculate event rates
        si_rate = (params['beta'] / N) * S * I * params['pSI']  # rate of S->I
        ir_rate = params['gamma'] * I * params['pIR']  # rate of I->R
        ih_rate = params['gamma'] * I * params['pIH']  # rate of I->H
        id_rate = params['gamma'] * I * params['pID']  # rate of I->D
        hr_rate = params['alpha'] * H * params['pHR']  # rate of H->R
        hd_rate = params['alpha'] * H * params['pHD']  # rate of H->D
        rs_rate = params['lambda'] * R * params['pRS']  # rate of R->S
        
        total_rate = si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate + rs_rate
        
        if total_rate == 0:
            break
        
        tau = -np.log(np.random.random()) / total_rate
        current_time = current_time + tau
        
        if current_time > params['T']:
            break
        
        # Determine which event occurs
        chance = np.random.random() * total_rate
        if chance < si_rate:
            # S->I transition
            if S > 0:
                S = S - 1
                I = I + 1
        elif chance < (si_rate + ir_rate):
            # I->R transition
            if I > 0:
                I = I - 1
                R = R + 1
        elif chance < (si_rate + ir_rate + ih_rate):
            # I->H transition
            if I > 0:
                I = I - 1
                H = H + 1
        elif chance < (si_rate + ir_rate + ih_rate + id_rate):
            # I->D transition
            if I > 0:
                I = I - 1
                D = D + 1
        elif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate):
            # H->R transition
            if H > 0:
                H = H - 1
                R = R + 1
        elif chance < (si_rate + ir_rate + ih_rate + id_rate + hr_rate + hd_rate):
            # H->D transition
            if H > 0:
                H = H - 1
                D = D + 1
        else:
            # R->S transition
            if R > 0:
                R = R - 1
                S = S + 1
        
        event_count = event_count + 1
        S_hist[event_count] = S
        I_hist[event_count] = I
        H_hist[event_count] = H
        R_hist[event_count] = R
        D_hist[event_count] = D
        time_pts[event_count] = current_time
    
    # Trim arrays
    S_hist = S_hist[:event_count]
    I_hist = I_hist[:event_count]
    H_hist = H_hist[:event_count]
    R_hist = R_hist[:event_count]
    D_hist = D_hist[:event_count]
    time_pts = time_pts[:event_count]
    
    # Check population conservation (allow small rounding errors)
    total_pop = S_hist + I_hist + H_hist + R_hist + D_hist
    if not np.allclose(total_pop, N, atol=1):
        print(f"Warning: Population conservation issue. Expected {N}, got {total_pop[-1]}")
        print(f"Final populations: S={S_hist[-1]}, I={I_hist[-1]}, H={H_hist[-1]}, R={R_hist[-1]}, D={D_hist[-1]}")
    
    return S_hist, I_hist, H_hist, R_hist, D_hist, time_pts

def interpolate_results(S_hist, I_hist, H_hist, R_hist, D_hist, time_pts, t, N):
    """Interpolate to fixed time grid using 'previous' method"""
    try:
        S_interp = np.interp(t, time_pts, S_hist) / N
        I_interp = np.interp(t, time_pts, I_hist) / N
        H_interp = np.interp(t, time_pts, H_hist) / N
        R_interp = np.interp(t, time_pts, R_hist) / N
        D_interp = np.interp(t, time_pts, D_hist) / N
        
        # Handle values beyond last event
        S_interp[t > np.max(time_pts)] = S_hist[-1] / N
        I_interp[t > np.max(time_pts)] = I_hist[-1] / N
        H_interp[t > np.max(time_pts)] = H_hist[-1] / N
        R_interp[t > np.max(time_pts)] = R_hist[-1] / N
        D_interp[t > np.max(time_pts)] = D_hist[-1] / N
        
        return S_interp, I_interp, H_interp, R_interp, D_interp
    except Exception as e:
        raise RuntimeError(f'Interpolation failed: {e}')

if __name__ == "__main__":
    generate_sihrs_csv()
