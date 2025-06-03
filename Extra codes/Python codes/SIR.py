import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import odeint
from typing import List, Tuple, Dict
from dataclasses import dataclass


@dataclass
class ModelParameters:
    """Parameters for the SIR model simulation."""
    beta: float = 1.1  # Infection rate
    gamma: float = 1.0  # Recovery rate
    p1: float = 0.5  # Probability of infection event occurring
    p2: float = 0.5  # Probability of recovery event occurring
    tend: float = 30  # Total time to simulate


def run_stochastic_simulation(
        N_init: int,
        params: ModelParameters
) -> Tuple[List[float], List[int], List[int], List[int]]:
    """
    Run a stochastic SIR simulation using the Gillespie algorithm.

    Args:
        N_init: Initial population size
        params: Model parameters

    Returns:
        Tuple containing (time points, S values, I values, R values)
    """
    # Initial conditions (96% S, 4% I, 0% R)
    S = [int(N_init * 0.96)]
    I = [int(N_init * 0.04)]
    R = [0]
    t = [0]

    while t[-1] < params.tend and (S[-1] + I[-1] >= 1):
        N = S[-1] + I[-1] + R[-1]
        props = [params.beta * I[-1] * S[-1] / N, params.gamma * I[-1]]
        prop_sum = sum(props)

        if prop_sum == 0:
            break

        tau = np.random.exponential(scale=1 / prop_sum)
        t.append(t[-1] + tau)

        rand = random.uniform(0, 1)

        if rand * prop_sum <= props[0]:  # Infection event
            if random.random() < params.p1:
                S.append(S[-1] - 1)
                I.append(I[-1] + 1)
                R.append(R[-1])
            else:
                S.append(S[-1])
                I.append(I[-1])
                R.append(R[-1])
        else:  # Recovery event
            if random.random() < params.p2:
                S.append(S[-1])
                I.append(I[-1] - 1)
                R.append(R[-1] + 1)
            else:
                S.append(S[-1])
                I.append(I[-1])
                R.append(R[-1])

    return t, S, I, R


def sim_deterministic(variables: List[float], t: float, params: List[float]) -> List[float]:
    """
    Deterministic SIR model differential equations.

    Args:
        variables: Current state [S, I, R]
        t: Current time point
        params: [beta, gamma]

    Returns:
        List of derivatives [dS/dt, dI/dt, dR/dt]
    """
    S, I, R = variables
    N = S + I + R
    beta, gamma = params

    dSdt = -beta * I * S / N
    dIdt = beta * I * S / N - gamma * I
    dRdt = gamma * I

    return [dSdt, dIdt, dRdt]


def plot_sir_comparison(params: ModelParameters):
    """
    Create and save comparison plots of stochastic and deterministic SIR models.
    """
    N_values = [316, 3162, 10000]
    colors = ['blue', 'orange', 'green']

    f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10, 12), sharex=True)
    peak_I_values = []
    peak_times = []

    # Stochastic simulations
    for N_init, color in zip(N_values, colors):
        t, S, I, R = run_stochastic_simulation(N_init, params)

        # Plot stochastic data
        ax1.plot(t, S, label=f'S (N={N_init}, Stochastic)', color=color, alpha=0.5)
        ax2.plot(t, I, label=f'I (N={N_init}, Stochastic)', color=color, alpha=0.5)
        ax3.plot(t, R, label=f'R (N={N_init}, Stochastic)', color=color, alpha=0.5)

        # Calculate peak information
        max_I = max(I)
        peak_idx = I.index(max_I)
        peak_I_values.append(max_I / N_init)
        peak_times.append(t[peak_idx])

    # Deterministic simulation
    t_det = np.linspace(0, params.tend, num=1000)
    y0 = [10000 * 0.96, 10000 * 0.04, 0]
    model_params = [params.beta, params.gamma]
    y = odeint(sim_deterministic, y0, t_det, args=(model_params,))

    # Plot deterministic results
    ax1.plot(t_det, y[:, 0], label='S (N=10000, Deterministic)', color='black', linestyle='--')
    ax2.plot(t_det, y[:, 1], label='I (N=10000, Deterministic)', color='black', linestyle='--')
    ax3.plot(t_det, y[:, 2], label='R (N=10000, Deterministic)', color='black', linestyle='--')

    # Add labels and formatting
    initial_labels = ["Initial: 0.96", "Initial: 0.04", "Initial: 0.00"]
    peak_labels = [f"N={N}: Peak: {peak:.3f} at t={t_i:.1f}"
                   for N, peak, t_i in zip(N_values, peak_I_values, peak_times)]

    for i, (label, ax) in enumerate(zip(initial_labels, [ax1, ax2, ax3])):
        ax.text(0.02, 0.9 - i * 0.1, label, transform=ax.transAxes, color='black', fontsize=10)

    for i, (label, ax) in enumerate(zip(peak_labels, [ax1, ax2, ax3])):
        ax.text(0.98, 0.9 - i * 0.1, label, transform=ax.transAxes, ha='right',
                color='black', fontsize=10)

    ax1.set_ylabel("S")
    ax2.set_ylabel("I")
    ax3.set_ylabel("R")
    ax3.set_xlabel("Time")
    ax1.legend()
    ax2.legend()
    ax3.legend()

    plt.tight_layout()
    plt.savefig('sir_simulation.png')
    plt.show()


if __name__ == "__main__":
    params = ModelParameters()
    plot_sir_comparison(params)