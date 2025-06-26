Demo at : https://grgaakash.github.io/SIHR%20Stochastic%20vs%20ODE/index.html

## Mathematical Background

The SIHR model extends the classic SIR model by adding a hospitalized compartment:

```
dS/dt = -p₁βSI
dI/dt = p₁βSI - p₂γI
dH/dt = p₂pₕγI - p₃αH
dR/dt = p₂(1-pₕ)γI + p₃αH
```

Where:
- S: Susceptible individuals
- I: Infected individuals  
- H: Hospitalized individuals
- R: Recovered individuals

