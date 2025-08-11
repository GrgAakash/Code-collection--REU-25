# SIHRS with Death Model - ODE System

## Mathematical Model

The SIHRS epidemic model is described by the following system of ordinary differential equations:

```
ds/dt = -β·S·I·pSI + λ·R·pRS
di/dt = β·S·I·pSI - γ·I·(1-pII)  
dh/dt = γ·I·pIH - α·H·(1-pHH)
dr/dt = γ·I·pIR + α·H·pHR - λ·R·pRS
dd/dt = γ·I·pID + α·H·pHD
```

## Parameter Definitions

### Rate Parameters
- **β** (beta): Infection rate - controls how quickly susceptible individuals become infected
- **γ** (gamma): I transition rate - controls how quickly infected individuals transition to other states
- **α** (alpha): H transition rate - controls how quickly hospitalized individuals transition to other states  
- **λ** (lambda): R transition rate - controls waning immunity (recovered individuals become susceptible again)

### Transition Probabilities
- **pSI**: Probability of S → I transition (susceptible to infected)
- **pII**: Probability of I → I transition (infected stays infected)
- **pIH**: Probability of I → H transition (infected to hospitalized)
- **pIR**: Probability of I → R transition (infected to recovered)
- **pID**: Probability of I → D transition (infected to dead)
- **pHH**: Probability of H → H transition (hospitalized stays hospitalized)
- **pHR**: Probability of H → R transition (hospitalized to recovered)
- **pHD**: Probability of H → D transition (hospitalized to dead)
- **pRR**: Probability of R → R transition (recovered stays recovered)
- **pRS**: Probability of R → S transition (recovered to susceptible)

### Model Constraints
- All transition probabilities must sum to 1 for each compartment
- All rates must be positive
- Initial conditions must sum to 1 (proportions)
