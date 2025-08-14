# JUSTIFICATION FOR SIHRS MODEL PARAMETERS

## EXPLANATION FOR R₀ (Basic Reproduction Number)

The most important parameter to get correct is the basic reproduction number:

```
R₀ = (β × pSI) / (γ × (1-pII))
```

From the paper [1], we have R₀ statistics for 1,151 counties:

### R₀ (Cases) Statistics:
- **Number of counties:** 1,151
- **Mean R₀:** 1.791268
- **Median R₀:** 1.663525
- **Standard Deviation:** 0.718795

### Range and Percentiles:
- **Minimum:** 0.381350
- **Maximum:** 12.436950
- **25th percentile:** 1.346918
- **75th percentile:** 2.111341
- **95th percentile:** 3.004858

**R₀ = 2.12** was chosen as it represents the 75th percentile, giving us a realistic but somewhat high transmission scenario that captures the variability seen in real-world COVID-19 data.

## EXPLANATION FOR TRANSITION PROBABILITIES

### pIR = 0.959
**Justification:** For COVID-19, a large majority of infected people recovered without requiring hospitalization. This reflects the fact that most cases were mild to moderate.

### pHR = 0.9882
**Justification:** Even among hospitalized COVID-19 patients, the majority eventually recovered. This accounts for the effectiveness of medical care and the body's natural recovery mechanisms.

### pRS = 0.98
**Justification:** The majority of the population loses immunity over time and becomes susceptible again. This reflects the waning immunity observed with COVID-19, where natural and vaccine-induced immunity decreases over several months.

## EXPLANATION FOR TRANSITION RATES

### λ (Lambda) = 0.0083
**Justification:** For COVID-19, the immunity period was on average 4 to 6 months. This rate parameter controls how quickly recovered individuals transition back to susceptible status, reflecting the natural decay of immunity.

## PARAMETER TUNING APPROACH

The parameters β, γ, pSI, and pII are tuned to achieve the target **R₀ = 2.12** within 4-5% accuracy. This ensures the model produces realistic epidemic dynamics while maintaining mathematical consistency.

## EPIDEMIC TIMESCALE CONSIDERATIONS

While multiple combinations of the parameters β, γ, pSI, and pII can yield the same basic reproductive number, R₀, the specific values are not arbitrary. They critically define the epidemic's timescale. For example, the parameter γ determines the average infectious period (as 1/γ).

A model with high rates (e.g., β=2.12, γ=1) implies a very rapid transmission and recovery cycle. In stochastic simulations, this can cause the infected population to "burn out" and drop to zero prematurely.

Therefore, to model more realistic dynamics that unfold over weeks or months, it is preferable to select parameter values that reflect a longer infectious period (e.g., β=0.212, γ=0.1), as this produces a more stable and prolonged epidemic curve.

## INITIAL CONDITIONS

The initial conditions (s0, i0, h0, r0, d0) are chosen for demonstration purposes to make plots more visually interpretable, though real-world scenarios typically have s0 ≈ 0.999 and i0 ≈ 0.001.

## REFERENCES

[1] Karla Therese L. Sy, Laura F. White, and Brooke E. Nichols, "Population density and basic reproductive number of covid-19 across united states counties," *PLOS ONE* 16, no. 4 (2021): 1–11.
