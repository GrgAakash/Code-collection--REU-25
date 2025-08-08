// SIHR-IPC Model Implementation
// Stochastic and Deterministic epidemic models

// Validate parameters for SIHRS
function validateParameters(params) {
    if ([params.beta, params.gamma, params.alpha, params.lambda].some(x => x <= 0)) {
        throw new Error('All rates must be positive');
    }
    if ([params.p_SI, params.p_IH, params.p_IR, params.p_ID, params.p_HR, params.p_HH, params.p_HD, params.p_RR, params.p_RS].some(x => x <= 0 || x > 1)) {
        throw new Error('All probabilities must be in (0,1]');
    }
    // Allow p_II to be 0 (no one remains infected)
    if (params.p_II < 0 || params.p_II > 1) {
        throw new Error('p_II must be in [0,1]');
    }
    if (Math.abs(params.p_IH + params.p_IR + params.p_II + params.p_ID - 1) > 1e-10) {
        throw new Error('p_IH + p_IR + p_II + p_ID must sum to 1');
    }
    if (Math.abs(params.p_HR + params.p_HH + params.p_HD - 1) > 1e-10) {
        throw new Error('p_HR + p_HH + p_HD must sum to 1');
    }
    if (Math.abs(params.p_RR + params.p_RS - 1) > 1e-10) {
        throw new Error('p_RR + p_RS must sum to 1');
    }
    if (Math.abs(params.s0 + params.i0 + params.h0 + params.r0 + params.d0 - 1) > 1e-10) {
        throw new Error('Initial conditions must sum to 1');
    }
}

// Stochastic SIHRS model using Gillespie algorithm
function sihrsStochasticModel(N, params, runNumber) {
    console.log(`🎲 Running stochastic simulation ${runNumber}/${totalRuns}`);

    // Initial conditions
    let s = Math.round(params.s0 * N);
    let i = Math.round(params.i0 * N);
    let h = Math.round(params.h0 * N);
    let r = Math.round(params.r0 * N);
    let d = Math.round(params.d0 * N);
    const total = s + i + h + r + d;
    if (total !== N) {
        s += N - total;
    }

    // Preallocate arrays like Julia implementation
    const maxEvents = Math.min(N * 10, 1000000); // Estimate max events, cap at 1M
    const T = new Array(maxEvents);
    const S_prop = new Array(maxEvents);
    const I_prop = new Array(maxEvents);
    const H_prop = new Array(maxEvents);
    const R_prop = new Array(maxEvents);
    const D_prop = new Array(maxEvents);
    
    // Initialize first values
    T[0] = 0;
    S_prop[0] = s / N;
    I_prop[0] = i / N;
    H_prop[0] = h / N;
    R_prop[0] = r / N;
    D_prop[0] = d / N;
    
    let t = 0;
    let eventCount = 0;

    // Main simulation loop
    while ((i > 0 || h > 0) && t < params.tmax && eventCount < maxEvents - 1)  {
        const nS = s, nI = i, nH = h, nR = r, nD = d;

        // Calculate event rates (SIHRS)
        const infectionRate = params.p_SI * params.beta * nS * nI / N;
        const toIIRate = params.gamma * params.p_II * nI; // I -> I (remains infected)
        const toHospitalRate = params.gamma * params.p_IH * nI; // I -> H
        const toRecoveredFromIRate = params.gamma * params.p_IR * nI; // I -> R
        const toDeadFromIRate = params.gamma * params.p_ID * nI; // I -> D
        const toRecoveredFromHRate = params.alpha * params.p_HR * nH; // H -> R
        const toDeadFromHRate = params.alpha * params.p_HD * nH; // H -> D
        const toSusceptibleFromR = params.p_RS * params.lambda * nR; // R -> S (reinfection)
        // H -> H (remains hospitalized) and R -> R (remains recovered) are not state changes
        const totalRate = infectionRate + toIIRate + toHospitalRate + toRecoveredFromIRate + toDeadFromIRate + toRecoveredFromHRate + toDeadFromHRate + toSusceptibleFromR;

        if (totalRate === 0) break;

        // Time to next event
        const dt = -Math.log(Math.random()) / totalRate;
        t += dt;

        if (t > params.tmax) {
            t = params.tmax;
            eventCount++;
            T[eventCount] = t;
            const total = s + i + h + r + d;
            S_prop[eventCount] = s / total;
            I_prop[eventCount] = i / total;
            H_prop[eventCount] = h / total;
            R_prop[eventCount] = r / total;
            D_prop[eventCount] = d / total;
            break;
        }

        // Determine which event occurs
        const chance = Math.random() * totalRate;
        let threshold = 0;
        if (chance < (threshold += infectionRate) && nS > 0) {
            s--; i++;
        } else if (chance < (threshold += toIIRate) && nI > 0) {
            // I -> I (remains infected), no state change
            // skip
        } else if (chance < (threshold += toHospitalRate) && nI > 0) {
            i--; h++;
        } else if (chance < (threshold += toRecoveredFromIRate) && nI > 0) {
            i--; r++;
        } else if (chance < (threshold += toDeadFromIRate) && nI > 0) {
            i--; d++;
        } else if (chance < (threshold += toRecoveredFromHRate) && nH > 0) {
            h--; r++;
        } else if (chance < (threshold += toDeadFromHRate) && nH > 0) {
            h--; d++;
        } else if (chance < (threshold += toSusceptibleFromR) && nR > 0) {
            r--; s++;
        }

        // Record every event using array indexing (Julia-style)
        eventCount++;
        T[eventCount] = t;
        const total = s + i + h + r + d;
        S_prop[eventCount] = s / total;
        I_prop[eventCount] = i / total;
        H_prop[eventCount] = h / total;
        R_prop[eventCount] = r / total;
        D_prop[eventCount] = d / total;
    }

    // Trim arrays to actual size (Julia-style)
    const finalT = T.slice(0, eventCount + 1);
    const finalS_prop = S_prop.slice(0, eventCount + 1);
    const finalI_prop = I_prop.slice(0, eventCount + 1);
    const finalH_prop = H_prop.slice(0, eventCount + 1);
    const finalR_prop = R_prop.slice(0, eventCount + 1);
    const finalD_prop = D_prop.slice(0, eventCount + 1);

    // Warn if we hit the event limit (rare for reasonable populations)
    if (eventCount >= maxEvents - 1) {
        console.warn(`⚠️ Simulation reached max events limit (${maxEvents}). Results may be incomplete.`);
    }

    return { 
        T: finalT, 
        S_prop: finalS_prop, 
        I_prop: finalI_prop, 
        H_prop: finalH_prop, 
        R_prop: finalR_prop, 
        D_prop: finalD_prop, 
        maxI: Math.max(...finalI_prop), 
        maxH: Math.max(...finalH_prop), 
        maxD: Math.max(...finalD_prop) 
    };
}

// Deterministic SIHRS model (Euler integration)
function solveDeterministicSIHRS(params) {
    console.log('📊 Solving deterministic model (SIHRS)');
    const dt = 0.1;
    const T = [];
    const Y = [[params.s0, params.i0, params.h0, params.r0, params.d0]];
    
    for (let t = 0; t <= params.tmax; t += dt) {
        T.push(t);
        if (t > 0) {
            const y = Y[Y.length - 1];
            const s = y[0], i = y[1], h = y[2], r = y[3], d = y[4];
            // ODE system (SIHRS)
            const dsdt = -params.beta * params.p_SI * s * i + params.p_RS * params.lambda * r;
            const didt = params.beta * params.p_SI * s * i - params.gamma * (1 - params.p_II) * i;
            const dhdt = params.p_IH * params.gamma * i - params.alpha * (1 - params.p_HH) * h;
            const drdt = params.p_IR * params.gamma * i + params.p_HR * params.alpha * h - params.p_RS * params.lambda * r;
            const ddot = params.p_ID * params.gamma * i + params.p_HD * params.alpha * h;
            // Euler method
            const yNext = [
                s + dt * dsdt,
                i + dt * didt,
                h + dt * dhdt,
                r + dt * drdt,
                d + dt * ddot
            ];
            Y.push(yNext);
        }
    }

    return { T, S_prop: Y.map(y => y[0]), I_prop: Y.map(y => y[1]), H_prop: Y.map(y => y[2]), R_prop: Y.map(y => y[3]), D_prop: Y.map(y => y[4]) };
}

// Function to integrate by Simpson's rule
function simpsonsRule(f, a, b, n) {
    if (n % 2 === 1) n++; // Ensure even number of intervals
    const h = (b - a) / n;
    let sum = f(a) + f(b);

    for (let i = 1; i < n; i++) {
      const x = a + i * h;
      sum += (i % 2 === 0 ? 2 : 4) * f(x);
    }

    return (h / 3) * sum;
}

// Function to calculate T(s): (-1/beta) * integral from s0 to s of dm/mI(m) 
function compute_T(s) {
  const I = m => params.s0 + params.i0 - m + (params.gamma * (1 - params.p_II)/(params.beta * params.p_SI)) * Math.log(m/params.s0);
  const f = m => {
    const Im = I(m);
    if (m <= 0 || !isFinite(Im) || Im === 0) return 0; // skip invalid points
    return (-1 / (params.beta * params.p_SI)) * 1 / (m * Im);
  };
  const a = params.s0;
  const b = s;
  const n = 50;
  // Check for valid integration range
  if (b <= 0 || b >= a) return null;
  return simpsonsRule(f, a, b, n);
}

function i_peak(){
  const s_p = params.gamma * (1 - params.p_II)/(params.beta * params.p_SI);
  const i_p =  params.s0 + params.i0 - s_p + (params.gamma * (1 - params.p_II)/(params.beta * params.p_SI)) * Math.log(s_p/params.s0);
  
  return i_p;
}

function compute_h_tpi(){
  const y = params.gamma * (1 - params.p_II) / (params.beta * params.p_SI); // R_0 ^ -1
  const tpi = compute_T(y);
  f = s => (1/s) * Math.exp(params.p_HR * params.alpha * compute_T(s) );
  
  const h_tpi = Math.exp(-params.p_HR * params.alpha * tpi) * ( params.h0 + params.p_IH * params.gamma * (-1/(params.beta * params.p_SI)) * simpsonsRule(f,params.s0,y,100) );
  
  return h_tpi;
}

// Function to calculate R₀ based on current parameter values
function calculateR0SIHRS(params) {
    return params.p_SI * params.beta / (params.gamma * (1 - params.p_II));
}

function calculateSIHRSThresholds(params) {
  const sigma = params.p_SI * params.beta * params.s0 / (params.gamma * (1 - params.p_II));
  const sigma1 = params.gamma * params.p_IH * params.i0 / (params.alpha * params.p_HR * params.h0);
 
  const i_p = i_peak();
  const h_tpi = compute_h_tpi();
  const sigma2 = params.gamma * params.p_IH * i_p / (params.alpha * params.p_HR * h_tpi);
  
  return [sigma, sigma1, sigma2];
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        validateParameters,
        sihrsStochasticModel,
        solveDeterministicSIHRS,
        simpsonsRule,
        compute_T,
        i_peak,
        compute_h_tpi,
        calculateR0SIHRS,
        calculateSIHRSThresholds
    };
} 