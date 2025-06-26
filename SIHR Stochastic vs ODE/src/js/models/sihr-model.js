// SIHR Model Implementation
// Stochastic and Deterministic epidemic models

// Validate parameters
function validateParameters(params) {
    if ([params.beta, params.gamma, params.alpha].some(x => x <= 0)) {
        throw new Error('All rates must be positive');
    }
    if ([params.p1, params.p2, params.p3, params.ph].some(x => x <= 0 || x > 1)) {
        throw new Error('All probabilities must be in (0,1]');
    }
    if (Math.abs(params.s0 + params.i0 + params.h0 + params.r0 - 1) > 1e-10) {
        throw new Error('Initial conditions must sum to 1');
    }
}

// Stochastic SIHR model using Gillespie algorithm
function sirAgentModel(N, params, runNumber) {
    console.log(`🎲 Running stochastic simulation ${runNumber}/${totalRuns}`);
    
    // Calculate R₀ for this run
    const R_0 = params.p1 * params.beta / (params.p2 * params.gamma);
    
    // Initial conditions
    let s = Math.round(params.s0 * N);
    let i = Math.round(params.i0 * N);
    let h = Math.round(params.h0 * N);
    let r = Math.round(params.r0 * N);
    const total = s + i + h + r;
    if (total !== N) {
        s += N - total;
    }

    const T = [0], I_prop = [i / N], H_prop = [h / N];
    let t = 0;

    // Main simulation loop
    while ((i > 0 || h > 0) && t < params.tmax)  {
        const nS = s, nI = i, nH = h;

        // Calculate event rates
        const infectionRate = params.p1 * params.beta * nS * nI / N;
        const toHospitalRate = params.p2 * params.ph * params.gamma * nI;
        const toRecoveredFromIRate = params.p2 * (1 - params.ph) * params.gamma * nI;
        const toRecoveredFromHRate = params.p3 * params.alpha * nH;
        const totalRate = infectionRate + toHospitalRate + toRecoveredFromIRate + toRecoveredFromHRate;

        if (totalRate === 0) break;

        // Time to next event
        const dt = -Math.log(Math.random()) / totalRate;
        t += dt;

        if (t > params.tmax) {
            t = params.tmax;
            T.push(t);
            I_prop.push(i / (s + i + h + r))
            H_prop.push(h / (s + i + h + r));
            break;
        }

        // Determine which event occurs
        const chance = Math.random() * totalRate;
        if (chance < infectionRate && nS > 0) {
            s--; i++;
        } else if (chance < infectionRate + toHospitalRate && nI > 0) {
            i--; h++;
        } else if (chance < infectionRate + toHospitalRate + toRecoveredFromIRate && nI > 0) {
            i--; r++;
        } else if (nH > 0) {
            h--; r++;
        }

        T.push(t);
        I_prop.push(i / (s + i + h + r));
        H_prop.push(h / (s + i + h + r));
    }

    return { T, I_prop, H_prop, maxI: Math.max(...I_prop), maxH: Math.max(...H_prop) };
}

// Deterministic SIHR model (RK4 integration)
function solveDeterministicSIR(params) {
    console.log('📊 Solving deterministic model');
    const dt = 0.1;
    const T = [];
    const Y = [[params.s0, params.i0, params.h0, params.r0]];
    
    for (let t = 0; t <= params.tmax; t += dt) {
        T.push(t);
        if (t > 0) {
            const y = Y[Y.length - 1];
            const s = y[0], i = y[1], h = y[2];
            
            // ODE system
            const dsdt = -params.p1 * params.beta * s * i;
            const didt = params.p1 * params.beta * s * i - params.p2 * params.gamma * i;
            const dhdt = params.p2 * params.ph * params.gamma * i - params.p3 * params.alpha * h;
            const drdt = params.p2 * (1 - params.ph) * params.gamma * i + params.p3 * params.alpha * h;
            
            // Euler method (simplified for performance)
            const yNext = [
                s + dt * dsdt,
                i + dt * didt,
                h + dt * dhdt,
                y[3] + dt * drdt
            ];
            Y.push(yNext);
        }
    }

    return { T, I_prop: Y.map(y => y[1]), H_prop: Y.map(y => y[2]) };
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
  const I = m => params.s0 + params.i0 - m + (params.p2 * params.gamma/(params.p1 * params.beta)) * Math.log(m/params.s0);
  const f = m => (-1 / (params.p1 * params.beta)) * 1 / (m * I(m) );
  const a = params.s0;
  const b = s;
  const n = 50;

  return simpsonsRule(f,a,b,n);
}

function i_peak(){
  const s_p = params.p2 * params.gamma/(params.p1 * params.beta);
  const i_p =  params.s0 + params.i0 - s_p + (params.p2 * params.gamma/(params.p1 * params.beta)) * Math.log(s_p/params.s0);
  
  return i_p;
}

function compute_h_tpi(){
  const y = params.p2 * params.gamma / (params.p1 * params.beta); // R_0 ^ -1
  const tpi = compute_T(y);
  f = s => (1/s) * Math.exp(params.p3 * params.alpha * compute_T(s) );
  
  const h_tpi = Math.exp(-params.p3 * params.alpha * tpi) * ( params.h0 + params.ph * params.p1 * params.gamma * (-1/(params.p1 * params.beta)) * simpsonsRule(f,params.s0,y,100) );
  
  return h_tpi;
}

// Function to calculate R₀ based on current parameter values
function calculateR0(params) {
    return params.p1 * params.beta / (params.p2 * params.gamma);
}

function calculate_thresholds(params) {
  const sigma = params.p1 * params.beta * params.s0/ (params.p2 * params.gamma);
  const sigma1 = params.gamma * params.ph * params.p2 * params.i0 / (params.alpha * params.p3 * params.h0);
 
  const i_p = i_peak();
  const h_tpi = compute_h_tpi();
  const sigma2 =  params.gamma * params.ph * params.p2 * i_p / (params.alpha * params.p3 * h_tpi);
  
  return [sigma, sigma1, sigma2];
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        validateParameters,
        sirAgentModel,
        solveDeterministicSIR,
        simpsonsRule,
        compute_T,
        i_peak,
        compute_h_tpi,
        calculateR0,
        calculate_thresholds
    };
} 