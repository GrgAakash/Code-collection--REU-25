// Main Application Entry Point
// Ties together all modules and initializes the application

// Global variables (shared across modules)
let params = {
    beta: 0.183,     // infection rate (adjusted for R₀ = 1.25)
    gamma: 0.126,   // I outflow rate (INCREASED for later peak)
    alpha: 0.1,    // H outflow rate
    lambda: 0.0083,   // reinfection rate (matching MATLAB)
    p_SI: 1.0,     // S->I probability
    p_II: 0.00,     // I->I (remains infected)
    p_IH: 0.04,     // I->H (matching MATLAB)
    p_IR: 0.959,     // I->R (matching MATLAB)
    p_ID: 0.001,     // I->D (death)
    p_HR: 0.9882,     // H->R
    p_HH: 0.01,     // H->H (remains hospitalized)
    p_HD: 0.0018,     // H->D (death)
    p_RR: 0.02,     // R->R (remains recovered)
    p_RS: 0.98,     // R->S (reinfection)
    tmax: 620,      // simulation end time (matching MATLAB)
    s0: 0.99,       // initial susceptible proportion
    i0: 0.01,       // initial infected proportion
    h0: 0.0,       // initial hospitalized proportion
    r0: 0.0,       // initial recovered proportion
    d0: 0.0,       // initial death proportion
    N: 3000,
    R_0_value: 0    // Add R₀ value to params
};

// Population size
let N = 3000;
let totalRuns = 50;

// Store results
const results = [];
let detResult = null;
let charts = {
    S: null,
    I: null,
    H: null,
    R: null,
    D: null
};

// Function to calculate dynamic y-axis range based on peak values (including real-world data)
function calculateDynamicYRange(compartment, detResult, results) {
    try {
        // Get all values for this compartment
        const allValues = [];
        
        // Add deterministic values
        if (detResult && detResult[`${compartment}_prop`] && Array.isArray(detResult[`${compartment}_prop`])) {
            allValues.push(...detResult[`${compartment}_prop`]);
        }
        
        // Add stochastic values from current run only (to avoid recursion)
        if (results && results.length > 0 && results[0] && results[0][`${compartment}_prop`] && Array.isArray(results[0][`${compartment}_prop`])) {
            allValues.push(...results[0][`${compartment}_prop`]);
        }
        
        // Add real-world data values if available
        if (typeof window.realDataHandler !== 'undefined' && 
            window.realDataHandler.showRealData && 
            window.realDataHandler.realData &&
            window.realDataHandler.realDataCompartment === compartment) {
            
            // Try to get real data from the chart if it exists and has been processed
            if (typeof charts !== 'undefined' && charts[compartment] && charts[compartment].data.datasets.length > 2) {
                const realDataDataset = charts[compartment].data.datasets[2]; // Real data is typically the 3rd dataset
                if (realDataDataset && realDataDataset.data && Array.isArray(realDataDataset.data)) {
                    const realDataValues = realDataDataset.data.map(point => point.y).filter(y => !isNaN(y) && y !== null);
                    if (realDataValues.length > 0) {
                        allValues.push(...realDataValues);
                        console.log(`📊 Dynamic Y-axis: Including real-world data peak (${Math.max(...realDataValues).toFixed(6)}) for compartment ${compartment}`);
                    }
                }
            }
        }
        
        if (allValues.length === 0) return { min: 0, max: 1 };
        
        const maxValue = Math.max(...allValues);
        const minValue = Math.min(...allValues);
        
        // Set y-axis range with padding
        const range = maxValue - minValue;
        
        // Special handling for infected (I) - use max(ODE peak, stochastic peak, real data peak) * 1.1
        if (compartment === 'I') {
            // Find the peak and min of deterministic values
            const deterministicValues = detResult && detResult[`${compartment}_prop`] ? detResult[`${compartment}_prop`] : [];
            const deterministicPeak = deterministicValues.length > 0 ? Math.max(...deterministicValues) : 0;
            const deterministicMin = deterministicValues.length > 0 ? Math.min(...deterministicValues) : 0;
            
            // Find the peak and min of stochastic values
            const stochasticValues = results && results.length > 0 && results[0] && results[0][`${compartment}_prop`] ? results[0][`${compartment}_prop`] : [];
            const stochasticPeak = stochasticValues.length > 0 ? Math.max(...stochasticValues) : 0;
            const stochasticMin = stochasticValues.length > 0 ? Math.min(...stochasticValues) : 0;
            
            // Find the peak and min of real-world data values
            let realDataPeak = 0;
            let realDataMin = Infinity;
            if (typeof window.realDataHandler !== 'undefined' && 
                window.realDataHandler.showRealData && 
                window.realDataHandler.realData &&
                window.realDataHandler.realDataCompartment === compartment &&
                typeof charts !== 'undefined' && charts[compartment] && charts[compartment].data.datasets.length > 2) {
                
                const realDataDataset = charts[compartment].data.datasets[2];
                if (realDataDataset && realDataDataset.data && Array.isArray(realDataDataset.data)) {
                    const realDataValues = realDataDataset.data.map(point => point.y).filter(y => !isNaN(y) && y !== null);
                    if (realDataValues.length > 0) {
                        realDataPeak = Math.max(...realDataValues);
                        realDataMin = Math.min(...realDataValues);
                    }
                }
            }
            
            // Use the highest peak among all three data sources and lowest min
            const maxPeak = Math.max(deterministicPeak, stochasticPeak, realDataPeak);
            const minValue = Math.min(deterministicMin, stochasticMin, realDataMin === Infinity ? 0 : realDataMin);
            const finalMax = maxPeak * 1.1;
            const finalMin = Math.max(0, minValue * 1.1);
            
            return {
                min: finalMin,
                max: finalMax
            };
        }
        
        // Special handling for hospitalization (H) - allow very small values, include real data peak
        if (compartment === 'H') {
            const padding = Math.max(range * 0.1, 0.000001); // Minimum padding of 0.000001 for H
            const calculatedMax = maxValue + padding;
            const finalMax = Math.min(calculatedMax, maxValue * 1.1); // Real data peak already included in maxValue
            
            return {
                min: Math.max(0, minValue - padding),
                max: finalMax
            };
        }
        
        // Default padding for other compartments (S, R, D) - include real data peak
        const padding = Math.max(range * 0.1, 0.01); // Minimum padding of 0.01
        
        // Allow Y-axis to go above 1.0 if needed (up to 1.1x peak) - real data peak already included in maxValue
        const calculatedMax = maxValue + padding;
        const finalMax = Math.min(calculatedMax, maxValue * 1.1); // Allow up to 1.1x the highest peak (including real data)
        
        return {
            min: Math.max(0, minValue - padding),
            max: finalMax
        };
    } catch (error) {
        console.warn('Error calculating dynamic range:', error);
        return { min: 0, max: 1 };
    }
}

// Check if Chart.js is loaded
function checkChartJS() {
    if (typeof Chart === 'undefined') {
        document.getElementById('errorMessage').style.display = 'block';
        document.getElementById('errorMessage').innerHTML = 
            'Error: Chart.js failed to load. Please check your internet connection and refresh the page.';
        return false;
    }
    return true;
}

// Initialize everything
async function initialize() {
    try {
        // Initialize audio context on first user interaction
        document.addEventListener('click', function initAudioOnClick() {
            if (!audioContext) {
                initAudio();
            }
            document.removeEventListener('click', initAudioOnClick);
        }, { once: true });

        // Check if Chart.js is loaded
        if (!checkChartJS()) return;
        
        validateParameters(params);
        
        // Run all simulations
        updateStatus('Running simulations...', 'running');
        for (let i = 0; i < totalRuns; i++) {
            results.push(sihrsStochasticModel(N, params, i + 1));
            if (i % 10 === 0) {
                updateStatus(`Running simulations... ${i + 1}/${totalRuns}`, 'running');
                await new Promise(resolve => setTimeout(resolve, 10)); // Allow UI to update
            }
        }
        
        detResult = solveDeterministicSIHRS(params);
        
        // Set up Charts.js
        updateStatus('Setting up visualization...', 'running');
        
        // Destroy existing charts if they exist
        Object.values(charts).forEach(chart => {
            if (chart) {
                chart.destroy();
            }
        });

        // Create separate chart for each compartment
        const compartments = [
            { id: 'S', color: '#FF6B6B', name: 'Susceptible' },
            { id: 'I', color: '#4ECDC4', name: 'Infected' },
            { id: 'H', color: '#45B7D1', name: 'Hospitalized' },
            { id: 'R', color: '#96CEB4', name: 'Recovered' },
            { id: 'D', color: '#DDA0DD', name: 'Death' }
        ];

        compartments.forEach(comp => {
            const ctx = document.getElementById(`chart${comp.id}`).getContext('2d');
            
            // Calculate dynamic y-axis range for this compartment
            const yRange = calculateDynamicYRange(comp.id, detResult, results);
            
            charts[comp.id] = new Chart(ctx, {
                type: 'line',
                data: {
                    datasets: [
                        {
                            label: 'Deterministic',
                            data: detResult.T.map((t, i) => ({ x: t, y: detResult[`${comp.id}_prop`][i] })),
                            borderColor: comp.color,
                            backgroundColor: `${comp.color}20`,
                            borderDash: [8, 4],
                            borderWidth: 3,
                            pointRadius: 0,
                            fill: false,
                            tension: 0.1
                        },
                        {
                            label: 'Stochastic',
                            data: results[0].T.map((t, i) => ({ x: t, y: results[0][`${comp.id}_prop`][i] })),
                            borderColor: comp.color,
                            backgroundColor: `${comp.color}20`,
                            borderWidth: 2,
                            pointRadius: 0,
                            fill: false,
                            tension: 0.2
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    animation: {
                        duration: 0
                    },
                    interaction: {
                        intersect: false,
                        mode: 'index'
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: `${comp.name} Proportion Over Time - Run 1/${totalRuns} (N=${N})`,
                            font: {
                                family: "'Press Start 2P', cursive",
                                size: 16,
                                weight: 'bold'
                            },
                            color: '#000'
                        },
                        legend: {
                            position: 'top',
                            labels: {
                                color: '#000',
                                font: {
                                    family: "'Press Start 2P', cursive",
                                    size: 10
                                }
                            }
                        },
                        tooltip: {
                            backgroundColor: 'rgba(0, 0, 0, 0.8)',
                            titleColor: '#fff',
                            bodyColor: '#fff',
                            borderColor: '#fff',
                            borderWidth: 1,
                            padding: 8,
                            boxPadding: 4,
                            usePointStyle: true,
                            callbacks: {
                                label: function(context) {
                                    return `${context.dataset.label}: ${context.parsed.y.toFixed(4)}`;
                                }
                            }
                        }
                    },
                    scales: {
                        x: {
                            type: 'linear',
                            title: {
                                display: true,
                                text: 'Time (days)',
                                font: {
                                    family: "'Press Start 2P', cursive",
                                    size: 10,
                                    weight: 'bold'
                                },
                                color: '#000'
                            },
                            min: 0, 
                            max: params.tmax,
                            ticks: { 
                                color: '#000',
                                font: {
                                    family: 'Inter',
                                    size: 10
                                }
                            },
                            grid: { 
                                color: 'rgba(0, 0, 0, 0.05)'
                            }
                        },
                        y: { 
                            title: {
                                display: true,
                                text: 'Proportion',
                                font: {
                                    family: "'Press Start 2P', cursive",
                                    size: 10,
                                    weight: 'bold'
                                },
                                color: '#000'
                            },
                            min: yRange.min,
                            max: yRange.max,
                            ticks: { 
                                color: '#000',
                                font: {
                                    family: 'Inter',
                                    size: 10
                                },
                                maxTicksLimit: 6, // Limit number of ticks to prevent duplicates
                                callback: function(value) {
                                    // Show more decimal places for small values
                                    if (value < 0.01) {
                                        return value.toFixed(4); // Show 4 decimal places for small values
                                    } else if (value < 0.1) {
                                        return value.toFixed(3); // Show 3 decimal places for medium values
                                    } else {
                                        return value.toFixed(2); // Show 2 decimal places for larger values
                                    }
                                }
                            },
                            grid: { 
                                color: 'rgba(0, 0, 0, 0.05)'
                            }
                        }
                    }
                }
            });
        });

        // Setup compartment buttons
        document.querySelectorAll('.comp-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                const selectedComp = this.getAttribute('data-comp');
                showCompartment(selectedComp);
            });
        });

        // Initialize real data handlers
        if (typeof initRealDataHandlers === 'function') {
            initRealDataHandlers();
        }
        
        // Start animation
        updateStatus(`Animation ready - Run 1/${totalRuns}`, 'running');
        setTimeout(() => animate(), 1000); // Start after 1 second
        
        console.log('✅ Initialization complete');
        
    } catch (error) {
        console.error('❌ Error:', error.message);
        updateStatus(`Error: ${error.message}`, 'stopped');
        document.getElementById('errorMessage').style.display = 'block';
        document.getElementById('errorMessage').innerHTML = `Error: ${error.message}`;
    }
}

// Function to show selected compartment chart
function showCompartment(compartment) {
    // Update button states
    document.querySelectorAll('.comp-btn').forEach(btn => {
        btn.classList.remove('active');
    });
    document.querySelector(`[data-comp="${compartment}"]`).classList.add('active');
    
    // Hide all chart containers
    document.querySelectorAll('.chart-container').forEach(container => {
        container.style.display = 'none';
    });
    
    // Show selected chart container
    document.getElementById(`chartContainer${compartment}`).style.display = 'block';
    
    // 🆕 Smart real-world data display based on current compartment view
    if (typeof window.realDataHandler !== 'undefined') {
        window.realDataHandler.updateRealDataForCompartment(compartment);
    }
    
    // Play sound effect
    playCoinSound();
}



// Set up event listeners for all parameters
['s0','i0','r0','h0','d0','beta', 'gamma', 'alpha', 'lambda', 'p_SI'].forEach(param => {
    const rangeInput = document.getElementById(param);
    const numberInput = document.getElementById(`${param}-number`);
    if (!rangeInput || !numberInput) return;
    rangeInput.addEventListener('input', (e) => {
        playCoinSound();
        updateParameter(param, e.target.value);
    });
    numberInput.addEventListener('input', (e) => {
        let value = e.target.value;
        if (value === '' || value === '.' || value === '-.') return;
        const numValue = parseFloat(value);
        if (!isNaN(numValue)) {
            const min = parseFloat(rangeInput.min);
            const max = parseFloat(rangeInput.max);
            if (numValue >= min && numValue <= max) {
                playCoinSound();
                updateParameter(param, numValue);
            }
        }
    });
    numberInput.addEventListener('blur', (e) => {
        let value = e.target.value;
        if (value === '' || value === '.' || value === '-.') {
            value = rangeInput.value;
            e.target.value = value;
        }
        const numValue = parseFloat(value);
        if (!isNaN(numValue)) {
            const min = parseFloat(rangeInput.min);
            const max = parseFloat(rangeInput.max);
            if (numValue >= min && numValue <= max) {
                updateParameter(param, numValue);
            }
        }
    });
});

// Auto-update logic for I outflow group (p_II, p_IH, p_IR, p_ID)
function updateIOutflow(changed) {
    let p_II = parseFloat(document.getElementById('p_II-number').value);
    let p_IH = parseFloat(document.getElementById('p_IH-number').value);
    let p_IR = parseFloat(document.getElementById('p_IR-number').value);
    let p_ID = parseFloat(document.getElementById('p_ID-number').value);
    if (changed === 'p_II') {
        p_IH = Math.max(0, Math.min(1, 1 - p_II - p_IR - p_ID));
        document.getElementById('p_IH-number').value = p_IH.toFixed(9);
    } else if (changed === 'p_IH') {
        p_IR = Math.max(0, Math.min(1, 1 - p_II - p_IH - p_ID));
        document.getElementById('p_IR-number').value = p_IR.toFixed(9);
    } else if (changed === 'p_IR') {
        p_ID = Math.max(0, Math.min(1, 1 - p_II - p_IH - p_IR));
        document.getElementById('p_ID-number').value = p_ID.toFixed(9);
    } else if (changed === 'p_ID') {
        p_II = Math.max(0, Math.min(1, 1 - p_IH - p_IR - p_ID));
        document.getElementById('p_II-number').value = p_II.toFixed(9);
    }
    params.p_II = parseFloat(document.getElementById('p_II-number').value);
    params.p_IH = parseFloat(document.getElementById('p_IH-number').value);
    params.p_IR = parseFloat(document.getElementById('p_IR-number').value);
    params.p_ID = parseFloat(document.getElementById('p_ID-number').value);
}
['p_II', 'p_IH', 'p_IR', 'p_ID'].forEach(param => {
    document.getElementById(`${param}-number`).addEventListener('input', (e) => {
        updateIOutflow(param);
    });
});

// Auto-update logic for H outflow group (p_HR, p_HH, p_HD)
function updateHOutflow(changed) {
    let p_HR = parseFloat(document.getElementById('p_HR-number').value);
    let p_HH = parseFloat(document.getElementById('p_HH-number').value);
    let p_HD = parseFloat(document.getElementById('p_HD-number').value);
    if (changed === 'p_HR') {
        p_HH = Math.max(0, Math.min(1, 1 - p_HR - p_HD));
        document.getElementById('p_HH-number').value = p_HH.toFixed(9);
    } else if (changed === 'p_HH') {
        p_HD = Math.max(0, Math.min(1, 1 - p_HR - p_HH));
        document.getElementById('p_HD-number').value = p_HD.toFixed(9);
    } else if (changed === 'p_HD') {
        p_HR = Math.max(0, Math.min(1, 1 - p_HH - p_HD));
        document.getElementById('p_HR-number').value = p_HR.toFixed(9);
    }
    params.p_HR = parseFloat(document.getElementById('p_HR-number').value);
    params.p_HH = parseFloat(document.getElementById('p_HH-number').value);
    params.p_HD = parseFloat(document.getElementById('p_HD-number').value);
}
['p_HR', 'p_HH', 'p_HD'].forEach(param => {
    document.getElementById(`${param}-number`).addEventListener('input', (e) => {
        updateHOutflow(param);
    });
});

// Auto-update logic for R outflow group (p_RR, p_RS)
function updateROutflow(changed) {
    let p_RR = parseFloat(document.getElementById('p_RR-number').value);
    let p_RS = parseFloat(document.getElementById('p_RS-number').value);
    if (changed === 'p_RR') {
        p_RS = Math.max(0, Math.min(1, 1 - p_RR));
        document.getElementById('p_RS-number').value = p_RS.toFixed(9);
    } else if (changed === 'p_RS') {
        p_RR = Math.max(0, Math.min(1, 1 - p_RS));
        document.getElementById('p_RR-number').value = p_RR.toFixed(9);
    }
    params.p_RR = parseFloat(document.getElementById('p_RR-number').value);
    params.p_RS = parseFloat(document.getElementById('p_RS-number').value);
}
['p_RR', 'p_RS'].forEach(param => {
    document.getElementById(`${param}-number`).addEventListener('input', (e) => {
        updateROutflow(param);
    });
});

// Population size input handling
const populationInput = document.getElementById('population');
const populationNumberInput = document.getElementById('population-number');
const populationValue = document.getElementById('population-value');

function updatePopulation(value) {
    value = Math.max(1, Math.min(100000, value)); // Increased max to 100,000
    populationInput.value = value;
    populationNumberInput.value = value;
    populationValue.textContent = value;
    params.N = value;
}

populationInput.addEventListener('input', (e) => {
    playCoinSound(); // Play coin sound for slider movement
    updatePopulation(parseInt(e.target.value));
});

populationNumberInput.addEventListener('input', (e) => {
    let value = parseInt(e.target.value);
    if (!isNaN(value)) {
        playCoinSound(); // Play coin sound for number input
        updatePopulation(value);
    }
});

// Number of runs input handling
const runsInput = document.getElementById('runs');
const runsNumberInput = document.getElementById('runs-number');
const runsValue = document.getElementById('runs-value');

function updateRuns(value) {
    value = Math.max(1, Math.min(1000, value));
    runsInput.value = value;
    runsNumberInput.value = value;
    runsValue.textContent = value;
    totalRuns = value;
}

runsInput.addEventListener('input', (e) => {
    playCoinSound(); // Play coin sound for slider movement
    updateRuns(parseInt(e.target.value));
});

runsNumberInput.addEventListener('input', (e) => {
    let value = parseInt(e.target.value);
    if (!isNaN(value)) {
        playCoinSound(); // Play coin sound for number input
        updateRuns(value);
    }
});

// Setup apply button
document.getElementById('applyParams').addEventListener('click', () => {
    playPowerUpSound();
    if (isRunning) {
        toggleAnimation();
    }
    params.beta = parseFloat(document.getElementById('beta').value);
    params.gamma = parseFloat(document.getElementById('gamma').value);
    params.alpha = parseFloat(document.getElementById('alpha').value);
    params.lambda = parseFloat(document.getElementById('lambda').value);
    params.p_SI = parseFloat(document.getElementById('p_SI').value);
    params.p_II = parseFloat(document.getElementById('p_II-number').value);
    params.p_IH = parseFloat(document.getElementById('p_IH-number').value);
    params.p_IR = parseFloat(document.getElementById('p_IR-number').value);
    params.p_ID = parseFloat(document.getElementById('p_ID-number').value);
    params.p_HR = parseFloat(document.getElementById('p_HR-number').value);
    params.p_HH = parseFloat(document.getElementById('p_HH-number').value);
    params.p_HD = parseFloat(document.getElementById('p_HD-number').value);
    params.p_RR = parseFloat(document.getElementById('p_RR-number').value);
    params.p_RS = parseFloat(document.getElementById('p_RS-number').value);
    params.s0 = parseFloat(document.getElementById('s0').value);
    params.i0 = parseFloat(document.getElementById('i0').value);
    params.r0 = parseFloat(document.getElementById('r0').value);
    params.h0 = parseFloat(document.getElementById('h0').value);
    params.d0 = parseFloat(document.getElementById('d0').value);
    N = parseInt(document.getElementById('population').value);
    totalRuns = parseInt(document.getElementById('runs').value);

    // Validate parameters
    try {
        validateParameters(params);
        
        // Calculate R₀ once when parameters are applied
        params.R_0_value = calculateR0SIHRS(params);
        
        // Clear existing results
        results.length = 0;
        currentRun = 0;

        // Run new simulations
        updateStatus('Running new simulations...', 'running');
        for (let i = 0; i < totalRuns; i++) {
            results.push(sihrsStochasticModel(N, params, i + 1));
            if (i % 10 === 0) {
                updateStatus(`Running simulations... ${i + 1}/${totalRuns}`, 'running');
            }
        }

        // Update deterministic solution
        detResult = solveDeterministicSIHRS(params);

                // Update all charts
        const compartments = ['S', 'I', 'H', 'R', 'D'];
        compartments.forEach(comp => {
            if (charts[comp]) {
                charts[comp].data.datasets[0].data = detResult.T.map((t, i) => ({ x: t, y: detResult[`${comp}_prop`][i] }));
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                charts[comp].options.scales.x.max = params.tmax;
                charts[comp].update('none');
            }
        });
        updateStatistics();

        updateStatus('Parameters applied successfully', 'running');
    } catch (error) {
        updateStatus(`Error: ${error.message}`, 'stopped');
    }
});

// Add click handler for analyze button
document.getElementById('analyzeBtn').addEventListener('click', function() {
    playPowerUpSound(); // Play power-up sound for analyze
    const patternSection = document.getElementById('patternAnalysis');
    if (patternSection.style.display === 'none') {
        // Get the currently active compartment
        const activeButton = document.querySelector('.comp-btn.active');
        const currentCompartment = activeButton ? activeButton.getAttribute('data-comp') : 'H';
        updatePatternAnalysis(results, currentCompartment);
    } else {
        patternSection.style.display = 'none';
    }
});

// Add click handlers for control buttons
document.getElementById('playBtn').addEventListener('click', toggleAnimation);
document.getElementById('resetBtn').addEventListener('click', resetAnimation);
document.getElementById('speedBtn').addEventListener('click', changeSpeed);

// Function to show pattern details in modal for a specific compartment
function showPatternDetails(patternIndex, runIndices, compartment = 'H') {
    const modal = document.getElementById('patternModal');
    const modalTitle = document.getElementById('modalTitle');
    const modalContent = document.getElementById('modalContent');
    
    modalTitle.textContent = `Pattern ${patternIndex + 1} - ${runIndices.length} Runs`;
    
    // Create content for modal
    let content = `
        <div class="pattern-runs-grid">
    `;
    
    runIndices.forEach((runIndex, i) => {
        const run = results[runIndex];
        const maxValue = Math.max(...run[`${compartment}_prop`]);
        const peakTime = run.T[run[`${compartment}_prop`].indexOf(maxValue)];
        
        content += `
            <div class="run-card">
                <div class="run-header">
                    <div class="run-selector">
                        <input type="radio" id="run${runIndex}" name="selectedRun" value="${runIndex}" class="run-radio" ${i === 0 ? 'checked' : ''}>
                        <label for="run${runIndex}" class="run-number">Run ${runIndex + 1}</label>
                    </div>
                    <span class="run-stats">Peak: ${maxValue.toFixed(4)} at ${peakTime.toFixed(1)}d</span>
                </div>
                <div class="run-chart">
                    <canvas id="modalRun${runIndex}Chart"></canvas>
                </div>
            </div>
        `;
    });
    
    content += `
        </div>
        <div class="download-section">
            <h3>Download Options</h3>
            <button class="download-btn mario-btn" onclick="downloadPatternAsPNG(${patternIndex}, ${JSON.stringify(runIndices)}, results, params, '${compartment}')">
                Download All Runs as PNG
            </button>
            <button class="download-btn mario-btn" onclick="downloadSelectedIndividualRun(${JSON.stringify(runIndices)}, results, params, '${compartment}')">
                Download Selected Individual Run
            </button>
        </div>
    `;
    
    modalContent.innerHTML = content;
    modal.style.display = 'block';
    
    // Create charts for each run with both deterministic and stochastic
    runIndices.forEach((runIndex, i) => {
        setTimeout(() => {
            const run = results[runIndex];
            const ctx = document.getElementById(`modalRun${runIndex}Chart`).getContext('2d');
            
            // Get deterministic solution for comparison (SIHRS version)
            const detResult = solveDeterministicSIHRS(params);
            
            // Create datasets for the specific compartment only
            const compartmentColors = {
                'S': '#FF6B6B',
                'I': '#4ECDC4', 
                'H': '#45B7D1',
                'R': '#96CEB4',
                'D': '#DDA0DD'
            };
            
            const datasets = [
                {
                    label: `Deterministic ${compartment.toLowerCase()}`,
                    data: detResult.T.map((t, j) => ({ x: t, y: detResult[`${compartment}_prop`][j] })),
                    borderColor: compartmentColors[compartment] || 'var(--deterministic_i-color)',
                    backgroundColor: `${compartmentColors[compartment] || 'var(--deterministic_i-color)'}20`,
                    borderDash: [8, 4],
                    borderWidth: 2,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.1
                },
                {
                    label: `Stochastic ${compartment.toLowerCase()}`,
                    data: run.T.map((t, j) => ({ x: t, y: run[`${compartment}_prop`][j] })),
                    borderColor: compartmentColors[compartment] || 'var(--stochastic-color)',
                    backgroundColor: `${compartmentColors[compartment] || 'var(--stochastic-color)'}20`,
                    borderWidth: 1,
                    pointRadius: 0,
                    fill: false,
                    tension: 0.2
                }
            ];
            
            new Chart(ctx, {
                type: 'line',
                data: {
                    datasets: datasets
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    layout: {
                        padding: {
                            bottom: 25,
                            left: 15,
                            right: 15
                        }
                    },
                    plugins: {
                        title: {
                            display: false
                        },
                        legend: { 
                            display: true,
                            position: 'top',
                            labels: {
                                usePointStyle: true,
                                padding: 12,
                                font: { size: 11 }
                            }
                        }
                    },
                    scales: {
                        x: {
                            type: 'linear',
                            title: { 
                                display: true, 
                                text: 'Time (days)',
                                font: { size: 13, weight: 'bold' }
                            },
                            ticks: { 
                                maxTicksLimit: 8,
                                font: { size: 11 }
                            },
                            grid: { 
                                color: 'rgba(0, 0, 0, 0.05)'
                            }
                        },
                        y: {
                            title: { 
                                display: true, 
                                text: 'Proportion',
                                font: { size: 13, weight: 'bold' }
                            },
                            beginAtZero: true,
                            ticks: { font: { size: 11 } },
                            grid: { 
                                color: 'rgba(0, 0, 0, 0.05)'
                            }
                        }
                    }
                }
            });
        }, i * 100); // Stagger chart creation
    });
}

// Function to close modal
function closePatternModal() {
    document.getElementById('patternModal').style.display = 'none';
}

// Close modal when clicking outside
window.onclick = function(event) {
    const modal = document.getElementById('patternModal');
    if (event.target === modal) {
        closePatternModal();
    }
}

// Modal logic for ODE info
const odeModal = document.getElementById('odeModal');
const closeOdeModal = document.getElementById('closeOdeModal');

// Attach click event only to the 1-star Dragon Ball
const dragonBall1 = document.querySelector('.dragon-ball-1');
if (odeModal && closeOdeModal && dragonBall1) {
    dragonBall1.addEventListener('click', () => {
        odeModal.style.display = 'flex';
        if (window.MathJax) {
            MathJax.typesetPromise();
        }
    });
    closeOdeModal.addEventListener('click', () => {
        odeModal.style.display = 'none';
    });
    // Close modal when clicking outside modal-content
    odeModal.addEventListener('click', (e) => {
        if (e.target === odeModal) {
            odeModal.style.display = 'none';
        }
    });
}

// Modal logic for Hospitalization Pattern Classification (2-star ball)
const patternClassModal = document.getElementById('patternClassModal');
const closePatternClassModal = document.getElementById('closePatternClassModal');
const dragonBall2 = document.querySelector('.dragon-ball-2');
if (patternClassModal && closePatternClassModal && dragonBall2) {
    dragonBall2.addEventListener('click', () => {
        patternClassModal.style.display = 'flex';
        if (window.MathJax) {
            MathJax.typesetPromise();
        }
    });
    closePatternClassModal.addEventListener('click', () => {
        patternClassModal.style.display = 'none';
    });
    patternClassModal.addEventListener('click', (e) => {
        if (e.target === patternClassModal) {
            patternClassModal.style.display = 'none';
        }
    });
}



// --- Dynamic initial conditions logic ---
function updateInitialConditions(changed) {
    let s0 = parseFloat(document.getElementById('s0-number').value);
    let i0 = parseFloat(document.getElementById('i0-number').value);
    let h0 = parseFloat(document.getElementById('h0-number').value);
    let r0 = parseFloat(document.getElementById('r0-number').value);
    let d0 = parseFloat(document.getElementById('d0-number').value);
    let total = s0 + i0 + h0 + r0 + d0;
    
    // If sum > 1, cap the changed value
    if (total > 1) {
        if (changed === 's0') {
            s0 = Math.max(0, 1 - i0 - h0 - r0 - d0);
            document.getElementById('s0').value = s0.toFixed(9);
            document.getElementById('s0-number').value = s0.toFixed(9);
        } else if (changed === 'i0') {
            i0 = Math.max(0, 1 - s0 - h0 - r0 - d0);
            document.getElementById('i0').value = i0.toFixed(9);
            document.getElementById('i0-number').value = i0.toFixed(9);
        } else if (changed === 'h0') {
            h0 = Math.max(0, 1 - s0 - i0 - r0 - d0);
            document.getElementById('h0').value = h0.toFixed(9);
            document.getElementById('h0-number').value = h0.toFixed(9);
        } else if (changed === 'r0') {
            r0 = Math.max(0, 1 - s0 - i0 - h0 - d0);
            document.getElementById('r0-number').value = r0.toFixed(9);
        } else if (changed === 'd0') {
            d0 = Math.max(0, 1 - s0 - i0 - h0 - r0);
            document.getElementById('d0-number').value = d0.toFixed(9);
        }
    }
    
    document.getElementById('s0-value').textContent = s0.toFixed(6);
    document.getElementById('i0-value').textContent = i0.toFixed(6);
    document.getElementById('h0-value').textContent = h0.toFixed(6);
    document.getElementById('r0-value').textContent = r0.toFixed(6);
    document.getElementById('d0-value').textContent = d0.toFixed(6);
    // Update params
    params.s0 = s0;
    params.i0 = i0;
    params.h0 = h0;
    params.r0 = r0;
    params.d0 = d0;
    // Update thresholds and chart if needed
    params.R_0_value = calculateR0SIHRS(params);
    const [sigma0, sigma1, sigma2] = calculateSIHRSThresholds(params);
    const tpi = compute_T(params.gamma * (1 - params.p_II) / (params.beta * params.p_SI));
    const h_tpi = compute_h_tpi();
    document.getElementById('R_0-value').textContent = params.R_0_value.toFixed(2);
    document.getElementById('sigma0-value').textContent = sigma0.toFixed(2);
    document.getElementById('sigma1-value').textContent = sigma1.toFixed(2);
    document.getElementById('sigma2-value').textContent = sigma2.toFixed(2);
    document.getElementById('tpi-value').textContent = isNaN(tpi) ? '--' : tpi.toFixed(2);
    document.getElementById('h_tpi-value').textContent = isNaN(h_tpi) ? '--' : h_tpi.toFixed(2);
}
['s0', 'i0', 'h0', 'r0', 'd0'].forEach(param => {
    document.getElementById(param).addEventListener('input', (e) => {
        updateInitialConditions(param);
    });
    document.getElementById(param + '-number').addEventListener('input', (e) => {
        updateInitialConditions(param);
    });
});

// Start when page loads
window.addEventListener('load', async () => {
    // Hide pattern details modal on load
    const patternModal = document.getElementById('patternModal');
    if (patternModal) patternModal.style.display = 'none';

    console.log('🔄 Page loaded, initializing...');
    initMatrixBackground(); // Initialize Matrix background
    await initialize();
    initDraggableCoins(); // Initialize draggable coins
    if (window.initDraggableDragonBalls) {
        window.initDraggableDragonBalls();
    } else if (typeof initDraggableDragonBalls === 'function') {
        initDraggableDragonBalls();
    }
});

function safeUpdateThresholds() {
    const sum = params.s0 + params.i0 + params.h0 + params.r0 + params.d0;
    if (Math.abs(sum - 1) > 1e-8 || params.s0 <= 0) {
        document.getElementById('tpi-value').textContent = '--';
        document.getElementById('sigma0-value').textContent = '--';
        document.getElementById('sigma1-value').textContent = '--';
        document.getElementById('sigma2-value').textContent = '--';
        document.getElementById('h_tpi-value').textContent = '--';
        document.getElementById('R_0-value').textContent = '--';
        return;
    }
    const s_p = params.gamma * (1 - params.p_II) / (params.beta * params.p_SI);
    if (s_p <= 0 || s_p >= params.s0) {
        document.getElementById('tpi-value').textContent = '--';
        document.getElementById('h_tpi-value').textContent = '--';
        return;
    }
    const [sigma0, sigma1, sigma2] = calculateSIHRSThresholds(params);
    const tpi = compute_T(s_p);
    const h_tpi = compute_h_tpi();
    params.R_0_value = calculateR0SIHRS(params);
    document.getElementById('R_0-value').textContent = isNaN(params.R_0_value) ? '--' : params.R_0_value.toFixed(2);
    document.getElementById('sigma0-value').textContent = isNaN(sigma0) ? '--' : sigma0.toFixed(2);
    document.getElementById('sigma1-value').textContent = isNaN(sigma1) ? '--' : sigma1.toFixed(2);
    document.getElementById('sigma2-value').textContent = isNaN(sigma2) ? '--' : sigma2.toFixed(2);
    document.getElementById('tpi-value').textContent = isNaN(tpi) ? '--' : tpi.toFixed(2);
    document.getElementById('h_tpi-value').textContent = isNaN(h_tpi) ? '--' : h_tpi.toFixed(2);
}

// Call safeUpdateThresholds after any parameter change
['s0', 'i0', 'h0', 'r0', 'd0', 'beta', 'gamma', 'alpha', 'lambda', 'p_SI', 'p_II', 'p_IH', 'p_IR', 'p_ID', 'p_HR', 'p_HH', 'p_HD', 'p_RR', 'p_RS'].forEach(param => {
    const rangeInput = document.getElementById(param);
    const numberInput = document.getElementById(`${param}-number`);
    if (rangeInput) rangeInput.addEventListener('input', safeUpdateThresholds);
    if (numberInput) numberInput.addEventListener('input', safeUpdateThresholds);
}); 