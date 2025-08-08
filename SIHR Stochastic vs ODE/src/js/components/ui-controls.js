// UI Controls Component
// Functions for handling user interface interactions and animations

// Animation controls
let isRunning = true;
let currentRun = 0;
let animationSpeed = 800; // milliseconds per frame
let animationTimeout;
const speeds = [400, 800, 1200]; // Fast, Normal, Slow
let speedIndex = 1;

// Update status
function updateStatus(message, className = 'running') {
    const statusEl = document.getElementById('status');
    statusEl.textContent = message;
    statusEl.className = className;
}

// Animation controls
function toggleAnimation() {
    playJumpSound(); // Play jump sound for play/pause
    isRunning = !isRunning;
    const btn = document.getElementById('playBtn');
    if (isRunning) {
        btn.innerHTML = '<span class="btn-text">Pause</span>';
        updateStatus(`Animation running - Run ${currentRun + 1}/${totalRuns}`, 'running');
        animate();
    } else {
        btn.innerHTML = '<span class="btn-text">Play</span>';
        updateStatus('Animation paused', 'paused');
        if (animationTimeout) clearTimeout(animationTimeout);
    }
}

function resetAnimation() {
    playGameOverSound(); // Play game over sound for reset
    // Stop any ongoing animation
    if (animationTimeout) {
        clearTimeout(animationTimeout);
    }
    
    // Reset current run
    currentRun = 0;

    // Hide pattern analysis section
    const patternSection = document.getElementById('patternAnalysis');
    if (patternSection) {
        patternSection.style.display = 'none';
    }

    // Clear existing results
    results.length = 0;

    // Rerun stochastic simulations with current parameters
    updateStatus('Running new simulations...', 'running');
    for (let i = 0; i < totalRuns; i++) {
        if (typeof sihrsStochasticModel !== 'undefined') {
            // SIHRS model
            results.push(sihrsStochasticModel(N, params, i + 1));
        } else if (typeof sirAgentModel !== 'undefined') {
            // Original SIHR model
            results.push(sirAgentModel(N, params, i + 1));
        }
    }

    // Update deterministic solution
    if (typeof solveDeterministicSIHRS !== 'undefined') {
        // SIHRS model
        detResult = solveDeterministicSIHRS(params);
    } else if (typeof solveDeterministicSIR !== 'undefined') {
        // Original SIHR model
        detResult = solveDeterministicSIR(params);
    }

    // Update all charts
    if (typeof charts !== 'undefined' && charts.S && charts.I && charts.H && charts.R && charts.D) {
        // SIHRS model with 5 compartments
        const compartments = ['S', 'I', 'H', 'R', 'D'];
        compartments.forEach(comp => {
            if (charts[comp]) {
                charts[comp].data.datasets[0].data = detResult.T.map((t, i) => ({ x: t, y: detResult[`${comp}_prop`][i] }));
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                charts[comp].options.scales.x.max = params.tmax;
                
                // Update dynamic y-axis range (simplified to avoid recursion)
                try {
                    if (typeof calculateDynamicYRange === 'function') {
                        const yRange = calculateDynamicYRange(comp, detResult, results);
                        if (yRange && typeof yRange.min === 'number' && typeof yRange.max === 'number') {
                            charts[comp].options.scales.y.min = yRange.min;
                            charts[comp].options.scales.y.max = yRange.max;
                        }
                    }
                } catch (error) {
                    console.warn('Error updating y-axis range:', error);
                }
                
                charts[comp].update('none');
            }
        });
        
        // Update charts with real data if available (for SIHRS model)
        if (typeof updateChartsWithRealData === 'function') {
            updateChartsWithRealData();
        }
        
        updateStatistics();
        updateStatus('Reset to Run 1 (new batch)', 'stopped');
    } else if (typeof charts !== 'undefined' && charts.I && charts.H) {
        // SIHR model with separate charts
        const compartments = ['I', 'H'];
        const compartmentNames = {
            'I': 'Infected',
            'H': 'Hospitalized'
        };
        
        compartments.forEach(comp => {
            if (charts[comp] && results[currentRun]) {
                // Update both deterministic and stochastic data
                charts[comp].data.datasets[0].data = detResult.T.map((t, i) => ({ x: t, y: detResult[`${comp}_prop`][i] }));
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                charts[comp].options.scales.x.max = params.tmax;
                
                // Update chart title with current run number
                charts[comp].options.plugins.title.text = `${compartmentNames[comp]} Proportion Over Time - Run 1/${totalRuns} (N=${N})`;
                
                // Update dynamic y-axis range for SIHR model (only if function exists)
                if (typeof calculateDynamicYRange === 'function' && detResult && results && results.length > 0) {
                    try {
                        const yRange = calculateDynamicYRange(detResult, results);
                        if (yRange && typeof yRange.min === 'number' && typeof yRange.max === 'number') {
                            charts[comp].options.scales.y.min = yRange.min;
                            charts[comp].options.scales.y.max = yRange.max;
                        }
                    } catch (error) {
                        console.warn('Error updating y-axis range:', error);
                    }
                }
                
                charts[comp].update('none');
            }
        });
        
        updateStatistics();
        updateStatus('Reset to Run 1 (new batch)', 'stopped');
    }
    
    // Reset UI and restart animation
    const playBtn = document.getElementById('playBtn');
    playBtn.innerHTML = '<span class="btn-text">Pause</span>';
    
    // Restart animation
    isRunning = true;
    animate();
}

function changeSpeed() {
    playCoinSound(); // Play coin sound for speed change
    speedIndex = (speedIndex + 1) % speeds.length;
    animationSpeed = speeds[speedIndex];
    const speedNames = ['Fast', 'Normal', 'Slow'];
    document.getElementById('speedBtn').innerHTML = 
        `<span class="btn-text">Speed: ${speedNames[speedIndex]}</span>`;
}

// Update chart with current run
function updateChart() {
    if (typeof charts !== 'undefined') {
        // SIHRS model with multiple charts
        const compartments = ['S', 'I', 'H', 'R', 'D'];
        compartments.forEach(comp => {
            if (charts[comp] && results[currentRun]) {
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                
                // Update dynamic y-axis range for current run (simplified)
                try {
                    if (typeof calculateDynamicYRange === 'function') {
                        const yRange = calculateDynamicYRange(comp, detResult, results);
                        if (yRange && typeof yRange.min === 'number' && typeof yRange.max === 'number') {
                            charts[comp].options.scales.y.min = yRange.min;
                            charts[comp].options.scales.y.max = yRange.max;
                        }
                    }
                } catch (error) {
                    console.warn('Error updating y-axis range:', error);
                }
                
                charts[comp].update('none');
            }
        });
    } else if (typeof charts !== 'undefined' && charts.I && charts.H && results[currentRun]) {
        // SIHR model with separate charts
        const compartments = ['I', 'H'];
        const compartmentNames = {
            'I': 'Infected',
            'H': 'Hospitalized'
        };
        
        compartments.forEach(comp => {
            if (charts[comp]) {
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                // Update chart title with current run number
                charts[comp].options.plugins.title.text = `${compartmentNames[comp]} Proportion Over Time - Run ${currentRun + 1}/${totalRuns} (N=${N})`;
                
                // Update dynamic y-axis range for SIHR model (only if function exists)
                if (typeof calculateDynamicYRange === 'function' && detResult && results && results.length > 0) {
                    try {
                        const yRange = calculateDynamicYRange(detResult, results);
                        if (yRange && typeof yRange.min === 'number' && typeof yRange.max === 'number') {
                            charts[comp].options.scales.y.min = yRange.min;
                            charts[comp].options.scales.y.max = yRange.max;
                        }
                    } catch (error) {
                        console.warn('Error updating y-axis range:', error);
                    }
                }
                
                charts[comp].update('none');
            }
        });
    }
}

// Update statistics display
function updateStatistics() {
    document.getElementById('currentRun').textContent = `${currentRun + 1}/${totalRuns}`;
    document.getElementById('populationSize').textContent = N;
    
    if (typeof calculateR0SIHRS !== 'undefined') {
        // SIHRS model
        params.R_0_value = calculateR0SIHRS(params);
    } else {
        // Original SIHR model
        params.R_0_value = calculateR0(params);
    }
    document.getElementById('R_0-value').textContent = params.R_0_value.toFixed(2);
    
    if (results[currentRun]) {
        const maxH = Math.max(...results[currentRun].H_prop);
        document.getElementById('peakHospitalization').textContent = maxH.toFixed(4);
    }
}

// Animation function just updates the display
function animate() {
    if (!isRunning || !results.length) return;
    
    // Update charts to current run (simplified to avoid recursion)
    if (typeof charts !== 'undefined') {
        // SIHRS model with multiple charts
        const compartments = ['S', 'I', 'H', 'R', 'D'];
        const compartmentNames = {
            'S': 'Susceptible',
            'I': 'Infected', 
            'H': 'Hospitalized',
            'R': 'Recovered',
            'D': 'Death'
        };
        compartments.forEach(comp => {
            if (charts[comp] && results[currentRun]) {
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                // Update chart title with current run number
                charts[comp].options.plugins.title.text = `${compartmentNames[comp]} Proportion Over Time - Run ${currentRun + 1}/${totalRuns} (N=${N})`;
                charts[comp].update('none');
            }
        });
    } else if (typeof charts !== 'undefined' && charts.I && charts.H && results[currentRun]) {
        // SIHR model with separate charts
        const compartments = ['I', 'H'];
        const compartmentNames = {
            'I': 'Infected',
            'H': 'Hospitalized'
        };
        
        compartments.forEach(comp => {
            if (charts[comp]) {
                charts[comp].data.datasets[1].data = results[currentRun].T.map((t, i) => ({ x: t, y: results[currentRun][`${comp}_prop`][i] }));
                // Update chart title with current run number
                charts[comp].options.plugins.title.text = `${compartmentNames[comp]} Proportion Over Time - Run ${currentRun + 1}/${totalRuns} (N=${N})`;
                charts[comp].update('none');
            }
        });
        
        // Update charts with real data if available (for SIHRS model)
        if (typeof updateChartsWithRealData === 'function') {
            updateChartsWithRealData();
        }
    }
    
    updateStatistics();
    updateStatus(`Animation running - Run ${currentRun + 1}/${totalRuns}`, 'running');
    
    // Move to next run with modulo to loop continuously
    currentRun = (currentRun + 1) % totalRuns;
    
    // Schedule next animation frame
    animationTimeout = setTimeout(animate, animationSpeed);
}

// Function to update parameter values and display
function updateParameter(paramName, value) {
    params[paramName] = parseFloat(value);
    document.getElementById(`${paramName}-value`).textContent = value;
    document.getElementById(`${paramName}`).value = value;
    document.getElementById(`${paramName}-number`).value = value;
    
    if (typeof calculateR0SIHRS !== 'undefined') {
        // SIHRS model
        params.R_0_value = calculateR0SIHRS(params);
        document.getElementById('R_0-value').textContent = params.R_0_value.toFixed(2);
    } else {
        // Original SIHR model
        params.R_0_value = calculateR0(params);
        const sigma0 = calculate_thresholds(params)[0];
        const sigma1 = calculate_thresholds(params)[1];
        const tpi = compute_T(params.p2 * params.gamma / (params.p1 * params.beta));
        const h_tpi = compute_h_tpi();
        const sigma2 = calculate_thresholds(params)[2];
        document.getElementById('R_0-value').textContent = params.R_0_value.toFixed(2);
        document.getElementById('sigma0-value').textContent = sigma0.toFixed(2);
        document.getElementById('sigma1-value').textContent = sigma1.toFixed(2);
        document.getElementById('sigma2-value').textContent = sigma2.toFixed(2);
        document.getElementById('tpi-value').textContent = tpi.toFixed(2);
        document.getElementById('h_tpi-value').textContent = h_tpi.toFixed(2);
    }
}

// Matrix-style falling code background
function initMatrixBackground() {
    const matrixBg = document.getElementById('matrixBg');
    const columns = Math.floor(window.innerWidth / 20); // One column every 20px
    const characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@#$%^&*()_+-=[]{}|;:,.<>?βγαρσμλπτφψωΔΣΠΩ∞∫∂∇√∑∏∐∈∉⊂⊃∪∩⊆⊇αβγδεζηθικλμνξοπρστυφχψωR₀SIRH';
    
    // Create matrix columns
    for (let i = 0; i < columns; i++) {
        const column = document.createElement('div');
        column.className = 'matrix-column';
        column.style.left = (i * 20) + 'px';
        column.style.animationDuration = (Math.random() * 15 + 12) + 's'; // 12-27 seconds (very slow)
        column.style.animationDelay = Math.random() * 5 + 's'; // Random delay
        
        // Generate random characters for this column
        let columnText = '';
        const columnLength = Math.floor(Math.random() * 20) + 10; // 10-30 characters
        for (let j = 0; j < columnLength; j++) {
            columnText += characters.charAt(Math.floor(Math.random() * characters.length)) + '<br>';
        }
        column.innerHTML = columnText;
        
        matrixBg.appendChild(column);
    }
    
    // Periodically refresh columns for continuous effect
    setInterval(() => {
        const existingColumns = matrixBg.querySelectorAll('.matrix-column');
        existingColumns.forEach(column => {
            // Randomly regenerate some columns
            if (Math.random() < 0.1) { // 10% chance each interval
                let columnText = '';
                const columnLength = Math.floor(Math.random() * 20) + 10;
                for (let j = 0; j < columnLength; j++) {
                    columnText += characters.charAt(Math.floor(Math.random() * characters.length)) + '<br>';
                }
                column.innerHTML = columnText;
                column.style.animationDuration = (Math.random() * 15 + 12) + 's';
                column.style.animationDelay = '0s';
            }
        });
    }, 15000); // Refresh every 15 seconds (much less frequent)
}

// Draggable Mario Coins
function initDraggableCoins() {
    const coins = document.querySelectorAll('.mario-coin');
    let isDragging = false;
    let currentCoin = null;
    let offsetX, offsetY;

    coins.forEach(coin => {
        coin.addEventListener('mousedown', startDrag);
        coin.addEventListener('touchstart', startDrag, { passive: false });
    });

    document.addEventListener('mousemove', drag);
    document.addEventListener('touchmove', drag, { passive: false });
    document.addEventListener('mouseup', endDrag);
    document.addEventListener('touchend', endDrag);

    function startDrag(e) {
        e.preventDefault();
        isDragging = true;
        currentCoin = e.target;
        
        const rect = currentCoin.getBoundingClientRect();
        const clientX = e.type === 'mousedown' ? e.clientX : e.touches[0].clientX;
        const clientY = e.type === 'mousedown' ? e.clientY : e.touches[0].clientY;
        
        offsetX = clientX - rect.left;
        offsetY = clientY - rect.top;
        
        currentCoin.style.zIndex = '1000';
        currentCoin.style.cursor = 'grabbing';
        
        // Play coin sound when starting to drag
        playCoinSound();
    }

    function drag(e) {
        if (!isDragging || !currentCoin) return;
        e.preventDefault();
        
        const clientX = e.type === 'mousemove' ? e.clientX : e.touches[0].clientX;
        const clientY = e.type === 'mousemove' ? e.clientY : e.touches[0].clientY;
        
        const x = clientX - offsetX;
        const y = clientY - offsetY;
        
        // Keep coin within viewport bounds
        const maxX = window.innerWidth - currentCoin.offsetWidth;
        const maxY = window.innerHeight - currentCoin.offsetHeight;
        
        currentCoin.style.left = Math.max(0, Math.min(x, maxX)) + 'px';
        currentCoin.style.top = Math.max(0, Math.min(y, maxY)) + 'px';
    }

    function endDrag(e) {
        if (!isDragging || !currentCoin) return;
        
        isDragging = false;
        currentCoin.style.zIndex = '10';
        currentCoin.style.cursor = 'grab';
        currentCoin = null;
        
        // Play coin sound when dropping
        playCoinSound();
    }
}

// Draggable Dragon Balls
function initDraggableDragonBalls() {
    const balls = document.querySelectorAll('.dragon-ball');
    let isDragging = false;
    let currentBall = null;
    let offsetX, offsetY;

    balls.forEach(ball => {
        ball.addEventListener('mousedown', startDrag);
        ball.addEventListener('touchstart', startDrag, { passive: false });
    });

    document.addEventListener('mousemove', drag);
    document.addEventListener('touchmove', drag, { passive: false });
    document.addEventListener('mouseup', endDrag);
    document.addEventListener('touchend', endDrag);

    function startDrag(e) {
        e.preventDefault();
        isDragging = true;
        currentBall = e.target.closest('.dragon-ball');
        if (!currentBall) return;
        const rect = currentBall.getBoundingClientRect();
        const clientX = e.type === 'mousedown' ? e.clientX : e.touches[0].clientX;
        const clientY = e.type === 'mousedown' ? e.clientY : e.touches[0].clientY;
        offsetX = clientX - rect.left;
        offsetY = clientY - rect.top;
        currentBall.style.zIndex = '1000';
        currentBall.style.cursor = 'grabbing';
        playCoinSound();
    }

    function drag(e) {
        if (!isDragging || !currentBall) return;
        e.preventDefault();
        const clientX = e.type === 'mousemove' ? e.clientX : e.touches[0].clientX;
        const clientY = e.type === 'mousemove' ? e.clientY : e.touches[0].clientY;
        const x = clientX - offsetX;
        const y = clientY - offsetY;
        // Keep ball within viewport bounds
        const maxX = window.innerWidth - currentBall.offsetWidth;
        const maxY = window.innerHeight - currentBall.offsetHeight;
        currentBall.style.left = Math.max(0, Math.min(x, maxX)) + 'px';
        currentBall.style.top = Math.max(0, Math.min(y, maxY)) + 'px';
    }

    function endDrag(e) {
        if (!isDragging || !currentBall) return;
        isDragging = false;
        currentBall.style.zIndex = '10';
        currentBall.style.cursor = 'grab';
        currentBall = null;
        playCoinSound();
    }
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        updateStatus,
        toggleAnimation,
        resetAnimation,
        changeSpeed,
        updateChart,
        updateStatistics,
        animate,
        updateParameter,
        initMatrixBackground,
        initDraggableCoins,
        initDraggableDragonBalls
    };
} 