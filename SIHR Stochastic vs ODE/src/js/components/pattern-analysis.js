// Pattern Analysis Component
// Functions for analyzing and grouping similar simulation runs

// Function to calculate key characteristics of a run for a specific compartment
function calculateRunCharacteristics(run, compartment = 'H') {
    const compartmentData = run[`${compartment}_prop`];
    const T = run.T;
    
    // Find peak
    const maxValue = Math.max(...compartmentData);
    const peakTime = T[compartmentData.indexOf(maxValue)];
    
    // Count major fluctuations (peaks and troughs)
    let fluctuations = 0;
    let prevDiff = 0;
    for (let i = 1; i < compartmentData.length; i++) {
        const diff = compartmentData[i] - compartmentData[i-1];
        if (prevDiff * diff < 0) { // Sign change indicates peak or trough
            fluctuations++;
        }
        prevDiff = diff;
    }
    
    return {
        maxValue,
        peakTime,
        fluctuations,
        curve: compartmentData, // Store the actual curve for shape comparison
        compartment: compartment
    };
}

// Function to group similar runs
function groupSimilarRuns(runCharacteristics) {
    const groups = [];
    
    runCharacteristics.forEach((char, index) => {
        let addedToGroup = false;
        
        // Try to add to existing group
        for (let group of groups) {
            if (areRunsSimilar(char, group.representative)) {
                group.runs.push(index);
                addedToGroup = true;
                break;
            }
        }
        
        // Create new group if no match found
        if (!addedToGroup) {
            groups.push({
                representative: char,
                runs: [index]
            });
        }
    });
    
    return groups;
}

// Function to compare two runs for similarity
function areRunsSimilar(run1, run2, thresholds = {
    peakHeight: 0.1,    // 10% difference allowed
    peakTime: 2,        // 2 days difference allowed
    fluctuations: 1,    // 1 fluctuation difference allowed
    shapeSimilarity: 0.8 // 80% shape similarity required
}) {
    // Compare peak heights
    const heightDiff = Math.abs(run1.maxH - run2.maxH) / Math.max(run1.maxH, run2.maxH);
    if (heightDiff > thresholds.peakHeight) return false;
    
    // Compare peak times
    const timeDiff = Math.abs(run1.peakTime - run2.peakTime);
    if (timeDiff > thresholds.peakTime) return false;
    
    // Compare number of fluctuations
    const fluctDiff = Math.abs(run1.fluctuations - run2.fluctuations);
    if (fluctDiff > thresholds.fluctuations) return false;
    
    // Compare overall shape using correlation
    const shapeSimilarity = calculateShapeSimilarity(run1.curve, run2.curve);
    if (shapeSimilarity < thresholds.shapeSimilarity) return false;
    
    return true;
}

// Function to calculate shape similarity between two curves
function calculateShapeSimilarity(curve1, curve2) {
    // Normalize curves to same length
    const length = Math.min(curve1.length, curve2.length);
    const norm1 = curve1.slice(0, length);
    const norm2 = curve2.slice(0, length);
    
    // Calculate correlation coefficient
    const mean1 = norm1.reduce((a, b) => a + b) / length;
    const mean2 = norm2.reduce((a, b) => a + b) / length;
    
    let numerator = 0;
    let denom1 = 0;
    let denom2 = 0;
    
    for (let i = 0; i < length; i++) {
        const diff1 = norm1[i] - mean1;
        const diff2 = norm2[i] - mean2;
        numerator += diff1 * diff2;
        denom1 += diff1 * diff1;
        denom2 += diff2 * diff2;
    }
    
    return Math.abs(numerator / Math.sqrt(denom1 * denom2));
}

// Function to create pattern cards
function createPatternCards(groups, results, params, compartment = 'H') {
    const container = document.getElementById('patternCards');
    container.innerHTML = '';

    groups.forEach((group, index) => {
        const card = document.createElement('div');
        card.className = 'pattern-card';
        card.innerHTML = `
            <div class="pattern-header">
                <div class="pattern-title">Pattern ${index + 1}</div>
                <div class="pattern-count">${group.runs.length} runs (${((group.runs.length / results.length) * 100).toFixed(1)}%)</div>
            </div>
            <div class="pattern-details">
                <div class="detail-item">
                    <span class="detail-label">Peak Height</span>
                    <span class="detail-value">${group.representative.maxValue.toFixed(4)}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">Peak Time</span>
                    <span class="detail-value">${group.representative.peakTime.toFixed(1)} days</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">Fluctuations</span>
                    <span class="detail-value">${group.representative.fluctuations}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">View Runs</span>
                    <span class="detail-value">
                        <a class="pattern-link" onclick="showPatternDetails(${index}, ${JSON.stringify(group.runs)}, '${compartment}')">
                            View ${group.runs.length} runs →
                        </a>
                    </span>
                </div>
            </div>
            <div class="pattern-chart">
                <canvas id="pattern${index}Chart"></canvas>
            </div>
        `;
        container.appendChild(card);

        // Create mini chart for this pattern
        const ctx = document.getElementById(`pattern${index}Chart`).getContext('2d');
        new Chart(ctx, {
            type: 'line',
            data: {
                labels: Array.from({length: group.representative.curve.length}, (_, i) => i),
                datasets: [{
                    data: group.representative.curve,
                    borderColor: `hsl(${index * 360 / groups.length}, 70%, 50%)`,
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0.2
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: { display: false },
                    tooltip: { enabled: false }
                },
                scales: {
                    x: { display: false },
                    y: { display: false }
                }
            }
        });
    });
}

// Function to create pattern distribution chart
function createPatternChart(groups) {
    const ctx = document.getElementById('patternChart').getContext('2d');
    
    // Properly destroy existing chart if it exists
    if (window.patternChart instanceof Chart) {
        window.patternChart.destroy();
    }
    
    // Create new chart
    window.patternChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: groups.map((_, i) => `Pattern ${i + 1}`),
            datasets: [{
                data: groups.map(g => g.runs.length),
                backgroundColor: groups.map((_, i) => 
                    `hsl(${i * 360 / groups.length}, 70%, 50%)`
                )
            }]
        },
        options: {
            responsive: true,
            plugins: {
                legend: { display: false },
                title: {
                    display: true,
                    text: 'Number of Runs per Pattern'
                }
            },
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Number of Runs'
                    }
                }
            }
        }
    });
}

// Function to update pattern analysis for a specific compartment
function updatePatternAnalysis(results, compartment = 'H') {
    const patternSection = document.getElementById('patternAnalysis');
    
    if (!results.length) {
        patternSection.style.display = 'none';
        return;
    }

    try {
        // Group the runs for the specific compartment
        const runCharacteristics = results.map(run => calculateRunCharacteristics(run, compartment));
        window.runGroups = groupSimilarRuns(runCharacteristics);
        
        // Update summary statistics
        const totalPatterns = window.runGroups.length;
        const mostCommonPattern = Math.max(...window.runGroups.map(g => g.runs.length)) / results.length * 100;
        const avgPeakHeight = window.runGroups.reduce((sum, g) => sum + g.representative.maxValue, 0) / totalPatterns;
        const avgPeakTime = window.runGroups.reduce((sum, g) => sum + g.representative.peakTime, 0) / totalPatterns;

        document.getElementById('totalPatterns').textContent = totalPatterns;
        document.getElementById('mostCommonPattern').textContent = `${mostCommonPattern.toFixed(1)}%`;
        document.getElementById('avgPeakHeight').textContent = avgPeakHeight.toFixed(4);
        document.getElementById('avgPeakTime').textContent = avgPeakTime.toFixed(1);

        // Create pattern cards with compartment info
        createPatternCards(window.runGroups, results, params, compartment);

        // Create pattern distribution chart
        createPatternChart(window.runGroups);

        patternSection.style.display = 'block';
    } catch (error) {
        console.error('Error in pattern analysis:', error);
        updateStatus(`Error in pattern analysis: ${error.message}`, 'stopped');
    }
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        calculateRunCharacteristics,
        groupSimilarRuns,
        areRunsSimilar,
        calculateShapeSimilarity,
        createPatternCards,
        createPatternChart,
        updatePatternAnalysis
    };
} 