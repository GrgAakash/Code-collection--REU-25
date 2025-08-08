// Download Utilities Component
// Functions for downloading charts and analysis results

// Function to download pattern as PNG for a specific compartment
function downloadPatternAsPNG(patternIndex, runIndices, results, params, compartment = 'H') {
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    canvas.width = 1200;
    canvas.height = 800;
    
    // Create a temporary chart for the combined view
    const tempDiv = document.createElement('div');
    tempDiv.style.position = 'absolute';
    tempDiv.style.left = '-9999px';
    document.body.appendChild(tempDiv);
    
    const tempCanvas = document.createElement('canvas');
    tempCanvas.width = 1200;
    tempCanvas.height = 600;
    tempDiv.appendChild(tempCanvas);
    
    const tempCtx = tempCanvas.getContext('2d');
    
    // Get deterministic solution for comparison
    const detResult = typeof solveDeterministicSIHRS !== 'undefined' ? solveDeterministicSIHRS(params) : solveDeterministicSIR(params);
    
    // Create combined chart for the specific compartment only
    const compartmentColors = {
        'S': '#FF6B6B',
        'I': '#4ECDC4', 
        'H': '#45B7D1',
        'R': '#96CEB4',
        'D': '#DDA0DD'
    };
    
    const datasets = [];
    
    // Add deterministic dataset for the specific compartment
    datasets.push({
        label: `Deterministic ${compartment.toLowerCase()}`,
        data: detResult.T.map((t, j) => ({ x: t, y: detResult[`${compartment}_prop`][j] })),
        borderColor: compartmentColors[compartment] || '#e53e3e',
        backgroundColor: `${compartmentColors[compartment] || '#e53e3e'}20`,
        borderDash: [8, 4],
        borderWidth: 3,
        pointRadius: 0,
        fill: false,
        tension: 0.1
    });
    
    // Add stochastic datasets for each run (only the specific compartment)
    runIndices.forEach((runIndex, i) => {
        datasets.push({
            label: `Stochastic ${compartment.toLowerCase()} (Run ${runIndex + 1})`,
            data: results[runIndex].T.map((t, j) => ({ x: t, y: results[runIndex][`${compartment}_prop`][j] })),
            borderColor: compartmentColors[compartment] || `hsl(${15 + i * 20}, 80%, 60%)`,
            backgroundColor: `${compartmentColors[compartment] || `hsl(${15 + i * 20}, 80%, 60%)`}20`,
            borderWidth: 2,
            pointRadius: 0,
            tension: 0.2,
            fill: false
        });
    });
    
    const combinedChart = new Chart(tempCtx, {
        type: 'line',
        data: {
            labels: results[0].T.map(t => t.toFixed(1)),
            datasets: datasets
        },
        options: {
            responsive: false,
            maintainAspectRatio: false,
            layout: {
                padding: {
                    bottom: 30,
                    left: 20,
                    right: 20
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: `Pattern ${patternIndex + 1} - All ${runIndices.length} Runs vs Deterministic`,
                    font: { size: 18, weight: 'bold' }
                },
                legend: {
                    display: true,
                    position: 'top',
                    labels: {
                        usePointStyle: true,
                        padding: 15,
                        font: { size: 12 }
                    }
                }
            },
            scales: {
                x: {
                    type: 'linear',
                    title: { 
                        display: true, 
                        text: 'Time (days)',
                        font: { size: 14, weight: 'bold' }
                    },
                    ticks: { font: { size: 12 } }
                },
                y: {
                    title: { 
                        display: true, 
                        text: 'Proportion',
                        font: { size: 14, weight: 'bold' }
                    },
                    beginAtZero: true,
                    ticks: { font: { size: 12 } }
                }
            }
        }
    });
    
    // Wait for chart to render then download
    setTimeout(() => {
        const link = document.createElement('a');
        link.download = `pattern_${patternIndex + 1}_${compartment}_runs_with_deterministic.png`;
        link.href = tempCanvas.toDataURL('image/png', 1.0);
        link.click();
        
        // Cleanup
        combinedChart.destroy();
        document.body.removeChild(tempDiv);
    }, 500);
}

// Function to download selected individual run (single selection)
function downloadSelectedIndividualRun(allRunIndices, results, params, compartment = 'H') {
    // Get the selected radio button
    const selectedRadio = document.querySelector('input[name="selectedRun"]:checked');
    
    if (!selectedRadio) {
        alert('Please select a run to download.');
        return;
    }
    
    const runIndex = parseInt(selectedRadio.value);
    const run = results[runIndex];
    
    // Create canvas and download single run
    const canvas = document.createElement('canvas');
    canvas.width = 800;
    canvas.height = 400;
    const ctx = canvas.getContext('2d');
    
    // Get deterministic solution for comparison
    const detResult = typeof solveDeterministicSIHRS !== 'undefined' ? solveDeterministicSIHRS(params) : solveDeterministicSIR(params);
    
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
            borderColor: compartmentColors[compartment] || '#e53e3e',
            backgroundColor: `${compartmentColors[compartment] || '#e53e3e'}20`,
            borderDash: [8, 4],
            borderWidth: 3,
            pointRadius: 0,
            fill: false,
            tension: 0.1
        },
        {
            label: `Stochastic ${compartment.toLowerCase()}`,
            data: run.T.map((t, j) => ({ x: t, y: run[`${compartment}_prop`][j] })),
            borderColor: compartmentColors[compartment] || '#ff6b6b',
            backgroundColor: `${compartmentColors[compartment] || '#ff6b6b'}20`,
            borderWidth: 2,
            pointRadius: 0,
            fill: false,
            tension: 0.2
        }
    ];
    
    const chart = new Chart(ctx, {
        type: 'line',
        data: {
            datasets: datasets
        },
        options: {
            responsive: false,
            maintainAspectRatio: false,
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
                    ticks: { font: { size: 11 } }
                },
                y: {
                    title: { 
                        display: true, 
                        text: 'Proportion',
                        font: { size: 13, weight: 'bold' }
                    },
                    beginAtZero: true,
                    ticks: { font: { size: 11 } }
                }
            }
        }
    });
    
    setTimeout(() => {
        const link = document.createElement('a');
        link.download = `run_${runIndex + 1}_${compartment}_deterministic_vs_stochastic.png`;
        link.href = canvas.toDataURL('image/png', 1.0);
        link.click();
        chart.destroy();
    }, 300);
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        downloadPatternAsPNG,
        downloadSelectedIndividualRun
    };
} 