// Real World Data Handler Component
// Functions for handling CSV data uploads and displaying real-world data on charts

// Global variables for real data
let realData = null;
let showRealData = false;
let realDataCompartment = 'I'; // Default to Infected
let cumulativeDeathsData = null; // Store cumulative deaths data (Julia-style)
let hospitalizationData = null; // Store hospitalization data (Julia-style)

// Function to parse CSV data
function parseCSV(csvText) {
    const lines = csvText.trim().split('\n');
    const headers = lines[0].split(',').map(h => h.trim());
    const data = [];
    
    for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => v.trim());
        const row = {};
        headers.forEach((header, index) => {
            row[header] = values[index];
        });
        data.push(row);
    }
    
    return { headers, data };
}

// Function to validate and process uploaded data
function processRealData(file) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        
        reader.onload = function(e) {
            try {
                const csvText = e.target.result;
                const { headers, data } = parseCSV(csvText);
                
                // Show format requirements if validation fails
                const formatInfo = `
📋 Supported CSV Formats:

1. **Simple Time Series** (like your simulation data):
   time,infected,hospitalized,recovered,deaths
   0,0.01,0.0,0.0,0.0
   1,0.015,0.001,0.0,0.0
   2,0.025,0.002,0.001,0.0

2. **Real-World Data** (like Carson City data):
   date,cases,deaths,daily_deaths,cumulative_deaths
   2020-03-25,1,0,0,0
   2020-03-26,2,0,0,0
   2020-03-27,3,0,0,0

3. **Active Cases Data**:
   date,active_cases,hospitalized,recovered
   2020-03-25,1,0,0
   2020-03-26,2,0,0

4. **Hospitalization Data** (like Carson City hospital data):
   collection_week,total_adult_and_pediatric_covid_patients
   8/2/20,7.9
   8/9/20,5.1
   8/16/20,8.3

Values can be:
- Proportions (0.0 to 1.0)
- Percentages (0 to 100)
- Raw counts (will be converted to proportions)
                `.trim();
                
                // Check for different data formats
                const hasDate = headers.some(h => h.toLowerCase().includes('date'));
                const hasTime = headers.some(h => h.toLowerCase().includes('time') || h.toLowerCase().includes('day'));
                const hasCases = headers.some(h => h.toLowerCase().includes('case'));
                const hasDeaths = headers.some(h => h.toLowerCase().includes('death'));
                const hasActive = headers.some(h => h.toLowerCase().includes('active'));
                const hasInfected = headers.some(h => h.toLowerCase().includes('infected'));
                const hasHospitalized = headers.some(h => h.toLowerCase().includes('hospitalized'));
                const hasRecovered = headers.some(h => h.toLowerCase().includes('recovered'));
                const hasSusceptible = headers.some(h => h.toLowerCase().includes('susceptible'));
                const hasCollectionWeek = headers.some(h => h.toLowerCase().includes('collection_week'));
                const hasCovidPatients = headers.some(h => h.toLowerCase().includes('total_adult_and_pediatric_covid_patients'));
                
                // Determine data format
                let dataFormat = 'simple';
                let timeColumn = null;
                let dataColumns = [];
                
                if (hasCollectionWeek && hasCovidPatients) {
                    // Hospitalization data format (like Julia implementation)
                    dataFormat = 'hospitalization';
                    timeColumn = headers.find(h => h.toLowerCase().includes('collection_week'));
                    dataColumns.push(headers.find(h => h.toLowerCase().includes('total_adult_and_pediatric_covid_patients')));
                } else if (hasDate) {
                    dataFormat = 'real_world';
                    timeColumn = headers.find(h => h.toLowerCase().includes('date'));
                    
                    // Find data columns based on what's available (like Julia implementation)
                    if (hasCases) dataColumns.push(headers.find(h => h.toLowerCase().includes('case')));
                    if (hasDeaths) dataColumns.push(headers.find(h => h.toLowerCase().includes('death')));
                    if (hasActive) dataColumns.push(headers.find(h => h.toLowerCase().includes('active')));
                    if (hasHospitalized) dataColumns.push(headers.find(h => h.toLowerCase().includes('hospitalized')));
                    if (hasRecovered) dataColumns.push(headers.find(h => h.toLowerCase().includes('recovered')));
                    
                    // Check for specific Carson City CSV columns like Julia implementation
                    const hasMovingAvg = headers.some(h => h.toLowerCase().includes('moving_avg'));
                    const hasCumulativeDeaths = headers.some(h => h.toLowerCase().includes('cumulative_deaths'));
                    const hasDailyDeaths = headers.some(h => h.toLowerCase().includes('daily_deaths'));
                    
                    if (hasMovingAvg) dataColumns.push(headers.find(h => h.toLowerCase().includes('moving_avg')));
                    if (hasCumulativeDeaths) dataColumns.push(headers.find(h => h.toLowerCase().includes('cumulative_deaths')));
                    if (hasDailyDeaths) dataColumns.push(headers.find(h => h.toLowerCase().includes('daily_deaths')));
                } else if (hasTime) {
                    dataFormat = 'simple';
                    timeColumn = headers.find(h => h.toLowerCase().includes('time') || h.toLowerCase().includes('day'));
                    
                    // Find all data columns
                    dataColumns = headers.filter(h => {
                        const lowerH = h.toLowerCase();
                        return lowerH.includes('infected') || lowerH.includes('hospitalized') || 
                               lowerH.includes('recovered') || lowerH.includes('death') || 
                               lowerH.includes('susceptible') || lowerH.includes('value');
                    });
                }
                
                if (!timeColumn) {
                    reject(new Error(`❌ No time/date column found!\n\n${formatInfo}`));
                    return;
                }
                
                if (dataColumns.length === 0) {
                    reject(new Error(`❌ No valid data columns found!\n\n${formatInfo}`));
                    return;
                }
                
                // Process data based on format
                let processedData = [];
                let compartment = 'I'; // Default to Infected
                let population = 56000; // Default population (Carson City)
                
                if (dataFormat === 'hospitalization') {
                    // Handle hospitalization data (like Julia implementation)
                    const valueCol = dataColumns[0];
                    
                    console.log(`🔍 Processing Carson City hospitalization data with column: ${valueCol}`);
                    
                    // Parse m/d/yy dates and filter out zero values (like Julia implementation)
                    const simulationStartDate = new Date('2020-03-25'); // March 25, 2020
                    const processedDataWithDates = data
                        .filter(row => row[timeColumn] && row[valueCol] && parseFloat(row[valueCol]) > 0) // Filter out zeros
                        .map(row => {
                            const dateStr = row[timeColumn];
                            
                            // Parse m/d/yy format and fix year to 20xx (like Julia implementation)
                            const [month, day, year] = dateStr.split('/');
                            const fullYear = parseInt(year) < 50 ? 2000 + parseInt(year) : 1900 + parseInt(year);
                            const date = new Date(fullYear, parseInt(month) - 1, parseInt(day));
                            
                            const daysFromStart = (date - simulationStartDate) / (1000 * 60 * 60 * 24);
                            const value = parseFloat(row[valueCol]);
                            
                            return {
                                time: daysFromStart,
                                value: value,
                                originalDate: dateStr,
                                parsedDate: date
                            };
                        })
                        .filter(point => !isNaN(point.time) && !isNaN(point.value) && point.time >= 0)
                        .sort((a, b) => a.time - b.time);
                    
                    processedData = processedDataWithDates;
                    compartment = 'H'; // Hospitalization compartment
                    
                    // Store hospitalization data globally (like Julia implementation)
                    storeHospitalizationData(processedDataWithDates, population);
                    
                    const maxHosp = Math.max(...processedData.map(p => p.value));
                    console.log(`✅ Peak hospitalization: ${maxHosp.toFixed(1)} (${(maxHosp/population).toFixed(6)} proportion)`);
                    
                } else if (dataFormat === 'real_world') {
                    // Handle real-world data with dates
                    const valueCol = dataColumns[0]; // Use first available data column
                    
                    console.log(`🔍 Processing Carson City data with column: ${valueCol}`);
                    
                    // Convert dates to days from start
                    const startDate = new Date(data[0][timeColumn]);
                    const processedDataWithDates = data
                        .filter(row => row[timeColumn] && row[valueCol])
                        .map(row => {
                            const date = new Date(row[timeColumn]);
                            const daysFromStart = (date - startDate) / (1000 * 60 * 60 * 24);
                            const value = parseFloat(row[valueCol]);
                            
                            // Create data point with all available columns (like Julia implementation)
                            const dataPoint = {
                                time: daysFromStart,
                                value: value,
                                originalDate: row[timeColumn]
                            };
                            
                            // Add specific columns if they exist (like Julia implementation)
                            if (row.cases !== undefined) dataPoint.cases = parseFloat(row.cases);
                            if (row.deaths !== undefined) dataPoint.deaths = parseFloat(row.deaths);
                            if (row.daily_deaths !== undefined) dataPoint.daily_deaths = parseFloat(row.daily_deaths);
                            if (row.cumulative_deaths !== undefined) dataPoint.cumulative_deaths = parseFloat(row.cumulative_deaths);
                            if (row.moving_avg_7day !== undefined) dataPoint.moving_avg_7day = parseFloat(row.moving_avg_7day);
                            if (row.active_cases !== undefined) dataPoint.active_cases = parseFloat(row.active_cases);
                            if (row.cumulative_cases !== undefined) dataPoint.cumulative_cases = parseFloat(row.cumulative_cases);
                            
                            return dataPoint;
                        })
                        .filter(point => !isNaN(point.time) && !isNaN(point.value))
                        .sort((a, b) => a.time - b.time);
                    
                    // 🔧 FIX: Calculate active cases if we're processing cumulative cases data
                    if (valueCol.toLowerCase().includes('case') && !valueCol.toLowerCase().includes('active')) {
                        console.log(`📊 Converting cumulative cases to active cases (Julia-style)`);
                        
                        // Extract cumulative cases
                        const cumulativeCases = processedDataWithDates.map(point => point.value);
                        
                        // Calculate active cases with 14-day recovery period (EXACTLY like Julia)
                        const recoveryDays = 14;
                        const cumulativeShifted = [];
                        
                        // Add zeros for first 14 days
                        for (let i = 0; i < recoveryDays; i++) {
                            cumulativeShifted.push(0);
                        }
                        // Add shifted cumulative cases
                        for (let i = 0; i < cumulativeCases.length - recoveryDays; i++) {
                            cumulativeShifted.push(cumulativeCases[i]);
                        }
                        
                        // Calculate active cases = cumulative - cumulative_shifted
                        processedData = processedDataWithDates.map((point, i) => {
                            const shifted = cumulativeShifted[i] || 0;
                            const activeCases = Math.max(point.value - shifted, 0);
                            return {
                                ...point,
                                value: activeCases, // Replace cumulative with active cases
                                originalCumulative: point.value, // Keep original for reference
                                activeCases: activeCases
                            };
                        });
                        
                        const maxActive = Math.max(...processedData.map(p => p.value));
                        console.log(`✅ Peak active cases: ${maxActive.toFixed(0)} (${(maxActive/population).toFixed(6)} proportion)`);
                        
                        // 🆕 JULIA-STYLE: Also process cumulative deaths if available in same CSV
                        if (processedDataWithDates.length > 0 && processedDataWithDates[0].deaths !== undefined) {
                            console.log(`📊 Also processing cumulative deaths from same CSV (Julia-style)`);
                            storeCumulativeDeathsData(processedDataWithDates, population);
                        }
                        
                    } else {
                        processedData = processedDataWithDates;
                    }
                    
                    // Determine compartment based on column name (like Julia implementation)
                    if (valueCol.toLowerCase().includes('moving_avg') || valueCol.toLowerCase().includes('daily_deaths')) {
                        compartment = 'D';
                    } else if (valueCol.toLowerCase().includes('cumulative_deaths')) {
                        compartment = 'D';
                    } else if (valueCol.toLowerCase().includes('death')) {
                        compartment = 'D';
                    } else if (valueCol.toLowerCase().includes('hospitalized') || valueCol.toLowerCase().includes('hosp')) {
                        compartment = 'H';
                    } else if (valueCol.toLowerCase().includes('recovered') || valueCol.toLowerCase().includes('rec')) {
                        compartment = 'R';
                    } else if (valueCol.toLowerCase().includes('susceptible') || valueCol.toLowerCase().includes('sus')) {
                        compartment = 'S';
                    } else if (valueCol.toLowerCase().includes('case') || valueCol.toLowerCase().includes('active')) {
                        compartment = 'I';
                    }
                    
                } else {
                    // Handle simple time series data
                    const valueCol = dataColumns[0];
                    
                    processedData = data
                        .filter(row => row[timeColumn] && row[valueCol])
                        .map(row => ({
                            time: parseFloat(row[timeColumn]),
                            value: parseFloat(row[valueCol])
                        }))
                        .filter(point => !isNaN(point.time) && !isNaN(point.value))
                        .sort((a, b) => a.time - b.time);
                    
                    // Determine compartment type based on column name
                    if (valueCol.toLowerCase().includes('hospitalized') || valueCol.toLowerCase().includes('hosp')) {
                        compartment = 'H';
                    } else if (valueCol.toLowerCase().includes('susceptible') || valueCol.toLowerCase().includes('sus')) {
                        compartment = 'S';
                    } else if (valueCol.toLowerCase().includes('recovered') || valueCol.toLowerCase().includes('rec')) {
                        compartment = 'R';
                    } else if (valueCol.toLowerCase().includes('death') || valueCol.toLowerCase().includes('dead')) {
                        compartment = 'D';
                    }
                }
                
                if (processedData.length === 0) {
                    reject(new Error('❌ No valid data points found. Check that your data contains numeric values only.'));
                    return;
                }
                
                // Store original values but don't convert to proportions yet (like Julia implementation)
                // Population scaling will be applied AFTER interpolation
                const maxValue = Math.max(...processedData.map(p => p.value));
                if (maxValue > 1.0 && maxValue <= 1000000) { // Likely raw counts
                    processedData = processedData.map(point => ({
                        ...point,
                        originalValue: point.value // Store original for reference
                    }));
                }
                
                resolve({
                    data: processedData,
                    compartment: compartment,
                    fileName: file.name,
                    totalPoints: processedData.length,
                    dataFormat: dataFormat,
                    timeRange: {
                        min: Math.min(...processedData.map(p => p.time)),
                        max: Math.max(...processedData.map(p => p.time))
                    },
                    valueRange: {
                        min: Math.min(...processedData.map(p => p.value)),
                        max: Math.max(...processedData.map(p => p.value))
                    },
                    originalValueRange: {
                        min: Math.min(...processedData.map(p => p.originalValue || p.value)),
                        max: Math.max(...processedData.map(p => p.originalValue || p.value))
                    }
                });
                
            } catch (error) {
                reject(new Error('Error parsing CSV file: ' + error.message));
            }
        };
        
        reader.onerror = function() {
            reject(new Error('Error reading file'));
        };
        
        reader.readAsText(file);
    });
}

// Function to upload and process real data
async function uploadRealData() {
    const fileInput = document.getElementById('realDataUpload');
    const file = fileInput.files[0];
    
    if (!file) {
        updateUploadStatus('Please select a CSV file', 'error');
        return;
    }
    
    if (!file.name.toLowerCase().endsWith('.csv')) {
        updateUploadStatus('Please select a CSV file', 'error');
        return;
    }
    
    updateUploadStatus('Processing data...', '');
    
    try {
        const result = await processRealData(file);
        realData = result;
        realDataCompartment = result.compartment;
        
        // Show success message and info
        updateUploadStatus(`✅ Data loaded successfully!`, 'success');
        showDataInfo(result);
        
        // Show control buttons
        document.getElementById('toggleRealDataBtn').style.display = 'inline-block';
        document.getElementById('clearRealDataBtn').style.display = 'inline-block';
        document.getElementById('updateInitialsBtn').style.display = 'inline-block';
        
        // Ensure simulation results are available for this compartment
        ensureSimulationResultsForCompartment(result.compartment);
        
        // Update charts with real data using Julia-style processing
        if (showRealData) {
            processRealDataWithJuliaAnalysis();
        }
        
        // Play success sound
        if (typeof playCoinSound === 'function') {
            playCoinSound();
        }
        
    } catch (error) {
        updateUploadStatus(`❌ ${error.message}`, 'error');
        console.error('Upload error:', error);
    }
}

// Function to update upload status
function updateUploadStatus(message, className) {
    const statusEl = document.getElementById('uploadStatus');
    statusEl.textContent = message;
    statusEl.className = `upload-status ${className}`;
}

// Function to show data information (concise version)
function showDataInfo(dataInfo) {
    const infoEl = document.getElementById('dataInfo');
    const textEl = document.getElementById('dataInfoText');
    
    const compartmentNames = {
        'S': 'Susceptible',
        'I': 'Infected', 
        'H': 'Hospitalized',
        'R': 'Recovered',
        'D': 'Deaths'
    };
    
    const compartmentName = compartmentNames[dataInfo.compartment] || dataInfo.compartment;
    
    // Concise single-line format
    const infoHTML = `📁 ${dataInfo.fileName} | ${dataInfo.totalPoints} points | ${compartmentName} (${dataInfo.compartment}) | ✅ Ready`;
    
    textEl.innerHTML = infoHTML;
    infoEl.style.display = 'block';
}

// Function to toggle real data visibility
function toggleRealData() {
    showRealData = !showRealData;
    const btn = document.getElementById('toggleRealDataBtn');
    
    if (showRealData) {
        btn.innerHTML = '<span class="btn-text">🙈 Hide Real Data</span>';
        updateUploadStatus('👁️ Real data is now visible on charts (Julia-style processing)', 'success');
    } else {
        btn.innerHTML = '<span class="btn-text">👁️ Show Real Data</span>';
        updateUploadStatus('🙈 Real data is now hidden', '');
    }
    
    // Use Julia-style processing instead of basic processing
    if (showRealData) {
        processRealDataWithJuliaAnalysis();
    } else {
        updateChartsWithRealData(); // This will remove the data
    }
    
    // Play sound
    if (typeof playCoinSound === 'function') {
        playCoinSound();
    }
}

// Function to clear real data
function clearRealData() {
    realData = null;
    showRealData = false;
    cumulativeDeathsData = null; // Clear cumulative deaths data too
    hospitalizationData = null; // Clear hospitalization data too
    
    // Hide buttons and info
    document.getElementById('toggleRealDataBtn').style.display = 'none';
    document.getElementById('clearRealDataBtn').style.display = 'none';
    document.getElementById('updateInitialsBtn').style.display = 'none';
    document.getElementById('dataInfo').style.display = 'none';
    
    // Clear status
    updateUploadStatus('', '');
    
    // Update charts without real data
    updateChartsWithRealData();
    
    // Play sound
    if (typeof playCoinSound === 'function') {
        playCoinSound();
    }
}

// Function to interpolate data to common time grid (like Julia implementation)
function interpolateToGrid(data, t_grid, extrapolate = true) {
    if (!data || data.length === 0) return new Array(t_grid.length).fill(0);
    
    // Create interpolation function
    const interpolate = (x, y, x_new) => {
        const result = [];
        for (let i = 0; i < x_new.length; i++) {
            const xi = x_new[i];
            
            // Find the two points that bracket xi
            let j = 0;
            while (j < x.length - 1 && x[j + 1] < xi) j++;
            
            if (j === x.length - 1 && xi > x[x.length - 1]) {
                // Beyond the end
                if (extrapolate) {
                    result.push(y[y.length - 1]);
                } else {
                    result.push(NaN); // Like Julia's extrapolation_bc=NaN
                }
            } else if (j === 0 && xi < x[0]) {
                // Before the beginning  
                if (extrapolate) {
                    result.push(y[0]);
                } else {
                    result.push(NaN); // Like Julia's extrapolation_bc=NaN
                }
            } else {
                // Linear interpolation
                const x0 = x[j];
                const x1 = x[j + 1];
                const y0 = y[j];
                const y1 = y[j + 1];
                const interpolated = y0 + (y1 - y0) * (xi - x0) / (x1 - x0);
                result.push(interpolated);
            }
        }
        return result;
    };
    
    return interpolate(data.map(p => p.time), data.map(p => p.value), t_grid);
}

// Function to compute rolling window (like Julia implementation)
function computeRollingWindow(data, windowSize) {
    const result = [...data];
    for (let t = windowSize; t < data.length; t++) {
        result[t] = data[t] - data[t - windowSize];
    }
    return result;
}

// Function to compute quantiles (like Julia implementation)
function computeQuantiles(data, qLower = 0.05, qUpper = 0.95) {
    const sortedData = [...data].sort((a, b) => a - b);
    const lowerIndex = Math.floor(qLower * sortedData.length);
    const upperIndex = Math.floor(qUpper * sortedData.length);
    
    const lower = sortedData[lowerIndex];
    const upper = sortedData[upperIndex];
    
    return [lower, upper];
}

// Function to update charts with real data using Julia-style processing
function updateChartsWithRealData() {
    if (!realData || !showRealData) {
        // Remove real data from all charts
        if (typeof charts !== 'undefined') {
            Object.values(charts).forEach(chart => {
                if (chart && chart.data.datasets.length > 2) {
                    chart.data.datasets = chart.data.datasets.slice(0, 2); // Keep only deterministic and stochastic
                    
                    // Recalculate y-axis range after removing real data
                    const compartment = getRealDataCompartment();
                    if (compartment) {
                        recalculateYAxisRange(chart, compartment);
                    }
                    
                    chart.update('none');
                }
            });
        }
        return;
    }
    
    // Get simulation parameters and results
    if (typeof params === 'undefined' || typeof results === 'undefined' || results.length === 0) {
        console.warn('Simulation data not available for Julia-style processing');
        return;
    }
    
    // Create common time grid (like Julia implementation)
    const tGrid = Array.from({length: Math.floor(params.tmax) + 1}, (_, i) => i);
    
    // Interpolate real data to the grid (like Julia implementation)
    const realDataInterpolated = interpolateToGrid(realData.data, tGrid);
    
    // Add real data to the appropriate chart with Julia-style processing
    if (typeof charts !== 'undefined' && charts[realDataCompartment]) {
        const chart = charts[realDataCompartment];
        
        // Remove existing real data dataset if it exists
        if (chart.data.datasets.length > 2) {
            chart.data.datasets = chart.data.datasets.slice(0, 2);
        }
        
        // Add real data dataset (interpolated to grid)
        chart.data.datasets.push({
            label: 'Real World Data (Interpolated)',
            data: tGrid.map((t, i) => ({ x: t, y: realDataInterpolated[i] })),
            borderColor: '#000000',
            backgroundColor: '#00000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#000000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
        
        // Recalculate y-axis range to include real data peak
        recalculateYAxisRange(chart, realDataCompartment);
        
        chart.update('none');
    }
    
    // If we have multiple stochastic runs, we could also compute uncertainty bands
    // This would require access to all simulation results, not just the current run
    console.log('Julia-style real data processing applied');
}



// Function to initialize real data handlers
function initRealDataHandlers() {
    // Upload button
    const uploadBtn = document.getElementById('uploadBtn');
    const fileInput = document.getElementById('realDataUpload');
    
    if (uploadBtn && fileInput) {
        uploadBtn.addEventListener('click', () => {
            fileInput.click();
        });
        
        fileInput.addEventListener('change', uploadRealData);
    }
    
    // Toggle button
    const toggleBtn = document.getElementById('toggleRealDataBtn');
    if (toggleBtn) {
        toggleBtn.addEventListener('click', toggleRealData);
    }
    
    // Clear button
    const clearBtn = document.getElementById('clearRealDataBtn');
    if (clearBtn) {
        clearBtn.addEventListener('click', clearRealData);
    }
    
    // Update initial conditions button
    const updateInitialsBtn = document.getElementById('updateInitialsBtn');
    if (updateInitialsBtn) {
        updateInitialsBtn.addEventListener('click', updateInitialConditionsFromRealData);
    }
    
    // Real data population input handlers
    const populationInput = document.getElementById('realDataPopulation');
    const populationNumberInput = document.getElementById('realDataPopulationNumber');
    const populationValue = document.getElementById('realDataPopulationValue');
    
    if (populationInput && populationNumberInput && populationValue) {
        // Update population value display
        function updateRealDataPopulation(value) {
            value = Math.max(1000, Math.min(1000000, value));
            populationInput.value = value;
            populationNumberInput.value = value;
            populationValue.textContent = value;
        }
        
        // Range slider
        populationInput.addEventListener('input', (e) => {
            updateRealDataPopulation(parseInt(e.target.value));
        });
        
        // Number input
        populationNumberInput.addEventListener('input', (e) => {
            let value = parseInt(e.target.value);
            if (!isNaN(value)) {
                updateRealDataPopulation(value);
            }
        });
    }
}

// Function to process real data with Julia-style analysis
function processRealDataWithJuliaAnalysis() {
    if (!realData || !showRealData || typeof results === 'undefined' || results.length === 0) {
        return;
    }
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    // 🔧 FIX: Use dedicated real data population slider, NOT simulation N
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    // Create common time grid (like Julia implementation)
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Process based on data format and compartment type
    if (realData.dataFormat === 'real_world') {
        // Handle real-world data like Julia implementation
        processRealWorldDataJuliaStyle(tGrid, N);
    } else {
        // Handle simple time series data
        const realDataInterpolated = interpolateToGrid(realData.data, tGrid);
        
        if (realDataCompartment === 'D') {
            // For deaths, compute rolling window like Julia implementation
            const window = 14; // 14-day rolling window like Julia
            const realDataRolling = computeRollingWindow(realDataInterpolated, window);
            
            // Update deaths chart with rolling window data
            if (charts.D) {
                updateChartWithJuliaData(charts.D, tGrid, realDataRolling, 'Real Deaths (14-day Rolling)');
            }
        } else {
            // For other compartments, use interpolated data directly
            if (charts[realDataCompartment]) {
                updateChartWithJuliaData(charts[realDataCompartment], tGrid, realDataInterpolated, 'Real World Data (Interpolated)');
            }
        }
    }
    
    // 🆕 JULIA-STYLE: Also update cumulative deaths if available 
    updateChartsWithCumulativeDeaths();
    
    // 🆕 JULIA-STYLE: Also update hospitalization data if available 
    updateChartsWithHospitalizationData();
    
    console.log('Julia-style analysis completed for', realDataCompartment, 'compartment');
}

// Function to update initial conditions based on real data
function updateInitialConditionsFromRealData() {
    if (!realData || !realData.data || realData.data.length === 0) {
        console.warn('No real data available for initial condition update');
        return;
    }
    
    // Get the real data population
    const realDataPopulation = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    // Get the first data point
    const firstDataPoint = realData.data[0];
    if (!firstDataPoint) {
        console.warn('No data points available');
        return;
    }
    
    // Use the actual initial value from the CSV file
    const initialValue = firstDataPoint.value;
    const initialInfected = initialValue / realDataPopulation;
    const initialSusceptible = 1.0 - initialInfected;
    const initialHospitalized = 0;
    const initialRecovered = 0;
    const initialDeaths = 0;
    
    // Update all input fields directly
    const inputs = {
        's0': initialSusceptible,
        'i0': initialInfected,
        'h0': initialHospitalized,
        'r0': initialRecovered,
        'd0': initialDeaths
    };
    
    // Update each parameter's UI elements
    Object.entries(inputs).forEach(([id, value]) => {
        const slider = document.getElementById(id);
        const numberInput = document.getElementById(`${id}-number`);
        const valueSpan = document.getElementById(`${id}-value`);
        
        if (slider) {
            slider.value = value.toFixed(9);
            // Trigger the input event to update the system
            slider.dispatchEvent(new Event('input', { bubbles: true }));
        }
        if (numberInput) {
            numberInput.value = value.toFixed(9);
            // Trigger the input event to update the system
            numberInput.dispatchEvent(new Event('input', { bubbles: true }));
        }
        if (valueSpan) {
            valueSpan.textContent = value.toFixed(9);
        }
    });
    
    // Update the params object
    if (typeof params !== 'undefined') {
        params.s0 = initialSusceptible;
        params.i0 = initialInfected;
        params.h0 = initialHospitalized;
        params.r0 = initialRecovered;
        params.d0 = initialDeaths;
    }
    
    console.log('Updated initial conditions from real data:', {
        s0: initialSusceptible,
        i0: initialInfected,
        h0: initialHospitalized,
        r0: initialRecovered,
        d0: initialDeaths
    });
}

// Function to process real-world data exactly like Julia implementation
function processRealWorldDataJuliaStyle(tGrid, N) {
    const startDateReal = new Date('2020-03-25'); // Default start date like Julia
    
    if (realDataCompartment === 'I') {
        // Process cases data like Julia implementation
        processCasesDataJuliaStyle(tGrid, N, startDateReal);
    } else if (realDataCompartment === 'D') {
        // Process deaths data like Julia implementation
        processDeathsDataJuliaStyle(tGrid, N, startDateReal);
    }
}

// Function to process cases data like Julia implementation  
function processCasesDataJuliaStyle(tGrid, N, startDateReal) {
    // 🔧 CHECK: If data was already converted to active cases during upload, use it directly
    if (realData.data[0].activeCases !== undefined) {
        console.log(`📊 Using pre-calculated active cases from upload (avoiding double processing)`);
        
        // Data was already converted during upload - use it directly
        const activeCasesCount = realData.data.map(point => point.activeCases);
        const dates = realData.data.map(point => new Date(point.originalDate));
        
        // Convert dates to days from start EXACTLY like Julia implementation
        const firstDate = dates[0]; 
        const carsonDays = dates.map(d => Math.floor((d - firstDate) / (1000 * 60 * 60 * 24)));
        
        // Interpolate active cases to simulation time grid (no need to recalculate)
        const activeCasesInterpolated = interpolateToGrid(
            carsonDays.map((day, i) => ({ time: day, value: activeCasesCount[i] })),
            tGrid
        );
        
        // Convert to proportions AFTER interpolation like Julia implementation
        const populationFactor = 1.0 / N;
        const activeCasesProportions = activeCasesInterpolated.map(val => val * populationFactor);
        
        // Update chart
        if (charts.I) {
            updateChartWithJuliaData(charts.I, tGrid, activeCasesProportions, 'Real Active Cases (14-day Recovery)');
        }
        return;
    }
    
    // FALLBACK: If data wasn't pre-processed, calculate active cases here  
    console.log(`📊 Calculating active cases from raw cumulative data`);
    const recoveryDays = 14; // Longer recovery period for rural areas (14 days)
    
    // Get cumulative cases from the data
    const cumulativeCases = realData.data.map(point => point.value);
    const dates = realData.data.map(point => new Date(point.originalDate));
    
    // Calculate active cases EXACTLY like Julia implementation
    // cumulative_shifted = [zeros(recovery_days); cumulative_cases_from_start[1:end-recovery_days]]
    const cumulativeShifted = [];
    for (let i = 0; i < recoveryDays; i++) {
        cumulativeShifted.push(0);
    }
    for (let i = 0; i < cumulativeCases.length - recoveryDays; i++) {
        cumulativeShifted.push(cumulativeCases[i]);
    }
    
    // active_cases_count = cumulative_cases_from_start - cumulative_shifted
    const activeCasesCount = cumulativeCases.map((val, i) => {
        const shifted = cumulativeShifted[i] || 0;
        return Math.max(val - shifted, 0); // Ensure no negative active cases
    });
    
    // Convert dates to days from start EXACTLY like Julia implementation
    // carson_days = [Dates.value(Day(d - dates_from_start[1])) for d in dates_from_start]
    const firstDate = dates[0];
    const carsonDays = dates.map(d => Math.floor((d - firstDate) / (1000 * 60 * 60 * 24)));
    
    // Interpolate active cases to simulation time grid
    const activeCasesInterpolated = interpolateToGrid(
        carsonDays.map((day, i) => ({ time: day, value: activeCasesCount[i] })),
        tGrid
    );
    
    // Convert to proportions AFTER interpolation like Julia implementation
    const populationFactor = 1.0 / N;
    const activeCasesProportions = activeCasesInterpolated.map(val => val * populationFactor);
    
    // Update chart
    if (charts.I) {
        updateChartWithJuliaData(charts.I, tGrid, activeCasesProportions, 'Real Active Cases (14-day Recovery)');
    }
}

// Function to process deaths data like Julia implementation
function processDeathsDataJuliaStyle(tGrid, N, startDateReal) {
    // For deaths, use the moving_avg_7day if available, otherwise use daily_deaths
    const dates = realData.data.map(point => new Date(point.originalDate));
    
    // Check if we have moving_avg_7day column (like Julia implementation)
    let deathValues;
    let columnType;
    
    if (realData.data[0].moving_avg_7day !== undefined) {
        // Use moving_avg_7day column like Julia implementation
        deathValues = realData.data.map(point => point.moving_avg_7day);
        columnType = 'moving_avg_7day';
    } else if (realData.data[0].daily_deaths !== undefined) {
        // Use daily_deaths column
        deathValues = realData.data.map(point => point.daily_deaths);
        columnType = 'daily_deaths';
    } else {
        // Fallback to value column
        deathValues = realData.data.map(point => point.value);
        columnType = 'value';
    }
    
    // Convert dates to days from start EXACTLY like Julia implementation
    // daily_death_days = [Dates.value(Day(d - start_date_real)) for d in daily_death_dates]
    const deathDays = dates.map(d => Math.floor((d - startDateReal) / (1000 * 60 * 60 * 24)));
    
    // Interpolate deaths to simulation time grid
    const deathsInterpolated = interpolateToGrid(
        deathDays.map((day, i) => ({ time: day, value: deathValues[i] })),
        tGrid
    );
    
    // Convert to proportions AFTER interpolation like Julia implementation
    const populationFactor = 1.0 / N;
    const deathsProportions = deathsInterpolated.map(val => val * populationFactor);
    
    // Update chart with appropriate label
    if (charts.D) {
        const label = columnType === 'moving_avg_7day' ? 'Real Deaths (7-day Moving Avg)' : 
                     columnType === 'daily_deaths' ? 'Real Daily Deaths' : 'Real Deaths';
        updateChartWithJuliaData(charts.D, tGrid, deathsProportions, label);
    }
}

// Function to update chart with Julia-style data
function updateChartWithJuliaData(chart, tGrid, data, label) {
    // Remove existing real data dataset if it exists
    if (chart.data.datasets.length > 2) {
        chart.data.datasets = chart.data.datasets.slice(0, 2);
    }
    
            // Add real data dataset (interpolated to grid)
        chart.data.datasets.push({
            label: label,
            data: tGrid.map((t, i) => ({ x: t, y: data[i] })),
            borderColor: '#000000',
            backgroundColor: '#00000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#000000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
    
    // Recalculate y-axis range to include real data peak
    const compartment = getRealDataCompartment();
    if (compartment) {
        recalculateYAxisRange(chart, compartment);
    }
    
    chart.update('none');
}

// Helper function to get the current real data compartment
function getRealDataCompartment() {
    return realDataCompartment;
}

// Function to recalculate y-axis range after real data is added
function recalculateYAxisRange(chart, compartment) {
    if (!chart || !chart.data.datasets) return;
    
    // Get all values from all datasets (deterministic, stochastic, real data)
    const allValues = [];
    
    chart.data.datasets.forEach(dataset => {
        if (dataset.data && Array.isArray(dataset.data)) {
            const values = dataset.data.map(point => point.y).filter(y => !isNaN(y) && y !== null);
            allValues.push(...values);
        }
    });
    
    if (allValues.length === 0) return;
    
    const maxValue = Math.max(...allValues);
    const minValue = Math.min(...allValues);
    const range = maxValue - minValue;
    
    let finalMax, finalMin;
    
    // Use the same logic as the main calculateDynamicYRange function
    if (compartment === 'I') {
        // Special handling for infected - use max * 1.1
        finalMax = maxValue * 1.1;
        finalMin = Math.max(0, minValue * 1.1);
    } else if (compartment === 'H') {
        // Special handling for hospitalization - minimum padding
        const padding = Math.max(range * 0.1, 0.000001);
        finalMax = Math.min(maxValue + padding, maxValue * 1.1);
        finalMin = Math.max(0, minValue - padding);
    } else {
        // Default handling for S, R, D
        const padding = Math.max(range * 0.1, 0.01);
        finalMax = Math.min(maxValue + padding, maxValue * 1.1);
        finalMin = Math.max(0, minValue - padding);
    }
    
    // Update the chart's y-axis scale
    if (chart.options && chart.options.scales && chart.options.scales.y) {
        chart.options.scales.y.min = finalMin;
        chart.options.scales.y.max = finalMax;
        
        console.log(`🎯 Recalculated Y-axis for ${compartment}: min=${finalMin.toFixed(6)}, max=${finalMax.toFixed(6)} (real data peak: ${maxValue.toFixed(6)})`);
    }
}

// Function to store cumulative deaths data (Julia-style) 
function storeCumulativeDeathsData(processedDataWithDates, population) {
    // Extract cumulative deaths data
    const cumulativeDeaths = processedDataWithDates.map(point => ({
        time: point.time,
        value: point.deaths, // Cumulative deaths count
        proportion: point.deaths / population, // Convert to proportion
        originalDate: point.originalDate
    }));
    
    // Store globally for chart rendering
    cumulativeDeathsData = {
        data: cumulativeDeaths,
        compartment: 'D',
        fileName: 'Cumulative Deaths (from combined CSV)',
        dataFormat: 'cumulative_deaths',
        totalPoints: cumulativeDeaths.length,
        valueRange: {
            min: Math.min(...cumulativeDeaths.map(p => p.value)),
            max: Math.max(...cumulativeDeaths.map(p => p.value))
        }
    };
    
    const maxDeaths = Math.max(...cumulativeDeaths.map(p => p.value));
    console.log(`✅ Peak cumulative deaths: ${maxDeaths.toFixed(0)} (${(maxDeaths/population).toFixed(6)} proportion)`);
}

// Function to store hospitalization data (Julia-style)
function storeHospitalizationData(processedDataWithDates, population) {
    // Store hospitalization data (7-day average of active hospitalizations)
    const hospitalization = processedDataWithDates.map(point => ({
        time: point.time,
        value: point.value, // Hospitalization count (already 7-day average)
        proportion: point.value / population, // Convert to proportion
        originalDate: point.originalDate
    }));
    
    // Store globally for chart rendering
    hospitalizationData = {
        data: hospitalization,
        compartment: 'H',
        fileName: 'Hospitalization Data (7-day Average)',
        dataFormat: 'hospitalization',
        totalPoints: hospitalization.length,
        valueRange: {
            min: Math.min(...hospitalization.map(p => p.value)),
            max: Math.max(...hospitalization.map(p => p.value))
        }
    };
    
    const maxHosp = Math.max(...hospitalization.map(p => p.value));
    console.log(`✅ Peak hospitalization (stored): ${maxHosp.toFixed(1)} (${(maxHosp/population).toFixed(6)} proportion)`);
}

// Function to update real-world data display based on current compartment view (Smart Display)
function updateRealDataForCompartment(currentCompartment) {
    if (!showRealData) {
        return; // Real data is hidden, do nothing
    }
    
    console.log(`🎯 Smart display: Updating real-world data for compartment ${currentCompartment}`);
    
    // Clear real-world data from ALL charts first
    clearRealDataFromAllCharts();
    
    // Show real-world data only for the current compartment
    if (currentCompartment === 'I' && realData && realDataCompartment === 'I') {
        // Show active cases data on Infected chart
        console.log(`📊 Showing active cases real-world data on Infected chart`);
        displayRealDataOnSpecificChart('I');
    } else if (currentCompartment === 'D' && cumulativeDeathsData) {
        // Show cumulative deaths data on Deaths chart
        console.log(`📊 Showing cumulative deaths real-world data on Deaths chart`);
        displayCumulativeDeathsOnChart();
    } else if (currentCompartment === 'H' && hospitalizationData) {
        // Show hospitalization data on Hospitalized chart
        console.log(`📊 Showing hospitalization real-world data on Hospitalized chart`);
        displayHospitalizationOnChart();
    }
    // For other compartments (R, S), show nothing unless we have specific data
}

// Function to clear real-world data from all charts
function clearRealDataFromAllCharts() {
    if (typeof charts !== 'undefined') {
        Object.values(charts).forEach(chart => {
            if (chart && chart.data.datasets.length > 2) {
                chart.data.datasets = chart.data.datasets.slice(0, 2); // Keep only deterministic and stochastic
                chart.update('none');
            }
        });
    }
}

// Function to display real-world data on a specific chart
function displayRealDataOnSpecificChart(compartment) {
    if (!realData || !showRealData) return;
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    // Create common time grid
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Process and display based on compartment
    if (compartment === 'I') {
        // Process cases data
        processRealWorldDataJuliaStyle(tGrid, N);
    }
}

// Function to display cumulative deaths on Deaths chart
function displayCumulativeDeathsOnChart() {
    if (!cumulativeDeathsData || !showRealData) return;
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    // Create common time grid
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Interpolate cumulative deaths to simulation time grid
    const cumulativeDeathsInterpolated = interpolateToGrid(cumulativeDeathsData.data, tGrid);
    
    // Convert to proportions
    const cumulativeDeathsProportions = cumulativeDeathsInterpolated.map(val => val / N);
    
    // Update deaths chart with cumulative deaths
    if (charts.D) {
        charts.D.data.datasets.push({
            label: 'Real Cumulative Deaths',
            data: tGrid.map((t, i) => ({ x: t, y: cumulativeDeathsProportions[i] })),
            borderColor: '#8B0000', // Dark red for cumulative deaths
            backgroundColor: '#8B000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#8B0000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
        
        // Recalculate y-axis range to include cumulative deaths
        recalculateYAxisRange(charts.D, 'D');
        
        charts.D.update('none');
        console.log(`📊 Added cumulative deaths to Deaths chart (smart display)`);
    }
}

// Function to display hospitalization data on Hospitalized chart
function displayHospitalizationOnChart() {
    if (!hospitalizationData || !showRealData) return;
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    console.log(`🔍 Processing hospitalization data: ${hospitalizationData.data.length} points, N=${N}`);
    
    // Create common time grid
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Interpolate hospitalization data to simulation time grid (no extrapolation like Julia)
    const hospitalizationInterpolated = interpolateToGrid(hospitalizationData.data, tGrid, false);
    
    // Convert to proportions (like Julia implementation)
    const hospitalizationProportions = hospitalizationInterpolated.map(val => val / N);
    
    // Filter out NaN values and create valid data points only (like Julia implementation)
    const validDataPoints = [];
    for (let i = 0; i < tGrid.length; i++) {
        if (!isNaN(hospitalizationProportions[i])) {
            validDataPoints.push({ x: tGrid[i], y: hospitalizationProportions[i] });
        }
    }
    
    // Debug: Show raw counts and proportions
    const maxCount = Math.max(...hospitalizationInterpolated.filter(v => !isNaN(v) && v > 0));
    const maxProportion = Math.max(...hospitalizationProportions.filter(v => !isNaN(v) && v > 0));
    console.log(`🏥 Hospitalization data: Max count=${maxCount.toFixed(1)}, Max proportion=${maxProportion.toFixed(6)}`);
    
    // Update hospitalized chart with hospitalization data
    if (charts.H) {
        charts.H.data.datasets.push({
            label: 'Real Hospitalization Data (7-day Avg)',
            data: tGrid.map((t, i) => ({ x: t, y: hospitalizationProportions[i] })),
            borderColor: '#000000', // Black solid line
            backgroundColor: '#00000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#000000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
        
        // 🔧 DYNAMIC: Calculate Y-axis based on ALL data sources (real-world + simulation)
        recalculateYAxisRange(charts.H, 'H');
        
        charts.H.update('none');
        console.log(`📊 Added hospitalization data to Hospitalized chart (smart display)`);
    }
}

// Function to update charts with cumulative deaths data (Julia-style)
function updateChartsWithCumulativeDeaths() {
    if (!cumulativeDeathsData || !showRealData) {
        return;
    }
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    // Create common time grid
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Interpolate cumulative deaths to simulation time grid
    const cumulativeDeathsInterpolated = interpolateToGrid(cumulativeDeathsData.data, tGrid);
    
    // Convert to proportions
    const cumulativeDeathsProportions = cumulativeDeathsInterpolated.map(val => val / N);
    
    // Update deaths chart with cumulative deaths
    if (charts.D) {
        // Remove existing real data dataset if it exists
        if (charts.D.data.datasets.length > 2) {
            charts.D.data.datasets = charts.D.data.datasets.slice(0, 2);
        }
        
        // Add cumulative deaths dataset
        charts.D.data.datasets.push({
            label: 'Real Cumulative Deaths',
            data: tGrid.map((t, i) => ({ x: t, y: cumulativeDeathsProportions[i] })),
            borderColor: '#8B0000', // Dark red for cumulative deaths
            backgroundColor: '#8B000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#8B0000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
        
        // Recalculate y-axis range to include cumulative deaths
        recalculateYAxisRange(charts.D, 'D');
        
        charts.D.update('none');
        console.log(`📊 Added cumulative deaths to Deaths chart (Julia-style)`);
    }
}

// Function to ensure simulation results are plotted on the target compartment chart
function ensureSimulationResultsForCompartment(compartment) {
    // Check if we have access to simulation results
    if (typeof results === 'undefined' || typeof detResult === 'undefined' || 
        !results || results.length === 0 || !detResult) {
        console.log(`⚠️  No simulation results available - need to run simulation for ${compartment} compartment`);
        return;
    }
    
    // Check if the target chart exists and has simulation data
    if (typeof charts !== 'undefined' && charts[compartment]) {
        const chart = charts[compartment];
        
        // Check if chart has at least 2 datasets (deterministic + stochastic)
        if (chart.data.datasets.length < 2) {
            console.log(`🔧 Adding missing simulation results to ${compartment} chart`);
            
            // Ensure we have at least 2 datasets
            while (chart.data.datasets.length < 2) {
                chart.data.datasets.push({
                    label: chart.data.datasets.length === 0 ? 'Deterministic (ODE)' : 'Stochastic (Current Run)',
                    data: [],
                    borderColor: chart.data.datasets.length === 0 ? '#0066CC' : '#00CC66',
                    backgroundColor: chart.data.datasets.length === 0 ? '#0066CC20' : '#00CC6620',
                    borderWidth: 2,
                    pointRadius: 0,
                    fill: false,
                    tension: 0
                });
            }
        }
        
        // Update datasets with current simulation results
        const currentRun = typeof window.currentRun !== 'undefined' ? window.currentRun : 0;
        
        // Deterministic data (dataset 0)
        if (detResult && detResult.T && detResult[`${compartment}_prop`]) {
            chart.data.datasets[0].data = detResult.T.map((t, i) => ({ 
                x: t, 
                y: detResult[`${compartment}_prop`][i] 
            }));
        }
        
        // Stochastic data (dataset 1)  
        if (results[currentRun] && results[currentRun].T && results[currentRun][`${compartment}_prop`]) {
            chart.data.datasets[1].data = results[currentRun].T.map((t, i) => ({ 
                x: t, 
                y: results[currentRun][`${compartment}_prop`][i] 
            }));
        }
        
        // Update chart
        chart.update('none');
        console.log(`✅ Added simulation results to ${compartment} chart`);
    }
}

// Function to update charts with hospitalization data (Julia-style)
function updateChartsWithHospitalizationData() {
    if (!hospitalizationData || !showRealData) {
        return;
    }
    
    // Get simulation parameters
    const tmax = typeof params !== 'undefined' ? params.tmax : 620;
    const N = parseInt(document.getElementById('realDataPopulation')?.value || 56000);
    
    console.log(`🔍 Processing hospitalization data: ${hospitalizationData.data.length} points, N=${N}`);
    
    // Create common time grid
    const tGrid = Array.from({length: Math.floor(tmax) + 1}, (_, i) => i);
    
    // Interpolate hospitalization data to simulation time grid (no extrapolation like Julia)
    const hospitalizationInterpolated = interpolateToGrid(hospitalizationData.data, tGrid, false);
    
    // Convert to proportions (like Julia implementation)
    const hospitalizationProportions = hospitalizationInterpolated.map(val => val / N);
    
    // Filter out NaN values and create valid data points only (like Julia implementation)
    const validDataPoints = [];
    for (let i = 0; i < tGrid.length; i++) {
        if (!isNaN(hospitalizationProportions[i])) {
            validDataPoints.push({ x: tGrid[i], y: hospitalizationProportions[i] });
        }
    }
    
    // Debug: Show raw counts and proportions
    const maxCount = Math.max(...hospitalizationInterpolated.filter(v => !isNaN(v) && v > 0));
    const maxProportion = Math.max(...hospitalizationProportions.filter(v => !isNaN(v) && v > 0));
    console.log(`🏥 Hospitalization data: Max count=${maxCount.toFixed(1)}, Max proportion=${maxProportion.toFixed(6)}`);
    
    // Update hospitalized chart with hospitalization data
    if (charts.H) {
        // Remove existing real data dataset if it exists
        if (charts.H.data.datasets.length > 2) {
            charts.H.data.datasets = charts.H.data.datasets.slice(0, 2);
        }
        
        // Add hospitalization dataset (only where data exists, like Julia)
        charts.H.data.datasets.push({
            label: 'Real Hospitalization Data (7-day Avg)',
            data: validDataPoints, // Only valid data points, no NaN extrapolation
            borderColor: '#000000', // Black solid line
            backgroundColor: '#00000020',
            borderWidth: 3,
            pointRadius: 2,
            pointBackgroundColor: '#000000',
            pointBorderColor: '#fff',
            pointBorderWidth: 1,
            fill: false,
            tension: 0,
            borderDash: []
        });
        
        // 🔧 DYNAMIC: Calculate Y-axis based on ALL data sources (real-world + simulation)
        recalculateYAxisRange(charts.H, 'H');
        
        charts.H.update('none');
        console.log(`📊 Added hospitalization data to Hospitalized chart (Julia-style)`);
    }
}

// Function to compute uncertainty bands like Julia implementation
function computeUncertaintyBands(allResults, compartment, tGrid) {
    if (!allResults || allResults.length === 0) return null;
    
    // Interpolate all simulations to common grid
    const allInterpolated = [];
    for (let i = 0; i < allResults.length; i++) {
        const result = allResults[i];
        if (result && result.T && result[`${compartment}_prop`]) {
            const interpolated = interpolateToGrid(
                result.T.map((t, j) => ({ time: t, value: result[`${compartment}_prop`][j] })),
                tGrid
            );
            allInterpolated.push(interpolated);
        }
    }
    
    if (allInterpolated.length === 0) return null;
    
    // Compute statistics like Julia implementation
    const mean = [];
    const lower = [];
    const upper = [];
    
    for (let t = 0; t < tGrid.length; t++) {
        const valuesAtTime = allInterpolated.map(sim => sim[t]).filter(v => !isNaN(v));
        if (valuesAtTime.length > 0) {
            mean.push(valuesAtTime.reduce((a, b) => a + b, 0) / valuesAtTime.length);
            const [lowerQuantile, upperQuantile] = computeQuantiles(valuesAtTime, 0.05, 0.95);
            lower.push(lowerQuantile);
            upper.push(upperQuantile);
        } else {
            mean.push(0);
            lower.push(0);
            upper.push(0);
        }
    }
    
    return { mean, lower, upper };
}

// Export functions for use in main files
window.realDataHandler = {
    initRealDataHandlers,
    updateChartsWithRealData,
    processRealDataWithJuliaAnalysis,
    computeUncertaintyBands,
    interpolateToGrid,
    computeRollingWindow,
    computeQuantiles,
    recalculateYAxisRange,
    getRealDataCompartment,
    updateChartsWithCumulativeDeaths,
    updateChartsWithHospitalizationData,
    storeCumulativeDeathsData,
    storeHospitalizationData,
    ensureSimulationResultsForCompartment,
    updateRealDataForCompartment,
    clearRealDataFromAllCharts,
    displayRealDataOnSpecificChart,
    displayCumulativeDeathsOnChart,
    displayHospitalizationOnChart,
    realData,
    showRealData,
    realDataCompartment,
    cumulativeDeathsData,
    hospitalizationData
}; 