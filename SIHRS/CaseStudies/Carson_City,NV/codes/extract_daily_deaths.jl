using CSV
using DataFrames
using Dates
using Printf
using Statistics

# Load the data
data_table = CSV.read("carson_city_combined.csv", DataFrame)

# Convert date column to Date type
data_table.date = Date.(data_table.date)

# Sort by date to ensure chronological order
sort!(data_table, :date)

println("Data loaded successfully!")
println("Date range: $(data_table.date[1]) to $(data_table.date[end])")
println("Total rows: $(nrow(data_table))")

# Calculate daily deaths using simple difference
# Daily Deaths = Cumulative(n) - Cumulative(n-1)
daily_deaths = zeros(Int, nrow(data_table))

# First day: daily deaths = cumulative deaths (no previous day)
daily_deaths[1] = data_table.deaths[1]

# Rest of the days: difference from previous day
for i in 2:nrow(data_table)
    daily_deaths[i] = data_table.deaths[i] - data_table.deaths[i-1]
end

# Create results DataFrame
results_df = DataFrame(
    date = data_table.date,
    cumulative_deaths = data_table.deaths,
    daily_deaths = daily_deaths
)

# Display key dates
println("\n=== Key Dates ===")
println("March 25, 2020:")
march25_idx = findfirst(results_df.date .== Date("2020-03-25"))
if !isnothing(march25_idx)
    println("  Cumulative deaths: $(results_df.cumulative_deaths[march25_idx])")
    println("  Daily deaths: $(results_df.daily_deaths[march25_idx])")
end

println("\nAugust 2, 2020:")
august2_idx = findfirst(results_df.date .== Date("2020-08-02"))
if !isnothing(august2_idx)
    println("  Cumulative deaths: $(results_df.cumulative_deaths[august2_idx])")
    println("  Daily deaths: $(results_df.daily_deaths[august2_idx])")
end

# Show a sample of the data
println("\n=== Sample Data (first 30 rows) ===")
println(results_df[1:30, :])

# Show statistics
println("\n=== Statistics ===")
println("Maximum daily deaths: $(maximum(results_df.daily_deaths))")
println("Date of max daily deaths: $(results_df.date[argmax(results_df.daily_deaths)])")
println("Total cumulative deaths: $(maximum(results_df.cumulative_deaths))")
println("Total days with deaths: $(count(x -> x > 0, results_df.daily_deaths))")

# Calculate 7-day moving average for smoothing
println("\n=== 7-Day Moving Average ===")
moving_avg_7 = zeros(Float64, nrow(results_df))
for i in 1:nrow(results_df)
    start_idx = max(1, i-6)
    moving_avg_7[i] = mean(results_df.daily_deaths[start_idx:i])
end

results_df.moving_avg_7day = moving_avg_7

# Show moving average for key dates
if !isnothing(august2_idx)
    println("August 2, 2020:")
    println("  Daily deaths: $(results_df.daily_deaths[august2_idx])")
    println("  7-day moving average: $(round(results_df.moving_avg_7day[august2_idx], digits=2))")
end

# Save results to CSV
CSV.write("carson_city_daily_deaths.csv", results_df)
println("\nResults saved to 'carson_city_daily_deaths.csv'")

# Show the calculation method
println("\n=== Daily Deaths Calculation Method ===")
println("1. Take cumulative deaths: [$(results_df.cumulative_deaths[1]), $(results_df.cumulative_deaths[2]), ..., $(results_df.cumulative_deaths[end])]")
println("2. Calculate: daily_deaths[i] = cumulative_deaths[i] - cumulative_deaths[i-1]")
println("3. For first day: daily_deaths[1] = cumulative_deaths[1]")
println("4. Result: daily deaths represent new deaths each day")
println("\nThis is much simpler and more logical than rolling averages for deaths!") 