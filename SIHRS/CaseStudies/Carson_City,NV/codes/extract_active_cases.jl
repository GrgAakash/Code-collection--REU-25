using CSV
using DataFrames
using Dates
using Printf

# Load the data
data_table = CSV.read("carson_city_combined.csv", DataFrame)

# Convert date column to Date type
data_table.date = Date.(data_table.date)

# Sort by date to ensure chronological order
sort!(data_table, :date)

println("Data loaded successfully!")
println("Date range: $(data_table.date[1]) to $(data_table.date[end])")
println("Total rows: $(nrow(data_table))")

# Calculate active cases using rolling average method
# This is the same method used in SIHRS_multiple_simulations_infected.jl

recovery_days = 14  # Recovery period for rural areas (14 days)

# Calculate active cases
cumulative_cases = data_table.cases
cumulative_shifted = [zeros(recovery_days); cumulative_cases[1:end-recovery_days]]
active_cases_count = cumulative_cases - cumulative_shifted

# Ensure no negative active cases due to data corrections
active_cases_count = max.(active_cases_count, 0)

# Create a new DataFrame with the results
results_df = DataFrame(
    date = data_table.date,
    cumulative_cases = data_table.cases,
    cumulative_deaths = data_table.deaths,
    active_cases = active_cases_count
)

# Display some key dates
println("\n=== Key Dates ===")
println("March 25, 2020:")
march25_idx = findfirst(results_df.date .== Date("2020-03-25"))
if !isnothing(march25_idx)
    println("  Cumulative cases: $(results_df.cumulative_cases[march25_idx])")
    println("  Cumulative deaths: $(results_df.cumulative_deaths[march25_idx])")
    println("  Active cases: $(results_df.active_cases[march25_idx])")
end

println("\nAugust 2, 2020:")
august2_idx = findfirst(results_df.date .== Date("2020-08-02"))
if !isnothing(august2_idx)
    println("  Cumulative cases: $(results_df.cumulative_cases[august2_idx])")
    println("  Cumulative deaths: $(results_df.cumulative_deaths[august2_idx])")
    println("  Active cases: $(results_df.active_cases[august2_idx])")
end

# Show a sample of the data
println("\n=== Sample Data (first 20 rows) ===")
println(results_df[1:20, :])

# Show some statistics
println("\n=== Statistics ===")
println("Maximum active cases: $(maximum(results_df.active_cases))")
println("Date of max active cases: $(results_df.date[argmax(results_df.active_cases)])")
println("Maximum cumulative cases: $(maximum(results_df.cumulative_cases))")
println("Maximum cumulative deaths: $(maximum(results_df.cumulative_deaths))")

# Save results to CSV
CSV.write("carson_city_active_cases.csv", results_df)
println("\nResults saved to 'carson_city_active_cases.csv'")

# Show the rolling average calculation method
println("\n=== Rolling Average Method ===")
println("1. Take cumulative cases: [$(cumulative_cases[1]), $(cumulative_cases[2]), ..., $(cumulative_cases[end])]")
println("2. Shift by recovery_days = $recovery_days days")
println("3. Calculate: active_cases = cumulative_cases - cumulative_shifted")
println("4. Ensure active_cases >= 0 (no negative values)")
println("5. Result: active cases represent currently infected people") 