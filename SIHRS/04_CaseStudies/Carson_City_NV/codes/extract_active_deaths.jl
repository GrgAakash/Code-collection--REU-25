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

# Calculate active deaths using rolling average method
# For deaths, we typically use a shorter window since death reporting is more immediate
# Let's try different windows: 7 days, 14 days, and 30 days

recovery_days_options = [7, 14, 30]

for recovery_days in recovery_days_options
    println("\n=== Active Deaths ($(recovery_days)-day window) ===")
    
    # Calculate active deaths
    cumulative_deaths = data_table.deaths
    cumulative_shifted = [zeros(recovery_days); cumulative_deaths[1:end-recovery_days]]
    active_deaths_count = cumulative_deaths - cumulative_shifted
    
    # Ensure no negative active deaths due to data corrections
    active_deaths_count = max.(active_deaths_count, 0)
    
    # Create a DataFrame with the results
    results_df = DataFrame(
        date = data_table.date,
        cumulative_deaths = data_table.deaths,
        active_deaths = active_deaths_count
    )
    
    # Display key dates
    println("March 25, 2020:")
    march25_idx = findfirst(results_df.date .== Date("2020-03-25"))
    if !isnothing(march25_idx)
        println("  Cumulative deaths: $(results_df.cumulative_deaths[march25_idx])")
        println("  Active deaths: $(results_df.active_deaths[march25_idx])")
    end
    
    println("August 2, 2020:")
    august2_idx = findfirst(results_df.date .== Date("2020-08-02"))
    if !isnothing(august2_idx)
        println("  Cumulative deaths: $(results_df.cumulative_deaths[august2_idx])")
        println("  Active deaths: $(results_df.active_deaths[august2_idx])")
    end
    
    # Show statistics
    println("Maximum active deaths: $(maximum(results_df.active_deaths))")
    println("Date of max active deaths: $(results_df.date[argmax(results_df.active_deaths)])")
    println("Total cumulative deaths: $(maximum(results_df.cumulative_deaths))")
    
    # Save results to CSV
    filename = "carson_city_active_deaths_$(recovery_days)day.csv"
    CSV.write(filename, results_df)
    println("Results saved to '$filename'")
end

# Let's also create a comprehensive file with all windows
println("\n=== Comprehensive Active Deaths Analysis ===")

# Calculate for all windows
cumulative_deaths = data_table.deaths
active_deaths_7 = max.(cumulative_deaths - [zeros(7); cumulative_deaths[1:end-7]], 0)
active_deaths_14 = max.(cumulative_deaths - [zeros(14); cumulative_deaths[1:end-14]], 0)
active_deaths_30 = max.(cumulative_deaths - [zeros(30); cumulative_deaths[1:end-30]], 0)

comprehensive_df = DataFrame(
    date = data_table.date,
    cumulative_deaths = data_table.deaths,
    active_deaths_7day = active_deaths_7,
    active_deaths_14day = active_deaths_14,
    active_deaths_30day = active_deaths_30
)

# Show comparison for August 2, 2020
august2_idx = findfirst(comprehensive_df.date .== Date("2020-08-02"))
if !isnothing(august2_idx)
    println("August 2, 2020 Comparison:")
    println("  Cumulative deaths: $(comprehensive_df.cumulative_deaths[august2_idx])")
    println("  Active deaths (7-day): $(comprehensive_df.active_deaths_7day[august2_idx])")
    println("  Active deaths (14-day): $(comprehensive_df.active_deaths_14day[august2_idx])")
    println("  Active deaths (30-day): $(comprehensive_df.active_deaths_30day[august2_idx])")
end

# Save comprehensive results
CSV.write("carson_city_active_deaths_comprehensive.csv", comprehensive_df)
println("Comprehensive results saved to 'carson_city_active_deaths_comprehensive.csv'")

# Show the rolling average calculation method for deaths
println("\n=== Rolling Average Method for Deaths ===")
println("1. Take cumulative deaths: [$(cumulative_deaths[1]), $(cumulative_deaths[2]), ..., $(cumulative_deaths[end])]")
println("2. Shift by recovery_days (7, 14, or 30 days)")
println("3. Calculate: active_deaths = cumulative_deaths - cumulative_shifted")
println("4. Ensure active_deaths >= 0 (no negative values)")
println("5. Result: active deaths represent recent deaths in the specified window")
println("\nNote: For deaths, shorter windows (7-14 days) are typically more appropriate")
println("since death reporting is more immediate than case recovery.") 