using CSV, DataFrames, Dates

# Read the data files
active_cases_df = CSV.read("carson_city_active_cases.csv", DataFrame)
daily_deaths_df = CSV.read("carson_city_daily_deaths.csv", DataFrame)

# Convert dates to Date objects
active_cases_df.date = Date.(active_cases_df.date)
daily_deaths_df.date = Date.(daily_deaths_df.date)

# Rename columns to avoid conflicts during join
rename!(daily_deaths_df, :cumulative_deaths => :cumulative_deaths_daily)

# Merge the dataframes on date, keeping only the essential columns
combined_df = leftjoin(active_cases_df, daily_deaths_df, on = :date)

# Select only the three columns needed for P(ID) calculation
result_df = DataFrame(
    date = combined_df.date,
    active_cases = combined_df.active_cases,
    daily_deaths = combined_df.daily_deaths
)

# Save the simplified data
CSV.write("carson_city_pid_simple.csv", result_df)

println("Simplified data saved to carson_city_pid_simple.csv")
println("Columns:")
println("- date: Date")
println("- active_cases: Active cases on that date")
println("- daily_deaths: Deaths on that specific date")

# Calculate some example P(ID) values for demonstration
println("\nExample P(ID) calculations:")
println("Using active cases vs daily deaths:")

# Find dates with reasonable case and death counts (avoiding early zeros)
count = 0
for i in 1:nrow(result_df)
    if result_df.active_cases[i] > 10 && result_df.daily_deaths[i] > 0 && count < 3
        pid_estimate = result_df.daily_deaths[i] / result_df.active_cases[i]
        println("Date: $(result_df.date[i])")
        println("  Active cases: $(result_df.active_cases[i])")
        println("  Daily deaths: $(result_df.daily_deaths[i])")
        println("  P(ID) estimate: $(round(pid_estimate, digits=4))")
        println()
        count += 1
    end
end
