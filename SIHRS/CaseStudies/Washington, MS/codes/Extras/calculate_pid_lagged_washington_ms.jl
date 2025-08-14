using CSV, DataFrames, Dates

# Read the simplified data file for Washington, MS
data_df = CSV.read("washington_ms_pid_simple.csv", DataFrame)

# Convert dates to Date objects
data_df.date = Date.(data_df.date)

# Create a new dataframe for the lagged calculation
result_df = DataFrame()

# We need at least 21 days of data to calculate the lag
# Start from day 21 (index 21) since we need T-20
for i in 21:nrow(data_df)
    # Day T (current day)
    day_t = data_df.date[i]
    active_deaths_t = data_df.daily_deaths[i]
    
    # Day T-20 (20 days earlier)
    day_t_minus_20 = data_df.date[i-20]
    active_cases_t_minus_20 = data_df.active_cases[i-20]
    
    # Calculate P(ID) = active_deaths_T / active_cases_T-20
    if active_cases_t_minus_20 > 0
        pid_value = active_deaths_t / active_cases_t_minus_20
    else
        pid_value = missing  # Handle division by zero
    end
    
    # Add to result dataframe
    push!(result_df, (
        date_T = day_t,
        date_T_minus_20 = day_t_minus_20,
        active_deaths_T = active_deaths_t,
        active_cases_T_minus_20 = active_cases_t_minus_20,
        P_ID_estimate = pid_value
    ))
end

# Save the lagged calculation data
CSV.write("washington_ms_pid_lagged_20days.csv", result_df)

println("Lagged P(ID) data saved to washington_ms_pid_lagged_20days.csv")
println("Columns:")
println("- date_T: Current date (day T)")
println("- date_T_minus_20: Date 20 days earlier (day T-20)")
println("- active_deaths_T: Deaths on day T")
println("- active_cases_T_minus_20: Active cases on day T-20")
println("- P_ID_estimate: P(ID) = active_deaths_T / active_cases_T_minus_20")

# Show some example calculations
println("\nExample P(ID) calculations using 20-day lag:")
println("Format: Date T | Deaths T | Cases T-20 | P(ID)")

example_count = 0
for i in 1:nrow(result_df)
    if result_df.active_deaths_T[i] > 0 && result_df.active_cases_T_minus_20[i] > 10 && example_count < 5
        println("$(result_df.date_T[i]) | $(result_df.active_deaths_T[i]) | $(result_df.active_cases_T_minus_20[i]) | $(round(result_df.P_ID_estimate[i], digits=4))")
        example_count += 1
    end
end
