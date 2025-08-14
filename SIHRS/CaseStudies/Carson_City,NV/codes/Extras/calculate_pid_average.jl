using CSV, DataFrames, Statistics

# Read the lagged P(ID) data
pid_df = CSV.read("carson_city_pid_lagged_20days.csv", DataFrame)

# Calculate the average P(ID) across all rows
# Filter out any missing or invalid values
valid_pid_values = filter(x -> !ismissing(x) && isfinite(x), pid_df.P_ID_estimate)

# Calculate statistics
total_rows = nrow(pid_df)
valid_rows = length(valid_pid_values)
average_pid = mean(valid_pid_values)
median_pid = median(valid_pid_values)
min_pid = minimum(valid_pid_values)
max_pid = maximum(valid_pid_values)

println("P(ID) Statistics from carson_city_pid_lagged_20days.csv:")
println("==================================================")
println("Total rows: $total_rows")
println("Valid P(ID) values: $valid_rows")
println()
println("Average P(ID): $(round(average_pid, digits=6))")
println("Median P(ID): $(round(median_pid, digits=6))")
println("Minimum P(ID): $(round(min_pid, digits=6))")
println("Maximum P(ID): $(round(max_pid, digits=6))")
println()
println("This average P(ID) value can be used as a realistic parameter")
println("for your SIHRS model instead of the current assumption of 0.001")
