using CSV, DataFrames, Dates, Statistics

# Read the combined hospitalization and active cases data
combined_df = CSV.read("washington_ms_hospitalization_active_cases_combined.csv", DataFrame)

# Function to convert hospitalization date format (M/d/yy) to Date
function parse_hosp_date(date_str)
    # Handle the M/d/yy format (e.g., "8/2/20" -> Date(2020, 8, 2))
    parts = split(date_str, "/")
    month = parse(Int, parts[1])
    day = parse(Int, parts[2])
    year = parse(Int, parts[3])
    
    # Convert 2-digit year to 4-digit year
    if year < 50
        year += 2000
    elseif year < 100
        year += 1900
    end
    
    return Date(year, month, day)
end

# Convert dates to Date objects for easier manipulation
combined_df.date_parsed = [parse_hosp_date(date_str) for date_str in combined_df.date]

# Sort by date to ensure proper ordering
sort!(combined_df, :date_parsed)

# Create a new dataframe for P(IH) calculations
pih_data = DataFrame()

# Calculate P(IH) with 14-day lag (using 2 rows back since data is weekly)
# We need at least 3 rows of data to calculate the lag (start from row 3)
for i in 3:nrow(combined_df)
    # Day T (current day)
    day_t = combined_df.date_parsed[i]
    hospitalized_t = combined_df.total_adult_and_pediatric_covid_patients[i]
    
    # Day T-14 (approximately 14 days earlier, 2 rows back in weekly data)
    day_t_minus_14 = combined_df.date_parsed[i-2]
    active_cases_t_minus_14 = combined_df.active_cases[i-2]
    
    # Check that the date difference is reasonable (not more than 15 days)
    date_diff = Dates.value(day_t - day_t_minus_14)
    
    if date_diff <= 15 && date_diff >= 10  # Between 10-15 days is acceptable
        # Calculate P(IH) = hospitalized_T / active_cases_T-14
        if active_cases_t_minus_14 > 0
            pih_value = hospitalized_t / active_cases_t_minus_14
        else
            pih_value = missing  # Handle division by zero
        end
        
        # Add to result dataframe
        push!(pih_data, (
            date_T = combined_df.date[i],
            date_T_minus_14 = combined_df.date[i-2],
            hospitalized_T = hospitalized_t,
            active_cases_T_minus_14 = active_cases_t_minus_14,
            date_difference_days = date_diff,
            P_IH_estimate = pih_value
        ))
    else
        println("Skipping date pair with unrealistic gap: $(combined_df.date[i]) -> $(combined_df.date[i-2]) ($(date_diff) days)")
    end
end

# Save the P(IH) calculation data
CSV.write("washington_ms_pih_lagged_14days.csv", pih_data)

println("P(IH) calculation data saved to washington_ms_pih_lagged_14days.csv")
println("Columns:")
println("- date_T: Current date (day T)")
println("- date_T_minus_14: Date ~14 days earlier (day T-14)")
println("- hospitalized_T: Hospitalized patients on day T")
println("- active_cases_T_minus_14: Active cases on day T-14")
println("- date_difference_days: Actual days between T and T-14")
println("- P_IH_estimate: P(IH) = hospitalized_T / active_cases_T_minus_14")

# Show some example calculations
println("\nExample P(IH) calculations using 14-day lag:")
println("Format: Date T | Hospitalized T | Cases T-14 | P(IH)")

example_count = 0
for i in 1:nrow(pih_data)
    global example_count
    if pih_data.hospitalized_T[i] > 0 && pih_data.active_cases_T_minus_14[i] > 10 && example_count < 10
        println("$(pih_data.date_T[i]) | $(pih_data.hospitalized_T[i]) | $(pih_data.active_cases_T_minus_14[i]) | $(round(pih_data.P_IH_estimate[i], digits=4))")
        example_count += 1
    end
end

# Calculate summary statistics
valid_pih = filter(row -> !ismissing(row.P_IH_estimate) && row.active_cases_T_minus_14 > 0, pih_data)

if nrow(valid_pih) > 0
    println("\n" * "="^60)
    println("P(IH) SUMMARY STATISTICS")
    println("="^60)
    println("Total P(IH) calculations: $(nrow(pih_data))")
    println("Valid P(IH) calculations: $(nrow(valid_pih))")
    println()
    println("Average P(IH): $(round(mean(skipmissing(valid_pih.P_IH_estimate)), digits=6))")
    println("Median P(IH): $(round(median(skipmissing(valid_pih.P_IH_estimate)), digits=6))")
    println("Minimum P(IH): $(round(minimum(skipmissing(valid_pih.P_IH_estimate)), digits=6))")
    println("Maximum P(IH): $(round(maximum(skipmissing(valid_pih.P_IH_estimate)), digits=6))")
    println()
    println("This average P(IH) value represents the probability that an")
    println("infected person becomes hospitalized after a 14-day progression period.")
    println("This can be used as the pIH parameter in your SIHRS model.")
    
    # Show distribution information
    avg_pih = mean(skipmissing(valid_pih.P_IH_estimate))
    println("\nRecommended parameter for SIHRS model:")
    println("pIH = $(round(avg_pih, digits=4))")
    
else
    println("\nNo valid P(IH) calculations could be performed.")
end
