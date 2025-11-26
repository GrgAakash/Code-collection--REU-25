using CSV, DataFrames, Dates, Statistics

# Read the active cases data
active_cases_df = CSV.read("carson_city_active_cases.csv", DataFrame)
active_cases_df.date = Date.(active_cases_df.date)

# Read the hospitalization data
hosp_df = CSV.read("hospitalization_Carson_filtered_new.csv", DataFrame)

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

# Convert hospitalization dates to Date objects
hosp_df.date_parsed = [parse_hosp_date(date_str) for date_str in hosp_df.collection_week]

# Sort both dataframes by date
sort!(active_cases_df, :date)
sort!(hosp_df, :date_parsed)

# Create a new dataframe for the combined data
combined_data = DataFrame()

# For each hospitalization date, find the closest active cases date
for i in 1:nrow(hosp_df)
    hosp_date = hosp_df.date_parsed[i]
    hospitalized_count = hosp_df.total_adult_and_pediatric_covid_patients[i]
    
    # Find the closest active cases date (within ±7 days)
    date_diffs = abs.(Dates.value.(active_cases_df.date .- hosp_date))
    closest_idx = argmin(date_diffs)
    closest_date = active_cases_df.date[closest_idx]
    date_diff = abs(Dates.value(closest_date - hosp_date))
    
    # Only include if dates are within 7 days of each other
    if date_diff <= 7
        active_cases_count = active_cases_df.active_cases[closest_idx]
        
        # Add to combined dataframe
        push!(combined_data, (
            date = hosp_df.collection_week[i],
            date_parsed = hosp_date,
            active_cases = active_cases_count,
            total_hospitalized = hospitalized_count,
            date_difference_days = date_diff
        ))
    end
end

# Sort by date
sort!(combined_data, :date_parsed)

# Save the combined data
CSV.write("carson_city_active_cases_hospitalization.csv", combined_data)

println("Combined data saved to carson_city_active_cases_hospitalization.csv")
println()
println("Data summary:")
println("Total rows: $(nrow(combined_data))")
println("Date range: $(Dates.format(combined_data.date_parsed[1], "yyyy-mm-dd")) to $(Dates.format(combined_data.date_parsed[end], "yyyy-mm-dd"))")
println()
println("Columns:")
println("- date: Hospitalization date (M/d/yy format)")
println("- date_parsed: Date in Date format")
println("- active_cases: Number of active cases on closest date")
println("- total_hospitalized: Total hospitalized patients")
println("- date_difference_days: Days between hospitalization and active cases dates")
println()

# Show some example data
println("Example data:")
println("Date | Active Cases | Hospitalized | Date Diff")
println("-" * "="^50)
for i in 1:min(10, nrow(combined_data))
    println("$(combined_data.date[i]) | $(combined_data.active_cases[i]) | $(combined_data.total_hospitalized[i]) | $(combined_data.date_difference_days[i])")
end

# Calculate some basic statistics
if nrow(combined_data) > 0
    println()
    println("Statistics:")
    println("Average active cases: $(round(mean(combined_data.active_cases), digits=1))")
    println("Average hospitalized: $(round(mean(combined_data.total_hospitalized), digits=1))")
    println("Maximum active cases: $(maximum(combined_data.active_cases))")
    println("Maximum hospitalized: $(maximum(combined_data.total_hospitalized))")
end
