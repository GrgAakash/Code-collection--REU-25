using CSV
using DataFrames
using Dates

# Test hospitalization data loading
println("Loading Washington, Mississippi hospitalization data...")
hosp_data_table = CSV.read("hospitalization_MS_filtered.csv", DataFrame)
println("Raw data loaded. First few rows:")
println(first(hosp_data_table, 5))

println("\nData types:")
println(typeof(hosp_data_table.collection_week))
println(typeof(hosp_data_table.total_adult_and_pediatric_covid_patients))

# Check for duplicates
println("\nChecking for duplicates...")
date_counts = combine(groupby(hosp_data_table, :collection_week), nrow => :count)
duplicate_dates = date_counts[date_counts.count .> 1, :]
println("Number of dates with duplicates: ", nrow(duplicate_dates))
if nrow(duplicate_dates) > 0
    println("Sample duplicate dates:")
    println(first(duplicate_dates, 3))
end

# Remove duplicates by keeping only non-zero values
println("\nRemoving duplicates by keeping only non-zero values...")
hosp_data_table_filtered = hosp_data_table[hosp_data_table.total_adult_and_pediatric_covid_patients .> 0, :]
println("After filtering: ", nrow(hosp_data_table_filtered), " rows")

# Convert dates
println("\nConverting dates...")
hosp_data_table_filtered.collection_week = Date.(hosp_data_table_filtered.collection_week, "m/d/yy")
println("Date conversion successful!")

# Fix the year to be 20xx instead of 00xx
for i in 1:length(hosp_data_table_filtered.collection_week)
    d = hosp_data_table_filtered.collection_week[i]
    hosp_data_table_filtered.collection_week[i] = Date(2000 + year(d), month(d), day(d))
end

# Sort by date
sort!(hosp_data_table_filtered, :collection_week)

println("\nDate range: $(hosp_data_table_filtered.collection_week[1]) to $(hosp_data_table_filtered.collection_week[end])")

# Test interpolation
simulation_start_date = Date("2020-03-25")
hosp_days = [Dates.value(Day(d - simulation_start_date)) for d in hosp_data_table_filtered.collection_week]
hospitalization_data = hosp_data_table_filtered.total_adult_and_pediatric_covid_patients

println("\nHospitalization data sample:")
for i in 1:min(10, length(hospitalization_data))
    println("Day $(hosp_days[i]): $(hospitalization_data[i]) patients")
end

println("\nMaximum hospitalization: $(maximum(hospitalization_data))")
println("Date of max: $(hosp_data_table_filtered.collection_week[argmax(hospitalization_data)])") 