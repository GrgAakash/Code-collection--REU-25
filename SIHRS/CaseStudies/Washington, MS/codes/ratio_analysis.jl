using CSV
using DataFrames
using Dates

# Read the data
active_cases_df = CSV.read("washington_mississippi_active_cases.csv", DataFrame)
hospitalization_df = CSV.read("hospitalization_MS_filtered.csv", DataFrame)

# Parse dates from hospitalization data (handle 2-digit years)
function parse_date_manual(date_str)
    parts = split(date_str, '/')
    month = parse(Int, parts[1])
    day = parse(Int, parts[2])
    year_2digit = parse(Int, parts[3])
    year = year_2digit + 2000
    return Date(year, month, day)
end

hospitalization_df.date = parse_date_manual.(hospitalization_df.collection_week)

# Convert active cases dates to Date type
active_cases_df.date = Date.(active_cases_df.date)

println("Finding exact date matches and calculating ratios...")

# Results storage
results = DataFrame(
    date = Date[],
    active_cases = Float64[],
    hospitalized = Float64[],
    ratio = Float64[]
)

# For each hospitalization date, find the exact matching active cases date
for hosp_row in eachrow(hospitalization_df)
    hosp_date = hosp_row.date
    hospitalized = hosp_row.total_adult_and_pediatric_covid_patients
    
    # Skip zero hospitalization dates to avoid infinite ratios
    if hospitalized == 0.0
        continue
    end
    
    # Find exact matching date in active cases
    active_match = filter(row -> row.date == hosp_date, active_cases_df)
    
    if nrow(active_match) > 0
        active_cases = active_match[1, :active_cases]
        ratio = active_cases / hospitalized
        
        push!(results, (hosp_date, active_cases, hospitalized, ratio))
    end
end

println("DATES WHERE RATIO IS AROUND 4 (between 3.5 and 4.5):")
println("="^60)

# Filter for ratios around 4
target_ratios = filter(row -> 3.5 <= row.ratio <= 4.5, results)
sort!(target_ratios, :ratio)

if nrow(target_ratios) > 0
    for row in eachrow(target_ratios)
        println("$(row.date): $(row.active_cases) active → $(row.hospitalized) hospitalized → Ratio: $(round(row.ratio, digits=2))")
    end
    println("\nTotal matches found: $(nrow(target_ratios))")
else
    println("No dates found with ratio around 4")
end

println("\nTOP 5 CLOSEST TO RATIO = 4.0:")
println("="^40)

# Sort by distance from 4.0
results_sorted = copy(results)
results_sorted.ratio_diff = abs.(results_sorted.ratio .- 4.0)
sort!(results_sorted, :ratio_diff)

for i in 1:min(5, nrow(results_sorted))
    row = results_sorted[i, :]
    println("$(row.date): Ratio $(round(row.ratio, digits=2)) ($(row.active_cases) active, $(row.hospitalized) hospitalized)")
end

# Save results
CSV.write("ratio_analysis_results.csv", results)
CSV.write("target_ratio_dates.csv", target_ratios)

println("\nResults saved to:")
println("- ratio_analysis_results.csv (all date matches)")
println("- target_ratio_dates.csv (ratios around 4)")
