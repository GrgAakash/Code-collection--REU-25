using CSV
using DataFrames
using Dates

# Read the data
active_cases_df = CSV.read("carson_city_active_cases.csv", DataFrame)
# Note: Carson City doesn't have hospitalization data in the same format as Washington MS
# We'll need to check what hospitalization data is available for Carson City

println("Checking available data files for Carson City...")

# Let's see what files we have available
try
    # Check if there's a hospitalization file for Carson City
    if isfile("hospitalization_Carson_filtered_new.csv")
        hospitalization_df = CSV.read("hospitalization_Carson_filtered_new.csv", DataFrame)
        println("Found hospitalization data: hospitalization_Carson_filtered_new.csv")
        println("Columns: ", names(hospitalization_df))
        println("First few rows:")
        println(first(hospitalization_df, 3))
        
        # Parse dates from hospitalization data (handle 2-digit years if needed)
        function parse_date_manual(date_str)
            # Handle different date formats
            if occursin("/", string(date_str))
                parts = split(string(date_str), '/')
                month = parse(Int, parts[1])
                day = parse(Int, parts[2])
                year_part = parts[3]
                
                # Handle 2-digit vs 4-digit years
                if length(year_part) == 2
                    year_2digit = parse(Int, year_part)
                    year = year_2digit + 2000
                else
                    year = parse(Int, year_part)
                end
                return Date(year, month, day)
            else
                # Try parsing as standard date format
                return Date(date_str)
            end
        end

        # Convert active cases dates to Date type
        active_cases_df.date = Date.(active_cases_df.date)
        
        # Find the hospitalization column (might have different name)
        hosp_columns = names(hospitalization_df)
        println("\nAvailable columns in hospitalization data:")
        for col in hosp_columns
            println("  - $col")
        end
        
        # Look for columns that might contain hospitalization data
        potential_hosp_cols = filter(col -> 
            occursin("hospital", lowercase(string(col))) || 
            occursin("covid", lowercase(string(col))) ||
            occursin("patient", lowercase(string(col))), hosp_columns)
            
        if !isempty(potential_hosp_cols)
            println("\nPotential hospitalization columns: $potential_hosp_cols")
            
            # Use the first potential column (you might need to adjust this)
            hosp_col = potential_hosp_cols[1]
            println("Using column: $hosp_col")
            
            # Try to parse dates
            if "collection_week" in names(hospitalization_df)
                hospitalization_df.date = parse_date_manual.(hospitalization_df.collection_week)
            elseif "date" in names(hospitalization_df)
                hospitalization_df.date = Date.(hospitalization_df.date)
            else
                println("Warning: Could not find date column in hospitalization data")
                # Use the first column that looks like dates
                date_col = names(hospitalization_df)[1]
                hospitalization_df.date = parse_date_manual.(hospitalization_df[!, date_col])
            end
            
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
                hospitalized = hosp_row[hosp_col]
                
                # Skip zero or missing hospitalization dates to avoid infinite ratios
                if ismissing(hospitalized) || hospitalized == 0.0
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

            println("DATES WHERE RATIO IS AROUND 8.37 (between 7.5 and 9.5):")
            println("="^60)

            # Filter for ratios around 8.37
            target_ratios = filter(row -> 7.5 <= row.ratio <= 9.5, results)
            sort!(target_ratios, :ratio)

            if nrow(target_ratios) > 0
                for row in eachrow(target_ratios)
                    println("$(row.date): $(row.active_cases) active → $(row.hospitalized) hospitalized → Ratio: $(round(row.ratio, digits=2))")
                end
                println("\nTotal matches found: $(nrow(target_ratios))")
            else
                println("No dates found with ratio around 8.37")
            end

            println("\nTOP 5 CLOSEST TO RATIO = 8.37:")
            println("="^40)

            # Sort by distance from 8.37
            if nrow(results) > 0
                results_sorted = copy(results)
                results_sorted.ratio_diff = abs.(results_sorted.ratio .- 8.37)
                sort!(results_sorted, :ratio_diff)

                for i in 1:min(5, nrow(results_sorted))
                    row = results_sorted[i, :]
                    println("$(row.date): Ratio $(round(row.ratio, digits=2)) ($(row.active_cases) active, $(row.hospitalized) hospitalized)")
                end

                # Save results
                CSV.write("carson_ratio_analysis_results.csv", results)
                CSV.write("carson_target_ratio_dates.csv", target_ratios)

                println("\nResults saved to:")
                println("- carson_ratio_analysis_results.csv (all date matches)")
                println("- carson_target_ratio_dates.csv (ratios around 4)")
            else
                println("No matching dates found between active cases and hospitalization data")
            end
        else
            println("Could not identify hospitalization column in the data")
        end
        
    else
        println("Hospitalization data file not found for Carson City")
        println("Available files in current directory:")
        for file in readdir(".")
            if endswith(file, ".csv")
                println("  - $file")
            end
        end
    end
    
catch e
    println("Error processing data: $e")
    println("\nAvailable CSV files in current directory:")
    for file in readdir(".")
        if endswith(file, ".csv")
            println("  - $file")
        end
    end
end
