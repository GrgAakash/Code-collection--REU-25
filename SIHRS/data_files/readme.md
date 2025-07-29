First we export the csv file from here:
https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/uqq2-txqb/about_data

This is also the same source that NY Times has used in order to plot the hospitalization curves.

The key difference in their plot is that they use Health Service Area (HSA) that intersects Tuscaloosa County but we will be just sticking with the Tuscaloosa county hospitals.
Due to this difference, sometimes the hospitalization peak might be a bit higher in the NYTimes graphs.


Now, once we download the csv file, we open it using Power Query in Excel
and then use the following query to filter out the data for hospitalization according the county we are interested is as follows:

Table.SelectRows(#"Changed column type", each ([collection_week] > #date(2020, 3, 11) and [collection_week] <= #date(2021, 12, 31)) and ([total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg] <> null and [total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg] <> -999999) and ([fips_code] = 1125))


1125 is the fips code for Tuscaloosa county which we are not interested in due to the high population.

Below are some of the good candidates based on the fact that they have R_0 between 1 and 1.5, Population size less than 60,000 and have enough hospitalization data.

By enough hospitalization data, I mean we have actual numbers instead of null or -999999 which means no data or no hospitalization data reported for that 7 day period.

Good candidates along with an arbitrary/subjective rating assigned :

R_0	  Population.  FIPS  Rating
1.047	41411	  54069	  7/10
1.167	43909	  28151	  9/10
1.218	46963	  40113	  8/10
1.18	39228	  28113	  9/10
1.222	55916	  32510	  10/10
1.399	41083	  51683	  8/10

