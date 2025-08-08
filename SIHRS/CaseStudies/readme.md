# 🏥 COVID-19 Hospitalization Data Extraction & Analysis

## 📄 Data Source

We use hospitalization data from the **[HealthData.gov COVID-19 Hospital Report](https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/uqq2-txqb/about_data)**. This dataset is also used by *The New York Times* to plot COVID-19 hospitalization trends.

> **Note:** While the NYT uses **Health Service Areas (HSAs)** that intersect a particular County(say Tuscaloosa), we focus **only on hospitals within that County**. As a result, hospitalization peaks may appear slightly lower in our analysis compared to the **[NYT Interactive: Tuscaloosa, Alabama COVID‑19 Cases](https://www.nytimes.com/interactive/2021/us/tuscaloosa-alabama-covid-cases.html)**.

---

## 🧾 Data Processing (Power Query)

After downloading the CSV file, we load it into **Power Query** in Excel and apply the following filter to extract relevant rows:

```powerquery
Table.SelectRows(
  #"Changed column type",
  each
    ([collection_week] > #date(2020, 3, 11) and
     [collection_week] <= #date(2021, 12, 31)) and
    ([total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg] <> null and
     [total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg] <> -999999) and
    ([fips_code] = 1125)
)
```

- `1125` is the **FIPS code** for **Tuscaloosa County**.
- The values `null` and `-999999` represent missing or unreported data for that 7-day period.

---

## ❌ Why Not Tuscaloosa?

Due to its **high population**, we chose not to use Tuscaloosa County. Instead, we selected smaller counties with:

- Population **less than 60,000**
- Estimated **$begin:math:text$ R_0 $end:math:text$** between **1.0 and 1.5**
- Sufficient hospitalization data (i.e., not missing or invalid)

---

## ✅ Selected Counties

| R₀     | Population | FIPS  | Rating  |
|--------|------------|-------|---------|
| 1.047  | 41,411     | 54069 | 7/10    |
| 1.167  | 43,909     | 28151 | 9/10    |
| 1.218  | 46,963     | 40113 | 8/10    |
| 1.180  | 39,228     | 28113 | 9/10    |
| 1.222  | 55,916     | 32510 | 10/10   |
| 1.399  | 41,083     | 51683 | 8/10    |

---

## 📌 Notes

- These ratings are **subjective**, based on data quality, completeness, and how well the counties match our target criteria.
- The selected counties can be used for simulations or comparative analysis across regions with different epidemic dynamics.
