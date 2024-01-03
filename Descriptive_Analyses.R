#load packages
library(tidyverse)
library(readxl)
library(writexl)
library(xlsx)
library(lubridate)
library(ciTools)
library(stats)
library(data.table)

Sys.setenv(TZ = "Europe/London")



#### LOAD DATA ----------------------------------------------------------------
### Set working directory
setwd("SELECT PATH")

### Load dataset
data <- read.csv2("dataset.csv",
                  header = TRUE)

# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"



# Define countries
Countries_Africa <- c("Algeria", "Angola", "Benin", "Botswana","Burkina Faso", "Cameroon", 
                      "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Kenya", 
                      "Madagascar", "Mali", "Morocco", "Mozambique", "Namibia", 
                      "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe")

### Load Our world in data dataset on reported SARS-CoV-2 cases
Reported_Cases <- read.csv("owid-covid-data_20-05-2022.csv") %>%
  # Rename Country variable
  rename(Country = location) %>%
  # Remove data on countries not represented in the PCR data
  filter(Country %in% Countries_Africa)

Reported_Cases$Country[Reported_Cases$Country == "Congo"] <- "Republic of the Congo"

#### TIDY DATA ----------------------------------------------------------------
# Recode Variant status
data_recoded <- data %>%
  mutate(Variant = ifelse(is.na(SARS2) & is.na(Omicron) & is.na(Delta), "Unclear",
                   ifelse(SARS2 == 0 & is.na(Omicron) & is.na(Delta), "Negative",
                   ifelse(is.na(Omicron) | is.na(Delta), "Unclear",
                   ifelse(Omicron == 1 & Delta == 0, "Omicron",
                   ifelse(Omicron == 0 & Delta == 1, "Delta",
                   ifelse(SARS2 == 1 & Omicron == 0 & Delta == 0, "Other",
                   ifelse(SARS2 == 0 & Omicron == 0 & Delta == 0, "Negative",
                          "Unclear"))))))),
         Sex = ifelse(Sex == "Male", "M",
               ifelse(Sex == "Female", "F", Sex))) %>%
  # Remove negative and unclear samples
  #filter(Variant %in% c("Delta", "Omicron", "Other")) %>%
  # Transform collection date
  mutate(Collection_date = as.POSIXct(Collection_date, format = "%d.%m.%Y")) %>%
  # Remove samples withour date
  filter(!is.na(Collection_date),
         # Remove samples collected before June 2021
         Collection_date >= as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
         # Remove samples collected after April 2022
         Collection_date < as.POSIXct("2022-05-01", format = "%Y-%m-%d")) %>%
  # Determine collection date on month level
  mutate(Collection_date_month = as.POSIXct(Collection_date, format = "%m_%Y"),
         Collection_date_month = format(Collection_date_month, "%Y-%m"))



#### SUMMARIZE DATA ----------------------------------------------------------
data_sum <- data_recoded %>%
  # group by Country
  group_by(Country) %>%
  # count toal sample number
  summarize(Total_Samples = n(),
            # count Delta positive samples
            Delta = sum(Variant == "Delta"),
            # count BA.1 positive samples
            BA.1 = sum(Variant == "Omicron"),
            # count samples positive for other variant
            Other_Variant = sum(Variant == "Other"),
            # count negative samples
            Negative = sum(Variant == "Negative"),
            # count samples with unclear result
            Unclear = sum(Variant == "Unclear")) %>%
  # Calculate samples included
  mutate(Samples_included = Total_Samples - Negative - Unclear)


### Identify when included samples were collected
Dates_Cases <- data_recoded %>%
  filter(Variant %in% c("Delta", "Omicron", "Other")) %>%
  group_by(Country) %>%
  summarize(First_collection_date = min(Collection_date),
            # Show last collection date for filtering cases
            Last_collection_date = max(Collection_date),
            # show study timeframe
            Collection_date_range =
              paste(min(Collection_date), "-", max(Collection_date), sep = " "))

### Merge summarized data with collection date range
data_sum <- merge(data_sum, Dates_Cases %>% select(Country, Collection_date_range), by = "Country", all.x = TRUE)



### Extract and summarize reported cases for studied timeframe
Cases_extracted <- merge(Reported_Cases %>% select(date, Country, new_cases),
                         Dates_Cases %>% select(-Collection_date_range), 
                         by = "Country", all.x = TRUE) %>%
  # filter data to match collection date range for each country
  filter(date >= First_collection_date,
         date <= Last_collection_date) %>%
  group_by(Country) %>%
  summarize(Total_Cases = sum(new_cases, na.rm = TRUE))


### Merge summarized data with reported cases in studied timeframe
data_sum_cases <- merge(data_sum, Cases_extracted, by = "Country", all.x = TRUE)


            
### Calculate fraction of reported cases included in this study
data_sum_cases <- data_sum_cases %>%
  mutate(Perncentage_of_cases_tested = round((Samples_included/Total_Cases)*100, digits = 2))




#### EXPORT SUMMARY TABLE -----------------------------------------------------
# Save results
write.csv(data_sum_cases,
          "~/PATH\\Cohort_Summary.csv")  



#### SUMMARIZE SEX AND AGE ----------------------------------------------------
data_sum_sex_age <- data_recoded %>%
  # Select only included samples
  filter(Variant %in% c("Delta", "Omicron", "Other")) %>%
  # group by Country
  group_by(Country) %>%
  # summarize data
  summarize(
    # calculate mean age
    Mean_Age = round(mean(as.numeric(Age), na.rm = TRUE), digits = 1),
    Female = round(sum(Sex == "F", na.rm = TRUE)/n(), digits = 3)*100,
    Male = round(sum(Sex == "M", na.rm = TRUE)/n(), digits = 3)*100,
    Unknown_Sex = round((n()-(sum(Sex == "F", na.rm = TRUE)+sum(Sex == "M", na.rm = TRUE)))/n(), digits = 3)*100)
    

### Identify Sampling sites for each country
Sampling_Sites <- data_recoded %>%
  # Select only included samples
  filter(Variant %in% c("Delta", "Omicron", "Other")) %>%
  # group by Country
  group_by(Country) %>%
  summarize(Locations = n_distinct(Location, na.rm = TRUE))


### Merge summarized data on sex and age with sampling sites
data_sum_sex_age <- merge(data_sum_sex_age, Sampling_Sites, 
                          by = "Country", all.x = TRUE)
  
# Save results
write.csv(data_sum_sex_age,
          "~/PATH\\Cohort_Sex_Age_Sites.csv")     



