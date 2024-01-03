#load packages
library(tidyverse)
library(readxl)
library(readr)
library(lubridate)
library(writexl)

################################################################################
# GLEAMviz exports individual tsv files for every day and every country within
# the simulation. This script was used to load and merge all those tsv files into
# one dataset. All original tsv files were stored in one folder/path.
################################################################################


### LOAD GLEAMviz results
# Set wd
setwd("SELECT PATH")


# Define function to read tsv files and include filename which contains day and class 
customized_read_tsv <- function(file){
  read_tsv(file) %>%
    mutate(fileName = file)
}
# load data
results_countries <- list.files(full.names = TRUE) %>% # list all the files
  lapply(customized_read_tsv) %>% # read them all in with our custom function
  reduce(bind_rows) 

### LOAD GLEAMviz CITY CODES
# Set wd
setwd("SELECT PATH")

# Load data
Country_Codes <- read_tsv("md_countries.tsv") %>% dplyr::select(2,3) %>% rename(Country_Id = 1) 


## MUTATE GLEAMviz simulation data
results_tidy <- results_countries %>%
  # Extract country ID
  mutate(Country_Id = sub("-.*", "", fileName),
         Country_Id = sub(".*/", "", Country_Id),
         # Extract data level ID
         Data_Level = sub(".tsv.*", "", fileName),
         Data_Level = sub(".*-", "", Data_Level)
         # Data Level: 2 = Infectious severe, 3 = Recovered, 4 = Infectious mild
  )     


# MERGE GLEAMviz RESULTS AND CITY NAMES/IDs
results_tidy_countries <- merge(results_tidy, Country_Codes, by = "Country_Id")

# Rename countries to match them in both files
results_tidy_countries[results_tidy_countries == "Egypt, Arab Rep."] <- "Egypt"
results_tidy_countries[results_tidy_countries == "Cape Verde"] <- "Cabo Verde"
results_tidy_countries[results_tidy_countries == "Congo, Dem. Rep."] <- "Republic of Congo"
results_tidy_countries[results_tidy_countries == "Congo, Rep."] <- "Republic of the Congo"
results_tidy_countries[results_tidy_countries == "Côte d'Ivoire"] <- "Cote dâ€™Ivoire"
results_tidy_countries[results_tidy_countries == "Sao Tome and Principe"] <- "São Tomé and Princip"
results_tidy_countries[results_tidy_countries == "eSwatini"] <- "Eswatini"



# Merge results and country codes
Country_Results <- results_tidy_countries %>%
  # Select and rename needed collumns
  dplyr::select(Day = 3,
                Cummulative_Median = 7,
                Data_Level = 11,
                Country = Name) %>%
  # Mutate data level
  mutate(Data_Level = ifelse(Data_Level == 0, "Susceptible",
                             ifelse(Data_Level == 1, "Exposed",
                                    ifelse(Data_Level == 2, "Infectious_Severe",
                                           ifelse(Data_Level == 3, "Recovered",
                                                  ifelse(Data_Level == 4, "Infectious_Mild", NA)))))) %>%
  pivot_wider(names_from = Data_Level, values_from = Cummulative_Median) %>%
  mutate(Infected = Infectious_Mild + Infectious_Severe,
         Infected_Recovered = Infectious_Mild + Infectious_Severe + Recovered) %>%
  # Add Days post June (Simulation was started at 2021-11-11)
  mutate(Days_post_June = as.numeric(round(difftime("2021-11-11", "2021-06-01", units= "days"), digits = 0)) + as.numeric(Day))



#### SAVE DATASET -------------------------------------------------------------
# set wd
setwd("SELECT PATH")

# save
#write.csv(Country_Results, "SimNew1_F1_I13_countries_tidy.csv", row.names=FALSE)
# ...



