#load packages
library(tidyverse)
library(readxl)
library(xlsx)
library(scico)
library(viridis)
library(rcompanion)
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


### Load population data
pop <- read_excel("WPP2022_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx") 
# DEfine countries
Countries_Africa <- c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Cameroon", 
                      "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Kenya", 
                      "Madagascar", "Mali", "Morocco", "Mozambique", "Namibia", 
                      "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe")

### Load Our world in data dataset on reported SARS-CoV-2 cases
Reported_Cases <- read.csv("owid-covid-data_20-05-2022.csv") %>%
  # Rename Country variable
  rename(Country = location) %>%
  # Remove data on countries not represented in the PCR data
  filter(Country %in% c("South Africa", "Botswana"))


# Mutate po data
pop_selection <- pop %>% rename(Country = 3) %>% 
  # Select countries
  filter(Country %in% Countries_Africa ) %>%
  # Select msot recent year (2021)
  filter(Year == max(Year)) %>%
  # Select collumns
  dplyr::select(Country, Year, 113) %>%
  # Show pop in millions
  mutate(Total = Total*1000)
# Rename Congo
pop_selection$Country[pop_selection$Country == "Congo"] <- "Republic of the Congo"




#### TIDY DATA ----------------------------------------------------------------
data_tidy <- data %>%
  # Change ages to completed life years
  mutate(Age = floor(as.numeric(Age))) %>% 
  # Remove SARS-CoV-2 negative samples
  filter(!is.na(Variant)) %>%
  # Change date and create binned dates
  mutate(Collection_date2 = as.POSIXct(Collection_date, format = "%d.%m.%Y")) %>%
  # Remove old samples
  filter(#year != 2020,
    Collection_date2 >= as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
    # Remove samples collected after April 2022
    Collection_date2 < as.POSIXct("2022-05-01", format = "%Y-%m-%d")) %>%
  mutate(Variant = ifelse(is.na(SARS2) & is.na(Omicron) & is.na(Delta), "Unclear",
                   ifelse(SARS2 == 0 & is.na(Omicron) & is.na(Delta), "Negative",
                   ifelse(is.na(Omicron) | is.na(Delta), "Unclear",
                   ifelse(Omicron == 1 & Delta == 0, "Omicron",
                   ifelse(Omicron == 0 & Delta == 1, "Delta",
                   ifelse(SARS2 == 1 & Omicron == 0 & Delta == 0, "Other",
                   ifelse(SARS2 == 0 & Omicron == 0 & Delta == 0, "Negative",
                          "Unclear")))))))) %>%
  # Remove negative and unclear samples
  filter(#Inbound_Travel == 0 | Inbound_Travel == 9 | is.na(Inbound_Travel),
         Variant %in% c("Delta", "Omicron", "Other")) %>%
  # Calculate days post first of June
  mutate(Days_post_June = as.numeric(round(difftime(Collection_date2, 
                                                    "2021-06-01", 
                                                    units= "days"), 
                                           digits = 0)),
         # Classify Omicron category based on all three PCR results 
         Omicron_reclass = ifelse(Variant == "Omicron", 1, 0)) #%>%
  # Remove traveler
  #filter(Inbound_Travel %in% c(9, 0, NA))



#### IDENTIFY DATE FOR OMICRON DOMINANCE --------------------------------------
### Define function to identify date of dominance
MODEL_DOMINANCE <- function(dataset, CountryX) {
  # Fit glm
  glmX = glm(Omicron ~ Days_post_June, data = dataset %>% filter(Country == CountryX), 
             family = binomial(link = "logit"))
  # Create empty dataset
  newdata = data.frame(Days_post_June = 0:320) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
  predicted <- as.data.frame(predict(glmX, newdata, #se.fit=TRUE,
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]   
  predicted <- predicted %>% 
    rename(Day_Difference = 1) %>%
    filter(Value > 0.5) %>%
    filter(Value == min(Value)) %>%
    mutate(Country = CountryX)
  # Return dataframe
  return(predicted)
}


## Aplly model function to the dataset
dom.Alg <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Algeria")
dom.Ang <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Angola")
dom.Ben <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Benin") 
dom.BuF <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Burkina Faso")
dom.Cam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Cameroon")
dom.Eth <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ethiopia")
#dom.Gab <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gabon") # No Omicron
dom.Gam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gambia")
dom.Gha <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ghana")
dom.Gui <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Guinea")
dom.Ken <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Kenya")
dom.Mad <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Madagascar")
dom.Mal <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mali") 
dom.Mor <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Morocco") 
dom.Moz <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mozambique") 
dom.Nam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Namibia") 
dom.Nig <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Niger") 
dom.Cog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Republic of the Congo") 
dom.Sen <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Senegal") 
dom.SAf <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "South Africa") 
dom.Tog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Togo") 
dom.Uga <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Uganda") 
dom.Zim <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Zimbabwe")

# Summarize takeover modelling results (Gabon excluded)
Dominance_Modelling <- rbind(dom.Alg, dom.Ang, dom.Ben, dom.BuF, dom.Cam, dom.Eth,
                             dom.Gam, dom.Gha, dom.Gui, dom.Ken, dom.Mal, dom.Mor, 
                             dom.Moz, dom.Nam, dom.Nig, dom.Cog, dom.Sen, dom.SAf, 
                             dom.Tog, dom.Uga, dom.Mad, dom.Zim) %>%
  # Reorder and select data
  dplyr::select(3, 1) 

# Merge modelling data and population
Dominance_Modelling_pop <- merge(Dominance_Modelling, pop_selection, by = "Country") %>%
  # Make days numeric
  mutate(Day_Difference = as.numeric(Day_Difference)) %>%
  # Calculate first case
  mutate(First_case = round(Day_Difference - (log10(0.5/(1/Total))/log10(2))*3), digits = 0) %>%
  # Define date
  mutate(First_Case_date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(First_case),
         First_Case_date_30 = First_Case_date - lubridate::days(30))

# Save results
#write.csv(Dominance_Modelling,
#          "~/PATH\\BA1_Dominance_Overview.csv")  


#### IDENTIFY FIRST CASES AND FIRST EXPECTED CASES BASED ON DOMINANCE ---------
## First Omicron
First_Omicron <- as.data.frame(data_tidy %>% 
                                 filter(Variant == "Omicron") %>%
                                 group_by(Country) %>%
                                 summarise(min(Collection_date2))) %>%
  rename(Collection_date = 2) 

## First expected Omicron
First_Expected_Omicron <- Dominance_Modelling_pop %>%
  # Select needed collumns
  dplyr::select(Country, First_Case_date)

## Merge data
First_observed_expected_Omicron <- merge(First_Omicron,
                                         First_Expected_Omicron,
                                         by = "Country") %>%
  # Rename collumns
  rename(First_BA1_tested = 2,
         First_BA1_Expected = 3)

## Export data on first observed and epected BA.1
#write.csv(First_observed_expected_Omicron,
#          "~/PATH\\First_Cases_25-04-2023.csv")  




#### REMOVE BA.1-POSITIVE SAMPLES COLLECTED BEFORE FIRST EXPECTED CASES -------
data_tidy_cleaned <- merge(data_tidy, 
                           First_Expected_Omicron,
                           by = "Country",
                           all.x = TRUE) %>%
  # Remove Gabon samples
  filter(Country != "Gabon") %>%
  # Remove samples positive for BA.1 before expected date
  filter(ifelse(Variant == "Omicron", 
                Collection_date2 >= First_Case_date,
                !is.na(Collection_date2)))
   


#### GROUP COUNTRIES BY REGION ------------------------------------------------
# Define regions
EA <- c("Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya",
                        "Madagascar", "Mauritius", "Rwanda", "Seychelles",
                        "Somalia", "South Sudan", "Sudan", "Uganda", "Tanzania")

WA <- c("Benin", "Burkina Faso", "Cape Verde", " Cote d'Ivoire", 
                        "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", 
                        "Mali", "Niger", "Nigeria", "Senegal", 
                        "Sierra Leone", "Togo")

CA <- c("Burundi", "Cameroon", "Central African Republic", 
                        "Chad", "Democratic Republic of Congo", "Equatorial Guinea", 
                        "Gabon", "Republic of the Congo", "Sao Tome and Principe")

NoA <- c("Algeria", "Egypt", "Libya", "Mauritania", "Morocco", "Tunisia")

SA <-  c("Angola", "Botswana", "Eswatini", "Lesotho", "Malawi",
                        "Mozambique", "Namibia", "Zambia", "Zimbabwe", "South Africa")


# Add region collumn to data
data_grouped <- data_tidy_cleaned %>%
  mutate(Region = ifelse(Country %in% WA, "Western Africa", 
                  ifelse(Country %in% NoA, "Northern Africa",
                  ifelse(Country %in% CA, "Central Africa",
                  ifelse(Country %in% EA, "Eastern Africa",
                  ifelse(Country %in% SA, "Southern Africa", NA)))))) %>%
  # remove MAdagascar to only show mainland EA
  filter(Country != "Madagascar")

#### EXCLUDE TRAVELER
#data_grouped <- data_grouped %>% filter(Inbound_Travel %in% c(9, 0, NA))





#### Identify cases on November 11, 2021 --------------------------------------
GLM_SA_NOV21 <- glm(Omicron ~ Days_post_June, 
                    data = data_grouped %>% filter(Region == "Southern Africa"), 
                    family = binomial(link = "logit"))
#GLM_SA_NOV21_fit <- as.data.frame(GLM_SA_NOV21$fitted.values)
  # Create empty dataset
GLM_SA_NOV21_newdata = data.frame(Days_post_June = 0:320) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
GLM_SA_NOV21_predicted <- as.data.frame(predict(GLM_SA_NOV21, 
                                                GLM_SA_NOV21_newdata, #se.fit=TRUE,
                                     type = "response")) %>%
    rename(Value = 1)
GLM_SA_NOV21_predicted <- setDT(GLM_SA_NOV21_predicted, keep.rownames = TRUE)[]   
GLM_SA_NOV21_predicted <- GLM_SA_NOV21_predicted %>% 
    rename(Days_post_June = 1) %>%
    mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
  filter(Date < "2021-11-12",  
         Date > "2021-11-06")

### Merge predicted BA.1 fraction and reported cases
SA_NOV11_2021_cases <- merge(x = GLM_SA_NOV21_predicted %>%
                               mutate(Date = as.character(Date)),
                             y = Reported_Cases %>% 
                               filter(Country == "South Africa") %>%
                               select(date, new_cases) %>%
                               mutate(date = as.character(date)),
                             by.x = "Date", 
                             by.y = "date",
                             all.x = TRUE) %>%
  mutate(BA1_cases = Value * new_cases) %>%
  summarize(BA1_Cases = sum(BA1_cases))

### Merge predicted BA.1 fraction and reported cases
BOTSW_NOV11_2021_cases <- merge(x = GLM_SA_NOV21_predicted %>%
                               mutate(Date = as.character(Date)),
                             y = Reported_Cases %>% 
                               filter(Country == "Botswana") %>%
                               select(date, new_cases) %>%
                               mutate(date = as.character(date)),
                             by.x = "Date", 
                             by.y = "date",
                             all.x = TRUE) %>%
  mutate(BA1_cases = Value * new_cases) %>%
  summarize(BA1_Cases = sum(BA1_cases))









  
  


