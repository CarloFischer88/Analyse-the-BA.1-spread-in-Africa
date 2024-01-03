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
library(EpiEstim)
library(incidence2)
library(matrixStats)

Sys.setenv(TZ = "Europe/London")



#### LOAD DATA ----------------------------------------------------------------
### Set working directory
setwd("SELECT PATH")
### Load dataset
data <- read.csv("dataset.csv",
                 header = TRUE,        # Whether to read the header or not
                 sep = ";",            # Separator of the values
                 dec = ".")

# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"

# Correct Morocco
data$Country[data$Country == "Marocco"] <- "Morocco"



## Define countries represented in the PCR data
## (Gabon is excluded as no BA.1 cases were observed in the sample)
Countries_Africa <- c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Cameroon", "Ethiopia", "Gabon",
                      "Gambia", "Ghana", "Kenya", "Madagascar", "Mali", "Morocco", "Mozambique",
                      "Namibia", "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe", "Guinea")
## Load Our world in data dataset on reported SARS-CoV-2 cases
Reported_Cases <- read.csv("owid-covid-data_20-05-2022.csv") %>%
  # Rename Country variable
  rename(Country = location) %>%
  # Remove data on countries not represented in the PCR data
  filter(Country %in% Countries_Africa)

# Correct Morocco
Reported_Cases$Country[Reported_Cases$Country == "Congo"] <- "Republic of the Congo"

### Load population data
pop <- read_excel("WPP2022_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx") 



#### TIDY AND MUTATE DATA -----------------------------------------------------
### Mutate population data
pop_selection <- pop %>% rename(Country = 3) %>% 
  # Select countries
  filter(Country %in% Countries_Africa ) %>%
  # Select msot recent year (2021)
  filter(Year == max(Year)) %>%
  # Select collumns
  select(Country, Year, 113) %>%
  # Show pop in millions
  mutate(Total = Total*1000)
# Rename Congo
pop_selection$Country[pop_selection$Country == "Congo"] <- "Republic of the Congo"


### Tidy PCR data
data_tidy <- data %>%
  # Change ages to completed life years
  mutate(Age = floor(as.numeric(Age))) %>% 
  # Remove SARS-CoV-2 negative samples
  filter(!is.na(Variant)) %>%
  # Change date and create binned dates
  mutate(Collection_date = as.POSIXct(Collection_date, format = "%d.%m.%Y"),
         Days_post_June = as.numeric(round(difftime(Collection_date, 
                                                    "2021-06-01", units= "days"), digits = 0))) %>%
  # Remove old samples
  filter(Collection_date >= as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
         Collection_date <= as.POSIXct("2022-03-01", format = "%Y-%m-%d")) %>%
  # Remove Zimbabwe as data are flawed
  #Country != "Zimbabwe"
  mutate(Variant = ifelse(is.na(SARS2) & is.na(Omicron) & is.na(Delta), "Unclear",
                          ifelse(SARS2 == 0 & is.na(Omicron) & is.na(Delta), "Negative",
                                 ifelse(is.na(Omicron) | is.na(Delta), "Unclear",
                                        ifelse(Omicron == 1 & Delta == 0, "Omicron",
                                               ifelse(Omicron == 0 & Delta == 1, "Delta",
                                                      ifelse(SARS2 == 1 & Omicron == 0 & Delta == 0, "Other",
                                                             ifelse(SARS2 == 0 & Omicron == 0 & Delta == 0, "Negative",
                                                                    "Unclear")))))))) %>%
  # Remove samples positive for both Delta and Omicron
  filter(Variant %in% c("Omicron", "Delta", "Other"))# %>%
# Remove traveler
#filter(Inbound_Travel == 0 | is.na(Inbound_Travel) | Inbound_Travel == 9)



#### MODEL OMICRON DOMINANCE --------------------------------------------------
### Define function to estimate date of dominance
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


### Apply model function
dom.Alg <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Algeria")
dom.Ang <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Angola") 
dom.Ben <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Benin") 
dom.Bot <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Botswana") 
dom.BuF <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Burkina Faso") 
dom.Cam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Cameroon") 
dom.Eth <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ethiopia") 
#dom.Gab <- MODEL_DOMINANCE(dataset = data_grouped, CountryX = "Gabon") # No BA.1
dom.Gam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gambia") 
dom.Gha <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ghana") 
dom.Gui <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Guinea") 
dom.Ken <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Kenya") 
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
dom.Mad <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Madagascar") 
dom.Zim <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Zimbabwe") 

# Summarize takeover modelling results
Dominance_Modelling <- rbind(dom.Alg, dom.Ang, dom.Ben, dom.Bot, dom.BuF, dom.Cam, 
                             dom.Eth, dom.Gam, dom.Gha, dom.Ken, dom.Mal, dom.Mor, 
                             dom.Moz, dom.Nam, dom.Nig, dom.Cog, dom.Sen, dom.SAf, 
                             dom.Tog, dom.Uga, dom.Mad, dom.Zim, dom.Gui) %>%
  # Reorder and select data
  select(3, 1)
# Rename Morocco
#Dominance_Modelling$Country[Dominance_Modelling$Country == "Marocco"] <- "Morocco"


### Merge modelling data and population
Dominance_Modelling_pop <- merge(Dominance_Modelling, pop_selection, by = "Country") %>%
  # Make days numeric
  mutate(Day_Difference = as.numeric(Day_Difference)) %>%
  # Show dominance date
  mutate(Date_Dominance = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Day_Difference)) %>%
  # Calculate first case
  mutate(First_expected_case = round(Day_Difference - (log10(0.5/(1/Total))/log10(2))*3), digits = 0) %>%
  # Define date
  mutate(First_expected_case_date = 
           as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(First_expected_case),
         First_expected_case_date_30 = Date_Dominance - lubridate::days(30))




#### Find first Omicron cases -------------------------------------------------
First_Omicron <- as.data.frame(data_tidy %>% 
                                 filter(Omicron == 1) %>%
                                 group_by(Country) %>%
                                 summarise(min(Collection_date))) %>%
  rename(First_Omicron_date = 2) 

#merge Dominance modelling and first Omicron case data
Dominance_Modelling_pop_first_case <- merge(Dominance_Modelling_pop, First_Omicron, by = "Country")



#### Remove BA.1 positive samples which were collected before the estimated ---
#### first case of the specific country ---------------------------------------
filter_dates <- Dominance_Modelling_pop_first_case %>% select(Country, First_expected_case_date)

# Merge testing data and earliest expected BA.1 date
data_filtered <- merge(data_tidy, filter_dates, by = "Country") %>%
  # Remove samples which were BA.1 positive before the first expected case
  filter(!(Omicron == 1 & Collection_date < First_expected_case_date))




#### Find first Omicron cases after removing too early detections -------------
First_Omicron_cleaned <- as.data.frame(data_filtered %>% 
                                         filter(Omicron == 1) %>%
                                         group_by(Country) %>%
                                         summarise(min(Collection_date))) %>%
  rename(First_Omicron_date = 2) 



##### ESTIMATE OMICRON FRACTION OF REPORTED SARS-CoV-2 CASES ------------------
### DEFINE FUNCTION
MODEL_OMICRON <- function(dataset, CountryX, FirstCase) {
  ## Select Omicron data
  data_loaded <- dataset %>%
    filter(Country == CountryX) %>%
    # Calculate difference between first June 2021 signal and collection date
    mutate(Diff_1st_Case = 
             round(as.numeric(difftime((as.POSIXct(Collection_date, format = "%Y-%m-%d")),
                                       as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
                                       units= "days")), digits = 0))
  ## Select data on reported cases
  data_cases <- Reported_Cases %>%
    filter(Country == CountryX) %>%
    # Select collumns with needed data
    select(3,4,13) %>% 
    # Transform dates into posixct format
    mutate(Collection_date = as.POSIXct(date, format = "%Y-%m-%d")) %>%
    # Calculate difference between June 2021 and reporting date 
    mutate(Diff_1st_Case = round(as.numeric(difftime(Collection_date,
                                                     as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
                                                     units= "days")), digits = 0)) %>%
    # Remove samples which were calculated more than 50 days before the first BA.1 signal
    filter(Diff_1st_Case >= -50) %>%
    # Select needed collumns
    select(Diff_1st_Case, new_cases_smoothed_per_million)

  
  #################### OMICRON
  glmO = glm(Omicron ~ Diff_1st_Case, data = data_loaded, 
             family = binomial(link = "logit"))
  ## Create dummy dataset 
  newdata = data.frame(Diff_1st_Case = 
                         as.numeric(0:600))
  ## Predict BA.1 fraction over time
  predicted <- as.data.frame(predict(glmO, newdata, 
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]      
  predicted <- predicted %>% 
    rename(Diff_1st_Case = 1) %>%
    mutate(Diff_1st_Case = as.numeric(Diff_1st_Case)) #%>%
  ## Merge data on cases and BA.1 fraction
  MERGED <- merge(predicted, data_cases, by = "Diff_1st_Case") %>%
    # Calculate BA.1 cases by multiplying BA.1 fraction with reported cases
    mutate(Omicron_Cases = round((Value*new_cases_smoothed_per_million
    ), 
    digits = 3))
  ## Identify dominance day
  Day_Dominance <- MERGED %>%
    filter(Value >= 0.5) %>%
    slice(1:1)
  ## Create days post/after dominance
  OMICRON <- MERGED %>%
    mutate(Diff_1st_Case = Diff_1st_Case - Day_Dominance$Diff_1st_Case)
  
  ################## DELTA
  glmD = glm(Delta ~ Diff_1st_Case, data = data_loaded, 
             family = binomial(link = "logit"))
  # Predict BA.1 fraction over time
  predicted_D <- as.data.frame(predict(glmD, newdata, 
                                       type = "response")) %>%
    rename(Value = 1)
  predicted_D <- setDT(predicted_D, keep.rownames = TRUE)[]      
  predicted_D <- predicted_D %>% 
    rename(Diff_1st_Case = 1) %>%
    # Make day difference numeric
    mutate(Diff_1st_Case = as.numeric(Diff_1st_Case)) #%>%
  # Merge data on cases and Omicron fraction
  DELTA <- merge(predicted_D, data_cases, by = "Diff_1st_Case") %>%
    mutate(Delta_Cases = round((Value*new_cases_smoothed_per_million
    ), 
    digits = 3)) %>%
    mutate(Diff_1st_Case = Diff_1st_Case - Day_Dominance$Diff_1st_Case)
  # Merge Omicron and Delta data
  MERGED <- merge(OMICRON, DELTA, by = "Diff_1st_Case") %>%
    # Select and rename collumns
    select(Days_post_Dom = 1, BA1_fraction = 2, 3, 4,
           Delta_fraction = 5, 7)
  ## Return dataframe
  return(MERGED)
}


##### APPLY FUNCTION TO ESTIMATE BA.1 CASES -----------------------------------
## Algeria
CALC_ALG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Algeria")
# Plot values for quality control
ggplot(data = CALC_ALG) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Angola
CALC_ANG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Angola")
# Plot values for quality control
ggplot(data = CALC_ANG) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Benin
CALC_BEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Benin")
# Plot values for quality control
ggplot(data = CALC_BEN) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Botswana
CALC_BOT <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Botswana")
# Plot values for quality control
ggplot(data = CALC_BOT) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Burkina Faso
CALC_BUF <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Burkina Faso")
# Plot values for quality control
ggplot(data = CALC_BUF) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Cameroon
CALC_CAM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Cameroon")
# Plot values for quality control
ggplot(data = CALC_CAM) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Ethiopia
CALC_ETH <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Ethiopia")
# Plot values for quality control
ggplot(data = CALC_ETH) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Gambia
CALC_GAM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Gambia")
# Plot values for quality control
ggplot(data = CALC_GAM) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Ghana
CALC_GHA <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Ghana")
# Plot values for quality control
ggplot(data = CALC_GHA) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Guinea
CALC_GUI <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Guinea")
# Plot values for quality control
ggplot(data = CALC_GUI) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Kenya
CALC_KEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Kenya")
# Plot values for quality control
ggplot(data = CALC_KEN) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")


## Madagascar
CALC_MAD <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Madagascar")

# Plot values for quality control
ggplot(data = CALC_MAD) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Mali
CALC_MAL <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Mali")
# Plot values for quality control
ggplot(data = CALC_MAL) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

## Morocco
CALC_MOR <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Morocco")
# Plot values for quality control
ggplot(data = CALC_MOR) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Mozambique
CALC_MOZ <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Mozambique")
# Plot values for quality control
ggplot(data = CALC_MOZ) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Namibia
CALC_NMB <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Namibia")
# Plot values for quality control
ggplot(data = CALC_NMB) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Niger
CALC_NIG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Niger")
# Plot values for quality control
ggplot(data = CALC_NIG) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Republic of the Congo
CALC_COG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Republic of the Congo")
# Plot values for quality control
ggplot(data = CALC_COG) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Senegal
CALC_SEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Senegal")
# Plot values for quality control
ggplot(data = CALC_SEN) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# South Africa
CALC_SAF <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "South Africa")
# Plot values for quality control
ggplot(data = CALC_SAF) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Togo
CALC_TOG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Togo")
# Plot values for quality control
ggplot(data = CALC_TOG) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Uganda
CALC_UGA <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Uganda")
# Plot values for quality control
ggplot(data = CALC_UGA) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")

# Zimbabwe
CALC_ZIM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Zimbabwe")
# Plot values for quality control
ggplot(data = CALC_ZIM) +
  geom_point(aes(x = Days_post_Dom, y = BA1_fraction), colour = "red")+
  geom_point(aes(x = Days_post_Dom, y = Delta_fraction), colour = "blue")




#### Merge predicted values for BA.1 ------------------------------------------
## Join predicted country data in one list to join them in a data frame
Predicted_list_BA1 <- list(CALC_ALG,
                           #CALC_ANG,
                           CALC_BEN,
                           #CALC_BOT,
                           CALC_BUF,
                           CALC_CAM,
                           CALC_ETH,
                           #CALC_GAM,
                           CALC_GHA,
                           CALC_GUI,
                           CALC_KEN,
                           CALC_MAL,
                           CALC_MOR,
                           #CALC_NMB,
                           #CALC_NIG,
                           CALC_SEN,
                           CALC_SAF,
                           CALC_TOG,
                           CALC_UGA,
                           CALC_MOZ,
                           CALC_COG,
                           CALC_MAD#,
                           #CALC_ZIM
)


## Merge all data frames of predicted cases
Predicted_merged <- Predicted_list_BA1 %>% reduce(full_join, by='Days_post_Dom') %>%
  # Select needed data
  select(#-contains("smoothed"),
    -contains("Value")) %>%
  arrange(Days_post_Dom)

## Extract Omicron cases
Predicted_merged_BA1_cases <- Predicted_merged %>%
  select(contains("Omicron_Cases") | contains("Days_post_Dom"))

## Calculate median of predicted BA.1 cases over time 
Predicted_merged_BA1_cases$Mean <- round(rowMeans(as.matrix(Predicted_merged_BA1_cases[,1:(ncol(Predicted_merged_BA1_cases)-1)]), na.rm = TRUE), digits = 3)

## Select needed collumns
Predicted_merged_BA1_cases <- Predicted_merged_BA1_cases %>%
  select(Days_post_Dom, Mean_BA1_Cases = Mean)



#### IDENTIFY COUNTRIES WITH VERY SLOW INCREASE IN BA.1 FRACTION --------------
### Define function to estimate time between <10% BA.1 and >80% BA.1 
### (3 doubling times)
### Define function to estimate date of dominance
MODEL_01_to_08 <- function(dataset) {
  # Identify date of roughly 10% BA.1
  BA1_01 <- dataset %>%
    filter(BA1_fraction <= 0.1) %>%
    filter(Days_post_Dom == max(Days_post_Dom)) %>%
    select(Days_post_Dom)
  # Identify date of roughly 80% BA.1
  BA1_08 <- dataset %>%
    filter(BA1_fraction >= 0.8) %>%
    filter(Days_post_Dom == min(Days_post_Dom)) %>%
    select(Days_post_Dom)
  BA1_01_to_08 <- as.data.frame(as.numeric(BA1_08$Days_post_Dom) - as.numeric(BA1_01$Days_post_Dom))
  # Add Country
  #BA1_01_to_08 <- BA1_01_to_08 %>% mutate(Country = print(names(i)))
  # Return dataframe
  return(BA1_01_to_08)
}

## APPLY FUNCTION TO DATASETS
MODEL_01_to_08_ALG <- as.data.frame(MODEL_01_to_08(CALC_ALG)) %>% mutate(Country = "ALG")
MODEL_01_to_08_ANG <- as.data.frame(MODEL_01_to_08(CALC_ANG)) %>% mutate(Country = "ANG")
MODEL_01_to_08_BEN <- as.data.frame(MODEL_01_to_08(CALC_BEN)) %>% mutate(Country = "BEN")
MODEL_01_to_08_BOT <- as.data.frame(MODEL_01_to_08(CALC_BOT)) %>% mutate(Country = "BOT")
MODEL_01_to_08_BUF <- as.data.frame(MODEL_01_to_08(CALC_BUF)) %>% mutate(Country = "BUF")
MODEL_01_to_08_CAM <- as.data.frame(MODEL_01_to_08(CALC_CAM)) %>% mutate(Country = "CAM")
MODEL_01_to_08_ETH <- as.data.frame(MODEL_01_to_08(CALC_ETH)) %>% mutate(Country = "ETH")
MODEL_01_to_08_GAM <- as.data.frame(MODEL_01_to_08(CALC_GAM)) %>% mutate(Country = "GAM")
MODEL_01_to_08_GHA <- as.data.frame(MODEL_01_to_08(CALC_GHA)) %>% mutate(Country = "GHA")
MODEL_01_to_08_GUI <- as.data.frame(MODEL_01_to_08(CALC_GUI)) %>% mutate(Country = "GUI")
MODEL_01_to_08_KEN <- as.data.frame(MODEL_01_to_08(CALC_KEN)) %>% mutate(Country = "KEN")
MODEL_01_to_08_MAL <- as.data.frame(MODEL_01_to_08(CALC_MAL)) %>% mutate(Country = "MAL")
MODEL_01_to_08_MOR <- as.data.frame(MODEL_01_to_08(CALC_MOR)) %>% mutate(Country = "MOR")
MODEL_01_to_08_NMB <- as.data.frame(MODEL_01_to_08(CALC_NMB)) %>% mutate(Country = "NMB")
MODEL_01_to_08_NIG <- as.data.frame(MODEL_01_to_08(CALC_NIG)) %>% mutate(Country = "NIG")
MODEL_01_to_08_SEN <- as.data.frame(MODEL_01_to_08(CALC_SEN)) %>% mutate(Country = "SEN")
MODEL_01_to_08_SAF <- as.data.frame(MODEL_01_to_08(CALC_SAF)) %>% mutate(Country = "SAF")
MODEL_01_to_08_TOG <- as.data.frame(MODEL_01_to_08(CALC_TOG)) %>% mutate(Country = "TOG")
MODEL_01_to_08_UGA <- as.data.frame(MODEL_01_to_08(CALC_UGA)) %>% mutate(Country = "UGA")
MODEL_01_to_08_MOZ <- as.data.frame(MODEL_01_to_08(CALC_MOZ)) %>% mutate(Country = "MOZ")
MODEL_01_to_08_COG <- as.data.frame(MODEL_01_to_08(CALC_COG)) %>% mutate(Country = "COG")
MODEL_01_to_08_MAD <- as.data.frame(MODEL_01_to_08(CALC_MAD)) %>% mutate(Country = "MAD")

## MERGE RESULTS
MODEL_01_to_08_All <- rbind(
  MODEL_01_to_08_ALG,
  MODEL_01_to_08_ANG,
  MODEL_01_to_08_BEN,
  MODEL_01_to_08_BOT,
  MODEL_01_to_08_BUF,
  MODEL_01_to_08_CAM,
  MODEL_01_to_08_ETH,
  MODEL_01_to_08_GAM,
  MODEL_01_to_08_GHA,
  MODEL_01_to_08_GUI,
  MODEL_01_to_08_KEN,
  MODEL_01_to_08_MAL,
  MODEL_01_to_08_MOR,
  MODEL_01_to_08_NMB,
  MODEL_01_to_08_NIG,
  MODEL_01_to_08_SEN,
  MODEL_01_to_08_SAF,
  MODEL_01_to_08_TOG,
  MODEL_01_to_08_UGA,
  MODEL_01_to_08_MOZ,
  MODEL_01_to_08_COG,
  MODEL_01_to_08_MAD
) %>%
  rename(Day_Difference = 1, Country = 2)




#### SUMMARIZE TOTAL NUMBER OF SARS2 CASES ------------------------------------
## Extract Total cases
Predicted_merged_All_cases <- Predicted_merged %>%
  select(contains("new_cases_smoothed") | contains("Days_post_Dom"))

## Calculate mean of total cases over time 
Predicted_merged_All_cases$Mean <- round(rowMeans(as.matrix(Predicted_merged_All_cases[,1:(ncol(Predicted_merged_All_cases)-1)]), na.rm = TRUE), digits = 3)

## Select needed collumns
Predicted_merged_All_cases <- Predicted_merged_All_cases %>%
  select(Days_post_Dom, Total_Cases = Mean)



#### Merge predicted values for Delta -----------------------------------------
## Join predicted country data in one list to join them in a data frame
Predicted_list_Delta <- list(CALC_ALG,
                             #CALC_ANG,
                             CALC_BEN,
                             #CALC_BOT,
                             #CALC_BUF,
                             #CALC_CAM,
                             #CALC_ETH,
                             CALC_GAM,
                             CALC_GHA,
                             CALC_GUI,
                             CALC_KEN,
                             CALC_MAL,
                             CALC_MOR,
                             #CALC_NMB,
                             #CALC_NIG,
                             CALC_SEN,
                             CALC_SAF,
                             CALC_TOG,
                             CALC_UGA,
                             #CALC_MOZ,
                             #CALC_COG,
                             CALC_MAD#,
                             #CALC_ZIM#
)


## Merge all data frames of predicted Delta cases
Predicted_merged_Delta <- Predicted_list_Delta %>% reduce(full_join, by='Days_post_Dom') %>%
  # Select needed data
  select(#-contains("smoothed"),
    -contains("Value")) %>%
  arrange(Days_post_Dom)

## Extract Omicron cases
Predicted_merged_Delta_cases <- Predicted_merged_Delta %>%
  select(contains("Delta_Cases") | contains("Days_post_Dom"))

## Calculate median of predicted BA.1 cases over time 
Predicted_merged_Delta_cases$Mean <- round(rowMeans(as.matrix(Predicted_merged_Delta_cases[,1:(ncol(Predicted_merged_Delta_cases)-1)]), na.rm = TRUE), digits = 3)

## Select needed collumns
Predicted_merged_Delta_cases <- Predicted_merged_Delta_cases %>%
  select(Days_post_Dom, Delta_Cases_Mean = Mean)




#### MERGE ALL PREDICTED AND REPORTED CASE NUMBERS-----------------------------
Predicted_Cases_merged <- merge(Predicted_merged_All_cases,
                                Predicted_merged_BA1_cases,
                                by = "Days_post_Dom")

Predicted_Cases_merged <- merge(Predicted_Cases_merged,
                                Predicted_merged_Delta_cases,
                                by = "Days_post_Dom")

## Calculate case numbers for non-Delta + non-Omicron variants
Predicted_Cases_merged <- Predicted_Cases_merged %>%
  mutate(Mean_Other_Cases = Total_Cases - Mean_BA1_Cases - Delta_Cases_Mean)

## Remove results >50 days before dominance
Predicted_Cases_merged_selected <- Predicted_Cases_merged %>%
  filter(Days_post_Dom > -51,
         !is.na(Total_Cases))

## Transform results to long for mat for plotting
Predicted_Cases_merged_selected_long <- Predicted_Cases_merged_selected %>%
  pivot_longer(cols = c("Mean_BA1_Cases", "Delta_Cases_Mean"),
               names_to = "Variant",
               values_to = "Cases")


##### ESTIMATE R0 FOR BA.1 ---------------------------------------------------
## Define the sliding window size for estimation of R0
t_start <- seq(2, nrow(Predicted_merged)-7)   
t_end <- t_start + 7

## Estimate R0
res_parametric_si <- estimate_R(Predicted_Cases_merged_selected$Mean_BA1_Cases, 
                                # Select prediction method
                                method="parametric_si",
                                config = make_config(list(
                                  # define serial interval length
                                  mean_si = 3.3, 
                                  # define serial interval deviation
                                  std_si = 2.4#, 
                                  #t_start = t_start, 
                                  #t_end = t_end
                                )))

## Extract needed results from R0 estimation results
R_results <- as.data.frame(res_parametric_si$R) %>%
  rename(Mean = 3, Q005 = 6, Q0025= 5, Q025 = 7, Median =8, Q075 = 9, 
         Q095 = 10, Q0975 = 11) %>%
  # Correct time as R0 prediction starts with day 0 instead of e.g., -50
  mutate(Time = t_end - 50)

I_results <- as.data.frame(res_parametric_si$I) %>%
  rename(Incidence = 1)
I_results$row_names <- row.names(I_results) 
I_results <- as.data.frame(I_results) %>%
  rename(Time = 2) %>%
  mutate(Time = as.numeric(Time)-50)

plot(res_parametric_si, legend = FALSE)





# Calculate max R
maxR <- R_results %>% filter(Time >-30 & Time < 50) %>%
  filter(Median == max(Median))


################### TO COUNTRY-specific R0s ----------------------------------

## Merge all data frames of predicted cases
R_GHA <- CALC_GHA %>% 
  # Select needed collumns
  select(Days_post_Dom, Omicron_Cases) %>%
  filter(Days_post_Dom > -50)

## Estimate R0
res_parametric_si_GHA <- estimate_R(R_GHA$Omicron_Cases, 
                                    # Select prediction method
                                    method="parametric_si",
                                    config = make_config(list(
                                      # define serial interval length
                                      mean_si = 3.3, 
                                      # define serial interval deviation
                                      std_si = 2.4#, 
                                      #t_start = t_start, 
                                      #t_end = t_end
                                    )))

## Extract needed results from R0 estimation results
R_results_GHA <- as.data.frame(res_parametric_si$R) %>%
  rename(Mean = 3, Q005 = 6, Q0025= 5, Q025 = 7, Median =8, Q075 = 9, 
         Q095 = 10, Q0975 = 11) %>%
  # Correct time as R0 prediction starts with day 0 instead of e.g., -50
  mutate(Time = t_end - 50)


##### FUNCTION
R_PRED_OMICRON <- function(dataset) {
  
  ## Merge all data frames of predicted cases
  Input_data <- dataset %>% 
    # Select needed collumns
    select(Days_post_Dom, Omicron_Cases) %>%
    filter(Days_post_Dom > -71,
           Days_post_Dom < 150,
           !is.na(Omicron_Cases))
  
  # Define sliding window size
  t_start <- seq(2, nrow(Input_data)-13)   
  t_end <- t_start + 13  
  ## Estimate R0
  res_parametric_si_FUN <- estimate_R(Input_data$Omicron_Cases, 
                                      # Select prediction method
                                      method="parametric_si",
                                      config = make_config(list(
                                        # define serial interval length
                                        mean_si = 3.3, 
                                        # define serial interval deviation
                                        std_si = 2.4, 
                                        t_start = t_start, 
                                        t_end = t_end
                                      )))
  
  ## Extract needed results from R0 estimation results
  R_PRED_results <- as.data.frame(res_parametric_si_FUN$R) %>%
    rename(Mean = 3, Q005 = 6, Q0025= 5, Q025 = 7, Median =8, Q075 = 9, 
           Q095 = 10, Q0975 = 11) %>%
    # Correct time as R0 prediction starts with day 0 instead of e.g., -50
    mutate(Time = t_end - 70) %>%
    # Select needed collumns
    select(Time, Median)
  # Define country name
  NAME = deparse(substitute(dataset))
  NAME = gsub(".*CALC_", "", NAME)
  # Smoothen output
  loessMod.2 <- loess(Median ~ Time, data = R_PRED_results, span = 0.4)
  Smoothed_Median_R <- as.data.frame(loessMod.2$x)
  Smoothed_Median_R$Median <- loessMod.2$fitted
  
  # Rename collumn
  R_PRED_results2 <- Smoothed_Median_R %>%
    rename_with(~NAME,
                .cols = Median)
  ## Return dataframe
  return(R_PRED_results2)
  #return(paste(deparse(substitute(dataset)), "Median_R", sep = "_"))
}

R_ALG <- R_PRED_OMICRON(CALC_ALG)
R_ANG <- R_PRED_OMICRON(CALC_ANG)
R_BEN <- R_PRED_OMICRON(CALC_BEN)
R_BOT <- R_PRED_OMICRON(CALC_BOT)
R_BUF <- R_PRED_OMICRON(CALC_BUF)
R_CAM <- R_PRED_OMICRON(CALC_CAM)
R_ETH <- R_PRED_OMICRON(CALC_ETH)
R_GAM <- R_PRED_OMICRON(CALC_GAM)
R_GHA <- R_PRED_OMICRON(CALC_GHA)
R_GUI <- R_PRED_OMICRON(CALC_GUI)
R_KEN <- R_PRED_OMICRON(CALC_KEN)
R_MAL <- R_PRED_OMICRON(CALC_MAL)
R_MOR <- R_PRED_OMICRON(CALC_MOR)
R_NMB <- R_PRED_OMICRON(CALC_NMB)
R_NIG <- R_PRED_OMICRON(CALC_NIG)
R_SEN <- R_PRED_OMICRON(CALC_SEN)
R_SAF <- R_PRED_OMICRON(CALC_SAF)
R_TOG <- R_PRED_OMICRON(CALC_TOG)
R_UGA <- R_PRED_OMICRON(CALC_UGA)
R_MOZ <- R_PRED_OMICRON(CALC_MOZ)
R_COG <- R_PRED_OMICRON(CALC_COG)
R_MAD <- R_PRED_OMICRON(CALC_MAD)
#R_ZIM <- R_PRED_OMICRON(CALC_ZIM)

# Create lsit from predicted Rs
R_MERGED_list <- list(R_ALG,
                      #R_ANG,
                      R_BEN,
                      #R_BOT,
                      R_BUF,
                      R_CAM,
                      R_ETH,
                      #R_GAM,
                      R_GHA,
                      R_GUI,
                      R_KEN,
                      R_MAL,
                      R_MOR,
                      #R_NMB,
                      #R_NIG,
                      R_SEN,
                      R_SAF,
                      R_TOG,
                      R_UGA,
                      R_MOZ,
                      R_COG,
                      R_MAD#,
                      #R_ZIM
)

R_MERGED <- R_MERGED_list %>% reduce(full_join, by = 'Time') 

## Calculate median of Rt
#R_MERGED$Median <- round(rowMedians(as.matrix(R_MERGED[,2:(ncol(R_MERGED))]), na.rm = TRUE), digits = 3)



## Prepare medians for plotting
R_MERGED_long <- R_MERGED %>%
  #select(-Median) %>%
  pivot_longer(cols = !Time,
               names_to = "Country") %>%
  # Remove entries with missing value (NA)
  filter(!is.na(value))

## Calculate Median R and 95% confidence intervals
R_MEDIAN_CIs95 <- groupwiseMedian(value ~ Time,
                                  data       = R_MERGED_long,
                                  conf       = 0.95,
                                  R          = 5000,
                                  percentile = TRUE,
                                  bca        = FALSE,
                                  basic      = FALSE,
                                  normal     = FALSE,
                                  wilcox     = FALSE,
                                  digits     = 3)
## Calculate Median R and 75% confidence intervals
R_MEDIAN_CIs75 <- groupwiseMedian(value ~ Time,
                                  data       = R_MERGED_long,
                                  conf       = 0.75,
                                  R          = 5000,
                                  percentile = TRUE,
                                  bca        = FALSE,
                                  basic      = FALSE,
                                  normal     = FALSE,
                                  wilcox     = FALSE,
                                  digits     = 3)

## Merge median R and CI data
R_MEDIAN_CIs <- merge(R_MEDIAN_CIs95 %>% select(Time, Median,
                                                CI95 = Percentile.upper,
                                                CI05 = Percentile.lower),
                      R_MEDIAN_CIs75 %>% select(Time,
                                                CI75 = Percentile.upper,
                                                CI25 = Percentile.lower),
                      by = "Time")




## Calculate median Rt
medRt <- R_MEDIAN_CIs %>%
  filter(Time %in% c(-30:0)) %>%
  #select(Median) %>%
  summarize(MedianR = median(Median))


#### PLOT R0 ESTIMATIONS AND CASES --------------------------------------------
ggplot(data = R_results, aes(x = Time, y = Median)) +
  # Plot daily Delta and BA.1 cases
  geom_bar(data = Predicted_Cases_merged_selected_long, 
           aes(x = Days_post_Dom, y = Cases/8, fill = Variant), stat = "identity",
           alpha = 0.3) +
  scale_fill_manual(values=c("#512DAD", "#3F9BBA"))+
  # add trend line for mean daily SARS-CoV-2 cases
  #geom_smooth(data = Predicted_Cases_merged_selected, 
  #            aes(x = Days_post_Dom, y = Total_Cases/5), 
  #           color = "orange",
  #           method = loess, fullrange = TRUE, se = FALSE) +
  geom_line(data = Predicted_Cases_merged_selected, 
            aes(x = Days_post_Dom, y = Total_Cases/8), 
            color = "orange",
            size = 1) +
  # add horizontal line to show R0/RT = 1
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.5) +
  geom_ribbon(data = R_MEDIAN_CIs, 
              aes(ymin = CI05,
                  ymax = CI95), 
              linetype = 0, alpha = 0.2) +
  geom_ribbon(data = R_MEDIAN_CIs, 
              aes(ymin = CI25,
                  ymax = CI75), 
              linetype = 0, alpha = 0.5) +
  #geom_line(data = R_MERGED_long, 
  #         aes(x = Time, y = value, group = Country), 
  #        color = "grey",
  #       size = 0.8,
  #      alpha = 0.8) +
  geom_line(data = R_MEDIAN_CIs, 
            aes(x = Time, y = Median), 
            color = "black",
            size = 1) +
  # Highlight max R0
  scale_x_continuous(limits = c(-30, 70), 
                     expand = c(0,0), 
                     breaks = c(-40, -20, 0, 20, 40, 60)) +
  scale_y_continuous(name="Median R0", sec.axis = sec_axis(~ 8*., name="BA.1 incidence/1 million"),
                     limits = c(0,7), 
                     expand = c(0,0)) +
  #scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6,8), expand =c(0,0))+
  theme_classic() +
  ggtitle("R0 of Omicron in Africa") +
  xlab("Days after BA.1 becoming the dominant variant") +
  #labs(tag = "maximum R0 = 4.86") +
  labs(tag = paste("Median Rt:", round(medRt$MedianR, digits = 1)), 
       sep = " ") +
  theme(plot.tag.position = c(0.1, 0.6),
        plot.tag = element_text(size = 12, hjust = 0, colour = "#B81135"),
        panel.spacing.y = unit(1.5, "lines")) 

ggsave("R0_InclTraveler.pdf", width = 20, height = 11, unit = "cm",
       dpi = 600)





###################################################################
#### SUMMARIZE Rt BY COUNTRY ---------------------------------------------------
# Calculate medians and CIs
Summarized_Rt_Country <- groupwiseMedian(value ~ Country,
                                  data       = R_MERGED_long %>% filter(Time < 0,
                                                                        Time > -31),
                                  conf       = 0.95,
                                  R          = 5000,
                                  percentile = TRUE,
                                  bca        = FALSE,
                                  basic      = FALSE,
                                  normal     = FALSE,
                                  wilcox     = FALSE,
                                  digits     = 3)

# Edit dataframe structure
Summarized_Rt_Country <- Summarized_Rt_Country %>%
  mutate(Median_Rt = paste(Median, " (", Percentile.lower, " - ", Percentile.upper, ")", sep = "")) %>%
  select(Country, Median_Rt)

# Export summarized Rt
write.csv(Summarized_Rt_Country, 
          "~/PATH\\Countrywise_Median_Rt_InclTraveller.csv", 
          row.names=FALSE)


