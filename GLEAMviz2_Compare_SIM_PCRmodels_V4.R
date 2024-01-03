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
library(egg)

Sys.setenv(TZ = "Europe/London")






#### LOAD DATA ----------------------------------------------------------------
### Set working directory
setwd("SELECT PATH")

### Load dataset
data <- read.csv("dataset.csv",
                 header = TRUE,        # Whether to rCAd the hCAder or not
                 sep = ";",            # Separator of the values
                 dec = ".")

# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"

# Correct Morocco
data$Country[data$Country == "Marocco"] <- "Morocco"


## Define countries represented in the PCR data
## (Gabon is excluded as no BA.1 cases were observed in the sample)
Countries_Africa <- c("Algeria", "Angola", "Benin", "Gabon", "Botswana", "Burkina Faso", "Cameroon", "Ethiopia", "Gabon",
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
  # Select most recent year (2021)
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
  # CrCAte empty dataset
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
dom.Alg <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Algeria") # 224
dom.Ang <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Angola") # 151
dom.Ben <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Benin") # 161  
dom.Bot <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Botswana") # 161  
dom.BuF <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Burkina Faso") # 205
dom.Cam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Cameroon") # 194
dom.Eth <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ethiopia") # 215
#dom.Gab <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gabon") # 
dom.Gam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Gambia") # 141
dom.Gha <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Ghana") # 166
dom.Gui <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Guinea") # 152
dom.Ken <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Kenya") # 165
dom.Mal <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mali") # 193
dom.Mor <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Morocco") # 204
dom.Moz <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Mozambique") # 173
dom.Nam <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Namibia") # 164
dom.Nig <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Niger") # 218
dom.Cog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Republic of the Congo") # 192
dom.Sen <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Senegal") # 194
dom.SAf <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "South Africa") # 163
dom.Tog <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Togo") # 189
dom.Uga <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Uganda") # 190
dom.Mad <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Madagascar") # 252
dom.Zim <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Zimbabwe") # 190

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
  # CAlculate first CAse
  mutate(First_expected_Case = round(Day_Difference - (log10(0.5/(1/Total))/log10(2))*3), digits = 0) %>%
  # Define date
  mutate(First_expected_Case_date = 
           as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(First_expected_Case),
         First_expected_Case_date_30 = Date_Dominance - lubridate::days(30))




#### Find first Omicron CAses -------------------------------------------------
First_Omicron <- as.data.frame(data_tidy %>% 
                                 filter(Omicron == 1) %>%
                                 group_by(Country) %>%
                                 summarise(min(Collection_date))) %>%
  rename(First_Omicron_date = 2) 

#merge Dominance modelling and first Omicron case data
Dominance_Modelling_pop_first_Case <- merge(Dominance_Modelling_pop, First_Omicron, by = "Country")



#### Remove BA.1 positive samples which were collected before the estimated ---
#### first case of the specific country ---------------------------------------
filter_dates <- Dominance_Modelling_pop_first_Case %>% select(Country, First_expected_Case_date)

# Merge testing data and earliest expected BA.1 date
data_filtered <- merge(data_tidy, filter_dates, by = "Country") %>%
  # Remove samples which were BA.1 positive before the first expected case
  filter(!(Omicron == 1 & Collection_date < First_expected_Case_date))




#### Find first Omicron CAses after removing too CArly detections -------------
First_Omicron_cleaned <- as.data.frame(data_filtered %>% 
                                         filter(Omicron == 1) %>%
                                         group_by(Country) %>%
                                         summarise(min(Collection_date))) %>%
  rename(First_Omicron_date = 2) 




##### ESTIMATE OMICRON FRACTION OF REPORTED SARS-CoV-2 CASES ------------------
### DEFINE FUNCTION
MODEL_OMICRON <- function(dataset, CountryX) {
  ## Select Omicron data
  data_loaded <- dataset %>%
    filter(Country == CountryX) %>%
    # CAlculate difference between first June 2021 signal and collection date
    mutate(Days_post_June = 
             round(as.numeric(difftime((as.POSIXct(Collection_date, format = "%Y-%m-%d")),
                                       as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
                                       units= "days")), digits = 0))
  ## Select data on reported cases
  data_Cases <- Reported_Cases %>%
    filter(Country == CountryX) %>%
    # Select collumns with needed data
    select(3,4,13) %>% 
    # Transform dates into posixct format
    mutate(Collection_date = as.POSIXct(date, format = "%Y-%m-%d")) %>%
    # CAlculate difference between June 2021 and reporting date 
    mutate(Days_post_June = round(as.numeric(difftime(Collection_date,
                                                      as.POSIXct("2021-06-01", format = "%Y-%m-%d"),
                                                      units= "days")), digits = 0)) %>%
    # Remove samples which were calculated more than 50 days before the first BA.1 signal
    filter(Days_post_June >= -50) %>%
    # Select needed collumns
    select(Days_post_June, new_cases_smoothed_per_million)

  
  #################### OMICRON
  glmO = glm(Omicron ~ Days_post_June, data = data_loaded, 
             family = binomial(link = "logit"))
  ## Create dummy dataset 
  newdata = data.frame(Days_post_June = 
                         as.numeric(0:600))
  ## Predict BA.1 fraction over time
  predicted <- as.data.frame(predict(glmO, newdata, 
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]      
  predicted <- predicted %>% 
    rename(Days_post_June = 1) %>%
    mutate(Days_post_June = as.numeric(Days_post_June)) #%>%
  ## Merge data on cases and BA.1 fraction
  MERGED <- merge(predicted, data_Cases, by = "Days_post_June") %>%
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
    mutate(Omicron_Cases_cumsum = cumsum(Omicron_Cases))
  
  ## Return dataframe
  return(OMICRON )
}



##### APPLY FUNCTION TO ESTIMATE BA.1 CASES -----------------------------------
## Algeria
CALC_ALG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Algeria")
# Plot values for quality control
ggplot(data = CALC_ALG) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Angola
CALC_ANG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Angola")
# Plot values for quality control
ggplot(data = CALC_ANG) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Benin
CALC_BEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Benin")
# Plot values for quality control
ggplot(data = CALC_BEN) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Botswana
CALC_BOT <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Botswana")
# Plot values for quality control
ggplot(data = CALC_BOT) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Burkina Faso
CALC_BUF <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Burkina Faso")
# Plot values for quality control
ggplot(data = CALC_BUF) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Cameroon
CALC_CAM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Cameroon")
# Plot values for quality control
ggplot(data = CALC_CAM) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Ethiopia
CALC_ETH <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Ethiopia")
# Plot values for quality control
ggplot(data = CALC_ETH) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Gambia
CALC_GAM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Gambia")
# Plot values for quality control
ggplot(data = CALC_GAM) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Ghana
CALC_GHA <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Ghana")
# Plot values for quality control
ggplot(data = CALC_GHA) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Guinea
CALC_GUI <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Guinea")
# Plot values for quality control
ggplot(data = CALC_GUI) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Kenya
CALC_KEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Kenya")
# Plot values for quality control
ggplot(data = CALC_KEN) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")


## Madagascar
CALC_MAD <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Madagascar")

# Plot values for quality control
ggplot(data = CALC_MAD) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Mali
CALC_MAL <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Mali")
# Plot values for quality control
ggplot(data = CALC_MAL) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

## Morocco
CALC_MOR <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Morocco")
# Plot values for quality control
ggplot(data = CALC_MOR) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Mozambique
CALC_MOZ <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Mozambique")
# Plot values for quality control
ggplot(data = CALC_MOZ) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Namibia
CALC_NMB <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Namibia")
# Plot values for quality control
ggplot(data = CALC_NMB) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Niger
CALC_NIG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Niger")
# Plot values for quality control
ggplot(data = CALC_NIG) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Republic of the Congo
CALC_COG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Republic of the Congo")
# Plot values for quality control
ggplot(data = CALC_COG) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Senegal
CALC_SEN <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Senegal")
# Plot values for quality control
ggplot(data = CALC_SEN) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# South Africa
CALC_SAF <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "South Africa")
# Plot values for quality control
ggplot(data = CALC_SAF) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Togo
CALC_TOG <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Togo")
# Plot values for quality control
ggplot(data = CALC_TOG) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Uganda
CALC_UGA <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Uganda")
# Plot values for quality control
ggplot(data = CALC_UGA) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")

# Zimbabwe
CALC_ZIM <- MODEL_OMICRON(dataset = data_filtered,
                          CountryX = "Zimbabwe")
# Plot values for quality control
ggplot(data = CALC_ZIM) +
  geom_point(aes(x = Days_post_June, y = Omicron_Cases), colour = "red")




#### Merge predicted values for BA.1 ------------------------------------------
## Join predicted country data in one list to join them in a data frame
Predicted_list_BA1 <- list(CALC_ALG,
                           CALC_ANG,#
                           CALC_BEN,
                           CALC_BOT,
                           CALC_BUF,
                           CALC_CAM,
                           CALC_ETH,
                           CALC_GAM,#
                           CALC_GHA,
                           CALC_GUI,
                           CALC_KEN,
                           CALC_MAL,#
                           CALC_MOR,
                           CALC_NMB,#
                           CALC_NIG,#
                           CALC_SEN,
                           CALC_SAF,
                           CALC_TOG,
                           CALC_UGA,
                           CALC_MOZ,
                           CALC_COG,
                           CALC_MAD#,
                           #CALC_ZIM
)


## Merge all data frames of predicted CAses
Predicted_merged <- Predicted_list_BA1 %>% reduce(full_join, by='Days_post_June') %>%
  # Select needed data
  select(#-contains("smoothed"),
    -contains("Value")) %>%
  arrange(Days_post_June)

## Extract Omicron cases
Predicted_merged_BA1_Cases <- Predicted_merged %>%
  select(contains("Omicron_Cases") | contains("Days_post_June"))

## Calculate median of predicted BA.1 cases over time 
Predicted_merged_BA1_Cases$mean <- round(rowMeans(as.matrix(Predicted_merged_BA1_Cases[,1:(ncol(Predicted_merged_BA1_Cases)-1)]), na.rm = TRUE), digits = 3)

## Select needed collumns
Predicted_merged_BA1_Cases <- Predicted_merged_BA1_Cases %>%
  select(Days_post_June, Mean_BA1_Cases = mean)





#### MERGE PCR-BASED GLM DATA FOR AFRICAN REGIONS ##############################

#### DEFINE FUNCTION ===========================================================
MERGE_GLM_PCR <- function(dataset) {
  ## Filter countries
  data_loaded <- dataset %>%
    reduce(full_join, by='Days_post_June') %>%
    # Select needed data
    select(-contains("Value")) %>%
    # sort data
    arrange(Days_post_June) %>%
    # select collumns
    select(contains("cumsum") | contains("Days_post_June"))
  
  ## Calculate mean of predicted BA.1 cases over time 
  data_loaded$mean<- round(rowMeans(as.matrix(data_loaded[,1:(ncol(data_loaded)-1)]), na.rm = TRUE), digits = 3)
  # Select collumns
  data_loaded <- data_loaded %>% select(Days_post_June, Mean_BA1_Cases = mean)
  
  ## Return dataframe
  return(data_loaded)
}



### Merge data for Western Africa ### 
# List dataframes to be merged
Predicted_list_BA1_WA <- list(CALC_BEN,
                              CALC_BUF,
                              CALC_GAM,#
                              CALC_GHA,
                              CALC_GUI,
                              CALC_MAL,#
                              CALC_NIG,#
                              CALC_SEN,
                              CALC_TOG)
# Merge dataframes
Predicted_merged_BA1_WA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_WA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_WA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




### Merge data for Central Africa ###
# List dataframes to be merged
Predicted_list_BA1_CA <- list(CALC_CAM,
                              CALC_COG)

# Merge dataframes
Predicted_merged_BA1_CA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_CA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_CA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




### Merge data for Central Africa ###
# List dataframes to be merged
Predicted_list_BA1_CA <- list(CALC_CAM,
                              CALC_COG)

# Merge dataframes
Predicted_merged_BA1_CA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_CA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_CA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




### Merge data for Eastern Africa ###
# List dataframes to be merged
Predicted_list_BA1_EA <- list(CALC_ETH,
                              CALC_KEN,
                              CALC_UGA)

# Merge dataframes
Predicted_merged_BA1_EA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_EA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_EA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




### Merge data for Northern Africa ###
# List dataframes to be merged
Predicted_list_BA1_NoA <- list(CALC_ALG,
                               CALC_MOR)

# Merge dataframes
Predicted_merged_BA1_NoA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_NoA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_NoA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




### Merge data for Southern Africa ###
# List dataframes to be merged
Predicted_list_BA1_SA <- list(CALC_ANG,
                              #CALC_BOT,
                              CALC_NMB,
                              CALC_SAF,
                              CALC_MOZ)

# Merge dataframes
Predicted_merged_BA1_SA_selected <- MERGE_GLM_PCR(dataset = Predicted_list_BA1_SA) %>%
  # filter specific period
  filter(Days_post_June >= 163, Days_post_June <= 272)

# Plot values for quality control
ggplot(data = Predicted_merged_BA1_SA_selected) +
  geom_point(aes(x = Days_post_June, y = Mean_BA1_Cases), colour = "red")




#### LOAD TIDY GLEAMviz SIMULATION DATA ========================================
# Define Countries of interest
Africa_GLEAM <- c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "CAmeroon", "Ethiopia", "Gabon",
                  "Gambia", "Ghana", "Kenya", "Madagascar", "Mali", "Morocco", "Mozambique",
                  "Namibia", "Niger", "Republic of the Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                  "Zimbabwe", "Guinea")


# set wd
setwd("C:/Users/Carlo Fischer/OneDrive/Dokumente/PhD/01_Projects/C_SARS/CF02_BA1_spread_Africa/Analyses/Analyses_Sep2023/GLEAMviz_Evaluation2/Tidy_Results")

# load data
SimNew1_F1_I13 <- read.csv("SimNew1_F1_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew1_F1_I31 <- read.csv("SimNew1_F1_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew1_F1_I49 <- read.csv("SimNew1_F1_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew2_F1_5_I13 <- read.csv("SimNew2_F1_5_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew2_F1_5_I31 <- read.csv("SimNew2_F1_5_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew2_F1_5_I49 <- read.csv("SimNew2_F1_5_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew3_F2_25_I13 <- read.csv("SimNew3_F2_25_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew3_F2_25_I31 <- read.csv("SimNew3_F2_25_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew3_F2_25_I49 <- read.csv("SimNew3_F2_25_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew4_F3_38_I13 <- read.csv("SimNew4_F3_38_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew4_F3_38_I31 <- read.csv("SimNew4_F3_38_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew4_F3_38_I49 <- read.csv("SimNew4_F3_38_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew5_F5_06_I13 <- read.csv("SimNew5_F5_06_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew5_F5_06_I31 <- read.csv("SimNew5_F5_06_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew5_F5_06_I49 <- read.csv("SimNew5_F5_06_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew6_F7_59_I13 <- read.csv("SimNew6_F7_59_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew6_F7_59_I31 <- read.csv("SimNew6_F7_59_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew6_F7_59_I49 <- read.csv("SimNew6_F7_59_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew7_F11_39_I13 <- read.csv("SimNew7_F11_39_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew7_F11_39_I31 <- read.csv("SimNew7_F11_39_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew7_F11_39_I49 <- read.csv("SimNew7_F11_39_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew8_F17_09_I13 <- read.csv("SimNew8_F17_09_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew8_F17_09_I31 <- read.csv("SimNew8_F17_09_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew8_F17_09_I49 <- read.csv("SimNew8_F17_09_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew9_F25_63_I13 <- read.csv("SimNew9_F25_63_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew9_F25_63_I31 <- read.csv("SimNew9_F25_63_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew9_F25_63_I49 <- read.csv("SimNew9_F25_63_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew10_F38_44_I13 <- read.csv("SimNew10_F38_44_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew10_F38_44_I31 <- read.csv("SimNew10_F38_44_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew10_F38_44_I49 <- read.csv("SimNew10_F38_44_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew11_F57_67_I13 <- read.csv("SimNew11_F57_67_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew11_F57_67_I31 <- read.csv("SimNew11_F57_67_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew11_F57_67_I49 <- read.csv("SimNew11_F57_67_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew12_F86_50_I13 <- read.csv("SimNew12_F86_50_I13_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew12_F86_50_I31 <- read.csv("SimNew12_F86_50_I31_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)
SimNew12_F86_50_I49 <- read.csv("SimNew12_F86_50_I49_countries_tidy.csv") %>%
  filter(Country %in% Africa_GLEAM)



#### SUMMARIZE SIMULATION DATA BY AFRICAN REGIONS ##############################

#### DEFINE FUNCTION ===========================================================
TIDY_SIM <- function(dataset, CountryVec, Mean_Name, Median_Name) {
  ## Filter countries
  data_loaded <- dataset %>%
    filter(Country %in% CountryVec) %>% 
    # select collumns
    select(Country, Infected_Recovered, Days_post_June) %>% 
    # Adapt values for 1,000,000 people
    mutate(Infected_Recovered = 1000*Infected_Recovered) %>%
    # Make table wide format
    pivot_wider(names_from = Country, values_from = Infected_Recovered)
  ## Calculate mean and median values
  data_loaded$Mean_Sim <- round(rowMeans(as.matrix(data_loaded[,2:(ncol(data_loaded))]), na.rm = TRUE), digits = 3)
  data_loaded$Median_Sim <- round(rowMedians(as.matrix(data_loaded[,2:(ncol(data_loaded)-1)]), na.rm = TRUE), digits = 3)
  # Select collumns
  data_loaded <- data_loaded %>%
    select(Days_post_June, Mean_Sim, Median_Sim) %>%
    # rename variables
    rename(!!sym(Mean_Name) := Mean_Sim,                                    
           !!sym(Median_Name) := Median_Sim)                               
  ## Return dataframe
  return(data_loaded)
}

# Define African regions in our datasets to group results
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

SA <-  c("Angola", #"Botswana", 
         "Eswatini", "Lesotho", "Malawi",
         "Mozambique", "Namibia", "Zambia", "Zimbabwe"#, "South Africa"
)



#### WESTERN AFRICA ============================================================
#### Apply function to extract simulation data ---------------------------------
# Simulation 1.1
SimNew1_F1_I13_WA <- TIDY_SIM(dataset = SimNew1_F1_I13, CountryVec = WA, 
                              Mean_Name = "Mean_SimNew1_F1_I13", Median_Name = "Median_SimNew1_F1_I13") 
# Simulation 1.2
SimNew1_F1_I31_WA <- TIDY_SIM(dataset = SimNew1_F1_I31, CountryVec = WA, 
                              Mean_Name = "Mean_SimNew1_F1_I31", Median_Name = "Median_SimNew1_F1_I31") 
# Simulation 1.3
SimNew1_F1_I49_WA <- TIDY_SIM(dataset = SimNew1_F1_I49, CountryVec = WA, 
                              Mean_Name = "Mean_SimNew1_F1_I49", Median_Name = "Median_SimNew1_F1_I49") 

# Simulation 2.1
SimNew2_F1_5_I13_WA <- TIDY_SIM(dataset = SimNew2_F1_5_I13, CountryVec = WA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I13", Median_Name = "Median_SimNew2_F1_5_I13") 
# Simulation 2.2
SimNew2_F1_5_I31_WA <- TIDY_SIM(dataset = SimNew2_F1_5_I31, CountryVec = WA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I31", Median_Name = "Median_SimNew2_F1_5_I31") 
# Simulation 2.3
SimNew2_F1_5_I49_WA <- TIDY_SIM(dataset = SimNew2_F1_5_I49, CountryVec = WA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I49", Median_Name = "Median_SimNew2_F1_5_I49")

# Simulation 3.1
SimNew3_F2_25_I13_WA <- TIDY_SIM(dataset = SimNew3_F2_25_I13, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I13", Median_Name = "Median_SimNew3_F2_25_I13") 
# Simulation 3.2
SimNew3_F2_25_I31_WA <- TIDY_SIM(dataset = SimNew3_F2_25_I31, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I31", Median_Name = "Median_SimNew3_F2_25_I31") 
# Simulation 3.3
SimNew3_F2_25_I49_WA <- TIDY_SIM(dataset = SimNew3_F2_25_I49, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I49", Median_Name = "Median_SimNew3_F2_25_I49")

# Simulation 4.1
SimNew4_F3_38_I13_WA <- TIDY_SIM(dataset = SimNew4_F3_38_I13, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I13", Median_Name = "Median_SimNew4_F3_38_I13") 
# Simulation 4.2
SimNew4_F3_38_I31_WA <- TIDY_SIM(dataset = SimNew4_F3_38_I31, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I31", Median_Name = "Median_SimNew4_F3_38_I31") 
# Simulation 4.3
SimNew4_F3_38_I49_WA <- TIDY_SIM(dataset = SimNew4_F3_38_I49, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I49", Median_Name = "Median_SimNew4_F3_38_I49")

# Simulation 5.1
SimNew5_F5_06_I13_WA <- TIDY_SIM(dataset = SimNew5_F5_06_I13, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I13", Median_Name = "Median_SimNew5_F5_06_I13") 
# Simulation 5.2
SimNew5_F5_06_I31_WA <- TIDY_SIM(dataset = SimNew5_F5_06_I31, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I31", Median_Name = "Median_SimNew5_F5_06_I31") 
# Simulation 5.3
SimNew5_F5_06_I49_WA <- TIDY_SIM(dataset = SimNew5_F5_06_I49, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I49", Median_Name = "Median_SimNew5_F5_06_I49")

# Simulation 6.1
SimNew6_F7_59_I13_WA <- TIDY_SIM(dataset = SimNew6_F7_59_I13, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I13", Median_Name = "Median_SimNew6_F7_59_I13") 
# Simulation 6.2
SimNew6_F7_59_I31_WA <- TIDY_SIM(dataset = SimNew6_F7_59_I31, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I31", Median_Name = "Median_SimNew6_F7_59_I31") 
# Simulation 6.3
SimNew6_F7_59_I49_WA <- TIDY_SIM(dataset = SimNew6_F7_59_I49, CountryVec = WA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I49", Median_Name = "Median_SimNew6_F7_59_I49")

# Simulation 7.1
SimNew7_F11_39_I13_WA <- TIDY_SIM(dataset = SimNew7_F11_39_I13, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I13", Median_Name = "Median_SimNew7_F11_39_I13") 
# Simulation 7.2
SimNew7_F11_39_I31_WA <- TIDY_SIM(dataset = SimNew7_F11_39_I31, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I31", Median_Name = "Median_SimNew7_F11_39_I31") 
# Simulation 7.3
SimNew7_F11_39_I49_WA <- TIDY_SIM(dataset = SimNew7_F11_39_I49, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I49", Median_Name = "Median_SimNew7_F11_39_I49")

# Simulation 8.1
SimNew8_F17_09_I13_WA <- TIDY_SIM(dataset = SimNew8_F17_09_I13, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I13", Median_Name = "Median_SimNew8_F17_09_I13") 
# Simulation 8.2
SimNew8_F17_09_I31_WA <- TIDY_SIM(dataset = SimNew8_F17_09_I31, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I31", Median_Name = "Median_SimNew8_F17_09_I31") 
# Simulation 8.3
SimNew8_F17_09_I49_WA <- TIDY_SIM(dataset = SimNew8_F17_09_I49, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I49", Median_Name = "Median_SimNew8_F17_09_I49")

# Simulation 9.1
SimNew9_F25_63_I13_WA <- TIDY_SIM(dataset = SimNew9_F25_63_I13, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I13", Median_Name = "Median_SimNew9_F25_63_I13") 
# Simulation 9.2
SimNew9_F25_63_I31_WA <- TIDY_SIM(dataset = SimNew9_F25_63_I31, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I31", Median_Name = "Median_SimNew9_F25_63_I31") 
# Simulation 9.3
SimNew9_F25_63_I49_WA <- TIDY_SIM(dataset = SimNew9_F25_63_I49, CountryVec = WA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I49", Median_Name = "Median_SimNew9_F25_63_I49")

# Simulation 10.1
SimNew10_F38_44_I13_WA <- TIDY_SIM(dataset = SimNew10_F38_44_I13, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I13", Median_Name = "Median_SimNew10_F38_44_I13") 
# Simulation 10.2
SimNew10_F38_44_I31_WA <- TIDY_SIM(dataset = SimNew10_F38_44_I31, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I31", Median_Name = "Median_SimNew10_F38_44_I31") 
# Simulation 10.3
SimNew10_F38_44_I49_WA <- TIDY_SIM(dataset = SimNew10_F38_44_I49, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I49", Median_Name = "Median_SimNew10_F38_44_I49")

# Simulation 11.1
SimNew11_F57_67_I13_WA <- TIDY_SIM(dataset = SimNew11_F57_67_I13, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I13", Median_Name = "Median_SimNew11_F57_67_I13") 
# Simulation 11.2
SimNew11_F57_67_I31_WA <- TIDY_SIM(dataset = SimNew11_F57_67_I31, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I31", Median_Name = "Median_SimNew11_F57_67_I31") 
# Simulation 11.3
SimNew11_F57_67_I49_WA <- TIDY_SIM(dataset = SimNew11_F57_67_I49, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I49", Median_Name = "Median_SimNew11_F57_67_I49")

# Simulation 12.1
SimNew12_F86_50_I13_WA <- TIDY_SIM(dataset = SimNew12_F86_50_I13, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I13", Median_Name = "Median_SimNew12_F86_50_I13") 
# Simulation 12.2
SimNew12_F86_50_I31_WA <- TIDY_SIM(dataset = SimNew12_F86_50_I31, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I31", Median_Name = "Median_SimNew12_F86_50_I31") 
# Simulation 12.3
SimNew12_F86_50_I49_WA <- TIDY_SIM(dataset = SimNew12_F86_50_I49, CountryVec = WA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I49", Median_Name = "Median_SimNew12_F86_50_I49")


#### Merge extracted simulation data with PCR-based glm data -------------------
#put all data frames into list
GLEAM_GLM_WA_merged <- list(Predicted_merged_BA1_WA_selected, 
                            SimNew1_F1_I13_WA,
                            SimNew1_F1_I31_WA,
                            SimNew1_F1_I49_WA,
                            SimNew2_F1_5_I13_WA,
                            SimNew2_F1_5_I31_WA,
                            SimNew2_F1_5_I49_WA,
                            SimNew3_F2_25_I13_WA,
                            SimNew3_F2_25_I31_WA,
                            SimNew3_F2_25_I49_WA,
                            SimNew4_F3_38_I13_WA,
                            SimNew4_F3_38_I31_WA,
                            SimNew4_F3_38_I49_WA,
                            SimNew5_F5_06_I13_WA,
                            SimNew5_F5_06_I31_WA,
                            SimNew5_F5_06_I49_WA,
                            SimNew6_F7_59_I13_WA,
                            SimNew6_F7_59_I31_WA,
                            SimNew6_F7_59_I49_WA,
                            SimNew7_F11_39_I13_WA,
                            SimNew7_F11_39_I31_WA,
                            SimNew7_F11_39_I49_WA,
                            SimNew8_F17_09_I13_WA,
                            SimNew8_F17_09_I31_WA,
                            SimNew8_F17_09_I49_WA,
                            SimNew9_F25_63_I13_WA,
                            SimNew9_F25_63_I31_WA,
                            SimNew9_F25_63_I49_WA,
                            SimNew10_F38_44_I13_WA,
                            SimNew10_F38_44_I31_WA,
                            SimNew10_F38_44_I49_WA,
                            SimNew11_F57_67_I13_WA,
                            SimNew11_F57_67_I31_WA,
                            SimNew11_F57_67_I49_WA,
                            SimNew12_F86_50_I13_WA,
                            SimNew12_F86_50_I31_WA,
                            SimNew12_F86_50_I49_WA)


#### Merge all data frames in list ---------------------------------------------
GLEAM_GLM_WA_merged <- GLEAM_GLM_WA_merged %>% 
  reduce(full_join, by = "Days_post_June") %>%
  #select only data before February
  filter(Days_post_June < 245)






#### CENTRAL AFRICA ============================================================
#### Apply function to extract simulation data ---------------------------------
# Simulation 1.1
SimNew1_F1_I13_CA <- TIDY_SIM(dataset = SimNew1_F1_I13, CountryVec = CA, 
                              Mean_Name = "Mean_SimNew1_F1_I13", Median_Name = "Median_SimNew1_F1_I13") 
# Simulation 1.2
SimNew1_F1_I31_CA <- TIDY_SIM(dataset = SimNew1_F1_I31, CountryVec = CA, 
                              Mean_Name = "Mean_SimNew1_F1_I31", Median_Name = "Median_SimNew1_F1_I31") 
# Simulation 1.3
SimNew1_F1_I49_CA <- TIDY_SIM(dataset = SimNew1_F1_I49, CountryVec = CA, 
                              Mean_Name = "Mean_SimNew1_F1_I49", Median_Name = "Median_SimNew1_F1_I49") 

# Simulation 2.1
SimNew2_F1_5_I13_CA <- TIDY_SIM(dataset = SimNew2_F1_5_I13, CountryVec = CA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I13", Median_Name = "Median_SimNew2_F1_5_I13") 
# Simulation 2.2
SimNew2_F1_5_I31_CA <- TIDY_SIM(dataset = SimNew2_F1_5_I31, CountryVec = CA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I31", Median_Name = "Median_SimNew2_F1_5_I31") 
# Simulation 2.3
SimNew2_F1_5_I49_CA <- TIDY_SIM(dataset = SimNew2_F1_5_I49, CountryVec = CA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I49", Median_Name = "Median_SimNew2_F1_5_I49")

# Simulation 3.1
SimNew3_F2_25_I13_CA <- TIDY_SIM(dataset = SimNew3_F2_25_I13, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I13", Median_Name = "Median_SimNew3_F2_25_I13") 
# Simulation 3.2
SimNew3_F2_25_I31_CA <- TIDY_SIM(dataset = SimNew3_F2_25_I31, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I31", Median_Name = "Median_SimNew3_F2_25_I31") 
# Simulation 3.3
SimNew3_F2_25_I49_CA <- TIDY_SIM(dataset = SimNew3_F2_25_I49, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I49", Median_Name = "Median_SimNew3_F2_25_I49")

# Simulation 4.1
SimNew4_F3_38_I13_CA <- TIDY_SIM(dataset = SimNew4_F3_38_I13, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I13", Median_Name = "Median_SimNew4_F3_38_I13") 
# Simulation 4.2
SimNew4_F3_38_I31_CA <- TIDY_SIM(dataset = SimNew4_F3_38_I31, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I31", Median_Name = "Median_SimNew4_F3_38_I31") 
# Simulation 4.3
SimNew4_F3_38_I49_CA <- TIDY_SIM(dataset = SimNew4_F3_38_I49, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I49", Median_Name = "Median_SimNew4_F3_38_I49")

# Simulation 5.1
SimNew5_F5_06_I13_CA <- TIDY_SIM(dataset = SimNew5_F5_06_I13, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I13", Median_Name = "Median_SimNew5_F5_06_I13") 
# Simulation 5.2
SimNew5_F5_06_I31_CA <- TIDY_SIM(dataset = SimNew5_F5_06_I31, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I31", Median_Name = "Median_SimNew5_F5_06_I31") 
# Simulation 5.3
SimNew5_F5_06_I49_CA <- TIDY_SIM(dataset = SimNew5_F5_06_I49, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I49", Median_Name = "Median_SimNew5_F5_06_I49")

# Simulation 6.1
SimNew6_F7_59_I13_CA <- TIDY_SIM(dataset = SimNew6_F7_59_I13, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I13", Median_Name = "Median_SimNew6_F7_59_I13") 
# Simulation 6.2
SimNew6_F7_59_I31_CA <- TIDY_SIM(dataset = SimNew6_F7_59_I31, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I31", Median_Name = "Median_SimNew6_F7_59_I31") 
# Simulation 6.3
SimNew6_F7_59_I49_CA <- TIDY_SIM(dataset = SimNew6_F7_59_I49, CountryVec = CA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I49", Median_Name = "Median_SimNew6_F7_59_I49")

# Simulation 7.1
SimNew7_F11_39_I13_CA <- TIDY_SIM(dataset = SimNew7_F11_39_I13, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I13", Median_Name = "Median_SimNew7_F11_39_I13") 
# Simulation 7.2
SimNew7_F11_39_I31_CA <- TIDY_SIM(dataset = SimNew7_F11_39_I31, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I31", Median_Name = "Median_SimNew7_F11_39_I31") 
# Simulation 7.3
SimNew7_F11_39_I49_CA <- TIDY_SIM(dataset = SimNew7_F11_39_I49, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I49", Median_Name = "Median_SimNew7_F11_39_I49")

# Simulation 8.1
SimNew8_F17_09_I13_CA <- TIDY_SIM(dataset = SimNew8_F17_09_I13, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I13", Median_Name = "Median_SimNew8_F17_09_I13") 
# Simulation 8.2
SimNew8_F17_09_I31_CA <- TIDY_SIM(dataset = SimNew8_F17_09_I31, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I31", Median_Name = "Median_SimNew8_F17_09_I31") 
# Simulation 8.3
SimNew8_F17_09_I49_CA <- TIDY_SIM(dataset = SimNew8_F17_09_I49, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I49", Median_Name = "Median_SimNew8_F17_09_I49")

# Simulation 9.1
SimNew9_F25_63_I13_CA <- TIDY_SIM(dataset = SimNew9_F25_63_I13, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I13", Median_Name = "Median_SimNew9_F25_63_I13") 
# Simulation 9.2
SimNew9_F25_63_I31_CA <- TIDY_SIM(dataset = SimNew9_F25_63_I31, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I31", Median_Name = "Median_SimNew9_F25_63_I31") 
# Simulation 9.3
SimNew9_F25_63_I49_CA <- TIDY_SIM(dataset = SimNew9_F25_63_I49, CountryVec = CA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I49", Median_Name = "Median_SimNew9_F25_63_I49")

# Simulation 10.1
SimNew10_F38_44_I13_CA <- TIDY_SIM(dataset = SimNew10_F38_44_I13, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I13", Median_Name = "Median_SimNew10_F38_44_I13") 
# Simulation 10.2
SimNew10_F38_44_I31_CA <- TIDY_SIM(dataset = SimNew10_F38_44_I31, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I31", Median_Name = "Median_SimNew10_F38_44_I31") 
# Simulation 10.3
SimNew10_F38_44_I49_CA <- TIDY_SIM(dataset = SimNew10_F38_44_I49, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I49", Median_Name = "Median_SimNew10_F38_44_I49")

# Simulation 11.1
SimNew11_F57_67_I13_CA <- TIDY_SIM(dataset = SimNew11_F57_67_I13, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I13", Median_Name = "Median_SimNew11_F57_67_I13") 
# Simulation 11.2
SimNew11_F57_67_I31_CA <- TIDY_SIM(dataset = SimNew11_F57_67_I31, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I31", Median_Name = "Median_SimNew11_F57_67_I31") 
# Simulation 11.3
SimNew11_F57_67_I49_CA <- TIDY_SIM(dataset = SimNew11_F57_67_I49, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I49", Median_Name = "Median_SimNew11_F57_67_I49")

# Simulation 12.1
SimNew12_F86_50_I13_CA <- TIDY_SIM(dataset = SimNew12_F86_50_I13, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I13", Median_Name = "Median_SimNew12_F86_50_I13") 
# Simulation 12.2
SimNew12_F86_50_I31_CA <- TIDY_SIM(dataset = SimNew12_F86_50_I31, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I31", Median_Name = "Median_SimNew12_F86_50_I31") 
# Simulation 12.3
SimNew12_F86_50_I49_CA <- TIDY_SIM(dataset = SimNew12_F86_50_I49, CountryVec = CA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I49", Median_Name = "Median_SimNew12_F86_50_I49")


#### Merge extracted simulation data with PCR-based glm data -------------------
#put all data frames into list
GLEAM_GLM_CA_merged <- list(Predicted_merged_BA1_CA_selected, 
                            SimNew1_F1_I13_CA,
                            SimNew1_F1_I31_CA,
                            SimNew1_F1_I49_CA,
                            SimNew2_F1_5_I13_CA,
                            SimNew2_F1_5_I31_CA,
                            SimNew2_F1_5_I49_CA,
                            SimNew3_F2_25_I13_CA,
                            SimNew3_F2_25_I31_CA,
                            SimNew3_F2_25_I49_CA,
                            SimNew4_F3_38_I13_CA,
                            SimNew4_F3_38_I31_CA,
                            SimNew4_F3_38_I49_CA,
                            SimNew5_F5_06_I13_CA,
                            SimNew5_F5_06_I31_CA,
                            SimNew5_F5_06_I49_CA,
                            SimNew6_F7_59_I13_CA,
                            SimNew6_F7_59_I31_CA,
                            SimNew6_F7_59_I49_CA,
                            SimNew7_F11_39_I13_CA,
                            SimNew7_F11_39_I31_CA,
                            SimNew7_F11_39_I49_CA,
                            SimNew8_F17_09_I13_CA,
                            SimNew8_F17_09_I31_CA,
                            SimNew8_F17_09_I49_CA,
                            SimNew9_F25_63_I13_CA,
                            SimNew9_F25_63_I31_CA,
                            SimNew9_F25_63_I49_CA,
                            SimNew10_F38_44_I13_CA,
                            SimNew10_F38_44_I31_CA,
                            SimNew10_F38_44_I49_CA,
                            SimNew11_F57_67_I13_CA,
                            SimNew11_F57_67_I31_CA,
                            SimNew11_F57_67_I49_CA,
                            SimNew12_F86_50_I13_CA,
                            SimNew12_F86_50_I31_CA,
                            SimNew12_F86_50_I49_CA)


#### Merge all data frames in list ---------------------------------------------
GLEAM_GLM_CA_merged <- GLEAM_GLM_CA_merged %>% 
  reduce(full_join, by = "Days_post_June") %>%
  #select only data before February
  filter(Days_post_June < 245)





#### EASTERN AFRICA ============================================================
#### Apply function to extract simulation data ---------------------------------
# Simulation 1.1
SimNew1_F1_I13_EA <- TIDY_SIM(dataset = SimNew1_F1_I13, CountryVec = EA, 
                              Mean_Name = "Mean_SimNew1_F1_I13", Median_Name = "Median_SimNew1_F1_I13") 
# Simulation 1.2
SimNew1_F1_I31_EA <- TIDY_SIM(dataset = SimNew1_F1_I31, CountryVec = EA, 
                              Mean_Name = "Mean_SimNew1_F1_I31", Median_Name = "Median_SimNew1_F1_I31") 
# Simulation 1.3
SimNew1_F1_I49_EA <- TIDY_SIM(dataset = SimNew1_F1_I49, CountryVec = EA, 
                              Mean_Name = "Mean_SimNew1_F1_I49", Median_Name = "Median_SimNew1_F1_I49") 

# Simulation 2.1
SimNew2_F1_5_I13_EA <- TIDY_SIM(dataset = SimNew2_F1_5_I13, CountryVec = EA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I13", Median_Name = "Median_SimNew2_F1_5_I13") 
# Simulation 2.2
SimNew2_F1_5_I31_EA <- TIDY_SIM(dataset = SimNew2_F1_5_I31, CountryVec = EA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I31", Median_Name = "Median_SimNew2_F1_5_I31") 
# Simulation 2.3
SimNew2_F1_5_I49_EA <- TIDY_SIM(dataset = SimNew2_F1_5_I49, CountryVec = EA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I49", Median_Name = "Median_SimNew2_F1_5_I49")

# Simulation 3.1
SimNew3_F2_25_I13_EA <- TIDY_SIM(dataset = SimNew3_F2_25_I13, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I13", Median_Name = "Median_SimNew3_F2_25_I13") 
# Simulation 3.2
SimNew3_F2_25_I31_EA <- TIDY_SIM(dataset = SimNew3_F2_25_I31, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I31", Median_Name = "Median_SimNew3_F2_25_I31") 
# Simulation 3.3
SimNew3_F2_25_I49_EA <- TIDY_SIM(dataset = SimNew3_F2_25_I49, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I49", Median_Name = "Median_SimNew3_F2_25_I49")

# Simulation 4.1
SimNew4_F3_38_I13_EA <- TIDY_SIM(dataset = SimNew4_F3_38_I13, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I13", Median_Name = "Median_SimNew4_F3_38_I13") 
# Simulation 4.2
SimNew4_F3_38_I31_EA <- TIDY_SIM(dataset = SimNew4_F3_38_I31, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I31", Median_Name = "Median_SimNew4_F3_38_I31") 
# Simulation 4.3
SimNew4_F3_38_I49_EA <- TIDY_SIM(dataset = SimNew4_F3_38_I49, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I49", Median_Name = "Median_SimNew4_F3_38_I49")

# Simulation 5.1
SimNew5_F5_06_I13_EA <- TIDY_SIM(dataset = SimNew5_F5_06_I13, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I13", Median_Name = "Median_SimNew5_F5_06_I13") 
# Simulation 5.2
SimNew5_F5_06_I31_EA <- TIDY_SIM(dataset = SimNew5_F5_06_I31, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I31", Median_Name = "Median_SimNew5_F5_06_I31") 
# Simulation 5.3
SimNew5_F5_06_I49_EA <- TIDY_SIM(dataset = SimNew5_F5_06_I49, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I49", Median_Name = "Median_SimNew5_F5_06_I49")

# Simulation 6.1
SimNew6_F7_59_I13_EA <- TIDY_SIM(dataset = SimNew6_F7_59_I13, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I13", Median_Name = "Median_SimNew6_F7_59_I13") 
# Simulation 6.2
SimNew6_F7_59_I31_EA <- TIDY_SIM(dataset = SimNew6_F7_59_I31, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I31", Median_Name = "Median_SimNew6_F7_59_I31") 
# Simulation 6.3
SimNew6_F7_59_I49_EA <- TIDY_SIM(dataset = SimNew6_F7_59_I49, CountryVec = EA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I49", Median_Name = "Median_SimNew6_F7_59_I49")

# Simulation 7.1
SimNew7_F11_39_I13_EA <- TIDY_SIM(dataset = SimNew7_F11_39_I13, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I13", Median_Name = "Median_SimNew7_F11_39_I13") 
# Simulation 7.2
SimNew7_F11_39_I31_EA <- TIDY_SIM(dataset = SimNew7_F11_39_I31, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I31", Median_Name = "Median_SimNew7_F11_39_I31") 
# Simulation 7.3
SimNew7_F11_39_I49_EA <- TIDY_SIM(dataset = SimNew7_F11_39_I49, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I49", Median_Name = "Median_SimNew7_F11_39_I49")

# Simulation 8.1
SimNew8_F17_09_I13_EA <- TIDY_SIM(dataset = SimNew8_F17_09_I13, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I13", Median_Name = "Median_SimNew8_F17_09_I13") 
# Simulation 8.2
SimNew8_F17_09_I31_EA <- TIDY_SIM(dataset = SimNew8_F17_09_I31, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I31", Median_Name = "Median_SimNew8_F17_09_I31") 
# Simulation 8.3
SimNew8_F17_09_I49_EA <- TIDY_SIM(dataset = SimNew8_F17_09_I49, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I49", Median_Name = "Median_SimNew8_F17_09_I49")

# Simulation 9.1
SimNew9_F25_63_I13_EA <- TIDY_SIM(dataset = SimNew9_F25_63_I13, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I13", Median_Name = "Median_SimNew9_F25_63_I13") 
# Simulation 9.2
SimNew9_F25_63_I31_EA <- TIDY_SIM(dataset = SimNew9_F25_63_I31, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I31", Median_Name = "Median_SimNew9_F25_63_I31") 
# Simulation 9.3
SimNew9_F25_63_I49_EA <- TIDY_SIM(dataset = SimNew9_F25_63_I49, CountryVec = EA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I49", Median_Name = "Median_SimNew9_F25_63_I49")

# Simulation 10.1
SimNew10_F38_44_I13_EA <- TIDY_SIM(dataset = SimNew10_F38_44_I13, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I13", Median_Name = "Median_SimNew10_F38_44_I13") 
# Simulation 10.2
SimNew10_F38_44_I31_EA <- TIDY_SIM(dataset = SimNew10_F38_44_I31, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I31", Median_Name = "Median_SimNew10_F38_44_I31") 
# Simulation 10.3
SimNew10_F38_44_I49_EA <- TIDY_SIM(dataset = SimNew10_F38_44_I49, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I49", Median_Name = "Median_SimNew10_F38_44_I49")

# Simulation 11.1
SimNew11_F57_67_I13_EA <- TIDY_SIM(dataset = SimNew11_F57_67_I13, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I13", Median_Name = "Median_SimNew11_F57_67_I13") 
# Simulation 11.2
SimNew11_F57_67_I31_EA <- TIDY_SIM(dataset = SimNew11_F57_67_I31, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I31", Median_Name = "Median_SimNew11_F57_67_I31") 
# Simulation 11.3
SimNew11_F57_67_I49_EA <- TIDY_SIM(dataset = SimNew11_F57_67_I49, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I49", Median_Name = "Median_SimNew11_F57_67_I49")

# Simulation 12.1
SimNew12_F86_50_I13_EA <- TIDY_SIM(dataset = SimNew12_F86_50_I13, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I13", Median_Name = "Median_SimNew12_F86_50_I13") 
# Simulation 12.2
SimNew12_F86_50_I31_EA <- TIDY_SIM(dataset = SimNew12_F86_50_I31, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I31", Median_Name = "Median_SimNew12_F86_50_I31") 
# Simulation 12.3
SimNew12_F86_50_I49_EA <- TIDY_SIM(dataset = SimNew12_F86_50_I49, CountryVec = EA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I49", Median_Name = "Median_SimNew12_F86_50_I49")


#### Merge extracted simulation data with PCR-based glm data -------------------
#put all data frames into list
GLEAM_GLM_EA_merged <- list(Predicted_merged_BA1_EA_selected, 
                            SimNew1_F1_I13_EA,
                            SimNew1_F1_I31_EA,
                            SimNew1_F1_I49_EA,
                            SimNew2_F1_5_I13_EA,
                            SimNew2_F1_5_I31_EA,
                            SimNew2_F1_5_I49_EA,
                            SimNew3_F2_25_I13_EA,
                            SimNew3_F2_25_I31_EA,
                            SimNew3_F2_25_I49_EA,
                            SimNew4_F3_38_I13_EA,
                            SimNew4_F3_38_I31_EA,
                            SimNew4_F3_38_I49_EA,
                            SimNew5_F5_06_I13_EA,
                            SimNew5_F5_06_I31_EA,
                            SimNew5_F5_06_I49_EA,
                            SimNew6_F7_59_I13_EA,
                            SimNew6_F7_59_I31_EA,
                            SimNew6_F7_59_I49_EA,
                            SimNew7_F11_39_I13_EA,
                            SimNew7_F11_39_I31_EA,
                            SimNew7_F11_39_I49_EA,
                            SimNew8_F17_09_I13_EA,
                            SimNew8_F17_09_I31_EA,
                            SimNew8_F17_09_I49_EA,
                            SimNew9_F25_63_I13_EA,
                            SimNew9_F25_63_I31_EA,
                            SimNew9_F25_63_I49_EA,
                            SimNew10_F38_44_I13_EA,
                            SimNew10_F38_44_I31_EA,
                            SimNew10_F38_44_I49_EA,
                            SimNew11_F57_67_I13_EA,
                            SimNew11_F57_67_I31_EA,
                            SimNew11_F57_67_I49_EA,
                            SimNew12_F86_50_I13_EA,
                            SimNew12_F86_50_I31_EA,
                            SimNew12_F86_50_I49_EA)


#### Merge all data frames in list ---------------------------------------------
GLEAM_GLM_EA_merged <- GLEAM_GLM_EA_merged %>% 
  reduce(full_join, by = "Days_post_June") %>%
  #select only data before February
  filter(Days_post_June < 245)





#### WESTERN AFRICA ============================================================
#### Apply function to extract simulation data ---------------------------------
# Simulation 1.1
SimNew1_F1_I13_SA <- TIDY_SIM(dataset = SimNew1_F1_I13, CountryVec = SA, 
                              Mean_Name = "Mean_SimNew1_F1_I13", Median_Name = "Median_SimNew1_F1_I13") 
# Simulation 1.2
SimNew1_F1_I31_SA <- TIDY_SIM(dataset = SimNew1_F1_I31, CountryVec = SA, 
                              Mean_Name = "Mean_SimNew1_F1_I31", Median_Name = "Median_SimNew1_F1_I31") 
# Simulation 1.3
SimNew1_F1_I49_SA <- TIDY_SIM(dataset = SimNew1_F1_I49, CountryVec = SA, 
                              Mean_Name = "Mean_SimNew1_F1_I49", Median_Name = "Median_SimNew1_F1_I49") 

# Simulation 2.1
SimNew2_F1_5_I13_SA <- TIDY_SIM(dataset = SimNew2_F1_5_I13, CountryVec = SA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I13", Median_Name = "Median_SimNew2_F1_5_I13") 
# Simulation 2.2
SimNew2_F1_5_I31_SA <- TIDY_SIM(dataset = SimNew2_F1_5_I31, CountryVec = SA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I31", Median_Name = "Median_SimNew2_F1_5_I31") 
# Simulation 2.3
SimNew2_F1_5_I49_SA <- TIDY_SIM(dataset = SimNew2_F1_5_I49, CountryVec = SA, 
                                Mean_Name = "Mean_SimNew2_F1_5_I49", Median_Name = "Median_SimNew2_F1_5_I49")

# Simulation 3.1
SimNew3_F2_25_I13_SA <- TIDY_SIM(dataset = SimNew3_F2_25_I13, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I13", Median_Name = "Median_SimNew3_F2_25_I13") 
# Simulation 3.2
SimNew3_F2_25_I31_SA <- TIDY_SIM(dataset = SimNew3_F2_25_I31, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I31", Median_Name = "Median_SimNew3_F2_25_I31") 
# Simulation 3.3
SimNew3_F2_25_I49_SA <- TIDY_SIM(dataset = SimNew3_F2_25_I49, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew3_F2_25_I49", Median_Name = "Median_SimNew3_F2_25_I49")

# Simulation 4.1
SimNew4_F3_38_I13_SA <- TIDY_SIM(dataset = SimNew4_F3_38_I13, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I13", Median_Name = "Median_SimNew4_F3_38_I13") 
# Simulation 4.2
SimNew4_F3_38_I31_SA <- TIDY_SIM(dataset = SimNew4_F3_38_I31, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I31", Median_Name = "Median_SimNew4_F3_38_I31") 
# Simulation 4.3
SimNew4_F3_38_I49_SA <- TIDY_SIM(dataset = SimNew4_F3_38_I49, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew4_F3_38_I49", Median_Name = "Median_SimNew4_F3_38_I49")

# Simulation 5.1
SimNew5_F5_06_I13_SA <- TIDY_SIM(dataset = SimNew5_F5_06_I13, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I13", Median_Name = "Median_SimNew5_F5_06_I13") 
# Simulation 5.2
SimNew5_F5_06_I31_SA <- TIDY_SIM(dataset = SimNew5_F5_06_I31, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I31", Median_Name = "Median_SimNew5_F5_06_I31") 
# Simulation 5.3
SimNew5_F5_06_I49_SA <- TIDY_SIM(dataset = SimNew5_F5_06_I49, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew5_F5_06_I49", Median_Name = "Median_SimNew5_F5_06_I49")

# Simulation 6.1
SimNew6_F7_59_I13_SA <- TIDY_SIM(dataset = SimNew6_F7_59_I13, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I13", Median_Name = "Median_SimNew6_F7_59_I13") 
# Simulation 6.2
SimNew6_F7_59_I31_SA <- TIDY_SIM(dataset = SimNew6_F7_59_I31, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I31", Median_Name = "Median_SimNew6_F7_59_I31") 
# Simulation 6.3
SimNew6_F7_59_I49_SA <- TIDY_SIM(dataset = SimNew6_F7_59_I49, CountryVec = SA, 
                                 Mean_Name = "Mean_SimNew6_F7_59_I49", Median_Name = "Median_SimNew6_F7_59_I49")

# Simulation 7.1
SimNew7_F11_39_I13_SA <- TIDY_SIM(dataset = SimNew7_F11_39_I13, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I13", Median_Name = "Median_SimNew7_F11_39_I13") 
# Simulation 7.2
SimNew7_F11_39_I31_SA <- TIDY_SIM(dataset = SimNew7_F11_39_I31, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I31", Median_Name = "Median_SimNew7_F11_39_I31") 
# Simulation 7.3
SimNew7_F11_39_I49_SA <- TIDY_SIM(dataset = SimNew7_F11_39_I49, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew7_F11_39_I49", Median_Name = "Median_SimNew7_F11_39_I49")

# Simulation 8.1
SimNew8_F17_09_I13_SA <- TIDY_SIM(dataset = SimNew8_F17_09_I13, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I13", Median_Name = "Median_SimNew8_F17_09_I13") 
# Simulation 8.2
SimNew8_F17_09_I31_SA <- TIDY_SIM(dataset = SimNew8_F17_09_I31, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I31", Median_Name = "Median_SimNew8_F17_09_I31") 
# Simulation 8.3
SimNew8_F17_09_I49_SA <- TIDY_SIM(dataset = SimNew8_F17_09_I49, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew8_F17_09_I49", Median_Name = "Median_SimNew8_F17_09_I49")

# Simulation 9.1
SimNew9_F25_63_I13_SA <- TIDY_SIM(dataset = SimNew9_F25_63_I13, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I13", Median_Name = "Median_SimNew9_F25_63_I13") 
# Simulation 9.2
SimNew9_F25_63_I31_SA <- TIDY_SIM(dataset = SimNew9_F25_63_I31, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I31", Median_Name = "Median_SimNew9_F25_63_I31") 
# Simulation 9.3
SimNew9_F25_63_I49_SA <- TIDY_SIM(dataset = SimNew9_F25_63_I49, CountryVec = SA, 
                                  Mean_Name = "Mean_SimNew9_F25_63_I49", Median_Name = "Median_SimNew9_F25_63_I49")

# Simulation 10.1
SimNew10_F38_44_I13_SA <- TIDY_SIM(dataset = SimNew10_F38_44_I13, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I13", Median_Name = "Median_SimNew10_F38_44_I13") 
# Simulation 10.2
SimNew10_F38_44_I31_SA <- TIDY_SIM(dataset = SimNew10_F38_44_I31, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I31", Median_Name = "Median_SimNew10_F38_44_I31") 
# Simulation 10.3
SimNew10_F38_44_I49_SA <- TIDY_SIM(dataset = SimNew10_F38_44_I49, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew10_F38_44_I49", Median_Name = "Median_SimNew10_F38_44_I49")

# Simulation 11.1
SimNew11_F57_67_I13_SA <- TIDY_SIM(dataset = SimNew11_F57_67_I13, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I13", Median_Name = "Median_SimNew11_F57_67_I13") 
# Simulation 11.2
SimNew11_F57_67_I31_SA <- TIDY_SIM(dataset = SimNew11_F57_67_I31, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I31", Median_Name = "Median_SimNew11_F57_67_I31") 
# Simulation 11.3
SimNew11_F57_67_I49_SA <- TIDY_SIM(dataset = SimNew11_F57_67_I49, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew11_F57_67_I49", Median_Name = "Median_SimNew11_F57_67_I49")

# Simulation 12.1
SimNew12_F86_50_I13_SA <- TIDY_SIM(dataset = SimNew12_F86_50_I13, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I13", Median_Name = "Median_SimNew12_F86_50_I13") 
# Simulation 12.2
SimNew12_F86_50_I31_SA <- TIDY_SIM(dataset = SimNew12_F86_50_I31, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I31", Median_Name = "Median_SimNew12_F86_50_I31") 
# Simulation 12.3
SimNew12_F86_50_I49_SA <- TIDY_SIM(dataset = SimNew12_F86_50_I49, CountryVec = SA, 
                                   Mean_Name = "Mean_SimNew12_F86_50_I49", Median_Name = "Median_SimNew12_F86_50_I49")


#### Merge extracted simulation data with PCR-based glm data -------------------
#put all data frames into list
GLEAM_GLM_SA_merged <- list(Predicted_merged_BA1_SA_selected, 
                            SimNew1_F1_I13_SA,
                            SimNew1_F1_I31_SA,
                            SimNew1_F1_I49_SA,
                            SimNew2_F1_5_I13_SA,
                            SimNew2_F1_5_I31_SA,
                            SimNew2_F1_5_I49_SA,
                            SimNew3_F2_25_I13_SA,
                            SimNew3_F2_25_I31_SA,
                            SimNew3_F2_25_I49_SA,
                            SimNew4_F3_38_I13_SA,
                            SimNew4_F3_38_I31_SA,
                            SimNew4_F3_38_I49_SA,
                            SimNew5_F5_06_I13_SA,
                            SimNew5_F5_06_I31_SA,
                            SimNew5_F5_06_I49_SA,
                            SimNew6_F7_59_I13_SA,
                            SimNew6_F7_59_I31_SA,
                            SimNew6_F7_59_I49_SA,
                            SimNew7_F11_39_I13_SA,
                            SimNew7_F11_39_I31_SA,
                            SimNew7_F11_39_I49_SA,
                            SimNew8_F17_09_I13_SA,
                            SimNew8_F17_09_I31_SA,
                            SimNew8_F17_09_I49_SA,
                            SimNew9_F25_63_I13_SA,
                            SimNew9_F25_63_I31_SA,
                            SimNew9_F25_63_I49_SA,
                            SimNew10_F38_44_I13_SA,
                            SimNew10_F38_44_I31_SA,
                            SimNew10_F38_44_I49_SA,
                            SimNew11_F57_67_I13_SA,
                            SimNew11_F57_67_I31_SA,
                            SimNew11_F57_67_I49_SA,
                            SimNew12_F86_50_I13_SA,
                            SimNew12_F86_50_I31_SA,
                            SimNew12_F86_50_I49_SA)


#### Merge all data frames in list ---------------------------------------------
GLEAM_GLM_SA_merged <- GLEAM_GLM_SA_merged %>% 
  reduce(full_join, by = "Days_post_June") %>%
  #select only data before February
  filter(Days_post_June < 245)





#### WESTERN AFRICA ============================================================
#### Apply function to extract simulation data ---------------------------------
# Simulation 1.1
SimNew1_F1_I13_NoA <- TIDY_SIM(dataset = SimNew1_F1_I13, CountryVec = NoA, 
                               Mean_Name = "Mean_SimNew1_F1_I13", Median_Name = "Median_SimNew1_F1_I13") 
# Simulation 1.2
SimNew1_F1_I31_NoA <- TIDY_SIM(dataset = SimNew1_F1_I31, CountryVec = NoA, 
                               Mean_Name = "Mean_SimNew1_F1_I31", Median_Name = "Median_SimNew1_F1_I31") 
# Simulation 1.3
SimNew1_F1_I49_NoA <- TIDY_SIM(dataset = SimNew1_F1_I49, CountryVec = NoA, 
                               Mean_Name = "Mean_SimNew1_F1_I49", Median_Name = "Median_SimNew1_F1_I49") 

# Simulation 2.1
SimNew2_F1_5_I13_NoA <- TIDY_SIM(dataset = SimNew2_F1_5_I13, CountryVec = NoA, 
                                 Mean_Name = "Mean_SimNew2_F1_5_I13", Median_Name = "Median_SimNew2_F1_5_I13") 
# Simulation 2.2
SimNew2_F1_5_I31_NoA <- TIDY_SIM(dataset = SimNew2_F1_5_I31, CountryVec = NoA, 
                                 Mean_Name = "Mean_SimNew2_F1_5_I31", Median_Name = "Median_SimNew2_F1_5_I31") 
# Simulation 2.3
SimNew2_F1_5_I49_NoA <- TIDY_SIM(dataset = SimNew2_F1_5_I49, CountryVec = NoA, 
                                 Mean_Name = "Mean_SimNew2_F1_5_I49", Median_Name = "Median_SimNew2_F1_5_I49")

# Simulation 3.1
SimNew3_F2_25_I13_NoA <- TIDY_SIM(dataset = SimNew3_F2_25_I13, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew3_F2_25_I13", Median_Name = "Median_SimNew3_F2_25_I13") 
# Simulation 3.2
SimNew3_F2_25_I31_NoA <- TIDY_SIM(dataset = SimNew3_F2_25_I31, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew3_F2_25_I31", Median_Name = "Median_SimNew3_F2_25_I31") 
# Simulation 3.3
SimNew3_F2_25_I49_NoA <- TIDY_SIM(dataset = SimNew3_F2_25_I49, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew3_F2_25_I49", Median_Name = "Median_SimNew3_F2_25_I49")

# Simulation 4.1
SimNew4_F3_38_I13_NoA <- TIDY_SIM(dataset = SimNew4_F3_38_I13, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew4_F3_38_I13", Median_Name = "Median_SimNew4_F3_38_I13") 
# Simulation 4.2
SimNew4_F3_38_I31_NoA <- TIDY_SIM(dataset = SimNew4_F3_38_I31, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew4_F3_38_I31", Median_Name = "Median_SimNew4_F3_38_I31") 
# Simulation 4.3
SimNew4_F3_38_I49_NoA <- TIDY_SIM(dataset = SimNew4_F3_38_I49, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew4_F3_38_I49", Median_Name = "Median_SimNew4_F3_38_I49")

# Simulation 5.1
SimNew5_F5_06_I13_NoA <- TIDY_SIM(dataset = SimNew5_F5_06_I13, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew5_F5_06_I13", Median_Name = "Median_SimNew5_F5_06_I13") 
# Simulation 5.2
SimNew5_F5_06_I31_NoA <- TIDY_SIM(dataset = SimNew5_F5_06_I31, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew5_F5_06_I31", Median_Name = "Median_SimNew5_F5_06_I31") 
# Simulation 5.3
SimNew5_F5_06_I49_NoA <- TIDY_SIM(dataset = SimNew5_F5_06_I49, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew5_F5_06_I49", Median_Name = "Median_SimNew5_F5_06_I49")

# Simulation 6.1
SimNew6_F7_59_I13_NoA <- TIDY_SIM(dataset = SimNew6_F7_59_I13, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew6_F7_59_I13", Median_Name = "Median_SimNew6_F7_59_I13") 
# Simulation 6.2
SimNew6_F7_59_I31_NoA <- TIDY_SIM(dataset = SimNew6_F7_59_I31, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew6_F7_59_I31", Median_Name = "Median_SimNew6_F7_59_I31") 
# Simulation 6.3
SimNew6_F7_59_I49_NoA <- TIDY_SIM(dataset = SimNew6_F7_59_I49, CountryVec = NoA, 
                                  Mean_Name = "Mean_SimNew6_F7_59_I49", Median_Name = "Median_SimNew6_F7_59_I49")

# Simulation 7.1
SimNew7_F11_39_I13_NoA <- TIDY_SIM(dataset = SimNew7_F11_39_I13, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew7_F11_39_I13", Median_Name = "Median_SimNew7_F11_39_I13") 
# Simulation 7.2
SimNew7_F11_39_I31_NoA <- TIDY_SIM(dataset = SimNew7_F11_39_I31, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew7_F11_39_I31", Median_Name = "Median_SimNew7_F11_39_I31") 
# Simulation 7.3
SimNew7_F11_39_I49_NoA <- TIDY_SIM(dataset = SimNew7_F11_39_I49, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew7_F11_39_I49", Median_Name = "Median_SimNew7_F11_39_I49")

# Simulation 8.1
SimNew8_F17_09_I13_NoA <- TIDY_SIM(dataset = SimNew8_F17_09_I13, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew8_F17_09_I13", Median_Name = "Median_SimNew8_F17_09_I13") 
# Simulation 8.2
SimNew8_F17_09_I31_NoA <- TIDY_SIM(dataset = SimNew8_F17_09_I31, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew8_F17_09_I31", Median_Name = "Median_SimNew8_F17_09_I31") 
# Simulation 8.3
SimNew8_F17_09_I49_NoA <- TIDY_SIM(dataset = SimNew8_F17_09_I49, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew8_F17_09_I49", Median_Name = "Median_SimNew8_F17_09_I49")

# Simulation 9.1
SimNew9_F25_63_I13_NoA <- TIDY_SIM(dataset = SimNew9_F25_63_I13, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew9_F25_63_I13", Median_Name = "Median_SimNew9_F25_63_I13") 
# Simulation 9.2
SimNew9_F25_63_I31_NoA <- TIDY_SIM(dataset = SimNew9_F25_63_I31, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew9_F25_63_I31", Median_Name = "Median_SimNew9_F25_63_I31") 
# Simulation 9.3
SimNew9_F25_63_I49_NoA <- TIDY_SIM(dataset = SimNew9_F25_63_I49, CountryVec = NoA, 
                                   Mean_Name = "Mean_SimNew9_F25_63_I49", Median_Name = "Median_SimNew9_F25_63_I49")

# Simulation 10.1
SimNew10_F38_44_I13_NoA <- TIDY_SIM(dataset = SimNew10_F38_44_I13, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew10_F38_44_I13", Median_Name = "Median_SimNew10_F38_44_I13") 
# Simulation 10.2
SimNew10_F38_44_I31_NoA <- TIDY_SIM(dataset = SimNew10_F38_44_I31, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew10_F38_44_I31", Median_Name = "Median_SimNew10_F38_44_I31") 
# Simulation 10.3
SimNew10_F38_44_I49_NoA <- TIDY_SIM(dataset = SimNew10_F38_44_I49, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew10_F38_44_I49", Median_Name = "Median_SimNew10_F38_44_I49")

# Simulation 11.1
SimNew11_F57_67_I13_NoA <- TIDY_SIM(dataset = SimNew11_F57_67_I13, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew11_F57_67_I13", Median_Name = "Median_SimNew11_F57_67_I13") 
# Simulation 11.2
SimNew11_F57_67_I31_NoA <- TIDY_SIM(dataset = SimNew11_F57_67_I31, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew11_F57_67_I31", Median_Name = "Median_SimNew11_F57_67_I31") 
# Simulation 11.3
SimNew11_F57_67_I49_NoA <- TIDY_SIM(dataset = SimNew11_F57_67_I49, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew11_F57_67_I49", Median_Name = "Median_SimNew11_F57_67_I49")

# Simulation 12.1
SimNew12_F86_50_I13_NoA <- TIDY_SIM(dataset = SimNew12_F86_50_I13, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew12_F86_50_I13", Median_Name = "Median_SimNew12_F86_50_I13") 
# Simulation 12.2
SimNew12_F86_50_I31_NoA <- TIDY_SIM(dataset = SimNew12_F86_50_I31, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew12_F86_50_I31", Median_Name = "Median_SimNew12_F86_50_I31") 
# Simulation 12.3
SimNew12_F86_50_I49_NoA <- TIDY_SIM(dataset = SimNew12_F86_50_I49, CountryVec = NoA, 
                                    Mean_Name = "Mean_SimNew12_F86_50_I49", Median_Name = "Median_SimNew12_F86_50_I49")


#### Merge extracted simulation data with PCR-based glm data -------------------
#put all data frames into list
GLEAM_GLM_NoA_merged <- list(Predicted_merged_BA1_NoA_selected, 
                             SimNew1_F1_I13_NoA,
                             SimNew1_F1_I31_NoA,
                             SimNew1_F1_I49_NoA,
                             SimNew2_F1_5_I13_NoA,
                             SimNew2_F1_5_I31_NoA,
                             SimNew2_F1_5_I49_NoA,
                             SimNew3_F2_25_I13_NoA,
                             SimNew3_F2_25_I31_NoA,
                             SimNew3_F2_25_I49_NoA,
                             SimNew4_F3_38_I13_NoA,
                             SimNew4_F3_38_I31_NoA,
                             SimNew4_F3_38_I49_NoA,
                             SimNew5_F5_06_I13_NoA,
                             SimNew5_F5_06_I31_NoA,
                             SimNew5_F5_06_I49_NoA,
                             SimNew6_F7_59_I13_NoA,
                             SimNew6_F7_59_I31_NoA,
                             SimNew6_F7_59_I49_NoA,
                             SimNew7_F11_39_I13_NoA,
                             SimNew7_F11_39_I31_NoA,
                             SimNew7_F11_39_I49_NoA,
                             SimNew8_F17_09_I13_NoA,
                             SimNew8_F17_09_I31_NoA,
                             SimNew8_F17_09_I49_NoA,
                             SimNew9_F25_63_I13_NoA,
                             SimNew9_F25_63_I31_NoA,
                             SimNew9_F25_63_I49_NoA,
                             SimNew10_F38_44_I13_NoA,
                             SimNew10_F38_44_I31_NoA,
                             SimNew10_F38_44_I49_NoA,
                             SimNew11_F57_67_I13_NoA,
                             SimNew11_F57_67_I31_NoA,
                             SimNew11_F57_67_I49_NoA,
                             SimNew12_F86_50_I13_NoA,
                             SimNew12_F86_50_I31_NoA,
                             SimNew12_F86_50_I49_NoA)


#### Merge all data frames in list ---------------------------------------------
GLEAM_GLM_NoA_merged <- GLEAM_GLM_NoA_merged %>% 
  reduce(full_join, by = "Days_post_June") %>%
  #select only data before February
  filter(Days_post_June < 245)




#### STATISTICALLY COMPARE SIMULATIONS WITH PCR DATA ###########################

#### DEFINE FUNCTION ===========================================================
SquareSum <- function(dataset, Region, DpJ) {
  ## Filter countries
  data_loaded <- dataset %>%
    # Filter based on days post June
    filter(Days_post_June <= DpJ) %>%
    # Calculate squared sifference between simulation and PCR based estimate
    mutate(Mean_SimNew1_F1_I13_Diff = (Mean_SimNew1_F1_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew1_F1_I31_Diff = (Mean_SimNew1_F1_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew1_F1_I49_Diff = (Mean_SimNew1_F1_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew2_F1_5_I13_Diff = (Mean_SimNew2_F1_5_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew2_F1_5_I31_Diff = (Mean_SimNew2_F1_5_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew2_F1_5_I49_Diff = (Mean_SimNew2_F1_5_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew3_F2_25_I13_Diff = (Mean_SimNew3_F2_25_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew3_F2_25_I31_Diff = (Mean_SimNew3_F2_25_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew3_F2_25_I49_Diff = (Mean_SimNew3_F2_25_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew4_F3_38_I13_Diff = (Mean_SimNew4_F3_38_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew4_F3_38_I31_Diff = (Mean_SimNew4_F3_38_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew4_F3_38_I49_Diff = (Mean_SimNew4_F3_38_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew5_F5_06_I13_Diff = (Mean_SimNew5_F5_06_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew5_F5_06_I31_Diff = (Mean_SimNew5_F5_06_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew5_F5_06_I49_Diff = (Mean_SimNew5_F5_06_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew6_F7_59_I13_Diff = (Mean_SimNew6_F7_59_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew6_F7_59_I31_Diff = (Mean_SimNew6_F7_59_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew6_F7_59_I49_Diff = (Mean_SimNew6_F7_59_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew7_F11_39_I13_Diff = (Mean_SimNew7_F11_39_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew7_F11_39_I31_Diff = (Mean_SimNew7_F11_39_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew7_F11_39_I49_Diff = (Mean_SimNew7_F11_39_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew8_F17_09_I13_Diff = (Mean_SimNew8_F17_09_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew8_F17_09_I31_Diff = (Mean_SimNew8_F17_09_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew8_F17_09_I49_Diff = (Mean_SimNew8_F17_09_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew9_F25_63_I13_Diff = (Mean_SimNew9_F25_63_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew9_F25_63_I31_Diff = (Mean_SimNew9_F25_63_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew9_F25_63_I49_Diff = (Mean_SimNew9_F25_63_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew10_F38_44_I13_Diff = (Mean_SimNew10_F38_44_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew10_F38_44_I31_Diff = (Mean_SimNew10_F38_44_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew10_F38_44_I49_Diff = (Mean_SimNew10_F38_44_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew11_F57_67_I13_Diff = (Mean_SimNew11_F57_67_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew11_F57_67_I31_Diff = (Mean_SimNew11_F57_67_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew11_F57_67_I49_Diff = (Mean_SimNew11_F57_67_I49 - Mean_BA1_Cases)^2,
           Mean_SimNew12_F86_50_I13_Diff = (Mean_SimNew12_F86_50_I13 - Mean_BA1_Cases)^2,
           Mean_SimNew12_F86_50_I31_Diff = (Mean_SimNew12_F86_50_I31 - Mean_BA1_Cases)^2,
           Mean_SimNew12_F86_50_I49_Diff = (Mean_SimNew12_F86_50_I49 - Mean_BA1_Cases)^2) %>%
    # Summarize squared difference over time
    summarize(SimNew1_F1_I13 = sum(Mean_SimNew1_F1_I13_Diff),
              SimNew1_F1_I31 = sum(Mean_SimNew1_F1_I31_Diff),
              SimNew1_F1_I49 = sum(Mean_SimNew1_F1_I49_Diff),
              SimNew2_F1_5_I13 = sum(Mean_SimNew2_F1_5_I13_Diff),
              SimNew2_F1_5_I31 = sum(Mean_SimNew2_F1_5_I31_Diff),
              SimNew2_F1_5_I49 = sum(Mean_SimNew2_F1_5_I49_Diff),
              SimNew3_F2_25_I13 = sum(Mean_SimNew3_F2_25_I13_Diff),
              SimNew3_F2_25_I31 = sum(Mean_SimNew3_F2_25_I31_Diff),
              SimNew3_F2_25_I49 = sum(Mean_SimNew3_F2_25_I49_Diff),
              SimNew4_F3_38_I13 = sum(Mean_SimNew4_F3_38_I13_Diff),
              SimNew4_F3_38_I31 = sum(Mean_SimNew4_F3_38_I31_Diff),
              SimNew4_F3_38_I49 = sum(Mean_SimNew4_F3_38_I49_Diff),
              SimNew5_F5_06_I13 = sum(Mean_SimNew5_F5_06_I13_Diff),
              SimNew5_F5_06_I31 = sum(Mean_SimNew5_F5_06_I31_Diff),
              SimNew5_F5_06_I49 = sum(Mean_SimNew5_F5_06_I49_Diff),
              SimNew6_F7_59_I13 = sum(Mean_SimNew6_F7_59_I13_Diff),
              SimNew6_F7_59_I31 = sum(Mean_SimNew6_F7_59_I31_Diff),
              SimNew6_F7_59_I49 = sum(Mean_SimNew6_F7_59_I49_Diff),
              SimNew7_F11_39_I13 = sum(Mean_SimNew7_F11_39_I13_Diff),
              SimNew7_F11_39_I31 = sum(Mean_SimNew7_F11_39_I31_Diff),
              SimNew7_F11_39_I49 = sum(Mean_SimNew7_F11_39_I49_Diff),
              SimNew8_F17_09_I13 = sum(Mean_SimNew8_F17_09_I13_Diff),
              SimNew8_F17_09_I31 = sum(Mean_SimNew8_F17_09_I31_Diff),
              SimNew8_F17_09_I49 = sum(Mean_SimNew8_F17_09_I49_Diff),
              SimNew9_F25_63_I13 = sum(Mean_SimNew9_F25_63_I13_Diff),
              SimNew9_F25_63_I31 = sum(Mean_SimNew9_F25_63_I31_Diff),
              SimNew9_F25_63_I49 = sum(Mean_SimNew9_F25_63_I49_Diff),
              SimNew10_F38_44_I13 = sum(Mean_SimNew10_F38_44_I13_Diff),
              SimNew10_F38_44_I31 = sum(Mean_SimNew10_F38_44_I31_Diff),
              SimNew10_F38_44_I49 = sum(Mean_SimNew10_F38_44_I49_Diff),
              SimNew11_F57_67_I13 = sum(Mean_SimNew11_F57_67_I13_Diff),
              SimNew11_F57_67_I31 = sum(Mean_SimNew11_F57_67_I31_Diff),
              SimNew11_F57_67_I49 = sum(Mean_SimNew11_F57_67_I49_Diff),
              SimNew12_F86_50_I13 = sum(Mean_SimNew12_F86_50_I13_Diff),
              SimNew12_F86_50_I31 = sum(Mean_SimNew12_F86_50_I31_Diff),
              SimNew12_F86_50_I49 = sum(Mean_SimNew12_F86_50_I49_Diff)) %>%
    mutate(Region = as.character(Region))
  ## Return dataframe
  return(data_loaded)
}


#### APPLY FUNCTION ============================================================
# Western Africa
GLEAM_GLM_WA_merged_Distance2 <- SquareSum(dataset = GLEAM_GLM_WA_merged,
                                           Region = "WesternAfrica",
                                           #DpJ = 228,
                                           DpJ = 213)

# Central Africa
GLEAM_GLM_CA_merged_Distance2 <- SquareSum(dataset = GLEAM_GLM_CA_merged,
                                           Region = "CentralAfrica",
                                           #DpJ = 228,
                                           DpJ = 213)

# Eastern Africa
GLEAM_GLM_EA_merged_Distance2 <- SquareSum(dataset = GLEAM_GLM_EA_merged,
                                           Region = "EasternAfrica",
                                           #DpJ = 228,
                                           DpJ = 213)

# Northern Africa
GLEAM_GLM_NoA_merged_Distance2 <- SquareSum(dataset = GLEAM_GLM_NoA_merged,
                                            Region = "NorthernAfrica",
                                            #DpJ = 228,
                                            DpJ = 213)

# Southern Africa
GLEAM_GLM_SA_merged_Distance2 <- SquareSum(dataset = GLEAM_GLM_SA_merged,
                                           Region = "SouthernAfrica",
                                           #DpJ = 228,
                                           DpJ = 213)



#### JOIN SQUARED SUM DATA =====================================================
GLEAM_GLM_merged_Distance2_joined <- rbind(GLEAM_GLM_WA_merged_Distance2,
                                           GLEAM_GLM_CA_merged_Distance2,
                                           GLEAM_GLM_EA_merged_Distance2,
                                           GLEAM_GLM_NoA_merged_Distance2,
                                           GLEAM_GLM_SA_merged_Distance2) %>%
  # Re-arrange dataset
  pivot_longer(cols = 1:36, names_to = "Simulation", values_to = "Squared_Sum") %>%
  pivot_wider(values_from = "Squared_Sum", names_from = "Region") %>%
  # Calculate squared sum ratios
  mutate(WA_ratio = WesternAfrica/min(WesternAfrica),
         CA_ratio = CentralAfrica/min(CentralAfrica),
         EA_ratio = EasternAfrica/min(EasternAfrica),
         NoA_ratio = NorthernAfrica/min(NorthernAfrica),
         SA_ratio = SouthernAfrica/min(SouthernAfrica))
# Summarize mean of squared sum ratios
# Sum
GLEAM_GLM_merged_Distance2_joined$Ratios_summed <- round(rowSums(GLEAM_GLM_merged_Distance2_joined[7:11]), digits = 4)
# Mean
GLEAM_GLM_merged_Distance2_joined$Ratios_mean <- round(rowMeans(GLEAM_GLM_merged_Distance2_joined[7:11]), digits = 4)
# Median
GLEAM_GLM_merged_Distance2_joined$Ratios_median <- round(rowMedians(as.matrix(GLEAM_GLM_merged_Distance2_joined[,c(7:11)])), digits = 4)



#### Export results as table ---------------------------------------------------
# Set wd
setwd("SELECT PATH")

# Export table 
write.csv(GLEAM_GLM_merged_Distance2_joined,
          "GLEAMviz2_vs_PCRestimate.csv")




#### PREPARE PLOTS #############################################################

#### HEAT MAP ==================================================================
#### Manipulate data -----------------------------------------------------------
Data_HeatMap <- GLEAM_GLM_merged_Distance2_joined %>%
  # select data
  select(Simulation, Ratios_mean) %>%
  # Unify names
  mutate(Simulation = ifelse(Simulation == "SimNew1_F1_I13", "SimNew1_F1_0_I13",
                      ifelse(Simulation == "SimNew1_F1_I31", "SimNew1_F1_0_I31",
                      ifelse(Simulation == "SimNew1_F1_I49", "SimNew1_F1_0_I49",
                             Simulation)))) %>%
  # split simualtion name into data collumn to access protective immunity
  separate(Simulation, c("A", "B", "C", "D")) %>%
  # Clean dataset
  select(Simulation = A, Immunity = D, Ratios_mean) %>%
  # Make variables numeric
  mutate(Simulation = as.numeric(str_remove(Simulation, "SimNew")),
         Immunity = as.numeric(str_remove(Immunity, "I"))) %>%
  # logarithmize values
  mutate(Ratios_mean_log = log10(Ratios_mean))
  

#### Plot data -----------------------------------------------------------------
# Heatmap 
plot.HeatMap <- ggplot(data = Data_HeatMap, aes(x = Immunity, y = Simulation, fill = Ratios_mean_log)) + 
  scale_y_continuous(expand = c(0,0),
                     breaks = c(1,2,3,4,5,6,7,8,9,10,11,12)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(13,31,49)) +
  scale_fill_gradient2(low = "#C70E46",
                      mid = "#FCE0B8",
                      high = "#4E7096",
                      midpoint = 0.85) +
  labs(x = "Protective Immunity (%)",
       fill = "Log10 of mean 
squared sum
ratio") +
  geom_tile() +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, face="bold"),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.ticks.length=unit(.12, "cm"),
        axis.line = element_line(size = 0.2),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 6.5))

# see plot
plot.HeatMap




#### REGIONAL CURVES ===========================================================
#### Tranform data for plotting ------------------------------------------------
# Transform days post June to Date
GLEAM_GLM_WA_merged <- GLEAM_GLM_WA_merged %>%
  mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Days_post_June))

GLEAM_GLM_CA_merged <- GLEAM_GLM_CA_merged %>%
  mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Days_post_June))

GLEAM_GLM_EA_merged <- GLEAM_GLM_EA_merged %>%
  mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Days_post_June))

GLEAM_GLM_NoA_merged <- GLEAM_GLM_NoA_merged %>%
  mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Days_post_June))

GLEAM_GLM_SA_merged <- GLEAM_GLM_SA_merged %>%
  mutate(Date = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + 
           lubridate::days(Days_post_June))


#### Define function to plot data ----------------------------------------------

PLOT_SIMS_V2<- function(dataset, Y_Limit, Region, bestFit) {
  # Plot values for quality control
  P1 <- ggplot(data = dataset) +
    geom_line(aes(x = Date, y = Mean_SimNew1_F1_I13), size = 0.3, colour = ifelse(bestFit == "SimNew1_F1_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew1_F1_I31), size = 0.3, colour = ifelse(bestFit == "SimNew1_F1_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew1_F1_I49), size = 0.3, colour = ifelse(bestFit == "SimNew1_F1_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew2_F1_5_I13), size = 0.3, colour = ifelse(bestFit == "SimNew2_F1_5_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew2_F1_5_I31), size = 0.3, colour = ifelse(bestFit == "SimNew2_F1_5_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew2_F1_5_I49), size = 0.3, colour = ifelse(bestFit == "SimNew2_F1_5_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew3_F2_25_I13), size = 0.3, colour = ifelse(bestFit == "SimNew3_F2_25_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew3_F2_25_I31), size = 0.3, colour = ifelse(bestFit == "SimNew3_F2_25_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew3_F2_25_I49), size = 0.3, colour = ifelse(bestFit == "SimNew3_F2_25_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew4_F3_38_I13), size = 0.3, colour = ifelse(bestFit == "SimNew4_F3_38_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew4_F3_38_I31), size = 0.3, colour = ifelse(bestFit == "SimNew4_F3_38_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew4_F3_38_I49), size = 0.3, colour = ifelse(bestFit == "SimNew4_F3_38_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew5_F5_06_I13), size = 0.3, colour = ifelse(bestFit == "SimNew5_F5_06_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew5_F5_06_I31), size = 0.3, colour = ifelse(bestFit == "SimNew5_F5_06_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew5_F5_06_I49), size = 0.3, colour = ifelse(bestFit == "SimNew5_F5_06_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew6_F7_59_I13), size = 0.3, colour = ifelse(bestFit == "SimNew6_F7_59_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew6_F7_59_I31), size = 0.3, colour = ifelse(bestFit == "SimNew6_F7_59_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew6_F7_59_I49), size = 0.3, colour = ifelse(bestFit == "SimNew6_F7_59_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew7_F11_39_I13), size = 0.3, colour = ifelse(bestFit == "SimNew7_F11_39_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew7_F11_39_I31), size = 0.3, colour = ifelse(bestFit == "SimNew7_F11_39_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew7_F11_39_I49), size = 0.3, colour = ifelse(bestFit == "SimNew7_F11_39_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew8_F17_09_I13), size = 0.3, colour = ifelse(bestFit == "SimNew8_F17_09_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew8_F17_09_I31), size = 0.3, colour = ifelse(bestFit == "SimNew8_F17_09_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew8_F17_09_I49), size = 0.3, colour = ifelse(bestFit == "SimNew8_F17_09_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew9_F25_63_I13), size = 0.3, colour = ifelse(bestFit == "SimNew9_F25_63_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew9_F25_63_I31), size = 0.3, colour = ifelse(bestFit == "SimNew9_F25_63_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew9_F25_63_I49), size = 0.3, colour = ifelse(bestFit == "SimNew9_F25_63_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew10_F38_44_I13), size = 0.3, colour = ifelse(bestFit == "SimNew10_F38_44_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew10_F38_44_I31), size = 0.3, colour = ifelse(bestFit == "SimNew10_F38_44_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew10_F38_44_I49), size = 0.3, colour = ifelse(bestFit == "SimNew10_F38_44_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew11_F57_67_I13), size = 0.3, colour = ifelse(bestFit == "SimNew11_F57_67_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew11_F57_67_I31), size = 0.3, colour = ifelse(bestFit == "SimNew11_F57_67_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew11_F57_67_I49), size = 0.3, colour = ifelse(bestFit == "SimNew11_F57_67_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew12_F86_50_I13), size = 0.3, colour = ifelse(bestFit == "SimNew12_F86_50_I13", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew12_F86_50_I31), size = 0.3, colour = ifelse(bestFit == "SimNew12_F86_50_I31", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_SimNew12_F86_50_I49), size = 0.3, colour = ifelse(bestFit == "SimNew12_F86_50_I49", "#CA4A6E", "#CBCBCB")) +
    geom_line(aes(x = Date, y = Mean_BA1_Cases), size = 0.3, colour = "orange") +
    geom_point(aes(x = Date, y = Mean_SimNew1_F1_I13), size = 0.5, colour = ifelse(bestFit == "SimNew1_F1_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew1_F1_I31), size = 0.5, colour = ifelse(bestFit == "SimNew1_F1_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew1_F1_I49), size = 0.5, colour = ifelse(bestFit == "SimNew1_F1_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew2_F1_5_I13), size = 0.5, colour = ifelse(bestFit == "SimNew2_F1_5_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew2_F1_5_I31), size = 0.5, colour = ifelse(bestFit == "SimNew2_F1_5_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew2_F1_5_I49), size = 0.5, colour = ifelse(bestFit == "SimNew2_F1_5_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew3_F2_25_I13), size = 0.5, colour = ifelse(bestFit == "SimNew3_F2_25_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew3_F2_25_I31), size = 0.5, colour = ifelse(bestFit == "SimNew3_F2_25_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew3_F2_25_I49), size = 0.5, colour = ifelse(bestFit == "SimNew3_F2_25_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew4_F3_38_I13), size = 0.5, colour = ifelse(bestFit == "SimNew4_F3_38_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew4_F3_38_I31), size = 0.5, colour = ifelse(bestFit == "SimNew4_F3_38_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew4_F3_38_I49), size = 0.5, colour = ifelse(bestFit == "SimNew4_F3_38_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew5_F5_06_I13), size = 0.5, colour = ifelse(bestFit == "SimNew5_F5_06_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew5_F5_06_I31), size = 0.5, colour = ifelse(bestFit == "SimNew5_F5_06_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew5_F5_06_I49), size = 0.5, colour = ifelse(bestFit == "SimNew5_F5_06_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew6_F7_59_I13), size = 0.5, colour = ifelse(bestFit == "SimNew6_F7_59_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew6_F7_59_I31), size = 0.5, colour = ifelse(bestFit == "SimNew6_F7_59_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew6_F7_59_I49), size = 0.5, colour = ifelse(bestFit == "SimNew6_F7_59_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew7_F11_39_I13), size = 0.5, colour = ifelse(bestFit == "SimNew7_F11_39_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew7_F11_39_I31), size = 0.5, colour = ifelse(bestFit == "SimNew7_F11_39_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew7_F11_39_I49), size = 0.5, colour = ifelse(bestFit == "SimNew7_F11_39_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew8_F17_09_I13), size = 0.5, colour = ifelse(bestFit == "SimNew8_F17_09_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew8_F17_09_I31), size = 0.5, colour = ifelse(bestFit == "SimNew8_F17_09_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew8_F17_09_I49), size = 0.5, colour = ifelse(bestFit == "SimNew8_F17_09_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew9_F25_63_I13), size = 0.5, colour = ifelse(bestFit == "SimNew9_F25_63_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew9_F25_63_I31), size = 0.5, colour = ifelse(bestFit == "SimNew9_F25_63_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew9_F25_63_I49), size = 0.5, colour = ifelse(bestFit == "SimNew9_F25_63_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew10_F38_44_I13), size = 0.5, colour = ifelse(bestFit == "SimNew10_F38_44_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew10_F38_44_I31), size = 0.5, colour = ifelse(bestFit == "SimNew10_F38_44_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew10_F38_44_I49), size = 0.5, colour = ifelse(bestFit == "SimNew10_F38_44_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew11_F57_67_I13), size = 0.5, colour = ifelse(bestFit == "SimNew11_F57_67_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew11_F57_67_I31), size = 0.5, colour = ifelse(bestFit == "SimNew11_F57_67_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew11_F57_67_I49), size = 0.5, colour = ifelse(bestFit == "SimNew11_F57_67_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew12_F86_50_I13), size = 0.5, colour = ifelse(bestFit == "SimNew12_F86_50_I13", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew12_F86_50_I31), size = 0.5, colour = ifelse(bestFit == "SimNew12_F86_50_I31", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_SimNew12_F86_50_I49), size = 0.5, colour = ifelse(bestFit == "SimNew12_F86_50_I49", "#CA4A6E", "#CBCBCB")) +
    geom_point(aes(x = Date, y = Mean_BA1_Cases), colour = "orange", size = 0.5) +
    #scale_x_continuous(limits = c(170, 228)) +
    scale_x_datetime(#expand=c(0,0),
      breaks= as.POSIXct(c("2021-11-15", "2021-12-01", "2021-12-15", "2021-12-31", "2022-01-15")),
      #breaks = seq(as.Date("2021-11-11"), as.Date("2022-01-15"), by="2 weeks"),
      #date_breaks= "14 days", 
      #date_minor_breaks = "2 weeks", 
      #date_labels = "%d/%m", 
      limits = as.POSIXct(c("2021-11-11", "2021-12-31"))) +
    scale_y_continuous(limits = c(0, Y_Limit)) +
    theme_classic() +
    ggtitle(Region) +
    theme(plot.title=element_text(margin=margin(t=30, b=-30))) +
    theme(plot.title=element_text(hjust = 0.05, vjust = 0.5)) +
    labs(y = "BA.1 cases/100,000 inhabitants", x = "Date") +
  theme(plot.title = element_text(size = 7, face="bold"),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, face="bold"),
        axis.ticks = element_line(colour="black", size=0.2),
        axis.ticks.length=unit(.12, "cm"),
        axis.line=element_line(size=0.2),
        axis.text.x=element_text(angle = 35, hjust = 1)
  )
        
        
  ## Return dataframe
  return(P1)
}




#### Apply function to plot data -----------------------------------------------

# Western Africa
WA_Plot <- PLOT_SIMS_V2(dataset = GLEAM_GLM_WA_merged, Y_Limit = 120, 
                        Region = "Western Africa", bestFit = "SimNew12_F86_50_I13")

# Central Africa
CA_Plot <- PLOT_SIMS_V2(dataset = GLEAM_GLM_CA_merged, Y_Limit = 150, 
                        Region = "Central Africa", bestFit = "SimNew11_F57_67_I31")

# Eastern Africa
EA_Plot <- PLOT_SIMS_V2(dataset = GLEAM_GLM_EA_merged, Y_Limit = 300, 
                        Region = "Eastern Africa", bestFit = "SimNew12_F86_50_I13")

# Northern Africa
NoA_Plot <- PLOT_SIMS_V2(dataset = GLEAM_GLM_NoA_merged, Y_Limit = 100, 
                         Region = "Northern Africa", bestFit = "SimNew12_F86_50_I13")

# Southern Africa
SA_Plot <- PLOT_SIMS_V2(dataset = GLEAM_GLM_SA_merged, Y_Limit = 3000, 
                        Region = "Southern Africa", bestFit = "SimNew11_F57_67_I13")




#### MERGE PLOTS ===============================================================
## Arrange plots
GLEAMviz_Plot_arranged <- grid.arrange(
  CA_Plot, EA_Plot, NoA_Plot,
  SA_Plot, WA_Plot, plot.HeatMap,
  nrow = 2
)

## Plot arrange dplots
GLEAMviz_Plot_arranged



#### Export merged plots -------------------------------------------------------
# Set wd
setwd("SELECT PATH")

ggsave("GLEAMviz_vs_PCRpredicted_arranged.pdf", GLEAMviz_Plot_arranged, width = 8, height = 5.5,
       limitsize = FALSE)



