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




#### LOAD DATA #################################################################
#### PCR DATA ------------------------------------------------------------------
## Set wd
setwd("SELECT PATH")

## Load dataset
data <- read.csv2("dataset.csv",
                 header = TRUE)

# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"


#### POPULATION DATA -----------------------------------------------------------
pop <- read_excel("WPP2022_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx") 
# Define countries
Countries_Africa <- c("Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Cameroon", 
                      "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Kenya", 
                      "Madagascar", "Mali", "Morocco", "Mozambique", "Namibia", 
                      "Niger", "Congo", "Senegal", "South Africa", "Togo", "Uganda", 
                      "Zimbabwe")


## Mutate population data
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
  # filter(Inbound_Travel %in% c(9, 0, NA))




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
dom.Bot <- MODEL_DOMINANCE(dataset = data_tidy, CountryX = "Botswana")
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
Dominance_Modelling <- rbind(dom.Alg, dom.Ang, dom.Ben, dom.Bot, dom.BuF, dom.Cam, dom.Eth,
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
write.csv(Dominance_Modelling,
          "~/PATH\\BA1_Dominance_Overview.csv")  


#### IDENTIFY FIRST CASES AND FIRST EXPECTED CASES BASED ON DOMINANCE ----------
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
###########          "~/PATH\\First_Cases.csv")  




#### REMOVE BA.1-POSITIVE SAMPLES COLLECTED BEFORE FIRST EXPECTED CASES --------
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
   


#### GROUP COUNTRIES BY REGION -------------------------------------------------
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
  # remove Madagascar to only show mainland EA
  filter(Country != "Madagascar")
#### EXCLUDE TRAVELER
#data_grouped <- data_grouped %>% filter(Inbound_Travel %in% c(9, 0, NA))


# Add region collumn to not-filered data
data_grouped_raw <- data_tidy %>%
  mutate(Region = ifelse(Country %in% WA, "Western Africa", 
                  ifelse(Country %in% NoA, "Northern Africa",
                  ifelse(Country %in% CA, "Central Africa",
                  ifelse(Country %in% EA, "Eastern Africa",
                  ifelse(Country %in% SA, "Southern Africa", NA)))))) %>%
  # remove MAdagascar to only show mainland EA
  filter(Country != "Madagascar")


#### MODELL TAKEOVER ----------------------------------------------------------
# create dummy dataframe to highlight first reported cases
first_case_point <- data.frame(Collection_date2  = as.POSIXct(c("2021-11-11")),
                               Omicron = c(1.04))
travel_restrictions_point <- data.frame(Collection_date2  = as.POSIXct(c("2021-11-27")),
                               Omicron = c(1.04))

# Plot glms
ggplot() +
  # Add line for travel restrictions
  geom_vline(xintercept = as.POSIXct(c("2021-11-27")),
             linetype="dotted", 
             color = "grey", size=0.8) +
  # Add point for travel restrictions
  geom_point(data = travel_restrictions_point,
             aes(x = Collection_date2,
                 y = Omicron), 
             color = "grey", size=3) +
  # Add line for first reported case
  geom_vline(xintercept = as.POSIXct(c("2021-11-11")),
             linetype="dotted", 
             color = "black", size=0.8) +
  # Add point for first reported case
  geom_point(data = first_case_point,
             aes(x = Collection_date2,
               y = Omicron), 
             color = "black", size=3) +
  # Add GLMs
  geom_smooth(data = data_grouped,
              aes(x = Collection_date2, y = Omicron, colour = Region, fill = Region),
              method.args=list(family="binomial"),
              method = glm,
              fullrange = TRUE, se = TRUE) +
  theme_classic() +
  # Define y axis limits
  scale_y_continuous(limits = c(0, 1.05),
                     expand = c(0,0)) +
  # Define x axis limits
  expand_limits(x = as.POSIXct(c("2021-04-01", "2022-08-01"))) +
  scale_x_datetime(breaks = as.POSIXct(c("2021-04-01", 
                                         "2021-07-01",
                                         "2021-10-01",
                                         "2022-01-01",
                                         "2022-04-01",
                                         "2022-07-01"))) +
  # relabel axes
  labs(y = "Omicron/BA.1 fraction",
       x = "Date")

# Save plot
ggsave("Takeover_Curves.pdf", height = 4, width = 7)
#ggsave("Takeover_Curves_NoTraveler.pdf", height = 4, width = 7)



#### CALCULATE DAYS BETWEEN FIST cASE AND DOMINANCE ---------------------------
days_to_dominance <- merge(First_Omicron,
                           Dominance_Modelling,
                           by = "Country") %>%
  mutate(date_dominance = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Day_Difference)) %>%
  dplyr::select(-3) %>%
  # Calculate days between first cases and modelled dominance
  mutate(Days_to_dom = as.numeric(round(difftime(date_dominance, 
                                                 Collection_date, 
                                                 units= "days"),
                                        digits = 0))) %>%
  # Remove countries with negative day difference as early cases were abviously missed
  filter(Days_to_dom > 0)

# Calculate median days between BA.1 emergence and dominance
Med_Dom_Days <- as.data.frame.matrix(groupwiseMedian(
  Days_to_dom ~ 1,
  data = days_to_dominance,
  conf = 0.95,
  bca = FALSE,
  R = 5000,
  exact = TRUE,
  digits = 2)) %>%
  rename(Lower = 5, Upper = 6)
           

## Plot days to dominance
ggplot() +
  geom_boxplot(data = days_to_dominance %>% filter(Days_to_dom > 0), 
               aes(x = 1, y = Days_to_dom)) +
  geom_jitter(data = days_to_dominance, 
              aes(x = 1, y = Days_to_dom), 
              width = 0.4, height = 0, size = 3, colour = "#FFA500") + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125),
                    limits = c(0, 135)) +
  theme_classic()

# Save plots
ggsave("Time_do_Dom.pdf", height = 5, width = 4.2)
#ggsave("Time_do_Dom_Notraveller.pdf", height = 5, width = 4.2)






#### IDENTIFY DATE FOR OMICRON DOMINANCE PER AFRICAN REGION ####################
#### Define function to identify date of dominance -----------------------------
MODEL_DOMINANCE_REGION <- function(dataset, RegionY) {
  # Fit glm
  glmY = glm(Omicron ~ Days_post_June, data = dataset %>% filter(Region == RegionY), 
             family = binomial(link = "logit"))
  # Create empty dataset
  newdata = data.frame(Days_post_June = 0:320) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
  predicted <- add_ci(newdata, 
                      glmY, 
                      names = c("lwr", "upr"), alpha = 0.1) 
  predicted <- setDT(predicted, keep.rownames = TRUE)[]   
  pred_pred <- predicted %>% filter(pred >= 0.5) %>% filter(pred == min(pred)) %>%
    mutate(Date_pred = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    dplyr::select(6)
  pred_lwr <- predicted %>% filter(lwr >= 0.5) %>% filter(lwr == min(lwr))%>%
    mutate(Date_lwr = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    dplyr::select(6)
  pred_upr <- predicted %>% filter(upr >= 0.5) %>% filter(upr == min(upr))%>%
    mutate(Date_upr = as.POSIXct("2021-06-01", format = "%Y-%m-%d") + lubridate::days(Days_post_June)) %>%
    dplyr::select(6)
  predicted_dates <- bind_cols(Region = RegionY, pred_pred, pred_lwr, pred_upr)
  # Return dataframe
  return(predicted_dates)
}



#### Apply function to identify dominance per Region ---------------------------
## For filtered data
Dom_WA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                 Region = "Western Africa")
Dom_CA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                 Region = "Central Africa")
Dom_EA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                 Region = "Eastern Africa")
Dom_SA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                 Region = "Southern Africa")
Dom_NA <- MODEL_DOMINANCE_REGION(dataset = data_grouped, 
                                 Region = "Northern Africa")

## For not-filtered data
Dom_WA_raw <- MODEL_DOMINANCE_REGION(dataset = data_grouped_raw, 
                                 Region = "Western Africa")
Dom_CA_raw <- MODEL_DOMINANCE_REGION(dataset = data_grouped_raw, 
                                 Region = "Central Africa")
Dom_EA_raw <- MODEL_DOMINANCE_REGION(dataset = data_grouped_raw, 
                                 Region = "Eastern Africa")
Dom_SA_raw <- MODEL_DOMINANCE_REGION(dataset = data_grouped_raw, 
                                 Region = "Southern Africa")
Dom_NA_raw <- MODEL_DOMINANCE_REGION(dataset = data_grouped_raw, 
                                 Region = "Northern Africa")

## Merge Modelling results
## Filtered data
Dom_Merged <- bind_rows(Dom_CA, Dom_EA, Dom_NA, Dom_SA, Dom_WA) %>%
  # Calculate difference in days
  mutate(lwr_days = difftime(Date_lwr, Date_pred, units= "days"),
         upr_days = difftime(Date_upr, Date_pred, units= "days"))

## Not-filtered data
Dom_Merged_raw <- bind_rows(Dom_CA_raw, Dom_EA_raw, Dom_NA_raw, Dom_SA_raw, 
                            Dom_WA_raw) %>%
  # Calculate difference in days
  mutate(lwr_days_raw = difftime(Date_lwr, Date_pred, units= "days"),
         upr_days_raw = difftime(Date_upr, Date_pred, units= "days")) %>%
  rename(Date_pred_raw = Date_pred,
         Date_upr_raw = Date_upr,
         Date_lwr_raw = Date_lwr)

## Merge filtered and not-filtered data
Dom_Merged_all <- merge(Dom_Merged, Dom_Merged_raw, by = "Region")




#### Export table --------------------------------------------------------------
## Set wd
setwd("SELECT PATH")

## Export data
write.csv(Dom_Merged_all,
         "Dom_Merged_all.csv")

  
  
  
  
  
  


