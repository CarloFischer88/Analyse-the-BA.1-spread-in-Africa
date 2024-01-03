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
library(matrixStats)
library(egg)
library(slider)



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
                  ifelse(Country %in% SA, "Southern Africa", NA)))))) #%>%
  # remove MAdagascar to only show mainland EA
  #filter(Country != "Madagascar")

#### EXCLUDE TRAVELER
#data_grouped <- data_grouped %>% filter(Inbound_Travel %in% c(9, 0, NA))



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




#### TEST DIFFERENT MODELS TO SELECT BEST FITTING MODEL -----------------------
#### Define function to identify best fitting model
AIC_TEST <- function(dataset, CountryX) {
  # filter data
  data_filtered <- dataset %>% filter(Country == CountryX) %>%
    filter(!is.na(Days_post_June)) %>%
    filter(Days_post_June < 251) %>%
    filter(!is.na(Omicron))
  # Calculate GLMs
  glm1 <- glm(Omicron ~ Days_post_June, data = data_filtered, 
                  family = binomial(link = "logit"))
  glm2 <- glm(Omicron ~ I(Days_post_June^2), data = data_filtered, 
                  family = binomial(link = "logit"))
  glm3 <- glm(Omicron ~ Days_post_June + I(Days_post_June^2), data = data_filtered, 
                  family = binomial(link = "logit"))
  glm4 <- glm(Omicron ~ I(Days_post_June^3), data = data_filtered, 
                  family = binomial(link = "logit"))
  glm5 <- glm(Omicron ~ Days_post_June + I(Days_post_June^3), data = data_filtered, 
              family = binomial(link = "logit"))
  glm6 <- glm(Omicron ~ I(Days_post_June/2), data = data_filtered, 
                  family = binomial(link = "logit"))
  glm7 <- glm(Omicron ~ I(Days_post_June*2), data = data_filtered, 
                  family = binomial(link = "logit"))
  # Extract AICs
  data_AIC <- as.data.frame(AIC(glm1,
                                glm2,
                                glm3,
                                glm4,
                                glm5,
                                glm6,
                                glm7))
  # Add model Name
  data_AIC$Model <- rownames(data_AIC)
  # Select data
  data_AIC <- data_AIC %>% dplyr::select(Model, AIC) %>%
    # Calculate ratio
    mutate(AIC_Ratio = round(as.numeric(AIC)/as.numeric(min(data_AIC$AIC)), digits = 4)) %>%
    # rename collumn
    rename(setNames("AIC", paste("AIC", CountryX, sep = "_"))) %>%
      rename(setNames("AIC_Ratio", paste("AIC_Ratio", CountryX, sep = "_")))

  # Return dataframe
  return(data_AIC)
}

## Aplly function
AIC_ALG <- AIC_TEST(data_grouped, CountryX = "Algeria")
AIC_ANG <- AIC_TEST(data_grouped, CountryX = "Angola")
AIC_BEN <- AIC_TEST(data_grouped, CountryX = "Benin")
AIC_BOT <- AIC_TEST(data_grouped, CountryX = "Botswana")
AIC_BuF <- AIC_TEST(data_grouped, CountryX = "Burkina Faso")
AIC_CAM <- AIC_TEST(data_grouped, CountryX = "Cameroon")
AIC_ETH <- AIC_TEST(data_grouped, CountryX = "Ethiopia")
AIC_GAM <- AIC_TEST(data_grouped, CountryX = "Gambia")
AIC_GHA <- AIC_TEST(data_grouped, CountryX = "Ghana")
AIC_GUI <- AIC_TEST(data_grouped, CountryX = "Guinea")
AIC_KEN <- AIC_TEST(data_grouped, CountryX = "Kenya")
AIC_MAD <- AIC_TEST(data_grouped, CountryX = "Madagascar")
AIC_MAL <- AIC_TEST(data_grouped, CountryX = "Mali")
AIC_MOR <- AIC_TEST(data_grouped, CountryX = "Morocco")
AIC_MOZ <- AIC_TEST(data_grouped, CountryX = "Mozambique")
AIC_NAM <- AIC_TEST(data_grouped, CountryX = "Namibia")
AIC_NIG <- AIC_TEST(data_grouped, CountryX = "Niger")
AIC_COG <- AIC_TEST(data_grouped, CountryX = "Republic of the Congo")
AIC_SEN <- AIC_TEST(data_grouped, CountryX = "Senegal")
AIC_SAf <- AIC_TEST(data_grouped, CountryX = "South Africa")
AIC_TOG <- AIC_TEST(data_grouped, CountryX = "Togo")
AIC_UGA <- AIC_TEST(data_grouped, CountryX = "Uganda")


## Merge AIC values
## Join predicted country data in one list to join them in a data frame
AIC_list <- list(AIC_ALG,
                 AIC_ANG,
                 AIC_BEN,
                 AIC_BOT,
                 AIC_BuF,
                 AIC_CAM,
                 AIC_ETH,
                 AIC_GAM,
                 AIC_GHA,
                 AIC_GUI,
                 AIC_KEN,
                 AIC_MAD,
                 AIC_MAL,
                 AIC_MOR,
                 AIC_MOZ,
                 AIC_NAM,
                 AIC_NIG,
                 AIC_COG,
                 AIC_SEN,
                 AIC_SAf,
                 AIC_TOG,
                 AIC_UGA)


## Merge all data frames of AICs
AICs_merged <- AIC_list %>% reduce(full_join, by='Model') 

## Sum all Ratios
AICs_summed <- AICs_merged %>%
  dplyr::select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45) 
# Sum
AICs_summed$AICs_summed <- round(rowSums(AICs_summed[2:23]), digits = 4)
# Mean
AICs_summed$AICs_mean <- round(rowMeans(AICs_summed[2:23]), digits = 4)
# Median
AICs_summed$AICs_median <- round(rowMedians(as.matrix(AICs_summed[,c(2:23)])), digits = 4)

## Add model descriptions
AICs_summed <- AICs_summed %>%
  mutate(Model_description = case_when(Model == "glm1" ~ "y ~ x",
                                       Model == "glm2" ~ "y ~ x^2",
                                       Model == "glm3" ~ "y ~ x + x^2",
                                       Model == "glm4" ~ "y ~ x^3",
                                       Model == "glm5" ~ "y ~ x + x^3",
                                       Model == "glm6" ~ "y ~ x/2",
                                       Model == "glm7" ~ "y ~ x*2"))

## Export AIC data
write.csv(AICs_summed,
          "~/PATH\\Countrywise_AICs.csv")




#### PLOT SELECTED GLMs TO VERIFY CORRECTNESS ---------------------------------
## Define function to modify nd plot data
MODEL_PLOT <- function(dataset, CountryX) {
  # filter data
  data_filtered <- dataset %>% 
    filter(Country == CountryX) %>% 
    filter(Days_post_June < 251)
  # Fit glm
  glmX = glm(#Omicron ~ Days_post_June + I(Days_post_June^2), 
    Omicron ~ Days_post_June, 
    data = data_filtered, 
    family = binomial(link = "logit"))
  # Create empty dataset
  newdata = data.frame(Days_post_June = 0:250) %>%
    mutate(Days_post_June = as.numeric(Days_post_June))
  # Predict values
  predicted <- as.data.frame(predict(glmX, newdata, #se.fit=TRUE,
                                     type = "response")) %>%
    rename(Value = 1)
  predicted <- setDT(predicted, keep.rownames = TRUE)[]   
  predicted <- predicted %>% 
    rename(Days_post_June = 1) %>%
    mutate(Days_post_June = as.numeric(Days_post_June),
           Value = round(Value, digits = 4))
  #filter(Days_post_June >= 50)
  # bin data
  data_bin <- data_filtered %>%
    mutate(bin = case_when(Days_post_June >= 30 & Days_post_June < 50 ~ 40,
                           Days_post_June >= 50 & Days_post_June < 70 ~ 60,
                           Days_post_June >= 70 & Days_post_June < 90 ~ 80,
                           Days_post_June >= 90 & Days_post_June < 110 ~ 100,
                           Days_post_June >= 110 & Days_post_June < 130 ~ 120,
                           Days_post_June >= 130 & Days_post_June < 150 ~ 140,
                           Days_post_June >= 150 & Days_post_June < 170 ~ 160,
                           Days_post_June >= 170 & Days_post_June < 190 ~ 180,
                           Days_post_June >= 190 & Days_post_June < 210 ~ 200,
                           Days_post_June >= 210 & Days_post_June < 230 ~ 220,
                           Days_post_June >= 230 & Days_post_June < 250 ~ 240)) %>%
    filter(!is.na(bin)) %>%
    group_by(bin) %>%
    summarize(mean_BA = mean(Omicron)) %>%
    mutate(mean_BA = round(mean_BA, digits = 3))
  data_rollingAverage <- data_grouped %>% 
    filter(Country == "Benin") %>%
    group_by(Days_post_June) %>%
    summarise(mean_Omicron = mean(Omicron, na.rm = TRUE)) %>%
    mutate(day20_mean = slide_dbl(mean_Omicron,                         # calculate on new_cases
                                  .f = ~mean(.x, na.rm = T),          # function is sum() with missing values removed
                                  .before = 10, 
                                  .after = 10))
  # rollng real mean
  data_BinnedMean_x <- data.frame()
  
  for(x in 0:250) {
    dataF <- as.data.frame(data_filtered %>% 
                             filter(Days_post_June > x-16,
                                    Days_post_June < x+16) %>%
                             summarize(Days_post_June = x,
                                       meanOmicron = mean(Omicron, na.rm = TRUE)))
    #print(dataF)
    #return(dataF)
    data_BinnedMean_x = rbind(data_BinnedMean_x, dataF)
  }
  
  data_BinnedMean_x <- data_BinnedMean_x %>%
    mutate(meanOmicron = as.numeric(meanOmicron))
  
  
  # Plot data
  plot1 <- ggplot() +
    geom_point(data = data_BinnedMean_x, 
               aes(x = Days_post_June, y = meanOmicron), colour = "black") +
    geom_point(data = predicted, 
               aes(x = Days_post_June, y = Value), colour = "red", alpha = 0.4) +
    #geom_point(data = data_bin, aes(x = bin, y = mean_BA), colour = "blue") +
    #geom_point(data = data_rollingAverage, aes(x = Days_post_June, y = day20_mean)) +
    ggtitle(CountryX) +
    theme(plot.title=element_text(margin=margin(t=30,b=-30))) +
    labs(y = "Omicron/BA.1 fraction", x = "Days post June")
  return(plot1)
}



## Apply function to plot GLMs
GLMplot_ALG <- MODEL_PLOT(dataset = data_grouped, CountryX = "Algeria")
GLMplot_ANG <- MODEL_PLOT(dataset = data_grouped, CountryX = "Angola")
GLMplot_BEN <- MODEL_PLOT(dataset = data_grouped, CountryX = "Benin")
GLMplot_BOT <- MODEL_PLOT(dataset = data_grouped, CountryX = "Botswana")
GLMplot_BUF <- MODEL_PLOT(dataset = data_grouped, CountryX = "Burkina Faso")
GLMplot_CAM <- MODEL_PLOT(dataset = data_grouped, CountryX = "Cameroon")
GLMplot_ETH <- MODEL_PLOT(dataset = data_grouped, CountryX = "Ethiopia")
GLMplot_GAM <- MODEL_PLOT(dataset = data_grouped, CountryX = "Gambia")
GLMplot_GHA <- MODEL_PLOT(dataset = data_grouped, CountryX = "Ghana")
GLMplot_GUI <- MODEL_PLOT(dataset = data_grouped, CountryX = "Guinea")
GLMplot_KEN <- MODEL_PLOT(dataset = data_grouped, CountryX = "Kenya")
GLMplot_MAD <- MODEL_PLOT(dataset = data_grouped, CountryX = "Madagascar")
GLMplot_MAL <- MODEL_PLOT(dataset = data_grouped, CountryX = "Mali")
GLMplot_MOR <- MODEL_PLOT(dataset = data_grouped, CountryX = "Morocco")
GLMplot_MOZ <- MODEL_PLOT(dataset = data_grouped, CountryX = "Mozambique")
GLMplot_NAM <- MODEL_PLOT(dataset = data_grouped, CountryX = "Namibia")
GLMplot_NIG <- MODEL_PLOT(dataset = data_grouped, CountryX = "Niger")
GLMplot_COG <- MODEL_PLOT(dataset = data_grouped, CountryX = "Republic of the Congo")
GLMplot_SEN <- MODEL_PLOT(dataset = data_grouped, CountryX = "Senegal")
GLMplot_SAf <- MODEL_PLOT(dataset = data_grouped, CountryX = "South Africa")
GLMplot_TOG <- MODEL_PLOT(dataset = data_grouped, CountryX = "Togo")
GLMplot_UGA <- MODEL_PLOT(dataset = data_grouped, CountryX = "Uganda")

# Arrange plots
GLM_Plot_arranged <- grid.arrange(
  GLMplot_ALG, GLMplot_ANG, GLMplot_BEN, GLMplot_BOT, GLMplot_BUF, GLMplot_CAM,
  GLMplot_ETH, GLMplot_GAM, GLMplot_GHA, GLMplot_GUI, GLMplot_KEN, GLMplot_MAD,
  GLMplot_MAL, GLMplot_MOR, GLMplot_MOZ, GLMplot_NAM, GLMplot_NIG, GLMplot_COG,
  GLMplot_SEN, GLMplot_SAf, GLMplot_TOG, GLMplot_UGA,
  nrow = 6,
  top = "GLM fit per country"
  )

# Save arranged plots
ggsave("GLM_Plot_arranged.pdf", GLM_Plot_arranged, width = 20, height = 25,
       limitsize = FALSE)

