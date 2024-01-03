#load packages
library(tidyverse)
library(readxl)
library(writexl)
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
### Set working directory
setwd("ENTER PATH")

### Load dataset
data <- read.csv2("dataset.csv",
                  header = TRUE)

# Change language settings to allow special letters for city names
Encoding(data$Location) <- "UTF-8"



#### TIDY DATA #################################################################
# Recode Variant status
data_recoded <- data %>%
  mutate(Variant = ifelse(is.na(SARS2) & is.na(Omicron) & is.na(Delta), "Unclear",
                   ifelse(SARS2 == 0 & is.na(Omicron) & is.na(Delta), "Negative",
                   ifelse(is.na(Omicron) | is.na(Delta), "Unclear",
                   ifelse(Omicron == 1 & Delta == 0, "Omicron",
                   ifelse(Omicron == 0 & Delta == 1, "Delta",
                   ifelse(SARS2 == 1 & Omicron == 0 & Delta == 0, "Other",
                   ifelse(SARS2 == 0 & Omicron == 0 & Delta == 0, "Negative",
                          "Unclear")))))))) %>%
  # Remove negative and unclear samples
  filter(Variant %in% c("Delta", "Omicron", "Other")) %>%
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



#### ANALYZE OVERALL SAMPLE DISTRIBUTION #######################################
#### COUNT SAMPLE ENTRIES PER MONTH --------------------------------------------
data_counted <- data_recoded %>%
  group_by(Country, Collection_date_month) %>%
  summarise(Samples = n())

#### PLOT SAMPLE DATA ----------------------------------------------------------
ggplot(data = data_counted, 
       aes(x = Collection_date_month, 
           y = factor(Country,
                      # Reorder y axis
                      levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = log10(Samples))) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#DAE8E8", 
                       mid = "#9788B3", 
                       high = "#5E0523", 
                       midpoint = 1.5,
                       labels = c("1", "10", "100", "1,000")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "Samples")

## Save plot
ggsave("Samples_overTime.pdf", height = 4, width = 7)



#### ANALYZE BA.1 FRACTION OVER TIME ###########################################
#### COUNT BA.1 FRACTION PER MONTH ---------------------------------------------
BA1_counted <- data_recoded %>%
  group_by(Country, Collection_date_month) %>%
  summarise(BA1_Fraction = mean(Omicron))

#### PLOT BA1 DATA -------------------------------------------------------------
ggplot(data = BA1_counted, 
       aes(x = Collection_date_month, 
           y = factor(Country,
                      # Reorder y axis
                      levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = BA1_Fraction)) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#DAE8E8", 
                       mid = "#9788B3", 
                       high = "#5E0523", 
                       midpoint = 0.5,
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "BA.1 fraction")

## Save plot
ggsave("BA1_overTime.pdf", height = 4, width = 7)



#### ANALYZE DELTA FRACTION OVER TIME ##########################################
#### COUNT DELTA FRACTION PER MONTH -------------------------------------------
Delta_counted <- data_recoded %>%
  group_by(Country, Collection_date_month) %>%
  summarise(Delta_Fraction = mean(Delta))


#### PLOT DELTA DATA ---------------------------------------------------------
ggplot(data = Delta_counted, 
       aes(x = Collection_date_month, 
           y = factor(Country,
                      # Reorder y axis
                      levels = rev(levels(factor(Country)))))) + 
  # Define heatmap
  geom_raster(aes(fill = Delta_Fraction)) + 
  theme_classic() +
  # Define colours
  scale_fill_gradient2(low ="#DAE8E8", 
                       mid = "#9788B3", 
                       high = "#5E0523", 
                       midpoint = 0.5,
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  xlab("Sampling period (by month)") + 
  ylab("Country")  +
  labs(fill = "Delta fraction")

## Save plot
ggsave("Delta_overTime.pdf", height = 4, width = 7)










