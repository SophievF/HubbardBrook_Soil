## Hubbard Brook archived soil samples project ##
## Data preparation ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)

# Load HBEF data (dataset 2)
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-11-07.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes)) 

head(HBEF_data)

# Rename mineral horizon names
HBEF_data$Horizon <- replace(HBEF_data$Horizon, HBEF_data$Horizon == "Mineral_0_10",
                             "Mineral")

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral"),
                            ordered = TRUE)

# Check mean values for 14C and SOC (w/o roots)
HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  group_by(Horizon) %>% 
  reframe(mean_14C = mean(Delta14C, na.rm = TRUE),
          mean_C = mean(`Measured_%_C`, na.rm = TRUE))

# Load earlier 14C data (dataset 1)
LitterData <- read_csv("./Data/LitterData_Driscoll.csv")

# fm * exp(lambda * (-obs_date_y + 1950)) - 1) * 1000
lambda <- 0.0001209681

# Modify litter data to match with other dataset
litter_data <- LitterData %>% 
  filter(Watershed == 6) %>% 
  mutate(Delta14C = (F14C * exp(lambda * (-Year + 1950)) -1) * 1000) %>% 
  mutate(Elevation = case_when(
    Plot < 70 ~ "High",
    Plot > 152 ~ "Low",
    TRUE ~ "Mid"
  )) %>% 
  #match Horizon names
  mutate(Horizon = case_when(
    Horizon == "Oie" ~ "Oi/Oe",
    Horizon == "Oa" ~ "Oa/A"
  )) 

# Count number of samples per year (field replicates)
litter_data %>%
  count(Year) 

# Load additional soil data from HBEF (bulk density)
soil_info <- read_csv("./Data/MassChemistryOrganicHorizonMineralSoil_WS6_1976_present/HubbardBrook_ForestFloor_SoilMass_W6.csv") %>%
  dplyr::select(-Watershed)

litter_all <- soil_info %>%
  #calculate SOC stocks (from kg/m2 to g/m2) and from OM to C
  mutate(SOC_g_m2 = OM_OM * 1000 * 0.58) %>%
  mutate(Horizon = case_when(
    Horizon == "Oie" ~ "Oi/Oe",
    Horizon == "Oa" ~ "Oa/A"
  )) %>%
  dplyr::select(Year, Plot, Horizon, SOC_g_m2) %>%
  right_join(litter_data)

#gap-fill missing 1969 data from 1978 and closest Plot
litter_all[31,4] <- 1235.4 #Plot 2
litter_all[32,4] <- 986 #Plot 129
litter_all[33,4] <- 2383.8 #Plot 114
litter_all[34,4] <- 777.2 #Plot 2
litter_all[35,4] <- 904.8 #Plot 129
litter_all[36,4] <- 1508 #Plot 114

# Summarize and merge data by horizon; remove roots for now
HBEF_all <- HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  filter(Year != 2020) %>% 
  mutate(DataSource = "Dataset1") %>% 
  mutate(SOC_g_m2 = (`Measured_%_C` * mean_BD_g_cm3 * mean_thick_cm) * 100) %>% 
  dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, DataSource, Elevation) %>% 
  rbind(litter_all %>% 
          dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, Elevation) %>% 
          mutate(DataSource = "Dataset2")) %>% 
  dplyr::select(Year:Delta14C, SOC_g_m2, DataSource, Elevation)

write_csv(HBEF_all,
          file = paste0("./Data/HBEF_ModelingDatasets_all_", Sys.Date(), ".csv"))
