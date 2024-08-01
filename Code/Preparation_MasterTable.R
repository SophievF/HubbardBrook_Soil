## Hubbard Brook archived soil samples project ##
## Prepare master table for data anlysis ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)

# Load master table
HBEF_master <- read_csv("./Data/HBBMaster.csv")

head(HBEF_master)

# Load bulk density data
BD_avg_data <- read_csv("./Data/soil_BD_avg_WS6.csv")

# Rename Horizon names to match with master table
BD_avg_data$Horizon[BD_avg_data$Horizon == "min"] <- "Mineral_0_10"
BD_avg_data$Horizon[BD_avg_data$Horizon == "Oa"] <- "Oa/A"
BD_avg_data$Horizon[BD_avg_data$Horizon == "Oie"] <- "Oi/Oe"

# Check duplicate samples
HBEF_dup <- HBEF_master %>% 
  filter(grepl("HBB_012|HBB_040|HBB_056|HBB_078", ID))

HBEF_dup_avg <- HBEF_dup %>% 
  mutate(ID_group = case_when(
    grepl("HBB_012", ID) ~ "HBB_012",
    grepl("HBB_040", ID) ~ "HBB_040",
    grepl("HBB_056", ID) ~ "HBB_056",
    grepl("HBB_078", ID) ~ "HBB_078"
  )) %>% 
  group_by(ID_group) %>% 
  mutate(across(c(where(is.numeric)), mean, na.rm = TRUE)) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  dplyr::select(-ID_group)

# Merge master file with averaged duplicate samples
# Add bulk density data and calculate D14C
#fm * exp(lambda * (-obs_date_y + 1950)) - 1) * 1000
lambda <- 0.00012097

HBEF_data <- HBEF_master %>% 
  filter(!grepl("HBB_012|HBB_040|HBB_056|HBB_078", ID)) %>% 
  rbind(HBEF_dup_avg) %>% 
  tibble() %>% 
  mutate(Delta14C = (F14C * exp(lambda * (-Year + 1950)) -1) * 1000) %>% 
  left_join(BD_avg_data)

## Once field notes are complete, calculate C stocks (Mg ha-1) = SOC% * BD (g/cm3) * horizon thickness

# Save file
write_csv(x = HBEF_data, file = paste0("./Data/HBEF_data_all_", Sys.Date(), ".csv"))

  
  
  
  
  
  

  
