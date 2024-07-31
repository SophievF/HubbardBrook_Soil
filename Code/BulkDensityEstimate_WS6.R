## Hubbard Brook archived soil samples project ##
## Calculate bulk density for WS6 ##
## Sophie von Fromm ##
## 2024-07-31 ##

# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-hbr.172.4

library(tidyverse)

## Grid and elevation data
# Load elevation data
grid_elevation <- read_csv("./Data/Grid_Elevation_WS6.csv") %>% 
  dplyr::select(Plot, Elevation_grp)

## Soil mass 
# Load soil mass data
soil_mass <- read_csv("./Data/MassChemistryOrganicHorizonMineralSoil_WS6_1976_present/HubbardBrook_ForestFloor_SoilMass_W6.csv")
head(soil_mass)

## Soil thickness
# Load soil thickness data
soil_thickness <- read_csv("./Data/MassChemistryOrganicHorizonMineralSoil_WS6_1976_present/HubbardBrook_ForestFloor_SiteInfo_W6.csv")
head(soil_thickness)

#Merge duplicate sample in year 2018 at plot 36
soil_2018_36 <- soil_thickness %>% 
  filter(Plot == 36.1| Plot == 36.2) %>% 
  dplyr::select(-c(E_Code:Sp_3)) %>% 
  mutate(across(FF_Thickness1:Core_4, ~na_if(., -9999.9))) %>% 
  summarise(across(Year:Core_4, mean))

soil_2018_36$Plot <- 36

# Replace -9999 with NA and join with elevation data
soil_thickness_NA <- soil_thickness %>% 
  # remove duplicate samples
  filter(Plot != 36.1, Plot != 36.2) %>% 
  mutate(across(FF_Thickness1:Core_4, ~na_if(., -9999.9))) %>% 
  dplyr::select(-c(E_Code:Sp_3)) %>%
  # add averaged values for duplicate sample
  rbind(soil_2018_36) %>% 
  #filter for year with measurements
  filter(Year > 1976) %>% 
  left_join(grid_elevation, by = "Plot")

soil_thickness_NA %>% 
  count(Elevation_grp)

soil_thick_sum <- soil_thickness_NA %>% 
  mutate(Oa_Thickness1 = FF_Thickness1 - Oie_Thickness1,
         Oa_Thickness2 = FF_Thickness2 - Oie_Thickness2,
         Oa_Thickness3 = FF_Thickness3 - Oie_Thickness3,
         Oa_Thickness4 = FF_Thickness4 - Oie_Thickness4,
         Oa_Thickness5 = FF_Thickness5 - Oie_Thickness5,
         Oa_Thickness6 = FF_Thickness6 - Oie_Thickness6,
         Oa_Thickness7 = FF_Thickness7 - Oie_Thickness7,
         Oa_Thickness8 = FF_Thickness8 - Oie_Thickness8,) %>% 
  rowwise() %>% 
  summarise(Year = Year,
            Plot = Plot,
            Elevation_grp = Elevation_grp,
            FF_thickness = mean(c(FF_Thickness1, FF_Thickness2, FF_Thickness3, 
                                  FF_Thickness4, FF_Thickness5, FF_Thickness6, 
                                  FF_Thickness7, FF_Thickness8), na.rm = TRUE),
            Oie_thickness = mean(c(Oie_Thickness1, Oie_Thickness2, Oie_Thickness3, 
                                   Oie_Thickness4, Oie_Thickness5, Oie_Thickness6, 
                                   Oie_Thickness7, Oie_Thickness8), na.rm = TRUE),
            Oa_thickness = mean(c(Oa_Thickness1, Oa_Thickness2, Oa_Thickness3, 
                                  Oa_Thickness4, Oa_Thickness5, Oa_Thickness6, 
                                  Oa_Thickness7, Oa_Thickness8), na.rm = TRUE),
            Min_thickness = mean(c(Core_1, Core_2, Core_3, Core_4), na.rm = TRUE)) 

soil_thick_sum %>% 
  group_by(Year, Elevation_grp) %>% 
  summarise(mean_FF = mean(FF_thickness, na.rm = TRUE),
            mean_Oie = mean(Oie_thickness, na.rm = TRUE),
            mean_Oa = mean(Oa_thickness, na.rm = TRUE),
            mean_min = mean(Min_thickness, na.rm = TRUE)) %>% 
  pivot_longer(!c(Year, Elevation_grp), names_to = "Horizon", values_to = "Thickness_cm") %>% 
  ggplot(aes(x = Year, y = Thickness_cm, color = Horizon)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Elevation_grp) +
  theme_bw()

soil_thick_sum %>% 
  pivot_longer(!c(Year, Elevation_grp, Plot), names_to = "Horizon", values_to = "Thickness_cm")

  
