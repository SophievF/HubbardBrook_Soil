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
soil_mass <- read_csv("./Data/MassChemistryOrganicHorizonMineralSoil_WS6_1976_present/HubbardBrook_ForestFloor_SoilMass_W6.csv") %>% 
  dplyr::select(-Watershed)
head(soil_mass)

#Merge duplicate sample in year 2013 at plot 156
soil_2013_156 <- soil_mass %>% 
  filter(grepl("2013-6-156", Sample_ID)) %>% 
  group_by(Horizon) %>% 
  summarise(across(OM_TM:OM_LOI, mean))

soil_2013_156$Plot <- 156
soil_2013_156$Year <- 2013

soil_mass_red <- soil_mass %>% 
  # remove duplicate samples
  filter(!grepl("2013-6-156", Sample_ID)) %>% 
  dplyr::select(-Sample_ID) %>% 
  # add averaged values for duplicate sample
  rbind(soil_2013_156) %>% 
  tibble()

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

## Merge thickness and mass data to calculate BD
soil_BD <- soil_thick_sum %>% 
  pivot_longer(!c(Year, Elevation_grp, Plot), names_to = "Horizon", values_to = "Thickness_cm") %>% 
  mutate(Thickness_cm = round(Thickness_cm, 2)) %>% 
  # set thickness values below 0.1 cm to 0
  mutate(Thickness_cm = case_when(
    Thickness_cm < 0.1 ~ 0,
    TRUE ~ Thickness_cm
  )) %>%
  mutate(Horizon = case_when(
    Horizon == "FF_thickness" ~ "Oie+a",
    Horizon == "Oie_thickness" ~ "Oie",
    Horizon == "Oa_thickness" ~ "Oa",
    Horizon == "Min_thickness" ~ "min",
  )) %>% 
  left_join(soil_mass_red, by = c("Year", "Plot", "Horizon")) %>% 
  #OM_TM: convert from kg/m2 in g/cm2
  mutate(BD_g_cm3 = (OM_TM * 0.1) / Thickness_cm) %>% 
  # Replace Inf values with NA
  mutate(BD_g_cm3 = na_if(BD_g_cm3, Inf))

summary(soil_BD$BD_g_cm3)

soil_BD %>% 
  ggplot(aes(x = BD_g_cm3, y = Thickness_cm)) +
  geom_point()

soil_BD %>% 
  drop_na(BD_g_cm3) %>% 
  ggplot(aes(x = Year, y = BD_g_cm3, color = Horizon)) +
  geom_point() +
  facet_wrap(~Elevation_grp) +
  theme_bw(base_size = 14) +
  geom_smooth(method = "lm")
ggsave(file = paste0("./Output/HBEF_BD_estimates_WS6_all_samples_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_BD %>% 
  group_by(Year, Horizon, Elevation_grp) %>% 
  reframe(mean_BD_g_cm3 = mean(BD_g_cm3, na.rm = TRUE),
          sd_BD_g_cm3 = sd(BD_g_cm3, na.rm = TRUE)) %>% 
  drop_na(mean_BD_g_cm3) %>% 
  ggplot(aes(x = Year, y = mean_BD_g_cm3, color = Horizon)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_BD_g_cm3 - sd_BD_g_cm3,
                    ymax = mean_BD_g_cm3 + sd_BD_g_cm3)) +
  facet_wrap(~Elevation_grp) +
  theme_bw(base_size = 14)
ggsave(file = paste0("./Output/HBEF_BD_estimates_WS6_averaged_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)

soil_BD %>% 
  ggplot(aes(x = Horizon, y = BD_g_cm3, color = Elevation_grp)) +
  geom_boxplot(notch = TRUE) +
  theme_bw(base_size = 14)

## Calculate averaged BD values for later use
# After visual inspection, calculate averaged values for each horizon
soil_BD_avg <- soil_BD %>% 
  group_by(Horizon) %>% 
  reframe(mean_BD_g_cm3 = mean(BD_g_cm3, na.rm = TRUE),
          sd_BD_g_cm3 = sd(BD_g_cm3, na.rm = TRUE),
          median_BD_g_cm3 = median(BD_g_cm3, na.rm = TRUE),
          mad_BD_g_cm3 = mad(BD_g_cm3, na.rm = TRUE)) %>% 
  #remove Oie+a for which we don't have samples
  filter(Horizon != "Oie+a")

soil_BD_avg

# Values seems reasonable, for European forests (0-10cm), median BD is 0.73 g/cm3
# https://www.sciencedirect.com/science/article/pii/S0167880924000252

# Same with organic layer data: Ol 0.03-0.12 g/cm3, Ofh: 0.07-0.27 g/cm3
# Findings based on a study from Poland: https://www.sciencedirect.com/science/article/pii/S0016706116304918

# Save file 
write_csv(x = soil_BD_avg, file = "./Data/soil_BD_avg_WS6.csv")
  
