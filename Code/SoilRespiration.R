## Hubbard Brook archived soil samples project ##
## Soil respiration data ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)

soil_resp <- read_csv("./Data/HubbardBrook_TraceGas_Watershed1_BearBrook_2002-2022.csv")

head(soil_resp)

soil_resp_sum <- soil_resp %>% 
  filter(Site == "Bear Brook",
         Elevation == "low") %>% 
  mutate(Year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(Year) %>% 
  summarise(mean_CO2 = mean(CO2_rate),
            sd_CO2 = sd(CO2_rate)) %>% 
  mutate(mean_gC_m2_yr = (mean_CO2*0.2729)*8760,
         sd_gC_m2_yr = (sd_CO2*0.2729)*8760) %>% 
  mutate(Breakpoint = case_when(
    Year > 2014 ~ "post-breakpoint",
    TRUE ~ "pre-breakpoint"
  )) 

soil_resp_sum %>% 
  ggplot(aes(x = Year, y = mean_gC_m2_yr, fill = Breakpoint)) +
  geom_point(shape = 21, size = 3) +
  geom_errorbar(aes(ymin =  mean_gC_m2_yr - sd_gC_m2_yr,
                    ymax =  mean_gC_m2_yr + sd_gC_m2_yr)) +
  theme_bw() +
  geom_smooth(aes(color = Breakpoint), method = "lm")

lm_post_break <- lm(mean_gC_m2_yr ~ Year, 
   data = soil_resp_sum %>% 
     filter(Breakpoint == "post-breakpoint"))
summary(lm_post_break)

summary(soil_resp_sum$mean_gC_m2_yr)

soil_resp_sum %>% 
  ggplot(aes(x = Year, y = mean_gC_m2_yr)) +
  geom_point(shape = 21, size = 3) +
  geom_errorbar(aes(ymin =  mean_gC_m2_yr - sd_gC_m2_yr,
                    ymax =  mean_gC_m2_yr + sd_gC_m2_yr)) +
  theme_bw() +
  geom_smooth() +
  scale_y_continuous("Soil respiration [gC/m2 yr]", expand = c(0,0)) +
  scale_x_continuous("Year", expand = c(0,0), breaks = seq(2002,2022,5)) +
  coord_cartesian(ylim = c(0,460))

ggsave("./Output/SoilRespiration_low_2002_2022.jpeg",
       width = 8, height = 6)         
