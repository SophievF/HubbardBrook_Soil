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
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  geom_smooth(aes(color = Breakpoint), method = "lm") +
  scale_y_continuous(expression(paste("Soil respiration [g C m"^-2,"yr]")), 
                     expand = c(0,0)) +
  scale_x_continuous("Year", expand = c(0,0), breaks = seq(2005,2020,5))
ggsave(paste0("./Output/HBEF_FigureS6_", Sys.Date(), ".jpeg"),
       width = 12, height = 6)

lm_post_break <- lm(mean_gC_m2_yr ~ Year, 
   data = soil_resp_sum %>% 
     filter(Breakpoint == "post-breakpoint"))
summary(lm_post_break)

summary(soil_resp_sum$mean_gC_m2_yr)

