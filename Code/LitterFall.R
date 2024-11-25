## Hubbard Brook archived soil samples project ##
## Litter fall data ##
## Sophie von Fromm ##
## 2024-10-24 ##

library(tidyverse)

#litter fall: Oven dry weight of sample (g); litter trap: 0.097m2

litter_fall <- read_csv("./Data/HBEF_fine_litter_1992-2022.csv")

litter_fall_sum <- litter_fall %>% 
  filter(SITE == "BB" & TRTMT == "noCA" & ELEV == "Low") %>%
  dplyr::select(-contains("COUNT")) %>% 
  #remove "pooled" tags
  filter(TAG != "pooled") %>% 
  # replace missing values with long-term average
  mutate(DRY_MASS = replace(DRY_MASS, DRY_MASS < 0, NA)) %>% 
  mutate(DRY_MASS = ifelse(is.na(DRY_MASS), mean(DRY_MASS, na.rm = TRUE), DRY_MASS)) %>% 
  #calculate sum for each tag for each year and convert from g/0.097m2 to g/m2
  group_by(YEAR, TAG) %>% 
  summarise(sum_dry_mass = sum(DRY_MASS) * (1/0.097)) %>% 
  ungroup(TAG) %>% 
  summarise(mean_annual_litter_g = mean(sum_dry_mass),
            sd_annual_litter_g = sd(sum_dry_mass))
  

litter_fall_sum %>% 
  filter(YEAR < 2006) %>% 
  summary()

summary(litter_fall_sum$mean_annual_litter_g*0.49)

litter_fall_sum %>% 
  # multiply by 0.49 (~C content in foliage at low HBEF)
  ggplot(aes(x = YEAR, y = mean_annual_litter_g*0.49)) +
  geom_point() +
  geom_errorbar(aes(ymax = mean_annual_litter_g*0.49 + sd_annual_litter_g*0.49,
                    ymin = mean_annual_litter_g*0.49 - sd_annual_litter_g*0.49)) +
  geom_smooth() +
  theme_bw() +
  scale_y_continuous("Litter fall [gC/m2 yr]", expand = c(0,0)) +
  scale_x_continuous("Year", expand = c(0,0), breaks = seq(1992,2022,5)) +
  coord_cartesian(ylim = c(0,460))
ggsave("./Output/LitterFall_low_1992_2022.jpeg",
       width = 8, height = 6)
  
  