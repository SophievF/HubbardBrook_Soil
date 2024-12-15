## Hubbard Brook archived soil samples project ##
## soil pH ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)

soil_data <- read_csv("./Data/HubbardBrook_microbialBiomass_1994-2023.csv")

pH_data <- soil_data %>% 
  filter(Treatment == "BearBrook") %>% 
  filter(El == "L") %>% 
  mutate(Year = as.numeric(format(Date, "%Y"))) %>% 
  mutate(pH = na_if(pH, -9999.99)) %>% 
  drop_na(pH) %>% 
  group_by(Year) %>% 
  mutate(pH_mean = mean(pH),
         pH_sd = sd(pH),
         breakpoint = case_when(
           Year < 2015 ~ "pre",
           TRUE ~ "post"
         )) %>% 
  ungroup()

pH_data %>%  
  ggplot(aes(x = Year, y = pH_mean, fill = breakpoint)) +
  geom_errorbar(aes(ymin = pH_mean - pH_sd, ymax = pH_mean + pH_sd)) +
  geom_point(shape = 21, size = 3,) +
  geom_smooth(aes(color = breakpoint), method = "lm", linewidth = 2, se = FALSE) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("pH", expand = c(0,0), limits = c(3.2,5.1))

ggsave("./Output/HBEF_pH.png", bg = "transparent",
       width = 6, height = 4.5)

lm_pH <- lm(pH ~ Year, 
            data = pH_data %>% 
              filter(breakpoint == "post"))
summary(lm_pH)

pH_data %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

lm_pH$coefficients[2]*8

pH_data_grp <- soil_data %>% 
  filter(Treatment == "BearBrook") %>% 
  filter(El == "L") %>% 
  mutate(Year = as.numeric(format(Date, "%Y"))) %>% 
  mutate(pH = na_if(pH, -9999.99)) %>% 
  drop_na(pH) %>% 
  group_by(Year, Hor) %>% 
  mutate(pH_mean = mean(pH),
         pH_sd = sd(pH),
         breakpoint = case_when(
           Year < 2015 ~ "pre",
           TRUE ~ "post"
         )) %>% 
  ungroup()

pH_data_grp$Hor <- factor(pH_data_grp$Hor, levels = c("Oi/Oe", "Oa/A", "Min"), 
                          ordered = TRUE)

pH_data_grp %>%  
  ggplot(aes(x = Year, y = pH_mean, fill = breakpoint)) +
  geom_errorbar(aes(ymin = pH_mean - pH_sd, ymax = pH_mean + pH_sd)) +
  geom_point(shape = 21, size = 3,) +
  facet_wrap(~Hor) +
  geom_smooth(aes(color = breakpoint), method = "lm", linewidth = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("pH", expand = c(0,0), limits = c(3.0,5.2))

ggsave("./Output/HBEF_pH_Horizon.png", bg = "transparent",
       width = 6, height = 4.5)

lm_pH_oie <- lm(pH ~ Year, 
            data = pH_data_grp %>% 
              filter(breakpoint == "post",
                     Hor == "Oi/Oe"))
summary(lm_pH_oie)

pH_data_grp %>% 
  filter(Hor == "Oi/Oe") %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

lm_pH_oie$coefficients[2]*8

lm_pH_oa <- lm(pH ~ Year, 
                data = pH_data_grp %>% 
                  filter(breakpoint == "post",
                         Hor == "Oa/A"))
summary(lm_pH_oa)

pH_data_grp %>% 
  filter(Hor == "Oi/Oe") %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

lm_pH_oa$coefficients[2]*8

lm_pH_oie$coefficients[2]*8

lm_pH_min <- lm(pH ~ Year, 
               data = pH_data_grp %>% 
                 filter(breakpoint == "post",
                        Hor == "Min"))
summary(lm_pH_min)

pH_data_grp %>% 
  filter(Hor == "Min") %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

lm_pH_min$coefficients[2]*8
