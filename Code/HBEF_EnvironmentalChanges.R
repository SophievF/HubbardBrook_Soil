## Hubbard Brook archived soil samples project ##
## Environmental data ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)
library(ggpubr)

### Precipitation ####
map_data <- read_csv("./Data/dailyWatershedPrecip1956-2024.csv")

map <- map_data %>% 
  filter(watershed == "W6") %>% 
  mutate(Year = as.numeric(format(DATE, "%Y"))) %>% 
  filter(Year < 2024) %>% 
  group_by(Year) %>% 
  mutate(AP = sum(Precip)) %>% 
  ungroup() %>% 
  distinct(Year, .keep_all = TRUE)

summary(map)
  
map_p <- map %>%  
  ggplot(aes(x = Year, y = AP)) +
  geom_line(linewidth = 0.5) +
  geom_smooth(method = "lm", linewidth = 2, fill = "#F7746B", color = "#F7746B", 
              se = FALSE) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Annual precipitation [mm]", limits = c(1000,2000), 
                     expand = c(0,0))

lm_map <- lm(AP ~ Year, data = map)
#Intercept is 3.6 (+3.6 mm per year)
summary(lm_map)

map %>% 
  count(Year)

# Increase since the mid-1990s (last 28 years)
lm_map$coefficients[2]*28

### MAT ####
temp_data <- read_csv("./Data/HBEF_air_temp_daily_1957-2024.csv")

temp <- temp_data %>% 
  mutate(Year = as.numeric(format(date, "%Y"))) %>% 
  #only select station 1 (most representative for sampling locations)
  filter(STA == "STA1") %>% 
  group_by(Year) %>% 
  filter(Year >= 1964 & Year < 2024) %>% 
  mutate(AVE = mean(AVE)) %>% 
  mutate(breakpoint = case_when(
    Year <= 1980 ~ "pre",
    TRUE ~ "post"
  ))

temp_data %>% 
  filter(STA == "STA1") %>% 
  mutate(Month = as.numeric(format(date, "%m"))) %>% 
  group_by(Month) %>% 
  summarise(MAT = mean(AVE))

temp_p <- temp %>%  
  ggplot(aes(x = Year, y = AVE)) +
  geom_line(linewidth = 0.5) +
  geom_smooth(aes(color = breakpoint),method = "lm", linewidth = 2, se = FALSE) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Mean annual temperature [°C]")

lm_temp <- lm(AVE ~ Year, 
              data = temp %>% 
                filter(breakpoint == "post"))

# Intercept is 0.05 (+0.05C per year)
summary(lm_temp)

temp %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

# Increase since the mid-1990s (last 28 years)
lm_temp$coefficients[2]*28

### soil pH ####
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

pH_p <- pH_data %>%  
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

lm_pH <- lm(pH ~ Year, 
            data = pH_data %>% 
              filter(breakpoint == "post"))

# Intercept is 0.07 (+0.07 per year)
summary(lm_pH)

pH_data %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

# Increase since break-point (~2015)
lm_pH$coefficients[2]*8

### Plot all data together ####
ggarrange(map_p, temp_p, pH_p, nrow = 1,
          labels = c("a)", "b)", "c)"))
ggsave(paste0("./Output/HBEF_Figure1_", Sys.Date(),
              ".jpeg"), width = 11, height = 5)

## pH by Horizon shows similar trends
# pH_data_grp <- soil_data %>% 
#   filter(Treatment == "BearBrook") %>% 
#   filter(El == "L") %>% 
#   mutate(Year = as.numeric(format(Date, "%Y"))) %>% 
#   mutate(pH = na_if(pH, -9999.99)) %>% 
#   drop_na(pH) %>% 
#   group_by(Year, Hor) %>% 
#   mutate(pH_mean = mean(pH),
#          pH_sd = sd(pH),
#          breakpoint = case_when(
#            Year < 2015 ~ "pre",
#            TRUE ~ "post"
#          )) %>% 
#   ungroup()
# 
# pH_data_grp$Hor <- factor(pH_data_grp$Hor, levels = c("Oi/Oe", "Oa/A", "Min"), 
#                           ordered = TRUE)
# 
# pH_data_grp %>%  
#   ggplot(aes(x = Year, y = pH_mean, fill = breakpoint)) +
#   geom_errorbar(aes(ymin = pH_mean - pH_sd, ymax = pH_mean + pH_sd)) +
#   geom_point(shape = 21, size = 3,) +
#   facet_wrap(~Hor) +
#   geom_smooth(aes(color = breakpoint), method = "lm", linewidth = 2) +
#   theme_classic(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA)) +
#   scale_x_continuous("Year", expand = c(0,0)) +
#   scale_y_continuous("pH", expand = c(0,0), limits = c(3.0,5.2))
# 
# lm_pH_oie <- lm(pH ~ Year, 
#                 data = pH_data_grp %>% 
#                   filter(breakpoint == "post",
#                          Hor == "Oi/Oe"))
# summary(lm_pH_oie)
# 
# pH_data_grp %>% 
#   filter(Hor == "Oi/Oe") %>% 
#   filter(breakpoint == "post") %>% 
#   count(Year)
# 
# lm_pH_oie$coefficients[2]*8
# 
# lm_pH_oa <- lm(pH ~ Year, 
#                data = pH_data_grp %>% 
#                  filter(breakpoint == "post",
#                         Hor == "Oa/A"))
# summary(lm_pH_oa)
# 
# pH_data_grp %>% 
#   filter(Hor == "Oi/Oe") %>% 
#   filter(breakpoint == "post") %>% 
#   count(Year)
# 
# lm_pH_oa$coefficients[2]*8
# 
# lm_pH_oie$coefficients[2]*8
# 
# lm_pH_min <- lm(pH ~ Year, 
#                 data = pH_data_grp %>% 
#                   filter(breakpoint == "post",
#                          Hor == "Min"))
# summary(lm_pH_min)
# 
# pH_data_grp %>% 
#   filter(Hor == "Min") %>% 
#   filter(breakpoint == "post") %>% 
#   count(Year)
# 
# lm_pH_min$coefficients[2]*8

