## Hubbard Brook archived soil samples project ##
## SOC stocks for dataset 2 ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)

HBEF_data <- read_csv("./Data/HBEF_data_all_2024-11-07.csv") 

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

# Summarize and merge data by horizon; remove roots for now
HBEF_all <- HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  filter(Year != 2020) %>% 
  mutate(SOC_g_m2 = (`Measured_%_C` * mean_BD_g_cm3 * mean_thick_cm) * 100) %>% 
  dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, Elevation) %>% 
  mutate(Breakpoint = case_when(
    Year < 2015 ~ "pre-breakpoint",
    TRUE ~ "post-breakpoint"
  ))

#Groffman data only
HBEF_all %>% 
  filter(Horizon == "Oi/Oe") %>% 
  group_by(Year) %>%
  summarise(mean_SOC_g_m2 = mean(SOC_g_m2),
            sd_SOC_g_m2 = sd(SOC_g_m2)) %>% 
  ggplot(aes(x = Year, y = mean_SOC_g_m2)) +
  geom_point(shape = 21, size = 3, fill = "#F7746B") +
  geom_errorbar(aes(ymin =  mean_SOC_g_m2 - sd_SOC_g_m2,
                    ymax =  mean_SOC_g_m2 + sd_SOC_g_m2)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  # facet_wrap(~Horizon) +
  geom_smooth(method = "lm", fill = "#F7746B", color = "#F7746B") +
  scale_y_continuous("Mean SOC stock [gC/m2]", expand = c(0,0),
                     limits = c(0,2000)) +
  scale_x_continuous(expand = c(0,0), limits = c(1997,2024))

#Calculate decline in SOC stocks in Oie
lm_oie_C <- lm(mean_SOC_g_m2 ~ Year,
               data = HBEF_all %>%
                 filter(DataSource == "Groffman") %>%
                 filter(Horizon == "Oi/Oe") %>%
                 # filter(Year >= 2015) %>% 
                 group_by(Year) %>%
                 summarise(mean_SOC_g_m2 = mean(SOC_g_m2)))
summary(lm_oie_C)

lm_oa_C <- lm(mean_SOC_g_m2 ~ Year,
               data = HBEF_all %>%
                 filter(DataSource == "Groffman") %>%
                 filter(Horizon == "Oa/A") %>%
                 group_by(Year) %>%
                 summarise(mean_SOC_g_m2 = mean(SOC_g_m2)))
summary(lm_oa_C)

lm_min_C <- lm(mean_SOC_g_m2 ~ Year,
              data = HBEF_all %>%
                filter(DataSource == "Groffman") %>%
                filter(Horizon == "Mineral_0_10") %>%
                group_by(Year) %>%
                summarise(mean_SOC_g_m2 = mean(SOC_g_m2)))
summary(lm_min_C)
