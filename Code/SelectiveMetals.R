## Hubbard Brook archived soil samples project ##
## Selective metal analysis ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-08-06.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes))

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)


HBEF_data %>% 
  drop_na(al_pyr_mg_g) %>% 
  ggplot(aes(x = Year, y = al_pyr_mg_g)) +
  geom_point(aes(color = Plot), size = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm") +
  scale_x_continuous("", limits = c(1998,2023)) +
  scale_y_continuous("Pyrophosphate extracted Al [mg/g]")

HBEF_data %>% 
  drop_na(al_pyr_mg_g, `Measured_%_C`) %>% 
  ggplot(aes(y = `Measured_%_C`*10, x = al_pyr_mg_g)) +
  geom_point(aes(color = Plot), size = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm") +
  scale_y_continuous("Total C [mg/g]") +
  scale_x_continuous("Pyrophosphate extracted Al [mg/g]")

HBEF_data %>% 
  drop_na(fe_pyr_mg_g) %>% 
  ggplot(aes(x = Year, y = fe_pyr_mg_g)) +
  geom_point(aes(color = Plot), size = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm") +
  scale_x_continuous("", limits = c(1998,2023)) +
  scale_y_continuous("Pyrophosphate extracted Fe [mg/g]")

HBEF_data %>% 
  drop_na(fe_pyr_mg_g, `Measured_%_C`) %>% 
  ggplot(aes(y = `Measured_%_C`*10, x = fe_pyr_mg_g)) +
  geom_point(aes(color = Plot), size = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm") +
  scale_y_continuous("Total C [mg/g]") +
  scale_x_continuous("Pyrophosphate extracted Fe [mg/g]")
