## Hubbard Brook archived soil samples project ##
## Selective metal analysis ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-09-.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes))

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

#Function to plot metal data over time by horizon
fun_metal_time <- function(metal){
  HBEF_data %>% 
    drop_na(.data[[metal]]) %>% 
    ggplot(aes(x = Year, y = .data[[metal]])) +
    geom_point(aes(color = Plot), size = 2) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~Horizon) +
    geom_smooth(method = "lm") +
    scale_x_continuous("", limits = c(1998,2023))
}

#pyro Al
fun_metal_time("al_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Al [mg/g]", limits = c(0,18),
                     expand = c(0,0))

#pyro Fe
fun_metal_time("fe_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Fe [mg/g]", limits = c(0,16),
                     expand = c(0,0))

#oxal Al
fun_metal_time("al_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Al [mg/g]", limits = c(0,18),
                     expand = c(0,0))

#oxal Fe
fun_metal_time("fe_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Fe [mg/g]", limits = c(0,16),
                     expand = c(0,0))

#Function to plot metal data against by horizon
fun_metal_oc <- function(metal){
  HBEF_data %>% 
    drop_na(.data[[metal]], `Measured_%_C`) %>% 
    ggplot(aes(y = `Measured_%_C`*10, x = .data[[metal]])) +
    geom_point(aes(color = Plot), size = 2) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~Horizon) +
    geom_smooth(method = "lm") +
    scale_y_continuous("Total C [mg/g]")
}

#pyro Al
fun_metal_oc("al_pyr_mg_g") +
  scale_x_continuous("Pyrophosphate extracted Al [mg/g]", limits = c(0,18),
                     expand = c(0,0))

#pyro Fe
fun_metal_oc("fe_pyr_mg_g") +
  scale_x_continuous("Pyrophosphate extracted Fe [mg/g]", limits = c(0,16),
                     expand = c(0,0))

#oxal Al
fun_metal_oc("al_ox_mg_g") +
  scale_x_continuous("Oxalate extracted Al [mg/g]", limits = c(0,18),
                     expand = c(0,0))

#oxal Fe
fun_metal_oc("fe_ox_mg_g") +
  scale_x_continuous("Oxalate extracted Fe [mg/g]", limits = c(0,16),
                     expand = c(0,0))




HBEF_data %>% 
  drop_na(fe_pyr_mg_g) %>% 
  group_by(Horizon, Year) %>% 
  summarise(mean = mean(fe_pyr_mg_g),
            sd = sd(fe_pyr_mg_g)) %>% 
  ggplot(aes(x = Year)) +
  geom_point(aes(y = mean), size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  facet_wrap(~Horizon) +
  scale_x_continuous("", limits = c(1998,2023)) +
  scale_y_continuous("Pyrophosphate extracted Fe [mg/g]")


