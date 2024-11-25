## Hubbard Brook archived soil samples project ##
## Selective metal analysis ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(ggpubr)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_selectiveAlFe_data_2024-09-20.csv") 

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
    geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
    scale_x_continuous("", limits = c(1998,2023))
}

#pyro Al
p1 <- fun_metal_time("al_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#pyro Fe
p2 <- fun_metal_time("fe_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Fe [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#oxal Al
p3 <- fun_metal_time("al_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#oxal Fe
p4 <- fun_metal_time("fe_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#dith Al
p5 <- fun_metal_time("al_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#dith Fe
p6 <- fun_metal_time("fe_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

ggarrange(p1, p3, p5, 
          p2, p4, p6, common.legend = TRUE)

ggsave(paste0("./Output/SelectiveMetall_all_", Sys.Date(),
              ".jpeg"), width = 18, height = 8)

##Check if quadratic function is significant

#pyro Al Oa/A
lm_quad_alpyr_oa <- lm(data = HBEF_data %>% 
                         drop_na(al_pyr_mg_g) %>% 
                         filter(Horizon == "Oa/A"), al_pyr_mg_g ~ Year + I(Year^2))
summary(lm_quad_alpyr_oa) #p > 0.05

lm_quad_alpyr_min <- lm(data = HBEF_data %>% 
                         drop_na(al_pyr_mg_g) %>% 
                         filter(Horizon == "Mineral_0_10"), 
                        al_pyr_mg_g ~ Year + I(Year^2))
summary(lm_quad_alpyr_min) #p > 0.05

#ox Al 
lm_quad_alox_oa <- lm(data = HBEF_data %>% 
                      drop_na(al_ox_mg_g) %>% 
                        filter(Horizon == "Oa/A"), al_ox_mg_g ~ Year + I(Year^2))
summary(lm_quad_alox_oa) #p > 0.05

lm_quad_alox_min <- lm(data = HBEF_data %>% 
                        drop_na(al_ox_mg_g) %>% 
                        filter(Horizon == "Mineral_0_10"), 
                       al_ox_mg_g ~ Year + I(Year^2))
summary(lm_quad_alox_min) #p > 0.05

#dith Al 
lm_quad_aldith_oa <- lm(data = HBEF_data %>% 
                     drop_na(al_dith_mg_g) %>% 
                       filter(Horizon == "Oa/A"), al_dith_mg_g ~ Year +I(Year^2))
summary(lm_quad_aldith_oa) #p > 0.05

lm_quad_aldith_min <- lm(data = HBEF_data %>% 
                          drop_na(al_dith_mg_g) %>% 
                          filter(Horizon == "Mineral_0_10"), 
                         al_dith_mg_g ~ Year +I(Year^2))
summary(lm_quad_aldith_min) #p > 0.05

#pyro Fe 
lm_quad_fepyr_oa <- lm(data = HBEF_data %>% 
                      drop_na(fe_pyr_mg_g) %>% 
                        filter(Horizon == "Oa/A"), fe_pyr_mg_g ~ Year + I(Year^2))
summary(lm_quad_fepyr_oa) #p > 0.05

lm_quad_fepyr_min <- lm(data = HBEF_data %>% 
                         drop_na(fe_pyr_mg_g) %>% 
                         filter(Horizon == "Mineral_0_10"), 
                        fe_pyr_mg_g ~ Year + I(Year^2))
summary(lm_quad_fepyr_min) #p > 0.05

#ox Fe 
lm_quad_feox_oa <- lm(data = HBEF_data %>% 
                     drop_na(fe_ox_mg_g) %>% 
                       filter(Horizon == "Oa/A"), fe_ox_mg_g ~ Year + I(Year^2))
summary(lm_quad_feox_oa) #p > 0.05

lm_quad_feox_min <- lm(data = HBEF_data %>% 
                        drop_na(fe_ox_mg_g) %>% 
                        filter(Horizon == "Mineral_0_10"), 
                       fe_ox_mg_g ~ Year + I(Year^2))
summary(lm_quad_feox_min) #p > 0.05

#dith Fe 
lm_quad_fedith_oa <- lm(data = HBEF_data %>% 
                       drop_na(fe_dith_mg_g) %>% 
                         filter(Horizon == "Oa/A"), fe_dith_mg_g ~ Year +I(Year^2))
summary(lm_quad_fedith_oa) #p > 0.05

lm_quad_fedith_min <- lm(data = HBEF_data %>% 
                          drop_na(fe_dith_mg_g) %>% 
                          filter(Horizon == "Mineral_0_10"), 
                         fe_dith_mg_g ~ Year +I(Year^2))
summary(lm_quad_fedith_min) #p > 0.05


#Function to plot metal data against C by horizon
fun_metal_oc <- function(metal){
  HBEF_data %>% 
    drop_na(.data[[metal]], `Measured_%_C`) %>% 
    ggplot(aes(y = `Measured_%_C`*10, x = .data[[metal]])) +
    geom_point(aes(color = Plot), size = 2) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~Horizon) +
    geom_smooth(method = "lm") +
    scale_y_continuous("Total C [mg/g]", transform = "log10")
}

#pyro Al
p7 <- fun_metal_oc("al_pyr_mg_g") +
  scale_x_continuous("Pyrophosphate extracted Al [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

#pyro Fe
p8 <-fun_metal_oc("fe_pyr_mg_g") +
  scale_x_continuous("Pyrophosphate extracted Fe [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

#oxal Al
p9 <-fun_metal_oc("al_ox_mg_g") +
  scale_x_continuous("Oxalate extracted Al [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

#oxal Fe
p10 <-fun_metal_oc("fe_ox_mg_g") +
  scale_x_continuous("Oxalate extracted Fe [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

#dith Al
p11 <-fun_metal_oc("al_dith_mg_g") +
  scale_x_continuous("Dithionite extracted Al [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

#dith Fe
p12 <-fun_metal_oc("fe_dith_mg_g") +
  scale_x_continuous("Dithionite extracted Fe [mg/g]", transform = "log10",
                     expand = c(0,0), limits = c(0.4,18))

ggarrange(p7, p9, p11, 
          p8, p10, p12, common.legend = TRUE)

ggsave(paste0("./Output/SelectiveMetall_SOC_all_log_", Sys.Date(),
              ".jpeg"), width = 18, height = 8)


#Function to plot summarized metal data over time by horizon
fun_metal_sum_time <- function(metal){
  HBEF_data %>% 
    drop_na(.data[[metal]]) %>% 
    group_by(Horizon, Year) %>% 
    summarise(mean = mean(.data[[metal]]),
              sd = sd(.data[[metal]])) %>% 
    ggplot(aes(x = Year, y = mean)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~Horizon) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
    scale_x_continuous("", limits = c(1998,2023))
}

#pyro Al
fun_metal_sum_time("al_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))

#pyro Fe
fun_metal_sum_time("fe_pyr_mg_g") +
  scale_y_continuous("Pyrophosphate extracted Fe [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))

#oxal Al
fun_metal_sum_time("al_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))

#oxal Fe
fun_metal_sum_time("fe_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))

#dith Al
fun_metal_sum_time("al_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))

#dith Fe
fun_metal_sum_time("fe_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,16)))


HBEF_data %>% 
  drop_na(al_ox_mg_g) %>% 
  ggplot(aes(x = Horizon, y = al_ox_mg_g)) +
  geom_boxplot(notch = TRUE)

HBEF_data %>% 
  drop_na(al_ox_mg_g) %>% 
  ggplot(aes(x = Horizon, y = fe_ox_mg_g)) +
  geom_boxplot(notch = TRUE)


