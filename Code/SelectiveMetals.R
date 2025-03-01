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
    group_by(Year, Horizon) %>% 
    summarise(mean_metal = mean(.data[[metal]]),
              sd_metal = sd(.data[[metal]])) %>% 
    ggplot(aes(x = Year, y = mean_metal, color = Horizon)) +
    geom_point(size = 2) +
    geom_path(linewidth = 1) +
    geom_errorbar(aes(ymax = mean_metal + sd_metal, ymin = mean_metal - sd_metal),
                  linewidth = 1) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          axis.title.x = element_blank()) +
    facet_wrap(~Horizon) +
    scale_x_continuous(limits = c(1998,2023)) +
    scale_color_manual(values = c("#b2df8a", "#a6cee3"))
}

#pyro Al
p1 <- fun_metal_time("al_pyr_mg_g") +
  scale_y_continuous("Pyroph. extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#pyro Fe
p2 <- fun_metal_time("fe_pyr_mg_g") +
  scale_y_continuous("Pyroph. extracted Fe [mg/g]",
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

ggarrange(p1, p3, p5, p2, p4, p6,
          labels = c("a)", "c)", "e)", "b)", "d)", "f)"))

ggsave(paste0("./Output/HBEF_FigureS7_", Sys.Date(),
              ".jpeg"), width = 18, height = 8)

#function to plot metal concentration by horizon
fun_metal <- function(metal){
  HBEF_data %>% 
    drop_na(.data[[metal]]) %>% 
    ggplot(aes(x = Horizon, y = .data[[metal]], fill = Horizon)) +
    geom_boxplot(notch = TRUE) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c("#b2df8a", "#a6cee3"))
}

#pyro Al
p7 <- fun_metal("al_pyr_mg_g") +
  scale_y_continuous("Pyroph. extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#pyro Fe
p8 <- fun_metal("fe_pyr_mg_g") +
  scale_y_continuous("Pyroph. extracted Fe [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#oxal Al
p9 <- fun_metal("al_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#oxal Fe
p10 <- fun_metal("fe_ox_mg_g") +
  scale_y_continuous("Oxalate extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#dith Al
p11 <- fun_metal("al_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Al [mg/g]",
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

#dith Fe
p12 <- fun_metal("fe_dith_mg_g") +
  scale_y_continuous("Dithionite extracted Fe [mg/g]", 
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(limits = c(0,18)))

ggarrange(p7, p9, p11, p8, p10, p12,
          labels = c("a)", "c)", "e)", "b)", "d)", "f)")) 

ggsave(paste0("./Output/HBEF_FigureS8_", Sys.Date(),
              ".jpeg"), width = 12, height = 8)





