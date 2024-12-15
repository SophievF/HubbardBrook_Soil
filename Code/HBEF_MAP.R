## Hubbard Brook archived soil samples project ##
## Annual precipitation ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)

map_data <- read_csv("./Data/dailyWatershedPrecip1956-2024.csv")

map <- map_data %>% 
  filter(watershed == "W6") %>% 
  mutate(Year = as.numeric(format(DATE, "%Y"))) %>% 
  filter(Year < 2024) %>% 
  group_by(Year) %>% 
  mutate(AP = sum(Precip)) %>% 
  ungroup() %>% 
  distinct(Year, .keep_all = TRUE)
  
map %>%  
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
ggsave("./Output/HBEF_AnnualPrecip.png", bg = "transparent",
       width = 6, height = 4.5)

lm_map <- lm(AP ~ Year, data = map)
summary(lm_map)

map %>% 
  count(Year)

lm_map$coefficients[2]*60
219