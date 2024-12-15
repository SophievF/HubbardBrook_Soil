## Hubbard Brook archived soil samples project ##
## Air temperature ##
## Sophie von Fromm ##
## 2024-11-15 ##

library(tidyverse)

temp_data <- read_csv("./Data/HBEF_air_temp_daily_1957-2024.csv")

temp <- temp_data %>% 
  mutate(Year = as.numeric(format(date, "%Y"))) %>% 
  filter(STA == "STA1") %>% 
  group_by(Year) %>% 
  filter(Year > 1955 & Year < 2024) %>% 
  mutate(AVE = mean(AVE)) %>% 
  mutate(breakpoint = case_when(
    Year <= 1980 ~ "pre",
    TRUE ~ "post"
  ))
  
temp %>%  
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
ggsave("./Output/HBEF_MAT_STA1.png", bg = "transparent",
       width = 6, height = 4.5)

lm_temp <- lm(AVE ~ Year, 
              data = temp %>% 
                filter(breakpoint == "post"))
summary(lm_temp)

temp %>% 
  filter(breakpoint == "post") %>% 
  count(Year)

lm_temp$coefficients[2]*43
