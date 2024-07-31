## Hubbard Brook archived soil samples project ##
## Grid and elevation map WS6 ##
## Sophie von Fromm ##
## 2024-07-31 ##

# https://hubbardbrook.org/experimental-watersheds-research-sites/watershed-6/

library(tidyverse)

# Load grid and elevation file
grid_elv <- read.table("./Data/Watershed6_grid_elevations.txt", header = TRUE,
                       sep = ",", dec = ".")

grid_elv_df <- grid_elv %>% 
  drop_na(Plot) %>% 
  tibble() %>% 
  mutate(Elevation_grp = case_when(
    Plot <= 86 ~ "Upper",
    Plot > 86 & Plot <= 158 ~ "Middle",
    Plot > 158 ~ "Lower"
  ))

# Check if elevation assignment makes sense
grid_elv_df %>% 
  ggplot(aes(x = Plot, y = Elevation_m, color = Elevation_grp)) +
  geom_point() +
  theme_bw()

grid_elv_df %>% 
  group_by(Elevation_grp) %>% 
  reframe(mean_elv = mean(Elevation_m))

# Save file 
write_csv(x = grid_elv_df, file = "./Data/Grid_Elevation_WS6.csv")
