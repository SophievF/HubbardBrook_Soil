## Hubbard Brook archived soil samples project ##
## Radiocarbon analysis - constraint model with 14C and SOC ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(SoilR)
library(forecast)
library(FME)
library(ggpubr)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-11-07.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes)) 

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  group_by(Horizon) %>% 
  reframe(mean_14C = mean(Delta14C, na.rm = TRUE),
          mean_C = mean(`Measured_%_C`, na.rm = TRUE))

HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  group_by(Year, Horizon) %>% 
  summarise(mean_C = mean(`Measured_%_C`, na.rm = TRUE),
            sd_C = sd(`Measured_%_C`, na.rm = TRUE)) %>% 
  ggplot(aes(x = Year, y = mean_C)) +
  geom_point() +
  geom_errorbar(aes(ymax = mean_C + sd_C, ymin = mean_C - sd_C)) +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm")

HBEF_data %>%
  drop_na(Delta14C) %>%
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_point() +
  facet_wrap(~Horizon) +
  geom_smooth(method = "lm") +
  theme_bw()

# Load 14C data from Charley Driscoll
# LitterData <- read_csv("./Data/LitterData_Driscoll.csv")
# 
# #fm * exp(lambda * (-obs_date_y + 1950)) - 1) * 1000
lambda <- 0.0001209681
# 
# litter_all <- LitterData %>% 
#   filter(Watershed == 6) %>% 
#   mutate(Delta14C = (F14C * exp(lambda * (-Year + 1950)) -1) * 1000) %>% 
#   mutate(Elevation = case_when(
#     Plot < 70 ~ "High",
#     Plot > 152 ~ "Low",
#     TRUE ~ "Mid"
#   )) %>% 
#   #match Horizon names
#   mutate(Horizon = case_when(
#     Horizon == "Oie" ~ "Oi/Oe",
#     Horizon == "Oa" ~ "Oa/A"
#   )) %>% 
#   #Add NA's for SOC stocks so model will not be contrained by SOC stocsk 
#   mutate(SOC_g_m2 = NA)
# 
# # Summarize and merge data by horizon; remove roots for now
# HBEF_all <- HBEF_data %>% 
#   filter(Plot != "all fine roots") %>% 
#   filter(Year != 2020) %>% 
#   mutate(DataSource = "Groffman") %>% 
#   mutate(SOC_g_m2 = (`Measured_%_C` * mean_BD_g_cm3 * mean_thick_cm) * 100) %>% 
#   dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, DataSource, Elevation) %>% 
#   rbind(litter_all %>% 
#           dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, Elevation) %>% 
#           mutate(DataSource = "Driscoll")) %>% 
#   dplyr::select(Year:Delta14C, SOC_g_m2, DataSource, Elevation)

HBEF_all <- HBEF_data %>%
  filter(Plot != "all fine roots") %>%
  filter(Year != 2020) %>%
  mutate(DataSource = "Groffman") %>%
  mutate(SOC_g_m2 = (`Measured_%_C` * mean_BD_g_cm3 * mean_thick_cm) * 100) %>%
  dplyr::select(Year, Horizon, Delta14C, SOC_g_m2, DataSource, Elevation)

min_data_14C <- HBEF_all %>% 
  filter(Horizon == "Mineral_0_10") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(min_14C = mean(Delta14C, na.rm = TRUE),
            sd = sd(Delta14C, na.rm = TRUE)) %>%
  data.frame()

min_data_C <- HBEF_all %>% 
  filter(Horizon == "Mineral_0_10") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(min_C = mean(SOC_g_m2, na.rm = TRUE),
            sd = sd(SOC_g_m2, na.rm = TRUE)) %>%
  data.frame()

oa_data_14C <- HBEF_all %>% 
  filter(Horizon == "Oa/A") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(oa_14C = mean(Delta14C),
            sd = sd(Delta14C)) %>% 
  data.frame()

oa_data_C <- HBEF_all %>% 
  filter(Horizon == "Oa/A") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  drop_na(SOC_g_m2) %>% 
  summarise(oa_C = mean(SOC_g_m2),
            sd = sd(SOC_g_m2)) %>% 
  data.frame()

oie_data_14C <- HBEF_all %>% 
  filter(Horizon == "Oi/Oe") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(oie_14C = mean(Delta14C),
            sd = sd(Delta14C)) %>%
  # Replace NA for 1998 sd with ~26 which is the mean sd for all years in the Oie
  replace(is.na(.), 26.66504) %>% 
  data.frame()

oie_data_C <- HBEF_all %>% 
  filter(Horizon == "Oi/Oe") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  drop_na(SOC_g_m2) %>% 
  summarise(oie_C = mean(SOC_g_m2),
            sd = sd(SOC_g_m2)) %>%
  # Replace NA for 1998 sd with ~334 which is the mean sd for all years in the Oie
  mutate_at(vars(sd), ~replace_na(., 334)) %>%
  data.frame()

## Set-up model at steady state
## set-up model
# atmospheric 14C
NHZone2 <- bind.C14curves(prebomb = IntCal20, postbomb = Hua2021$NHZone2,
                          time.scale = "AD")

# atmospheric 14C forecast
y0 <- 1975 #initial year used for detecting trend in time series
q <- 4 #quarterly time-scale
y <- seq(y0, 2019.375, by = 1/q)

# Function to transform data series fo F' from F
FpfromF <- function(X){ # X must be a data.frame with first column time and second column Fraction Modern
  y = X[,1]
  fM = X[,2]
  Fp = fM*exp((1950-y)*lambda)
  return(data.frame(year = y, Fp = Fp))
}

FpNZ2 <- FpfromF(Hua2021$NHZone2[,c(1,4)]) #Absolute fraction modern F'
qNZ2 <- spline(FpNZ2, xout = y) #Spline interpolation of the NH_Zone 2 data set at a quarterly basis
NZ2 <- ts(qNZ2$y, start = y0, freq = q) #Transformation into a time-series object

mNZ2 <- ets(NZ2) #Fits an exponential smoothing state space model to the time series

foryrs <- 10 # Number of years to forecast
fNZ2 <- forecast(mNZ2, h = foryrs*q, level = c(69,90)) #Uses the fitted model to forecast 10 years into the future

NHZone2_2023 <- data.frame(Year = c(NHZone2[-dim(NHZone2)[1],1],
                                    seq(tsp(fNZ2$mean)[1], tsp(fNZ2$mean)[2], 
                                        by = 1/tsp(fNZ2$mean)[3])),
                           Delta14C = c(NHZone2[-dim(NHZone2)[1],2], 
                                        as.numeric((fNZ2$mean)-1)*1000))

NHZone2_2023 %>% 
  filter(Year >= 1949) %>% 
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() 

# time interval for model
# years <- seq(-53042, 2025, by = 0.5)
# years <- seq(-10000, 2025, by = 0.5)
years <- seq(1997, 2023, by = 0.5)

# initial C stocks in each pool (based 3-yr average)
# C0 <- c(mean(oie_data_C[1:3,2]),
#         mean(oa_data_C[1:3,2]), mean(min_data_C[1:3,2]))
#Average of first three years (including data from Driscoll)
# C0 <- c(1409, 1615, 2143)

## initial Delta14C in each pool (based on steady-state 3p model)
sens_14c <- read_csv("./Output/HBEF_3ps_steady_long_sens_14C_3_2024-12-14.csv") %>% 
  filter(Year == 1997)
init14C <- c(sens_14c[1,8], sens_14c[2,8], sens_14c[3,8])
init14C <- as.numeric(init14C)

sens_c <- read_csv("./Output/HBEF_3ps_steady_long_sens_C_3_2024-12-14.csv") %>% 
  filter(Year == 1997)
C0 <- c(sens_c[1,8], sens_c[2,8], sens_c[3,8])
C0 <- as.numeric(C0)
rm(sens_14c, sens_c)

# lag-time before C enters soils: based on communication with Josh
lag_time <- 3

# Number of model iterations
itr <- 15000

# C inputs
# In <- data.frame(year = years, Inputs = rep(210, length(years)))
In <- 210

## Set-up 14C pool (three pools in series)
# initial values are based on current C budget, assuming three pools (no roots)
ThreePSeriesModel_fun <- function(pars){
  mod = ThreepSeriesModel14(
    t = years,
    ks = pars[1:3],
    C0 = as.numeric(C0),
    F0_Delta14C = init14C,
    In = In,
    a21 = pars[4]*pars[1],
    a32 = pars[5]*pars[2],
    inputFc = NHZone2_2023,
    lag = lag_time
  )
  res_14C = getF14(mod)
  res_C = getC(mod)
  return(data.frame(time = years,
                    oie_14C = res_14C[,1],
                    oie_C = res_C[,1],
                    oa_14C = res_14C[,2],
                    oa_C = res_C[,2],
                    min_14C = res_14C[,3],
                    min_C = res_C[,3]))
}

tpsCost <- function(pars){
  funccall = ThreePSeriesModel_fun(pars)
  cost1 = modCost(model = funccall, obs = oie_data_14C, err = "sd")
  cost2 = modCost(model = funccall, obs = oa_data_14C, err = "sd", cost = cost1)
  cost3 = modCost(model = funccall, obs = min_data_14C, err = "sd", cost = cost2)
  cost4 = modCost(model = funccall, obs = oie_data_C, err = "sd", cost = cost3)
  cost5 = modCost(model = funccall, obs = oa_data_C, err = "sd", cost = cost4)
  cost6 = modCost(model = funccall, obs = min_data_C, err = "sd", cost = cost5)
  return(cost6)
}

#values based on stocks/fluxes; alternative use values from 3 pool model at steady state
init_pars <- c(k1 = 1/8, k2 = 1/19, k3 = 1/59, 
               alpha21 = 100/(100 + 110), alpha32 = 39/(39 + 61))

#double-check lower/upper again
tpsModelShortFit <- FME::modFit(f = tpsCost, p = init_pars, method = "Marq", 
                                upper = c(3, rep(1,4)), lower = rep(0,5)) 

#sum squared residuals
tpsModelShortFit$ssr

#mean squared residuals
tpsModelShortFit$ms

#AIC
(2*length(tpsModelShortFit$par))-(2*log(tpsModelShortFit$ms))

#Mean squared residuals per variable/horizon
sqrt(tpsModelShortFit$var_ms)

model_summary <- data.frame(ssr = tpsModelShortFit$ssr,
                            msr = tpsModelShortFit$ms,
                            aic = (2*length(tpsModelShortFit$par))-(2*log(tpsModelShortFit$ms)))

write.csv(model_summary, row.names = TRUE, quote = FALSE,
          file = paste0("./Output/HBEF_3ps_short_14C_C_summary_stats_", lag_time, "_",
                        Sys.Date(), ".csv"))

tpsVar <- tpsModelShortFit$var_ms_unweighted

#double-check lower/upper again
tpsShortMcmcFits <- FME::modMCMC(f = tpsCost, p = tpsModelShortFit$par, niter = itr, ntrydr = 5,
                                 updatecov = 50, var0 = tpsVar, upper = c(3, rep(1,4)),
                                 lower = rep(0,5)) #Create a new object to record fit stats

tpsModelShortOutput <- ThreePSeriesModel_fun(pars = as.numeric(summary(tpsShortMcmcFits)[1,1:5]))

# Save output
save(tpsModelShortFit, tpsShortMcmcFits, tpsModelShortOutput, 
     file = paste0("./Output/HBEF_3ps_short_14C_C_", lag_time, "_", Sys.Date(), ".Rdata"))
write_csv(summary(tpsShortMcmcFits), 
          file = paste0("./Output/HBEF_3ps_short_14C_C_summary_", lag_time, "_",
                        Sys.Date(), ".csv"))

# Create long dataframe
tpsModelOutput_df <- tpsModelShortOutput %>% 
  # filter(time > 1945) %>% 
  pivot_longer(!time,
               cols_vary = "slowest",
               names_to = c("Horizon", ".value"),
               names_pattern = "(.*)_(.*)") %>% 
  rename(Delta14C = `14C`,
         SOC_Stock = C,
         Year = time)

tpsModelOutput_df$Horizon <- factor(tpsModelOutput_df$Horizon,
                                    levels = c("oie", "oa", "min"),
                                    ordered = TRUE)

# Summarise HBEF data
HBEF_data_14C_C_sum <- HBEF_all %>% 
  filter(Year != 2020) %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C, na.rm = TRUE),
            Delta14C_sd = sd(Delta14C, na.rm = TRUE),
            C_mean = mean(SOC_g_m2, na.rm = TRUE),
            C_sd = sd(SOC_g_m2, na.rm = TRUE)) %>% 
  mutate(Horizon = case_when(
    Horizon == "Oi/Oe" ~ "oie",
    Horizon == "Oa/A" ~ "oa",
    Horizon == "Mineral_0_10" ~ "min"
  ))

HBEF_data_14C_C_sum$Horizon <- factor(HBEF_data_14C_C_sum$Horizon,
                                      levels = c("oie", "oa", "min"),
                                      ordered = TRUE)

tpsShortMcmcFits$bestpar

summary(tpsShortMcmcFits)

#Check for convergence: if model is converged, there should be no visible drift
jpeg(paste0("./Output/HBEF_3ps_short_14C_C_converg_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
plot(tpsShortMcmcFits)
dev.off()

jpeg(paste0("./Output/HBEF_3ps_short_14C_C_pairs_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
pairs(tpsShortMcmcFits)
dev.off()

#### Uncertainty analysis
pars <- tpsShortMcmcFits$pars

num <- 1000

sens_oie_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_14C"))) %>% 
  mutate(Horizon = "oie")

sens_oa_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                             sensvar = c("oa_14C"))) %>% 
  mutate(Horizon = "oa")

sens_min_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                              sensvar = c("min_14C"))) %>% 
  mutate(Horizon = "min")

sens_all_14C <- rbind(sens_oie_14C, sens_oa_14C, sens_min_14C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_14C, 
          file = paste0("./Output/HBEF_3ps_short_sens_14C_", lag_time, "_",
                        Sys.Date(), ".csv"))

sens_oie_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_C"))) %>% 
  mutate(Horizon = "oie")

sens_oa_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                 sensvar = c("oa_C"))) %>% 
  mutate(Horizon = "oa")

sens_min_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("min_C"))) %>% 
  mutate(Horizon = "min")

sens_all_C <- rbind(sens_oie_C, sens_oa_C, sens_min_C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_C, 
          file = paste0("./Output/HBEF_3ps_short_sens_C_", lag_time, "_",
                        Sys.Date(), ".csv"))

atm_mod <- data.frame(Year = NHZone2_2023$Year,
                      Mean = NA,
                      Sd = NA,
                      Min = NA,
                      Max = NA,
                      q05 = NA,
                      q25 = NA,
                      q50 = NHZone2_2023$Delta14C,
                      q75 = NA,
                      q95 = NA,
                      Horizon = "Atmosphere") %>% 
  dplyr::filter(Year > 1968)

sens_all_14C$Horizon <- factor(sens_all_14C$Horizon,
                               levels = c("oie", "oa", "min"), ordered = TRUE)

sens_all_C$Horizon <- factor(sens_all_C$Horizon,
                               levels = c("oie", "oa", "min"), ordered = TRUE)

sens_all_14C_p <- sens_all_14C %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = Horizon), alpha = 0.4) +
  geom_line(data = NHZone2_2023 %>% 
              filter(Year > 1997), aes(y = Delta14C)) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_C_sum,
                aes(y = Delta14C_mean, ymin = Delta14C_mean - Delta14C_sd,
                    ymax = Delta14C_mean + Delta14C_sd,
                    group = Horizon),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_C_sum, aes(y = Delta14C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1997,2024), expand = c(0,0),
                     breaks = seq(1969,2023,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), limits = c(-175,1000),
                     expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled", label = c("Oie", "Oa", "0-10 cm"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oie", "Oa/A", "0-10 cm"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) 

ggsave(file = paste0("./Output/HBEF_3ps_short_14C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6) 

sens_all_C_p <- sens_all_C %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = Horizon), alpha = 0.4) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_C_sum,
                aes(y = C_mean, ymin = C_mean - C_sd,
                    ymax = C_mean + C_sd,
                    group = Horizon),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_C_sum, aes(y = C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1997,2024), expand = c(0,0),
                     breaks = seq(1997,2023,10)) +
  scale_y_continuous("SOC stocks") +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled", label = c("Oie", "Oa", "0-10 cm"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oie", "Oa/A", "0-10 cm"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon, ncol = 1)

ggsave(file = paste0("./Output/HBEF_3ps_short_C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6)

ggarrange(sens_all_14C_p, sens_all_C_p, common.legend = TRUE)
ggsave(file = paste0("./Output/HBEF_3ps_short_14C_C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

sens_all_C %>% 
  filter(Horizon == "oie") %>% 
  ggplot(aes(x = Year)) +
  #add regression line based on measured values
  geom_smooth(data = HBEF_data_14C_C_sum %>% 
                filter(Horizon == "oie"), aes(y = C_mean),
              method = "lm", color = "black", linetype = "dashed") +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = Horizon), alpha = 0.4) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_C_sum %>% 
                  filter(Horizon == "oie"),
                aes(y = C_mean, ymin = C_mean - C_sd,
                    ymax = C_mean + C_sd),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_C_sum %>% 
               filter(Horizon == "oie"), aes(y = C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", expand = c(0,0), limits = c(1997,2024)) +
  scale_y_continuous("SOC stocks [g C/m2]", limits = c(750,1750), expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual("Modeled", label = c("Oie", "Oa", "0-10 cm"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oie", "Oa/A", "0-10 cm"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon) 
ggsave(file = paste0("./Output/HBEF_3ps_short_14C_C_Sensitivity_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 7, height = 5)

#Plot predicted vs observed (mean + SD) for each Horizon and compute linear regression
model_14C_pred_obs <- sens_all_14C %>%
  right_join(HBEF_data_14C_C_sum)

model_14C_pred_obs_res_oie <- model_14C_pred_obs %>% 
  filter(Horizon == "oie") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - Delta14C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_14C_pred_obs_res_oie %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_14C_pred_obs_res_oie),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_14C_pred_obs_res_oie$Delta14C_mean,
    x = model_14C_pred_obs_res_oie$Mean,
    method = "pearson")

model_14C_pred_obs_res_oa <- model_14C_pred_obs %>% 
  filter(Horizon == "oa") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - Delta14C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_14C_pred_obs_res_oa %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_14C_pred_obs_res_oa),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_14C_pred_obs_res_oa$Delta14C_mean,
    x = model_14C_pred_obs_res_oa$Mean,
    method = "pearson")

model_14C_pred_obs_res_min <- model_14C_pred_obs %>% 
  filter(Horizon == "min") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - Delta14C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_14C_pred_obs_res_min %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_14C_pred_obs_res_min),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_14C_pred_obs_res_min$Delta14C_mean,
    x = model_14C_pred_obs_res_min$Mean,
    method = "pearson")

fun_pred_obs_14C <- function(x){
  model_14C_pred_obs %>%
    filter(Horizon == x) %>% 
    ggplot(aes(x = Mean, y = Delta14C_mean, color = Horizon)) +
    geom_abline(intercept = 1, linetype = "dashed") +
    geom_point() +
    geom_errorbar(aes(ymin = Delta14C_mean - Delta14C_sd, 
                      ymax = Delta14C_mean + Delta14C_sd)) +
    geom_errorbar(aes(xmin = Mean - Sd, 
                      xmax = Mean + Sd)) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none") +
    facet_wrap(~Horizon) 
}

oie_14C <- fun_pred_obs_14C(x = "oie") +
  #Only plot the one that is significant
  geom_smooth(data = model_14C_pred_obs %>%
                filter(Horizon == "oie"),
              method = "lm") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(0,475), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(0,525), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))

oa_14C <- fun_pred_obs_14C(x = "oa") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-20,155), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(-20,230), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_14C <- fun_pred_obs_14C(x = "min") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-35,-10), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(-110,20), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_14C, oa_14C, min_14C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_3ps_short_14C_Obs_Pred_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

model_C_pred_obs <- sens_all_C %>%
  right_join(HBEF_data_14C_C_sum) %>% 
  drop_na(C_mean)

# manually calculate residuals and RMSE
model_C_pred_obs_res_oie <- model_C_pred_obs %>% 
  filter(Horizon == "oie") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_oie %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_oie),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_C_pred_obs_res_oie$C_mean,
    x = model_C_pred_obs_res_oie$Mean,
    method = "pearson")

model_C_pred_obs_res_oa <- model_C_pred_obs %>% 
  filter(Horizon == "oa") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_oa %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_oa),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_C_pred_obs_res_oa$C_mean,
    x = model_C_pred_obs_res_oa$Mean,
    method = "pearson")

model_C_pred_obs_res_min <- model_C_pred_obs %>% 
  filter(Horizon == "min") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_min %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_min),
            #take the square root
            rmse = sqrt(msr))

cor(y = model_C_pred_obs_res_min$C_mean,
    x = model_C_pred_obs_res_min$Mean,
    method = "pearson")

fun_pred_obs_C <- function(x){
  model_C_pred_obs %>%
    filter(Horizon == x) %>% 
    ggplot(aes(x = Mean, y = C_mean, color = Horizon)) +
    geom_abline(intercept = 1, linetype = "dashed") +
    geom_point() +
    geom_errorbar(aes(ymin = C_mean - C_sd, 
                      ymax = C_mean + C_sd)) +
    geom_errorbar(aes(xmin = Mean - Sd, 
                      xmax = Mean + Sd)) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none") +
    facet_wrap(~Horizon) 
}

oie_C <- fun_pred_obs_C(x = "oie") +
  scale_x_continuous("Predicted SOC stocks [g/m2]",
                     limits = c(1280,1580), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(750,1750), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))
ggsave(file = paste0("./Output/HBEF_3ps_short_C_Obs_Pred_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 5, height = 6)

oa_C <- fun_pred_obs_C(x = "oa") +
  scale_x_continuous("Predicted SOC stocks [g/m2]",
                     limits = c(1650,2010), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(750,3200), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_C <- fun_pred_obs_C(x = "min") +
  scale_x_continuous("Predicted SOC stocks [g/m2]",
                     limits = c(2175,2420), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(1000,4100), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_C, oa_C, min_C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_3ps_short_C_Obs_Pred_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

### Transit time and C age ###

## Calculate transit time and system age based on best parameter fits
#https://github.com/MPIBGC-TEE/Lanna/blob/v1.0.0/code/modelFitsLanna.R
# tau <- seq(0,500)
# # A1 <- meanModel@mat@map(1970)
# parsMCMC <- summary(tpsShortMcmcFits)
# 
# 
# A1min <- -1*diag(parsMCMC[3,1:3])
# A1min[2,1] <- abs(parsMCMC[3,4])*parsMCMC[3,1] 
# A1min[3,2] <- abs(parsMCMC[3,5])*parsMCMC[3,2]
# A1max <- -1*diag(parsMCMC[4,1:3])
# A1max[2,1] <- parsMCMC[4,4]*parsMCMC[4,1] 
# A1max[2,1] <- parsMCMC[4,5]*parsMCMC[4,2] 
# 
# 
# 
# SA1=systemAge(A=A1,u=c(1,0),a=tau)
# TT1=transitTime(A=A1,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
# TT1min<-transitTime(A=A1min,u=c(1,0),q=c(0.1,0.5,0.9))
# TT1max<-transitTime(A=A1max,u=c(1,0),q=c(0.1,0.5,0.9))
# 
# SA2=systemAge(A=A2,u=c(1,0),a=tau)
# TT2=transitTime(A=A2,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
# TT2min<-transitTime(A=A2min,u=c(1,0),q=c(0.1,0.5,0.9))
# TT2max<-transitTime(A=A2max,u=c(1,0),q=c(0.1,0.5,0.9))
# 
# 
# # see Stoner et al 2021 and Gonzalez-Sosa et al. 2024
# propagation_fun <- function(data, num_iter, input_vector){
#   
#   pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
#   
#   n_iter <- num_iter
#   
#   iter <- numeric(n_iter)
#   Iter_number <- numeric(n_iter)
#   System_age <- numeric(n_iter)
#   oie_age <- numeric(n_iter)
#   oa_age <- numeric(n_iter)
#   min_age <- numeric(n_iter)
#   Transit_time <- numeric(n_iter)
#   
#   
#   for (i in 1:n_iter){
#     subset_par <- data.frame(k1 = sample(data$k1, 1),
#                              k2 = sample(data$k2, 1),
#                              k3 = sample(data$k3, 1),
#                              alpha21 = sample(data$alpha21, 1),
#                              alpha31 = sample(data$alpha32, 1))
#     
#     #-------------- A matrix and inputs
#     ks <- subset_par[1:3]
#     A <- -1 * diag(ks) 
#     
#     ## Add transfer coefficients to A matrix:
#     alpha_2_1 <- subset_par[4]
#     alpha_3_2 <- subset_par[5]
#     
#     A[2,1] <- as.numeric(alpha_2_1*subset_par[1])
#     A[3,2] <- as.numeric(alpha_3_2*subset_par[2])
#     
#     u <- matrix((input_vector), ncol = 1)
#     
#     #---------------- Age and transit time ----------------------
#     ages <- c(0,500)
#     Sist_age <- systemAge(A = A, u = u, a = ages)
#     Trans_time <- transitTime(A = A, u = u, a = ages)
#     
#     Iter_number[i] <- i
#     System_age[i] <- as.numeric(Sist_age$meanSystemAge)
#     oie_age[i] <- as.numeric(Sist_age$meanPoolAge[1])
#     oa_age[i] <- as.numeric(Sist_age$meanPoolAge[2])
#     min_age[i] <- as.numeric(Sist_age$meanPoolAge[3])
#     Transit_time[i] <- as.numeric(Trans_time$meanTransitTime)
#     
#     
#     age_results <- as.data.frame(cbind(
#       Iter_number, System_age, oie_age, oa_age, min_age, Transit_time
#     ))
#     
#     setTxtProgressBar(pb, i)
#     
#   }
#   
#   return(age_results)
# }
# 
# # exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
# pars_df <- as.data.frame(tpsShortMcmcFits$pars[-(1:1000), ])
# 
# age_transit_dist <- propagation_fun(pars_df, 10000, C0)
# 
# write_csv(age_transit_dist, 
#           file = paste0("./Output/HBEF_3ps_short_Age_Transit_Distribution_", lag_time, "_",
#                         Sys.Date(), ".csv"))
# 
# age_transit_dist_sum <- age_transit_dist %>% 
#   pivot_longer(!Iter_number, values_to = "age_yr", names_to = "Pools") %>% 
#   group_by(Pools) %>% 
#   summarise(mean_age = mean(age_yr),
#             sd_age = sd(age_yr),
#             median_age = median(age_yr),
#             mad_age = mad(age_yr))
# 
# age_transit_dist %>% 
#   pivot_longer(!Iter_number, values_to = "age_yr", names_to = "Pools") %>% 
#   filter(Pools != "System_age", Pools != "Transit_time") %>%
#   mutate(Pools = factor(Pools, levels = c("oie_age", "oa_age", "min_age"),
#                         ordered = TRUE)) %>% 
#   ggplot(aes(x = age_yr, color = Pools)) +
#   geom_density(linewidth = 1) +
#   facet_wrap(~Pools, scales = "free") +
#   geom_vline(aes(xintercept = mean_age), 
#              data = age_transit_dist_sum %>% 
#                filter(Pools != "System_age",
#                       Pools != "Transit_time"),
#              linetype = "dashed") +
#   theme_bw() +
#   theme_classic(base_size = 16) +
#   theme(axis.text = element_text(color = "black"),
#         legend.position = "none") +
#   scale_color_manual(values = c("#33a02c", "#b2df8a", "#a6cee3")) +
#   scale_x_continuous("Age [yr]", expand = c(0,0))
# 
# ggsave(file = paste0("./Output/HBEF_3ps_short_Age_Distribution_", 
#                      lag_time, "_", Sys.Date(), ".jpeg"), width = 10, height = 6) 
# 
# 
# 
# 
