## Hubbard Brook archived soil samples project ##
## Radiocarbon analysis - constraint model with 14C and SOC ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(SoilR)
library(forecast)
library(FME)
library(ggpubr)
library(bayestestR)

# Load data
HBEF_all <- read_csv("./Data/HBEF_ModelingDatasets_all_2025-01-02.csv")

HBEF_all$Horizon <- factor(HBEF_all$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral"),
                            ordered = TRUE)

## Prepare each horizon for model
min_data_14C <- HBEF_all %>% 
  filter(Horizon == "Mineral") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(min_14C = mean(Delta14C, na.rm = TRUE),
            sd = sd(Delta14C, na.rm = TRUE)) %>%
  data.frame()

min_data_C <- HBEF_all %>% 
  filter(Horizon == "Mineral") %>% 
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
  data.frame()

# Replace NA for 1998 sd with the mean sd for all years in the Oie from dataset 2
oie_data_14C$sd <- replace(oie_data_14C$sd, is.na(oie_data_14C$sd),
                           mean(oie_data_14C$sd[c(6:14)]))

oie_data_C <- HBEF_all %>% 
  filter(Horizon == "Oi/Oe") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  drop_na(SOC_g_m2) %>% 
  summarise(oie_C = mean(SOC_g_m2),
            sd = sd(SOC_g_m2)) %>%
  data.frame()

# Replace NA for 1998 sd with the mean sd for all years in the Oie from dataset 2
oie_data_C$sd <- replace(oie_data_C$sd, is.na(oie_data_C$sd),
                         mean(oie_data_C$sd[c(6:14)]))

#### Model set-up ####

## Define lambda 
lambda <- 0.0001209681

## Create atmospheric 14C curve
NHZone2 <- bind.C14curves(prebomb = IntCal20, postbomb = Hua2021$NHZone2,
                          time.scale = "AD")

# Make atmospheric 14C forecast (Hua 2021 only goes until 2019)
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
# Check atmospheric curve
NHZone2_2023 %>% 
  filter(Year >= 1949) %>% 
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() 

## Define initial values for model 
# Define time interval for model
years <- seq(-53042, 2023, by = 0.5)

# initial C stocks based on average values
C0 <- c(mean(oie_data_C[,2]),
        mean(oa_data_C[,2]), mean(min_data_C[,2]))

## initial Delta14C in each pool 
init14C <- c(0, 0, -23)

# lag-time before C enters soils: based on communication with Josh
lag_time <- 3

# Number of model iterations
itr <- 15000

# C inputs (based on values from Fahey et al. 2005)
In <- 210

# Initial k and a values (k = stocks/flux; modified from Fahey et al 2005)
init_pars <- c(k1 = 1/7, k2 = 1/19, k3 = 1/59,
               alpha21 = 100/(100 + 110), alpha32 = 39/(39 + 61))

## Set-up model (three pools in series)
ThreePSeriesModel_fun <- function(pars){
  mod = SoilR::ThreepSeriesModel14(
    # time interval for model to find solutions
    t = years,
    # decomposition rates (k)
    ks = pars[1:3],
    # initial SOC stocks
    C0 = as.numeric(C0),
    # Initial 14C values
    F0_Delta14C = init14C,
    # C inputs (constant over time based on litter fall data)
    In = In,
    # transfer rate from pool 1 to 2
    a21 = pars[4]*pars[1],
    # transfer rate from pool 2 to 3
    a32 = pars[5]*pars[2],
    # atmospheric 14C for each time point
    inputFc = NHZone2_2023,
    # lag time for radiocarbon when it enters the system
    lag = lag_time
  )
  # get 14C for each pool and each time point
  res_14C = SoilR::getF14(mod)
  # get C for each pool and each time point
  res_C = SoilR::getC(mod)
  # return dataframe with years and 14C and C values for each pool
  return(data.frame(time = years,
                    oie_14C = res_14C[,1],
                    oie_C = res_C[,1],
                    oa_14C = res_14C[,2],
                    oa_C = res_C[,2],
                    min_14C = res_14C[,3],
                    min_C = res_C[,3]))
}

# Cost function to constrain model with observational data
tpsCost <- function(pars){
  funccall = ThreePSeriesModel_fun(pars)
  cost1 = FME::modCost(model = funccall, obs = oie_data_14C, err = "sd")
  cost2 = FME::modCost(model = funccall, obs = oa_data_14C, err = "sd", cost = cost1)
  cost3 = FME::modCost(model = funccall, obs = min_data_14C, err = "sd", cost = cost2)
  cost4 = FME::modCost(model = funccall, obs = oie_data_C, err = "sd", cost = cost3)
  cost5 = FME::modCost(model = funccall, obs = oa_data_C, err = "sd", cost = cost4)
  cost6 = FME::modCost(model = funccall, obs = min_data_C, err = "sd", cost = cost5)
  return(cost6)
}

#### NO NEED TO RUN THE FOLLOWING CODE (OBJECTS HAVE BEEN SAVED) #####
## Un-comment the code if you wish to run it ##

## Fitting the model
# tpsModelFit <- FME::modFit(f = tpsCost, p = init_pars, method = "Marq", 
#                            upper = c(3, rep(1,4)), lower = rep(0,5))

## Model evaluation
# sum squared residuals
# tpsModelFit$ssr

# mean squared residuals
# tpsModelFit$ms

# AIC
# (2*length(tpsModelFit$par))-(2*log(tpsModelFit$ms))

# Mean squared residuals per variable/horizon
# sqrt(tpsModelFit$var_ms)

# model_summary <- data.frame(ssr = tpsModelFit$ssr,
#                             msr = tpsModelFit$ms,
#                             aic = (2*length(tpsModelFit$par))-(2*log(tpsModelFit$ms)))

# write.csv(model_summary, row.names = TRUE, quote = FALSE,
#           file = paste0("./Output/HBEF_3ps_steady_long_14C_C_summary_stats_", lag_time, "_",
#                         Sys.Date(), ".csv"))

## Perform Markov Chain Monte Carlo simulation
# tpsVar <- tpsModelFit$var_ms_unweighted
# tpsMcmcFits <- FME::modMCMC(f = tpsCost, p = tpsModelFit$par, niter = itr, ntrydr = 5,
#                             updatecov = 50, var0 = tpsVar, upper = c(3, rep(1,4)),
#                             lower = rep(0,5)) #Create a new object to record fit stats

# Fit model with mean parameter values
# tpsModelOutput <- ThreePSeriesModel_fun(pars = as.numeric(summary(tpsMcmcFits)[1,1:6]))

# write.csv(summary(tpsMcmcFits), 
#           file = paste0("./Output/HBEF_3ps_steady_long_14C_C_summary_", lag_time, "_",
#                         Sys.Date(), ".csv"))

# Check for convergence: if model is converged, there should be no visible drift
# plot(tpsMcmcFits)
# 
# pairs(tpsMcmcFits)

#### Uncertainty analysis
# exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
# pars <- tpsMcmcFits$pars[-(1:1000), ]

# Number of times the model will be run
# num <- 1000

# Estimate uncertainties for each horizon
# sens_oie_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                                   sensvar = c("oie_14C"))) %>% 
#   mutate(Horizon = "oie")
 
# sens_oa_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                              sensvar = c("oa_14C"))) %>% 
#   mutate(Horizon = "oa")
 
# sens_min_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                               sensvar = c("min_14C"))) %>% 
#   mutate(Horizon = "min")

# Combine results for each horizon into one data frame and save as csv file
# sens_all_14C <- rbind(sens_oie_14C, sens_oa_14C, sens_min_14C) %>% 
#   tibble() %>% 
#   dplyr::rename(Year = x) %>% 
#   filter(Year > 1968)

# write_csv(sens_all_14C, 
#           file = paste0("./Output/HBEF_3ps_steady_long_sens_14C_", lag_time, "_",
#                         Sys.Date(), ".csv"))

# Estimate uncertainties for each horizon
# sens_oie_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                                   sensvar = c("oie_C"))) %>% 
#   mutate(Horizon = "oie")

# sens_oa_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                                  sensvar = c("oa_C"))) %>% 
#   mutate(Horizon = "oa")

# sens_min_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
#                                   sensvar = c("min_C"))) %>% 
#   mutate(Horizon = "min")

# Combine results for each horizon into one data frame and save as csv file
# sens_all_C <- rbind(sens_oie_C, sens_oa_C, sens_min_C) %>% 
#   tibble() %>% 
#   dplyr::rename(Year = x) %>% 
#   filter(Year > 1968)

# write_csv(sens_all_C, 
#           file = paste0("./Output/HBEF_3ps_steady_long_sens_C_", lag_time, "_",
#                         Sys.Date(), ".csv"))

# Save all model outputs
# save(tpsModelFit, tpsMcmcFits, tpsModelOutput, sens_all_C, sens_all_14C,
#      file = paste0("./Output/HBEF_3ps_steady_long_14C_C_", lag_time, "_",
#                    Sys.Date(), ".Rdata"))

##### Plot and evaluate model results ######
# Load model results (from code above)
load("./Output/HBEF_3ps_steady_long_14C_C_3_2025-01-02.RData")

# Create long dataframe from model output
tpsModelOutput_df <- tpsModelOutput %>% 
  pivot_longer(!time,
               cols_vary = "slowest",
               names_to = c("Horizon", ".value"),
               names_pattern = "(.*)_(.*)") %>% 
  rename(Delta14C = `14C`,
         SOC_Stock = C,
         Year = time) %>% 
  mutate(Horizon = case_when(
    Horizon == "oie" ~ "Oi/Oe",
    Horizon == "oa" ~ "Oa/A",
    Horizon == "min" ~ "Mineral"))

tpsModelOutput_df$Horizon <- factor(tpsModelOutput_df$Horizon,
                                    levels = c("Oi/Oe", "Oa/A", "Mineral"),
                                    ordered = TRUE)

# Check best parameter values
tpsMcmcFits$bestpar

# Check summary statistics for estimated parameter values
summary(tpsMcmcFits, quantile.type = c(0.005,0.5,0.95))[,1:5]

tpsMcmcFits$pars %>% 
  as.data.frame() %>%
  summarise(median_k1 = median(k1),
            q05_k1 = quantile(k1, 0.05),
            q95_k1 = quantile(k1, 0.95),
            median_a21 = median(alpha21),
            q05_a21 = quantile(alpha21, 0.05),
            q95_a21 = quantile(alpha21, 0.95),
            median_k2 = median(k2),
            q05_k2 = quantile(k2, 0.05),
            q95_k2 = quantile(k2, 0.95),
            median_k3 = median(k3),
            q05_k3 = quantile(k3, 0.05),
            q95_k3 = quantile(k3, 0.95),
            median_a32 = median(alpha32),
            q05_a32 = quantile(alpha32, 0.05),
            q95_a32 = quantile(alpha32, 0.95))

## Prepare data for plotting
# Summarise HBEF data for plotting and merging with model results
HBEF_data_14C_C_sum <- HBEF_all %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C, na.rm = TRUE),
            Delta14C_sd = sd(Delta14C, na.rm = TRUE),
            C_mean = mean(SOC_g_m2, na.rm = TRUE),
            C_sd = sd(SOC_g_m2, na.rm = TRUE)) 

# Rename horizon names from sensitivity analysis
sens_all_14C <- sens_all_14C %>% 
  mutate(Horizon = case_when(
    Horizon == "oie" ~ "Oi/Oe",
    Horizon == "oa" ~ "Oa/A",
    Horizon == "min" ~ "Mineral"
  ))

sens_all_14C$Horizon <- factor(sens_all_14C$Horizon,
                               levels = c("Oi/Oe", "Oa/A", "Mineral"),
                               ordered = TRUE)

sens_all_C <- sens_all_C %>% 
  mutate(Horizon = case_when(
    Horizon == "oie" ~ "Oi/Oe",
    Horizon == "oa" ~ "Oa/A",
    Horizon == "min" ~ "Mineral"
  ))

sens_all_C$Horizon <- factor(sens_all_C$Horizon,
                             levels = c("Oi/Oe", "Oa/A", "Mineral"),
                             ordered = TRUE)

## Plot 14C results
sens_all_14C_p <- sens_all_14C %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = Horizon), alpha = 0.4) +
  geom_line(data = NHZone2_2023 %>% 
              filter(Year > 1968), aes(y = Delta14C)) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_C_sum,
                aes(y = Delta14C_mean, ymin = Delta14C_mean - Delta14C_sd,
                    ymax = Delta14C_mean + Delta14C_sd,
                    group = Horizon),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_C_sum, aes(y = Delta14C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1968.5,2024), expand = c(0,0),
                     breaks = seq(1969,2023,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), limits = c(-175,1000),
                     expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual("Modeled", label = c("Oi/Oe", "Oa/A", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oi/Oe", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) 

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
  scale_x_continuous("Year", limits = c(1968,2024), expand = c(0,0),
                     breaks = seq(1969,2023,10)) +
  scale_y_continuous(expression(paste("SOC stocks [g C m"^-2,"]")),
                     limits = c(500,5000), expand = c(0,0),
                     breaks = seq(500,5000,1500)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = "transparent")) +
  scale_color_manual("Modeled", label = c("Oi/Oe", "Oa/A", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oi/Oe", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon, ncol = 1) 

ggarrange(sens_all_14C_p, sens_all_C_p, common.legend = TRUE, widths = c(1.5,1),
          labels = c("a)", "b)"))
ggsave(file = paste0("./Output/HBEF_Figure2_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

### Plot predicted vs observed (mean + SD) for each Horizon and calculate RMSE
## 14C Data
model_14C_pred_obs <- sens_all_14C %>%
  right_join(HBEF_data_14C_C_sum)

model_14C_pred_obs_res_oie <- model_14C_pred_obs %>% 
  filter(Horizon == "Oi/Oe") %>% 
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
  filter(Horizon == "Oa/A") %>% 
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
  filter(Horizon == "Mineral") %>% 
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

# Function to plot data
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

# Plot data
oie_14C <- fun_pred_obs_14C(x = "Oi/Oe") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(0,475), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(0,525), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))

oa_14C <- fun_pred_obs_14C(x = "Oa/A") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-20,110), expand = c(0,0)) +
  scale_y_continuous("", limits = c(-20,230), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_14C <- fun_pred_obs_14C(x = "Mineral") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-30,-11), expand = c(0,0)) +
  scale_y_continuous("", limits = c(-110,20), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_14C, oa_14C, min_14C, nrow = 1, labels = c("a)", "b)", "c)"))
ggsave(file = paste0("./Output/HBEF_FigureS3_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)
## SOC data
model_C_pred_obs <- sens_all_C %>%
  right_join(HBEF_data_14C_C_sum) %>% 
  drop_na(C_mean)

#Cannot make correlation estimates for SOC stocks as the standard deviation is 0 for predicted values

# Calculate residuals and RMSE
model_C_pred_obs_res_oie <- model_C_pred_obs %>% 
  filter(Horizon == "Oi/Oe") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_oie %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_oie),
            #take the square root
            rmse = sqrt(msr))

model_C_pred_obs_res_oa <- model_C_pred_obs %>% 
  filter(Horizon == "Oa/A") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_oa %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_oa),
            #take the square root
            rmse = sqrt(msr))

model_C_pred_obs_res_min <- model_C_pred_obs %>% 
  filter(Horizon == "Mineral") %>% 
  #subtract the predicted value from the actual value (residual)
  mutate(res = Mean - C_mean,
         #Square each of the calculated residuals
         sq_res = (res)^2)

model_C_pred_obs_res_min %>% 
  #Calculate the mean of the squared differences
  summarise(msr = sum(sq_res)/length(model_C_pred_obs_res_min),
            #take the square root
            rmse = sqrt(msr))

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

oie_C <- fun_pred_obs_C(x = "Oi/Oe") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(1400,1735), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed SOC stocks [g C m"^-2,"]")),
                     limits = c(550,3750), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))

oa_C <- fun_pred_obs_C(x = "Oa/A") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(1700,2035), expand = c(0,0)) +
  scale_y_continuous("", limits = c(500,5200), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_C <- fun_pred_obs_C(x = "Mineral") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(2200,2365), expand = c(0,0)) +
  scale_y_continuous("", limits = c(1100,4100), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_C, oa_C, min_C, nrow = 1, labels = c("a)", "b)", "c)"))
ggsave(file = paste0("./Output/HBEF_FigureS4_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

#### Calculate transit time/ system age 
## see Stoner et al 2021 and Gonzalez-Sosa et al. 2024
propagation_fun <- function(data, num_iter, input_vector){
  pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
  
  n_iter <- num_iter
  iter <- numeric(n_iter)
  Iter_number <- numeric(n_iter)
  System_age <- numeric(n_iter)
  oie_age <- numeric(n_iter)
  oa_age <- numeric(n_iter)
  min_age <- numeric(n_iter)
  Transit_time <- numeric(n_iter)

    for (i in 1:n_iter){
        subset_par <- data.frame(k1 = sample(data$k1, 1),
                                 k2 = sample(data$k2, 1),
                                 k3 = sample(data$k3, 1),
                                 alpha21 = sample(data$alpha21, 1),
                                 alpha31 = sample(data$alpha32, 1))

    #-------------- A matrix and inputs
    ks <- subset_par[1:3]
    A <- -1 * diag(ks)

    ## Add transfer coefficients to A matrix:
    alpha_2_1 <- subset_par[4]
    alpha_3_2 <- subset_par[5]

    A[2,1] <- as.numeric(alpha_2_1*subset_par[1])
    A[3,2] <- as.numeric(alpha_3_2*subset_par[2])

    u <- matrix((input_vector), ncol = 1)

    #---------------- Age and transit time ----------------------

    Sys_age <- systemAge(A = A, u = u)
    Trans_time <- transitTime(A = A, u = u)

    Iter_number[i] <- i
    System_age[i] <- as.numeric(Sys_age$meanSystemAge)
    oie_age[i] <- as.numeric(Sys_age$meanPoolAge[1])
    oa_age[i] <- as.numeric(Sys_age$meanPoolAge[2])
    min_age[i] <- as.numeric(Sys_age$meanPoolAge[3])
    Transit_time[i] <- as.numeric(Trans_time$meanTransitTime)


    age_results <- as.data.frame(cbind(
      Iter_number, System_age, oie_age, oa_age, min_age, Transit_time
    ))

    setTxtProgressBar(pb, i)

  }

  return(age_results)
}

### NO NEED TO RUN THE FOLLOWING CODE (OBJECT HAS BEEN SAVED) ###
## Un-comment to run code ##
# exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
# pars_df <- as.data.frame(tpsMcmcFits$pars[-(1:1000), ])

# age_transit_dist <- propagation_fun(data = pars_df, num_iter = 10000, 
#                                     input_vector = c(210,0,0))

# write_csv(age_transit_dist,
#           file = paste0("./Output/HBEF_3ps_steady_long_MeanAge_Transit_Distribution_", lag_time, "_",
#                         Sys.Date(), ".csv"))

age_transit_dist <- read_csv("./Output/HBEF_3ps_steady_long_MeanAge_Transit_Distribution_3_2024-12-31.csv")

## Calculate system age and transit time based on median par fits
u <- matrix(c(210,0,0))
pars <- c(median(tpsMcmcFits$pars[-(1:1000),1]),
          median(tpsMcmcFits$pars[-(1:1000),2]),
          median(tpsMcmcFits$pars[-(1:1000),3]),
          median(tpsMcmcFits$pars[-(1:1000),4]),
          median(tpsMcmcFits$pars[-(1:1000),5]))

k <- pars[1:3]
A <- -1 * diag(k)
A[2,1] <- k[1]*pars[4]
A[3,2] <- k[2]*pars[5]

ua <- seq(0,1600)

# System Age
SA <- systemAge(A = A, u = u, a = ua)
SA$meanSystemAge
SA$meanPoolAge
SA$quantilesSystemAge

# Transit age
TT <- transitTime(A = A, u = u, a = ua)
TT$meanTransitTime

sa_tt <- data.frame(
  age = ua,
  system_age = SA$systemAgeDensity,
  oie = SA$poolAgeDensity[,1],
  oa = SA$poolAgeDensity[,2],
  min = SA$poolAgeDensity[,3],
  transit_time = TT$transitTimeDensity
)

sum_age_tt_fun <- function(vector){
  data.frame(mean_age = mean(vector),
             mean_age_sd = sd(vector),
             q1 = as.numeric(quantile(vector, probs = c(0.05))),
             median = median(vector),
             q3 = as.numeric(quantile(vector, probs = c(0.95))),
             ci_low = as.numeric(ci(vector))[2],
             ci_high = as.numeric(ci(vector))[3] 
             )
}

SA_TT_sum <- data.frame(
  rbind(
    cbind(Variable = "system_age", sum_age_tt_fun(age_transit_dist$System_age)),
    cbind(Variable = "Oi/Oe", sum_age_tt_fun(age_transit_dist$oie_age)),
    cbind(Variable = "Oa/A", sum_age_tt_fun(age_transit_dist$oa_age)),
    cbind(Variable = "Mineral", sum_age_tt_fun(age_transit_dist$min_age)),
    cbind(Variable = "transit_time", sum_age_tt_fun(age_transit_dist$Transit_time))
  )
)

## Plot results
pool_age_dens <- sa_tt %>% 
  dplyr::select(age, oie, oa, min) %>% 
  pivot_longer(!age, values_to = "age_dens", names_to = "Horizon") %>% 
  mutate(Horizon = case_when(
    Horizon == "oie" ~"Oi/Oe",
    Horizon == "oa" ~ "Oa/A",
    Horizon == "min" ~ "Mineral"
  ))

pool_age_dens$Horizon <- factor(pool_age_dens$Horizon, 
                                levels = c("Oi/Oe", "Oa/A", "Mineral"),
                                ordered = TRUE)

pool_age_fun <- function(x, ci_low, ci_high, mean_age){
  pool_age_dens %>%
    filter(Horizon == x) %>%
    ggplot(aes(y = age_dens, color = Horizon, x = age)) +
    annotate("rect", xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf,
             alpha = 0.1) +
    geom_line(linewidth = 1) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none") +
    facet_wrap(~Horizon, scales = "free") +
    geom_vline(xintercept = mean_age, linetype = "dotdash")
}

oie_age_p <- pool_age_fun(x = "Oi/Oe", mean_age = SA_TT_sum$mean_age[2],
                          ci_low = SA_TT_sum$ci_low[2], ci_high = SA_TT_sum$ci_high[2]) +
  scale_x_continuous("", expand = c(0,0)) +
  coord_cartesian(xlim = c(0,50)) +
  scale_color_manual(values = c("#33a02c")) +
  theme(plot.margin = margin(5,10,5,5, "pt")) +
  annotate(geom = "text", x = 20, y = 0.12,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[2],0), " yrs")) +
  scale_y_continuous("Density function", expand = c(0,0), breaks = seq(0,0.14,0.04),
                     limits = c(0,0.14))

oa_age_p <- pool_age_fun(x = "Oa/A", mean_age = SA_TT_sum$mean_age[3],
                         ci_low = SA_TT_sum$ci_low[3], ci_high = SA_TT_sum$ci_high[3]) +
  scale_x_continuous("Pool age [yr]", expand = c(0,0)) +
  coord_cartesian(xlim = c(0,650)) +
  scale_color_manual(values = c("#b2df8a")) +
  theme(plot.margin = margin(5,10,5,5, "pt")) +
  annotate(geom = "text", x = 300, y = 0.0085,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[3],0), " yrs")) +
  scale_y_continuous("", expand = c(0,0), breaks = seq(0,0.01, 0.0025),
                     limits = c(0,0.01))

min_age_p <- pool_age_fun(x = "Mineral", mean_age = SA_TT_sum$mean_age[4],
                          ci_low = SA_TT_sum$ci_low[4], ci_high = SA_TT_sum$ci_high[4]) +
  scale_x_continuous("", expand = c(0,0)) +
  coord_cartesian(xlim = c(0,1600)) +
  scale_color_manual(values = c("#a6cee3")) +
  theme(plot.margin = margin(5,15,5,0, "pt")) +
  annotate(geom = "text", x = 800, y = 0.0025,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[4],0), " yrs")) +
  scale_y_continuous("", expand = c(0,0), breaks = seq(0,0.003, 0.001),
                     limits = c(0,0.003))

ggarrange(oie_age_p, oa_age_p, min_age_p, nrow = 1, labels = c("a)", "b)", "c)"))
ggsave(file = paste0("./Output/HBEF_FigureS5_",
                     Sys.Date(), ".jpeg"), width = 12, height = 5)

# Calculate and plot how much C is cycling on different timescales
pool_age_perc <- pool_age_dens %>%
  group_by(Horizon) %>%
  mutate(age_perc = (age_dens/sum(age_dens))*100,
         perc_sum = cumsum(age_perc))

pool_age_cum_fun <- function(x, ci_low, ci_high, mean_age){
  pool_age_perc %>%
    filter(Horizon == x) %>%
    ggplot(aes(x = age, y = perc_sum, color = Horizon)) +
    geom_vline(xintercept = mean_age, linetype = "dotdash") +
    annotate("rect", xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf,
             alpha = 0.1) +
    geom_path(linewidth = 2) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          panel.grid.major = element_line(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_rect(fill = "transparent")) +
    facet_wrap(~Horizon)
}

oie_age_cum_p <- pool_age_cum_fun(x = "Oi/Oe", mean_age = SA_TT_sum$mean_age[2],
                                  ci_low = SA_TT_sum$ci_low[2], 
                                  ci_high = SA_TT_sum$ci_high[2]) +
  scale_y_continuous("Cumulative distribution [%]", limits = c(0,100), expand = c(0,0)) +
  scale_x_continuous("", expand = c(0,0)) +
  coord_cartesian(xlim = c(0,50)) +
  scale_color_manual(values = c("#33a02c")) +
  annotate(geom = "text", x = 20, y = 10,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[2],0), " yrs"))

oa_age_cum_p <- pool_age_cum_fun(x = "Oa/A", mean_age = SA_TT_sum$mean_age[3],
                                 ci_low = SA_TT_sum$ci_low[3], 
                                 ci_high = SA_TT_sum$ci_high[3]) +
  scale_y_continuous("", expand = c(0,0), limits = c(0,100)) +
  scale_x_continuous("Age [yr]", expand = c(0,0), breaks = seq(0,1600,100)) +
  coord_cartesian(xlim = c(0,650)) +
  scale_color_manual(values = c("#b2df8a")) +
  annotate(geom = "text", x = 300, y = 10,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[3],0), " yrs"))

min_age_cum_p <- pool_age_cum_fun(x = "Mineral", mean_age = SA_TT_sum$mean_age[4],
                                  ci_low = SA_TT_sum$ci_low[4], 
                                  ci_high = SA_TT_sum$ci_high[4]) +
  scale_y_continuous("", expand = c(0,0), limits = c(0,100)) +
  scale_x_continuous("", expand = c(0,0), breaks = seq(0,1600,250)) +
  coord_cartesian(xlim = c(0,1600)) +
  scale_color_manual(values = c("#a6cee3")) +
  annotate(geom = "text", x = 800, y = 10,
           label = paste0("Mean age = ", round(SA_TT_sum$mean_age[4],0), " yrs"))

ggarrange(oie_age_cum_p, oa_age_cum_p, min_age_cum_p, nrow = 1,
          labels = c("a)", "b)", "c)"))
ggsave(file = paste0("./Output/HBEF_Figure3_",
                     Sys.Date(), ".jpeg"), width = 12, height = 5)

## Calculate fluxes and stocks for all pools and transfer rates (Figure 4)

# Parameter estimates for each model run
pf_3ps <- tpsMcmcFits$pars %>% 
  as.data.frame() %>% 
  # add estimated C stocks for each pool
  mutate(C1 = mean(sens_all_C$q50[sens_all_C$Horizon == "oie"]),
         C2 = mean(sens_all_C$q50[sens_all_C$Horizon == "oa"]),
         C3 = mean(sens_all_C$q50[sens_all_C$Horizon == "min"])) %>% 
  mutate(
    #Output flux from pool 1
    out1 = (1 - alpha21) * k1 * C1,
    #transfer rate from pool 1 to pool 2
    trans1 = alpha21 * k1 * C1,
    out2 = (1- alpha32) * k2 * C2,
    trans2 = alpha32 * k2 * C2,
    out3 = k3 * C3
  )

pf_3ps %>% 
  summarise(median_out1 = round(median(out1), 0),
            lci_out1 = round(quantile(out1, 0.05), 0),
            uci_out1 = round(quantile(out1, 0.95), 0),
            median_trans1 = round(median(trans1), 0),
            lci_trans1 = round(quantile(trans1, 0.05), 0),
            uci_trans1 = round(quantile(trans1, 0.95), 0),
            median_out2 = round(median(out2), 0),
            lci_out2 = round(quantile(out2, 0.05), 0),
            uci_out2 = round(quantile(out2, 0.95), 0),
            median_trans2 = round(median(trans2), 0),
            lci_trans2 = round(quantile(trans2, 0.05), 0),
            uci_trans2 = round(quantile(trans2, 0.95), 0),
            median_out3 = round(median(out3), 0),
            lci_out3 = round(quantile(out3, 0.05), 0),
            uci_out3 = round(quantile(out3, 0.95), 0))
