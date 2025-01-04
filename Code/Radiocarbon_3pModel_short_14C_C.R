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
HBEF_all <- read_csv("./Data/HBEF_ModelingDatasets_all_2025-01-02.csv") %>% 
  #only select dataset 2 (1998 to 2023)
  filter(DataSource == "Dataset2")

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
                           mean(oie_data_14C$sd, na.rm = TRUE))

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
                           mean(oie_data_C$sd, na.rm = TRUE))

#### Model set-up ####

# Define lambda
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
years <- seq(1997, 2023, by = 0.5)

## initial Delta14C in each pool (based on steady-state 3p model)
sens_14c <- read_csv("./Output/HBEF_3ps_steady_long_sens_14C_3_2024-12-30.csv") %>% 
  filter(Year == 1997)
init14C <- c(sens_14c[1,8], sens_14c[2,8], sens_14c[3,8])
init14C <- as.numeric(init14C)

## initial SOC in each pool (based on steady-state 3p model)
sens_c <- read_csv("./Output/HBEF_3ps_steady_long_sens_C_3_2024-12-31.csv") %>% 
  filter(Year == 1997)
C0 <- c(sens_c[1,8], sens_c[2,8], sens_c[3,8])
C0 <- as.numeric(C0)
rm(sens_14c, sens_c)

# lag-time before C enters soils: based on communication with Josh
lag_time <- 3

# Number of model iterations
itr <- 15000

# C inputs (based on values from Fahey et al. 2005)
In <- 210

# Initial k and a values (k = stocks/flux; modified from Fahey et al 2005)
init_pars <- c(k1 = 1/8, k2 = 1/19, k3 = 1/59, 
               alpha21 = 100/(100 + 110), alpha32 = 39/(39 + 61))

## Set-up model (three pools in series)
ThreePSeriesModel_fun <- function(pars){
  mod = ThreepSeriesModel14(
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
  res_14C = getF14(mod)
  # get C for each pool and each time point
  res_C = getC(mod)
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
  cost1 = modCost(model = funccall, obs = oie_data_14C, err = "sd")
  cost2 = modCost(model = funccall, obs = oa_data_14C, err = "sd", cost = cost1)
  cost3 = modCost(model = funccall, obs = min_data_14C, err = "sd", cost = cost2)
  cost4 = modCost(model = funccall, obs = oie_data_C, err = "sd", cost = cost3)
  cost5 = modCost(model = funccall, obs = oa_data_C, err = "sd", cost = cost4)
  cost6 = modCost(model = funccall, obs = min_data_C, err = "sd", cost = cost5)
  return(cost6)
}

#### NO NEED TO RUN THE FOLLOWING CODE (OBJECTS HAVE BEEN SAVED) #####
## Un-comment the code if you wish to run it ##

## Fitting the model
tpsModelShortFit <- FME::modFit(f = tpsCost, p = init_pars, method = "Marq", 
                                upper = c(3, rep(1,4)), lower = rep(0,5)) 

## Model evaluation
# sum squared residuals
tpsModelShortFit$ssr

# mean squared residuals
tpsModelShortFit$ms

# AIC
(2*length(tpsModelShortFit$par))-(2*log(tpsModelShortFit$ms))

# Mean squared residuals per variable/horizon
sqrt(tpsModelShortFit$var_ms)

model_summary <- data.frame(ssr = tpsModelShortFit$ssr,
                            msr = tpsModelShortFit$ms,
                            aic = (2*length(tpsModelShortFit$par))-(2*log(tpsModelShortFit$ms)))

write.csv(model_summary, row.names = TRUE, quote = FALSE,
          file = paste0("./Output/HBEF_3ps_short_14C_C_summary_stats_", lag_time, "_",
                        Sys.Date(), ".csv"))

## Perform Markov Chain Monte Carlo simulation
tpsVar <- tpsModelShortFit$var_ms_unweighted

tpsShortMcmcFits <- FME::modMCMC(f = tpsCost, p = tpsModelShortFit$par, niter = itr, ntrydr = 5,
                                 updatecov = 50, var0 = tpsVar, upper = c(3, rep(1,4)),
                                 lower = rep(0,5)) #Create a new object to record fit stats

# Fit model with mean parameter values
tpsModelShortOutput <- ThreePSeriesModel_fun(pars = as.numeric(summary(tpsShortMcmcFits)[1,1:5]))

write.csv(summary(tpsShortMcmcFits), 
          file = paste0("./Output/HBEF_3ps_short_14C_C_summary_", lag_time, "_",
                        Sys.Date(), ".csv"))

# Check for convergence: if model is converged, there should be no visible drift
jpeg(paste0("./Output/HBEF_3ps_short_14C_C_converg_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
plot(tpsShortMcmcFits)
dev.off()

jpeg(paste0("./Output/HBEF_3ps_short_14C_C_pairs_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
pairs(tpsShortMcmcFits)
dev.off()

#### Uncertainty analysis ####
# exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
pars <- tpsShortMcmcFits$pars[-(1:1000), ]

# Number of times the model will be run
num <- 1000

### Estimate uncertainties for each horizon
## Radiocarbon
sens_oie_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_14C"))) %>% 
  mutate(Horizon = "Oi/Oe")

sens_oa_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                             sensvar = c("oa_14C"))) %>% 
  mutate(Horizon = "Oa/A")

sens_min_14C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                              sensvar = c("min_14C"))) %>% 
  mutate(Horizon = "Mineral")

# Combine results for each horizon into one data frame and save as csv file
sens_all_14C <- rbind(sens_oie_14C, sens_oa_14C, sens_min_14C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_14C, 
          file = paste0("./Output/HBEF_3ps_short_sens_14C_", lag_time, "_",
                        Sys.Date(), ".csv"))

## Carbon
sens_oie_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_C"))) %>% 
  mutate(Horizon = "Oi/Oe")

sens_oa_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                 sensvar = c("oa_C"))) %>% 
  mutate(Horizon = "Oa/A")

sens_min_C <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                                  sensvar = c("min_C"))) %>% 
  mutate(Horizon = "Mineral")

# Combine results for each horizon into one data frame and save as csv file
sens_all_C <- rbind(sens_oie_C, sens_oa_C, sens_min_C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_C, 
          file = paste0("./Output/HBEF_3ps_short_sens_C_", lag_time, "_",
                        Sys.Date(), ".csv"))

### Save all model outputs
save(tpsModelShortFit, tpsShortMcmcFits, tpsModelShortOutput, sens_all_C, sens_all_14C, 
     file = paste0("./Output/HBEF_3ps_short_14C_C_", lag_time, "_", Sys.Date(), ".Rdata"))

##### Plot and evaluate model results ######
# Load model results (from code above)
# load("./Output/HBEF_3ps_short_14C_C_3_2025-01-04.RData")

# Create long dataframe
tpsModelOutput_df <- tpsModelShortOutput %>% 
  # filter(time > 1945) %>% 
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
tpsShortMcmcFits$bestpar

# Check summary statistics for estimated parameter values
summary(tpsShortMcmcFits)[,1:5]

## Prepare data for plotting
# Summarise HBEF data for plotting and merging with model results
HBEF_data_14C_C_sum <- HBEF_all %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C, na.rm = TRUE),
            Delta14C_sd = sd(Delta14C, na.rm = TRUE),
            C_mean = mean(SOC_g_m2, na.rm = TRUE),
            C_sd = sd(SOC_g_m2, na.rm = TRUE))

sens_all_14C$Horizon <- factor(sens_all_14C$Horizon,
                               levels = c("Oi/Oe", "Oa/A", "Mineral"),
                               ordered = TRUE)

sens_all_C$Horizon <- factor(sens_all_C$Horizon,
                             levels = c("Oi/Oe", "Oa/A", "Mineral"),
                             ordered = TRUE)

## Plot 14C results
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
  scale_color_manual("Modeled", label = c("Oi/Oe", "Oa/A", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oi/Oe", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) 

ggsave(file = paste0("./Output/HBEF_3ps_short_14C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6) 

## Plot SOC results
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
                     breaks = seq(1969,2023,10)) +
  scale_y_continuous(expression(paste("SOC stocks [g C m"^-2,"]")),
                     limits = c(500,5000), expand = c(0,0),
                     breaks = seq(500,5000,1500)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled", label = c("Oi/Oe", "Oa/A", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oi/Oe", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon, ncol = 1)

ggsave(file = paste0("./Output/HBEF_3ps_short_C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6)

ggarrange(sens_all_14C_p, sens_all_C_p, common.legend = TRUE, widths = c(1.5,1))
ggsave(file = paste0("./Output/HBEF_3ps_short_14C_C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

## Plot SOC Oie results only
sens_all_C %>% 
  filter(Horizon == "Oi/Oe") %>% 
  ggplot(aes(x = Year)) +
  #add regression line based on measured values
  geom_smooth(data = HBEF_data_14C_C_sum %>% 
                filter(Horizon == "Oi/Oe"), aes(y = C_mean),
              method = "lm", color = "black", linetype = "dashed") +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = Horizon), alpha = 0.4) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_C_sum %>% 
                  filter(Horizon == "Oi/Oe"),
                aes(y = C_mean, ymin = C_mean - C_sd,
                    ymax = C_mean + C_sd),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_C_sum %>% 
               filter(Horizon == "Oi/Oe"), aes(y = C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", expand = c(0,0), limits = c(1997.9,2023.1)) +
  scale_y_continuous(expression(paste("SOC stocks [g C m"^-2,"]")), 
                     limits = c(750,1750), expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = c("#33a02c")) +
  scale_fill_manual(values = c("#33a02c")) +
  facet_wrap(~Horizon) 

ggsave(file = paste0("./Output/HBEF_3ps_short_14C_C_Sensitivity_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 7, height = 5)

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
                     limits = c(0,225), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(0,225), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))

oa_14C <- fun_pred_obs_14C(x = "Oa/A") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(75,90), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(0,130), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_14C <- fun_pred_obs_14C(x = "Mineral") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-28,-14), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(-110,20), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_14C, oa_14C, min_14C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_3ps_short_14C_Obs_Pred_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

## SOC data
model_C_pred_obs <- sens_all_C %>%
  right_join(HBEF_data_14C_C_sum) %>% 
  drop_na(C_mean)

# calculate residuals and RMSE
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

cor(y = model_C_pred_obs_res_oie$C_mean,
    x = model_C_pred_obs_res_oie$Mean,
    method = "pearson")

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

cor(y = model_C_pred_obs_res_oa$C_mean,
    x = model_C_pred_obs_res_oa$Mean,
    method = "pearson")

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

cor(y = model_C_pred_obs_res_min$C_mean,
    x = model_C_pred_obs_res_min$Mean,
    method = "pearson")

# Function to plot data
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

# Plot data
oie_C <- fun_pred_obs_C(x = "Oi/Oe") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(1170,1610), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed SOC stocks [g C m"^-2,"]")),
                     limits = c(750,1750), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))
ggsave(file = paste0("./Output/HBEF_3ps_short_C_Obs_Pred_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 5, height = 6)

oa_C <- fun_pred_obs_C(x = "Oa/A") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(1490,1905), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed SOC stocks [g C m"^-2,"]")),
                     limits = c(750,3200), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_C <- fun_pred_obs_C(x = "Mineral") +
  scale_x_continuous(expression(paste("Predicted SOC stocks [g C m"^-2,"]")),
                     limits = c(2185,2430), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed SOC stocks [g C m"^-2,"]")),
                     limits = c(1000,4100), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_C, oa_C, min_C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_3ps_short_C_Obs_Pred_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

