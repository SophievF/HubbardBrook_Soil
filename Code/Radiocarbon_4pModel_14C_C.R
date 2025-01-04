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

# Make atmospheric 14C forecast (Hua 2021 only goes until 2019
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

## Define initial values for model 
# Define time interval for model
years <- seq(1997, 2023, by = 0.5)

## initial Delta14C in each pool (based on steady-state 3p model)
sens_14c <- read_csv("./Output/HBEF_3ps_steady_long_sens_14C_3_2024-12-30.csv") %>% 
  filter(Year == 1997)
init14C <- ConstFc(as.numeric(c(250, 150, sens_14c[2,8], sens_14c[3,8])),
                   format = "Delta14C")

sens_c <- read_csv("./Output/HBEF_3ps_steady_long_sens_C_3_2024-12-31.csv") %>% 
  filter(Year == 1997)
#Split between Oi and Oe based on data from Chris Johnson
C0 <- c(sens_c[1,8]*0.25, sens_c[1,8]*0.75, sens_c[2,8], sens_c[3,8])
C0 <- as.numeric(C0)
rm(sens_14c, sens_c)

# lag-time before C enters soils: based on communication with Josh
lag_time <- 3

# Number of model iterations
itr <- 15000

# Initial k and a values (k = stocks/flux; modified from Fahey et al 2005)
init_pars <- c(k1 = 1/2, k2 = 1/7, k3 = 1/19, k4 = 1/59, 
               alpha21 = 180/(180 + 30), alpha32 = 100/(100 + 80),
               alpha43 = 39/(39 + 61)) 

## Set-up model (four pools in series)
FourPSeriesModel_fun <- function(pars){
  # Define and construct matrix A
  ks <- pars[1:4]
  A <- -1 * diag(ks)
  A[2,1] <- pars[5]*pars[1]
  A[3,2] <- pars[6]*pars[2]
  A[4,3] <- pars[7]*pars[3]
  
  mod = Model_14(
    # time interval for model to find solutions
    t = years,
    # matrix that contains k and a values
    A = A, 
    ivList = C0,
    # Need to define inputs into each pool separately; only Oi recieves inputs
    inputFluxes = BoundInFluxes(map = function(t){matrix(c(210,0,0,0), 
                                                       ncol = 1, nrow = 4)}),
    # Initial 14C values
    initialValF = init14C,
    # atmospheric 14C for each time point, including lag-time
    inputFc = BoundFc(map = NHZone2_2023,  format = "Delta14C", lag = lag_time)
  )
  # get 14C for each pool and each time point
  res_14C = getF14(mod)
  # get C for each pool and each time point
  res_C = getC(mod)
  # return dataframe with years and 14C and C values for each pool
  return(data.frame(time = years,
                    # Calculate bulk value for Oie (C-weighted)
                    oie_14C = (res_14C[,1]*(res_C[,1]/rowSums(res_C[,1:2]))) + 
                      (res_14C[,2]*(res_C[,2]/rowSums(res_C[,1:2]))),
                    oie_C = rowSums(res_C[,1:2]),
                    oi_14C = res_14C[,1],
                    oi_C = res_C[,1],
                    oe_14C = res_14C[,2],
                    oe_C = res_C[,2],
                    oa_14C = res_14C[,3],
                    oa_C = res_C[,3],
                    min_14C = res_14C[,4],
                    min_C = res_C[,4]))
}

# Cost function to constrain model with observational data
fpsCost <- function(pars){
  funccall = FourPSeriesModel_fun(pars)
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
fpsModelFit <- FME::modFit(f = fpsCost, p = init_pars, method = "Marq", 
                           upper = c(3, rep(1,6)), lower = rep(0,7)) 

## Model evaluation
# sum squared residuals
fpsModelFit$ssr

# mean squared residuals
fpsModelFit$ms

# AIC
(2*length(fpsModelFit$par))-(2*log(fpsModelFit$ms))

# Mean squared residuals per variable/horizon
sqrt(fpsModelFit$var_ms)

model_summary <- data.frame(ssr = fpsModelFit$ssr,
                            msr = fpsModelFit$ms,
                            aic = (2*length(fpsModelFit$par))-(2*log(fpsModelFit$ms)))

write.csv(model_summary, row.names = TRUE, quote = FALSE,
          file = paste0("./Output/HBEF_4ps_short_14C_C_summary_stats_", 
                        Sys.Date(), ".csv"))

## Perform Markov Chain Monte Carlo simulation
fpsVar <- fpsModelFit$var_ms_unweighted

fpsMcmcFits <- FME::modMCMC(f = fpsCost, p = fpsModelFit$par, niter = itr, ntrydr = 5,
                            updatecov = 50, var0 = fpsVar, upper = c(3, rep(1,6)),
                            lower = rep(0,7)) #Create a new object to record fit stats

# Fit model with mean parameter values
fpsModelOutput <- FourPSeriesModel_fun(pars = as.numeric(summary(fpsMcmcFits)[1,1:7]))

write.csv(summary(fpsMcmcFits), 
          file = paste0("./Output/HBEF_4ps_short_14C_C_summary_", 
                        Sys.Date(), ".csv"))

# Check for convergence: if model is converged, there should be no visible drift
jpeg(paste0("./Output/HBEF_4ps_short_14C_C_converg_", 
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
plot(fpsMcmcFits)
dev.off()

jpeg(paste0("./Output/HBEF_4ps_short_14C_C_pairs_", 
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
pairs(fpsMcmcFits)
dev.off()

#### Uncertainty analysis ####
# exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
pars <- fpsMcmcFits$pars[-(1:1000), ]

# Number of times the model will be run
num <- 1000

### Estimate uncertainties for each horizon
## Radiocarbon
sens_oi_14C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                  sensvar = "oi_14C")) %>% 
  mutate(Horizon = "Oi")

sens_oe_14C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                 sensvar = c("oe_14C"))) %>% 
  mutate(Horizon = "Oe")

sens_oie_14C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_14C"))) %>% 
  mutate(Horizon = "Oi/Oe")

sens_oa_14C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                             sensvar = c("oa_14C"))) %>% 
  mutate(Horizon = "Oa/A")

sens_min_14C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                              sensvar = c("min_14C"))) %>% 
  mutate(Horizon = "Mineral")

# Combine results for each horizon into one data frame and save as csv file
sens_all_14C <- rbind(sens_oi_14C, sens_oe_14C, sens_oie_14C, sens_oa_14C, sens_min_14C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_14C, 
          file = paste0("./Output/HBEF_4ps_short_SensitivityAnalysis_14C_", 
                        Sys.Date(), ".csv"))

## Carbon
sens_oi_C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                sensvar = c("oi_C"))) %>% 
  mutate(Horizon = "Oi")

sens_oe_C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                               sensvar = c("oe_C"))) %>% 
  mutate(Horizon = "Oe")

sens_oie_C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                  sensvar = c("oie_C"))) %>% 
  mutate(Horizon = "Oi/Oe")

sens_oa_C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                 sensvar = c("oa_C"))) %>% 
  mutate(Horizon = "Oa/A")

sens_min_C <- summary(sensRange(num = num, func = FourPSeriesModel_fun, parInput = pars, 
                                  sensvar = c("min_C"))) %>% 
  mutate(Horizon = "Mineral")

# Combine results for each horizon into one data frame and save as csv file
sens_all_C <- rbind(sens_oi_C, sens_oe_C, sens_oie_C, sens_oa_C, sens_min_C) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1968)

write_csv(sens_all_C, 
          file = paste0("./Output/HBEF_4ps_short_SensitivityAnalysis_C_",
                        Sys.Date(), ".csv"))

### Save all model outputs
save(fpsModelFit, fpsMcmcFits, fpsModelOutput, sens_all_C, sens_all_14C,
     file = paste0("./Output/HBEF_4ps_short_14C_C_", Sys.Date(), ".Rdata"))

##### Plot and evaluate model results ######
# Load model results (from code above)
# load("./Output/HBEF_4ps_short_14C_C_3_2025-01-04.RData")

# Create long dataframe
fpsModelOutput_df <- fpsModelOutput %>% 
  # filter(time > 1945) %>% 
  pivot_longer(!time,
               cols_vary = "slowest",
               names_to = c("Horizon", ".value"),
               names_pattern = "(.*)_(.*)") %>% 
  rename(Delta14C = `14C`,
         SOC_Stock = C,
         Year = time) %>% 
  mutate(Horizon = case_when(
    Horizon == "oi" ~ "Oi",
    Horizon == "oe" ~ "Oe",
    Horizon == "oie" ~ "Oi/Oe",
    Horizon == "oa" ~ "Oa/A",
    Horizon == "min" ~ "Mineral"))

fpsModelOutput_df$Horizon <- factor(fpsModelOutput_df$Horizon,
                                    levels = c("Oi", "Oe", "Oi/Oe", "Oa/A", "Mineral"),
                                    ordered = TRUE)

# Check best parameter values
fpsMcmcFits$bestpar

# Check summary statistics for estimated parameter values
summary(fpsMcmcFits)

## Prepare data for plotting
# Summarise HBEF data for plotting and merging with model results
HBEF_data_14C_C_sum <- HBEF_all %>% 
  filter(Year != 2020) %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C, na.rm = TRUE),
            Delta14C_sd = sd(Delta14C, na.rm = TRUE),
            C_mean = mean(SOC_g_m2, na.rm = TRUE),
            C_sd = sd(SOC_g_m2, na.rm = TRUE)) 

sens_all_14C$Horizon <- factor(sens_all_14C$Horizon,
                               levels = c("Oi", "Oe", "Oi/Oe", "Oa/A", "Mineral"))

sens_all_C$Horizon <- factor(sens_all_C$Horizon,
                               levels = c("Oi", "Oe", "Oi/Oe", "Oa/A", "Mineral"))

## Plot 14C results
sens_all_14C_p <- sens_all_14C %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, group = Horizon), alpha = 0.2) +
  geom_line(aes(y = q50, color = Horizon, linetype = Horizon), linewidth = 1) +
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
  scale_linetype_manual("Modeled", label = c("Oi", "Oe", "Oi/Oe", "Oa/A", "Mineral"),
                        values = c(2, 3, 1, 1, 1)) +
  scale_color_manual("Modeled", label = c("Oi", "Oe", "Oi/Oe", "Oa/A", "Mineral"),
                     values = c("#33a02c", "#33a02c", "#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", label = c("Oi/Oe", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) 

ggsave(file = paste0("./Output/HBEF_4ps_short_14C_Sensitivity_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6) 

## Plot SOC results
sens_all_C_p <- sens_all_C %>% 
  ggplot(aes(x = Year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, group = Horizon), alpha = 0.2) +
  geom_line(aes(y = q50, color = Horizon), linewidth = 1) +
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
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_linetype_manual("Modeled", values = c(2, 3, 1, 1, 1)) +
  scale_color_manual("Modeled", values = c("#33a02c", "#33a02c", "#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon)

ggsave(file = paste0("./Output/HBEF_4ps_short_C_Sensitivity_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6)

ggarrange(sens_all_14C_p, sens_all_C_p, common.legend = TRUE)
ggsave(file = paste0("./Output/HBEF_4ps_short_14_C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)

## Plot SOC Oie results only
sens_all_C %>% 
  filter(Horizon == "Oi/Oe") %>% 
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
  scale_x_continuous("Year", expand = c(0,0), limits = c(1997.9,2023.1)) +
  scale_y_continuous(expression(paste("SOC stocks [g C m"^-2,"]")), 
                     limits = c(750,1750), expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual("Modeled", values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured", values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  facet_wrap(~Horizon) 
ggsave(file = paste0("./Output/HBEF_4ps_short_14_C_Sensitivity_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 7, height = 5)

#Plot predicted vs observed (mean + SD) for each Horizon and compute residuals
model_14C_pred_obs <- sens_all_14C %>%
  filter(Horizon == "oie"| Horizon == "oa"| Horizon == "min") %>% 
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
                     limits = c(0,410), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(0,525), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))

oa_14C <- fun_pred_obs_14C(x = "oa") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(89,105), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(-20,150), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_14C <- fun_pred_obs_14C(x = "min") +
  scale_x_continuous(expression(paste("Predicted ", Delta^14, "C [‰]")),
                     limits = c(-29,-14), expand = c(0,0)) +
  scale_y_continuous(expression(paste("Observed ", Delta^14, "C [‰]")),
                     limits = c(-110,25), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_14C, oa_14C, min_14C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_4ps_short_14C_Obs_Pred_", lag_time, "_",
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
                     limits = c(1170,1610), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(750,1750), expand = c(0,0)) +
  scale_color_manual(values = c("#33a02c"))
ggsave(file = paste0("./Output/HBEF_4ps_C_Obs_Pred_oie_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 5, height = 6)

oa_C <- fun_pred_obs_C(x = "oa") +
  scale_x_continuous("Predicted SOC stocks [g/m2]",
                     limits = c(1490,1890), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(750,3250), expand = c(0,0)) +
  scale_color_manual(values = c("#b2df8a"))

min_C <- fun_pred_obs_C(x = "min") +
  scale_x_continuous("Predicted SOC stocks [g/m2]",
                     limits = c(2185,2430), expand = c(0,0)) +
  scale_y_continuous("Observed SOC stocks [g/m2]",
                     limits = c(1000,4000), expand = c(0,0)) +
  scale_color_manual(values = c("#a6cee3"))

ggarrange(oie_C, oa_C, min_C, nrow = 1)
ggsave(file = paste0("./Output/HBEF_4ps_short_C_Obs_Pred_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 12, height = 6)
