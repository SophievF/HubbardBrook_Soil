## Hubbard Brook archived soil samples project ##
## Radiocarbon analysis - constraint model with 14C ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(SoilR)
library(forecast)
library(FME)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-09-20.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes))

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  group_by(Horizon) %>% 
  reframe(mean_14C = mean(Delta14C, na.rm = TRUE))

HBEF_data %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Oa/A") %>% 
  ggplot(aes(x = `Measured_%_C`, y = Delta14C)) +
  geom_point(aes(color = Plot)) +
  facet_wrap(vars(Horizon, Year)) +
  # geom_smooth(method = "lm") +
  theme_bw()

# Load 14C data from Charley Driscoll
LitterData <- read_csv("./Data/LitterData_Driscoll.csv")

#fm * exp(lambda * (-obs_date_y + 1950)) - 1) * 1000
lambda <- 0.00012097

litter_data <- LitterData %>% 
  filter(Watershed == 6) %>% 
  mutate(Delta14C = (F14C * exp(lambda * (-Year + 1950)) -1) * 1000) %>% 
  mutate(Elevation = case_when(
    Plot < 70 ~ "high",
    Plot > 152 ~ "low",
    TRUE ~ "mid"
  )) %>% 
  #match Horizon names
  mutate(Horizon = case_when(
    Horizon == "Oie" ~ "Oi/Oe",
    Horizon == "Oa" ~ "Oa/A"
  ))

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
years <- seq(-10000, 2025, by = 0.5)

# initial C stocks in each pool
C0 <- c(1394, 1484, 3170)

# initial Delta14C in each pool: 0 for all fast cycling pools; -23 for mineral (= average value over all years)
init14C <- c(0, 0, -23)

# lag-time before C enters soils: based on communication with Josh
lag_time <- 3

# Number of model iterations
itr <- 15000

# C inputs
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
  model_result = getF14(mod)
  return(data.frame(time = years, oie = model_result[,1],
                    oa = model_result[,2], min = model_result[,3]))
}

# model_df <- ThreePSeriesModel_fun(init_pars)
# model_df %>% filter(time == 1998)

# Summarize and merge data by horizon; remove roots for now
HBEF_all <- HBEF_data %>% 
  filter(Plot != "all fine roots") %>% 
  dplyr::select(Year, Horizon, Delta14C) %>% 
  rbind(litter_data %>% 
          dplyr::select(Year, Horizon, Delta14C))

min_data <- HBEF_all %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Mineral_0_10") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(min = mean(Delta14C),
            sd = sd(Delta14C)) %>%
  data.frame()

oa_data <- HBEF_all %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Oa/A") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(oa = mean(Delta14C),
            sd = sd(Delta14C)) %>% 
  data.frame()

oie_data <- HBEF_all %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Oi/Oe") %>% 
  rename(time = Year) %>% 
  group_by(time) %>% 
  summarise(oie = mean(Delta14C),
            sd = sd(Delta14C)) %>%
  # Replace NA for 1998 sd with ~26 which is the mean sd for all years in the Oie
  replace(is.na(.), 26.66504) %>% 
  data.frame()

tpsCost <- function(pars){
  funccall = ThreePSeriesModel_fun(pars)
  cost1 = modCost(model = funccall, obs = oie_data, err = "sd")
  cost2 = modCost(model = funccall, obs = oa_data, err = "sd", cost = cost1)
  cost3 = modCost(model = funccall, obs = min_data, err = "sd", cost = cost2)
  return(cost3)
}

init_pars <- c(k1 = 1/6, k2 = 1/14, k3 = 1/81, 
               alpha21 = 100/(100 + 110), alpha32 = 39/(39 + 61))

tpsModelFit <- FME::modFit(f = tpsCost, p = init_pars, method = "Marq", 
                           upper = c(3, rep(1,4)), lower = rep(0,5)) 

tpsVar <- tpsModelFit$var_ms_unweighted

tpsMcmcFits <- FME::modMCMC(f = tpsCost, p = tpsModelFit$par, niter = itr, ntrydr = 5,
                            updatecov = 50, var0 = tpsVar, upper = c(3, rep(1,4)),
                            lower = rep(0,5)) #Create a new object to record fit stats

# tpsMcmcFits <- FME::modMCMC(f = tpsCost, p = tpsModelFit$par, niter = itr, ntrydr = 5, 
#                             updatecov = 50, var0 = tpsVar, upper = c(rep(1,5)), 
#                             lower = rep(0,5)) #Create a new object to record fit stats

tpsModelOutput <- ThreePSeriesModel_fun(pars = as.numeric(summary(tpsMcmcFits)[1,1:6]))

# Save output
save(tpsMcmcFits, tpsModelOutput, 
     file = paste0("./Output/ThreePoolSeriesModel_", lag_time, "_", Sys.Date(), ".Rdata"))
write_csv(summary(tpsMcmcFits), 
          file = paste0("./Output/ThreePoolSeriesModel_summary_", lag_time, "_",
                        Sys.Date(), ".csv"))

# Create long dataframe
tpsModelOutput_df <- tpsModelOutput %>% 
  filter(time > 1945) %>% 
  pivot_longer(!time, values_to = "Delta14C", names_to = "Horizon") %>% 
  rename(Year = time)

tpsModelOutput_df$Horizon <- factor(tpsModelOutput_df$Horizon,
                                    levels = c("oie", "oa", "min"),
                                    ordered = TRUE)

# Summarise HBEF data
HBEF_data_14C_sum <- HBEF_all %>% 
  drop_na(Delta14C) %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C),
            Delta14C_sd = sd(Delta14C)) %>% 
  mutate(Horizon = case_when(
    Horizon == "Oi/Oe" ~ "oie",
    Horizon == "Oa/A" ~ "oa",
    Horizon == "Mineral_0_10" ~ "min"
  ))

HBEF_data_14C_sum$Horizon <- factor(HBEF_data_14C_sum$Horizon,
                                    levels = c("oie", "oa", "min"),
                                    ordered = TRUE)

tpsMcmcFits$bestpar

summary(tpsMcmcFits)

#Check for convergence: if model is converged, there should be no visible drift
jpeg(paste0("./Output/HBEFall_SteadyStateModel_tpsModelFit_14C_converg_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
plot(tpsMcmcFits)
dev.off()

jpeg(paste0("./Output/HBEFall_SteadyStateModel_tpsModelFit_14C_pairs_", lag_time, "_",
            Sys.Date(), ".jpeg"), width = 1550, height = 1000)
pairs(tpsMcmcFits)
dev.off()

## Plot measured and modeled data together  
NHZone2_2023 %>%  
  filter(Year > 1945) %>% 
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() +
  geom_line(data = tpsModelOutput_df,
            aes(color = Horizon), linewidth = 1) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_sum,
                aes(y = Delta14C_mean, ymin = Delta14C_mean - Delta14C_sd,
                    ymax = Delta14C_mean + Delta14C_sd,
                    group = Horizon),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_sum, aes(y = Delta14C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1955,2025), expand = c(0,0),
                     breaks = seq(1955,2025,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), limits = c(-175,1000),
                     expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled\nhorizon data", label = c("Oie", "Oa", "0-10 cm"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured\nhorizon data", label = c("Oie", "Oa/A", "0-10 cm"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3"))

ggsave(file = paste0("./Output/HBEFall_SteadyStateModel_tpsModelFit_14C_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6)

#### Uncertainty analysis
pars <- tpsMcmcFits$pars

num <- 1000

# sens_all <- summary(FME::sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars))

sens_oie <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                              sensvar = c("oie"))) %>% 
  mutate(Horizon = "oie")

sens_oa <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                             sensvar = c("oa"))) %>% 
  mutate(Horizon = "oa")

sens_min <- summary(sensRange(num = num, func = ThreePSeriesModel_fun, parInput = pars, 
                              sensvar = c("min"))) %>% 
  mutate(Horizon = "min")

sens_all <- rbind(sens_oie, sens_oa, sens_min) %>% 
  tibble() %>% 
  dplyr::rename(Year = x) %>% 
  filter(Year > 1945)

write_csv(sens_all , 
          file = paste0("./Output/ThreePoolSeriesModel_SensitivityAnalysis_", lag_time, "_",
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
  dplyr::filter(Year > 1945)

sens_all$Horizon <- factor(sens_all$Horizon,
                           levels = c("oie", "oa", "min"), ordered = TRUE)


sens_all %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean, color = Horizon), linewidth = 1) +
  geom_ribbon(aes(ymin = Mean - Sd, ymax = Mean + Sd, fill = Horizon), alpha = 0.4) +
  geom_line(data = NHZone2_2023 %>% 
              filter(Year > 1945), aes(y = Delta14C)) +
  # Add measured data points
  geom_errorbar(data = HBEF_data_14C_sum,
                aes(y = Delta14C_mean, ymin = Delta14C_mean - Delta14C_sd,
                    ymax = Delta14C_mean + Delta14C_sd,
                    group = Horizon),
                width = 0.3) +
  geom_point(data = HBEF_data_14C_sum, aes(y = Delta14C_mean, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1955,2025), expand = c(0,0),
                     breaks = seq(1955,2025,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), limits = c(-175,1000),
                     expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled\nhorizon data", label = c("Oie", "Oa", "0-10 cm"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured\nhorizon data", label = c("Oie", "Oa/A", "0-10 cm"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3")) 

ggsave(file = paste0("./Output/HBEFall_SteadyStateModel_tpsModelFit_14C_Sensitivity_", lag_time, "_",
                     Sys.Date(), ".jpeg"), width = 10, height = 6) 

#### Calculate system C age and turnover time
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
    
    Sist_age <- systemAge(A = A, u = u)
    Trans_time <- transitTime(A = A, u = u)
    
    Iter_number[i] <- i
    System_age[i] <- as.numeric(Sist_age$meanSystemAge)
    oie_age[i] <- as.numeric(Sist_age$meanPoolAge[1])
    oa_age[i] <- as.numeric(Sist_age$meanPoolAge[2])
    min_age[i] <- as.numeric(Sist_age$meanPoolAge[3])
    Transit_time[i] <- as.numeric(Trans_time$meanTransitTime)
    
    
    age_results <- as.data.frame(cbind(
      Iter_number, System_age, oie_age, oa_age, min_age, Transit_time
    ))
    
    setTxtProgressBar(pb, i)
    
  }
  
  return(age_results)
}

# tpsMcmcFits$pars <- tpsMcmcFits$pars[-(1:1000), ]

# exclude first 1000 rows: MCMC algorithms are sensitive to their starting point
pars_df <- as.data.frame(tpsMcmcFits$pars[-(1:1000), ])

age_transit_dist <- propagation_fun(pars_df, 10000, C0)

write_csv(age_transit_dist, 
          file = paste0("./Output/ThreePoolSeriesModel_Age_Transit_Distribution_", lag_time, "_",
                        Sys.Date(), ".csv"))

age_transit_dist_sum <- age_transit_dist %>% 
  pivot_longer(!Iter_number, values_to = "age_yr", names_to = "Pools") %>% 
  group_by(Pools) %>% 
  summarise(mean_age = mean(age_yr),
            sd_age = sd(age_yr),
            median_age = median(age_yr),
            mad_age = mad(age_yr))
  
age_transit_dist %>% 
  pivot_longer(!Iter_number, values_to = "age_yr", names_to = "Pools") %>% 
  filter(Pools != "System_age", Pools != "Transit_time") %>%
  mutate(Pools = factor(Pools, levels = c("oie_age", "oa_age", "min_age"),
                        ordered = TRUE)) %>% 
  ggplot(aes(x = age_yr, color = Pools)) +
  geom_density(linewidth = 1) +
  facet_wrap(~Pools, scales = "free_y") +
  geom_vline(aes(xintercept = mean_age), 
             data = age_transit_dist_sum %>% 
               filter(Pools != "System_age",
                      Pools != "Transit_time"),
             linetype = "dashed") +
  theme_bw() +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_x_continuous("Age [yr]", limits = c(0,650), expand = c(0,0))
ggsave(file = paste0("./Output/HBEF_SteadyStateModel_tpsModelFit_14C_Age_Distribution_", 
                     lag_time, "_", Sys.Date(), ".jpeg"), width = 10, height = 6) 
  

 



