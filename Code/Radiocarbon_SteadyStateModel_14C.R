## Hubbard Brook archived soil samples project ##
## Radiocarbon analysis - contraint model with 14C ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(SoilR)
library(forecast)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-08-01.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes))

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

HBEF_data %>% 
  group_by(Horizon) %>% 
  reframe(mean_14C = mean(Delta14C, na.rm = TRUE))

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
lambda <- 0.00012097
FpfromF <- function(X){ # X must be a data.frame with first column time and second column Fraction Modern
  y = X[,1]
  fM = X[,2]
  Fp = fM*exp((1950-y)*lambda)
  return(data.frame(year = y, Fp = Fp))
}

FpNZ2 <- FpfromF(Hua2021$NHZone2[,c(1,4)]) #Absolute fraction modern F'
qNZ2<- spline(FpNZ2, xout = y) #Spline interpolation of the NH_Zone 2 data set at a quarterly basis
NZ2 <- ts(qNZ2$y, start = y0, freq = q) #Transformation into a time-series object

mNZ2 <- ets(NZ2) #Fits an exponential smoothing state space model to the time series

foryrs <- 10 # Number of years to forecast
fNZ2 <- forecast(mNZ2, h = foryrs*q, level = c(69,90)) #Uses the fitted model to forecast 10 years into the future

NHZone2_2023 <- data.frame(Year = c(NHZone2[-dim(NHZone2)[1],1],
                                    seq(tsp(fNZ2$mean)[1], tsp(fNZ2$mean)[2], by = 1/tsp(fNZ2$mean)[3])),
                           Delta14C = c(NHZone2[-dim(NHZone2)[1],2], as.numeric((fNZ2$mean)-1)*1000))

NHZone2_2023 %>% 
  filter(Year >= 1949) %>% 
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() 

# time interval for model
Year <- seq(-53050, 2025, by = 1)

# initial C stocks in each pool
C0 <- c(1394, 1484, 3170)

# initial Delta14C in each pool: 0 for all fast cycling pools; -26 for mineral (= average value over all years)
init14C <- c(0, 0, -25)

# lag-time before C enters soils: ask Josh again
lag_time <- 8

# Number of model iterations: set to 15000 later
itr <- 500

# C inputs
In <- 210

###NOT DONE YET####

## Set-up 14C pool (three pools in series)
ThreePSeriesModel <- function(pars){
  mod=ThreepSeriesModel14(
    t = years,
    ks = pars[1:3],
    C0 = as.numeric(C0), 
    F0_Delta14C = iniD14C,
    In = In,
    a21= pars[4]*pars[1], #0.5*pars[1], #change if do not want to estimate alpha
    a32= pars[5]*pars[2], #0.1*pars[2], #change if do not want to estimate alpha
    inputFc = NHZone2_2023,
    lag = lag_time
  )
  model_result = getF14(mod)
  return(data.frame(time = years, Oie = model_result[,1],
                    Oa = model_result[,2], min = model_result[,3]))
}

min_data <- HBEF_data %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Mineral_0_10") %>% 
  dplyr::select(Year, Delta14C)

oa_data <- HBEF_data %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Oa/A") %>% 
  dplyr::select(Year, Delta14C)

oie_data <- HBEF_data %>% 
  drop_na(Delta14C) %>% 
  filter(Horizon == "Oi/Oe") %>% 
  dplyr::select(Year, Delta14C)

costF1 <- function(pars){
  funccall = ThreePSeriesModel(pars)
  cost1 = modCost(model = funccall, obs = flf,  err = "sd")
  cost2 = modCost(model = funccall, obs = oa_data,  err = "sd", cost = cost1)
  cost3 = modCost(model = funccall, obs = min_data,   err = "sd", cost = cost2)
  return(cost3)
}



