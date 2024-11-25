## Hubbard Brook archived soil samples project ##
## Radiocarbon analysis ##
## Sophie von Fromm ##
## 2024-08-01 ##

library(tidyverse)
library(SoilR)
library(forecast)

# Load HBEF data
HBEF_data <- read_csv("./Data/HBEF_data_all_2024-11-07.csv") %>% 
  # not needed once field notes are entered
  dplyr::select(-c(Plot_1_Horizons:Notes))

head(HBEF_data)

HBEF_data$Horizon <- factor(HBEF_data$Horizon,
                            levels = c("Oi/Oe", "Oa/A", "Mineral_0_10"),
                            ordered = TRUE)

# Quick look at the radiocarbon data
HBEF_data %>% 
  drop_na(Delta14C) %>% 
  ggplot(aes(x = Year, y = Delta14C, color = Plot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 14) +
  facet_wrap(~Horizon)

### Set-up model at steady state based on current C budget from HBEF
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
  # filter(Year >= 1949) %>% 
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() 

# time interval for model
Year <- seq(-53050, 2025, by = 1)

# initial C stocks in each pool
C0 <- c(115, 1394, 595, 1484, 1480, 12770)

# construct decomposition matrix
ks <- c(1/1.5, 1/6, 1/4, 1/14, 1/6, 1/81)

A <- -abs(diag(ks))
A[2,1] <- ks[1] * 41/(41+35)
A[4,2] <- ks[2] * 66/(66+175+10)
A[4,3] <- ks[3] * 37/(37+100)
A[6,4] <- ks[4] * 28/(28+75)
A[6,5] <- ks[5] * 119/(119+125)
A[6,2] <- ks[2] * 10/(10+66+175)

# C inputs
LI <- 210
RI1 <- 76
RI3 <- 137
RI5 <- 244

#In <- matrix(nrow = 6, ncol = 1, c(RI1, LI, RI3, 0, RI5, 0))
In <- c(RI1, LI, RI3, 0, RI5, 0)

# initial F values; mineral soil: avg Delta14C value
F0 <- ConstFc(values = c(0, 0, 0, 0, 0, 0), format = "Delta14C")

HBEF_model <- Model_14(
  t = Year,
  A = A,
  ivList = C0,
  initialValF = F0,
  inputFluxes = In,
  inputFc = BoundFc(map = NHZone2_2023,  format = "Delta14C")
)

## Get 14C data
C14pools <- getF14(HBEF_model) # 14C of each pools

# Convert into dataframe
C14pools_df <- C14pools %>% 
  data.frame(Year) %>% 
  rename_with(~gsub("X", "pool_", .x, fixed = TRUE)) %>% 
  pivot_longer(!Year, names_to = "pool", values_to = "Delta14C")

# Plot: only the three main pools (Oie: Pool 2, Oa/A: Pool 4; Mineral: Pool 6)
HBEF_data_14C_sum <- HBEF_data %>% 
  drop_na(Delta14C) %>% 
  #remove roots
  filter(Plot != "all fine roots") %>% 
  group_by(Year, Horizon) %>% 
  summarise(Delta14C_mean = mean(Delta14C),
            Delta14C_sd = sd(Delta14C))

NHZone2_2023 %>%  
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() +
  #filter for pool 2, 4, and 6
  geom_line(data = C14pools_df %>%
              filter(pool == "pool_2"| pool == "pool_4"| pool == "pool_6"),
            aes(color = pool), linewidth = 1) +
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
  scale_color_manual("Modeled\nhorizon data", label = c("Oie", "Oa", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured\nhorizon data", label = c("Oie", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3"))

ggsave(file = paste0("./Output/HBEF_SteadyStateModel_CBudget_14C_", Sys.Date(),
                     ".jpeg"), width = 10, height = 6)

# Plot: only the root pools (Oie: Pool 1, Oa/A: Pool 3; Mineral: Pool 5)
HBEF_root_14C_sum <- HBEF_data %>% 
  drop_na(Delta14C) %>% 
  #only roots
  filter(Plot == "all fine roots") 

NHZone2_2023 %>%  
  ggplot(aes(x = Year, y = Delta14C)) +
  geom_line() +
  #filter for pool 2, 4, and 6
  geom_line(data = C14pools_df %>%
              filter(pool == "pool_1"| pool == "pool_3"| pool == "pool_5"),
            aes(color = pool), linewidth = 1) +
  # Add measured data points
  geom_point(data = HBEF_root_14C_sum, aes(y = Delta14C, fill = Horizon),
             shape = 21, size = 2) +
  scale_x_continuous("Year", limits = c(1995,2025), expand = c(0,0),
                     breaks = seq(1995,2025,10)) +
  scale_y_continuous(expression(paste(Delta^14, "C [‰]")), limits = c(0,500),
                     expand = c(0,0)) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual("Modeled\nhorizon data", label = c("Oie", "Oa", "Mineral"),
                     values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  scale_fill_manual("Measured\nhorizon data", label = c("Oie", "Oa/A", "Mineral"),
                    values = c("#33a02c", "#b2df8a", "#a6cee3"))

## System age and transit time
ages <- seq(0,200)
SA <- systemAge(A = A, u = In, a = ages)
TT <- transitTime(A = A, u = In, a = ages)

SA$meanSystemAge
SA$meanPoolAge

TT$meanTransitTime

# Plot density distribution
par(mfrow=c(3,2))
pools <- c("Roots1", "Oi", "Root2", "Oa/A", "Roots 3", "Mineral")
for(i in 1:6){
  plot(ages, SA$poolAgeDensity[,i], type = "l", main = pools[i], 
       ylab = "Probability density", xlab = "Age", bty = "n")
}

#maybe add mean pool age and transit time for each pool to the plot

R14t <- getF14R(HBEF_model) # Average 14C of total respiration

C14t <- getF14C(HBEF_model) # Average 14C of all pools in bulk


