#Load packages
library(tidyverse)
# install.packages("writexl")
library(writexl)

#Create primary df
HBBMaster <- read.csv("./Data/HBBMaster.csv")
HBBclean <- HBBMaster %>% filter(!is.na(C))

HBBHistoric <- read.csv("./Data/Historic.csv") %>% 
  filter(Year>=1997) %>% 
  filter(PerCentN > 0) %>%
  filter(Plot>=159)


#Divide Spreadsheet by horizon
HBB_A <- HBBclean %>% filter(Horizon == "Oa/A")
HBB_i <- HBBclean %>% filter(Horizon == "Oi/Oe")
HBB_m <- HBBclean %>% filter(Horizon == "Mineral_0_10")

Historic_A <- HBBHistoric %>% filter(Horizon == "Oa")
Historic_i <- HBBHistoric %>% filter(Horizon == "Oie")
Historic_m <- HBBHistoric %>% filter(Horizon == "min")

#Graph C
Cmain <- ggplot() +
  geom_point(data=HBBclean, aes(x=Year, y=C, color=Horizon)) +
  labs(y="Percentage Carbon")
Cmain

#Take mean for each year across the sample sites - not currently used in a later step
HBByearavgs_A <- HBB_A %>%
  group_by(Year) %>%
  summarize(avgC = mean(C), avgN = mean(N))

HBByearavgs_i <- HBB_i %>%
  group_by(Year) %>%
  summarize(avgC = mean(C), avgN = mean(N))

HBByearavgs_m <- HBB_m %>%
  group_by(Year) %>%
  summarize(avgC = mean(C), avgN = mean(N))

Historicyearavgs_A <- Historic_A %>%
  group_by(Year) %>%
  summarize(avgC = mean(PerCentC), avgN = mean(PerCentN))

Historicyearavgs_i <- Historic_i %>%
  group_by(Year) %>%
  summarize(avgC = mean(PerCentC), avgN = mean(PerCentN))

Historicyearavgs_m <- Historic_m %>%
  group_by(Year) %>%
  summarize(avgC = mean(PerCentC), avgN = mean(PerCentN))

#Graph average values - These are bad, ignore
avgC_A <- ggplot() +
  geom_smooth(data=HBByearavgs_A, aes(x=Year, y=avgC), color="lightcoral", fill="pink", method="lm")+
  geom_point(data=HBByearavgs_A, aes(x=Year, y=avgC), color="tomato3") +
  labs(title="Mean Percentage Carbon in Oa/A Horizon
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")

avgN_A <- ggplot() +
  geom_smooth(data=HBByearavgs_A, aes(x=Year, y=avgN), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBByearavgs_A, aes(x=Year, y=avgN), color="skyblue4")+
  labs(title="Mean Percentage Nitrogen in Oa/A Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")

avgC_i <- ggplot() +
  geom_smooth(data=HBByearavgs_i, aes(x=Year, y=avgC), color="lightcoral", fill="pink", method="lm")+
  geom_point(data=HBByearavgs_i, aes(x=Year, y=avgC), color="tomato3")+
  labs(title="Mean Percentage Carbon in Oi/Oe Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")

avgN_i <- ggplot() +
  geom_smooth(data=HBByearavgs_i, aes(x=Year, y=avgN), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBByearavgs_i, aes(x=Year, y=avgN), color="skyblue4")+
  labs(title="Mean Percentage Nitrogen in Oi/Oe Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")

avgC_m <- ggplot() +
  geom_smooth(data=HBByearavgs_m, aes(x=Year, y=avgC), color="lightcoral", fill="pink", method="lm")+
  geom_point(data=HBByearavgs_m, aes(x=Year, y=avgC), color="tomato3")+
  labs(title="Mean Percentage Carbon in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")

avgN_m <- ggplot() +
  geom_smooth(data=HBByearavgs_m, aes(x=Year, y=avgN), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBByearavgs_m, aes(x=Year, y=avgN), color="skyblue4") +
  labs(title="Mean Percentage Nitrogen in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")

#Same thing but with the og data set to reduce the ribbon size
C_A <- ggplot() +
  geom_smooth(data=HBB_A, aes(x=Year, y=C), color="lightcoral", fill="pink", method="lm")+
  geom_point(data=HBB_A, aes(x=Year, y=C), color="tomato3")+
  labs(title="Mean Percentage Carbon in Oa/A Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")

N_A <- ggplot() +
  geom_smooth(data=HBB_A, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBB_A, aes(x=Year, y=N), color="skyblue4")+
  labs(title="Mean Percentage Nitrogen in Oa/A Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")

C_i <- ggplot() +
  geom_smooth(data=HBB_i, aes(x=Year, y=C), color="darkorchid4", fill="orchid2", method="lm")+
  geom_point(data=HBB_i, aes(x=Year, y=C), color="darkorchid1")+
  labs(title="Mean Percentage Carbon in Oi/Oe Horizon 
       at Hubbard Brook West of Watershed Six Low Elevation", y="Percentage C")

N_i <- ggplot() +
  geom_smooth(data=HBB_i, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBB_i, aes(x=Year, y=N), color="skyblue3")+
  labs(title="Mean Percentage Nitrogen in Oi/Oe Horizon 
       at Hubbard Brook West of Watershed Six Low Elevation", y="Percentage N")

C_m <- ggplot() +
  geom_smooth(data=HBB_m, aes(x=Year, y=C), color="lightcoral", fill="pink", method="lm")+
  geom_point(data=HBB_m, aes(x=Year, y=C), color="tomato3")+
  labs(title="Mean Percentage Carbon in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")

N_m <- ggplot() +
  geom_smooth(data=HBB_m, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue", method="lm")+
  geom_point(data=HBB_m, aes(x=Year, y=N), color="skyblue4")+
  labs(title="Mean Percentage Nitrogen in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")

#graph comparisons with historic data set
compareC_A <- ggplot() +
  geom_smooth(data=HBB_A, aes(x=Year, y=C), color="lightcoral", fill="pink")+
  geom_point(data=HBB_A, aes(x=Year, y=C), color="tomato3")+
  geom_smooth(data=Historic_A, aes(x=Year, y=PerCentC), color="mediumorchid2") +
  geom_point(data=Historic_A, aes(x=Year, y=PerCentC), color="darkorchid4") +
  labs(title="Percentage Carbon in Oa/A Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")
compareN_A <- ggplot() +
  geom_smooth(data=HBB_A, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue")+
  geom_point(data=HBB_A, aes(x=Year, y=N), color="skyblue4")+
  geom_smooth(data=Historic_A, aes(x=Year, y=PerCentN), color="mediumorchid2") +
  geom_point(data=Historic_A, aes(x=Year, y=PerCentN), color="darkorchid4") +
  labs(title="Percentage Nitrogen in Oa/A Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")
compareC_i <- ggplot() +
  geom_smooth(data=HBB_i, aes(x=Year, y=C), color="lightcoral", fill="pink")+
  geom_point(data=HBB_i, aes(x=Year, y=C), color="tomato3")+
  geom_smooth(data=Historic_i, aes(x=Year, y=PerCentC), color="mediumorchid2") +
  geom_point(data=Historic_i, aes(x=Year, y=PerCentC), color="darkorchid4") +
  labs(title="Percentage Carbon in Oi/Oe Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")
compareN_i <- ggplot() +
  geom_smooth(data=HBB_i, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue")+
  geom_point(data=HBB_i, aes(x=Year, y=N), color="skyblue4")+
  geom_smooth(data=Historic_i, aes(x=Year, y=PerCentN), color="mediumorchid2") +
  geom_point(data=Historic_i, aes(x=Year, y=PerCentN), color="darkorchid4") +
  labs(title="Percentage Nitrogen in Oi/Oe Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")
compareC_m <- ggplot() +
  geom_smooth(data=HBB_m, aes(x=Year, y=C), color="lightcoral", fill="pink")+
  geom_point(data=HBB_m, aes(x=Year, y=C), color="tomato3")+
  geom_smooth(data=Historic_m, aes(x=Year, y=PerCentC), color="mediumorchid2") +
  geom_point(data=Historic_m, aes(x=Year, y=PerCentC), color="darkorchid4") +
  labs(title="Percentage Carbon in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage C")
compareN_m <- ggplot() +
  geom_smooth(data=HBB_m, aes(x=Year, y=N), color="midnightblue", fill="lightskyblue")+
  geom_point(data=HBB_m, aes(x=Year, y=N), color="skyblue4")+
  geom_smooth(data=Historic_m, aes(x=Year, y=PerCentN), color="mediumorchid2") +
  geom_point(data=Historic_m, aes(x=Year, y=PerCentN), color="darkorchid4") +
  labs(title="Percentage Nitrogen in Mineral Horizon 
       at Hubbard Brook Watershed Six Low Elevation", y="Percentage N")
  
#Write data frame to excel file
write_xlsx(HBB_i, "C:HBB_i.xlsx")
