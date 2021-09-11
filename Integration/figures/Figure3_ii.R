library(readr) #read data
library(dplyr) #organise df
library(tidyr) #clean df
library(ggplot2) #Visualization
library(ggthemes)
library(ggpubr)
require(sjstats)
library(forcats)
library(RColorBrewer)
setwd("/Volumes/NATHAN/SDEP/code_and_output/") #Path to the folder with .csv
wide <- read_csv('FCR_behav_SD.csv') #df wide format


ggplot(wide, aes(x=Accuracy_WR_SD, y=FCR_WR_SD)) + 
  geom_point(size=3)+
  geom_smooth(method=lm)+
  ylim(-0.05,0.2)+xlim(-55,40)+
  theme_classic()

ggplot(wide, aes(x=RT_WR_SD, y=FCR_WR_SD)) + 
  geom_point(size=3)+
  geom_smooth(method=lm)+
  ylim(-0.05,0.2)+xlim(-400,620)+
  theme_classic()

ggplot(wide, aes(x=Accuracy_SD_PRN, y=FCR_SD_PRN)) + 
  geom_point(size=3)+
  geom_smooth(method=lm)+
  theme_classic()

ggplot(wide, aes(x=RT_SD_PRN, y=FCR_SD_PRN)) + 
  geom_point(size=3)+
  geom_smooth(method=lm)+
  theme_classic()

