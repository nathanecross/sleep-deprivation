library(readr) #read data
library(dplyr) #organise df
library(tidyr) #clean df
library(ggplot2) #Visualization
library(ggthemes)
library(ggpubr)
require(sjstats)
library(forcats)
library(RColorBrewer)
setwd(".../integration/data/integration") # <- Set as path to the folder with .csv
wide <- read_csv('FCR_behav_SD.csv') 


### GSF
gsf <- read_csv('gsf_Int.csv') #df wide format
ggplot(gsf, aes(x=gsf_WR, y=I_WR)) + 
  geom_point(size=3,colour ='#0072bd' )+
  geom_smooth(method=lm,colour='darkgrey')+
  ylim(600,900)+xlim(1,3.5)+
  theme_classic()

ggplot(gsf, aes(x=gsf_SD, y=I_SD)) + 
  geom_point(size=3,colour='#ce1c31')+
  geom_smooth(method=lm,colour='darkgrey')+
  ylim(600,900)+xlim(1,5)+
  theme_classic()

ggplot(gsf, aes(x=gsf_PRN, y=I_PRN)) + 
  geom_point(size=3,colour='#edb020')+
  geom_smooth(method=lm,colour='darkgrey')+
  ylim(600,900)+xlim(1,5)+
  theme_classic()

