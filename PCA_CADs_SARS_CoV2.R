
install_github('sinhrks/ggfortify')
install.packages("readr")
# You need to load a package (like magrittr or dplyr) that defines the function first, then it should work.
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(tidyverse)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(readr)

library(devtools)
library(ggfortify)
library(ggplot2)

CADs_SARS_screens <- read_csv("OneDrive-v3/Data/CADs/v2/CADs_SARS-CoV-2_CRISPRscreens.csv")
df <- CADs_SARS_screens 
df['sars_5'] <- df['lfc_MOI1']
df['sars_4'] <- df['Cas9-v1 Avg.']
df0 <- df[ , c("verap", "amiod", "sert", "chloro", "nortrip", "fluox", 'sars_5', 'sars_4')]
df0_na <- na.omit(df0)

maxs <- apply(df0_na, 2, max)    
mins <- apply(df0_na, 2, min)
#scale(data, center = mins, scale = maxs - mins)

dat.scale <- scale(df0_na, center = mins, scale = maxs - mins)

df0_na_scaled <- as.data.frame(dat.scale)

tt2_t <- t(df0_na_scaled)

autoplot(prcomp(tt2_t), data = tt2_t, label = TRUE)

CADs_1_t_na <- na.omit(CADs_1_t)
tt2_t_na <- na.omit(tt2_t)


