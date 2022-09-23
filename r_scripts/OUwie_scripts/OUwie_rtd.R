# This script includes code for building different OU models using the R package
# OUwie, to determine if root tissue density (RTD) varies among functional groups:
# Nitrogen fixers (N), Forbs (F), Graminoids (G)

#remove any objects in the environment
rm(list=ls(all=TRUE))

#load package OUwie
library(OUwie)

#load package phytools
library(phytools)

#load dplyr
library(dplyr)

library(picante)

#include root traits file
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\OUwie_data\\OUwie_data_rtd.csv")

# read in phylogenetic tree:
phy <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\prairieTreeDichotomous.tre")
# this tree already has the selective regimes coded in (F, N, G)

# build the models using OUwie
#brown_OUwie <- OUwie(phy = phy, data= rtd_data, model='BM1')
#OU1 <- OUwie(phy = phy, data= rtd_data, model='OU1')
OU2 <- OUwie(phy = phy, data= roots, model='OUM')
#brown_OUwie # AIC = 465.7124
#OU1 # AIC = 193.1146
OU2 # AIC = 172.1931

#my_aic <- c(brown_OUwie$AICc, OU1$AICc, OU2$AICc)
#names(my_aic)<-c("fitBM1", "fitOU1", "fitOUM")
#aic.w(my_aic)

#Selective optima:
# F: -1.19474918
# G: -1.0629458 
# N: -0.2062934
boxplot(roots$RTD ~ roots$func_group, ylab = "Root Tissue Density (RTD)", 
        xlab="Functional Group", main="RTD vs Functional Group")
text(c(-1.19474918, -1.0629458, -0.2062934), "*", cex=2)
