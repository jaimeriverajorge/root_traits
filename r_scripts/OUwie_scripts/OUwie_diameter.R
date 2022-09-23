#This script includes code to run different OU models using the OUwie package for
# the average diameter of 101 prairie species


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
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\OUwie_data\\OUwie_data_diameter.csv")

# read in phylogenetic tree:
phy <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\prairieTreeDichotomous.tre")
# this tree already has the selective regimes coded in (F, N, G)
# and it has been ran through the method multi2di to clean up any polytomies
roots$diameter <- log(roots$diameter + 0.001)

#run the three different models

#brown_OUwie <- OUwie(phy = phy, data = diameter_data, model='BM1')
#OU1 <- OUwie(phy = phy, data = diameter_data, model='OU1')
OU2 <- OUwie(phy = phy, data = roots, model = 'OUM')
#brown_OUwie # AICc = 61.51637
#OU1 # -162.6002
OU2 # 7.704128

#my_aic <- c(brown_OUwie$AICc, OU1$AICc, OU2$AICc)
#names(my_aic)<-c("fitBM1", "fitOU1", "fitOUM")
#aic.w(my_aic)

# Selective Optima:
# F: -1.03508596 
# G: -1.52696514 
# N: -1.00129532

boxplot(roots$diameter ~ roots$func_group, ylab = "Diameter (mm)", 
        xlab="Functional Group", main="Diameter vs Functional Group")
text(c(-1.03508596, -1.52696514, -1.00129532), "*", cex = 2)


#calculate selective optima (reverse log)
f_optima = exp(-1.03508596)
g_optima = exp(-1.52696514)
n_optima = exp(-1.00129532)
boxplot(roots$diameter ~ roots)