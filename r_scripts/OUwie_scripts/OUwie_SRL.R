# R code for running an OUwie model with phylogenetic data from the prairie

#remove any objects in the environment
rm(list=ls(all=TRUE))

#load package OUwie
library(OUwie)

#load package phytools
library(phytools)

#load dplyr
library(dplyr)

library(picante)

#take in data file
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\OUwie_data\\root_trait_data_OUwie_SRL.csv")
roots$SRL <- log(roots$SRL)
View(roots)

# trying with the dichotomous tree
phy1 <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\prairieTreeDichotomous.tre")

# trying the BM model
brown_OUwie <- OUwie(phy = phy1, data = roots, model = 'BM1')

# the argument scaleHeight=1 was also giving the error, although this was not the source of it before
# since I had tried running it without that argument and still could not get it to work 
OU1 <- OUwie(phy = phy1, data = roots, model='OU1')
OU2 <- OUwie(phy = phy1, data = roots, model='OUM')
brown_OUwie # AIC = 670.5203
OU1 # AIC = 322.6895
OU2 # AIC = 312.7813 

#compare the fit of the models using weights
my_aic <- c(brown_OUwie$AICc, OU1$AICc, OU2$AICc)
names(my_aic)<-c("fitBM1", "fitOU1", "fitOUM")
aic.w(my_aic)

# Now we will plot the selective optima of each functional group
# for F: 8.1973958 
# for G: 9.0249427 
# for N: 7.3312831
boxplot(roots$SRL ~ roots$func_group, ylab = "Specific Root Length", 
        xlab="Functional Group", main="SRL vs Functional Group")
text(c(8.1973958, 9.0249427, 7.3312831), "*", cex=2)




# normal to have AIC values around 300