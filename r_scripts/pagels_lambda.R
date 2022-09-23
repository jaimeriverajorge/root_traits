# CLEAN UP
rm(list=ls(all=TRUE))

# LOAD PACKAGES
# section to import any packages from library
if(TRUE){
  library(vegan)
  library(picante)
  library(ape)
  library(dplyr)
  library(phytools)
  library(PhylogeneticEM)
}

phy <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre")
plot(phy, cex=0.35)
#Ntip(phy) # 127
#phy$tip.label[1:5]

#get the trait data
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized.csv", header=TRUE, row.names=1)
n_data <- roots[c('N_roots')]
n_data <- na.omit(n_data)

#match the n_data with the phy tree
combined_n <- match.phylo.data(phy, n_data)
phy <- combined_n$phy
n_data <- combined_n$data
n_vec <- log(n_data$`data[res$phy$tip.label, ]`)
names(n_vec) <- phy$tip.label
n_lam <- phylosig(phy, n_vec, method = 'lambda', test = TRUE)
n_lam
#use match.phylo.data that takes a data object and phylo object and reports
# any species that are not present in the both datasets, and
# outputs a version of each object in the same order containing
# the same species
combined <- match.phylo.data(phy, roots)
#View(combined)
# replace original data with the matching data
phy <- combined$phy

roots <- combined$data
fixed_roots <- roots

N_conc <- log(as.double(roots$N_roots) + 0.0001)
names(N_conc) <- phy$tip.label
n_lam <- phylosig(phy, N_conc, method = 'lambda', test = TRUE)
n_lam # lambda: 0.560594 #p-val < 0.0001
# combined$data changes our variables into characters, below is a for-loop that
# goes through the columns of our dataframe and changes them back into doubles
for(i in 3:ncol(fixed_roots)){
  fixed_roots[,i] <- as.double(fixed_roots[,i])
}

n_data <- roots[c('N_roots')]
n_data <- na.omit(n_data)
n_vec <- log(as.double(n_data$N_roots) + 0.0001)


# transform diameter and then run phylosig and try to get p-val
logD <- log(fixed_roots$diameter + 0.0001) 
names(logD) <- phy$tip.label
diam_Lam <- phylosig(phy, logD, method = "lambda", test = TRUE)
diam_Lam #lambda: 0.80691, p-val: 1.34601e-11

#transform SRL and run phylosig
logSRL <- log(fixed_roots$SRL + 0.0001)
names(logSRL) <- phy$tip.label
srl_Lam <- phylosig(phy, logSRL, method = "lambda", test = TRUE)
srl_Lam #lambda: 0.611645, p-val = 1.11822e-06 

#transform RTD and run phylosig

#make a loop instead
fixed_roots <- fixed_roots[,3:13]
cols = colnames(fixed_roots)
lambda_values <- matrix(, nrow = 11, ncol = 3)
lambda_values

cols

rowCount <- 1
colCount <- 1
# a for loop to calculate the lambda values in our data
for(i in cols){
  colCount <- 1
  if(i == 'SLA' || i == 'RDMC'){
    lambda_values[rowCount,colCount] <- i
    colCount <- colCount+1
    myTrait <- sqrt(fixed_roots[,i])
    names(myTrait) <- phy$tip.label
    lam <- phylosig(phy, myTrait, method="lambda", test=TRUE)
    lambda <- lam$lambda
    p <- lam$P
    lambda_values[rowCount,colCount] <- lambda
    colCount <- colCount + 1
    lambda_values[rowCount,colCount] <- p
  }
  else {
    lambda_values[rowCount,colCount] <- i
    colCount <- colCount+1
    myTrait <- log(fixed_roots[,i] + 0.001)
    names(myTrait) <- phy$tip.label
    lam <- phylosig(phy, myTrait, method="lambda", test=TRUE)
    lambda <- lam$lambda
    p <- lam$P
    lambda_values[rowCount,colCount] <- lambda
    colCount <- colCount + 1
    lambda_values[rowCount,colCount] <- p
  }
  rowCount <- rowCount + 1
}
lambda_values
lam
library(MASS)
write.matrix(lambda_values, file="lambda.csv")
#TODO: Clean up this mess of code
