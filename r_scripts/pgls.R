# This R script includes code to run a phylogenetic generalized least squares (PGLS) regression
# on the root trait data of the prairie plants

# CLEAN UP
rm(list=ls(all=TRUE))

# LOAD PACKAGES
# section to import any packages from library
if(TRUE){
  library(vegan)
  
  library(picante)
  library(phytools)
  library(ape)
  library(dplyr)
}

# read in phylogenetic tree of the prairie plants
# this is the path to the tree if at working computer
tree_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre"
#this is the file path to the phylo tree if at home computer
#tree_path = "C:\\Users\\Jorge\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre"
phy <- read.tree(tree_path)
#phy$tip.label
plot(phy, cex=0.35)

#get the trait data
# this is the file path to the root trait data if at working computer
data_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv"
# file path to the root trait data if at home computer
#data_path = "C:\\Users\\Jorge\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv"
roots <- read.csv(data_path, header=TRUE, row.names=1)
#roots$diameter

# use match.phylo to match up the data with the phylogeny and get rid of any species
# that are not present in both
combined <- match.phylo.data(phy, roots)
phy <-combined$phy
roots <- combined$data
roots$diameter
#assigning the ordered roots dataframe to the previous existing one changes all
# numbers to characters, so we must change them back to doubles before proceeding
# might as well log transform the data at the same time
roots$diameter <- log(as.double(roots$diameter) + 0.0001)
roots$SRL <- log(as.double(roots$SRL) + 0.0001)
roots$RTD <- log(as.double(roots$RTD) + 0.0001)
roots$branching_intensity <- log(as.double(roots$branching_intensity) + 0.0001)

# add a new column that has the species names, then move it back to the first place
roots$species <- row.names(roots)
row.names(roots) <- 1:nrow(roots)

roots <- roots %>% select(species, everything())
#View(roots)

# had not called caper before because it masked the 'select' function of the dplyr package
library(caper)

# create a comparative.data object to have the correct format for PGLS

# set the node labels to null to prevent the error: 
# Labels duplicated between tips and nodes in phylogeny

phy$node.label <- NULL
prairie <- comparative.data(data=roots, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# get our variables for our regression
diameter <- roots$diameter
srl <- roots$SRL
branching <- roots$branching_intensity

diam_SRL_lm <- lm(srl ~ diameter)
summary(diam_SRL_lm)
# Multiple R-squared:  0.5015,	Adjusted R-squared:  0.4966 
# p-val < 2.2e-16
# create our first PGLS model, comparing diameter to SRL
model.pgls <- pgls(srl ~ diameter, data = prairie, lambda = "ML")
summary(model.pgls) # got a p-value of 1.063e-15 
# Multiple R-squared: 0.4692,	Adjusted R-squared: 0.464 

plot(srl ~ diameter, data = roots)
abline(model.pgls)
abline(diam_SRL_lm, lty="dashed")
#plotting the likelihood surface of lambda
profile_lambda=pgls.profile(model.pgls, which="lambda")
plot(profile_lambda)

# linear model and pgls model comparing SRL to branching intensity for all species
SRL_BI_lm <- lm(srl ~ branching)
summary(SRL_BI_lm)
#Multiple R-squared:  0.0854,	Adjusted R-squared:  0.07644 
# p-val = 0.002611
mod1.pgls <- pgls(srl ~ branching, data = prairie, lambda="ML")
summary(mod1.pgls) # p-val = 0.01119
# Multiple R-squared: 0.06142,	Adjusted R-squared: 0.05222 
#plotting srl ~ branching and adding line from pgls regression
plot(srl ~ branching, data=roots)
abline(mod1.pgls)
abline(SRL_BI_lm, lty="dashed")

# next: make linear models and pgls models for each of the following:
# 1. For Forbs, compare:
#   a. SRL vs diameter
#     i. linear model
#     ii. pgls
#   b. SRL vs branching intensity
#     i. linear model
#     ii. pgls
# 2. For Graminoids, compare:
#   a. SRL vs diameter
#     i. linear model
#     ii. pgls
#   b. SRL vs branching intensity
#     i. linear model
#     ii. pgls
# 3. For N-fixers, compare:
#   a. SRL vs diameter
#     i. linear model
#     ii. pgls
#   b. SRL vs branching intensity
#     i. linear model
#     ii. pgls
# 4. For monocots, compare:
#   a. SRL vs diameter
#     i. linear model
#     ii. pgls
#   b. SRL vs branching intensity
#     i. linear model
#     ii. pgls
# 5. For eudicots (dicots), compare:
#   a. SRL vs diameter
#     i. linear model
#     ii. pgls
#   b. SRL vs branching intensity
#     i. linear model
#     ii. pgls

#1. Forbs
# make our data only include Forbs
roots_F.df <- subset(roots, func_group=='F')
# have to go through the process of matching our phylogenetic tree to it once again,
# and then moving around the species column and making everything an integer,
# but first lets try to just
# create a new comparative data object for this using our current phylogeny
prairie_F <- comparative.data(data=roots_F.df, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
# just got a warning that data was dropped in compiling comparative data object, 
# so that may have taken care of our mismatching phylogenetic tree
# it looks like the comparative.data object did the matching, as the number of tip
# labels match up with the amount of rows we have in roots_F.df (72)

# now we can create our variables once again
diameter_F <- roots_F.df$diameter
srl_F <- roots_F.df$SRL
branching_F <- roots_F.df$branching_intensity

#a.i
#create our linear model for SRL vs diameter
srl_diam_F_lm <- lm(srl_F ~ diameter_F)
summary(srl_diam_F_lm)
# Multiple R-squared:  0.5673,	Adjusted R-squared:  0.5611 
# p-value: 2.313e-14

#a.ii: creating PGLS model for SRL vs diameter for forbs
srl_diam_F_pg.pgls <- pgls(srl_F ~ diameter_F, data = prairie_F, lambda = "ML")
summary(srl_diam_F_pg.pgls)
#Multiple R-squared: 0.6063,	Adjusted R-squared: 0.6006 
# p-value: 8.227e-16 

# lets plot this out
plot(srl_F ~ diameter_F, data=roots_F.df)
abline(srl_diam_F_pg.pgls)
abline(srl_diam_F_lm, lty="dashed")

#b.i creating linear model for SRL vs branching intensity of Forbs
srl_branching_F_lm <- lm(srl_F ~ branching_F)
summary(srl_branching_F_lm)
#Multiple R-squared:  0.06179,	Adjusted R-squared:  0.04839 
# p-value: 0.03525

#b.ii creating pgls model for SRL vs branching intensity of Forbs
srl_branching_F_pg.pgls <- pgls(srl_F ~ branching_F, data = prairie_F, lambda = "ML")
summary(srl_branching_F_pg.pgls)
# Multiple R-squared: 0.1136,	Adjusted R-squared: 0.1009 
# p-value: 0.003789 

# plot this out:
plot(srl_F ~ branching_F, data = roots_F.df)
abline(srl_branching_F_pg.pgls)
abline(srl_branching_F_lm, lty="dashed")

#2. Graminoids:
# create a data set with just the graminoids
roots_G <- subset(roots, func_group=='G')
# create another comparative data object for the graminoids
prairie_G <- comparative.data(data=roots_G, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# create our variables again
diameter_G <- roots_G$diameter
srl_G <- roots_G$SRL
branch_G <- roots_G$branching_intensity

#2.a.i : linear model for SRL vs Diameter of G
srl_diam_G_lm <- lm(srl_G ~ diameter_G)
summary(srl_diam_G_lm)
# Multiple R-squared:  0.5534,	Adjusted R-squared:  0.5236 
# p-value: 0.0006178

#2.a.ii: creating PGLS model for SRL vs diameter of G
srl_diam_G_pg.pgls <- pgls(srl_G ~ diameter_G, data = prairie_G, lambda = "ML")
summary(srl_diam_G_pg.pgls)
# Multiple R-squared: 0.5534,	Adjusted R-squared: 0.5236 
# p-value: 0.0006178 

# plot these:
plot(srl_G ~ diameter_G, data = roots_G)
abline(srl_diam_G_pg.pgls)
abline(srl_diam_G_lm, lty="dashed") # the lines are pretty much on top of each other

#2.b.i: create linear model comparing SRL vs branching intensity for G
srl_branching_G_lm <- lm(srl_G ~ branch_G)
summary(srl_branching_G_lm)
#Multiple R-squared:  0.01208,	Adjusted R-squared:  -0.05378
# p-value: 0.6745

#2.b.ii: create pgls model comparing SRL vs branching intensity for G
srl_branching_G_pg <- pgls(srl_G ~ branch_G, data = prairie_G, lambda = "ML")
summary(srl_branching_G_pg)
# Multiple R-squared: 0.07215,	Adjusted R-squared: 0.01029 
# p-value: 0.2972 

# plot these
plot(srl_G ~ branch_G, data = roots_G)
abline(srl_branching_G_pg)
abline(srl_branching_G_lm, lty="dashed")

#3. same analyses but for N-fixers this time
# create our dataset that has only n-fixers:
roots_N <- subset(roots, func_group=='N')

# create new comparative data object:
prairie_N <- comparative.data(data=roots_N, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# create our variables again:
diameter_N <- roots_N$diameter
srl_N <- roots_N$SRL
branch_N <- roots_N$branching_intensity

#3.a.i create linear model comparing SRL vs Diameter for N-fixers
srl_diam_N_lm <- lm(srl_N ~ diameter_N)
summary(srl_diam_N_lm)
# Multiple R-squared:  0.1496,	Adjusted R-squared:  0.06451 
# p-value: 0.2143

#3.a.ii create pgls model comparing SRL vs diameter for N-fixers
srl_diam_N_pg <- pgls(srl_N ~ diameter_N, data = prairie_N, lambda = "ML")
summary(srl_diam_N_pg)
# Multiple R-squared: 0.1496,	Adjusted R-squared: 0.06451 
# p-value: 0.2143 

# plot these
plot(srl_N ~ diameter_N, data = roots_N)
abline(srl_diam_N_pg)
abline(srl_diam_N_lm, lty="dashed")

#3.b.i create linear model comparing SRL vs branching intensity of N-fixers
srl_branching_N_lm <- lm(srl_N ~ branch_N)
summary(srl_branching_N_lm)
# Multiple R-squared:  0.1861,	Adjusted R-squared:  0.1047 
# p-value: 0.1614

#3.b.ii create pgls model comparing SRL vs branching intensity of N-fixers
srl_branching_N_pg <- pgls(srl_N ~ branch_N, data = prairie_N, lambda = "ML")
summary(srl_branching_N_pg)
# Multiple R-squared: 0.1861,	Adjusted R-squared: 0.1047 
# p-value: 0.1614 

# plot these:
plot(srl_N ~ branch_N, data = roots_N)
abline(srl_branching_N_pg)
abline(srl_branching_N_lm, lyt="dashed")

#4. Monocots!
#create our dataset that only has monocots
roots_M <- subset(roots, cot=="M")
# create a new comparative data object
prairie_M <- comparative.data(data=roots_M, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# create our variables again:
diameter_M <- roots_M$diameter
srl_M <- roots_M$SRL
branching_M <- roots_M$branching_intensity

#4.a.i: create a linear model comparing SRL vs Diameter for monocots
srl_diam_M_lm <- lm(srl_M ~ diameter_M)
summary(srl_diam_M_lm)
# Multiple R-squared:  0.4231,	Adjusted R-squared:  0.391 
# p-value: 0.001902

#4.a.ii: create a pgls model comparing SRL vs Diameter for monocots
srl_diam_M_pg <- pgls(srl_M ~ diameter_M, data = prairie_M, lambda = "ML")
summary(srl_diam_M_pg)
# Multiple R-squared: 0.448,	Adjusted R-squared: 0.4173
# p-value: 0.001249 

# plot these:
plot(srl_M ~ diameter_M, data = roots_M)
abline(srl_diam_M_pg)
abline(srl_diam_M_lm, lty="dashed")

#4.b.i: create a linear model comparing SRL vs Branching intensity for monocots
srl_branching_M_lm <- lm(srl_M ~ branching_M)
summary(srl_branching_M_lm)
# Multiple R-squared:  0.0933,	Adjusted R-squared:  0.04292
# p-value: 0.1903

#4.b.ii: create a pgls model comparing SRL vs Branching intensity for monocots
srl_branching_M_pg <- pgls(srl_M ~ branching_M, data = prairie_M, lambda = "ML")
summary(srl_branching_M_pg)
# Multiple R-squared: 0.0933,	Adjusted R-squared: 0.04292 
# p-value: 0.1903 

# plot these:
plot(srl_M ~ branching_M, data = roots_M)
abline(srl_branching_M_pg)
abline(srl_branching_M_lm, lty="dashed")

#5: Dicots!

#create a dataset of just dicots:
roots_D <- subset(roots, cot=='D')
#create new comparative data object of dicots
prairie_D <- comparative.data(data=roots_D, phy=phy, names.col="species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# create our variables again:
srl_D <- roots_D$SRL
diameter_D <- roots_D$diameter
branching_D <- roots_D$branching_intensity

#5.a.i: create a linear model comparing SRL vs Diameter for dicots
srl_diam_D_lm <- lm(srl_D ~ diameter_D)
summary(srl_diam_D_lm)
#Multiple R-squared:  0.4789,	Adjusted R-squared:  0.4725 
# p-value: 3.107e-13

#5.a.ii: create a pgls model comparing SRL vs Diameter for dicots
srl_diam_D_pg <- pgls(srl_D ~ diameter_D, data = prairie_D, lambda = "ML")
summary(srl_diam_D_pg)
# Multiple R-squared: 0.5075,	Adjusted R-squared: 0.5015 
# p-value: 2.977e-14 

# plot these:
plot(srl_D ~ diameter_D, data = roots_D)
abline(srl_diam_D_pg)
abline(srl_diam_D_lm, lty="dashed")

#5.b.i: create a linear model comparing SRL vs branching intensity of dicots
srl_branching_D_lm <- lm(srl_D ~ branching_D)
summary(srl_branching_D_lm)
# Multiple R-squared:  0.006521,	Adjusted R-squared:  -0.005595
# p-value: 0.4653

#5.b.ii: create a pgls model comparing SRL vs branching intensity of dicots 
srl_branching_D_pg <- pgls(srl_D ~ branching_D, data = prairie_D, lambda = "ML")
summary(srl_branching_D_pg)
# Multiple R-squared: 0.05947,	Adjusted R-squared: 0.048 
# p-value: 0.02539 

# plot these:
plot(srl_D ~ branching_D, data = roots_D)
abline(srl_branching_D_pg)
abline(srl_branching_D_lm, lty="dashed")

#TODO: put all of these R-squared values and p-values into a table