# the following R script conducts a PERMANOVA analysis on the root trait data for
# diameter, specific root length, and root tissue density
library(tidyverse)
library(vegan)

# CLEAN UP
rm(list=ls(all=TRUE))

roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv")
# we are going to analyze by functional group, so we want to remove the shrubs
roots_short.df <- roots.df[c("species","func_group", "diameter", "SRL", "RTD")]
roots_short.df <- subset(roots_short.df, func_group!='SS' & func_group!='SSN')


# do a db-RDA with distance matrix as D, SRL, RTD and then constrain
# the ordination by functional group
# would look like fxn(dist_matrix ~ func_group)

# constrained vs unconstrained? one functional group at a time or the column
# of functional group? should still get rid of shrubs? Yes

# for distance matrices, we only care about the trait data
my_dist <- vegdist(roots_short.df[3:5], method="euclidean")

# run adonis
my_test <- adonis2(formula = my_dist ~ func_group, data=roots_short.df, permutations = 10000)
my_test # P = 2e-4 
# what does this tell us? which part do we care about?
