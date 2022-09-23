# Statistics using the phylogenetic tree of the prarie plants
# taken from this tutorial: http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
# and from Liam Revell's 2013 paper: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12066

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
Ntip(phy) # 127
phy$tip.label[1:5]

#get the trait data
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv", header=TRUE, row.names=1)
head(roots)
#pairs(roots)
#use match.phylo.data that takes a data object and phylo object and reports
# any species that are not present in the both datasets, and
# outputs a version of each object in the same order containing
# the same species
combined <- match.phylo.data(phy, roots)
#View(combined)
# replace original data with the matching data
phy <- combined$phy

plot(phy, cex=0.5)

#use phytools to map trait values on tree

roots <- combined$data
fixed_roots <- roots
# combined$data changes our variables into characters, below is a for-loop that
# goes through the columns of our dataframe and changes them back into doubles
for(i in 3:ncol(fixed_roots)){
  fixed_roots[,i] <- as.double(fixed_roots[,i])
}
# transform diameter and then run phylosig and try to get p-val
logD <- log(fixed_roots$diameter + 0.0001) 
names(logD) <- phy$tip.label
diam_Lam <- phylosig(phy, logD, method = "lambda")
diam_Lam$
View(fixed_roots)
# mapping the trait values of average diameter
diameter <- roots[, c('diameter')]
length(diameter)
names(diameter) = phy$tip.label
diameter_map <- contMap(phy, diameter, fsize = 0.5)
title(main="Diameter", cex.main = 1)

#change diameter into doubles
diameter <- as.double(diameter)
logD <- log(diameter + 0.01)
# calculating the phylogenetic signal in diameter
diam_l <- phylosig(phy, logD, method = "lambda") #lambda : 0.792489
print(diam_l)

# code to plot to png
png(file="C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\diam_map.png",
    height=2.5, width=3.5, units="in", res=600)
plot(diameter_map)
dev.off()

# mapping the trait values of SRL
srl <- roots[, c('SRL')]
srl <- as.double(srl)
srl
length(srl)
names(srl) = phy$tip.label
srl_map <- contMap(phy, srl, fsize = 0.5)
logSRL <- log(srl + 0.001)
phylosig(phy, logSRL, method = "lambda") # 0.498835

# mapping the trait values of RTD
rtd <- roots[, c('RTD')]
length(rtd)
names(rtd) = phy$tip.label
rtd_map <- contMap(phy, rtd, fsize = 0.5)

phylosig(phy, rtd, method = "lambda") # lambda : 0.803513 

# mapping the trait values of SLA
SLA <- roots[, c('SLA')]
length(SLA)
names(SLA) = phy$tip.label
SLA_map <- contMap(phy, SLA, fsize = 0.5)

phylosig(phy, SLA, method = "lambda") # lambda : 0.584393 

#mapping the trait values of RDMC
RDMC <- roots[, c('RDMC')]
length(RDMC)
names(RDMC) = phy$tip.label
RDMC_map <- contMap(phy, RDMC, fsize = 0.5) 

phylosig(phy, RDMC, method = "lambda") # lambda : 0.71332 

# mapping the trait values of "total_length"
total_length <- roots[, c('total_length')]
length(total_length)
names(total_length) = phy$tip.label
rtd_map <- contMap(phy, total_length, fsize = 0.5)

phylosig(phy, total_length, method = "lambda") # lambda : 0.615003 

# mapping the trait values of "branching_intensity"
branching_intensity <- roots[, c('branching_intensity')]
length(branching_intensity)
names(branching_intensity) = phy$tip.label
branching_map <- contMap(phy, branching_intensity, fsize = 0.5)

phylosig(phy, branching_intensity, method = "lambda") # lambda : 0.941081 


# mapping the trait values of "total_volume"
total_volume <- roots[, c('total_volume')]
length(total_volume)
names(total_volume) = phy$tip.label
volume_map <- contMap(phy, total_volume, fsize = 0.5)

phylosig(phy, total_volume, method = "lambda") # lambda : 6.61092e-05 


# creating a matrix of SRL and total average diameter to make
# a phylomorphospace
srl_diameter <- roots[c('SRL','total_avg_diameter')]
head(srl_diameter)
phylomorphospace(phy, srl_diameter,xlab="diameter",ylab="SRL", fsize=0.75, axes=1)

# will it work for more than two traits?
# it can only map two traits, unless you use phylomorphospace3d
srl_diam_rtd <- roots[c('SRL', 'total_avg_diameter', 'RTD')]
phylomorphospace3d(phy, srl_diam_rtd, angle=-30, label='off')

# creating an ordination of SRL and branching intensity
srl_branching <- roots[c('SRL', 'branching_intensity')]
phylomorphospace(phy, srl_branching)


# creating a phyloEM plot from the average diameter and SRL
View(t(srl_diameter))
mySim <- PhyloEM(phy, t(srl_diameter))

plot(mySim)
#summary <- roots %>% group_by(species) %>% summarise_each(funs(mean))
#View(summary)
#write.csv(summary, "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\root_traits_summarized.csv")
