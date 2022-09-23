# this is an R-script to create a figure of a phylogenetic tree with the heatmap
# of the following traits: root diameter, specific root length, 
# root tissue density, and %N

# CLEAN UP
rm(list=ls(all=TRUE))

# load packages
library(phytools)
library(picante)
library(RColorBrewer)

# load in phylogenetic tree
tree_path = "C:\\Users\\Jorge\\Dropbox\\PC (3)\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre"
phy <- read.tree(tree_path)

file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\phylogenetic_tree.png"
png(file=file_name,height=5,width=3,units="in",res=800)
plot(phy, cex = .66)
dev.off()

#get the trait data
data_path <- "C:\\Users\\Jorge\\Dropbox\\PC (3)\\Documents\\data\\root_traits_summarized.csv"
roots <- read.csv(data_path, header = TRUE, row.names = 1)
#shorten data to just include the variables we want
roots <- roots[c('diameter', 'SRL', 'RTD', 'N_roots')]
roots <- na.omit(roots)
#use match.phylo to put the phylogenetic tree and dataframe in the same order
combined <- match.phylo.data(phy, roots)
phy <- combined$phy
roots <- combined$data

roots <- log(roots + 0.0001)

myPal <- brewer.pal(9, "PuBuGn")

tree_name = "C:\\Users\\Jorge\\Dropbox\\PC (3)\\Documents\\figures\\phylo_heatmap_nenew.png"
png(file=tree_name,height=12,width=8,units="in",res=1600)
phylo.heatmap(phy, roots, standardize = TRUE, labels = FALSE,
              colors = myPal, 
              fsize = c(0.5, 1, 0.6), split=c(0.70,0.3), pts=FALSE)
dev.off()
# tips for beautifying this figure?
# what is the best way to edit it while preserving the quality?
# I want to add a few pictures along the branches and also the trait labels 
# myself so that I can add the lambda values and p-values, but I am also 
# saving it as a PDF to preserve the quality which makes it hard to edit in powerpoint
# but saving it as an image in the export function really degrades the quality


