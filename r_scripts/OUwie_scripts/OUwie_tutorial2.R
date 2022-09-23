# following this tutorial: http://www.phytools.org/Cordoba2017/ex/10/Multi-regime.html

library(phytools)
library(OUwie)

rm(list=ls(all=TRUE))

anoleTree <- read.tree("Anolis.tre")
anoleData <- read.csv("svl.csv")

crownGiantSpecies <- c("equestris", "luteogularis", "smallwoodi", 
                       "noblei", "baracoae", "baleatus", "barahonae", "ricordii", "cuvieri", "garmani")

isCrownGiant <- anoleData[,1] %in% crownGiantSpecies

ecomorph <- sapply(isCrownGiant, function(x) if(x) "cg" else "other")
ecomorph <- as.factor(ecomorph)
names(ecomorph) <- anoleData[,1]

ecomorphTree <- make.simmap(anoleTree, ecomorph, model="ER")
plot(ecomorphTree)

svl <- anoleData[,2]
names(svl)<-anoleData[,1]
svl

fitBrownie<-brownie.lite(ecomorphTree, svl)
fitBrownie
ecomorph

svl_df<-data.frame(anoleData[,1], ecomorph, svl)
fitOUM<-OUwie(ecomorphTree, svl_df,model="OUM",simmap.tree=TRUE)


roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\fixed_data.csv", row.names = 1)
phy <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\tr.analysis.tre")
combined <- match.phylo.data(phy, roots)
phy <- combined$phy
roots<- combined$data
roots$species <- row.names(roots)
# change the row names back to numbers:
row.names(roots) <- 1:nrow(roots)
# our dataframe currently has the species name as the last column, but we want it
# in the first column, willl use dplyr to make this change

fixed_data <- roots %>% select(species, everything())
View(fixed_data)
states <- c("F", "F", "F", "F", "F", "F", "F" ,"F" ,"F" ,"F",
            "F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
            "F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
            "F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
            "F", "F", "F", "F", "F", "F", "F", "F", "F", "F",
            "F", "F", "F", "F", "F", "F", "F", "N", "N", "N",
            "N", "N", "N", "N", "N", "N", "N", "N", "F", "F",
            "F", "F", "F", "F", "F", "F", "G", "G", "G", "G",
            "G", "G", "G", "G", "G", "G", "G", "G", "G", "G",
            "G", "G", "F")
phy$node.label <- states
srl <- roots[,3]
names(srl)<-roots[,1]

myBrownie <- brownie.lite(phy, srl)
