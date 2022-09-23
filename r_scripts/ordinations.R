# the following script conducts ordinations on the root trait data, 
# which includes using a db-RDA, and afterwards a regular PCA
# functional group as a constraining factor
#db-RDA: distance based redundancy analysis

library(vegan)

# CLEAN UP
rm(list=ls(all=TRUE))

roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv")
# we are going to analyze by functional group, so we want to remove the shrubs
roots_short.df <- roots.df[c("species","func_group", "diameter", "SRL", "RTD")]
roots_short.df <- subset(roots_short.df, func_group!='SS' & func_group!='SSN')

# for distance matrices, we only care about the trait data
my_dist <- vegdist(roots_short.df[3:5], method="euclidean")

#construct db-RDA
dbRDA = capscale(my_dist~func_group, data=roots_short.df)
dbRDA

#what does this all mean?
plot(dbRDA, display=c('bp','sites'), scaling='species')
#nmds <- metaMDS(my_dist, k = 2, try = 999)
plot(nmds)

# below we conduct a principal component analysis of our trait data (diameter, SRL, RTD)
pca <- princomp(roots.df[4:6])
loadings<-pca$loadings[,1:2]
loadings
biplot(pca)

#TODO: read about ANCOVA, try doing full dataset (nmds first ,if not insufficient data then db-RDA)
#trying with the full dataset (2 samples per species)
all_roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_temp_Copy.csv")
all_roots.df <- all_roots.df[c("Name", "Species", "func_group", "total_avg_diameter",
                               "SRL", "RTD", "RDMC", "branching_intensity")]
#take out shrubs
all_roots.df <- subset(all_roots.df, func_group!='SS' & func_group!='SSN')
#make distance matrix
all_dist <- vegdist(all_roots.df[4:8], method = "euclidean")
# run nmds
#nmds <- metaMDS(all_dist, k = 10, try = 999)
#plot(nmds)

pca1 <- princomp(all_roots.df[4:6])
biplot(pca1)
pcoa <- cmdscale(all_dist, k=10, eig=FALSE)
ordiplot(pcoa, display = "sites")

library(ape)

m_dist <- vegdist(all_roots.df[4:8], method = "mahalanobis")
pcoa_m <- cmdscale(m_dist, k=10, eig=FALSE)
ordiplot(pcoa_m, display = "sites")
pca2<- princomp(all_roots.df[4:8])
biplot(pca2)

#pcoa using ape
new_pcoa <- pcoa(m_dist)
biplot.pcoa(new_pcoa)
biplot(new_pcoa, all_roots.df[4:8])

#TODO: take the averaged dataset and log transform the variables then construct a PCOA and PCA
roots_short.df$diameter <- log(roots_short.df$diameter + 0.0001)
roots_short.df$SRL <- log(roots_short.df$SRL + 0.0001)
roots_short.df$RTD <- log(roots_short.df$RTD + 0.0001)

#create distance matrix
mah_dist <- vegdist(roots_short.df[3:5], method = "mahalanobis")
pcoa_log <- pcoa(mah_dist)
biplot(pcoa_log, roots_short.df[3:5])

# pca with the log transformed data and no species names
new_pca <- princomp(roots_short.df[3:5])
# taking the log transformation of the variables helped a lot, now on to plotting it with 
# a different color for each functional group
library(ggfortify)
autoplot(new_pca, data=roots_short.df, colour="func_group", loadings=TRUE)

#plot(pcoa_log, display=c('bp','wa'),scaling='species')

library(factoextra)
fviz_pca_biplot(new_pca, col.ind = roots_short.df$func_group)

roots_more_traits.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv")
roots_more_traits.df <- roots_more_traits.df[c("species","func_group", "diameter",
                                               "SRL", "RTD", "RDMC", "branching_intensity",
                                               "SLA","root_shoot_ratio")]
#remove row 3 because it has NA for SLA
test.df <- roots_more_traits.df[-c(3),]
#remove the shrubs
test.df <- subset(test.df, func_group!='SS' & func_group!='SSN')

#log transform variables first
cols <- colnames(test.df[3:9])
test.df[cols] <- log(test.df[cols])

#create pca
pca_more_traits <- princomp(test.df[3:9])
fviz_pca_biplot(pca_more_traits, col.ind = test.df$func_group)

# for figures: take out the centroid, take out axes and gridlines, but put
# a box around the plot, and get rid of the title, put the legend outside of the
# figure, change the errors to black 

# right now these are just PCAs of the main three (D, SRL, RTD) and of most of the
# traits, which ordination should I do next?

#constrained ordination with functional groups as a constraint (db-RDA)
# going to answer the same question that a PERMANOVA would: 
#to what extent do these groups drive variation in these traits?
