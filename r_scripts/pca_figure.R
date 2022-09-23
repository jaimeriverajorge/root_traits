# this R script creates a nice clean figure of an ordination 

# remove variables from environment
rm(list=ls(all=TRUE))

# load in packages
library(vegan)
library(ggfortify)
library(factoextra)

# read in data
roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized.csv")
#remove the shrubs
roots.df <- subset(roots.df, func_group != 'SS' & func_group != 'SSN')
#just get the columns we want (species, func_group, diameter, SRL, RTD)
roots.df <- roots.df[c('species', 'func_group', 'diameter', 'SRL', 'RTD', 'N_roots')]
roots.df <- na.omit(roots.df)
#log transform our variables
cols <- colnames(roots.df[3:6])
roots.df[cols] <- log(roots.df[cols] + 0.0001)

#create pca
my_pca <- princomp(roots.df[3:6])

#plot with fviz_pca_biplot
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\pca_4_5.png"
png(file=file_name,height=4,width=5,units="in",res=1200)
fviz_pca_biplot(my_pca, palette =c("#E69F00", "#009E73", "#56B4E9"), 
                habillage = roots.df$func_group, addEllipses = TRUE,
                title = '', col.var = 'black', label = 'none', 
               invisible = c('axes'))
dev.off()
# problem with the above plot is I cannot take out the centroid or the 
# grid lines or axes or legend
# TODO: save figure into png with high resolution, put into inkscape
#try with autoplot
autoplot(my_pca, data = roots.df, colour = 'func_group', loadings=TRUE, 
         loadings.colour='black', loadings.label=TRUE, loadings.label.colour='black',
         loading.label.size= 5, frame=TRUE, frame.type='norm', legendLabs='none')
# cant figure out how to set the colors for each functional group




# for figures: take out the centroid, take out axes and gridlines, but put
# a box around the plot, and get rid of the title, put the legend outside of the
# figure, change the errors to black 