# R-script to create scatter plots of 8 regressions:
# SRL x Diameter for  all  plants, then Forbs, Graminoids, and N-fixers separately
# RTD x N for all plants, then Forbs, Graminoids, and N-fixers separately

# CLEAN UP
rm(list=ls(all=TRUE))

#Load packages
library(vegan)
library(ggplot2)
library(caper)

#read in phylogenetic tree
tree_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre"
phy <- read.tree(tree_path)
#plot(phy)
# set node labels to NULL
phy$node.label <- NULL

# get trait data
data_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized.csv"
# file path to the root trait data if at home computer
#data_path = "C:\\Users\\Jorge\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv"
roots.df <- read.csv(data_path)

#shorten roots.df to only include species, func_group, diameter, SRL, RTD (and in the future %N)
roots.df <- roots.df[c('species', 'func_group', 'diameter', 'SRL', 'RTD')]

#log transform data
roots.df$diameter <- log(roots.df$diameter + 0.0001)
roots.df$SRL <- log(roots.df$SRL + 0.0001)
roots.df$RTD <- log(roots.df$RTD + 0.0001)

# create dataframe of just forbs
roots_F.df <- subset(roots.df, func_group=='F')
#create comparative data object for forbs
prairie_forbs <- comparative.data(data = roots_F.df, phy=phy, names.col = 'species', vcv = TRUE, na.omit = FALSE)

#recreate variables for plotting
diameter_F <- roots_F.df$diameter
srl_F <- roots_F.df$SRL
rtd_F <- roots_F.df$RTD

#forbs SRL vs D linear model
srl_forbs.lm <- lm(srl_F ~ diameter_F)
summary(srl_forbs.lm) 
#Multiple R-squared:  0.5673,	Adjusted R-squared:  0.5611
#p-val = 2.31e-14

#create pgls model
forbs.pgls <- pgls(srl_F ~ diameter_F, data = prairie_forbs, lambda = 'ML')
summary(forbs.pgls)
# Multiple R-squared: 0.5673,	Adjusted R-squared: 0.5611
# p-value: 2.31e-14 

# plotting the logged transformed data

plot(srl_F ~ diameter_F, data = roots_F.df)
abline(forbs.pgls)

# reverse log transform data to get the raw data that we will plot
diam_F_raw <- exp(diameter_F)
srl_F_raw <- exp(srl_F)

#make the dataframe include the raw data instead so we can plot it

roots_F.df$diameter <- log(roots_F.df$diameter + 0.0001)
roots_F.df$SRL <- log(roots_F.df$SRL + 0.0001)
# the PGLS was calculated with the log transformed data but we want to plot the 
# raw data on a log scale so that it is in the same units as the box plot
# I plotted it but I am not sure what I should do to do the PGLS abline to make it
# fit on to the raw data, or how to plot it
# I just reverse logged the intercept from to get the value for the raw data, is 
# this correct?

forbs_plot <- ggplot(data = roots_F.df, aes(x = diameter, y = SRL))+ 
  geom_point(color='#E69F00')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Diameter", paste(~log(mm)))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
forbs_plot
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\srl_forbs_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(forbs_plot)
dev.off()

#make plot of all plants SRL vs D
srl_plot_all <- ggplot(data = roots.df, aes(x = diameter, y = SRL))+ 
  geom_point()+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Diameter", paste(~log(mm)))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
srl_plot_all
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\srl_all_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(srl_plot_all)
dev.off()

# linear model of all plants
srl_all.lm <- lm(roots.df$SRL ~ roots.df$diameter)
summary(srl_all.lm)
# Multiple R-squared:  0.5015,	Adjusted R-squared:  0.4966 
# p-value: < 2.2e-16

#create subset of just N-fixers
roots_N.df <- subset(roots.df, func_group=='N')
srl_N.lm <- lm(roots_N.df$SRL ~ roots_N.df$diameter)
summary(srl_N.lm)
# Multiple R-squared:  0.1496,	Adjusted R-squared:  0.06458 
# p-value: 0.2142

#create plot of SRL vs D for N-fixers
srl_plot_N <- ggplot(data = roots_N.df, aes(x = diameter, y = SRL))+ 
  geom_point(color='#56B4E9')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Diameter", paste(~log(mm)))))+
  theme_classic()
srl_plot_N
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\srl_N_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(srl_plot_N)
dev.off()

#create subset of just graminoids
roots_G.df <- subset(roots.df, func_group=='G')
srl_G.lm <- lm(roots_G.df$SRL ~ roots_G.df$diameter)
summary(srl_G.lm)
# Multiple R-squared:  0.5533,	Adjusted R-squared:  0.5235
# p-value: 0.000619
#create plot of SRL vs D for graminoids
srl_plot_G <- ggplot(data = roots_G.df, aes(x = diameter, y = SRL))+ 
  geom_point(color='#009E73')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Diameter", paste(~log(mm)))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
srl_plot_G
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\srl_G_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(srl_plot_G)
dev.off()


#Now we create plots of a regression of RTD vs %N
# load in data sheet that has the %N data
n_data <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized.csv")
# we just want the following columns: species, functional group, RTD, %N
n_data <- n_data[c('species', 'func_group', 'RTD', 'N_roots')]
# get rid of NA values
n_data <- na.omit(n_data)

#linear model of RTD vs N for all plants
rtd_all.lm <- lm(n_data$RTD ~ n_data$N_roots)
summary(rtd_all.lm)
# Multiple R-squared:  2.727e-05,	Adjusted R-squared:  -0.009873 
# p-value: 0.9582

#Create plot of RTD vs N for all plants
rtd_plot_all <- ggplot(data = n_data, aes(x = N_roots, y = RTD))+ 
  geom_point()+
  labs(y=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))),
       x=expression(atop("Root [N]", paste(~log('%')))))+
  theme_classic()
rtd_plot_all
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\rtd_all_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(rtd_plot_all)
dev.off()

#create subset of n_data for just Forbs
n_data_f <- subset(n_data, func_group=='F')

#linear model of RTD vs N for forbs
rtd_F.lm <- lm(n_data_f$RTD ~ n_data_f$N_roots)
summary(rtd_F.lm)
# Multiple R-squared:  0.07512,	Adjusted R-squared:  0.06171 
# p-value: 0.02073
# create plot of RTD vs N for Forbs
rtd_plot_F <- ggplot(data = n_data_f, aes(x = N_roots, y = RTD))+ 
  geom_point(color='#E69F00')+
  labs(y=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))),
       x=expression(atop("Root [N]", paste(~log('%')))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
rtd_plot_F
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\rtd_forb_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(rtd_plot_F)
dev.off()

#create subset of n_data for just N-fixers
n_data_N <- subset(n_data, func_group=='N')

#linear model of RTD vs N for N-fixers
rtd_N.lm <- lm(n_data_N$RTD ~ n_data_N$N_roots)
summary(rtd_N.lm)
# p-value: 0.7186
#create plot of RTD vs N for N-fixers
rtd_plot_N <- ggplot(data = n_data_N, aes(x = N_roots, y = RTD))+ 
  geom_point(color='#56B4E9')+
  labs(y=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))),
       x=expression(atop("Root [N]", paste(~log('%')))))+
  theme_classic()
rtd_plot_N
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\rtd_N_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(rtd_plot_N)
dev.off()

#create subset of n_data for just graminoids
n_data_G <- subset(n_data, func_group=='G')

# create linear model of RTD vs N for graminoids
rtd_G.lm <- lm(n_data_G$RTD ~ n_data_G$N_roots)

summary(rtd_G.lm)
# p-value: 0.7565

#create plot of RTD vs N for graminoids
rtd_plot_G <- ggplot(data = n_data_G, aes(x = N_roots, y = RTD))+ 
  geom_point(color='#009E73')+
  labs(y=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))),
       x=expression(atop("Root [N]", paste(~log('%')))))+
  theme_classic()
rtd_plot_G
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\rtd_G_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(rtd_plot_G)
dev.off()




