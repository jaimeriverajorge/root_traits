# R-script to create scatter plots of 8 regressions:
# SRL x Diameter for  all  plants, then Forbs, Graminoids, and N-fixers separately
# RTD x N for all plants, then Forbs, Graminoids, and N-fixers separately

# added code to make regressions of:
# SRL x RTD for all plants, as well as forbs, graminoids, and n-fixers (expect a negative correlation)
# SRL x SLA for all plants, as well as forbs, graminoids, and n-fixers (expect a positive correlation)

# CLEAN UP
rm(list=ls(all=TRUE))

#Load packages
library(vegan)
library(ggplot2)
library(caper)

# originally this script was to run a PGLS on each of our other traits, but there was almost always
# no difference between a PGLS and a regular regression, so we continued with just regular 
# regressions and linear models to get the R-squared and p-values, much of the commented out code
# corresponds to the PGLS analysis

#read in phylogenetic tree
#tree_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre"
#phy <- read.tree(tree_path)
#plot(phy)
# set node labels to NULL
#phy$node.label <- NULL

# get trait data
data_path = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\data\\root_traits_summarized.csv"
# file path to the root trait data if at home computer
#data_path = "C:\\Users\\Jorge\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv"
roots.df <- read.csv(data_path)

#shorten roots.df to only include species, func_group, diameter, SRL, RTD
roots.df <- roots.df[c('species', 'func_group', 'diameter', 'SRL', 'RTD')]

#log transform data
roots.df$diameter <- log(roots.df$diameter + 0.0001)
roots.df$SRL <- log(roots.df$SRL + 0.0001)
roots.df$RTD <- log(roots.df$RTD + 0.0001)

# create dataframe of just forbs
roots_F.df <- subset(roots.df, func_group=='F')
#create comparative data object for forbs
#prairie_forbs <- comparative.data(data = roots_F.df, phy=phy, names.col = 'species', vcv = TRUE, na.omit = FALSE)

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
#forbs.pgls <- pgls(srl_F ~ diameter_F, data = prairie_forbs, lambda = 'ML')
#summary(forbs.pgls)
# Multiple R-squared: 0.5673,	Adjusted R-squared: 0.5611
# p-value: 2.31e-14 

#make the nice plot

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

# make linear model of SRL vs RTD for Forbs
srl_rtd_forbs.lm <- lm(srl_F ~ rtd_F)
summary(srl_rtd_forbs.lm)
# Multiple R-squared:  0.6465,	Adjusted R-squared:  0.6414 
# p-value: < 2.2e-16

#make plot of SRL vs RTD for forbs
forbs_plot_srl_rtd <- ggplot(data = roots_F.df, aes(x = RTD, y = SRL))+ 
  geom_point(color='#E69F00')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
forbs_plot_srl_rtd
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\srl_rtd_forb_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(forbs_plot_srl_rtd)
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

# linear model of all plants (SRL vs Diameter)
srl_all.lm <- lm(roots.df$SRL ~ roots.df$diameter)
summary(srl_all.lm)
# Multiple R-squared:  0.5015,	Adjusted R-squared:  0.4966 
# p-value: < 2.2e-16


#make plot of all plants SRL vs RTD
all_plot_srl_rtd <- ggplot(data = roots.df, aes(x = RTD, y = SRL))+ 
  geom_point(color='black')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
all_plot_srl_rtd
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\srl_rtd_all_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(all_plot_srl_rtd)
dev.off()

# make linear model of SRL vs RTD for all plants
srl_rtd_all.lm <- lm(roots.df$SRL ~ roots.df$RTD)
summary(srl_rtd_all.lm)
#Multiple R-squared:  0.5623,	Adjusted R-squared:  0.558 
# p-value: < 2.2e-16

#create subset of just N-fixers
roots_N.df <- subset(roots.df, func_group=='N')

# linear model of SRL x D for N-fixers
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

#create plot of SRL vs RTD for N-fixers
N_plot_srl_rtd <- ggplot(data = roots_N.df, aes(x = RTD, y = SRL))+ 
  geom_point(color='#56B4E9')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
N_plot_srl_rtd
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\srl_rtd_nfix_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(N_plot_srl_rtd)
dev.off()

#linear model of SRL vs RTD for N-fixers
srl_rtd_n.lm <- lm(roots_N.df$SRL ~ roots_N.df$RTD)
summary(srl_rtd_n.lm)
# Multiple R-squared:  0.6307,	Adjusted R-squared:  0.5938 
# p-value: 0.002036


#create subset of just graminoids
roots_G.df <- subset(roots.df, func_group=='G')

# create linear model of SRL x D for Graminoids
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

#create plot of SRL vs RTD for graminoids
G_plot_srl_rtd <- ggplot(data = roots_G.df, aes(x = RTD, y = SRL))+ 
  geom_point(color='#009E73')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Root Tissue Density", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
G_plot_srl_rtd
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\srl_rtd_G_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(G_plot_srl_rtd)
dev.off()

#linear model of SRL vs RTD for graminoids
gram_srl_rtd.lm <- lm(roots_G.df$SRL ~ roots_G.df$RTD)
summary(gram_srl_rtd.lm)
# Multiple R-squared:  0.5257,	Adjusted R-squared:  0.4941 
# p-value: 0.0009898


#Now we create plots of a regression of RTD vs %N, then of SRL vs SLA
# load in data sheet that has the %N data
n_data <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\data\\root_traits_summarized.csv")
# we just want the following columns: species, functional group, RTD, %N
n_data <- n_data[c('species', 'func_group', 'RTD', 'N_roots', 'SRL', 'SLA')]
# get rid of NA values
n_data <- na.omit(n_data)

# log transform data
n_data$RTD <- log(n_data$RTD + 0.0001)
n_data$N_roots <- log(n_data$N_roots + 0.0001)
n_data$SRL <- log(n_data$SRL + 0.0001)
n_data$SLA <- sqrt(n_data$SLA)

# the same specimen that has no N data has no SLA data, so we can still use the above 
# dataframe when making regressions of SRL and SLA

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

#linear model of SRL vs SLA for all plants
sla_all.lm <- lm(n_data$SRL ~ n_data$SLA)
summary(sla_all.lm)
#Multiple R-squared:  0.003494,	Adjusted R-squared:  -0.006372 
#p-value: 0.5531
sla_plot_all <- ggplot(data = n_data, aes(x = SLA, y = SRL))+ 
  geom_point(color='black')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Specific Leaf Area", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
sla_plot_all
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\sla_all_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(G_plot_srl_rtd)
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

#linear model of SRL vs SLA for forbs
sla_F.lm <- lm(n_data_f$SRL ~ n_data_f$SLA)
summary(sla_F.lm)
#Multiple R-squared:  0.01487,	Adjusted R-squared:  0.0005929 
#p-value: 0.311


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

#linear model of SRL vs SLA for n-fixers
sla_N.lm <- lm(n_data_N$SRL ~ n_data_N$SLA)
summary(sla_N.lm)
#Multiple R-squared:  0.3685,	Adjusted R-squared:  0.3054
# p-value: 0.03633

# plot of SRL vs SLA for n-fixers
sla_plot_N <- ggplot(data = n_data_N, aes(x = SLA, y = SRL))+ 
  geom_point(color='#56B4E9')+
  labs(y=expression(atop("Specific Root Length", paste(~log(cm/g)))),
       x=expression(atop("Specific Leaf Area", paste(~log(g/cm^{3})))))+
  geom_smooth(method = 'lm', formula = y~x, color='black', size=0.5)+
  theme_classic()
sla_plot_N
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC (3)\\Documents\\figures\\sla_N_regression.png"
png(file=file_name,height=3.5,width=4.5,units="in",res=1200)
plot(sla_plot_N)
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

# create linear model of SRL vs SLA for graminoids
sla_G.lm <- lm(n_data_G$SRL ~ n_data_G$SLA)
summary(sla_G.lm)
# Multiple R-squared:  0.01051,	Adjusted R-squared:  -0.05546
# p-value: 0.6954

