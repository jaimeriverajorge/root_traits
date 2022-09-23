# This R-script is used to look at the data for the nitrogen concentration
# of the roots of prairie plants, and see its normality along with any outliers

#CLEAN UP
rm(list=ls(all=TRUE))

#read in root traits file
file_path <- "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_belowground_organs.csv"
roots.df <- read.csv(file_path)

# get rid of NA values
N_conc <- na.omit(roots.df$N_roots)

hist(N_conc)
shapiro.test(N_conc) # p-val = 0.0149, not normal but not terrible

# try log transforming
N_log <- log(N_conc + 0.001)
shapiro.test(N_log) #p-val = 0.1343
hist(N_log)

#log transforming seems to make our data normal, now lets make some box plots
# by functional group

#create a subset of the data with just func_group, and N_conc, remove shrubs,
# then remove NA
roots_N <- roots.df[c('species', 'func_group', 'N_roots')]
roots_N <- subset(roots_N, func_group != 'SS' & func_group != 'SSN')
roots_N <- na.omit(roots_N)
#log transform N data
roots_N$N_roots <- log(roots_N$N_roots + 0.0001)
# create boxplot of functional group
boxplot(roots_N$N_roots ~ roots_N$func_group)

#create ANOVA for N_root by functional group
N_func.aov <- aov(roots_N$N_roots ~ roots_N$func_group)
qqnorm(residuals(N_func.aov))
qqline(residuals(N_func.aov))
plot(residuals(N_func.aov) ~ fitted(N_func.aov))
summary(N_func.aov) #p-val < 0.0001
TukeyHSD(N_func.aov)
# G < F, p-val = 0
# N = F, p-val = 0.9637
# G < N, p-val = 0

# conduct a regression of RTD by N_root
#create a data sheet of just species, RTD, and N_root
roots_rtd <- roots.df[c('species', 'RTD', 'N_roots')]
# get rid of NA
roots_rtd <- na.omit(roots_rtd)

#create a linear model of rtd by N
rtd <- log(roots_rtd$RTD + 0.0001)
n <- log(roots_rtd$N_roots + 0.0001)
rtd_N.lm <- lm(rtd ~ n)
summary(rtd_N.lm)

plot(rtd ~ n)
abline(rtd_N.lm)

# look within functional groups
