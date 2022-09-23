# this is the cleaned up R script with all of the statistics on the root traits
# 2 ANOVAs were run for each of the following traits:
# diameter, SRL, RTD, SLA, RDMC, branching_intensity, root_shoot_ratio, and above_to_below_scan
# One ANOVA was to compare functional groups (Forbs, N-fixers, Graminoids), whereas the second was to
# compare monocots and eudicots

# CLEAN UP
rm(list=ls(all=TRUE))


# load the root traits file
roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\all_root_traits_summarized.csv")


# to create the anovas of functional groups, we do not want the species that are shrubs,
# we need to create a new dataframe without those species
roots_noSS.df <- subset(roots.df, func_group!='SS' & func_group!='SSN')
View(roots_noSS.df)

# for diameter, SRL, and RTD, root_shoot_ratio, above_below_scan, branching_intensity
#the log transformation gave the most normal data
if(TRUE){
  diam <- log(roots.df$diameter + 0.01)
  SRL <- log(roots.df$SRL +0.01)
  rtd <- log(roots.df$RTD + 0.01)
  root_shoot <- log(roots.df$root_shoot_ratio + 0.01)
  above_below_scan <- log(roots.df$above_to_below_scan + 0.01)
  branching <- log(roots.df$branching_intensity + 0.01)
  # for SLA and RDMC, use the square root transformation
  sla <- sqrt(roots.df$SLA)
  rdmc <- sqrt(roots.df$RDMC)
  
  #creating the same variables but using the dataframe without the shrubs,
  # to use in the ANOVA comparing functional groups
  diam_noSS <- log(roots_noSS.df$diameter + 0.01)
  SRL_noSS <- log(roots_noSS.df$SRL + 0.01)
  RTD_noSS <- log(roots_noSS.df$RTD + 0.01)
  root_shoot_noSS <- log(roots_noSS.df$root_shoot_ratio + 0.01)
  above_below_noSS <- log(roots_noSS.df$above_to_below_scan + 0.01)
  branching_noSS <- log(roots_noSS.df$branching_intensity + 0.01)
  
  sla_noSS <- sqrt(roots_noSS.df$SLA)
  rdmc_noSS <- sqrt(roots_noSS.df$RDMC)
  
}

#creating the cots (monocot or dicot) variable
cots <- roots.df$cot

#creating the functional group variable
group <- roots_noSS.df$func_group

# creating ANOVA for diameter by cot
diameter_cot.aov <- aov(diam ~ cots)
qqnorm(residuals(diameter_cot.aov))
qqline(residuals(diameter_cot.aov))
plot(residuals(diameter_cot.aov) ~ fitted(diameter_cot.aov))
boxplot(diam ~ cots) # D > M
summary(diameter_cot.aov) # p-val = 6.47e-07

#creating ANOVA for diameter by functional group
diameter_func.aov <- aov(diam_noSS ~ group)
qqnorm(residuals(diameter_func.aov))
qqline(residuals(diameter_func.aov))
plot(residuals(diameter_func.aov) ~ fitted(diameter_func.aov))
boxplot(diam_noSS ~ group)
summary(diameter_func.aov)
TukeyHSD(diameter_func.aov)
# G < F p-val = 0
# N > F, p-val = 0.7474
# N > G, p-val = 0.000001

#create ANOVA for SRL by cot
SRL_cot.aov <- aov(SRL ~ cots)
qqnorm(residuals(SRL_cot.aov))
qqline(residuals(SRL_cot.aov))
plot(residuals(SRL_cot.aov) ~ fitted(SRL_cot.aov))
boxplot(SRL ~ cots) # M > D
summary(SRL_cot.aov) # p-val = 5.3e-05

# create ANOVA for SRL by func group
SRL_func.aov <- aov(SRL_noSS ~ group)
qqnorm(residuals(SRL_func.aov))
qqline(residuals(SRL_func.aov))
plot(residuals(SRL_func.aov) ~ fitted(SRL_func.aov))
boxplot(SRL_noSS ~ group)
summary(SRL_func.aov)
TukeyHSD(SRL_func.aov)
# G > F p-val = 0.0010400
# N < F p-val = 0.00335564
# N < G p-val = 0.00000013


#create ANOVA for RTD by cot
RTD_cot.aov <- aov(rtd ~ cots)
qqnorm(residuals(RTD_cot.aov))
qqline(residuals(RTD_cot.aov))
plot(residuals(RTD_cot.aov) ~ fitted(RTD_cot.aov))
boxplot(rtd ~ cots) # no sig dif
summary(RTD_cot.aov) # p-val = 0.15

# create ANOVA for RTD by func group
RTD_func.aov <- aov(RTD_noSS ~ group)
qqnorm(residuals(RTD_func.aov))
qqline(residuals(RTD_func.aov))
plot(residuals(RTD_func.aov) ~ fitted(RTD_func.aov))
boxplot(RTD_noSS ~ group)
summary(RTD_func.aov)
TukeyHSD(RTD_func.aov)
# N > G p-val = 0.0000336
# N > F p-val = 0
# G = F, p-val = 0.569

#create ANOVA for root_shoot by cot
root_shoot_cot.aov <- aov(root_shoot ~ cots)
qqnorm(residuals(root_shoot_cot.aov))
qqline(residuals(root_shoot_cot.aov))
plot(residuals(root_shoot_cot.aov) ~ fitted(root_shoot_cot.aov))
boxplot(root_shoot ~ cots) # M < D
summary(root_shoot_cot.aov) # p-val = 0.0286

# create ANOVA for root_shoot by func group
root_shoot_func.aov <- aov(root_shoot_noSS ~ group)
qqnorm(residuals(root_shoot_func.aov))
qqline(residuals(root_shoot_func.aov))
plot(residuals(root_shoot_func.aov) ~ fitted(root_shoot_func.aov))
boxplot(root_shoot_noSS ~ group)
summary(root_shoot_func.aov)
TukeyHSD(root_shoot_func.aov)
# N > F p-val = 0.0002379
# N > G p-val = 0.0000067
# G < F p-val = 0.0648626

#create ANOVA for above_below_scan by cot
above_below_cot.aov <- aov(above_below_scan ~ cots)
qqnorm(residuals(above_below_cot.aov))
qqline(residuals(above_below_cot.aov))
plot(residuals(above_below_cot.aov) ~ fitted(above_below_cot.aov))
boxplot(above_below_scan ~ cots) # M > D (no sig dif)
summary(above_below_cot.aov) # p-val = 0.0672

# create ANOVA for above_below_scan by func group
above_below_func.aov <- aov(above_below_noSS ~ group)
qqnorm(residuals(above_below_func.aov))
qqline(residuals(above_below_func.aov))
plot(residuals(above_below_func.aov) ~ fitted(above_below_func.aov))
boxplot(above_below_noSS ~ group)
summary(above_below_func.aov)
TukeyHSD(above_below_func.aov)
# G > F p-val = 0.0279878
# N < F p-val = 0.1266468
# N < G p-val = 0.0020714

#create ANOVA for branching intensity by cot
branching_cot.aov <- aov(branching ~ cots)
qqnorm(residuals(branching_cot.aov))
qqline(residuals(branching_cot.aov))
plot(residuals(branching_cot.aov) ~ fitted(branching_cot.aov))
boxplot(branching ~ cots) # M > D
summary(branching_cot.aov) # p-val = 2.49e-09

# create ANOVA for branching intensity by func group
branching_func.aov <- aov(branching_noSS ~ group)
qqnorm(residuals(branching_func.aov))
qqline(residuals(branching_func.aov))
plot(residuals(branching_func.aov) ~ fitted(branching_func.aov))
boxplot(branching_noSS ~ group)
summary(branching_func.aov)
TukeyHSD(branching_func.aov)
# G > F p-val = 0
# N > F p-val = 0.0220035
# N < G p-val = 0.0000168

#create ANOVA for SLA by cot
SLA_cot.aov <- aov(sla ~ cots)
qqnorm(residuals(SLA_cot.aov))
qqline(residuals(SLA_cot.aov))
plot(residuals(SLA_cot.aov) ~ fitted(SLA_cot.aov))
boxplot(sla ~ cots) # D > M
summary(SLA_cot.aov) # p-val = 3.92e-06

# create ANOVA for SLA by func group
sla_func.aov <- aov(sla_noSS ~ group)
qqnorm(residuals(sla_func.aov))
qqline(residuals(sla_func.aov))
plot(residuals(sla_func.aov) ~ fitted(sla_func.aov))
boxplot(sla_noSS ~ group)
summary(sla_func.aov)
TukeyHSD(sla_func.aov)
# G < F p-val = 0.000387
# N > F p-val = 0.5761417
# N > G p-val = 0.0011002

#create ANOVA for RDMC by cot
rdmc_cot.aov <- aov(rdmc ~ cots)
qqnorm(residuals(rdmc_cot.aov))
qqline(residuals(rdmc_cot.aov))
plot(residuals(rdmc_cot.aov) ~ fitted(rdmc_cot.aov))
boxplot(rdmc ~ cots) # no sig dif
summary(rdmc_cot.aov) # p-val 0.348

# create ANOVA for RDMC by func group
rdmc_func.aov <- aov(rdmc_noSS ~ group)
qqnorm(residuals(rdmc_func.aov))
qqline(residuals(rdmc_func.aov))
plot(residuals(rdmc_func.aov) ~ fitted(rdmc_func.aov))
boxplot(rdmc_noSS ~ group)
summary(rdmc_func.aov)
TukeyHSD(rdmc_func.aov)
# G > F p-val = 0.0000157
# N > F p-val = 0
# N > G p-val = 0.0113185


