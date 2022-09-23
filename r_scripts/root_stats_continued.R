# this file continues statistics from the root traits data sheet
# making histograms and anovas of: Specific leaf area (SLA), root dry matter content
# (RDMC), root to shoot ratio (root_shoot), aboveground to belowground ratio
# for the scanned indiv (above_below_scan), and aboveground to belowground ratio
# for all individuals (above_below_total), and branching intensity (branching)

# CLEAN UP
rm(list=ls(all=TRUE))

#root traits data frame
roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\root_traits_master_sheet_fixed.csv")
sla <- roots.df$SLA
#remove zero values from SLA
sla <- sla[sla != 0]
hist(sla)

# find Q1, Q3 and interquartile range for values in SLA
Q1 <- quantile(roots.df$SLA, .25)
Q3 <- quantile(roots.df$SLA, .75)
IQR <- IQR(roots.df$SLA)

# remove outliers for SRL, creates a new dataframe so need to assign SRL, cots, and func again
#ARNPLA4 is the biggest outlier, how did it have such high leaf area but low weight?
noOut <- subset(roots.df, roots.df$SLA> (Q1 - 1.5*IQR) & roots.df$SLA< (Q3 + 1.5*IQR))

#reassigning sla
sla <- noOut$SLA
#remove zero values from sla
sla <- sla[sla != 0]
hist(sla)
# testing for normality of SLA
shapiro.test(sla) #p-value = 0.002766
# see what the residuals look like before transforming

# doing the same as above but for root dry matter content (RDMC)
rdmc <- na.omit(roots.df$RDMC)
# removing zero values
rdmc <- rdmc[rdmc != 0]
hist(rdmc)
shapiro.test(rdmc) # p-value = 4.964e-08

#root_shoot, above_below_scan, above_below_total, branching
root_shoot <- roots.df$root_shoot_ratio
root_shoot <- root_shoot[root_shoot != 0] #getting rid of zero values
hist(root_shoot)
shapiro.test(root_shoot) # p-value < 2.2e-16

above_below_scan <- roots.df$above_to_below_scan
above_below_scan <- above_below_scan[above_below_scan !=0]
hist(above_below_scan)
shapiro.test(above_below_scan) # p-value < 2.2e-16

above_below_total <- roots.df$above_to_below_total
above_below_total <- above_below_total[above_below_total != 0]
hist(above_below_total)
shapiro.test(above_below_total) # p-value = 1.472e-15

branching <- roots.df$branching_intensity
hist(branching)
shapiro.test(branching) # p-value = 2.566e-13

#taking square root of data then testing for normality
if(TRUE) {
  sqSLA <- sqrt(sla)
  sqRDMC <- sqrt(rdmc)
  sqRootShoot <- sqrt(root_shoot)
  sqAboveBelowScan <- sqrt(above_below_scan)
  sqAboveBelowTotal <- sqrt(above_below_total)
  sqBranching <- sqrt(branching)
}
hist(sqSLA) # use this one for the ANOVA
hist(sqRDMC)# use this one
hist(sqRootShoot)
hist(sqAboveBelowScan)
hist(sqAboveBelowTotal)
hist(sqBranching)
shapiro.test(sqSLA) # p-value = 0.2198
shapiro.test(sqRDMC) # p-value = 0.009243
shapiro.test(sqRootShoot) # p-value = 1.876e-12
shapiro.test(sqAboveBelowScan) # p-value = 3.687e-08
shapiro.test(sqAboveBelowTotal) # p-value = 2.961e-05
shapiro.test(sqBranching) # p-value = 1.573e-06


#taking log of data then testing for normality
if(TRUE) {
  logSLA <- log(sla + 0.001)
  logRDMC <- log(rdmc + 0.001)
  logRootShoot <- log(root_shoot + 0.001)
  logAboveBelowScan <- log(above_below_scan + 0.001)
  logAboveBelowTotal <- log(above_below_total + 0.001)
  logBranching <- log(branching + 0.001)
}
hist(logSLA)
hist(logRDMC)
hist(logRootShoot)
hist(logAboveBelowScan)
hist(logAboveBelowTotal)
hist(logBranching)
shapiro.test(logSLA) # p-value = 1.908e-07
shapiro.test(logRDMC) # p-value = 0.001803
shapiro.test(logRootShoot) # p-value = 0.01369
shapiro.test(logAboveBelowScan) # p-value = 0.01693
shapiro.test(logAboveBelowTotal) # p-value = 7.028e-05
shapiro.test(logBranching) # p-value = 0.2023

# for SLA and RDMC, use square rooted version for ANOVAs
# for RootShoot, AboveBelowScan, Branching, use log transformed for ANOVAS
# could not get a normal distribution for above_below_total, may have to 
# remove outliers first

cots <- roots.df$cot
group <- roots.df$func_group

#ANOVAS for SLA vs monocot
slaCots.aov <- aov(sqSLA ~ cots) #this gives an error since cots has a
# different length than SLA because of getting rid of the 0 values,
# will need to create a subset of the datasheet instead
noZeroSLA.df <- subset(roots.df, SLA!=0)
cotSLA <- noZeroSLA.df$cot
sla <- sqrt(noZeroSLA.df$SLA)
slaCots.aov <- aov(sla ~ cotSLA)
# plotting residuals of diameter 
qqnorm(residuals(slaCots.aov))
qqline(residuals(slaCots.aov))
plot(residuals(slaCots.aov) ~ fitted(slaCots.aov))
#printing output from the sla Cot ANOVA
boxplot(sla ~ cotSLA)
summary(slaCots.aov) # p-val = 3.92e-07 
# dicots have a larger SLA than monocots 
#SLA vs LMA?

#ANOVA for SLA and functional groups
# using the previously made subset that does not have SLA with zero values,
# but removing the SS and SSN functional groups 
noSS.df <- subset(noZeroSLA.df, func_group != 'SS')
noSS.df <- subset(noSS.df, func_group != 'SSN')
View(noSS.df)
func_group_SLA <- noSS.df$func_group
sla <- sqrt(noSS.df$SLA)
slaFunc.aov <- aov(sla ~ func_group_SLA)
# plotting residuals of diameter 
qqnorm(residuals(slaFunc.aov))
qqline(residuals(slaFunc.aov))
plot(residuals(slaFunc.aov) ~ fitted(slaFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(sla ~ func_group_SLA) # why do SS and SSN still appear if I made a 
# subset without them before declaring func_group_SLA ?
summary(slaFunc.aov)
TukeyHSD(slaFunc.aov) #  for 
#G-F, p < 0.0001 (F > G) 
#N-G, p < 0.0001 (N > G)


#ANOVAS for root dry matter content
noZeroRDMC.df <- subset(roots.df, RDMC!=0)
cotRDMC <- noZeroRDMC.df$cot
rdmc <- sqrt(noZeroRDMC.df$RDMC)
rdmcCots.aov <- aov(rdmc ~ cotRDMC)
# plotting residuals of diameter 
qqnorm(residuals(rdmcCots.aov))
qqline(residuals(rdmcCots.aov))
plot(residuals(rdmcCots.aov) ~ fitted(rdmcCots.aov))
#printing output from the rdmc Cot ANOVA
boxplot(rdmc ~ cotRDMC)
summary(rdmcCots.aov) #p-val = 0.104, no sig difference

#ANOVA for RDMC and functional group
noSSRDMC.df <- subset(noZeroRDMC.df, func_group!='SS' & func_group!='SSN')
funcRDMC <- noSSRDMC.df$func_group
rdmcNoSS <- sqrt(noSSRDMC.df$RDMC)
rdmcFunc.aov <- aov(rdmcNoSS ~ funcRDMC)
qqnorm(residuals(rdmcFunc.aov))
qqline(residuals(rdmcFunc.aov))
plot(residuals(rdmcFunc.aov) ~ fitted(rdmcFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(rdmcNoSS ~ funcRDMC)
summary(rdmcFunc.aov) #p-val = <2e-16
TukeyHSD(rdmcFunc.aov)
# G-F p-val < 0.001, G > F
# N-F p-val = 0, N > F
# N-G p-val = 0.006, N > G

#ANOVAS for root to shoot ratio
noZeroRootShoot.df <- subset(roots.df, roots.df$root_shoot_ratio!=0)
rootShoot <- log(noZeroRootShoot.df$root_shoot_ratio + 0.001)
rsCot <- noZeroRootShoot.df$cot
rootShootCot.aov <- aov(rootShoot ~ rsCot)
qqnorm(residuals(rootShootCot.aov))
qqline(residuals(rootShootCot.aov))
plot(residuals(rootShootCot.aov) ~ fitted(rootShootCot.aov))
#printing output from the sla Cot ANOVA
boxplot(rootShoot ~ rsCot)
summary(rootShootCot.aov) #p-val = 0.00135, D > M

#need to remove the shrubs from the dataset
noSSRootShoot.df <- subset(noZeroRootShoot.df, func_group!='SS' & func_group!='SSN')
rsNoSS <- log(noSSRootShoot.df$root_shoot_ratio + 0.001)
funcRootShoot <- noSSRootShoot.df$func_group
rootShootFunc.aov <- aov(rsNoSS ~ funcRootShoot)
qqnorm(residuals(rootShootFunc.aov))
qqline(residuals(rootShootFunc.aov))
plot(residuals(rootShootFunc.aov) ~ fitted(rootShootFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(rsNoSS ~ funcRootShoot)
summary(rootShootFunc.aov)
TukeyHSD(rootShootFunc.aov)
#N-F p-val < 0.001, N > F
#N-G p-val < 0.001, N > G
#G-F p-val = 0.015, F > G

#ANOVAS for above to below for the scanned individual
noZeroAboveScan.df <- subset(roots.df, roots.df$above_to_below_scan != 0)
aboveBelowScan <- log(noZeroAboveScan.df$above_to_below_scan + 0.001)
aboveScanCot <- noZeroAboveScan.df$cot
aboveBelowScanCot.aov <- aov(aboveBelowScan ~ aboveScanCot)
qqnorm(residuals(aboveBelowScanCot.aov))
qqline(residuals(aboveBelowScanCot.aov))
plot(residuals(aboveBelowScanCot.aov) ~ fitted(aboveBelowScanCot.aov))
#printing output from the sla Cot ANOVA
boxplot(aboveBelowScan ~ aboveScanCot)
summary(aboveBelowScanCot.aov) #p-val < 0.001, M > D

noSSAboveScan.df <- subset(noZeroAboveScan.df, func_group!='SS' & func_group!='SSN')
aboveScanNoSS <- log(noSSAboveScan.df$above_to_below_scan + 0.001)
noSSABSFunc <- noSSAboveScan.df$func_group
aboveScanFunc.aov <- aov(aboveScanNoSS ~ noSSABSFunc)
qqnorm(residuals(aboveScanFunc.aov))
qqline(residuals(aboveScanFunc.aov))
plot(residuals(aboveScanFunc.aov) ~ fitted(aboveScanFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(aboveScanNoSS ~ noSSABSFunc)
summary(aboveScanFunc.aov)
TukeyHSD(aboveScanFunc.aov)
#G-F p-val = 0.0012149, G > F
#N-F p-val = 0.0065308 N < F
#N-G p-val < 0.001, N < G

#TODO: create ANOVAs for aboveToBelowTotal and for Branching
noZeroABTotal.df <- subset(roots.df, roots.df$above_to_below_total!=0)
aboveBelowTotal <- log(noZeroABTotal.df$above_to_below_total + 0.001)
abTotalCot <- noZeroABTotal.df$cot
aboveBelowTotalCot.aov <- aov(aboveBelowTotal ~ abTotalCot)
qqnorm(residuals(aboveBelowTotalCot.aov))
qqline(residuals(aboveBelowTotalCot.aov))
plot(residuals(aboveBelowTotalCot.aov) ~ fitted(aboveBelowTotalCot.aov))
#printing output from the sla Cot ANOVA
boxplot(aboveBelowTotal ~ abTotalCot)
summary(aboveBelowTotalCot.aov) #p-val = 0.000976, M > D

noSSABTotal.df <- subset(noZeroABTotal.df, func_group!='SS' & func_group!='SSN')
abTotalNoSS <- log(noSSABTotal.df$above_to_below_total + 0.001)
abTotalFunc <- noSSABTotal.df$func_group
aboveBelowTotalFunc.aov <- aov(abTotalNoSS ~ abTotalFunc)
qqnorm(residuals(aboveBelowTotalFunc.aov))
qqline(residuals(aboveBelowTotalFunc.aov))
plot(residuals(aboveBelowTotalFunc.aov) ~ fitted(aboveBelowTotalFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(abTotalNoSS ~ abTotalFunc)
summary(aboveBelowTotalFunc.aov)
TukeyHSD(aboveBelowTotalFunc.aov)
#G-F p-val < 0.001, G > F
#N-F p-val = 0.00135, N < F
#N-G p-val < 0.001, N < G

#ANOVAs for branching intensity
#branching intensity did not have any zero values, so we can use the original
#roots.df dataframe, where cot is already declared and so is logBranching
branchingCot.aov <- aov(logBranching ~ cots)
qqnorm(residuals(branchingCot.aov))
qqline(residuals(branchingCot.aov))
plot(residuals(branchingCot.aov) ~ fitted(branchingCot.aov))
#printing output from the sla Cot ANOVA
boxplot(logBranching ~ cots)
summary(branchingCot.aov) # p-val < 0.001, M > D

# we do still have to create a subset to get rid of the shrubs
noShrubs.df <- subset(roots.df, func_group!='SS' & func_group!='SSN')
noSSBranching <- log(noShrubs.df$branching_intensity + 0.001)
noSSFunc <- noShrubs.df$func_group
branchingFunc.aov <- aov(noSSBranching ~ noSSFunc)
qqnorm(residuals(branchingFunc.aov))
qqline(residuals(branchingFunc.aov))
plot(residuals(branchingFunc.aov) ~ fitted(branchingFunc.aov))
#printing output from the sla Cot ANOVA
boxplot(noSSBranching ~ noSSFunc)
summary(branchingFunc.aov)
TukeyHSD(branchingFunc.aov)
#G-F p-val = 0, G > F
#N-F p-val < 0.001, N > F
#N-G p-val = 0, N < G