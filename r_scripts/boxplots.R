# This script is used to make box plots of the root traits, specifically of SRL,
# RTD, and diameter compared by functional group (Forbs, Graminoids, N-fixers),
# each plot includes the selective optimum of each functional group

# Selective Optima for SRL (log transformed on the left, normal on the right):
# for F: 8.1973958, 3631.481
# for G: 9.0249427, 8307.738
# for N: 7.3312831, 1527.34

# Selective Optima for RTD:
# F: -1.19474918, 0.3027799
# G: -1.0629458, 0.3454367
# N: -0.2062934 0.8135943

# Selective Optima for Diameter:
# F: -1.03508596, 0.3551958
# G: -1.52696514, 0.2171938
# N: -1.00129532, 0.3674032

# read in csv data file
roots.df <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv")
roots.df <- subset(roots.df, func_group!='SS' & func_group!='SSN')

# plotting SRL vs functional group:
boxplot(roots.df$SRL ~ roots.df$func_group, ylab = "Specific Root Length (cm/g)", 
        xlab="Functional Group", xaxt="n")
axis(1, at=1:3, labels = c("Forbs", "Graminoids", "N-Fixers"))
text(c(3631.481, 8307.738, 1527.34), "*", cex=2)

# plotting RTD vs functional group:
boxplot(roots.df$RTD ~ roots.df$func_group, ylab = "Root Tissue Density (g/cm3)", 
        xlab="Functional Group", xaxt='n')
axis(1, at=1:3, labels = c("Forbs", "Graminoids", "N-Fixers"))
text(c(0.3027799, 0.3454367, 0.8135943), "*", cex=2)

# plotting diameter vs functional group:
boxplot(roots.df$diameter ~ roots.df$func_group, ylab = "Root Diameter (mm)", 
        xlab="Functional Group", xaxt='n')
axis(1, at=1:3, labels = c("Forbs", "Graminoids", "N-Fixers"))
text(c(0.3551958, 0.2171938, 0.3674032), "*", cex = 2)


