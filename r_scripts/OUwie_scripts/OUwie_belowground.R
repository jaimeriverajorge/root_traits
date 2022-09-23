# this script uses the R packag OUwie to determine trait optima for different
# growth forms related to non-root belowground organs:
# R: species with rhizomes
# B: species with storage organs (bulbs, corms, and tubers)
# P: perennial species with none of the above organs
# we would have two more, but the dataset only had one of each so we remove these
# from our data:
# A: annual species with none of the above organs (chamaecrista_fasciculata)
# S: species with stolons (rudbeckia_fulgida)

# there is one species: Helianthus strumosus (HELSTR) that has both rhizomes
# and tubers, so we will run the analysis once with it classified as R and once
# again with it classified as B

library(phytools)
library(OUwie)
library(picante)
library(purrr)

# CLEAN UP
rm(list=ls(all=TRUE))

# load in dataframe and phylogenetic tree
roots <- read.csv("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized.csv", header = TRUE, row.names = 1)
phy <- read.tree("C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\r_scripts\\phylogenetic_trees\\tr.analysis.tre")

#remove chamaecrista_fasciculata and rudbeckia_fulgida from our dataset
rows_to_remove <- c('chamaecrista_fasciculata','rudbeckia_fulgida')
roots <- roots[!(row.names(roots) %in% rows_to_remove),]

# make phy rooted and fully dichotomous by using multi2di
phy <- multi2di(phy)
# put data in the same order as phy tree
combined <- match.phylo.data(phy, roots)
phy <- combined$phy
roots <- combined$data

# change diameter, SRL, and RTD back into numbers, log transform at the same time
roots$diameter <- log(as.double(roots$diameter) + 0.001)
roots$SRL <- log(as.double(roots$SRL) + 0.0001)
roots$RTD <- log(as.double(roots$RTD) + 0.0001)

# remake the species name column
roots$species <- row.names(roots)



#change the row names back to numbers
row.names(roots) <- 1:nrow(roots)

#make separate dataframes for the diameter, SRL, and RTD
roots_diam <- roots[, c('species', 'organ', 'diameter')]
roots_SRL <- roots[, c('species', 'organ', 'SRL')]
roots_rtd <- roots[, c('species', 'organ', 'RTD')]

# starting with just the diameter data, we will begin the stocahastic mapping
# process

# set up data by making a vector of tip states and name them with species names
states_OU2 <- setNames(roots_diam[,2], roots_diam[,1])
states_OU2 <- states_OU2[phy$tip.label]
names(states_OU2)

# compare the AIC of two models: ER: equal rates, ARD: all rates different
AIC(ace(states_OU2, phy, type = 'discrete', model='ER')) # 221
AIC(ace(states_OU2, phy, type = 'discrete', model='ARD')) # 171

# set up data for OU3, will use three different models: ER, ARD, and 
# SYM: symmetric rates (rates are the same between any pairs of states)
OU3_er <- ace(states_OU2, phy, type = 'discrete', model = 'ER')
OU3_sym <- ace(states_OU2, phy, type = 'discrete', model = 'SYM')
OU3_ard <- ace(states_OU2, phy, type = 'discrete', model = 'ARD')

#print out the AIC values
map_dbl(list(OU3_er, OU3_sym, OU3_ard), AIC)
# 221.4397 183.9225 171.7899
# lowest value from ARD: 171

# conduct our SIMMAP analysis, using the ER model for both
nsims <- 200
#take time before
time_before_stoch = Sys.time()
simmap_OU2 <- make.simmap(phy, states_OU2, model = 'ER', nsim = nsims)
time_after_stoch = Sys.time()
OU2_OUwie_sm <- lapply(simmap_OU2, OUwie, data= roots_diam, model= 'OUM', simmap.tree = T)
time_after_looping = Sys.time()

# calculate time differences
time_to_simulate = time_after_stoch - time_before_stoch
print(time_to_simulate) # 2.387183 secs
time_to_loop = time_after_looping - time_after_stoch
print(time_to_loop) # 1.196404 mins for 10 trees, #23.40393 mins for 200 trees

# Conduct SIMMAP analysis using the SYM model
time_before_sim2 = Sys.time()
simmap_OU2_sym <- make.simmap(phy, states_OU2, model = 'SYM', nsim = nsims)
time_after_sim2 = Sys.time()
OU2_OUwie_sm_sym <- lapply(simmap_OU2_sym, OUwie, data= roots_diam, model= 'OUM', simmap.tree = T)
time_after_looping_sym = Sys.time()

time_to_simulate_SYM = time_after_sim2 - time_before_sim2
print(time_to_simulate_SYM) # 6.942846 secs
time_to_loop_SYM = time_after_looping_sym - time_after_sim2
print(time_to_loop_SYM) # 2.376897 mins