# root_traits
Data and R code from the root traits project at the Morton Arboretum. September 2021 - September 2022

## Data_files: 
This folder has all of my spreadsheets that I used for my analyses. Below is a description of each file/folder.

OUwie_data: This folder has our trait data for a specific trait, organized in the following format which is required to be the input for OUwie: species, func_group, trait
	The species are also in the same order as our phylogenetic tree, and there is a file for the three traits: diameter, SRL, and RTD

Phylogenetic_trees: has four .tre files which were used in our analysis
	prairieTreeDichotomous.tre: the phylogenetic tree of our prairie species, with the shrubs taken out. This tree has the “selective regimes” coded into it on each node (F for forbs, N for n-fixers, G for graminoids) at each node in the tree. It was ran through the function “multi2di” from the “ape” package which “collapse or resolve multichotomies in phylogenetic trees” since OUwie needs a fully dichotomous tree.
	prairieTreeNotDichotomous.tre: same as the above tree, just without resolving the multichotomies.
	tr.acceptedNames.tre & tr.analysis.tre: both of these are the full phylogenetic trees of our prairie plants that Andrew sent me, with some slight naming differences between them. For example, in the acceptedNames tree, we have "Pycnanthemum_verticillatum_var._pilosum" , while in the analysis tree, we have "pycnanthemum_pilosum". For most of our data analyses I used the tr.analysis.tre file.

AMF_data: a CSV file with the mycorrhizal colonization data. The only three species missing are Allium canadense, Senna marilandica, and Senna hebecarpa. For the Allium, we did not have any roots leftover to make a slide with, while the two Senna species had roots that were too dark and thick to properly stain, even with the root lab’s methods of soaking in KOH for a week.

Pagels_lambda_values: a CSV file with pagel’s lambda and p-values for most of our traits.

PRISM_greenhouse_tissueCN: the elementar data with the N and C concentrations of our samples.

Root Trait Tin Labels and Weights:  a CSV file with the weights of all of our dried samples, written down straight off of the notebook. I later used this in the master sheet to organize them better, but this was the first data entry of the weights.

Root_traits_master_sheet: all of our root trait data, with the calculated values such as SRL, RTD, RDMC, and the weights, and notes on each individual, and the winrhizo data at the end of it as well. This contains all of the indivuals.

Root_traits_summarized: this CSV file has the “summarized” root trait data, meaning that there is only one row for each species. So the average was calculated when a species had multiple individuals harvested.

Roots_winRhizo: this CSV file has the root trait data taken straight from WinRhizo, organized by the separate scans. The master sheet took this data and compiled the different scans from the same individual and calculated the average values.

Species_list: this CSV file is a list of all of the individuals that I actually took data on.

Species_with_tap_roots: a list of species in our study that had tap roots.

## R_scripts
OUWie_scripts: This folder contains the different scripts that were used to run OUwie on the functional groups.

Boxplots: a script to make box plots of our different traits by functional group (SRL, RTD, diameter).

N_data: This R-script is used to look at the data for the nitrogen concentration of the roots of prairie plants, and see its normality along with any outliers.

Ordinations: A script that conducts ordinations on the root trait data, which includes using a db-RDA, and afterwards a regular PCA.

Pagels_lambda: this script runs the pagel’s lambda test on our different root traits.

Pca_figure: This script creates a nice clean figure of an ordination with our root trait data, which was used in my presentation and ESA poster.

Permanova: This script conducts a PERMANOVA analysis on the root trait data for diameter, specific root length, and root tissue density.

PGLS: includes code to run a PGLS regression on the root trait data. 

Pgls_figures: script to create scatter plots of 8 regressions: SRL x Diameter for all plants, then forbs, graminoids and n-fixers separately; RTD x N for all plants, then forbs, graminoids, and n-fixers separately.

Phylo_heatmap_figure: script to create a figure of a phylogenetic tree with a heatmap of the following traits: root diameter, specific root length, root tissue density, and root %N.

Phylo_stats: Preliminary statistics using the phylogenetic tree of the prarie plants taken from this tutorial: http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
and from Liam Revell's 2013 paper: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12066

Regression_figures: R-script to create scatter plots of 8 regressions:
SRL x Diameter for  all  plants, then Forbs, Graminoids, and N-fixers separately
RTD x N for all plants, then Forbs, Graminoids, and N-fixers separately
Later added code to make regressions of:
SRL x RTD for all plants, as well as forbs, graminoids, and n-fixers (expect a negative correlation)
SRL x SLA for all plants, as well as forbs, graminoids, and n-fixers (expect a positive correlation)

root_stats_ANOVAs: this is the cleaned up (as the root_stats_continued script was the first time running those statistics and included many unnecessary lines of code) R script with all of the statistics on the root traits, 2 ANOVAs were run for each of the following traits:
diameter, SRL, RTD, SLA, RDMC, branching_intensity, root_shoot_ratio, and above_to_below_scan
One ANOVA was to compare functional groups (Forbs, N-fixers, Graminoids), whereas the second was to compare monocots and eudicots.

Root_stats_continued: this file continues statistics from the root traits data sheet
making histograms and anovas of: Specific leaf area (SLA), root dry matter content
(RDMC), root to shoot ratio (root_shoot), aboveground to belowground ratio
for the scanned indiv (above_below_scan), and aboveground to belowground ratio
for all individuals (above_below_total), and branching intensity (branching)
