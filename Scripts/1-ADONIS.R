#############################################
# Analysis of Plover preen wax compositions #
#############################################
rm(list = ls())
dev.off()
# A. Load Data, source Packages & run GCalignR
# --------------------------------------------
library(GCalignR)
library(vegan)
library(ggplot2)
source("R/nmds_plotter.R")
source("R/nmds_calculator.R")

# The original text file (same name) was transformed to a RData file before
# load("data/charadrius_peaks.RData")

# "w41" is a septum-sample that was intended to control for effects of those on the profiles.
# Therefore a septum remove from sample 65 was washed several times in dichlormethane before
# a volume of dichlormethane after incubating the septum in it was taken as gc-sample. The
# high number of peaks (126) makes it unlikely to be a proper control, hence exclude from here onwards.
# charadrius_peaks[["w41"]] <- NULL # remove w41



#check_input(data = charadrius_peaks,show_peaks = T) # check format, check peak distribution
# charadrius_peaks_aligned <- align_chromatograms(data =charadrius_peaks ,
#                                   rt_col_name = "rt", # retention time
#                                   conc_col_name = "area", # peak abundance
#                                   reference = "w62", # 170 peak, highest count in the sample
#                                   write_output = c("rt","area"),
#                                   blanks = c("w17","w37","w47","w57","w67","w77"),
#                                   delete_single_peak = T,
#                                   min_diff_peak2peak = 0.03,
#                                   max_diff_peak2mean = 0.03,
#                                   rt_cutoff_low = 8 # peaks before the solvent are treated as uncertain
#                                   )
# save(charadrius_peaks_aligned,file = "data/charadrius_peaks_aligned.RData")

# load data and get an quick idea of the data
# -------------------------------------------
load("data/charadrius_peaks_aligned.RData") # GCalignR output
class(charadrius_peaks_aligned) # object of class GCalign
names(charadrius_peaks_aligned) # includes three lists
plot(x = charadrius_peaks_aligned) # summarizes linear adjustments and output variation 
summary(charadrius_peaks_aligned) # summary of the alignment procedure

# get the area, normalise & log-tranform 
# --------------------------------------
scent <- norm_peaks(charadrius_peaks_aligned,conc_col_name = "area",rt_col_name = "rt",out = "data.frame") # normalise are
scent <- log(scent + 1)
head(scent[1:6]) # Scent Matrix, Individuals in rows, peaks in columns

# get the corresponding factors
# -----------------------------
factors <- read.csv("data/factors.csv",sep = ";",skip = 1) # Comprises all available  covariates for every sample, excl. blanks
row.names(factors) <- as.character(tolower(factors$GC_Sample)) # for reference to the scent data
factors <- factors[,c("Species","Ring","Nest_Brood","Brood_Status","ID","Sex","Age","GC_Sample")]
head(factors) 


# Nonmetric Multidimensional Scaling (NMDS) Ordination using vegan::
# http://cran.r-project.org/web/packages/vegan/vegan.pdf
# ------------------------------------------------------------------

scent_nmds <- nmds_calculator(scent = scent,factors = factors) 
p1 <- nmds_plotter(nmds = scent_nmds,main = "All Individuals") 

scent_nmds <- nmds_calculator(scent = scent,factors = factors,sub = list(Age="A")) 
p2 <- nmds_plotter(nmds = scent_nmds,main = "Adults") 
p2

scent_nmds <- nmds_calculator(scent = scent,factors = factors,sub = list(Age="A",Brood_Status=c("Brood","Nest"))) 
p3 <- nmds_plotter(nmds = scent_nmds,main = "Adults") 
p3


#?adonis #vegan::adonis
# Permutational Multivariate Analysis of Variance using Distance Matrices
# Thought to less sensitive to dispersion effects than anosim
# More robust than permanova, i.e. can take continous variables
# betadisper can be used to study disperion within the same framework

# ADONIS
#########
stats <- adonis(scent~factors$Species*factors$Brood_Status, permutations=999,method="bray")
                #strata=factors$age # groups within which to constrain permuations: For nested design

stats
R2 <- paste("R² =",round(stats$aov.tab$R2[1],2)) #label for plots
Pval <- paste("p =",round(stats$aov.tab$`Pr(>F)`[1],3)) # label for plots

# BETADISPER
############
anova(betadisper(vegdist(scent,method = "bray"), factors$Species)) # test for dispersion differences



