###################################################################################
# Analysis of preen wax compositions of KIP, MP & WFP during parental care stages #
###################################################################################

rm(list = ls())
###############################################
# A. Load Data, source Packages & run GCalignR
# #############################################
library(GCalignR)
library(vegan)
library(ggplot2)
source("R/nmds_plotter.R") # little for nmds plots using ggplot
source("R/nmds_calculator.R") # does nmds scaling and makes is ready for plotting

################################################################################################
# Loading the raw data, gas-chromatography peaks called with Xcalibur (Thermo Fisher Scientific)
################################################################################################
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
#factors <- factors[-c("w43","w53","w64")]
scent <- scent[match(row.names(factors),row.names(scent)),] # same order is crucial 


# Nonmetric Multidimensional Scaling (NMDS) Ordination using vegan::
# http://cran.r-project.org/web/packages/vegan/vegan.pdf
# ------------------------------------------------------------------

# With all samples
# ----------------
m1 <- nmds_calculator(scent = scent,factors = factors)
nmds_m1 <- m2$nmds
scent_m1 <- m1$scent
factors_m1 <- m1$factors
p1 <- nmds_plotter(nmds = nmds_m1,main = "All Individuals") 
anova(vegan::betadisper(vegdist(scent_m1,method = "bray"), factors_m1$Species)) # May be okay
adonis(scent_m1~factors_m1$Species*factors_m1$Brood_Status, permutations=999,method="bray")
#strata=factors$age # groups within which to constrain permuations: For nested design

# Only adults
# ----------------
#factors$family <- paste0(factors$Nest_Brood,factors$Species)
m2 <- nmds_calculator(scent = scent,factors = factors,sub = list(Age="A",Brood_Status="Nest",Species="KIP"))
nmds_m2 <- m2$nmds
scent_m2 <- m2$scent
factors_m2 <- m2$factors
p2 <- nmds_plotter(nmds = nmds_m2,main = "Adults") 
beta_2 <- vegan::betadisper(vegdist(scent_m2,method = "bray"), factors_m2$Species) # Dispersion effect
anova(beta_2)
TukeyHSD(beta_2)# KIP & WFP differ, meaining KIP´s are more similar to conspecifics
adonis(scent_m2~factors_m2$Species*factors_m2$Brood_Status,strata = factors_m2$Nest_Brood, permutations=999,method="bray")
#strata=factors$age # groups within which to constrain permuations: For nested design

# SIMPER Analysis

# getting 15 best substances and their contribution to Species dissimilarity
simp_Species <- vegan::simper(scent_m2, factors_m2$Species)

simp_Species_names <- rownames(summary(simp_Species, ordered = TRUE)[[1]])[1:15]
contribution <- summary(simp_Species, ordered = TRUE)[[1]]$contr[1:15]
ind_col <- paste(which(names(scent_m2)%in%simp_Species_names), collapse = ",")
# connect to data frame and compute contribution in percent
col_simp <- data.frame(comp = simp_Species_names, contrib = contribution*100, stringsAsfactors_m2 = FALSE)
col_simp
Species dissimilarity based on 15 compounds.
# overall (number of permutations is 1000 instead of 10,000 in the paper)
anosim(dat = scent_m2[col_simp$comp], grouping = factors_m2$Species,
       distance = "bray", permutations = 1000)


#?adonis #vegan::adonis
# Permutational Multivariate Analysis of Variance using Distance Matrices
# Thought to less sensitive to dispersion effects than anosim
# More robust than permanova, i.e. can take continous variables
# betadisper can be used to study disperion within the same framework
# One important assumption, but often violoated, is that beta-dispersion among groups should be the
# same! (see http://thebiobucket.blogspot.de/2011/04/assumptions-for-permanova-with-adonis.html)
# 

# Lets check the beta-dispersion
anova(vegan::betadisper(vegdist(scent,method = "bray"), factors$Species)) # May be okay

# ADONIS Permutational Multivariate Analysis of Variance Using Distance Matrices
# Use strata for nested design, i.e. blocking factor
model <- adonis(scent~factors$Species*factors$Brood_Status, permutations=999,method="bray")
                #strata=factors$age # groups within which to constrain permuations: For nested design
print(model)


# stats 
# R2 <- paste("R² =",round(stats$aov.tab$R2[1],2)) #label for plots
# Pval <- paste("p =",round(stats$aov.tab$`Pr(>F)`[1],3)) # label for plots

# Posthoc tests for adonis ....
#  If anybody figures out how to do post hoc test for adonis() we may incorporate submitted code into vegan (though personally I think that the whole idea of post hoc tests is weird, to put it mildly). J. Oksane, vegan creator

