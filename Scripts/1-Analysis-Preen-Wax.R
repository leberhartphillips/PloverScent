#############################################
# Analysis of Plover preen wax compositions #
#############################################
rm(list = ls())

# A. Load Data, source Packages & run GCalignR
# --------------------------------------------
library(GCalignR)
# check_input(data = "data/CharadriusPreen.txt",plot=T) # check format, check peak distribution
# CharadriusPreen <- align_chromatograms(data = "data/CharadriusPreen.txt",
#                                   rt_col_name = "rt", # retention time
#                                   conc_col_name = "area", # peak abundance
#                                   reference = "w62", # 170 peak, highest count in the sample
#                                   write_output = c("rt","area"),
#                                   blanks = c("w17","w37","w47","w57","w67","w77"),
#                                   delete_single_peak = T,
#                                   min_diff_peak2peak = 0.03,
#                                   max_diff_peak2mean = 0.02,
#                                   rt_cutoff_low = 8 # peaks before the solvent are treated as uncertain
#                                   )
# save(CharadriusPreen,file = "data/CharadriusPreen.RData")

# load data and get an quick idea of the data
# -------------------------------------------
load("data/CharadriusPreen.RData") # GCalignR output
# source("R/ChromaSimFunctions.R")
class(CharadriusPreen) # object of class GCalign
names(CharadriusPreen) # includes three lists
plot(x = CharadriusPreen) # summarizes linear adjustments and output variation 
summary(CharadriusPreen) # summary of the alignment procedure

# get the area, normalise & log-tranform 
# --------------------------------------
scent <- norm_peaks(CharadriusPreen,conc_col_name = "area",rt_col_name = "rt",out = "data.frame") # normalise are
scent <- log(scent + 1)
head(scent[1:6]) # Scent Matrix, Individuals in rows, peaks in columns

# get the corresponding factors
# -----------------------------
factors <- read.csv("data/factors.csv",sep = ";",skip = 1) # Comprises all available  covariates for every sample, excl. blanks
#factors <- factors[,c("Species","Ring","Sex","Age","Nest_Brood","Brood_Status","GC_Sample")] # subsetting
row.names(factors) <- as.character(tolower(factors$GC_Sample)) # for reference to the scent data
head(factors) 

scent_mds <- MASS::isoMDS(vegan:: vegdist(scent)) # Bray-curtis similarity matrix
vegan::ordiplot(scent_mds, type = "t", ylab = "", xlab = "",axes=FALSE, frame.plot=TRUE) #NMDS plots


scent <- scent[match(row.names(factors),row.names(scent)),] # same order of rows is convenient



# For Cara in Primer 6
factors <- subset(factors,factors$Age=="A" & factors$Nest_Brood!="Unknown") # only adults with know brood
scent <- scent[match(row.names(factors),row.names(scent)),]

write.csv2(scent,file = "ScentPrimer6.csv")
write.csv2(factors,file = "FactorsPrimer6.csv")


# Anosim
# -----------------
vegan::anosim(dat = scent, grouping = factors$Species, distance = "bray", permutations = 1000)

