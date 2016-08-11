nmds_calculator <- function(scent,factors,sub=NULL){
# Using the package vegan, this function allows to do a non metric multidimensional scaling using bray-curtis dissimilarity matrices. The maximum number of random starts is set to 100, which is more the enough to find a optimal configuration, which is usually done within less than 10 runs. Using the argument "sub" it is possible to subset factors based on their levels
# scent is a data frame of log-transformed normalised peaks
# factors a data frame containing variables, 
# rownames of those df´s are expected to be corresponding to each other, i.e. allow to link vars to the values
# sub is a list of factor names and their arguments, they must correspond to the factor df

  # subsetting, if sub is specified
  if(!is.null(sub)){
    index <- which(names(factors) %in% names(sub))[match(names(factors),names(sub))[!is.na(match(names(factors),names(sub)))]]
    vars <- (sub)
    
    for(x in 1:length(sub)){
      factors <- subset(factors,factors[,index[x]] %in% vars[[x]]) # extract factors by sub specifications
      factors
    }
    
  }
  
  scent <- scent[match(row.names(factors),row.names(scent)),] # same order is crucial 
  scent_nmds <- vegan::metaMDS(scent,
                               distance = "bray", # Bray-Curtis dissimilarity
                               k=2, # Number of dimension, where n > 2*k+1
                               trymax=100, # Max. number of random starts
                               autotransform=FALSE, # Use simple heuristics for possible data transformation of typical community data if TRUE
                               expand=FALSE, 
                               plot=F)
  nmds <- as.data.frame(scent_nmds$points)
  nmds <- cbind(nmds,factors)
  
  
  return(list(nmds = nmds,factors = factors,scent = scent))
}

