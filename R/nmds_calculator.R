nmds_calculator <- function(scent,factors,sub=NULL){
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
  
  # subsetting 
  if(!is.null(sub)){
    index <- which(names(nmds)==names(sub))[match(names(nmds),names(sub))[!is.na(match(names(nmds),names(sub)))]]
    vars <- (sub)
    
    for(x in 1:length(sub)){
    nmds <- subset(nmds,nmds[,index[x]] %in% vars[[x]])
    nmds
    }
    
  }
  
  return(nmds)
}

