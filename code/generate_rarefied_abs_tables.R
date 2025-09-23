require(abind)
require(extraDistr)

rarefy_sample <- function(otu_vector, iterations, depth, abs.count, return_rel=FALSE) { 
  
  otu_counts <- rmvhyper(iterations, otu_vector, depth)
  
  if(return_rel){
    rare_abs_counts <- otu_counts/depth
  }else{
  rare_abs_counts <- round((otu_counts / depth) * abs.count)
  }

  return(rare_abs_counts)
  }


generated_rarefied_abs_tables <- function(physeq, iterations, rare.depth, abs.counts, seed = 1, return_rel = FALSE){
  
  if(taxa_are_rows(physeq)){
    
    df <-
      physeq %>%
      otu_table() %>%
      data.frame()
    
  }else{
    
    df <-
      physeq %>%
      otu_table() %>%
      t() %>%
      data.frame()
  }
  
  if(any(colnames(df) != names(abs.counts))){
    stop("Colnames of OTU table don't match IDs associated with cell counts; normalization could be inaccurate", call. = FALSE)
  }
  
  set.seed(seed)
  rare_tables <- pmap(list(df, abs.counts), \(x,y)rarefy_sample(otu_vector = x, iterations = iterations, depth = rare.depth, 
                                                                abs.count = y, return_rel = return_rel))
  
  bound_together <- abind(rare_tables, along = 0)
  
  dimnames(bound_together)[[2]] <- paste0("iter_",c(1:iterations))
  
  dimnames(bound_together)[[3]] <- row.names(df)
  
  return(bound_together)
}
