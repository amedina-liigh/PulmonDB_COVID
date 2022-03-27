## Functions for Regulons Selection
## deciding FDR threshold

### make a table of number of regulons retrievng when using a range of FDR adj pvalue thrsholds
regsPerFDRthr <- function(pvalues.df, fdr_lowth = 0.1, fdr_upthr = 0.2){
  # list with number of regulons per experiment when changing the FDR threshold 
  fdrs <- list()
  i <- 1
  for (fdr_thr in seq(fdr_lowth,fdr_upthr, by=0.01)){
    fdrs[[i]] <- apply(pvalues.df, 2, function(x) length(which(x < fdr_thr)))
    i=i+1
  }
  names(fdrs) <- paste0("FDR_",seq(fdr_lowth,fdr_upthr, by=0.01))
  
  # turn to dataframe of thresholds
  thrs <- t(data.frame(matrix(unlist(fdrs), nrow=length(fdrs), byrow=TRUE)))
  colnames(thrs) <- paste0("FDR_",seq(fdr_lowth,fdr_upthr, by=0.01))
  rownames(thrs) <- colnames(pvalues.df)
  thrs <- as.data.frame(thrs) %>% arrange_all(desc) # arrange rows in descending order
  
  return(thrs)
}

### get the FDRs threshold to get a min number of regulons
###   limit.range.at : if at a given FDR no regulons are retrieve, just take whichever FDR that first retrieve at least 1 regulon
###                     specify col in thresholds table 
getChosenFDRth <- function(FDRthresholds, min.regs = c(100,1), limit.range.at = nrow(FDRthresholds)){
  
  # at which minimal fdr are at least n regulons
  chosen <- apply(FDRthresholds,1,function(x) which(as.double(x[1:limit.range.at]) >= min.regs[1])[1])
  # for those that not the min was found, change min regs to 1
  lowns <- which(is.na(chosen))
  chosen[lowns] <- apply(FDRthresholds[lowns,],1,function(x) which(as.double(x) >= min.regs[2])[1])
  # if not the min was found, take the minimal val founnd
  lowns <- which(is.na(chosen))
  chosen[lowns] <- apply(FDRthresholds[lowns,],1,function(x) which(as.double(x) >= min(FDRthresholds))[1])
  chosen <- as.data.frame(chosen)
  
  # get corresponding threshold
  chosen['FDR_thr'] <- colnames(FDRthresholds)[chosen$chosen]
  # add number of regulons (also equal to the number of accepted pvalues adj)
  nregs <- c()
  for (i in 1:nrow(chosen)) {
    nregs[i] <- FDRthresholds[i,chosen$chosen[i]]
  }
  chosen['nregulons'] <- nregs
  # adjust format
  chosen <- chosen %>% mutate(FDR_thr = as.double(str_remove(FDR_thr,"FDR_"))) %>% arrange(FDR_thr)
  
  return(chosen)
}
