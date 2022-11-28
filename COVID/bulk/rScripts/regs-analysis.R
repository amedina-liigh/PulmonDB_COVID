## regs-analysis.R
## mpadilla
## 2022
## Functions for regulons analysis - COVID-19 Project
## used in regulons-selection and regulons-exploration reports 

## Load libraries
library(tidyverse) # dplyr, stringr, ggplot2
library(clusterProfiler) # bitr, enrichGO
library(org.Hs.eg.db) # bitr
library(enrichplot) # dotplot

## volcano_per_exp
## to restore, used in regulons-selection reports

## distinctlySharedRegs
## Get distinctly shared sets from a list as a table
distinctlySharedRegs <- function(tfs_per_exp_list){
  
  n = length(names(tfs_per_exp_list))
  
  ## Get shared elements (regulons) between each pair of conditions (each set)
  sharedRegs <- list()
  i = 1
  j = 2
  universeShared <- c()
  for (exp_i in names(tfs_per_exp_list)[i:n-1]) {
    tmp_exp_list <- list()
    for (exp_j in names(tfs_per_exp_list)[j:n]) {
      tmp_exp_list[[exp_j]] <- intersect(tfs_per_exp_list[[exp_i]], tfs_per_exp_list[[exp_j]])
      tmp_exp_list[[exp_j]] <- str_remove(tmp_exp_list[[exp_j]], "_\\(\\+\\)")
      ## make universe of shared elements
      universeShared <- c(universeShared,tmp_exp_list[[exp_j]])
    }
    sharedRegs[[exp_i]] <- tmp_exp_list
    i = i+1
    j = j+1
  }
  
  ## Initiate table
  table <- as.data.frame(matrix(nrow = n, ncol = n))
  rownames(table) <- names(tfs_per_exp_list)
  colnames(table) <- names(tfs_per_exp_list)
  ## list of sets that were only shared once
  firstSeen <- match(unique(universeShared),universeShared)
  uniquelyShared <- setdiff(universeShared[firstSeen], universeShared[-firstSeen])
  
  ## refine list of shared elements to keep only 
  ##   elements shared in a given pair of experiments
  i = 1
  j = 2
  for (exp_i in names(tfs_per_exp_list)[i:n-1]) {
    for (exp_j in names(tfs_per_exp_list)[j:n]) {
      ## concatenate elements to add to table
      table[exp_i, exp_j] <- paste0(intersect(sharedRegs[[exp_i]][[exp_j]], uniquelyShared), collapse = ",")
    }
    # also report elements only seen in a given condition
    table[exp_i, exp_i] <- paste0(setdiff(str_remove(tfs_per_exp_list[[exp_i]], "_\\(\\+\\)"),
                                       unique(universeShared)), collapse = ",")
    i = i+1
    j = j+1
  }
  ## add last condition unique set
  table[n, n] <- paste0(setdiff(str_remove(unname(unlist(tfs_per_exp_list[n])), "_\\(\\+\\)"),
                                        unique(universeShared)), collapse = ",")
  
  return(table)
} # EOF

## getDistinct
## Get the complement from two groups in a list
getDistinct <- function(tfs_per_exp_list, indexInterestGroup){
  
  ## Make list of groups (group1 = interest group)
  group1 <- c()
  group2 <- c()
  i = 1
  for (exp in names(tfs_per_exp_list)) {
    if(i %in% indexInterestGroup){
      group1 <- c(unlist(tfs_per_exp_list[[exp]]))
    } else{
      group2 <- c(unlist(tfs_per_exp_list[[exp]]))
    }
    i = i+1
  }
  
  ## Return complement
  complement <- setdiff(group1, group2)
  complement <- sort(unique(str_remove(complement, "_\\(\\+\\)")))
  return(complement)
}

## selectTFs
## Get the elements from a list
selectTFs <- function(tfs_per_exp_list){
  all_tfs <- c()
  for (exp in names(tfs_per_exp_list)){
    all_tfs <- c(all_tfs, str_remove(unlist(tfs_per_exp_list[[exp]]), "_\\(\\+\\)"))
  }
  return(unique(all_tfs))
}

## family_annot
## Get a df with Family and DBD annotation for a set of TFs
## Requires file params$TF_Family_DBD_hsapiens_db
family_annot <- function(tfs_interest, tfs_db_csv = params$TF_Family_DBD_hsapiens_db, verbosity = 1){
  # Read DB
  db <- read_csv(file = tfs_db_csv, col_names = T, quote = "")
  
  # Search TFs of interest in db
  found <- tfs_interest[tfs_interest %in% db$TF_Name]
  notfound <- tfs_interest[tfs_interest %in% db$TF_Name == FALSE]
  if(verbosity > 1){print(paste0("There were ",length(found)," TFs found in the TF annotation database"))}
  if(length(notfound) > 0){print(paste0("Warning! ",length(notfound), " TFs were not found in the database. Annotating them as UNKNOWN"))}

  # Extract info from db
  found_db_annot <- db[db$TF_Name %in% tfs_interest,]
  ## add missing tfs
  notfound_db_annot <- data.frame(TF_Name = notfound,
                               Family_Name = rep("Unknown",length(notfound)),
                               DBDs = rep("UNKNOWN", length(notfound)))
  tfs_annot <- dplyr::full_join(found_db_annot, notfound_db_annot)
  
  return(tfs_annot)
} # EOF

## regsGenesPerExp
## Retrieve all regulon target genes from the list of regulons per experiment
## requires list of gene targets per regulon, return genes as entrez ids if TRUE
regsGenesPerExp <- function(tfs_per_exp_list, tfs_targets, entrezids = F){
  
  ## Make list per experiment
  regs_genes_per_exp <- list()
  for (exp in names(tfs_per_exp_list)){
    reg_genes <- c()
    ## search targets for each regulon
    for (regulon in tfs_per_exp_list[[exp]]) {
                                      ## TF                             ## target genes
      reg_genes <- c(reg_genes, str_remove(regulon, "_\\(\\+\\)"), tfs_targets[[regulon]])
    }
    ## array of all genes in all regulons in a given experiment
    reg_genes <- unique(reg_genes)
    ## return entrezids if specified
    if(entrezids == TRUE){
      reg_genes <- bitr(reg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    }
    regs_genes_per_exp[[exp]] <- reg_genes
  }
  
  return(regs_genes_per_exp)
}

## iterDotplot
## Generate Enriched Terms Dotplots from ClusterProfiler
iterDotplot <- function(ego_perexp_da, db = "GO:BP", savepngs = FALSE, outdir = params$outdir){
  
  if(savepngs==TRUE){db.save = str_remove(db, ":")}
  
  ## Generate dotplots per ego enrichement analysis
  plots <- list()
  for (exp in names(ego_perexp_da)){
    ego <- ego_perexp_da[[exp]]
    ## check that enriched terms were found
    if(identical(ego[,1], character(0))){
      next
    }
    plots[[exp]] <- enrichplot::dotplot(ego, showCategory=20, label_format = 50, font.size=7) + 
      ggtitle(paste0(db," - ", exp))
    ## Save png if specified
    if(savepngs == TRUE){
      ggsave(paste0(outdir, "dotplot-",db.save,"_",exp,"_regulons.png"))
    }
  }
  
  return(plots)
} # EOF
