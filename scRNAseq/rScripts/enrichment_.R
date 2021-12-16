### FUNCTIONAL ENRICHMNET ANALYSIS
library(gprofiler2)
library(GSEABase)
library(dplyr)
library(tidyverse)

# Read data
TF_targets <- readRDS("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/TF_targets.rds")

genesets_h <- read.gmt("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/resources/h.all.v7.4.symbols.gmt")

diffexpgenes <- readRDS("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/diffexpgenesMAST.rds")

# GSEA on DGE
for (i in 1:9) {
  gsea_ranks <- diffexpgenes[[i]]$avg_log2FC
  names(gsea_ranks) <- as.character(rownames(diffexpgenes[[i]]))
  gsea_ranks <- sort(gsea_ranks, decreasing = T)
  
  gsea_res <- GSEA(gsea_ranks, TERM2GENE = genesets_h, verbose = F)
  
  gsea_res_df <- as_tibble(gsea_res@result)
  gsea_res_df <- gsea_res_df %>%
    mutate(Disease = case_when(
      NES > 0 ~ "Y",
      NES < 0 ~ "N"))
  gsea_res_df <- gsea_res_df[order(gsea_res_df$NES, decreasing = T),]
  write.table(gsea_res_df, paste0("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/GSEA/gsea_",names(diffexpgenes)[i],".tsv"),
              quote = F, row.names = F, sep = "\t")
}

# GO analysis on regulons by cell type and condition
go_paths <- list()
for (i in 1:9){
  go_paths_ct <- list()
  print(names(TF_targets)[i])
  for (j in 1:length(TF_targets[[i]])) {
    print(names(TF_targets[[i]][j]))
    genes <- TF_targets[[i]][j]
    gost_res <- gost(genes, organis = "hsapiens", correction_method = "fdr")
    gost_res_sub <- gost_res$result[c('term_name','term_id','p_value','source')][grep("GO",gost_res$result[c('term_name','term_id','p_value','source')]$term_id),]
    gost_res_sub <- gost_res_sub[order(gost_res_sub$p_value),]
    go_paths_ct <- c(go_paths_ct,list(gost_res_sub[grep("BP",gost_res_sub$source),][,1]))
  }
  names(go_paths_ct) <- names(TF_targets[[i]])
  go_paths <- c(go_paths,list(go_paths_ct))
}
names(go_paths) <- names(TF_targets)

saveRDS(go_paths,"/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/go_paths.rds")


