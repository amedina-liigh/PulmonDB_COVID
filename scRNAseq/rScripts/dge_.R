### DIFFERENTIAL EXPRESSIÓN ANALYSIS
### Seurat
library(Seurat)
library(dplyr)

# Read data and add information
nCoV <- readRDS("/home/amedina/amedinalab/larteaga/COVID19/R/nCoV.rds")
metadata <- read.delim2("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/all.cell.annotation.meta.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")

nCoV <- nCoV[,nCoV$ID %in% metadata$ID]
nCoV$group2 <- metadata$group
nCoV$group2 <- as.factor(nCoV$group2)
nCoV$celltype <- metadata$celltype
nCoV$ctc <- paste0(metadata$celltype,"_",metadata$disease)
nCoV$cts <- paste0(metadata$celltype,"_",metadata$group)

levels(nCoV)
nCoV2 <- nCoV
Idents(nCoV2) <- nCoV2$cts

diffexpgenes <- list()
celltypes = sort(unique(nCoV2$celltype))[sort(unique(nCoV2$celltype)) != 'Neutrophil']
for (i in celltypes) {
  print(i)
  diffexpgenes <- c(diffexpgenes,list(FindMarkers(nCoV2, ident.1 = paste0(i,"_Y"), ident.2 = paste0(i,"_N"),test.use = "MAST")))
}
names(diffexpgenes) <- celltypes

saveRDS(diffexpgenes, "/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/diffexpgenesMAST.rds")
diffexpgenes <- readRDS("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/reticulate/diffexpgenesMAST.rds")



