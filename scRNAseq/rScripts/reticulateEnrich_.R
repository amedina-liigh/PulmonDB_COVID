## Load libraries
library(reticulate) # use python within R
library(dplyr)
library(tidyverse)
library(UpSetR) # graphics
library(ComplexHeatmap)
### SET CONDA ENV TO USE WITH RETICULATE
use_condaenv(condaenv = "scanpy", conda = "/cm/shared/apps/anaconda3/2021.05/envs/scanpy", 
             required = TRUE)
py_config()
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/signif-regulons/FDRthresholds_functions.R")

############
# Run script "difregs_.py" to use results in R
## Pass variables with input files

## Leo's single cell data
py_run_string('path_scenic_loom = "/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/SCENIC_final_multi.loom"')
# to run with changed difregs_.py script, it would require some changes: (cols are sample|cell, source|cell_type and disease_YorN)
py_run_string('path_meta = "/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/all.cell.annotation.meta.txt"')
py_run_file("/mnt/Citosina/amedina/mpadilla/COVID19/PulmonDB_COVID/scRNAseq/pyScripts/difregs_.py")

## 1. celllines experiment
##  scenic multiruns output loom
py_run_string('path_scenic_loom = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/celllines/out/scenic/celllines100/SCENIC_SCope_output.loom"')
##  metadata csv (cols are sample|cell, source|cell_type and disease_YorN)
py_run_string('path_meta = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease-cellline_subCelllines.csv"')
##  flag variable: filter auc mtx, to filter samples not included in metadata file
py_run_string('filter_auc = 1')
##  flag variable: mock samples are controls for several cellline experiments
py_run_string('celllines = 1')
##  run script to compute encircled and differentialy activated regulons
py_run_file("/mnt/Citosina/amedina/mpadilla/COVID19/PulmonDB_COVID/scRNAseq/pyScripts/difregs_.py")
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease-cellline_subCelllines.csv", sep = ","))
celllines <- 1

## 2. lung bladesdel experiment
py_run_string('path_scenic_loom = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/lung-bladesdel/out/scenic/lung-bdd100/SCENIC_SCope_output.loom"')
py_run_string('path_meta = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/blanco2020-desai2020-delorey2021_allLungSamples_ID-source-disease.csv"')
py_run_file("/mnt/Citosina/amedina/mpadilla/COVID19/PulmonDB_COVID/scRNAseq/pyScripts/difregs_.py")
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/blanco2020-desai2020-delorey2021_allLungSamples_ID-source-disease.csv", sep = ","))

## 3. multiorgans desai-gtex experiment
py_run_string('path_scenic_loom = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/organs-desgtex/out/scenic/org_desgtex100/SCENIC_SCope_output.loom"')
py_run_string('path_meta = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/desai2020_gtexv8-subTissues-sample_source_disease.csv"')
py_run_file("/mnt/Citosina/amedina/mpadilla/COVID19/PulmonDB_COVID/scRNAseq/pyScripts/difregs_.py")

############
# Save variables from python script to R
keys_ct = names(py$pvalues_adjusted) # cases names
n_cases = n_distinct(keys_ct) 
rss_both_ctc = r_to_py(py$rss_group) # rsscores for each group of source|celltype and disease state Y|N

# define cases and controls names
if(celllines == 1){
  cases.lines <- as.data.frame(meta %>% filter(disease == "Y") %>% select(source,cellline) %>% unique())
  cases.lines <- cases.lines[match(keys_ct,cases$source),]
}
controls <- meta %>% filter(disease == "N") %>% select(source) %>% unique()
controls <- sort(as.array(controls[,1]), decreasing=T)

# get TFs (leading regulons) names
tfs = str_subset(names(py$auc_control), "_\\(\\+\\)")

############
# Get common DA from both tests using the pvals from the 1stQT
# Take regulons significaly differentialy activated from both tests in all cases, one test is suffice though
tfs_mwu_ks = list()
for (i in 1:n_cases) {
  # assess which passed the mwu or ks tests
  seto1 = tfs[unname(unlist(py$regsPass_1stQ[[i]]))]
  seto2 = tfs[unname(unlist(py$regsPass_1stQ2[[i]]))]
  # get tfs that passed both tests
  tfs_mwu_ks = c(tfs_mwu_ks,list(intersect(seto1,seto2)))
}
names(tfs_mwu_ks) <- names(py$regsPass_1stQ)

############
### Filter adjusted pvalues for a given FDR threshold
# retrieve adjusted pvalues from difregs script
pvalues_adj <- py$pvalues_adjusted
pvalues_adj2 <- py$pvalues_adjusted2

# turn to dataframes
pv.adj.mwu <- t(data.frame(matrix(unlist(pvalues_adj), nrow=length(pvalues_adj), byrow=TRUE)))
colnames(pv.adj.mwu) <- keys_ct
pv.adj.klm <- t(data.frame(matrix(unlist(pvalues_adj2), nrow=length(pvalues_adj2), byrow=TRUE)))
colnames(pv.adj.klm) <- keys_ct
if(celllines==1){
  colnames(pv.adj.mwu) <- str_replace_all(colnames(pv.adj.mwu), "-", "")
  colnames(pv.adj.mwu) <- str_replace(colnames(pv.adj.mwu), "w_vector_expressing_", "")
  colnames(pv.adj.mwu) <- str_replace(colnames(pv.adj.mwu), "1hr_Ruxolitinib_pretreatment", "pt")
  colnames(pv.adj.klm) <- colnames(pv.adj.mwu)
}

# Find FDR thresholds (within a 0.1-0.2) per experiment that retrieve ~100 regulons 
thrs.mwu <- regsPerFDRthr(pv.adj.mwu)
fdrs.mwu <- getChosenFDRth(thrs.mwu)
fdrs.mwu
# for celllines case in klm, fdr thr extended to 0.3
thrs.klm <- regsPerFDRthr(pv.adj.klm, fdr_upthr = 0.3)
fdrs.klm <- getChosenFDRth(thrs.klm, min.regs = 82, limit.range.at = 11)
fdrs.klm

# Apply FDR filter to adjusted pvalues
fdrs.mwu <- fdrs.mwu[match(colnames(pv.adj.mwu),rownames(fdrs.mwu)),] # arrange rows as pvalues df
regsPassFDR.mwu <- as.data.frame(pv.adj.mwu <= fdrs.mwu$FDR_thr)
fdrs.klm <- fdrs.klm[match(colnames(pv.adj.klm),rownames(fdrs.klm)),] # arrange rows as pvalues df
regsPassFDR.klm <- as.data.frame(pv.adj.klm <= fdrs.klm$FDR_thr)

# apply less strict filter
regsPassFDR.mwu <- as.data.frame(pv.adj.mwu <= 0.2)
regsPassFDR.klm <- as.data.frame(pv.adj.klm <= 0.25)

############
# Get common DA from both tests using the pvals from the FDR cutoffs
tfs_mwu_ks_fdr = list()
for (i in 1:n_cases) {
  # assess which passed the mwu or ks tests
  seto1 = tfs[regsPassFDR.mwu[,i]]
  seto2 = tfs[regsPassFDR.klm[,i]]
  # get tfs that passed both tests
  tfs_mwu_ks_fdr = c(tfs_mwu_ks_fdr,list(intersect(seto1,seto2)))
}
names(tfs_mwu_ks_fdr) <- names(regsPassFDR.mwu)

## check how many regulons were kept
t(as.data.frame(lapply(tfs_mwu_ks_fdr,length)))

## 
firstqt <- 0 # in my case im using the fdr threshold
if(firstqt==0){
  tfs_mwu_ks <- tfs_mwu_ks_fdr 
}

############
# df of cases and controls names in RSS matrix

namesInRSS <- data.frame(
  cases = character(),
  controls = character()
)
for (i in 1:n_cases) {
  # this is how cases and controls are written in the rss table
  caseY = paste0(keys_ct[i],'Y')
  if(celllines == 1){
    if(cases.lines$cellline[i] == "A549_w_vector_expressing_hACE2_1hr_Ruxolitinib_pre-treatment"){
      control = "Mock_treated_A549"
    } else {
      control = str_subset(controls, paste0(cases.lines$cellline[i],"$")) # get corresponding control
    }
    controlN = paste0(control,'N')
  } else {
    controlN = paste0(keys_ct[i],'N')
  }
  namesInRSS <- add_row(namesInRSS, cases= caseY, controls=controlN)
}

############
# AQUÍ YA PODRÍAS TENER LOS REGULONES QUE ESTÁN DIFERENCIALMENTE ACTIVOS (DA)
# PERO LO QUE HICE AQUÍ FUE DE LOS QUE SALIERON DA, SELECCIÓN LOS QUE ESTABAN MÁS ACTIVOS EN CADA TIPO CELULAR
# SEGÚN EL RSS
ccmk = list()
for (i in 1:n_cases) {
  # get ordered regulon specificity scores for case
  top_rss_cov_Y = py_to_r(rss_both_ctc$T$nlargest(length(tfs),namesInRSS$cases[i])[namesInRSS$cases[i]])
  # get regulon specificity scores for its control
  top_rss_cov_N = py_to_r(rss_both_ctc$T$nlargest(length(tfs),namesInRSS$cases[i])[namesInRSS$controls[i]])
  # filter for scores > 1st quartile
  temp <- top_rss_cov_Y > summary(top_rss_cov_Y)[2]
  top_rss_cov_Y = top_rss_cov_Y[temp]
  top_rss_cov_N = top_rss_cov_N[temp]
  top_rss_cov_Y = top_rss_cov_Y[top_rss_cov_Y > top_rss_cov_N] # only take regulons more specific to cases
  # take regulons that are both enriched and differentially active
  ccmk = c(ccmk,list(intersect(names(top_rss_cov_Y),tfs_mwu_ks[[i]])))
}
names(ccmk) <- names(tfs_mwu_ks)

## check how many regulons were kept
t(as.data.frame(lapply(ccmk,length)))

############
filter.by.lfc <- 1 # filter by an optimal high LFC and min regulons
min.regulons <- 10 # min regulons to retrieve 
lfc.scores <- 0 # store LFC used

# TAMBIÉN PODRÍAS HACER UN PSEUDO-FOLD CHANGE PA VER CUALES SALEN
lfc <- list()
firstqt <- 0 # first QT cutoff is 0 and FDR adj pval cutoff is 1
min.lfcn <- 0.01 # LFC min val

for (i in 1:n_cases) {
  if(firstqt == 1){
    # regulons that passed 1st qu
    p <- tfs[unlist(py$regsPass_1stQ[keys_ct[i]])]
  } else {
    # regulons that passed fdr thresholds
    p <- tfs_mwu_ks_fdr[[i]]
  }
  lfc_Y = py_to_r(rss_both_ctc$T$nlargest(length(tfs),namesInRSS$cases[i])[namesInRSS$cases[i]])
  lfc_N = py_to_r(rss_both_ctc$T$nlargest(length(tfs),namesInRSS$cases[i])[namesInRSS$controls[i]])
  q <- log2(lfc_Y/lfc_N) # lfc foreach regulon
  
  if(leo == 1){
    if (i == 1 | i == 8 | i == 9) {
      lfc <- c(lfc,list(names(q[q > 0.05 & names(q) %in% p]))) 
    } else if (i == 4){
      lfc <- c(lfc,list(names(q[q > 0.01 & names(q) %in% p]))) 
    } else {
      lfc <- c(lfc,list(names(q[q > 0.1 & names(q) %in% p])))
    }
  }
  
  ## use a LFC that allows a min number of regulons
   if(filter.by.lfc == 1){
    lfcn <- 1.5
    if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
      lfcn <- 1
      if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
        lfcn <- 0.5
        if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
          lfcn <- 0.1
          if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
            lfcn <- 0.01
            if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
              lfcn <- 0.0005
              if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
                lfcn <- 0.00005

              }
            }
          }
        }
      }
    }
    # for (i in seq(0,0.1,by=0.0005)) {
    #   if(length(q[q > lfcn & names(q) %in% p]) < min.regulons){
    #     next
    #   } else {
    #     break
    #   }
    # }
    # keep track of the lfc used for each condition
    lfc.scores[i] <- lfcn
  } else {
    lfcn <- min.lfcn
  }
  
  # keep regulons with lfc > lfcn and with good pval adjusted of diff activation
  lfc <- c(lfc,list(names(q[q > lfcn & names(q) %in% p])))
}
names(lfc) <- keys_ct

# check how many regulons were kept
t(as.data.frame(lapply(lfc,length)))

############
## make a list of the TFs (or regulons) that resulted after filtering
selectedTFs <- c()
groupInterest <- "SARS" # variable to take only TFs enriched in a certain group

for(i in 1:n_cases){
  #selectedTFs <- c(selectedTFs,lfc[[i]])
  ## take only those relevant to samples of more interest
  if(str_detect(names(lfc)[i],groupInterest)){
    selectedTFs <- c(selectedTFs,lfc[[i]])
  }
}
selectedTFs <- unique(sort(selectedTFs))

# make df of adj pvalues from DA analysis for this subset of TFs
indexTFs <- match(selectedTFs, tfs)
pv.adj.mwu.sel <- pv.adj.mwu[indexTFs,]
rownames(pv.adj.mwu.sel) <- str_remove(selectedTFs,"_\\(\\+\\)")
write.csv(pv.adj.mwu.sel, file = paste0(outdir,"selectedRegs_pvals-adj.csv"), quote = F)

# make df of enrichment from RSS for this subset of TFs
RSS <- py$rss_group
RSS.sel <- RSS[match(namesInRSS$cases,rownames(RSS)),indexTFs] # fltr for selected tfs and cases 
colnames(RSS.sel) <- str_remove(colnames(RSS.sel),"_\\(\\+\\)")
rownames(RSS.sel) <- str_remove(rownames(RSS.sel), "[YN]$")
if(celllines==1){
  rownames(RSS.sel) <- str_replace_all(rownames(RSS.sel), "-", "")
  rownames(RSS.sel) <- str_replace(rownames(RSS.sel), "w_vector_expressing_", "")
  rownames(RSS.sel) <- str_replace(rownames(RSS.sel), "1hr_Ruxolitinib_pretreatment", "pt")
}
write.csv(t(RSS.sel), file = paste0(outdir,"selectedRegs_RSSs.csv"), quote = F)

############
## Leo's data
py_run_string('path_regulons_targets_pkl = "/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/regulons_multiruns2.pkl"')
py_run_string('path_regulons_targets_csv = "/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/multi_runs_features_mtf.csv"')
## celllines
py_run_string('path_regulons_targets_pkl = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/celllines/out/scenic/celllines100/multi_runs_cistarget/multi_runs_regulons_mtf.pkl.gz"')
py_run_string('path_regulons_targets_csv = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/celllines/out/scenic/celllines100/multi_runs_cistarget/multi_runs_features_mtf.csv"')4
py_run_string('pkl_format = 0')
# CORRO ESTE SCRIPT PA AGARRAR LOS REGULONES Y SUS TARGETS
py_run_file("/mnt/Citosina/amedina/mpadilla/COVID19/PulmonDB_COVID/scRNAseq/pyScripts/TFsTargetGenes_ccmk.py")

# AQUÍ RECUPERO LOS TFS Y SUS TARGETS
TF_targets <- py$TF_targets
names(TF_targets)
for (i in 1:n_cases) {
  for (j in 1:length(names(TF_targets[[i]]))) {
    TF_targets[[i]][[j]] <- unlist(TF_targets[[i]][[j]])
  }
}


