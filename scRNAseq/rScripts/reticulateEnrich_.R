### SET CONDA ENV TO USE WITH RETICULATE
library(reticulate)
library(dplyr)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
use_condaenv(condaenv = "scanpy", conda = "/cm/shared/apps/anaconda3/2021.05/envs/scanpy", required = TRUE)
py_config()

# CORRO EL SCRIPT "difregs_.py" PA USAR LOS RESULTADOS EN R
py_run_file("difregs_.py")

# GUARDO LAS VARIABLES DE PYTHON A R
keys_ct = names(py$pvalues_adjusted)
rss_both_ctc = r_to_py(py$rss_both_ctc)

# AQUÍ AGARRO LOS QUE SALIERON SIGNIFACTIVOS TANTO CON WILCOXON COMO CON KOLMOGOROV, PERO POS NO ES NECESARIO
# CON UNO DE LOS TEST BASTA
tfs_mwu_ks = list()
for (i in 1:9) {
  seto1 = names(py$auc_control)[names(py$auc_control) != 'celltype'][unname(unlist(py$regsPass_1stQ[[i]]))]
  seto2 = names(py$auc_control)[names(py$auc_control) != 'celltype'][unname(unlist(py$regsPass_1stQ2[[i]]))]
  tfs_mwu_ks = c(tfs_mwu_ks,list(intersect(seto1,seto2)))
}
names(tfs_mwu_ks) <- names(py$regsPass_1stQ)

# AQUÍ YA PODRÍAS TENER LOS REGULONES QUE ESTÁN DIFERENCIALMENTE ACTIVOS (DA)
# PERO LO QUE HICE AQUÍ FUE DE LOS QUE SALIERON DA, SELECCIÓN LOS QUE ESTABAN MÁS ACTIVOS EN CADA TIPO CELULAR
# SEGÚN EL RSS
ccmk = list()
for (i in 1:9) {
  top_rss_cov_Y = py_to_r(rss_both_ctc$T$nlargest(length(names(py$auc_control)),paste0(keys_ct[i],'_Y'))[paste0(keys_ct[i],'_Y')])
  top_rss_cov_N = py_to_r(rss_both_ctc$T$nlargest(length(names(py$auc_control)),paste0(keys_ct[i],'_Y'))[paste0(keys_ct[i],'_N')])
  temp <- top_rss_cov_Y > summary(top_rss_cov_Y)[2]
  top_rss_cov_Y = top_rss_cov_Y[temp]
  top_rss_cov_N = top_rss_cov_N[temp]
  top_rss_cov_Y = top_rss_cov_Y[top_rss_cov_Y > top_rss_cov_N]
  ccmk = c(ccmk,list(intersect(names(top_rss_cov_Y),tfs_mwu_ks[[i]])))
}
names(ccmk) <- names(tfs_mwu_ks)

# TAMBIÉN PODRÍAS HACER UN PSEUDO-FOLD CHANGE PA VER CUALES SALEN
lfc <- list()
for (i in 1:9) {
  p <- names(py$auc_control)[unlist(py$regsPass_1stQ[keys_ct[i]])]
  lfc_Y = py_to_r(rss_both_ctc$T$nlargest(length(names(py$auc_control)),paste0(keys_ct[i],'_Y'))[paste0(keys_ct[i],'_Y')])
  lfc_N = py_to_r(rss_both_ctc$T$nlargest(length(names(py$auc_control)),paste0(keys_ct[i],'_Y'))[paste0(keys_ct[i],'_N')])
  q <- log2(lfc_Y/lfc_N)
  if (i == 1 | i == 8 | i == 9) {
    lfc <- c(lfc,list(names(q[q > 0.05 & names(q) %in% p]))) 
  } else if (i == 4){
    lfc <- c(lfc,list(names(q[q > 0.01 & names(q) %in% p]))) 
  } else {
    lfc <- c(lfc,list(names(q[q > 0.1 & names(q) %in% p])))
  }
}
names(lfc) <- keys_ct
str(lfc)

# CORRO ESTE SCRIPT PA AGARRAR LOS REGULONES Y SUS TARGETS
py_run_file("TFsTargetGenes_ccmk.py")

# AQUÍ RECUPERO LOS TFS Y SUS TARGETS
TF_targets <- py$TF_targets
names(TF_targets)
for (i in 1:9) {
  for (j in 1:length(names(TF_targets[[i]]))) {
    TF_targets[[i]][[j]] <- unlist(TF_targets[[i]][[j]])
  }
}

