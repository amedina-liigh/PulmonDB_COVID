import math
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
from pyscenic.plotting import plot_rss
from pyscenic.rss import regulon_specificity_scores
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests as multest
from itertools import compress

# AQUÍ LEO LOS DATOS Y LES AGREGO INFO, TAMBIÉN LOS FILTRO PORQUE EL DATASET TIENE MÁS CÉLULAS QUE LAS QUE VIENEN EN LA METADATA
exp_mat_both = sc.read_loom("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/exp_mat_both.loom", validate=False)

cell_annotations_both = pd.read_csv("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/all.cell.annotation.meta.txt",sep="\t")

exp_mat_both = exp_mat_both[exp_mat_both.obs.index.isin(cell_annotations_both['ID'])]

exp_mat_both.obs['celltype'] = list(cell_annotations_both[['celltype']].iloc[:,0])
exp_mat_both.obs['ctc'] = [i+'_'+str(j) for i,j in zip(list(cell_annotations_both[['celltype']].iloc[:,0]),list(cell_annotations_both[['disease']].iloc[:,0]))]

lf = lp.connect("/home/amedina/amedinalab/larteaga/COVID19/pyscenic/authors/recuperacion/avg3/SCENIC_final_multi.loom", mode='r+', validate=False)
auc_both_multi = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

auc_both_multi = auc_both_multi[auc_both_multi.index.isin(cell_annotations_both['ID'])]

# AQUÍ SACO LOS RSS
rss_both = regulon_specificity_scores(auc_both_multi, exp_mat_both.obs.celltype)
rss_both_ctc = regulon_specificity_scores(auc_both_multi, exp_mat_both.obs.ctc)

# HAGO LAS MATRICES DE CONTROL Y COVID A PARTIR DE LA DE LOS DATOS INTEGRADOS
cell_annotations_control = cell_annotations_both[cell_annotations_both['disease']=='N']
cell_annotations_cov = cell_annotations_both[cell_annotations_both['disease']=='Y']

auc_control = auc_both_multi[auc_both_multi.index.isin(cell_annotations_control['ID'])]
auc_cov = auc_both_multi[auc_both_multi.index.isin(cell_annotations_cov['ID'])]

# Sort regs alphabetically
auc_control = auc_control.reindex(sorted(auc_control.columns), axis=1)
auc_cov = auc_cov.reindex(sorted(auc_cov.columns), axis=1)

# Add cell type info
auc_control['celltype'] = list(cell_annotations_control['celltype'])
auc_cov['celltype'] = list(cell_annotations_cov['celltype'])

# Remove the Neutrophil cell type form cov patients samples because it is not present in control samples
auc_cov = auc_cov[auc_cov['celltype'] != 'Neutrophil']

# AQUÍ HAGO UNA FUNCIÓN QUE VA A IR HACIENDO EL WILCOXON TEST CONTRL VS COVID EN CADA TIPO CELULAR
# USANDO LOS DATOS DE LAS CÉLULAS DE LA MATRIX DE ACTIVIDADES DE LOS REGULONES, LA ÚLTIMA QUE ARROJA SCENIC
def mwuAUC(auc1,auc2):
  celltypes = sorted(list(auc1['celltype'].unique()))
  pvalues = {}
  for i in celltypes:
    print(i)
    cellTypePvalueVector = []
    for j in range((len(auc1.columns)-1)):
      #print(j)
      res = mwu(auc1[auc1['celltype'] == i].iloc[:,j],auc2[auc2['celltype'] == i].iloc[:,j], method = "asymptotic")
      cellTypePvalueVector.append(res[1])
    pvalues[i] = cellTypePvalueVector
  return pvalues

pvalues = mwuAUC(auc_control,auc_cov)

# AJUSTO LOS PVALUES Y TAMBIÉN LOS HAGO -LOG10 Y LOS QUE PASEN EL LÍMITE DEL PRIMER CUARTIL

pvalues_adjusted = {}
for i in pvalues:
  pvalues_adjusted[i] = multest(pvals = pvalues[i],method = 'fdr_bh')[1]

# -LOG10
pvalues_adjusted_minusLog10 = {}
for i in pvalues_adjusted:
  pvalues_adjusted_minusLog10[i] = list(map(lambda x: -1*math.log10(x) if x !=0 else x,pvalues_adjusted[i]))
  
# LÍMITE DEL PRIMER CUARTIL
regsPass_1stQ = {}
for i in pvalues_adjusted:
  regsPass_1stQ[i] = list(pvalues_adjusted[i] <= np.quantile(pvalues_adjusted[i],.25))

# AQUÍ EN VEZ DEL WILCONXON, USO EL TEST DE KOLMOGOROV. NO HACE FALTA QUE USES LOS DOS, CON UNO BASTA
def ksAUC(auc1,auc2):
  celltypes = sorted(list(auc1['celltype'].unique()))
  pvalues = {}
  for i in celltypes:
    cellTypePvalueVector = []
    for j in range((len(auc1.columns)-1)):
      res = ks_2samp(auc1[auc1['celltype'] == i].iloc[:,j],auc2[auc2['celltype'] == i].iloc[:,j])
      cellTypePvalueVector.append(res[1])
    pvalues[i] = cellTypePvalueVector
  return pvalues

pvalues2 = ksAUC(auc_control,auc_cov)

pvalues_adjusted2 = {}
for i in pvalues2:
  pvalues_adjusted2[i] = multest(pvals = pvalues2[i],method = 'fdr_bh')[1]

pvalues_adjusted_minusLog102 = {}
for i in pvalues_adjusted2:
  pvalues_adjusted_minusLog102[i] = list(map(lambda x: -1*math.log10(x) if x !=0 else x,pvalues_adjusted2[i]))

regsPass_1stQ2 = {}
for i in pvalues_adjusted2:
  regsPass_1stQ2[i] = list(pvalues_adjusted2[i] <= np.quantile(pvalues_adjusted2[i],.25))

