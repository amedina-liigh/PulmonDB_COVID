import re
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

## WILCOXON TEST CONTRL VS COVID
# Performs wilcoxon u test of regulon activity (using auc matrix) between control and covid samples
#  for each celltype or source tissue
def mwuAUC(auc1,auc2):
    print("Performing Mann Whitney Test on ...")
    celltypes = sorted(list(auc1['source'].unique()))
    pvalues = {}
    for i in celltypes:
        print(i)
        cellTypePvalueVector = []
        for j in range((len(auc1.columns)-1)):
          #print(j)
          res = mwu(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], method = "asymptotic")
          cellTypePvalueVector.append(res[1])
        pvalues[i] = cellTypePvalueVector
    return pvalues


## KOLMOGOROV TEST CONTROL VS COVID
def ksAUC(auc1,auc2):
    print("Performing Kolmogrov Test on ...")
    celltypes = sorted(list(auc1['source'].unique()))
    pvalues = {}
    for i in celltypes:
     print(i)
     cellTypePvalueVector = []
     for j in range((len(auc1.columns)-1)):
       res = ks_2samp(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j])
       cellTypePvalueVector.append(res[1])
     pvalues[i] = cellTypePvalueVector
    return pvalues

# Functions adapted for celllines samples
# Celllines mocks are used for several experiments (different infections)
def mwuAUCcelllines(auc_cntrl,auc_case):
    print("Performing Mann Whitney Test on celllines ...")
    cases = sorted(list(auc_case['source'].unique()))
    celllines = sorted(list(auc_cntrl['cellline'].unique()))
    celllines.reverse() # more especific names should go first to pair with cases correctly
    pvalues = {}
    for i in cases:
        print(i)
        # check which cellline is to call corresponding controls
        for c in celllines:
         pttrn = re.compile(c)
         if pttrn.search(i):
             break
        print("\tcorresponding control : " + c)
        regulonsPvalueVector = []
        for j in range((len(auc_case.columns)-1)):
            #print(j)
            res = mwu(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j], method = "asymptotic")
            regulonsPvalueVector.append(res[1])
        pvalues[i] = casesPvalueVector
    return pvalues

def ksAUCcelllines(auc_cntrl,auc_case):
    print("Performing Kolmogrov Test on celllines ...")
    cases = sorted(list(auc_case['source'].unique()))
    celllines = sorted(list(auc_cntrl['cellline'].unique()))
    celllines.reverse()
    pvalues = {}
    for i in cases:
        print(i)
        regulonsPvalueVector = []
        # check which cellline is to call corresponding controls
        for c in celllines:
         pttrn = re.compile(c)
         if pttrn.search(i):
             break
        for j in range((len(auc_case.columns)-1)):
            res = ks_2samp(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j])
            regulonsPvalueVector.append(res[1])
        pvalues[i] = casesPvalueVector
    return pvalues

#def main():
#############################
## Read input and clean data
# input vars: path_meta, path_scenic_loom, filter_auc, celllines

# Read metadata or annotations ## first col is samples, second is source (of sample), third is disease flag
meta = pd.read_csv( path_meta, sep = ",")
if celllines :
    metadf = pd.DataFrame(meta).set_axis(['ID','source','disease','cellline'],axis=1)
else :
    metadf = pd.DataFrame(meta).set_axis(['ID','source','disease'],axis=1)
# dictionaries of different ways to group samples
meta_dis = pd.Series(data=dict(metadf[['ID','disease']].values), index=list(metadf['ID'])) # format as series
groups = list(metadf['source'])
    # concatenate source and disease state
for a in range(len(groups)):
    groups[a] = groups[a]+list(metadf['disease'])[a]

metadf['groups'] = groups
meta_group = pd.Series(data=dict(metadf[['ID','groups']].values), index=list(metadf['ID'])) # format as series

# Retrive regulon activity matrix (auc_mtx) from pyscenic
lf = lp.connect(path_scenic_loom, mode='r+', validate=False)
auc_mtx_multi = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

# in case filtering of auc_mtrx is necessary
if filter_auc:
    # filter out cells not included in metadata file
    auc_mtx_multi = auc_mtx_multi[auc_mtx_multi.index.isin(meta['ID'])]

#############################
## Compute RSS

print("Computing RSS ...")

# RSScores for each source tissue or celltype
rss_group = regulon_specificity_scores(auc_mtx_multi, meta_group)
# RSScores for control and covid
rss_dis = regulon_specificity_scores(auc_mtx_multi, meta_dis)

#############################
## Format data prior computation of Differentially Activated Regulons (DAR)

# Split annotations for control and cases
meta_control = meta[meta['disease']=='N']
meta_cov = meta[meta['disease']=='Y']

# Split auc matrix into controls and cases
auc_control = auc_mtx_multi[auc_mtx_multi.index.isin(meta_control['ID'])] ## auc matrix subset for control cells
auc_cov = auc_mtx_multi[auc_mtx_multi.index.isin(meta_cov['ID'])] ## auc matrix subset for covid cells

# Sort regulons alphabetically
auc_control = auc_control.reindex(sorted(auc_control.columns), axis=1)
auc_cov = auc_cov.reindex(sorted(auc_cov.columns), axis=1)

# Add cell type or tissue or experiment info
auc_control['source'] = list(meta_control['source'])
auc_cov['source'] = list(meta_cov['source'])

# Add cellline info
if celllines:
    auc_control['cellline'] = list(meta_control['cellline'])

#############################
## WILCOXON TEST CONTRL VS COVID, adjust pvalues

# Get pvalues from test
if celllines:
    pvalues = mwuAUCcelllines(auc_control,auc_cov)
else :
    pvalues = mwuAUC(auc_control,auc_cov)

print("Adjusting pvalues ...")

# AJUSTO LOS PVALUES Y TAMBIÉN LOS HAGO -LOG10 Y LOS QUE PASEN EL LÍMITE DEL PRIMER CUARTIL
# Adjust pvalues for multiple testing
pvalues_adjusted = {}
for i in pvalues:
  pvalues_adjusted[i] = multest(pvals = pvalues[i],method = 'fdr_bh')[1]

# -LOG10 of adjusted pvalues
pvalues_adjusted_minusLog10 = {}
for i in pvalues_adjusted:
  pvalues_adjusted_minusLog10[i] = list(map(lambda x: -1*math.log10(x) if x !=0 else x,pvalues_adjusted[i]))

# LÍMITE DEL PRIMER CUARTIL
regsPass_1stQ = {}
for i in pvalues_adjusted:
  regsPass_1stQ[i] = list(pvalues_adjusted[i] <= np.quantile(pvalues_adjusted[i],.25))

#############################
## KOLMOGOROV TEST CONTROL VS COVID

# Get pvalues from test
if celllines:
    pvalues2 = ksAUCcelllines(auc_control,auc_cov)
else :
    pvalues2 = ksAUC(auc_control,auc_cov)

print("Adjusting pvalues ...")

pvalues_adjusted2 = {}
for i in pvalues2:
  pvalues_adjusted2[i] = multest(pvals = pvalues2[i],method = 'fdr_bh')[1]

pvalues_adjusted_minusLog102 = {}
for i in pvalues_adjusted2:
  pvalues_adjusted_minusLog102[i] = list(map(lambda x: -1*math.log10(x) if x !=0 else x,pvalues_adjusted2[i]))

regsPass_1stQ2 = {}
for i in pvalues_adjusted2:
  regsPass_1stQ2[i] = list(pvalues_adjusted2[i] <= np.quantile(pvalues_adjusted2[i],.25))

#if __name__ == '__main__':
#    main()
