import os
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

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

## WILCOXON TEST CONTRL VS COVID
# Performs wilcoxon u test of regulon activity (using auc matrix) between control and covid samples
#  for each celltype or source tissue
def mwuAUC(auc1,auc2, directionals=0):
    print("Performing Mann Whitney Test on ...")
    celltypes = sorted(list(auc1['source'].unique()))
    pvalues = {}
    for i in celltypes:
        print(i)
        # make lists for different tests if directional alternative hypothesis is chosen
        if directionals :
            regulonsPvalueVectorUp = []
            regulonsPvalueVectorDown = []
        else :
            cellTypePvalueVector = [] # non-directional case
        for j in range((len(auc1.columns)-1)):
          #print(j)
          if directionals :
              res = mwu(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], method = "asymptotic", alternative = "greater")
              regulonsPvalueVectorUp.append(res[1])
              res = mwu(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], method = "asymptotic", alternative = "less")
              regulonsPvalueVectorDown.append(res[1])
          else :
              res = mwu(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], method = "asymptotic")
              cellTypePvalueVector.append(res[1])
        if directionals :
            pvalues[i] = [regulonsPvalueVectorUp,regulonsPvalueVectorDown]
        else :
            pvalues[i] = cellTypePvalueVector
    return pvalues


## KOLMOGOROV TEST CONTROL VS COVID
def ksAUC(auc1,auc2, directionals=0):
    print("Performing Kolmogrov Test on ...")
    celltypes = sorted(list(auc1['source'].unique()))
    pvalues = {}
    for i in celltypes:
     print(i)
     if directionals :
         regulonsPvalueVectorUp = []
         regulonsPvalueVectorDown = []
     else :
         cellTypePvalueVector = [] # non-directional case
     for j in range((len(auc1.columns)-1)):
        if directionals :
            res = ks_2samp(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], alternative = "greater")
            regulonsPvalueVectorUp.append(res[1])
            res = ks_2samp(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j], alternative = "less")
            regulonsPvalueVectorDown.append(res[1])
        else :
            res = ks_2samp(auc1[auc1['source'] == i].iloc[:,j],auc2[auc2['source'] == i].iloc[:,j])
            cellTypePvalueVector.append(res[1])
     if directionals :
         pvalues[i] = [regulonsPvalueVectorUp,regulonsPvalueVectorDown]
     else :
         pvalues[i] = cellTypePvalueVector
    return pvalues

# Functions adapted for celllines samples
# Celllines mocks are used for several experiments (different infections)
def mwuAUCcelllines(auc_cntrl,auc_case, directionals=0):
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
        # make lists for different tests if directional alternative hypothesis is chosen
        if directionals :
            regulonsPvalueVectorUp = []
            regulonsPvalueVectorDown = []
        else :
            regulonsPvalueVector = [] # non-directional case
        for j in range((len(auc_case.columns)-1)):
            #print(j)
            if directionals :
                res = mwu(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j], alternative = "greater")
                regulonsPvalueVectorUp.append(res[1])
                res = mwu(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j], alternative = "less")
                regulonsPvalueVectorDown.append(res[1])
            else :
                res = mwu(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j])
                regulonsPvalueVector.append(res[1])
        if directionals :
            pvalues[i] = [regulonsPvalueVectorUp,regulonsPvalueVectorDown]
        else :
            pvalues[i] = regulonsPvalueVector
    return pvalues

def ksAUCcelllines(auc_cntrl,auc_case, directionals=0):
    print("Performing Kolmogrov Test on celllines ...")
    cases = sorted(list(auc_case['source'].unique()))
    celllines = sorted(list(auc_cntrl['cellline'].unique()))
    celllines.reverse()
    pvalues = {}
    for i in cases:
        print(i)
        # make lists for different tests if directional alternative hypothesis is chosen
        if directionals :
            regulonsPvalueVectorUp = []
            regulonsPvalueVectorDown = []
        else :
            regulonsPvalueVector = [] # non-directional case
        # check which cellline is to call corresponding controls
        for c in celllines:
         pttrn = re.compile(c)
         if pttrn.search(i):
             break
        for j in range((len(auc_case.columns)-1)):
            if directionals :
                res = ks_2samp(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j], alternative = "greater")
                regulonsPvalueVectorUp.append(res[1])
                res = ks_2samp(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j], alternative="less")
                regulonsPvalueVectorDown.append(res[1])
            else :
                res = ks_2samp(auc_case[auc_case['source'] == i].iloc[:,j],auc_cntrl[auc_cntrl['cellline'] == c].iloc[:,j])
                regulonsPvalueVector.append(res[1])
        if directionals :
            pvalues[i] = [regulonsPvalueVectorUp,regulonsPvalueVectorDown]
        else :
            pvalues[i] = regulonsPvalueVector
    return pvalues

def adj_pvals(pvalues) :
    print("Adjusting pvalues ...")
    pvalues_adjusted = {}
    for i in pvalues:
      pvalues_adjusted[i] = multest(pvals = pvalues[i],method = 'fdr_bh')[1]
    return pvalues_adjusted

# Log2 fold change calculation
def l2fc_(auc1,auc2):
  celltypes = sorted(list(set(auc1['celltype'].unique()).intersection(set(auc2['celltype'].unique()))))
  l2fc = {}
  for i in celltypes:
    cellTypeL2fcVector = []
    for j in range((len(auc1.columns)-1)):
      mean1 = (np.mean(auc1[auc1['celltype'] == i].iloc[:,j])+1)
      mean2 = (np.mean(auc2[auc2['celltype'] == i].iloc[:,j])+1)
      res = np.log2(mean1/mean2)
      cellTypeL2fcVector.append(res)
    l2fc[i] = cellTypeL2fcVector
  return l2fc

# Build dataframe containing P values, adjusted P values, log2 fold changes, and regulones
def dfs_(auc1,auc2,pvalues,pvalues_adjusted,l2fc):
  celltypes = sorted(list(set(auc1['celltype'].unique()).intersection(set(auc2['celltype'].unique()))))
  dfs = {}
  for i in celltypes:
    df = pd.DataFrame({"reg":auc1.columns[0:(len(auc1.columns)-1)],"pval":pvalues[i],"adj_pval":pvalues_adjusted[i],"l2fc":l2fc[i]})
    dfs[i] = df
  return dfs

#def main():
#############################
## Read input and clean data
# input vars: path_meta, path_scenic_loom, filter_auc, filter_auc_regs, path_regulons_counts, min_reg_counts, test_sep, min_gene_occ

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

# create a dictionary of regulons (no filter):
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

# regulons gene occurance
regGeneOcc = pd.DataFrame(lf.ra.RegulonGeneOccurrences, columns = list(auc_mtx_multi.columns), index= list(lf.ra.Gene))

lf.close()

# in case filtering of auc_mtrx is necessary
if filter_auc_samples:
    # filter out cells not included in metadata file
    auc_mtx_multi = auc_mtx_multi[auc_mtx_multi.index.isin(meta['ID'])]

# apply regulons counts (from scenic iterations) cutoff
if filter_auc_regs :
    # read regulons.tsv multiruns output file
    reg_counts = pd.read_csv( path_regulons_counts, sep = "\t").set_axis(['regulon','count'],axis=1)
    # choose min_reg_counts
    if min_reg_counts is None:
        min_reg_counts = reg_counts['count'][0]/10
    # regs that pass cutoff
    regs = list(reg_counts[reg_counts['count'] >= min_reg_counts]['regulon'])
    regs_ = []
    for reg in regs: regs_.append(re.sub("\(\+\)","_(+)",reg))
    # apply regs filter
    auc_mtx_multi = auc_mtx_multi.T[auc_mtx_multi.columns.isin(regs_)].T
    # apply regs filter in regulons list
    regGeneOcc = regGeneOcc.T[regGeneOcc.T.index.isin(regs_)].T
    regulons_fltr = {}
    for tf in regGeneOcc:
        regulons_fltr[tf] = list(regGeneOcc[tf][regGeneOcc[tf] > min_gene_occ].index)

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

# Get pvalues from test for DA down and upregulated separetely
if test_sep :
    if celllines:
        pvalues = mwuAUCcelllines(auc_control,auc_cov, directionals=1)
    else :
        pvalues = mwuAUC(auc_cov, auc_control, directionals=1)
else :
    # Get pvalues from test for DA down or upregulated
    if celllines:
        pvalues = mwuAUCcelllines(auc_control,auc_cov)
    else :
        pvalues = mwuAUC(auc_control,auc_cov)

# split dictionary if directional test were performed to facilitate the downstream pipeline
if test_sep :
    pvalues_up = {}
    pvalues_down = {}
    for i in pvalues:
        pvalues_up[i] = pvalues[i][0]
        pvalues_down[i] = pvalues[i][1]

# AJUSTO LOS PVALUES Y TAMBIÉN LOS HAGO -LOG10 Y LOS QUE PASEN EL LÍMITE DEL PRIMER CUARTIL
# Adjust pvalues for multiple testing
if test_sep :
    pvalues_up_adjusted = adj_pvals(pvalues_up)
    pvalues_down_adjusted = adj_pvals(pvalues_down)
else :
    pvalues_adjusted = adj_pvals(pvalues)

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
if test_sep :
    if celllines:
        pvalues2 = ksAUCcelllines(auc_control,auc_cov, directionals=1)
    else :
        pvalues2 = ksAUC(auc_cov,auc_control, directionals=1)
else :
    if celllines:
        pvalues2 = ksAUCcelllines(auc_control,auc_cov)
    else :
        pvalues2 = ksAUC(auc_control, auc_cov)

# split dictionary if directional test were performed to facilitate the downstream pipeline
if test_sep :
    pvalues2_up = {}
    pvalues2_down = {}
    for i in pvalues2:
        pvalues2_up[i] = pvalues[i][0]
        pvalues2_down[i] = pvalues[i][1]

if test_sep :
    pvalues_up_adjusted2 = adj_pvals(pvalues2_up)
    pvalues_down_adjusted2 = adj_pvals(pvalues2_down)
else :
    pvalues_adjusted2 = adj_pvals(pvalues2)

    pvalues_adjusted_minusLog102 = {}
    for i in pvalues_adjusted2:
      pvalues_adjusted_minusLog102[i] = list(map(lambda x: -1*math.log10(x) if x !=0 else x,pvalues_adjusted2[i]))

    regsPass_1stQ2 = {}
    for i in pvalues_adjusted2:
      regsPass_1stQ2[i] = list(pvalues_adjusted2[i] <= np.quantile(pvalues_adjusted2[i],.25))

#############################
## LOG2 FOLD CHANGE TEST CONTROL VS COVID
l2fc = l2fc_(auc_cov,auc_control)

#############################
## BUILD DATAFRAME CONTAINING P VALUES, ADJUSTED P VALUES, LOG2 FOLD CHANGES, AND REGULONS
dfs_balf = dfs_(auc_cov,auc_control,pvalues,pvalues_adjusted,l2fc)

#if __name__ == '__main__':
#    main()
