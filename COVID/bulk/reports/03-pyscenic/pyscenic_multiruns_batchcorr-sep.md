PulmonDB - COVID19.
================
Monica Padilla

  

## **Multiruns pyscenic on Blanco 2020 + Desai 2020 + Delorey 2021 data**

This report describes how the scenic multiruns vsn-pipeline was run for
a merged dataset from human autopsies or cell lines samples of datasets
by Blanco-Melo et al 2020, Desai et al, 2020 and bulk samples from
Delorey et al., 2021. 3 analysis runs: celllines, organs-desgtex and
lung-bladesdel. Data download and pre-processing described in another
report. Does not include Blanco’s dataset Ferret samples.

  

### **Prepare terminal**

``` bash
ssh -Y mpadilla@dna.liigh.unam.mx
screen -S bladesdel
qlogin
cd /mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/
module load nextflow/20.10.0
module load singularity/3.7.0
module unload r/3.2.1
module load r/4.0.2
umask 2
```

  

### **Prepare Input data**

  

**Convert csv to loom file**

``` bash
mkdir data
module load anaconda3/2021.05
source activate scanpy # loompy 3.0.6
python3
```

Multiruns pipeline requires a loom file.

  

at =
“/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns3/bladesdel”

  

#### Create loom file with loompy version 3.0.6 (install with pip)

-   Adapted from :
    <https://linnarssonlab.org/loompy/apiwalkthrough/index.html>
-   All count matrixes were corrected by batch effect

``` python
import pandas as pd
import numpy as np
import loompy # verify version >= 3
# 

## lung-bladesdel 
file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs-gfltr.csv"
## celllines all bc dates
file =  "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.csv"
## celllines all bc dates no SARS-CoV-2-infected-A549 +ACE2 +pt since no control
file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE_nopt-counts-fltrtrxs.csv"
## multiorgans covid healthy desai gtex
file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_countsfltrtrxs.csv"

######
## lung-blades batch correction on dates
file = "blanco2020_adj-desai2020_adj-Lung-counts-fltrtrxs.csv.gz"
## lung-blades-2bc 2 batch corrections: date and dataset
file = "blanco2020-desai2020-mer-adj-counts-Lung-fltrtrxs.csv.gz"


## celllines without hINFB treated NHBE
file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_noINFBNHBE-counts-fltrtrxs.csv"
## multiorgans covid - desai
file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020-counts-adj-selCovidOrgans_fltrtrxs.csv.gz"
## multiorgans healthy - gtex


# read counts file 
data = pd.read_csv(file) # this csv is cols genes and rows samples
loomfile = file.replace(".csv",".loom") # outfile name

# genes rows and cells cols (samples in this case)
matrix = np.array(data.T.values)
row_attrs = { "Gene" : list(data.columns) }
col_attrs = { "CellID" : list(data.index) }
loompy.create(loomfile, matrix, row_attrs, col_attrs)
```

------------------------------------------------------------------------

  

### **Multiruns pipeline**

  

Pipeline described at:
<https://vsn-pipelines-examples.readthedocs.io/en/latest/PBMC10k_multiruns.html>
and
<https://vsn-pipelines.readthedocs.io/en/latest/pipelines.html#scenic-scenic>

  

#### **Load pipeline and create config file**

Load pipelines of `nextflow v0.25.0` :

``` bash
nextflow pull vib-singlecell-nf/vsn-pipelines -r v0.25.0
```

  

Use `scenic` entry pipeline, specify `scenic_multiruns` to load
multiruns parameters, specify `cistarget` to load params for scenic
input files and dbs, genome and `singularity` container type instead of
docker.

``` bash
nextflow config vib-singlecell-nf/vsn-pipelines \
    -profile loom,scenic,scenic_multiruns,scenic_use_cistarget_motifs,scenic_use_cistarget_tracks,hg38,singularity \
    > bladesdel-scenic.config
    
nano blades-scenic-261.config
```

  

#### **Make changes to** `.config` **file:**

Change general parameters:

-   in params.global add : `seed = 777`
-   in params.global change : `project_name = 'bladesdel-l3'`

Specify input loom file :

-   at params.sc.scenic :
    `filteredLoom = '/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/blanco2020_minusFerr-desai2020-delorey2021-counts-l3.loom'`

Change `params.sc.scenic.grn` and `params.sc.scenic.cisTarget`
parameters to take motif dbs and other resources from the ones
previously saved in the cluster.

-   at `params.sc.scenic.grn` change :
    `tfs = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/demo/resources/allTFs_hg38.txt'`

-   at `params.sc.scenic.cisTarget` change :

<!-- -->

    motifsDb = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/demo/resources/hg38__refseq-r80__*feather'
    motifsAnnotation = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/demo/resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
    tracksDb = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/demo/resources/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__*feather'
    tracksAnnotation = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/demo/resources/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv'

Multiruns-specific parameters:

-   number of runs of scenic : `params.sc.scenic.numRuns = 2`
-   minimum genes per regulon :
    `params.sc.scenic.aucell.min_genes_regulon = 5`
-   minimum ocurrences of genes at a given regulon through runs :
    `params.sc.scenic.aucell.min_regulon_gene_occurrence = 2`

Change input file path and suffix

-   At `params.data.csv` :

<!-- -->

             file_paths = '/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/blanco2020_minusFerr-desai2020-delorey2021-counts-l3.loom'
             suffix = '-counts-l3.loom'

Change executor params for grn since it requires more resources:

-   At `process.withLabel:compute_resources__scenic_grn` :
    `executor = 'qsub'` and `cpus = 19`

At last:

-   at singularity add :
    `cacheDir = '/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns3/bladesdel/work/singularity'`
-   at singularity change :
    `runOptions = '--cleanenv -H $PWD -B ${HOME}'`

  

#### **Run pipeline**

``` bash
############
# with loom file 3.0.2 and scenic 10.4
# bladesdel
nextflow -C lung_bladesdel-scenic10.config \
    run vib-singlecell-nf/vsn-pipelines \
    -entry scenic \
    -r v0.25.0
```

  

If an error occurs, fix it and then resume:

``` bash
nextflow -C bladesdel-scenic.config \
    run vib-singlecell-nf/vsn-pipelines \
    -entry scenic \
    -resume -r v0.25.0
```

##### **Jobs**

``` bash
#!/bin/bash
# Use current working directory
#$ -cwd
#
# Join stdout and stderr
#$ -j y
#
# Run job through bash shell
#$ -S /bin/bash
#
#You can edit the scriptsince this line
#
# Your job name
#$ -N celllines100
#
# Send an email after the job has finished
#$ -m e
#$ -M monicapadilla905@gmail.com
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
module load nextflow/20.10.0
module load singularity/3.7.0
#
# Write your commands in the next line
nextflow -C celllines-scenic100.config \
    run vib-singlecell-nf/vsn-pipelines \
    -entry scenic \
    -r v0.25.0
```

stdout example:

    [[8b/bb6102] process > scenic:SCENIC:ARBORETO_WITH_MULTIPROCESSING (10)        [100%] 10 of 10 ✔
    [21/3bc0f1] process > scenic:SCENIC:ADD_PEARSON_CORRELATION (10)              [100%] 10 of 10 ✔
    [c8/cf4778] process > scenic:SCENIC:CISTARGET__MOTIF (10)                     [100%] 10 of 10 ✔
    [47/82a008] process > scenic:SCENIC:AUCELL__MOTIF (10)                        [100%] 10 of 10 ✔
    [f4/c09e5b] process > scenic:SCENIC:MULTI_RUNS_TO_LOOM__MOTIF:AGGR_FEATURE... [100%] 1 of 1 ✔
    [86/dfb96d] process > scenic:SCENIC:MULTI_RUNS_TO_LOOM__MOTIF:AGGR_REGULON... [100%] 1 of 1 ✔
    [95/24ec26] process > scenic:SCENIC:MULTI_RUNS_TO_LOOM__MOTIF:AUCELL (1)      [100%] 1 of 1 ✔
    [b0/e6b540] process > scenic:SCENIC:MULTI_RUNS_TO_LOOM__MOTIF:FEATURES_TO_... [100%] 1 of 1 ✔            [f8/c99860] process > scenic:SCENIC:MULTI_RUNS_TO_LOOM__MOTIF:SAVE_TO_LOOM... [100%] 1 of 1 ✔            [28/51c123] process > scenic:SCENIC:VISUALIZE (1)                             [100%] 1 of 1 ✔            [2b/3e76de] process > scenic:SCENIC:PUBLISH_LOOM (1)                          [100%] 1 of 1 ✔            [85/4474e3] process > scenic:PUBLISH_SCENIC:COMPRESS_HDF5 (1)                 [100%] 1 of 1 ✔            [e6/28985a] process > scenic:PUBLISH_SCENIC:SC__PUBLISH (1)                   [100%] 1 of 1 ✔            Pulling Singularity image docker://vibsinglecellnf/hdf5:1.10.5-r2 [cache /mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep/organs-desgtex/work/singularity/vibsinglecellnf-hdf5-1.10.5-r2.img]
    Completed at: 16-Dec-2021 09:55:28
    Duration    : 1d 3m 22s
    CPU hours   : 117.4
    Succeeded   : 49

For 10 runs

-   lung blades

<!-- -->

    Completed at: 21-Dec-2021 15:25:11
    Duration    : 6d 3h 21m 43s
    CPU hours   : 732.9
    Succeeded   : 49

-   organs-desgtex

<!-- -->

    Completed at: 16-Dec-2021 09:55:28
    Duration    : 1d 3m 22s
    CPU hours   : 117.4
    Succeeded   : 49

-   celllines

<!-- -->

    Completed at: 16-Dec-2021 15:33:34
    Duration    : 1d 5h 46m 47s
    CPU hours   : 160.0
    Succeeded   : 49

  

#### **Configs, results files and run info**

This pipeline was run 100 times for each independent analysis:

Main dir =
`/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-bc-sep`

-   `celllines`
    -   A549, NHBE and Calu-3 mock and infected samples for virus: RSV,
        IAV, IAVdNS1, HPIV3, SARS-CoV-2 and expressing hACE2 or
        inhibiting INF
    -   data from Blanco 2020
    -   config file = `celllines/celllines-scenic100.config`
    -   input data =
        `celllines/data/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.loom`
    -   SCENIC output loom =
        `celllines/out/scenic/celllines100/SCENIC_SCope_output.loom`
    -   run info : `Duration    : 11d 11h 19m 14s`,
        `CPU hours   : 1'525.3`
-   `lung-bladesdel`
    -   Lung biopsies healthy and COVID
    -   data from Blanco 2020, Desai 2020, Delorey 2021
    -   config file = `lung-bladesdel/lung_bladesdel-scenic100.config`
    -   input data =
        `lung-bladesdel/data/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs-gfltr.loom`
    -   SCENIC output loom =
        `lung-bladesdel/out/scenic/lung-bdd100/SCENIC_SCope_output.loom`
    -   run info : `Duration    : 10d 15h 35m 30s`,
        `CPU hours   : 1'303.1`
-   `organs-desgtex`
    -   Lung, Heart, Liver, Bowel and Kidney biopies for healthy or
        COVID-19 diseased patients
    -   data from Desai 2020 for COVID-19 and GTEx 2020 for healthy
    -   config file = `organs-desgtex/organs_desgtex-scenic100.config`
    -   input data =
        `organs-desgtex/data/desai2020_gtexv8-subTissues-adj_countsfltrtrxs.loom`
    -   SCENIC output loom =
        `organs-desgtex/out/scenic/org_desgtex100/SCENIC_SCope_output.loom`
    -   run info : `Duration    : 6d 12h 57m 31s`, `CPU hours   : 791.2`

------------------------------------------------------------------------

### **Get enriched regulons**

  

pySCENIC pipeline outputs regulons enriched per sample, but we are
interested in regulons enriched in given biological conditions, i.e. in
a group of samples.

Two steps are used here to subset enriched regulons activated in certain
conditions:

1.  Regulon Specificity Score (RSS) (described at \[ref\])

To assign scores to regulons more activated in certain groups of
samples.

2.  Compute Wilcoxon test (according to RSS)

To select differentialy activated regulons accroding to 1.

  

#### **Regulon Specificity Score**

Adapted or copied from:
<https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_downstream-analysis.ipynb>

  

**Load packages and data**

``` python
import pandas as pd
import loompy as lp
from pyscenic.rss import regulon_specificity_scores

## scenic mutiruns output loom
path_scenic_loom = "/mnt/Citosina/amedina/mpadilla/COVID19/pyscenic/multiruns-sep/celllines-IN/out/scenic/celllines_wIN50/SCENIC_SCope_output.loom"
## metadata csv
path_meta = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease_subCelllines.csv"
meta = pd.read_csv( path_meta, sep = ",")
meta_df = pd.DataFrame(meta, columns = ['ID','source','disease'])
meta_se = pd.Series(data=dict(meta_df[['ID','source']].values), index=list(meta_df['ID'])) # format as series

## get regulon activity matrix, final scenic output
lf = lp.connect(path_scenic_loom, mode='r+', validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
```

  

**Compute RSS**

``` python
# Compute RSS
rss_exp = regulon_specificity_scores( auc_mtx2, meta )
rss_exp
```

Top 5 regulons per sample

``` python
cats = sorted(list(set(meta)))
cats_sc2 = sorted(list(set(meta_sc2)))

topreg5 = []
for i,c in enumerate(cats):
    topreg5.extend(
        list(rss_exp.T[c].sort_values(ascending=False)[:5].index)
    )


topreg5 = list(set(topreg5))

## top 2
topreg2 = []
for i,c in enumerate(cats):
    topreg2.extend(
        list(rss_exp.T[c].sort_values(ascending=False)[:2].index)
    )
    
    
topreg2 = list(set(topreg2))


auc_mtx_sc2_Z = pd.DataFrame( index=auc_mtx_sc2.index )
for col in list(auc_mtx_sc2.columns):
    auc_mtx_sc2_Z[ col ] = ( auc_mtx_sc2[col] - auc_mtx_sc2[col].mean()) / auc_mtx_sc2[col].std(ddof=0)
```

RSSs plot

``` python
## for all blades 
matplotlib.use('png')
fig = plt.figure(figsize=(15, 8))
for c,num in zip(cats_sc2, range(1,len(cats_sc2)+1)):
    x=rss_exp_sc2.T[c]
    ax = fig.add_subplot(2,5,num)
    plot_rss(rss_exp_sc2, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig("blades-RSS-top5.png", dpi=600, bbox_inches = "tight")
```

Generate heatmap

``` python
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f
  
import matplotlib
matplotlib.use('png')
import matplotlib.pyplot as plt


colors = sns.color_palette('bright',n_colors=len(cats_sc2) )
colorsd = dict( zip( cats_sh_sc2, colors ))
colormap = [ colorsd[x] for x in meta_sc2 ]

plt.figure(figsize=(20,20))
sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats_sh_sc2, size=1.0)
plt.savefig("blades-heatmap-legend-top2-sc2.png", dpi=600)


cats_sh = ['HPIV3 A549', 'IAV A549', 'IAV NHBE', 'IAVd NHBE', "LungHealthy b20","LungCOVID b20", "Mock A549", "Mock A549.1", "Mock Calu3", "Mock NHBE", 'RSV A549', 'SC2 A549', 'SC2 A549.1', 'SC2 A549.2', 'SC2 Calu-3', 'SC2 NHBE', 'Bowel', 'Fat', 'Heart', 'IFNB NHBE', 'Jejunum', 'Kidney', 'Liver', 'LungCOVID d20', 'LungHealthy d20', 'Marrow', 'Placenta', 'Skin']
```

Heatmap

``` python
plt.figure(figsize=(20,20))
sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx2_Z[topreg2], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig("blades-heatmap-top2.png", dpi=600, bbox_inches = "tight")

## sc2
plt.figure(figsize=(20,20))
sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_sc2_Z[topregsc2],  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=2, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig("blades-heatmap-top2-sc2.png", dpi=600, bbox_inches = "tight")
```

### **Refs**

  

In collaboration with Leonardo Arteaga, Ana Betty Villaseñor Altamirano,
Karen Nuñez Reza and Alejandra Medina Rivera.

  

**References:**

Datasets original publication:

\[1\] Blanco-Melo, D., Nilsson-Payant, B. E., Liu, W. C., Uhl, S.,
Hoagland, D., Møller, R., … & Albrecht, R. A. (2020). Imbalanced host
response to SARS-CoV-2 drives development of COVID-19. Cell, 181(5),
1036-1045.

\[2\] Desai, N., Neyaz, A., Szabolcs, A., Shih, A. R., Chen, J. H.,
Thapar, V., … & Deshpande, V. (2020). Temporal and spatial heterogeneity
of host response to SARS-CoV-2 pulmonary infection. Nature
communications, 11(1), 1-15.

\[3\] Delorey, T. M., Ziegler, C. G., Heimberg, G., Normand, R., Yang,
Y., Segerstolpe, Å., … & Regev, A. (2021). COVID-19 tissue atlases
reveal SARS-CoV-2 pathology and cellular targets. Nature, 1-8.

Methods:

-   loompy
-   pyscenic
-   vsn-pipelines multiruns
-   rss Revealing the Critical Regulators of Cell Identity in the Mouse
    Cell Atlas

  

**Programs and versions:**

-   R 4.0.2
-   python3.8
-   nextflow vsn-pipelines v0.25.0
-   singularity

packages:

-   tidyverse
-   loompy 3.0.2

  

*Last update: March 23rd, 2022*
