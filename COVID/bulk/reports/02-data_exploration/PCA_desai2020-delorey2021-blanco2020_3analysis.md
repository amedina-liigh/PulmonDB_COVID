PulmonDB - COVID19.
================
Monica Padilla

  

## **PCA of 3 parallel analysis**

`celllines`: Infected vs mock cell lines (Blanco 2020), `lung-covid`:
Lung COVID vs Healthy (Blanco 2020, Desai 2020, Delorey 2021) and
`organs-desgtex`: COVID tissues vs healthy (Desai 2020, GTEx v8)

Original sources: Delorey et al, 2021 \[1\], Desai et al. 2020 \[2\],
Blanco, GTEx Data preprocessing and matrixes for each experiment
generation is described in another report.

Explore overall transcriptomic profile of count matrixes that will be
used for further analysis.

### **Prepare terminal**

``` bash
ssh -Y mpadilla@dna.liigh.unam.mx
screen -S search
qlogin
cd /mnt/Citosina/amedina/mpadilla/COVID19/pca
module unload r/3.2.1
module load r/4.0.2
umask 2

mkdir plots_3analysis-bc
```

------------------------------------------------------------------------

  

### **Cell lines**

  

#### **All celllines experiments**

  

##### **Load and clean data**

``` r
library(tidyverse)
library(ggrepel)

## axes titles PC % explained var
pc.axis <- function(PCnum){
  paste0("PC",PCnum," (",
         summary(pca)$importance[2,][PCnum]*100,
         "% explained var.)")
}

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

# counts matrix - batch corrected
celllines <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.csv", sep = ","))
dim(celllines) # 74 20299

sub.celllines <- str_detect(meta$source, "Mock") | str_detect(meta$source, "infected") | str_detect(meta$source, "IFNB_treated_NHBE")
# subset for celllines samples
meta.sub <- meta[sub.celllines,]

## filter genes with zero counts and replace 0s for NAs
celllines[is.na(celllines)] <- 0
celllines <- celllines[,which(apply(celllines, 2, var) != 0)] #  74 19035
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(celllines),meta.sub$sample)
meta.sub <- meta.sub[order,]

# Create vector of experiment labels
exprmts <- rev(unique(meta.sub$source))
exprmts.short <- rev(c("hIFNB NHBE", "SC2-inf-NHBE", "SC2-inf-Calu3", "SC2-inf-A549 + ACE2 + pt",
                       "SC2-inf-A549 + ACE2", "SC2-inf-A549", "RSV-inf-A549", "Mock NHBE",
                       "Mock Calu3", "Mock A549 + ACE2", "Mock A549", "IAVdNS1-inf-NHBE",
                       "IAV-inf-NHBE", "IAV-inf-A549", "HPIV3-inf-A549"))
experiments <- meta.sub$source
for (i in 1:length(exprmts)) {
  experiments <- str_replace(experiments, exprmts[i], exprmts.short[i])
  i = i + 1
}
```

PCA computation and visualization

``` r
pca <- prcomp(celllines, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/celllines/pca_blanco2020_all-celllines-exp_PC56.png", width = 1000, height = 800)
ggplot(pca$x, aes(x=PC5, y=PC6, group=experiments)) +
  geom_point(size=4, aes(shape=experiments, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  #geom_mark_ellipse(aes(label=experiments,color=experiments), show.legend = F) +
  scale_shape_manual(values=c(0:14), name="Experiment") +
  scale_colour_discrete(name="Experiment") +
  labs(x = pc.axis(5), y = pc.axis(6)) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text( aes(label=experiments,color=experiments),
             nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T)
dev.off()
```

  

-   `Mock Calu-3` : Mock treated Calu-3 cell lines samples

-   `Mock NHBE` : Mock treated NHBE cell lines samples

-   `Mock A549` : Mock treated A459 cell lines samples

-   `Mock A549+ACE2` : Mock treated A459 cell lines with vector
    expression of ACE2 samples

-   `SC2-inf-A549` : SARS-CoV-2 infected A459 cell lines samples

-   `SC2-inf-A549+ACE2` : SARS-CoV-2 infected A459 cell lines with
    vector expression of ACE2 samples

-   `SC2-inf-A549+ACE2+pt` : SARS-CoV-2 infected A459 cell lines with
    vector expression of ACE2 with 1hr Ruxolitinib pre-treatment (IFN-I
    signaling blocked) samples

-   `SC2-inf-Calu3` : SARS-CoV-2 infected Calu-3 cell lines samples

-   `SC2-inf-NHBE` : SARS-CoV-2 infected NHBE cell lines samples

-   `RSV-inf-A549` : RSV infected A459 cell lines samples

-   `IAVdNS1-inf-NHBE` : IAVdNS1 infected NHBE cell lines samples

-   `IAV-inf-NHBE` : IAV infected NHBE cell lines samples

-   `IAV-inf-A549` : IAV infected A459 cell lines samples

-   `HPIV3-inf-A549` : HPIV3 infected A459 cell lines samples

-   `hIFNB NHBE` : NHBE cell lines samples treated with human
    Interferon-β (IFN-β)

-   PC1&2

![](imgs/pca/plots_3analysis-bc/celllines/pca_blanco2020_all-celllines-exp_PC12.png)

I think that what we’re seeing here is that the `PC1` is separating
experiments by expression of Interferon and its associated genes or the
lack of it. Because at the right is the experiment
`SARS-CoV-2-infected-A549 + ACE2 expression + Ruxolitinib pretreatment`
(purple stars) and `Mock NHBE` (green squares), whereas at the left is
the experiment `hINFB NHBE` which should be expressing and exacerbated
INF response in comparision to its control.

`SARS-CoV-2 A549 infected` samples wihtout pre-treatment are clustering
togetherat the uppper left position, joined by
`SARS-CoV-2 NHBE infected`, `hINFB NHBE` and close together to
`ACE2 expressing samples` is the `Mock A549 + ACE2`.

Infection to other viruses such as `HPIv3` and `IAV` seems to diverge a
lot from that of `SARS-CoV-2`, and have a high INF response following
the logic described above.

-   PC3&4

![](imgs/pca/plots_3analysis-bc/celllines/pca_blanco2020_all-celllines-exp_PC34.png)

Response to `IAV` and `RSV` infection seem to be very different from the
rest.

-   PC5&6

![](imgs/pca/plots_3analysis-bc/celllines/pca_blanco2020_all-celllines-exp_PC56.png)

-   PC7&8

![](imgs/pca/plots_3analysis-bc/celllines/pca_blanco2020_all-celllines-exp_PC78.png)

------------------------------------------------------------------------

  

#### **All celllines experiments except hIFNB treated NHBE**

  

##### **Load and clean data**

``` r
# counts matrix - batch corrected
celllines <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_noINFBNHBE-counts-fltrtrxs.csv", sep = ","))
dim(celllines) # 68 20299

# subset metadata
sub.celllines <- str_detect(meta$source, "Mock") | str_detect(meta$source, "infected")
# subset for celllines samples
meta.sub <- meta[sub.celllines,]

## filter genes with zero counts and replace 0s for NAs
celllines[is.na(celllines)] <- 0
celllines <- celllines[,which(apply(celllines, 2, var) != 0)] #  68 18982
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(celllines),meta.sub$sample)
meta.sub <- meta.sub[order,]

# Create vector of experiment labels
exprmts <- rev(unique(meta.sub$source))
exprmts.short <- rev(c("SC2-inf-NHBE", "SC2-inf-Calu3", "SC2-inf-A549 + ACE2 + pt",
                       "SC2-inf-A549 + ACE2", "SC2-inf-A549", "RSV-inf-A549", "Mock NHBE",
                       "Mock Calu3", "Mock A549 + ACE2", "Mock A549", "IAVdNS1-inf-NHBE",
                       "IAV-inf-NHBE", "IAV-inf-A549", "HPIV3-inf-A549"))
experiments <- meta.sub$source
for (i in 1:length(exprmts)) {
  experiments <- str_replace(experiments, exprmts[i], exprmts.short[i])
  i = i + 1
}
```

PCA computation and visualization

``` r
pca <- prcomp(celllines, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/celllines/pca_blanco2020_celllines-exp-noINFBNHBE_PC12.png", width = 1000, height = 800)
ggplot(pca$x, aes(x=PC1, y=PC2, group=experiments)) +
  geom_point(size=4, aes(shape=experiments, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(0:13), name="Experiment") +
  scale_colour_discrete(name="Experiment") +
  labs(x = pc.axis(1), y = pc.axis(2)) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text( aes(label=experiments, color=experiments), nudge_x = 8, nudge_y = 8, show.legend = F)
dev.off()
```

  

-   PC1&2
-   PC3&4
-   PC5&6
-   PC7&8

------------------------------------------------------------------------

  

### **Tissues COVID vs Healthy**

  

This experiment was causing a lot of bias in the exploratory analysis
where all count matrixes were merged together. In general, it was
expected since there where a lot of different sources. To minimize that
technical bias, the analysis was split as shown in this report. For this
experiment of **Tissues COVID vs Healthy** or `organs-desgtex`, there is
still expected to have a technical bias since we could only get a good
number of cases and control samples each from different sources, that
is, Desai COVID-19 samples and GTEx healthy samples. Batch correction
was performed jointly or separately taking into account dataset and
experiment technical variation sources.

  

#### **Desai and GTEx: Lung, Liver, Heart, Kidney, Bowel - Corrected separately**

  

##### **Load, join and clean data**

``` r
# counts matrix - batch corrected
gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-fltrensg.csv", sep = ","))
dim(gtex) #   50 13697

desai <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020/kallisto_gen/desai2020-counts-adj-selCovidOrgans_fltrtrxs.csv", sep = ","))
dim(desai) #  72 20299

# subset metadata
meta.sub <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/desai2020_gtexv8-subTissues-sample_source.csv", sep = ","))

## Join datasets
desai.g <- rownames_to_column(as.data.frame(t(desai)), var = "gene")
gtex.g <- rownames_to_column(as.data.frame(t(gtex)), var = "gene")
multiorg <- full_join(desai.g, gtex.g, by = "gene")
multiorg <- t(column_to_rownames(multiorg, var = "gene"))
dim(multiorg) #  122 20381

## filter genes with zero counts and replace 0s for NAs
multiorg[is.na(multiorg)] <- 0
multiorg <- multiorg[,which(apply(multiorg, 2, var) != 0)] #   122 19613
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(multiorg),meta.sub$sample)
meta.sub <- meta.sub[order,]

# Create vector of experiment labels
experiments <- meta.sub$source
dataset <- meta.sub$dataset
```

PCA computation and visualization

``` r
pca <- prcomp(multiorg, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)
# pca$x <- add_column(as.data.frame(pca$x), dataset)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_PC78-labhealthy.png", 
    width = 800, height = 800)
ggplot(pca$x, aes(x=PC7, y=PC8, group=experiments)) +
  geom_point(size=4, aes(shape=experiments, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(0:9), name="Tissue") +
  scale_colour_discrete(name="Tissue") +
  labs(x = pc.axis(7), y = pc.axis(8)) + #xlim(57.017922,57.017925) + ylim(3.9374615,3.937462) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text(aes(label=experiments,color=experiments), 
            nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T) +
  geom_label_repel(data = pca$x %>% select(PC7,PC8,experiments) %>%
                     slice(c(73,75,76,79,84)), aes(label=experiments,color=experiments),
             nudge_x = 55, nudge_y = 55, show.legend = F)
dev.off()
```

  

-   PC1&2

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_PC12-labhealthy.png)

PC1 is separating samples by source dataset.

-   PC3&4

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_PC34-labhealthy.png)

We can see now that samples from different datasets merge together,
however GTEx healthy samples are still clustered together at one point.

-   PC5&6

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_PC56-labhealthy.png)

-   PC7&8

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_PC78-labhealthy.png)

------------------------------------------------------------------------

  

#### **Desai and GTEx: Lung, Liver, Heart, Kidney, Bowel - Corrected together**

  

##### **Load, join and clean data**

``` r
# counts matrix - batch corrected
multiorg <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_countsfltrtrxs.csv", sep = ","))
dim(multiorgans) #   122 20299

# subset metadata
meta.sub <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/desai2020_gtexv8-subTissues-sample_source.csv", sep = ","))

## filter genes with zero counts and replace 0s for NAs
multiorg[is.na(multiorg)] <- 0
multiorg <- multiorg[,which(apply(multiorg, 2, var) != 0)] #   122 20299
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(multiorg),meta.sub$sample)
meta.sub <- meta.sub[order,]

# Create vector of experiment labels
experiments <- meta.sub$source
dataset <- meta.sub$dataset
```

PCA computation and visualization

``` r
pca <- prcomp(multiorg, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)
pca$x <- add_column(as.data.frame(pca$x), dataset)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_mer-adj_PC78-labhealthy.png", 
    width = 800, height = 800)
ggplot(pca$x, aes(x=PC7, y=PC8, group=experiments)) +
  geom_point(size=4, aes(shape=experiments, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(0:9), name="Tissue") +
  scale_colour_discrete(name="Tissue") +
  labs(x = pc.axis(7), y = pc.axis(8)) + #xlim(57.017922,57.017925) + ylim(3.9374615,3.937462) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text(aes(label=experiments,color=experiments), 
            nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T) +
  geom_label_repel(data = pca$x %>% dplyr::select(PC7,PC8,experiments) %>%
                     slice(c(73,75,76,79,84)), aes(label=experiments,color=experiments),
             nudge_x = 55, nudge_y = 55, show.legend = F)
dev.off()
```

  

-   PC1&2

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_mer-adj_PC12-labhealthy.png)

Performing batch correction with the datasets joined improved and
lessened the bias. As before, PC1 separates samples by source dataset,
but here `Healhty GTEx` samples are not clustered together at one point,
but have dispersed a little bit and there seems to be congruence
according to the tissue.

For COVID-19 samples, `Heart COVID` and `Bowel COVID` samples simulated
distributions are more constrained and less disperse than those of other
organs. `Liver COVID` sampels locate at the middle and `Lung COVID`
samples are the most disperse, these are also the most numerous.

-   PC3&4

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_mer-adj_PC34-labhealthy.png)

In PC3 and PC4, we see that samples from different datasets merge
together and again we can see that `healthy GTEx` samples have improved
in terms of dispersion.

As a conclusion, we will use this data matrix for further downstream
analysis.

-   PC5&6

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_mer-adj_PC56-labhealthy.png)

-   PC7&8

![](imgs/pca/plots_3analysis-bc/multiorgans/pca_desai2020-gtexv8_subTissues_mer-adj_PC78-labhealthy.png)

  

### **Lung COVID vs Healthy**

  

This dataset was also problematic since Blanco’s Lung samples are only 2
for cases and 2 for controls and Delorey’s data was taken from the
counts data matrix, so the pre-analysis wasn’t uniform with the other
datasets, also it has no controls. Desai’s samples are equilibrated.

To improve the congruence between datasets and lessen the technical
variance, we decided to also perform the **Lung COVID** analysis
seprated from the others.

#### **Metadata Lung samples**

``` r
## metadata bladesdel
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

## Subset metadata for lung samples
meta.lung <- str_detect(meta$source, "[Ll]ung")
meta.sub.lung <- meta[meta.lung,]

## Change tags across dataset to match biological condition
tissue <- rev(unique(meta.sub.lung$source))
bio.cond <- meta.sub.lung$source
tissue.covid <- c("Lung_COVID","Lung_healthy","Lung_COVID","Lung_COVID","Lung_healthy")
for (i in 1:length(tissue)) {
  bio.cond <- str_replace(bio.cond, tissue[i], tissue.covid[i])
  i = i + 1
}
meta.sub.lung$source <- bio.cond

## dataset 
meta.des <- str_detect(meta.sub.lung$sample, "^SAMN")
meta.del <- str_detect(meta.sub.lung$sample, "^\\d+")
dataset <- rep("Blanco_2020",length(meta.sub.lung$sample))
dataset[meta.des] <- "Desai_2020"
dataset[meta.del] <- "Delorey_2021"
meta.sub.lung[,"dataset"] <- dataset

## format delorey samples

## write out that file
write.table(x = meta.sub.lung, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/blanco2020-desai2020-delorey2021_allLungSamples_sample-source-dataset.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = F, col.names = T, quote = F)
```

  

#### **Blanco + Desai (batch correction separately)**

  

##### **Load, join and clean data**

``` r
# counts matrix - batch corrected
lung <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020_adj-desai2020_adj-Lung-counts-fltrtrxs.csv", sep = ","))
dim(lung) #    61 20299

# subset metadata
meta.lungdel <- meta.sub.lung$dataset == 3
# subset for lung samples from desai and blanco only
meta.sub <- meta.sub.lung[!meta.lungdel,]

## filter genes with zero counts and replace 0s for NAs
lung[is.na(lung)] <- 0
lung <- lung[,which(apply(lung, 2, var) != 0)] #   61 19618
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(lung),meta.sub.lung$sample)
meta.sub.lung <- meta.sub.lung[order,]

# Create vector of experiment labels
experiments <- meta.sub.lung$source
dataset <- as.factor(meta.sub.lung$dataset)
```

PCA computation and visualization

``` r
library(ghibli)
pca <- prcomp(lung, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)
pca$x <- add_column(as.data.frame(pca$x), as.factor(dataset))

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/lung/pca_blanco2020_adj-desai2020_adj-Lung-PC78.png",
    width = 800, height = 800)
ggplot(pca$x, aes(x=PC7, y=PC8, group=experiments)) +
  geom_point(size=4, aes(shape=dataset, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(15,1), name="Dataset") +
  scale_colour_manual(values=c("#f38989","#7289da"), name="Tissue") +
  labs(x = pc.axis(7), y = pc.axis(8)) + 
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text(aes(label=experiments,color=experiments),
            nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T)
dev.off()
```

  

-   PC1&2

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_adj-desai2020_adj-Lung-PC12.png)

-   PC3&4

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_adj-desai2020_adj-Lung-PC34.png)

-   PC5&6

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_adj-desai2020_adj-Lung-PC56.png)

-   PC7&8

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_adj-desai2020_adj-Lung-PC78.png)

------------------------------------------------------------------------

  

#### **Blanco + Desai (batch correction joined)**

  

##### **Load, join and clean data**

``` r
# counts matrix - batch corrected
lung <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-Lung-fltrtrxs.csv", sep = ","))
dim(lung) #    61 20299

# subset metadata
meta.lungdel <- str_detect(meta.sub.lung$dataset, "Delorey_2021")
# subset for lung samples from desai and blanco only
meta.sub <- meta.sub.lung[!meta.lungdel,]

## filter genes with zero counts and replace 0s for NAs
lung[is.na(lung)] <- 0
lung <- lung[,which(apply(lung, 2, var) != 0)] #   61 19618
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(lung),meta.sub.lung$sample)
meta.sub.lung <- meta.sub.lung[order,]

# Create vector of experiment labels
experiments <- meta.sub.lung$source
dataset <- as.factor(meta.sub.lung$dataset)
```

PCA computation and visualization

``` r
pca <- prcomp(lung, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)
pca$x <- add_column(as.data.frame(pca$x), dataset)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/lung/pca_blanco2020_desai2020-mer_adj-Lung-PC56.png",
    width = 800, height = 800)
ggplot(pca$x, aes(x=PC5, y=PC6, group=experiments)) +
  geom_point(size=4, aes(shape=dataset, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(15,1), name="Dataset") +
  scale_colour_manual(values=c("#f38989","#7289da"), name="Tissue") +
  labs(x = pc.axis(5), y = pc.axis(6)) + 
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text(aes(label=experiments,color=experiments),
            nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T)
dev.off()
```

  

-   PC1&2

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_desai2020-mer_adj-Lung-PC12.png)

The good news is that we can see two groups: healthy samples and COVID
samples. Some COVID samples merge together, whereas others drive
variation in PC1. However, there is no clear separation between datasets
in a condition wise manner.

-   PC3&4

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_desai2020-mer_adj-Lung-PC34.png)

-   PC5&6

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_desai2020-mer_adj-Lung-PC56.png)

Here we can finally see healthy and COVID samples separated regardless
of source dataset.

-   PC7&8

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020_desai2020-mer_adj-Lung-PC78.png)

------------------------------------------------------------------------

  

#### **Blanco + Desai + Delorey (batch correction joined)**

  

##### **Load, join and clean data**

``` r
# counts matrix - batch corrected
lung <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs.csv", sep = ","))
dim(lung) #   78 20299

## format delorey samples names
rownames(lung) <- str_remove(rownames(lung),"^X")
rownames(lung) <- str_replace_all(rownames(lung),"\\.","-")

## filter genes with zero counts and replace 0s for NAs
lung[is.na(lung)] <- 0
lung <- lung[,which(apply(lung, 2, var) != 0)] #    78 19919
```

##### **PCA**

Label vectors:

``` r
# order meta
order <- match(rownames(lung),meta.sub.lung$sample)
meta.sub.lung <- meta.sub.lung[order,]

# Create vector of experiment labels
experiments <- meta.sub.lung$source
dataset <- as.factor(meta.sub.lung$dataset)
```

PCA computation and visualization

``` r
pca <- prcomp(lung, scale = T)

## add labels to pca obj
pca$x <- add_column(as.data.frame(pca$x), experiments)
pca$x <- add_column(as.data.frame(pca$x), dataset)

## ggplot2 col-SC2inf  shape -experiments
png(file = "./plots_3analysis-bc/lung/pca_blanco2020-desai2020-delorey2021_mer-adj_Lung_PC34.png",
    width = 800, height = 800)
ggplot(pca$x, aes(x=PC3, y=PC4, group=experiments)) +
  geom_point(size=4, aes(shape=dataset, color=experiments)) +
  stat_ellipse(aes(color=experiments), show.legend = F) +
  scale_shape_manual(values=c(15,1,9), name="Dataset") +
  scale_colour_manual(values=c("#f38989","#7289da"), name="Tissue") +
  labs(x = pc.axis(3), y = pc.axis(4)) + 
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        axis.title.x=element_text(face = "bold"), axis.title.y=element_text(face = "bold")) +
  geom_text(aes(label=experiments,color=experiments),
            nudge_x = 8, nudge_y = 8, show.legend = F, check_overlap = T)
dev.off()
```

  

-   PC1&2

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020-desai2020-delorey2021_mer-adj_Lung_PC12.png)

I think PC1 is mostly separating Delorey and Desai+Blancos’s datasets.
The trend we saw before in the PCA’s of Desai + Blancos data is surely
repeating itself here but at the right part of this plot. We see the
group of healthy samples together.

Delorey’s data is very disperse. But the definition of groups has
definitely improved in comparision from when the analysis wasnt split.

-   PC3&4

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020-desai2020-delorey2021_mer-adj_Lung_PC34.png)

There is still the group of `Healthy Lung` samples. Also, `Lung COVID`
samples from different datasets are merging together, but with still a
bias.

-   PC5&6

![](imgs/pca/plots_3analysis-bc/lung/pca_blanco2020-desai2020-delorey2021_mer-adj_Lung_PC56.png)

------------------------------------------------------------------------
