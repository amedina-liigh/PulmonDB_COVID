PulmonDB - COVID19.
================
Monica Padilla

  

## **Data pre-processing additional datasets: Desai2020 (Lung + other organs), Delorey2020 (Lung COVID) and GTEx (some tissues)**

Original sources: Delorey et al, 2021 \[1\], and Desai et al. 2020
\[2\]. In this report is also reported a re-run of experiment count
matrix generation and batch correction of blanco’s data.

### **Prepare terminal**

``` bash
ssh -Y mpadilla@dna.liigh.unam.mx
screen -S search
qlogin
cd /mnt/Citosina/amedina/mpadilla/COVID19/
module unload r/3.2.1
module load r/4.0.2
module load e-utilities/27abr20
umask 2
#mkdir deloney2021
cd Data/deloney2021
```

------------------------------------------------------------------------

  

### **Search and Retrieve data**

#### **Search Data**

We need RNA-seq Lung samples with controls and infected with SARS-CoV-2
or from patients with COVID-19.

at web geo datasets:
`(rna-seq[All Fields] AND ("lung"[MeSH Terms] OR lung[All Fields]) AND (covid[All Fields] OR ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR sars-cov-2[All Fields]))) NOT scRNA-seq[All Fields]) AND "Homo sapiens"[orgn]`

``` bash
einfo -dbs # view dbs, gds = geo data sets
einfo -db gds -fields # search db specific fields
esearch -db gds -query "rna-seq[ALL] AND (\"lung\"[MESH] OR lung[ALL]) AND (covid[ALL] OR \"Severe acute respiratory syndrome coronavirus 2\"[ALL] OR sars-cov-2[ALL])) AND \"Homo sapiens\"[ORGN]" | efetch -format uid | esummary -db gds | perl -pe 's/\n/,/g'


esearch -db gds -query "rna-seq[ALL] AND (covid[ALL] OR \"Severe acute respiratory syndrome coronavirus 2\"[ALL] OR sars-cov-2[ALL])) AND \"Homo sapiens\"[ORGN] NOT (scRNA-seq[ALL] OR \"single-cell RNA-seq\"[ALL])"
```

\*resulting UIDs are shown below but were also saved at
`uid_gds_lungSC2rnaseq.txt`

``` bash
esummary -db gds -id 200176393,200167403,200167402,200176201,200182297,200181827,200179184,200157490,200178824,200151513,200164013,200163529,200162911,200171668,200155223,200167131,200163547,200161934,200157344,200161262,200161382,200155518,200160351,200158127,200152586,200159191,200155286,200156988,200154613,200151161,200147903,200150819,200150316,200148829,305471867,305471866,305471865,305471864,305471863,305471862,305471861,305471860,305471859,304796155,304796154,304796153,304796152,304796119,304796118,304796117,304796116,304796115,304796114,304796113,304766939,304766938,305230128,305230124,305230121,305230115,305230114,305230059,305230057,305230044,305230040,305230037,305230031,305230030,305229975,305229973,305230112,305230111,305230110,305230107,305230106,305230098,305230092,305230091,305230090,305230087,305230085,305230080,305230075,305230068,305230066,305230064,305230063,305230028,305230027,305230026,305230023,305230022,305230014,305230008,305230007,305230006,305230003,305230001,305229996,305229991,305229984,305229982,305229980,305229979,304822344,304822347,304822346,304822345,304698499,304698498,304698497,304698496,304698495,304698494,304698493,304698492,304698491,304698490,304698489,304698488,304698487 > esummary_gds_lungSC2rnaseq.xml
```

Saved file located at =
`mpadilla/COVID19/Data/deloney2021/metadata/esummary_gds_lungSC2rnaseq.xml`

  

  

------------------------------------------------------------------------

  

#### **Data selected**

-   [GSE171668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171668),
    that comes from this publication [COVID-19 tissue atlases reveal
    SARS-CoV-2 pathology and cellular
    targets](https://www.nature.com/articles/s41586-021-03570-8) by Toni
    M Delorey, …, Aviv Regev; April 2021.

-   [GSE150316](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150316),
    also has raw data available:
    [SRP261138](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP261138),
    PRJNA631753. That comes from this publication [Temporal and spatial
    heterogeneity of host response to SARS-CoV-2 pulmonary
    infection.](https://www.nature.com/articles/s41467-020-20139-7) by
    Desai N et al, 2020.

------------------------------------------------------------------------

  

#### **Data retrieval:GSE171668 - delorey2021**

Changed r version in cluster:

``` bash
module load r/3.6.1
```

Use GEOquery:

``` r
library(GEOquery)
gse <- getGEO("GSE171668", destdir = ".")
# supp <- getGEOSuppFiles("GSE171668", fetch_files = TRUE, baseDir = ".") # i didnt run this bc i didnt want all suppfiles
```

downloaded files are in = `delorey2021/GEOquery/`

  

Get bulk counts and metadata:

``` bash
mkdir counts
mkdir metadata
cd metadata
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171668/suppl/GSE171668_bulk_metadata.csv.gz
cd ../counts
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171668/suppl/GSE171668_rsem.genes.counts.matrix.txt.gz
```

Checking that samples correspond to bulk samples:

``` bash
# next file has the sample name of the counts matrix
head -n1 GSE171668_rsem.genes.counts.matrix.txt | perl -pe 's/\s/\n/g' | sort | sed 1d > GSE171668_rsem_samples.txt
# next file has bulk samples names
cut -d, -f1 ../metadata/GSE171668_bulk_metadata.csv | sort | head -n20 | perl -pe 's/-/_/g' > bulk_samples.txt
# report identical files
diff -s GSE171668_rsem_samples.txt bulk_samples.txt
```

`Files GSE171668_rsem_samples.txt and bulk_samples.txt are identical`

------------------------------------------------------------------------

  

#### **Data retrieval: SRP261138 - desai2020**

  

Metadata was manually downloaded from [SRA Run
Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA631753&o=acc_s%3Aa)
and lives here:
`mpadilla/COVID19/Data/desai2020/metadata/SraRunTable_SRP261138_desai2020.txt`

  

The below `fastq-dump` arguments were selected with guidance of this
[manual](https://edwards.sdsu.edu/research/fastq-dump/).

`fastq-dump` automatically stores cache in `~` dir, which ended up
consuming all my home space and I got the `disk quota exceeded`. To
disable this option I did the following:

``` bash
qlogin
module load sra/2.9.6-1
vdb-config -i
```

That opens an interface, press `2` (disable local file caching) &gt; `6`
(save) &gt; `7` (exit).

  

Code below was run in sge file
`mpadilla/COVID19/Data/desai2020/fastq/fastqdump.sge`.

``` bash
module load e-utilities/27abr20
module load sra/2.9.6-1

# Get SRR accessions, then download fastq files
esearch -db sra -query "SRP261138" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" | xargs fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files
```

job log file = `less fastqdump.o192951`

SOme files were missing, run with `fastqdump2.sge`

``` bash
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772358
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772359
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772360
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772361
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772362
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772363
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772364
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772365
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772366
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR11772367
fastq-dump --outdir . --gzip --skip-technical --read-filter pass --dumpbase --split-files SRR12340064
```

------------------------------------------------------------------------

  

#### **Data retrieval: GTEx v8 RNA-seq**

  

Datasets retrieved from: <https://gtexportal.org/home/datasets>

##### **RNA-seq transcript TPM counts**

Gene quantification done with RSEM, as reported here:
<https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods>

``` bash
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
```

location:
`/mnt/Citosina/amedina/mpadilla/resources/gtex/data_v8/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz`

  

##### **Sample annotation**

``` bash
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
```

location:
`/mnt/Citosina/amedina/mpadilla/resources/gtex/data_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`

  

##### **GENCODE v26 annotation**

``` bash
wget https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf
```

------------------------------------------------------------------------

  

### **Data pre-processing**

#### **Delorey 2021**

  

##### **Batch note**

  

Delorey’s dataset bulk RNA-seq samples all have the same batch,
i.e. samples were prepared the same way, the same library preparation
and sequenced in the same day. Thus, no within dataset correction is
necessary.

##### **Normalization to CPMs**

  

**Count per Million CPM**

For each sample:

1.  Sum counts for all genes
2.  Divide genes\_counts by sum
3.  Multiply genes by 10^6

``` r
# Load data of delorey 2021
del <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/counts/GSE171668_rsem.genes.counts.matrix.txt", sep = "\t"))
dim(del) # 23686    20

del <- as.matrix(del)

# 1
sum.counts.samp <- apply(del,2,sum)
del.cpm <- t((t(del)/sum.counts.samp) * 10^6) # it is the same
```

------------------------------------------------------------------------

  

##### **Filter samples for Lung**

  

**Select relevant samples : Lung COVID-19 biopsies**

Taken and adapted from:
<https://github.com/klarman-cell-observatory/covid19-autopsy/blob/main/LungBulk/deconvolution_of_lung_bulk_samples.R>

``` r
meta.del <- read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/metadata/GSE171668_bulk_metadata-sub.csv", sep = ",")

# Samples that were present in the input data but are not part of bulk lung samples
#                     D6:brain              D3:heart              D10:Lung                  
exclude.samples <- c("X02_P334354_S106_R01","X02_P166169_S027_R01","X04_P054921_S029_R01")
# D3 and D6 have 3 lung samples c("02-P166169-S038-R01","02-P334354-S088-R01","02-P334354-S066-R01")
# but these were not excluded in their deconvolution analysis as shown in the link
exclude.samples.donor <- c("D6","D3")

# Remove samples not used here from the bulk
del <- del[,!(colnames(del) %in% exclude.samples)]
dim(del) #  23686    17

# Write file with only bulk lung samples
write.table(x = del, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/counts/GSE171668-Lung_rsem.genes.cpm_counts.matrix.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

------------------------------------------------------------------------

  

##### **Convert gene\_names to ENST**

This step is necessary to perform batch correction between datasets with
a different gene definition and notation.

Convert gene names of Delorey to ENST ids, where the ENST ids to
synonyms are the ones that were chosen as “best” previously, with
Desai’s and Blanco’s transcripts, since these are expected to be ones
assesed by other alignment programs.

``` r
# Load data of delorey 2021
del <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/counts/GSE171668-Lung_rsem.genes.cpm_counts.matrix.csv", sep = ","))

## 1. load table of conversion between enstids, ensgids and gene names after filtering for best trxs
# load script with matrix summarizing functions
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")

## Table of conversion

# ensembl_ids that result from running kallisto with gencode transcripts
ensembl_ids.kall <- read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/ensembl_ids.txt", header = F, sep = "|") # 106,143 transcripts
# gencode adds a series of annotations
colnames(ensembl_ids.kall) <- c("ensembl_transcript_id_version","ensembl_gene_id_version",
                                "havana_gene_id", "havana_transcript_id",
                                "transcript_name","gene_name",c(7:11))

# conversion table of ENST ENSG and gene names from Desai's count matrix
res.des <- Id2GeneName(ensembl_ids.kall$ensembl_transcript_id_version, 
                   filter_transcripts = T, gencode_filter = T, get_gene_ids = T)

## 2. Check which genes of delorey's notation also appear in kallisto run
genes2change <- res.des$external_gene_name %in% rownames(del)
res.des <- res.des[genes2change,]
order <- match(rownames(del),res.des$external_gene_name)
res.des <- res.des[order,]

#  change their names to enstids
genes2changedel <- match(na.omit(res.des$external_gene_name), rownames(del))
del.enst <- del[genes2changedel,]
del.genename <- setdiff(del, del.enst)
rownames(del.enst) <- res.des$target_id[genes2changedel]

## for those not found, just keep the gene name. conversion of common genes was necessary for accurate comparison
del.all <- union(del.enst,del.genename)

## save matrix with conversion of gene names to enstids
write.table(x = del.all, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/counts/GSE171668-Lung_rsem.enst-genes.cpm_counts.matrix.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

------------------------------------------------------------------------

  

#### **Desai 2020**

  

Repeat the whole data processing done with Blanco’s dataset.

  

##### **Make Fastqcs**

``` bash
umask 2
mkdir fastqc
ls fastq/ > fastqc/makefastqc.sge
cd fastqc
perl -pe 's/^/fastqc ..\//' makefastqc.sge | perl -pe 's/fastq\.gz/fastq\.gz -f fastq -o . -d \.\/tmp/' > fastqc.sge
qsub fastqc.sge
```

------------------------------------------------------------------------

  

##### **Trimming**

``` bash
mkdir trimmed
module load trimmomatic/0.33
for file in ../fastq/*1.fastq.gz; do trimmomatic PE -phred33 -basein $file -baseout ${file//_1.fastq/_trmd_1.fastq} ILLUMINACLIP:/mnt/Citosina/amedina/mpadilla/resources/trimmomatic/adapters/TruSeq-SE.fa:2:30:10 SLIDINGWINDOW:5:30 MINLEN:40; done
```

Order reads

``` bash
mv fastq/*trmd* trimmed/
cd trimmed
cut -d, -f1,2 ../metadata/subset_SRP261138_desai2020.csv | perl -pe 's/(SRR.*),(SAMN.*)/mkdir $2; mv $1\* $2\//'
```

  

------------------------------------------------------------------------

##### **kallisto + gencode**

Download [GENCODE v38](https://www.gencodegenes.org/human/) ref
transcriptome :

``` bash
### reference Protein-coding transcript sequences
### gencode.v38.pc_transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.pc_transcripts.fa.gz
```

Annotation references:

-   EMBL
-   CCDS : [Consensus Coding Sequence
    Project](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi)
-   RefSeq\_dna, \_peptide
-   Uniprot
-   SWISSPROT

  

Make gencode transcriptome index for kalisto :

``` bash
### index_kallisto45_gencode38_hsapiens_cds
kallisto index -i ./kallisto/index_kallisto45_gencode38_hsapiens_cds gencode/v38/gencode.v38.pc_transcripts.fa.gz
```

  

To make sge file with `kallisto quant` command, I did:

dir =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020/kallisto_gen/`

``` bash
cut -d, -f1,2 ../metadata/subset_SRP261138_desai2020.csv | perl -pe 's/.*,(SAMN.*)/mkdir $1; kallisto quant -i \/mnt\/Citosina\/amedina\/mpadilla\/resources\/kallisto\/index_kallisto45_gencode38_hsapiens_cds -o $1 ..\/trimmed\/$1\/\*P.fastq.gz/' > kquant_gen.sge
```

The pipeline below shows how `kallisto quant` was run.

``` bash
mkdir kallisto
module load kallisto/0.45.0
mkdir sample
kallisto quant -i /mnt/Citosina/amedina/mpadilla/resources/kallisto/index_kallisto45_ensemblv102_hsapiens_cdna -o sample ../trimmed/sample/*P.fastq.gz
```

  

`abundance.tsv` files contain counts for `106143` (protein coding)
transcripts.

------------------------------------------------------------------------

  

##### **Gene annotation gencodev38**

To reduce the number of transcripts to approx 1 per gene, we used
GENCODE’s gene basic annotation for protein coding transcripts.

The table with ensembl ENSG ids, ENT ids and corresponding gene names
are located in =
`/mnt/Citosina/amedina/mpadilla/resources/gencode/v38/fltr_transcripts/gencode.v38.basic.pc.transcripts.ensgid_enstid_genename.csv`

This file has `61067` transcripts.

Some genes still have many transcripts as shown below.

``` r
counts <- read.csv("transcript_enstids_counts.csv", header = T, sep = ",")
mt30 <- which(counts$count_enst_ids > 30) # more than 30 transcripts
png(file = "./transcript_enstids_30counts.png", width = 1000, height = 400)
counts %>% 
  mutate(name = fct_reorder(gene_name, count_enst_ids)) %>%
  slice(1:22) %>%
  ggplot(aes(x=name, y=count_enst_ids)) +
  geom_col(fill="#825CA6") +
  labs(x="", y="# transcripts ids")
dev.off()
```

![](imgs/transcript_enstids_30counts.png)

For the `13626` genes that have more than one transcript assigned, we
try to select for the best one with trasncript flags as before.

------------------------------------------------------------------------

  

##### **TPM Count matrix all transcripts, Blanco - Desai**

TPMs, also dont filter genes:

For Blanco’s data:

``` r
library(tidyverse)

# get file paths of count matrixes from blanco2020 and desai2020
bla.samples.dirs <- str_subset(str_subset(list.dirs("."), "counts3"), "/.+/.+/.+/.+")
bla.samples.dirs <- str_subset(bla.samples.dirs, "Ferret", negate=T)
bla.files <- paste0(bla.samples.dirs, "/abundance.tsv")

# regex to identy sample names
sample_regex <- "[SAMNG]{3,4}\\d+"

# Join all samples counts into one matrix
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R") # code make experiment count matrix

# blanco
MakeExpCountsMatrix(files = bla.files, outfile_name = "blanco2020/counts3/blanco2020-counts_gencv38-alltrxs.csv", 
                    sample_regex = sample_regex, filter_transcripts = F, 
                    host_ensembl = "https://may2021.archive.ensembl.org",
                    gencode_filter = F, output_matrix_transcripts_ids = T, 
                    use_gencode = T)
```

  

For Desai’s data:

``` r
des.samples.dirs <- str_subset(str_subset(list.dirs("."), "kallisto_gen"), "/.+/.+/.+")
des.files <- paste0(des.samples.dirs, "/abundance.tsv")

# desai
MakeExpCountsMatrix(files = des.files, outfile_name = "desai2020/kallisto_genc/desai2020-counts-alltrxs.csv", 
                    sample_regex = sample_regex, filter_transcripts = F, 
                    host_ensembl = "https://may2021.archive.ensembl.org",
                    gencode_filter = F, output_matrix_transcripts_ids = T, 
                    use_gencode = T)
```

  

------------------------------------------------------------------------

##### **Batch Correction on Blanco - Desai TPM matrix all trxs**

  

**Define batches and run combat-seq**

###### **Per Dataset (dates and platform)**

  

**For Blanco’s data**

``` r
# Load blanco's experiment count matrix
blanco <- as.data.frame(read.csv("blanco2020/counts3/blanco2020-counts_gencv38-alltrxs_trxs.csv", sep = ","))
dim(blanco) #    78 106143
# Load metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/date_sample_source-Human.tsv", sep = "\t"))
colnames(meta) <- c("date","sample","source")

# Filter for Blanco's data
meta.bla <- str_detect(meta$sample, "^GSM")
meta.bla <- meta[meta.bla,]

# Put metadata in the same order
index <- match(rownames(blanco), meta.bla$sample)
meta.bla <- meta.bla[index,]

# Define batches (in this case by date)
dates <- unique(meta.bla$date)
batch <- meta.bla$date
for (i in 1:length(dates)) {
  batch <- str_replace(batch, dates[i], paste0(i))
  i = i + 1
}
batch <- as.numeric(batch)

# Define bio condition (co variate), to preserve signal in adj data
exp <- rev(unique(meta.bla$source)) # rev order so that it doesnt match other experiments
bio.cond <- meta.bla$source
for (i in 1:length(exp)) {
  bio.cond <- str_replace(bio.cond, exp[i], paste0(i))
  i = i + 1
}
bio.cond <- as.numeric(bio.cond)

# transpose and turn to matrix
blanco <- t(as.matrix(blanco)) # combat seq requires rows genes

# Correction of batch effects with ComBat-seq
library(sva)
adj.bla <- ComBat_seq(blanco, batch = batch, group = bio.cond)
```

Check created object

``` r
class(adj.bla) # "matrix" "array"
# adj.bla <- adj.bla[which(apply(adj.bla,1,var) != 0),] # quit counts that transform to 0
dim(adj.bla) # 106143     78
object.size(adj.bla)/1e6 # 75.6 MB
```

*Note:* Object saved as and can be loaded as following.

``` r
# write.csv(adj.bla, file = "blanco2020/counts3/blanco2020-counts_gencv38-alltrxs_trxs-adj.csv", quote = F)
```

  

**For Desai’s data**

``` r
# Load desai's experiment count matrix
desai <- as.data.frame(read.csv("desai2020/kallisto_gen/desai2020-counts-alltrxs_trxs.csv", sep = ","))
dim(desai)
# Load metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/date_sample_source-Human.tsv", sep = "\t"))
colnames(meta) <- c("date","sample","source")

# Filter for Blanco's data
meta.des <- str_detect(meta$sample, "^SAMN")
meta.des <- meta[meta.des,]

# Put metadata in the same order
index <- match(rownames(desai), meta.des$sample)
meta.des <- meta.des[index,]

# tag lung controls of desai 2020
lung.control.d20 <- meta.des$sample %in% c("SAMN14891157","SAMN14891156","SAMN14891155","SAMN14891154","SAMN14891153")
meta.des$source[lung.control.d20] <- "lung_healthy_d20"

# Define batches (in this case by date)
dates <- unique(meta.des$date)
batch <- meta.des$date
for (i in 1:length(dates)) {
  batch <- str_replace(batch, dates[i], paste0(i))
  i = i + 1
}
batch <- as.numeric(batch)

# there's also a batch taking the platform used into account (all other samples were done using nextseq 500)
illumina.miseq <- c("SAMN15660687","SAMN15660692","SAMN15660698")
miseq.samples <- which(rownames(desai) %in% illumina.miseq)
batch[miseq.samples] <- 3 # add batch
  
# Define bio condition (co variate), to preserve signal in adj data
exp <- rev(unique(meta.des$source)) # rev order so that it doesnt match other experiments
bio.cond <- meta.des$source
for (i in 1:length(exp)) {
  bio.cond <- str_replace(bio.cond, exp[i], paste0(i))
  i = i + 1
}
bio.cond <- as.numeric(bio.cond)

# transpose and turn to matrix
desai <- t(as.matrix(desai)) # combat seq requires rows genes

# Correction of batch effects with ComBat-seq
library(sva)
adj.des <- ComBat_seq(desai, batch = batch, group = bio.cond)
```

    Found 3 batches
    Using full model in ComBat-seq.
    Adjusting for 10 covariate(s) or covariate level(s)
    Estimating dispersions
    Fitting the GLM model
    Shrinkage off - using GLM estimates for parameters
    Adjusting the data

Check created object

``` r
class(adj.des) # "matrix" "array"
# adj.des <- adj.des[which(apply(adj.des,1,var) != 0),] # quit counts that transform to 0
dim(adj.des) # 106143     88
object.size(adj.des)/1e6 # 84.1 MB
```

*Note:* Object saved as and can be loaded as following.

``` r
# write.csv(adj.des, file = "desai2020/kallisto_gen/desai2020-counts-alltrxs-adj.csv", quote = F)
```

###### **Second Correction: between datasets**

  

``` r
## adj.bla, adj.des
## Join adjusted counts from each
adj.bla.g <- rownames_to_column(as.data.frame(adj.bla), var = "trx_id")
adj.des.g <- rownames_to_column(as.data.frame(adj.des), var = "trx_id")
adj.both <- full_join(adj.bla.g, adj.des.g, by = "trx_id")
adj.both <- column_to_rownames(adj.both, var = "trx_id")
adj.both <- as.matrix(adj.both)

# subset for blades samples
meta.bla <- str_detect(meta$sample, "^GSM")
meta.des <- str_detect(meta$sample, "^SAMN")
meta <- meta[meta.bla | meta.des,]

# Put metadata in the same order
index <- match(colnames(adj.both), meta$sample)
meta <- meta[index,]

# Define batches (in this case by origin of dataset)
samp.bla <- which(str_detect(meta$sample, "^GSM"))
samp.des <- which(str_detect(meta$sample, "^SAMN"))
batch <- rep(0,length(meta$sample))
batch[samp.bla] <- 1 # blanco batch
batch[samp.des] <- 2 # desai batch

# tag lung controls of desai 2020
lung.control.d20 <- meta$sample %in% c("SAMN14891157","SAMN14891156","SAMN14891155","SAMN14891154","SAMN14891153")
meta$source[lung.control.d20] <- "Lung_Healthy"
# tag equally lung covid samples across datasets
lung.covid <- c("lung","Lung_sample_from_postmortem_COVID-19_patient")
lung.covid.samples <- which(meta$source %in% lung.covid)
meta$source[lung.covid.samples] <- "Lung_COVID"
# tag equally lung healthy samples across datasets
lung.control <- c("Lung_Healthy","Lung_biopsy_for_heatly_negative_control")
lung.control.samples <- which(meta$source %in% lung.control)
meta$source[lung.control.samples] <- "Lung_Healthy"

# Define bio condition (co variate), to preserve signal in adj data
exp <- rev(unique(meta$source)) # rev order so that it doesnt match other experiments
bio.cond <- meta$source
for (i in 1:length(exp)) {
  bio.cond <- str_replace(bio.cond, exp[i], paste0(i))
  i = i + 1
}
bio.cond <- as.numeric(bio.cond)

# Correction of batch effects with ComBat-seq
library(sva)
adj.both <- ComBat_seq(adj.both, batch = batch, group = bio.cond)
```

    Found 2 batches
    Using full model in ComBat-seq.
    Adjusting for 25 covariate(s) or covariate level(s)
    Estimating dispersions

Check created object

``` r
class(adj.both) # "matrix" "array"
# adj.des <- adj.des[which(apply(adj.des,1,var) != 0),] # quit counts that transform to 0
dim(adj.both) # 106143    166
object.size(adj.both)/1e6 # 150.9 MB
```

*Note:* Object saved as and can be loaded as following.

``` r
# write.csv(adj.both, file = "multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-alltrxs-adj.csv", quote = F)
```

------------------------------------------------------------------------

#### **GTEx v8 RNA-seq**

  

##### **Subset for data of interest**

Applied filters and criteria:

-   RNA-seq data

-   Tissues of interest:

-   Bowel (intestine) = `Small Intestine - Terminal Ileum`

-   Kidney = `Kidney - Cortex`, !`Kidney - Medulla` because onlu miRNA
    data

-   Liver = `Liver`

-   Lung = `Lung`

-   Heart = `Heart - Left Ventricle`, `Heart - Atrial Appendage`

-   Batch repetition within and across tissues, to correct for that
    effect

-   Good reported quality

  

Filter sample annotation by experiment (rnaseq) and tissues:

``` r
library(dplyr)

# read sample annotation file
meta <- as.data.frame(read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep = "\t"))

# filter by experiment type of interest (mrna rnaseq)
#        and by tissues of interest
meta <- meta %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMNABTCHT , SMGEBTCH, SMGEBTCHT, SMAFRZE) %>%
  filter(SMNABTCHT == "RNA Extraction from Paxgene-derived Lysate Plate Based", 
         SMAFRZE == "RNASEQ", SMGEBTCHT == "TruSeq.v1",
         SMTS == "Heart" | SMTS == "Lung" | SMTS == "Small Intestine" | SMTS == "Liver"
         | SMTS == "Kidney" ) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH)
```

  

##### **Check which batches to use**

  

We want \~10 samples per tissue. Criteria for choosing them:

-   Only `Heart` samples are split into two categories:
    `Atrial appendage` and `Left Ventricle`, choose same number of
    samples for each.
-   Have the minimum number of batches, provided they are repeated
    within and across tissues.

``` r
# repetition of batches
meta %>%
  group_by(SMNABTCH, SMGEBTCH) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count))

# repetition of tissues and batches
batch_rep_perts <- meta %>%
  group_by(SMTS, SMNABTCH, SMGEBTCH) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count > 1)

# repetition of batches between diff tissues
batch_rep_manyts <- meta %>%
  group_by(SMTS, SMNABTCH, SMGEBTCH) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count > 1) %>% 
  group_by(SMNABTCH, SMGEBTCH) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count > 1)

# batches repeated across 3 or 4 organs
batch_rep_ts34 <- batch_rep_manyts$SMNABTCH[1:9]
# batches repeated across 2 organs
batch_rep_ts2 <- batch_rep_manyts$SMNABTCH[-c(1:9)]
```

Batches that are most repeated across diff organs (counts of organs
present shown):

    SMNABTCH SMGEBTCH   Count  +selected
       <chr>    <chr>      <int>
     1 BP-48576 LCSET-4955     4 +
     2 BP-75925 LCSET-9767     4 +
     3 BP-43881 LCSET-4793     3
     4 BP-44261 LCSET-4796     3 +
     5 BP-45390 LCSET-4633     3
     6 BP-47696 LCSET-4907     3 +
     7 BP-67833 LCSET-8132     3 +
     8 BP-69094 LCSET-8723     3
     9 BP-73249 LCSET-9410     3
    10 BP-42087 LCSET-4038     2

Check for batches repeated across 4 or 3 tissues to start selection. If
not all tissues are present in those samples (the case for Kidney and
Intestine), take batches repeated across 2 tissues to complete number
wanted of samples.

``` r
# limiting organ is Kidney and Small Intestine, take whichever batches that are repeated in each
# Kidney
batch_rep_perts %>%
  filter(SMNABTCH %in% batch_rep_ts34) %>%
  filter(SMTS == "Kidney") # 2:BP-67833, 2:BP-75925, 2:BP-49102, 1:BP-47696, 1:BP-44261 = 8
# Bowel
batch_rep_perts %>%
  filter(SMNABTCH %in% batch_rep_ts34) %>%
  filter(SMTS == "Small Intestine") 
                        # 2:BP-44261, 2:BP-48576, 2:BP-56446, 1:BP-47696, 1:BP-75925 = 8
# Liver
batch_rep_perts %>%
  filter(SMNABTCH %in% batch_rep_ts34) %>%
  filter(SMTS == "Liver") # 3:BP-75925, 2:BP-48576, 2:BP-56446, 3: BP-47696 = 10
                          # after filtering
                          # 3:BP-75925, 2:BP-48576, 2:BP-56446, 3: BP-47696 = 10
# Heart -atrial,left ventricle
meta %>%
  group_by(SMTS, SMTSD, SMNABTCH, SMGEBTCH) %>%
  filter(SMTS == "Heart") %>%
  summarize(Count=n()) %>%
  arrange(desc(Count)) %>%
  filter(Count > 1) %>%
  filter(SMNABTCH %in% batch_rep_ts34) 
  # Atrial = 2: BP-67833, 4: BP-44261 = 6 | 3: BP-49102, 3: BP-47696 = 12 
  # Left = 2: BP-75925, 2: BP-48576 = 4 | 2: BP-49102, 3: BP-47696 = 9

# subset for heart = 12  as follows:
  # Atrial = 1: BP-67833, 1: BP-44261 = 2 | 2: BP-49102, 2: BP-47696 = 6
  # Left = 1: BP-75925, 1: BP-48576 = 2   | 2: BP-49102, 2: BP-47696 = 6

# after filtering
# Atrial = 1: BP-67833, 4: BP-44261, 1: BP-48576 = 6 | 3: BP-49102, 3: BP-47696 = 12
  # Left = 2: BP-75925, 2: BP-48576, 1: BP-44261, 1: BP-67833 =6 | 3: BP-49102, 3: BP-47696 = 12

# filter for all  BP-48576(atrial), BP-44261(left), BP-67833(atrial), 
#   3->BP-44261(atrial), 1-> BP-75925(left), 1-> BP-48576(left)
#   2-> BP-49102 (1left, 1atrial), 2-> BP-47696 (1left, 1atrial)

# Lung
batch_rep_perts %>%
  filter(SMNABTCH %in% batch_rep_ts34) %>%
  filter(SMTS == "Lung") 
# batches counts:
# 3: BP-44261, 3:BP-75925, 2: BP-48576, 2:BP-67833, 2: BP-49102, 6:BP-47696 = 18

# subset for lung as follows:
# 2: BP-44261, 2:BP-75925, 2: BP-48576, 2:BP-67833, 2: BP-49102, 2:BP-47696 = 12

# after filtering
# 3: BP-44261, 3:BP-75925, 2: BP-48576, 2:BP-67833, 3: BP-49102, 6:BP-47696 = 16

# filter for all BP-56446 
#     1-> BP-49102-->( LCSET-5305)
#     1-> BP-44261, 1-> BP-75925, 4->BP-47696


# to complete required sample number, take 2 from kidney and intestine
# from batches that are repeated 2 times across organs
batch_rep_perts %>%
  filter(SMNABTCH %in% batch_rep_ts2) %>%
  filter(SMTS == "Kidney") # 2:BP-49102 = 6
```

7 SMNA Batches selected =
`BP-67833, BP-75925, BP-49102, BP-44261, BP-48576, BP-BP-56446, BP-47696`

8 batches in combination with SMGE batch:

    SMNABTCH SMGEBTCH #diff.tissues
    1 BP-48576 LCSET-4955 4
    2 BP-75925 LCSET-9767 5
    3 BP-44261 LCSET-4796 4
    4 BP-47696 LCSET-4907 5
    5 BP-67833 LCSET-8132 3
    6 BP-56446 LCSET-6445 3
    7 BP-49102 LCSET-5305 2
    8 BP-49102 LCSET-5304 2

  

##### **Apply batches filters**

``` r
##########
## Apply filters for all samples over a well-checked subset of batches

# SMNABTCHes chosen, repeated within and across organs
smnabatches_selected <- c("BP-67833","BP-75925","BP-49102","BP-44261",
                      "BP-48576","BP-56446","BP-47696")
# paired SMGEBTCHes
gebatches_selected <- c("LCSET-4955","LCSET-9767","LCSET-4796","LCSET-4907",
                        "LCSET-8132","LCSET-5305","LCSET-5304","LCSET-6445")

# subset for selected batches
meta <- meta %>%
  filter(SMNABTCH %in% smnabatches_selected, SMGEBTCH %in% gebatches_selected)

##########
## Apply further filter for Heart and Lung

# Filter for specific batches on certain organs
#       (the ones that were magically added after filtering for batches)
#   "BP-48576"-atrial
#   "BP-44261"-left
#   "BP-67833"-atrial
#   "BP-56446"-lung
meta <- meta %>%
  filter(
    !(
      (SMTSD == "Heart - Atrial Appendage" & SMNABTCH %in% c("BP-48576","BP-67833")) |
      (SMTSD == "Heart - Left Ventricle" & SMNABTCH == "BP-44261") |
      (SMTS == "Lung" & SMNABTCH == "BP-56446")
      )
    )
# Filter for a subset of samples from a certain batch:
# "BP-44261" 3-atrial,        1-lung +
# "BP-47696" 1-left,1-atrial  4-lung +
# "BP-48576" 1-left                  +
# "BP-49102" 1-left,1-atrial  1-lung-LCSET-5305 +
# "BP-75925" 1-left           1-lung +
exclude.samples <- meta %>%
  filter( (SMTS == "Heart" & 
            SMNABTCH %in% c("BP-44261","BP-75925","BP-48576","BP-49102","BP-47696")) |
          (SMTS == "Lung" & 
             SMNABTCH %in% c("BP-44261","BP-75925","BP-49102","BP-47696"))
          ) %>%
  arrange(SMNABTCH) %>%
  # took the first that met criteria for subsetting
  slice(1:3,6,8:11,13:14,20,22:23,25,31:32) %>%
  select(SAMPID)
exclude.samples <- as.array(exclude.samples$SAMPID)

meta <- filter(meta,!(SAMPID %in% exclude.samples))
##########

# Table with samples to use, with annotations of tissue and batches
write.table(x = meta, 
  file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8_SampleAttributesDS-fltr.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = F, col.names = T, quote = F)
```

    SAMPID  SMTS    SMTSD   SMNABTCH    SMGEBTCH
    GTEX-145MO-2326-SM-5NQ9K    Liver   Liver   BP-47696    LCSET-4907
    GTEX-146FH-1126-SM-5NQAT    Heart   Heart - Atrial Appendage    BP-47696    LCSET-4907
    GTEX-146FH-1526-SM-5NQBU    Liver   Liver   BP-47696    LCSET-4907
    GTEX-14753-1626-SM-5NQ9L    Liver   Liver   BP-47696    LCSET-4907
    GTEX-147F3-0726-SM-5NQ9U    Lung    Lung    BP-47696    LCSET-4907
    GTEX-147F4-1226-SM-5NQAY    Heart   Heart - Left Ventricle  BP-47696    LCSET-4907
    GTEX-147JS-1126-SM-5RQIW    Liver   Liver   BP-48576    LCSET-4955
    ...

------------------------------------------------------------------------

  

##### **Subset samples GTEx v8 gene count matrix**

  

Subset for samples of interest:

``` r
librar(tidyverse)

# gtex data
gtex <- as.data.frame(read.csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads2.gct", sep = "\t"), header = F)
dim(gtex)
head(gtex)[,1:5]
# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8_SampleAttributesDS-fltr.csv", sep = ","))
dim(meta)

samples <- as.array(meta$SAMPID)
samples <- str_replace_all(samples, "-",".")
smpl.fltr <- colnames(gtex) %in% samples
# first two cols are ensg_id and gene_name, I am working with the latter
smpl.fltr[1] = TRUE
smpl.fltr[2] = TRUE
gtex <- gtex[,smpl.fltr]

#write that table
write.table(x = meta, 
  file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-counts.csv", sep = ",", eol = "\n", na = "NA", row.names = F, col.names = T, quote = F)

# Filter for unexpressed genes across tissues
# gtex_fltr <- gtex[which(apply(gtex,1,var) != 0),]
```

------------------------------------------------------------------------

  

##### **Normalization of gene counts**

  

**Count per Million CPM**

For each sample:

1.  Sum counts for all genes
2.  Divide genes\_counts by sum
3.  Multiply genes by 10^6

``` r
# gtex samples of interest with gene ids rownames
gtex.ensg <- column_to_rownames(gtex[,-2], var = "Name")
gtex.ensg <- as.matrix(gtex.ensg)

# 1
sum.counts.samp <- apply(gtex.ensg,2,sum)
gtex.ensg.cpm <- t((t(gtex.ensg)/sum.counts.samp) * 10^6)
```

  

**Transcripts per Million TPM**

1.  Get A. For each gene:

-   A = (sum gene\_counts \* 10^3)/ gene\_length

2.  A/sum(A) \* 10^6

------------------------------------------------------------------------

  

##### **Batch correction**

``` r
# gtex samples of interest with gene ids rownames
#gtex.ensg <- column_to_rownames(gtex[,-2], var = "Name")
# gtex.ensg <- gtex.ensg.cpm

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8_SampleAttributesDS-fltr.csv", sep = ","))
dim(meta)

# order meta as in gtex.ensg
meta$SAMPID <- str_replace_all(meta$SAMPID, "-",".") # slight diff sample notation correction
order.meta <- match(colnames(gtex.ensg), meta$SAMPID)
meta <- meta[order.meta,]

# Define batches (SMNA)
SMNA <- unique(meta$SMNABTCH)
batch <- meta$SMNABTCH
for (i in 1:length(SMNA)) {
  batch <- str_replace(batch, SMNA[i], paste0(i))
  i = i + 1
}
batch <- as.numeric(batch)
# Add GE batches combination (only two diff GE batches in BP-49102)
gebtch.49102 <- meta$SMNABTCH == "BP-49102" & meta$SMGEBTCH == "LCSET-5305"
batch[gebtch.49102] <- 8 

# Define batches according to biological condition
tissue.esp <- rev(unique(meta$SMTSD))
bio.cond <- rev(meta$SMTSD)
for (i in 1:length(tissue.esp)) {
  bio.cond <- str_replace(bio.cond, tissue.esp[i], paste0(i))
  i = i + 1
}
bio.cond <- as.numeric(bio.cond)

# Correction of batch effects with ComBat-seq
gtex.ensg <- as.matrix(gtex.ensg)
library(sva)
adj.gtex <- ComBat_seq(gtex.ensg, batch = batch, group = bio.cond)
```

    Found 8 batches
    Using full model in ComBat-seq.
    Adjusting for 5 covariate(s) or covariate level(s)
    Estimating dispersions
    Fitting the GLM model
    Shrinkage off - using GLM estimates for parameters
    Adjusting the data

Check created object

``` r
class(adj.gtex) # "matrix" "array"
dim(adj.gtex) # 56200    50
object.size(adj.gtex)/1e6 # 27.4 MB
```

Not filtering for genes with zero counts because it follows a joining
with desai’s genes.

*Note:* Object saved as and can be loaded as following.

``` r
# write.csv(adj.gtex, file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts.csv", quote = F)
adj.gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts.csv", sep = ","))
```

------------------------------------------------------------------------

  

##### **Get “best” ENST for each ENSG**

Necessary for batch correction between datasets gtex and desai.
(transcripts in desai have ENST, this matrix cannot be shriked as
gtex’s, but we can convert the names from gtex)

  

``` r
ensg.gtex <- rownames(adj.gtex)

# load script with matrix summarizing functions
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")

## load table of conversion of all posible transcripts for all ensgs and gene names of gencode v26
gencodev26.b <- as.data.frame(read.csv(
          "/mnt/Citosina/amedina/mpadilla/resources/gencode/v26/fltr_transcripts/gencode.v26.basic.annotation.gene_id-trx_id-genename.tsv",
          header = T, sep = "\t"))
common.ensg <- which(gencodev26.b$ensembl_gene_id_version %in% ensg.gtex)
ensembl_ids <- gencodev26.b[common.ensg,"ensembl_transcript_id_version"]
all.ensg.gtex.enst <- gencodev26.b[common.ensg,] ## all posible transcripts in gtex engs ids

length(ensg.gtex %in% gencodev26.b$ensembl_gene_id_version)

## 2. Search best transcipt with Id2GeneName changed so that it searches 1 trx per ensg instead of 1trx per gene
res.gen26 <- mulENST2oneperENSG(ensembl_ids,filter_transcripts = T, 
                                gencode_filter = T, get_gene_ids = T)

## Take enst from first match of the ensg id
  # format
  all.ensg.gtex.enst <- relocate(all.ensg.gtex.enst, ensembl_gene_id_version, 
                                 .after = ensembl_transcript_id_version)
  colnames(all.ensg.gtex.enst) <- colnames(res.gen26)
dif <- setdiff(all.ensg.gtex.enst, res.gen26)
lacking.ensgs <- match(dif$ensembl_gene_id_version, all.ensg.gtex.enst$ensembl_gene_id_version)
rest.ensg.gtex.enst <- distinct(all.ensg.gtex.enst[lacking.ensgs,])


## 3. Join tables of conversion between ids and gene_name
ensg.gtex.enst <- union(res.gen26,rest.ensg.gtex.enst)
```

  

Change ENSGids for ENSTids to enable comparison between datasets:

``` r
## check which ensgs to retrieve and match them
index <- match(rownames(adj.gtex), ensg.gtex.enst$ensembl_gene_id_version)
ensg.gtex.enst2 <- ensg.gtex.enst[index,]

## save the conversion matrix just in case
write.csv(ensg.gtex.enst2, file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/conversion_enst_ensg_genename.csv", quote = F, row.names = F, col.names = T)

# Rename genes as their "best" enst for further comparison
rownames(adj.gtex) <- ensg.gtex.enst2$target_id
```

Write file out

``` r
write.csv(adj.gtex, file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-enst.csv", quote = F, row.names = T, col.names = T)

# adj.gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-enst.csv", sep = ","))
```

##### **Filter genes**

  

pyscenic takes gene names as input, not ids. We need gene names but only
one per gene.

  

``` r
# Filter by ENSTs found before with function that searches best transrcipt fpr a gene name
ensembl_ids <- rownames(adj.gtex)
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T, get_gene_ids = T)
dim(res)

# Format count matrix to rows gene names and cols samples
gtex <- FormatExpCountsMatrix(adj.gtex, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-fltrensg.csv")
# gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-fltrensg.csv", sep = ","))
```

  

  

### **Matrixes per experiment (3)**

  

We’ve decided to split the analysis in three, as shown below.

Runs that include batch correction: `multiruns-bc-sep`

  

#### **Cell lines**

  

##### **Subset for Blanco’s all Cell lines experiments**

``` r
# blanco's data with batch correction
blanco <- as.data.frame(read.csv("blanco2020/counts3/blanco2020-counts_gencv38-alltrxs_trxs-adj.csv", sep = ","))
dim(blanco) # 106143     79

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.mock <- str_detect(meta$source, "Mock")
meta.inf <- str_detect(meta$source, "infected")
meta.ifnb <- str_detect(meta$source, "IFNB_treated_NHBE")
# subset for celllines samples
meta <- meta[meta.mock | meta.inf | meta.ifnb,]

celllines <- blanco[,which(colnames(blanco) %in% meta$sample)]
dim(celllines)# 106143     74

# Write out filtered matrix
write.table(x = celllines, 
            file = "blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-alltrx.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

**Filter trasncripts and get gene names (one per trx)**

``` r
celllines <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-alltrx.csv", sep = ","))
## celllines without infb nhbe
ensembl_ids <- rownames(celllines)

source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")

# get gene names and filteres trxs
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
celllines <- FormatExpCountsMatrix(celllines, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.csv`

Metadata =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease_subCelllines.csv`

-   tpm counts from kallisto
-   Batch corrected with all blanco’s data
-   all celllines experiments
-   (20,299) gene names instead of ensembl ids (one per transcript)

------------------------------------------------------------------------

  

##### **Subset for Blanco’s Cell lines without SARS-CoV-2\_inf\_A549 + ACE2 + pretreatment**

  

``` r
celllines <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE-counts-fltrtrxs.csv", sep = ","))

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease_subCelllines-nopt.csv", sep =",", header = T))

## filter samples
celllines <- celllines[which(rownames(celllines) %in% meta$ID),]
dim(celllines)#   71 20299

## write to file
write.table(x = celllines, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE_nopt-counts-fltrtrxs.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_wINFBNHBE_nopt-counts-fltrtrxs.csv`

Metadata =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source-disease_subCelllines-nopt.csv`

-   tpm counts from kallisto
-   Batch corrected with all blanco’s data
-   all celllines experiments
-   (20,299) gene names instead of ensembl ids (one per transcript)
-   SARS-CoV-2\_inf\_A549 + ACE2 + pretreatment excluded because they
    had no mock control with pretreatment

------------------------------------------------------------------------

  

##### **Subset for Blanco’s Cell lines without hIFNB treated NHBE**

``` r
# blanco's data with batch correction
blanco <- as.data.frame(read.csv("blanco2020/counts3/blanco2020-counts_gencv38-alltrxs_trxs-adj.csv", sep = ","))
dim(blanco) # 106143     79

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.mock <- str_detect(meta$source, "Mock")
meta.inf <- str_detect(meta$source, "infected")
# subset for celllines samples
meta <- meta[meta.mock | meta.inf,]

celllines <- blanco[,which(colnames(blanco) %in% meta$sample)]
dim(celllines) # 106143     68

# Write out filtered matrix
write.table(x = celllines, 
            file = "blanco2020/counts3/blanco2020_adj-celllines_noIFNBHNBE-counts-alltrxs.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

**Filter trasncripts and get gene names (one per trx)**

``` r
celllines <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_noIFNBHNBE-counts-alltrxs.csv", sep = ","))
# res is the same from before

# format matrix with gene names found
celllines <- FormatExpCountsMatrix(celllines, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_noINFBNHBE-counts-fltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/blanco2020/counts3/blanco2020_adj-celllines_noINFBNHBE-counts-fltrtrxs.csv`

-   tpm counts from kallisto
-   Batch corrected with all blanco’s data
-   all celllines experiments, minus hINFB treated NHBE
-   (20,299) gene names instead of ensembl ids (one per transcript)

------------------------------------------------------------------------

#### **COVID organs**

  

##### **Subset for Desai’s organs (liver, bowel, kidney, lung, heart)**

``` r
# desai's data with batch correction
desai <- as.data.frame(read.csv("desai2020/kallisto_gen/desai2020-counts-alltrxs-adj.csv", sep = ","))
dim(desai) # 106143     88

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.lung <- str_detect(meta$source, "lung")
meta.lungh <- str_detect(meta$source, "lung_healthy")
meta.liver <- str_detect(meta$source, "liver")
meta.bowel <- str_detect(meta$source, "bowel")
meta.heart <- str_detect(meta$source, "heart")
meta.kidney <- str_detect(meta$source, "kidney")
# subset for lung samples from desai and blanco only
meta <- meta[(meta.lung & !meta.lungh) | meta.liver | meta.bowel | meta.heart | meta.kidney,]

des.organs <- desai[,which(colnames(desai) %in% meta$sample)]
dim(des.organs) # 106143     77

# Write out filtered matrix
write.table(x = des.organs, 
            file = "desai2020/kallisto_gen/desai2020-counts-alltrxs-adj-selCovidOrgans.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

**Filter transcripts and get gene names (one per trx)**

``` r
# desai's data with batch correction
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")
desai <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020/kallisto_gen/desai2020-counts-alltrxs-adj-selCovidOrgans.csv", sep = ","))
dim(desai) # 106143     72

## 
ensembl_ids <- rownames(desai)
# get gene names and filteres trxs
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
data <- FormatExpCountsMatrix(desai, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020/kallisto_gen/desai2020-counts-adj-selCovidOrgans_fltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/desai2020/kallisto_gen/desai2020-counts-adj-selCovidOrgans_fltrtrxs.csv`

-   tpm counts from kallisto
-   data was corrected by batch effect (date)
-   desai’s covid organs: bowel, lung, liver, heart, kidney
-   (20,299) gene names instead of ensembl ids (one per transcript)

------------------------------------------------------------------------

  

##### **GTEx + Desai’s organs (liver, bowel, kidney, lung, heart)**

  

###### **Load and Join matrixes and metadata**

  

**Count matrixes**

``` r
# desai's data with batch correction and subset for organs of interest
desai <- as.data.frame(read.csv("../desai2020/kallisto_gen/desai2020-counts-alltrxs-adj-selCovidOrgans.csv", sep = ","))
dim(desai) # 106143     72

# gtex with subset of samples, batch correction and enst ids
gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-enst.csv", sep = ","))
dim(gtex) #  56200    50

# join data
desai.g <- rownames_to_column(as.data.frame(desai), var = "enstid")
gtex.g <- rownames_to_column(as.data.frame(gtex), var = "enstid")
multiorgans <- full_join(desai.g, gtex.g, by="enstid")
multiorgans[is.na(multiorgans)] <- 0 # genes that were not present in a given datset are marked as NA
multiorgans <- column_to_rownames(multiorgans,var = "enstid")
multiorgans <- as.matrix(multiorgans)
```

**Metadata**

``` r
# meta for desai done in section COVID organs > Desai subset
meta.gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8_SampleAttributesDS-fltr.csv", sep = ","))
meta.gtex <- select(meta.gtex,SAMPID,SMTS)
colnames(meta.gtex) <- colnames(meta)
dim(meta.gtex)
# join metadata
meta <- union(meta,meta.gtex)

# Change source(tissue) tags to distinguish between healthy and covid
# Define batches according to biological condition
tissue <- rev(unique(meta$source))
bio.cond <- meta$source
tissue.covid <- c("Kidney_healthy","Lung_healthy","Heart_healthy","Liver_healthy","Bowel_healthy",
                  "Heart_COVID","Kidney_COVID","Bowel_COVID","Lung_COVID","Liver_COVID")
for (i in 1:length(tissue)) {
  bio.cond <- str_replace(bio.cond, tissue[i], tissue.covid[i])
  i = i + 1
}
meta$source <- bio.cond

# dataset
dataset <- rep(c(1,2),c(72,50))
meta.sub[,"dataset"] <- dataset

# format
meta$sample <- str_replace_all(meta$sample,"-",".")

# write joined meta
write.table(x = meta, 
            file = "../Data/multi-datasets/metadata/desai2020_gtexv8-subTissues-sample_source.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = F, col.names = T, quote = F)
```

------------------------------------------------------------------------

  

###### **Batch correction by dataset**

In this case, covariates are not used since they coincide with the
batch, this step potentially masks biological signal

``` r
# Define batches (SMNA)
batch <- rep(c(1,2),c(72,50))

library(sva)
adj.multiorg <- ComBat_seq(multiorgans, batch = batch)
```

Check created object

``` r
class(adj.multiorg) # "matrix" "array"
dim(adj.multiorg) # 154598    122
object.size(adj.multiorg)/1e6 # 164.5 MB
```

Not filtering for genes with zero counts because it follows a joining
with desai’s genes.

*Note:* Object saved as and can be loaded as following.

``` r
# write.csv(adj.multiorg, file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_counts-alltrxs.csv", quote = F, row.names = T, col.names = T)
adj.multiorg <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_counts-alltrxs.csv", sep = ","))
```

------------------------------------------------------------------------

  

###### **Filter transcripts and get gene names (one per trx)**

``` r
source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")
## 
ensembl_ids <- rownames(adj.multiorg)

# get gene names and filteres trxs
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
data <- FormatExpCountsMatrix(as.data.frame(adj.multiorg), res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_countsfltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/desai2020_gtexv8-subTissues-adj_countsfltrtrxs.csv`

-   desai : tpm counts from kallisto, gtex: gene counts from STAR

-   data was corrected by batch effect (date) on each dataset

-       then by dataset at merge

-   desai’s covid organs: bowel, lung, liver, heart, kidney

-     + gtex healthy organs: bowel (small intestine), lung, liver, heart (atrial,ventricle), kidney (cortex)

-   (20,299) gene names instead of ensembl ids (one per transcript)

##### **GTEx organs (liver, bowel, kidney, lung, heart)**

  

``` r
gtex <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-fltrensg.csv", sep = ","))
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/gtex_v8-fltr/GTEx_v8-RNASEQ-subTissues-adj-counts-fltrensg.csv`

-   gtex: gene counts from STAR normalized to CPM
-   data was corrected by batch effect SMNA and GE
-   gtex healthy organs: bowel (small intestine), lung, liver, heart
    (atrial,ventricle), kidney (cortex)
-   (13,697) gene names instead of ensembl ids (one per transcript)

#### **Lung COVID**

  

##### **Lung COVID and Healthy (Blanco + Desai, one batch correction)**

``` r
# blanco's data with batch correction
blanco <- as.data.frame(read.csv("blanco2020/counts3/blanco2020-counts_gencv38-alltrxs_trxs-adj.csv", sep = ","))
dim(blanco) # 106143     78

# desai's data with batch correction
desai <- as.data.frame(read.csv("desai2020/kallisto_gen/desai2020-counts-alltrxs-adj.csv", sep = ","))
dim(desai) # 106143     88

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.lung <- str_detect(meta$source, "[Ll]ung")
meta.lungdel <- str_detect(meta$source, "Lung_d1")
# subset for lung samples from desai and blanco only
meta <- meta[meta.lung & !meta.lungdel,]

lung.bla <- blanco[,which(colnames(blanco) %in% meta$sample)]
lung.des <- desai[,which(colnames(desai) %in% meta$sample)]
dim(lung)# 61 19629

## Join samples from blanco and desai
# since both datasets have been through the exact same preprocessing, their rows are equal
# length(which(rownames(lung.bla) == rownames(lung.des))) # so we can do this:
lung <- bind_cols(lung.bla,lung.des)

# Write out filtered matrix
write.table(x = lung, 
            file = "multi-datasets/counts/blanco2020_adj-desai2020_adj-Lung-counts-alltrxs.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

**Filter trasncripts and get gene names (one per trx)**

``` r
lung <- as.data.frame(read.csv(file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020_adj-desai2020_adj-Lung-counts-alltrxs.csv", sep = ","))

## 
ensembl_ids <- rownames(lung)
# get gene names and filteres trxs
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
lung <- FormatExpCountsMatrix(lung, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020_adj-desai2020_adj-Lung-counts-fltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020_adj-desai2020_adj-Lung-counts-fltrtrxs.csv`

-   tpm counts from kallisto
-   data from diff datasets was corrected separatedly
-   all lung covid and healthy biopsies from autopsy
-   (20,299) gene names instead of ensembl ids (one per transcript)

------------------------------------------------------------------------

  

##### **Lung COVID and Healthy (Blanco + Desai, two batch corrections)**

``` r
# Blades data with correction by dataset
blades <- as.data.frame(read.csv("multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-alltrxs-adj.csv", sep = ","))
dim(blades) # 106143    166

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.lung <- str_detect(meta$source, "[Ll]ung")
meta.lungdel <- str_detect(meta$source, "Lung_d1")
# subset for lung samples from desai and blanco only
meta <- meta[meta.lung & !meta.lungdel,]

lung <- blades[,which(colnames(blades) %in% meta$sample)]

# Write out filtered matrix
write.table(x = lung, 
            file = "multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-Lung-alltrxs.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

**Filter trasncripts and get gene names (one per trx)**

``` r
lung <- as.data.frame(read.csv(file="/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-Lung-alltrxs.csv", sep = ","))
## 
# ensembl_ids <- rownames(lung)
# 
# # get gene names and filteres trxs
# res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
lung <- FormatExpCountsMatrix(lung, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-Lung-fltrtrxs.csv")
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-Lung-fltrtrxs.csv`

-   tpm counts from kallisto
-   data from diff datasets was corrected separatedly (by dates) and
    then joined (by dataset)
-   all lung covid and healthy biopsies from autopsy
-   (20,299) gene names instead of ensembl ids (one per transcript)

------------------------------------------------------------------------

  

##### **Lung COVID and Healthy (Desai + Blanco + Delorey)**

  

-   Join matrixes blades with dates and dataset correction and delorey
    cpm

Lung COVID and Healthy (Desai + Blanco + Delorey)

``` r
# unfiltered dataset-batch-corrected Blanco+Desai matrix with gencode
# Blades data with correction by dataset
blades.un <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-mer-adj-counts-alltrxs-adj.csv", sep = ","))
dim(blades.un) # 106143    166

# Delorey lung samples cpm counts
del <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/delorey2021/counts/GSE171668-Lung_rsem.enst-genes.cpm_counts.matrix.csv", sep = ","))
dim(del) # 19840    17

# metadata
meta <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv", sep ="\t", header = F))
colnames(meta) <- c("sample", "source")

meta.lung <- str_detect(meta$source, "[Ll]ung")
# subset for lung samples
meta <- meta[meta.lung,]

lung.blades <- blades.un[,which(colnames(blades.un) %in% meta$sample)]
dim(lung.blades) #  106143     61

## Join blades lung and delorey lung data
lung.blades.g <- rownames_to_column(as.data.frame(lung.blades), var = "enstid")
del.g <- rownames_to_column(del, var = "enstid")
bladesdel.lung <- full_join(lung.blades.g, del.g, by = "enstid")
bladesdel.lung <- column_to_rownames(bladesdel.lung, var = "enstid")
dim(bladesdel.lung) # 109700     78
```

-   Batch correction between datasests ( desai-blanco and delorey)

``` r
# match meta and counts matrix
colnames(bladesdel.lung) <- str_replace(colnames(bladesdel.lung),"X","")
colnames(bladesdel.lung) <- str_replace_all(colnames(bladesdel.lung),"_","-")
order <- match(colnames(bladesdel.lung), meta$sample)
meta <- meta[order,]

# Define batches according to dataset
blades.samp <- str_detect(meta$sample, "^[GSMAN]")
batch <- rep(1,length(meta$sample))
batch[blades.samp] <- 2

# change tags across dataset to match biological condition
tissue <- rev(unique(meta$source))
bio.cond <- meta$source
tissue.covid <- c("Lung_COVID","Lung_healthy","Lung_COVID","Lung_COVID","Lung_healthy")
for (i in 1:length(tissue)) {
  bio.cond <- str_replace(bio.cond, tissue[i], tissue.covid[i])
  i = i + 1
}
meta$source <- bio.cond

# Define covariates for biological condition
tissue.esp <- rev(unique(meta$source))
bio.cond <- meta$source
for (i in 1:length(tissue.esp)) {
  bio.cond <- str_replace(bio.cond, tissue.esp[i], paste0(i))
  i = i + 1
}
bio.cond <- as.numeric(bio.cond)

# Correction of batch effects with ComBat-seq
bladesdel.lung <- as.matrix(bladesdel.lung)

library(sva)
adj.bladesdel.lung <- ComBat_seq(bladesdel.lung, batch = batch, group = bio.cond)
```

    Found 2 batches
    Using full model in ComBat-seq.
    Adjusting for 1 covariate(s) or covariate level(s)
    Estimating dispersions
    Fitting the GLM model
    Shrinkage off - using GLM estimates for parameters
    Adjusting the data

Check created object

``` r
class(adj.bladesdel.lung) # "matrix" "array"
dim(adj.bladesdel.lung) # 109700     78
object.size(adj.bladesdel.lung)/1e6 # 78 MB
```

Save to file:

``` r
write.table(x = adj.bladesdel.lung, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-alltrxs.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

  

Change ensts for gene names, perform transcript filtering.

``` r
adj.bladesdel.lung <- as.data.frame(read.csv("/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-alltrxs.csv", header = T, sep = ",", quote = ""))

source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")

# get gene names and filteres trxs
ensembl_ids <- rownames(adj.bladesdel.lung)
res <- Id2GeneName(ensembl_ids, filter_transcripts = T, gencode_filter = T)

# format matix with gene names found
adj.bladesdel.lung <- as.data.frame(adj.bladesdel.lung)
adj.bladesdel.lung2 <- FormatExpCountsMatrix(adj.bladesdel.lung, res, filter_transcripts = T, outfile_name = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs.csv")
```

  

**Filter genes with 0 counts and NAs**

  

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
## save gene-filtered matrix
write.table(x = lung, 
            file = "/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs-gfltr.csv", 
            sep = ",", eol = "\n", na = "NA", row.names = T, col.names = T, quote = F)
```

Matrix =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/counts/blanco2020-desai2020-delorey2021-mer-adj-counts-Lung-fltrtrxs-gfltr.csv`

-   tpm counts from kallisto
-   data from diff datasets was corrected separatedly (by dates) and
    then joined (by dataset)
-   all lung covid and healthy biopsies from autopsy
-   genes with zero counts filtered
-   (19,919) gene names instead of ensembl ids (one per transcript)

metadata =
`/mnt/Citosina/amedina/mpadilla/COVID19/Data/multi-datasets/metadata/meta_sample-source.tsv`

  

------------------------------------------------------------------------

### **Refs**

  

In collaboration with Leonardo Arteaga, Ana Betty Villaseñor Altamirano,
Karen Nuñez Reza and Alejandra Medina Rivera.

  

**References:**

\[1\] Delorey, T. M., Ziegler, C. G., Heimberg, G., Normand, R., Yang,
Y., Segerstolpe, Å., … & Regev, A. (2021). COVID-19 tissue atlases
reveal SARS-CoV-2 pathology and cellular targets. Nature, 1-8.

\[2\] Desai, N., Neyaz, A., Szabolcs, A., Shih, A. R., Chen, J. H.,
Thapar, V., … & Deshpande, V. (2020). Temporal and spatial heterogeneity
of host response to SARS-CoV-2 pulmonary infection. Nature
communications, 11(1), 1-15.

\[3\] Li, B., & Dewey, C. N. (2011). RSEM: accurate transcript
quantification from RNA-Seq data with or without a reference genome. BMC
bioinformatics, 12(1), 1-16.

  

**Programs and versions:**

-   R 4.0.2

    packages:

    -   tidyverse
        -   dplyr
        -   ggplot2
    -     

*Last update: March 23rd,2022*