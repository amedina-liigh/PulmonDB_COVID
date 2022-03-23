## MakeExperimentCountMatrixes.R
## mpadilla@lcgej.unam.mx
## Last update = Nov 26, 2021
# source("/mnt/Citosina/amedina/mpadilla/COVID19/scripts/MakeExperimentCountMatrixes.R")

#################################################################################################
## Function: Join counts matrixes from several samples (rna-seq kallisto quant runs)
##           Output: Data frame. Experiments count matrix
##                  Format: Cols : target_id, sample1_counts, ..., sampleN_counts
##                  count type: est_counts 4 and tpm 5

JoinSamplesCountsMatrixes <- function(files, counts_col = 5, sample_regex){

  print("In function: JoinSamplesCountsMatrixes")

  library(stringr)
  i = 1

  for (file in files) {
    ## Get id and tpm columns from counts file
    counts_file = as.data.frame(read.csv(file = file, header = TRUE, sep = "\t"))
    sample_name = str_extract(file, sample_regex)
    colnames(counts_file)[counts_col] = sample_name # replace counts_col (tpm counts col) for name of sample

    ## Sort table by transcript id
    sort_index = match(sort(counts_file[,1]), counts_file[,1])
    counts_file = counts_file[sort_index,]
    # Note: ids are not repeated

    ## Add new counts to experiment counts table
    if (i == 1){
      # first counts table
      exp_counts = counts_file[,c("target_id", sample_name)]
      i = i+1
    } else {
      # # 2-n counts table

      # merge matrixes
      exp_counts = full_join(x = exp_counts, y = counts_file[,c("target_id", sample_name)], by = "target_id")

      # sort exp_counts by target_id again
      sort_index = match(sort(exp_counts[,1]), exp_counts[,1])
      exp_counts = exp_counts[sort_index,]
    } # else 2-n table
  } # for file

  print(head(exp_counts))
  print(paste0("Generated matrix is dim ", dim(exp_counts)))

  return(exp_counts)
}

#################################################################################################
## Function: Convert transcript ids (w/ version) to gene_names with ensemblv103 and filter

Id2GeneName <- function(ensembl_ids, dataset = "hsapiens_gene_ensembl",
                        filter_transcripts = F, gencode_filter = F, get_gene_ids = F,
                        host_ensembl = "https://feb2021.archive.ensembl.org"){

  print("In function: Id2GeneName")

  if( dataset == "ferret"){
    dataset = "mpfuro_gene_ensembl"
  }

  # Use biomart package from Bioconductor
  library(biomaRt)
  library(stringr)
  library(dplyr)

  ## Check ensembl version to use
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                 dataset = dataset,
                 host = host_ensembl) # check listEnsemblArchives()

  print(paste("   biomart ready with dataset", dataset))


  if(filter_transcripts){

    if(gencode_filter == T){
      ## gencode transcript ids are bizarre so...
      ensembl_ids <- str_extract(ensembl_ids,"ENST\\d+\\.\\d+")

      # Use this when reference index was gencode v38
      ## Filter 1 according to GENCODE v38 basic protein coding genes annotation (differs from GENCODE basic here)
        gencodev38 <- as.data.frame(read.csv(
          "/mnt/Citosina/amedina/mpadilla/resources/gencode/v38/fltr_transcripts/gencode.v38.basic.pc.transcripts.ensgid_enstid_genename.csv",
          header = T, sep = ","))

      # apply filter for kallisto's matrix transcripts
      gencode.fltr <- which(gencodev38$ensembl_transcript_id_version %in% ensembl_ids)
      gencodev38 <- gencodev38[gencode.fltr,]

      # Submit duplicated transcripts to the next filters
        # First subset ensembl ids and match them with filtered gencode annotation table
        gencode.fltr <- ensembl_ids %in% gencodev38$ensembl_transcript_id_version
        ensembl_ids <- ensembl_ids[gencode.fltr]
        ensembl_ids <- ensembl_ids[match(ensembl_ids, gencodev38$ensembl_transcript_id_version)]
      dup <- duplicated(gencodev38$gene_name) | duplicated(gencodev38$gene_name, fromLast=TRUE)
      ensembl_ids <- ensembl_ids[dup]

      # Filter table for genes with one transcript
      gencodev38 <- gencodev38[!dup,] # 6324 transcripts = genes

      ## Note: this filter was applied first bc when submitting transripts to biomaRT
      ##       some genes where removed (580)
    }

    # Create a query and send it to biomart server, takes 1 min
    res = getBM(attributes = c('ensembl_transcript_id_version',
                               'ensembl_gene_id_version',
                               'external_gene_name',
                               'transcript_tsl',
                               'transcript_appris',
                               'transcript_gencode_basic'),
                filters = 'ensembl_transcript_id_version',
                values = ensembl_ids,
                mart = mart)
      # 50 transcripts lost from previous step

    ## Filter 1 : according to a list of preferable transcript flags (without removing genes)
    fltr <- res$transcript_gencode_basic == "GENCODE basic" |
        (res$transcript_tsl == "tsl1" | res$transcript_tsl == "tsl2" | res$transcript_tsl == "tsl3") |
        (res$transcript_appris == "principal1" | res$transcript_appris == "principal2" |
           res$transcript_appris == "principal3")
    res <- res[fltr,]

    ############################
    ## Filter 2 : search for best transcript for genes with multiple
    dup <- duplicated(res$external_gene_name) | duplicated(res$external_gene_name, fromLast=TRUE)

    # Keep transcripts with best supporting evidence
    genesdup <- length(unique(res[dup,]$external_gene_name)) # duplicated genes
    # filter genes with this priority
    appris <- c("principal1", "principal2", "principal3", "principal4")
    # all categories from this method so no genes are eliminated, but selected in this priority
    tsl_fltr_regex <- c("tsl1","tsl2","tsl3","tsl4", "tsl5", "tslNA", "")
    gencode_basic <- c("GENCODE basic", "")
    # only check for duplicated genes hereafter
    transcripts_checked <- !dup # results in true for checked ones
    i = 1

    # select according to appris and tsl flags
    for(principal in appris){
      for(tsl in tsl_fltr_regex){
        # transcripts to keep, remove those for which genes have already been selected
        fltrdup <- dup & res$transcript_appris == principal & # appris is the principal filter here
          res$transcript_gencode_basic == "GENCODE basic" &
          str_detect(res$transcript_tsl, pattern = tsl ) & !transcripts_checked

        # keep track of transcripts to keep
        if (i == 1){
          fltrdup_final <- fltrdup
          genesrm <- ""
          i = i+1
        } else{
          fltrdup_final <- fltrdup_final | fltrdup
        }

        # ignore those transcripts' genes for next round
        genesfltr <- unique(res[which(fltrdup),]$external_gene_name) # genes from previous filter
        genesrm <- c(genesrm, genesfltr)
        transcripts_checked <- res$external_gene_name %in% genesrm | !dup

        #keep track of genes left
        genesdup <- genesdup - length(genesfltr)
      }
    }

    # In case there are still some genes left
    # select according to appris and gencode basic flags
    if( genesdup > 1 ){
      for(principal in appris){
        for(basic in gencode_basic){
          fltrdup <- dup & res$transcript_appris == principal & # appris is the principal filter here
            res$transcript_gencode_basic == basic & !transcripts_checked

          fltrdup_final <- fltrdup_final | fltrdup

          # ignore those transcripts' genes for next round
          genesfltr <- unique(res[which(fltrdup),]$external_gene_name) # genes from previous filter
          genesrm <- c(genesrm, genesfltr)
          transcripts_checked <- res$external_gene_name %in% genesrm | !dup

          #keep track of genes left
          genesdup <- genesdup - length(genesfltr)
        }
      }
    }

    # In case there are still some genes left
    # select according to tsl and gencode flags
    if( genesdup > 1 ){
      for(tsl in tsl_fltr_regex){
        for(basic in gencode_basic){
          fltrdup <- dup & str_detect(res$transcript_tsl, pattern = tsl ) & # appris is the principal filter here
            res$transcript_gencode_basic == basic & !transcripts_checked

          fltrdup_final <- fltrdup_final | fltrdup

          # ignore those transcripts' genes for next round
          genesfltr <- unique(res[which(fltrdup),]$external_gene_name) # genes from previous filter
          genesrm <- c(genesrm, genesfltr)
          transcripts_checked <- res$external_gene_name %in% genesrm | !dup

          #keep track of genes left
          genesdup <- genesdup - length(genesfltr)
        }
      }

      fltrdup <- dup & res$transcript_gencode_basic == "GENCODE basic" & !transcripts_checked
      fltrdup_final <- fltrdup_final | fltrdup
      genesfltr <- unique(res[which(fltrdup),]$external_gene_name)
      genesdup <- genesdup - length(genesfltr)
    } ## end of Filter 2
    ############################

    ## Keep transcripts that passed filters
    res <- res[!dup | fltrdup_final,]

    ## Some genes are repited because they have the same id and others bc they have the same flags
    res <- distinct(res[,c(1:3)])

    ## Add transcripts table and gene names of gencode filter was used first
    if(gencode_filter == T){
      gencodev38 <- relocate(gencodev38, ensembl_gene_id_version, .after = ensembl_transcript_id_version)
      colnames(gencodev38)[3] <- "external_gene_name"
      gencodev38 <- union(gencodev38,res)

      # remove duplicates from ensembl
      dup <- duplicated(gencodev38$external_gene_name)
      gencodev38 <- gencodev38[!dup,]

      # format as before
      res <- gencodev38[,c(1,2,3)]
    } else{
      res <- res[,c(1,2,3)]
    }

    print(paste0(" Transcripts # after filter = ",dim(res)[1])) # 41791     5
    print(paste0(" Genes # after filter = ", length(unique(res$external_gene_name)) ))

  } else {
    ## No filter to apply. Keep all transcripts
    res = getBM(attributes = c('ensembl_transcript_id_version',
                               'ensembl_gene_id_version',
                               'external_gene_name'),
                filters = 'ensembl_transcript_id_version',
                values = ensembl_ids,
                mart = mart)
  }

  # Sort results (transcript id and gene name df) by transcript_id
  sort_index = match(sort(res$ensembl_transcript_id_version),
                     res$ensembl_transcript_id_version)
  res = res[sort_index,]
  # Note: several transcripts can correspond to the same gene
  colnames(res)[1] = "target_id" # so it matches with counts matrix col

  if(get_gene_ids){
    res <- res[,c("target_id","ensembl_gene_id_version","external_gene_name")]
  } else {
    res <- res[,c("target_id","external_gene_name")]
  }

  return(res)
}


#################################################################################################
## Function: Replace target_ids for gene_names in experiment counts matrix
##          Output: Experiment counts matrix
##                  Format: df: samples (rows) x genes (columns)

FormatExpCountsMatrix <- function(exp_counts, res, transposed = T, filter_transcripts = F,
                                  output_matrix_transcripts_ids = F, gencode_filter = F,
                                  outfile_name = outfile_name){

  print("In function: FormatExpCountsMatrix")

  # In case the input matrix does not come from this pipeline
  if(!"target_id" %in% names(exp_counts)){
    print("Adjusting format\nChanging transcript row to column and calling it target_id")

    exp_counts <- rownames_to_column(exp_counts, var = "target_id")
  }

  ## Index of ids for which an external_gene_name was found
  index_genenames = match(res$target_id, exp_counts$target_id)

  if(filter_transcripts){
    # Filter matrix by transcripts that passed the filter
    exp_counts <- exp_counts[index_genenames,]

    # Add gene names
    exp_counts[, "external_gene_name"] = res$external_gene_name
  } else {
    # Add gene names where it was found
    # exp_counts[index_genenames, "external_gene_name"] = res$external_gene_name[index_genenames]

    # Replace na's in gene_name col with ensembl_ids
    exp_counts$external_gene_name[-index_genenames] = exp_counts$target_id[-index_genenames]
  }


  # Make new experiments counts matrix: samples (rows) x gene names (columns)
  if(transposed){
    only_counts = exp_counts %>% select_if(is.numeric)
    exp_counts_t = as.data.frame(t(as.matrix( only_counts )))
    colnames(exp_counts_t) = exp_counts$external_gene_name

    # output experiment counts matrix
    write.table(x = exp_counts_t, file = outfile_name, sep = ",", eol = "\n",
                na = "NA", row.names = T, col.names = T, quote = F)

    if(output_matrix_transcripts_ids==T){
      exp_counts_t.trx <- exp_counts_t
      colnames(exp_counts_t.trx) = exp_counts$target_id
      # output experiment counts matrix
      write.table(x = exp_counts_t.trx, file = str_replace(outfile_name, ".csv","_trxs.csv"), sep = ",", eol = "\n",
                  na = "NA", row.names = T, col.names = T, quote = F)
    }

  }
  # else {
  #   exp_counts_t = exp_counts %>% relocate(external_gene_name, .before = target_id)
  #   exp_counts_t = exp_counts_t[,-2] # erase transcript_id col
  # }

  return(exp_counts_t)
}

#################################################################################################
## Function: Join the actions of the functions before
##            to make a single matrix joining the gene counts of several samples
##            coming from a single experiment and output that to csv
## sample_regex : regex to identify samples names in file name, like "GSM\\d+" or "SAMN\\d+"

MakeExpCountsMatrix <- function(files, sample_regex = "GSM\\d+", outfile_name = "experiment_counts.csv",
                                transposed = T, counts_col = 5, dataset = "hsapiens_gene_ensembl",
                                host_ensembl = "https://feb2021.archive.ensembl.org",
                                filter_transcripts = F, gencode_filter = F, use_gencode = F,
                                output_matrix_transcripts_ids = F){

  print("In function: MakeExpCountsMatrix")

  library(dplyr)

  # 1. Join counts matrixes from all biological replicates/ samples
  exp_counts = JoinSamplesCountsMatrixes(files, sample_regex = sample_regex)

  # 2. Get conversion of transcript ids to gene_names
  ensembl_ids = as.array(exp_counts[, "target_id"]) # transcript_ids from all samples
  if(gencode_filter){
    use_gencode = T
  }
  if(use_gencode){
    ensembl_ids <- str_extract(ensembl_ids,"ENST\\d+\\.\\d+")
    exp_counts$target_id <- ensembl_ids
  }
  id_gene_df =  Id2GeneName(ensembl_ids, dataset = dataset, filter_transcripts = filter_transcripts,
                            gencode_filter = gencode_filter, host_ensembl = host_ensembl)

  # 3. Adjust format of experiment matrix: samples x genes
  exp_counts_t = FormatExpCountsMatrix(exp_counts = exp_counts, res = id_gene_df, transposed = transposed,
                                       filter_transcripts = filter_transcripts,
                                       output_matrix_transcripts_ids = output_matrix_transcripts_ids,
                                       outfile_name = outfile_name, counts_col = counts_col)

}


#################################################################################################
## Function: Convert transcript ids (w/ version) to gene_names with ensemblv103 and filter
## change: filter 1 trx per 1 eng

mulENST2oneperENSG <- function(ensembl_ids, dataset = "hsapiens_gene_ensembl",
                               filter_transcripts = F, gencode_filter = F, get_gene_ids = F,
                               host_ensembl = "https://feb2021.archive.ensembl.org"){

  print("In function: Id2GeneName")

  # Use biomart package from Bioconductor
  library(biomaRt)
  library(stringr)

  ## Check ensembl version to use
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                 dataset = dataset,
                 host = host_ensembl) # check listEnsemblArchives()

  print(paste("   biomart ready with dataset", dataset))


  if(filter_transcripts){

    if(gencode_filter == T){
      ## gencode transcript ids are bizarre so...
      #ensembl_ids <- str_extract(ensembl_ids,"ENST\\d+\\.\\d+")

      # Use this when reference index was gencode v38
      ## Filter 1 according to GENCODE v38 basic protein coding genes annotation (differs from GENCODE basic here)
      gencodev38 <- as.data.frame(read.csv(
        "/mnt/Citosina/amedina/mpadilla/resources/gencode/v26/fltr_transcripts/gencode.v26.basic.annotation.gene_id-trx_id-genename.tsv",
        header = T, sep = "\t"))

      # apply filter for kallisto's matrix transcripts
      gencode.fltr <- which(gencodev38$ensembl_transcript_id_version %in% ensembl_ids)
      gencodev38 <- gencodev38[gencode.fltr,]

      # Submit duplicated transcripts to the next filters
      # First subset ensembl ids and match them with filtered gencode annotation table
      gencode.fltr <- ensembl_ids %in% gencodev38$ensembl_transcript_id_version
      ensembl_ids <- ensembl_ids[gencode.fltr]
      ensembl_ids <- ensembl_ids[match(ensembl_ids, gencodev38$ensembl_transcript_id_version)]
      dup <- duplicated(gencodev38$ensembl_gene_id_version) |
        duplicated(gencodev38$ensembl_gene_id_version, fromLast=TRUE)
      ensembl_ids <- ensembl_ids[dup]

      # Filter table for genes with one transcript
      gencodev38 <- gencodev38[!dup,] # 6324 transcripts = genes

      ## Note: this filter was applied first bc when submitting transripts to biomaRT
      ##       some genes where removed (580)
    }

    # Create a query and send it to biomart server, takes 1 min
    res = getBM(attributes = c('ensembl_transcript_id_version',
                               'ensembl_gene_id_version',
                               'external_gene_name',
                               'transcript_tsl',
                               'transcript_appris',
                               'transcript_gencode_basic'),
                filters = 'ensembl_transcript_id_version',
                values = ensembl_ids,
                mart = mart)
    # 50 transcripts lost from previous step

    ## Filter 1 : according to a list of preferable transcript flags (without removing genes)
    fltr <- res$transcript_gencode_basic == "GENCODE basic" |
      (res$transcript_tsl == "tsl1" | res$transcript_tsl == "tsl2" | res$transcript_tsl == "tsl3") |
      (res$transcript_appris == "principal1" | res$transcript_appris == "principal2" |
         res$transcript_appris == "principal3")
    res <- res[fltr,]

    ############################
    ## Filter 2 : search for best transcript for genes with multiple
    dup <- duplicated(res$ensembl_gene_id_version) | duplicated(res$ensembl_gene_id_version, fromLast=TRUE)

    # Keep transcripts with best supporting evidence
    genesdup <- length(unique(res[dup,]$ensembl_gene_id_version)) # duplicated genes
    # filter genes with this priority
    appris <- c("principal1", "principal2", "principal3", "principal4")
    # all categories from this method so no genes are eliminated, but selected in this priority
    tsl_fltr_regex <- c("tsl1","tsl2","tsl3","tsl4", "tsl5", "tslNA", "")
    gencode_basic <- c("GENCODE basic", "")
    # only check for duplicated genes hereafter
    transcripts_checked <- !dup # results in true for checked ones
    i = 1

    # select according to appris and tsl flags
    for(principal in appris){
      for(tsl in tsl_fltr_regex){
        # transcripts to keep, remove those for which genes have already been selected
        fltrdup <- dup & res$transcript_appris == principal & # appris is the principal filter here
          res$transcript_gencode_basic == "GENCODE basic" &
          str_detect(res$transcript_tsl, pattern = tsl ) & !transcripts_checked

        # keep track of transcripts to keep
        if (i == 1){
          fltrdup_final <- fltrdup
          genesrm <- ""
          i = i+1
        } else{
          fltrdup_final <- fltrdup_final | fltrdup
        }

        # ignore those transcripts' genes for next round
        genesfltr <- unique(res[which(fltrdup),]$ensembl_gene_id_version) # genes from previous filter
        genesrm <- c(genesrm, genesfltr)
        transcripts_checked <- res$ensembl_gene_id_version %in% genesrm | !dup

        #keep track of genes left
        genesdup <- genesdup - length(genesfltr)
      }
    }

    # In case there are still some genes left
    # select according to appris and gencode basic flags
    if( genesdup > 1 ){
      for(principal in appris){
        for(basic in gencode_basic){
          fltrdup <- dup & res$transcript_appris == principal & # appris is the principal filter here
            res$transcript_gencode_basic == basic & !transcripts_checked

          fltrdup_final <- fltrdup_final | fltrdup

          # ignore those transcripts' genes for next round
          genesfltr <- unique(res[which(fltrdup),]$ensembl_gene_id_version) # genes from previous filter
          genesrm <- c(genesrm, genesfltr)
          transcripts_checked <- res$ensembl_gene_id_version %in% genesrm | !dup

          #keep track of genes left
          genesdup <- genesdup - length(genesfltr)
        }
      }
    }

    # In case there are still some genes left
    # select according to tsl and gencode flags
    if( genesdup > 1 ){
      for(tsl in tsl_fltr_regex){
        for(basic in gencode_basic){
          fltrdup <- dup & str_detect(res$transcript_tsl, pattern = tsl ) & # appris is the principal filter here
            res$transcript_gencode_basic == basic & !transcripts_checked

          fltrdup_final <- fltrdup_final | fltrdup

          # ignore those transcripts' genes for next round
          genesfltr <- unique(res[which(fltrdup),]$ensembl_gene_id_version) # genes from previous filter
          genesrm <- c(genesrm, genesfltr)
          transcripts_checked <- res$ensembl_gene_id_version %in% genesrm | !dup

          #keep track of genes left
          genesdup <- genesdup - length(genesfltr)
        }
      }

      fltrdup <- dup & res$transcript_gencode_basic == "GENCODE basic" & !transcripts_checked
      fltrdup_final <- fltrdup_final | fltrdup
      genesfltr <- unique(res[which(fltrdup),]$ensembl_gene_id_version)
      genesdup <- genesdup - length(genesfltr)
    } ## end of Filter 2
    ############################

    ## Keep transcripts that passed filters
    res <- res[!dup | fltrdup_final,]

    ## Some genes are repited because they have the same id and others bc they have the same flags
    res <- distinct(res[,c(1:3)])

    ## Add transcripts table and gene names of gencode filter was used first
    if(gencode_filter == T){
      gencodev38 <- relocate(gencodev38, ensembl_gene_id_version, .after = ensembl_transcript_id_version)
      colnames(gencodev38)[3] <- "external_gene_name"
      gencodev38 <- union(gencodev38,res)

      # remove duplicates from ensembl
      dup <- duplicated(gencodev38$ensembl_gene_id_version)
      gencodev38 <- gencodev38[!dup,]

      # format as before
      res <- gencodev38[,c(1,2,3)]
    } else{
      res <- res[,c(1,2,3)]
    }

    print(paste0(" Transcripts # after filter = ",dim(res)[1])) # 41791     5
    print(paste0(" Genes # after filter = ", length(unique(res$ensembl_gene_id_version)) ))

  } else {
    ## No filter to apply. Keep all transcripts
    res = getBM(attributes = c('ensembl_transcript_id_version',
                               'ensembl_gene_id_version',
                               'external_gene_name'),
                filters = 'ensembl_transcript_id_version',
                values = ensembl_ids,
                mart = mart)
  }

  # Sort results (transcript id and gene name df) by transcript_id
  sort_index = match(sort(res$ensembl_transcript_id_version),
                     res$ensembl_transcript_id_version)
  res = res[sort_index,]
  # Note: several transcripts can correspond to the same gene
  colnames(res)[1] = "target_id" # so it matches with counts matrix col

  if(get_gene_ids){
    res <- res[,c("target_id","ensembl_gene_id_version","external_gene_name")]
  } else {
    res <- res[,c("target_id","external_gene_name")]
  }

  return(res)
}
