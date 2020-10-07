salmon_load = function (data_table=NULL, salmon_dir=NULL, DESeq_design=NULL, control=NULL, make_tx2gene=F, annotation=NULL ,returnMode = "all",
                        batch_correction=F, batch_group=NULL, k=1)
{
  library("tximport")
  library("GenomicFeatures")
  library("DESeq2")
  library("RUVSeq")
  #library(rtracklayer)
  #library(readr)
  
  ## This function can load salmon quantification files and perform DESeq2 differential expression analysis.
  ## Batch correction is supported. The function can output various corrected and uncorrected matrices, such as raw counts, TPMs, vst
  ## transformed counts and DESeq results.
  ##
  ## Arguments
  ## data_table                                                A delimiter separated table with data of one sample on each line
  ##                                                           sample names as the first line entry (as row names wih no header for the column).
  ##                                                           All other columns must have a header/colname. One of the columns must contain 
  ##                                                           the name of salmon_dir subfolders that hold the salmon quant file for a sample.
  ## salmon_dir                                                Directory that holds subfolders with salmon output.
  ## DESeq_design                                              The variable that DESeq will test the effect of (must be a column of data_table).
  ## control                                                   Control samples from the column that is indicated in DESeq_design.
  ## make_tx2gene                                              Boolean option to create a tx2gene database, will use a loaded object called tx2gene
  ##                                                           in the operation if false. tx2gene is saved as an object after the first run.
  ##                                                           You may only specify this option and annotation to only create the tx2gene object.
  ## annotation                                                Annotation GTF or GFF3 files that are used to convert transcript ids to gene ids.
  ##                                                           This should be a whole genome annotation file from the same supplier as the fasta 
  ##                                                           file that was used to create the salmon index.
  ## returnMode                                                Character vector of datasets to return (all, counts, TPMs, vst, dds). If combined 
  ##                                                           with batch_correction, corrected versions of objects are saved. Note that counts saved
  ##                                                           are raw counts, which are filtered later for use in vst and dds. Use counts(dds) or
  ##                                                           counts(vst) to get the filtered counts
  ## batch_correction
  ## batch_group
  ## k
  ################################################################################################################################################
  ## FUNCTIONS
  RUVs_rep_matrix = function(meta_table=NULL, similarity=NULL){
    # Make a matrix where every row contains indices of replicate samples, based on the values in column "similarity" of "meta_table"
    samples = read.table(file.path(meta_table), header = TRUE)     # Read the table indicated in "meta_table"
    if (isFALSE(similarity %in% colnames(samples))){
      stop(paste0("similarity '", similarity, "' must be a column of '", meta_table,"'."))
    }
    longest = 0                                                    # 
    sim_levels = levels(samples[[similarity]])                     # The levels in column "similarity"
    for (lvl in sim_levels){                                       # Loop over the levels in the column from option "similarity"
      if (sum(samples[similarity] == lvl, na.rm = T) > longest ){  # Determine the length of the most occurring level from column "similarity"
        longest = sum(samples[similarity] == lvl, na.rm = T)       # Remember the longest length
      }
    }
    n_levels = length(sim_levels)                                  # Calculate the amount of levels in "similarity"
    scIdx = matrix(nrow = n_levels, ncol = longest)                # Set up an empty matrix to output into
    count = 1                                                      #
    for (lvl in sim_levels){                                       # Loop over the levels in the column from option "similarity"
      lvl_i = which(samples[similarity]==lvl)                      # The indices for all values matching the level
      lvl_len = sum(samples[similarity]==lvl, na.rm = T)           # The number of values matching the level
      scIdx[count,] = c(lvl_i, seq(from=-1, to=-1,                 # Create the row from the indices, appending -1 values if the row is
                                   length.out =longest-lvl_len))   # shorter than the longest row
      count = count+1                                              # Keep count of the rows
    }
    row.names(scIdx) = sim_levels                                  # Label the matrix (for user reference only)
    return(scIdx)
  }
  
  count2TPM = function(counts, length){
    # Use a matrix with counts and a matrix with transcript lengths to calculate TPM
    RPK = counts / length              # calculate RPK for all genes in all samples
    TPM = RPK                          # keep the same column and row names as the input matrices
    for (col in 1:NCOL(RPK)){
      colsum = sum(RPK[,col])/1000000  # calculate total sample RPK, divide by one million to calculate the scaling factor
      TPM[,col] = TPM[,col]/colsum     # divide the values in a sample by the scaling factor
    }
    return(TPM)
  }
  ################################################################################################################################################
  ## PRE-CHECKS
  # Check if all requirements are met to generate tx2gene
  if (!isFALSE(make_tx2gene) & is.null(annotation) ){
    stop("When making tx2gene, you must specify a GTF or GFF3 genome annotation file.")
  }
  # Check if returnMode is set to an allowed options, or a character string with allowed options
  if ( FALSE %in% sapply(returnMode, function(x) x %in% c("all","counts","TPMs","vst","dds")) ){
    stop("returnMode may only be one or multiple of: \"all\",\"counts\",\"TPMs\",\"vst\", \"dds\".")
  }
  # Check if all requirments are met to pass to DESeq2
  if ( isTRUE(is.null(data_table) | is.null(salmon_dir) | is.null(DESeq_design) | is.null(control)) & isFALSE(make_tx2gene) ){
    stop("Please specify values for the options \"data_table\", \"salmon_dir\", \"DESeq_design\", \"control\".")
  }
  ################################################################################################################################################
  ## TX2GENE CREATION
  # Create a database for converting transcript to gene ids
  if ( !isFALSE(make_tx2gene)){                                                                   
    txdb = makeTxDbFromGFF(annotation, dataSource=deparse(substitute(annotation)),          # Make TxDb from the annotation file
                           organism="Homo sapiens", taxonomyId=NA , chrominfo=NULL,         # Add metadata
                           circ_seqs=DEFAULT_CIRC_SEQS, miRBaseBuild=NA, format=c("auto"))  # Automatically choose GFF or GTF formatting
    k = keys(txdb, keytype = "TXNAME")                                                      # Load a list with transcript ids
    tx2gene = select(txdb, k, "GENEID", "TXNAME")                                           # A matrix with tx-ids and gene-ids in column one and two
    assign("tx2gene", tx2gene, envir = .GlobalEnv)                                          # Save to the global environment (to save time later)
  }
  # Stop if after making tx2gene if no other options are provided
  if ( is.null(data_table) & is.null(salmon_dir) & is.null(DESeq_design) & is.null(control) & !isFALSE(make_tx2gene) ){
    return(paste0("Stopped after creating 'tx2gene' due to 'data_table', 'salmon_dir', 'DESeq_design', and 'control' not being set.
                 Please ignore the error message as the operation should have been completed normally,
                 and there is no clean way to gracefully exit from an R script"))
  }
  ################################################################################################################################################
  ## CHECK IF DATA_TABLE MEETS REQUIREMENTS
  # Load data_table
  samples = read.table(file.path(data_table), header = TRUE)
  # Check if there is a column called path in "data_table"
  if (isFALSE("path" %in% colnames(samples))){ 
    stop(paste0("data_table \"",data_table,"\" must have a column with \"path\" (without quotes) as a header that contains subfolders of salmon_dir."))
  }
  # Check is the option "DESeq_design" is a column of data_table
  if (isFALSE(DESeq_design %in% colnames(samples))){
    stop(paste0("Option DESeq_design \"",DESeq_design,"\" must refer to a column of data_table \"", data_table,"\"."))
  }
  # Check if the option "control" is an entry in the column of data_table indicated in DESeq_design
  if (isFALSE(control %in% samples[[DESeq_design]])){
    stop(paste0("Option control \"",control,"\" must be a value present in column \"",DESeq_design,"\" of \"", data_table,"\"."))
  }
  ################################################################################################################################################
  ## LOAD AND PROCESS SALMON QUANT FILES
  files = file.path(salmon_dir, samples$path , "quant.sf")          # Load file paths
  names(files) = paste0(row.names(samples))                         # Link sample names to file paths
  txi = tximport(files, type = "salmon", tx2gene = tx2gene)         # Import quant.sf files and sum to gene
  if (isFALSE(batch_correction)){
    # Save raw counts
    if ("counts" %in% returnMode | "all" %in% returnMode ){
      assign("counts", txi$counts, envir = .GlobalEnv)                # Save counts as an object in the global environment
      message("Saved counts the object 'counts'")
    }
    # Save TPMs
    if ("TPMs" %in% returnMode | "all" %in% returnMode ){
      assign("TPMs", txi$abundance, envir = .GlobalEnv)               # Save TPMs as an object in the global environment
      message("Saved TPMs in the object 'TPMs'")
    }
  }
  ################################################################################################################################################
  ## BATCH CORRECTION ON COUNTS
  if (!isFALSE(batch_correction)){
    scIdx = RUVs_rep_matrix(meta_table = data_table,                            # Create a matrix of replicates based on column "batch group"
                          similarity = batch_group)                             # from "data_table"
    message("Starting batch correction using below similarity matrix")
    print(scIdx)                                                                # Print the replicates matrix
    m = as.matrix(round(txi$counts))                                            # Get counts
    genes = rownames(m)                                                         # Take all genenames as control
    low_var = matrix(as.integer(RUVs(m, genes, k=k,scIdx)$normalizedCounts),    # 
                     nrow = nrow(m), ncol = ncol(m))                            #
    low_var[is.na(low_var)] = 0                                                 #
    row.names(low_var) = row.names(m)                                           #
    colnames(low_var) = colnames(m)                                             #
    message("Finished batch correction")
    # Save corrected raw counts
    if ("counts" %in% returnMode | "all" %in% returnMode ){
      assign("counts", low_var, envir = .GlobalEnv)                             # Save counts as an object in the global environment
      message("Saved corrected counts in the object 'counts'")                  #
    }
    # Generate and save corrected TPMs 
    if ("TPMs" %in% returnMode | "all" %in% returnMode ){
      message("Converting batch corrected counts to TPM")
      TPM_low_var = count2TPM(counts = low_var ,length = txi$length)            # Convert batch corrected counts to TPMs
      assign("TPMs", TPM_low_var, envir = .GlobalEnv)                           # Save TPMs as an object in the global environment
      message("Saved TPMs in the object 'TPMs'")
    }
    samples$path = NULL                                                         # Remove the path column this is only used for reading files
    samples[[DESeq_design]] = relevel(samples[[DESeq_design]], control)         # Define which samples are the control
    dds = DESeqDataSetFromMatrix(low_var, samples, 
                                   design = as.formula(paste("~",DESeq_design)))#
    dds = dds[rowMeans2(counts(dds)) >= 10,]                                    #
    # Perform a variance stabilising tranformation and save vst (for visualisation)
    if ("vst" %in% returnMode | "all" %in% returnMode ){
      vst = vst(dds, blind = FALSE)                                             # Perform a variance stabilising transformation
      assign("vst", vst, envir = .GlobalEnv)                                    # Save vst as an object in the global environment
    }
    # Run the DESeq2 algorithms and save dds (for analysis)
    if ("dds" %in% returnMode | "all" %in% returnMode ){
      dds = DESeq(dds)                                                          # run the DESeq algorithm
      assign("dds", dds, envir = .GlobalEnv)                                    # Save dds as an object in the global environment
    }
  }
  ################################################################################################################################################
  ## PROCESS COUNT MATRICES WITH DESEQ2
  if (isFALSE(batch_correction)){
  samples$path = NULL                                                         # Remove the path column this is only used for reading files
  samples[[DESeq_design]] = relevel(samples[[DESeq_design]], control)         # Define which samples are the control
  dds = DESeqDataSetFromTximport(txi, samples, 
                                 design = as.formula(paste("~",DESeq_design)))# Make DESeqDataSet from counts, with DESeq_design as design formula
  dds = dds[rowMeans2(counts(dds)) >= 10,]                                    # Filter out genes with a mean of < 10 counts across samples
  # Perform a variance stabilising tranformation and save vst (for visualisation)
  if ("vst" %in% returnMode | "all" %in% returnMode ){
    vst = vst(dds, blind = FALSE)                                             # Perform a variance stabilising transformation
    assign("vst", vst, envir = .GlobalEnv)                                    # Save vst as an object in the global environment
    message("Saved variance stabilised counts in the object 'vst'")
  }
  # Run the DESeq2 algorithms and save dds (for analysis)
  if ("dds" %in% returnMode | "all" %in% returnMode ){
    dds = DESeq(dds)                                                          # run the DESeq algorithm
    assign("dds", dds, envir = .GlobalEnv)                                    # Save dds as an object in the global environment
    message("Saved the DESeqDataSet in the object 'dds'")
  }
  }
}
