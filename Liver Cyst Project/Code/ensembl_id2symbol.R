ensembl_id2symbol = function(x, keepids = F)
{
  ## A function that changes the rownames of a data frame from ENSEMBL ids to gene symbols
  ##
  ## Arguments
  ## x                                      Matrix that has Ensembl ids as row names 

  library("org.Hs.eg.db")
  library("AnnotationDbi")
  
  row.names(x) = sapply(strsplit(row.names(x), ".", fixed=T), function(x) x[1])         # Split rownames to turn gencode ids into ENSMBL ids
  if (isTRUE(keepids)){
    x$symbol = mapIds(org.Hs.eg.db, keys=row.names(x), column="SYMBOL",                 # Make an extra column with the gene sybmol matching the ENSMBL id
                      keytype="ENSEMBL", multiVals="first")                             # 
    return(x)
  }
  m = as.matrix(mapIds(org.Hs.eg.db, keys=row.names(x), column="SYMBOL",                # Map ENSMBL ids to gene symbols
                       keytype="ENSEMBL", multiVals="first"))                           #
  m[which(is.na(m))] = rownames(m)[which(is.na(m))]                                     # Keep ids without a symbol as ids
  row.names(x) = as.vector(m)                                                           # Set the new rownames as gene symbols
  return(x)
}
