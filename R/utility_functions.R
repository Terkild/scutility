#' Calculate barcode rank
#'
#' @return sorted vector of colnames (barcodes)
#' @importFrom Matrix colSums
#' @export
barcode_rank <- function(data){
  data.colsum <- Matrix::colSums(data)
  data.rank <- names(data.colsum)[order(data.colsum, decreasing=TRUE)]
}

#' Subset a matrix
#'
#' If includeAll is TRUE (default) and if a barcode is not present in the
#' matrix, set its values to 0 to allow incorporation into a common object
#'
#' @return A subsetted sparse Matrix
#' @importFrom Matrix Matrix cbind2
#' @export
subset_matrix <- function(data, features=c(), barcodes=c(), includeAll=TRUE, na.value=0){
  if(length(barcodes) < 1) barcodes <- colnames(data)
  cols <- intersect(colnames(data),barcodes)
  cols.diff <- setdiff(barcodes,colnames(data))

  if(length(features) < 1) features <- rownames(data)
  rows <- intersect(rownames(data),features)
  rows.diff <- setdiff(features,rownames(data))

  newmatrix <- data[,cols]

  if(includeAll == TRUE & length(cols.diff) > 0){
    newmatrix.diff <- Matrix::Matrix(data=na.value,
                                nrow=nrow(data),
                                ncol=length(cols.diff),
                                dimnames=list(rownames(data),cols.diff))

    newmatrix <- Matrix::cbind2(newmatrix, newmatrix.diff)
  }

  if(includeAll == TRUE & length(rows.diff) > 0){
    newmatrix.diff <- Matrix::Matrix(data=na.value,
                                     nrow=length(rows.diff),
                                     ncol=length(barcodes),
                                     dimnames=list(rows.diff,barcodes))

    newmatrix <- Matrix::rbind2(newmatrix, newmatrix.diff)
  }

  return(newmatrix[,barcodes])
}

#' Translate ENSG annotation to gene symbol
#'
#' @return Count Matrix with translated rownames
#' @export
translate_ensg_to_symbol <- function(raw_mtx, t2g.file){
  # LIBRARY SHOULD NOT BE LOADED INSIDE FUNCTION
  library("Matrix")
  t2g <- unique(read.csv(t2g.file, sep = '\t', header=F)[,2:3]) # load t2g file
  t2g <- data.frame(t2g[,2], row.names = t2g[,1])
  gene_sym <- t2g[as.character(rownames(raw_mtx)),1] # get symbols for gene ids

  # Which rows have same gene symbol (but different Ensembl gene id)
  gene_sym.duplicated <- which(gene_sym %in% gene_sym[which(duplicated(gene_sym))])

  # Which genes are have duplicated entries
  gene_sym.duplicated.unique <- unique(gene_sym[gene_sym.duplicated])

  # Make placeholder matrix for duplicate gene symbols
  raw_mtx_dedup <- Matrix(data=0,nrow=length(gene_sym.duplicated.unique),ncol=ncol(raw_mtx))
  rownames(raw_mtx_dedup) <- gene_sym.duplicated.unique
  colnames(raw_mtx_dedup) <- colnames(raw_mtx)

  # Combine counts from genes with same gene symbol (but different Ensembl gene id)
  for(i in seq_along(gene_sym.duplicated)){
    curGene <- gene_sym[gene_sym.duplicated[i]]
    curRow <- gene_sym.duplicated.unique == curGene
    raw_mtx_dedup[curRow,] <- raw_mtx_dedup[curRow,] + raw_mtx[gene_sym.duplicated[i],]
  }

  # Merged combined counts duplicate gene symbol with matrix of unique gene symbol counts
  raw_mtx <- raw_mtx[-gene_sym.duplicated,]
  rownames(raw_mtx) <- gene_sym[-gene_sym.duplicated]
  raw_mtx <- rbind(raw_mtx,raw_mtx_dedup)

  return(raw_mtx)
}


#' Biexponential transformation (inspired by flowJo)
#'
#' @export
#' @import scales
biexp_trans <- function(lim = 5, decade.size = lim){
  trans <- function(x){
    ifelse(x <= lim,
           x,
           lim + decade.size * (suppressWarnings(log(x, 10)) -
                                  log(lim, 10)))
  }
  inv <- function(x) {
    ifelse(x <= lim,
           x,
           10^(((x-lim)/decade.size) + log(lim,10)))
  }
  breaks <- function(x) {
    if (all(x <= lim)) {
      scales::pretty_breaks()(x)
    } else if (all(x > lim)) {
      scales::breaks_log(10)(x)
    } else {
      unique(c(scales::pretty_breaks()(c(x[1],lim)),
               scales::breaks_log(10)(c(lim, x[2]))))
    }
  }
  scales::trans_new(paste0("biexp-",format(lim)), trans, inv, breaks)
}
