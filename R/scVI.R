
#' Load scVI output from adata object into SingleCellExperiment
#'
#' @param adata AnnData object containing scVI output
#' @param sce Existing sce object. If supplied, coldata is added to this object and scVI normalized expression is added as an altExp
#' @param obs_ignore Which columns in adata$obs should be ignored before importing?
#' @param obs_include Which columns in adata$obs should be included (default "all")
#' @param obs_overwrite Should columns existing in the sce be overwritten?
#' @param obs_prefix Prefix obs columns before importing
#' @param barcodes Function to handle if different barcodes are included in either SCE or adata objects? Intersect to include only barcodes present in both objects, union to include all barcodes - will fill missing values with NA and missing expression with 0 (default: intersect)
#' @param RNA_layer Name of adata$layers that contains the scVI normalized expression matrix
#' @param altexp_name Name of altExp to place scVI normalized SCE object if existing sce is supplied
#' @param reducedDim_include Which reducedDims should be transferred from adata$obsm
#' @param reducedDim_prefix Prefix reducedDim names before adding to sce
#' @param reducedDimNames_replace String to replace in obsm names before adding to the sce (removes "^X_" by default)
#' @param rownames_prefix Prefix for rownames in scVI normalized SCE
#'
#' @return Returns a SingleCellExperiment (SCE). If existing sce is not given the returned SCE contains scVI normalized expression, else scVI normalized expression is added as an altExp
#' @export
#'
scVI_load <- function(adata, sce=NULL, obs_ignore=c("_scvi_labels", "_scvi_batch"), obs_include="all", obs_overwrite=FALSE, obs_prefix="", barcodes=intersect, reducedDim_prefix="", reducedDim_include=c("X_scVI", "X_umap"), reducedDimNames_replace="^X_", RNA_layer="scvi_normalized", altexp_name="VI_RNA", rownames_prefix="VI_"){

  # Get scVI normalized expression into an SCE
  if(RNA_layer %in% names(adata$layers)){
    RNA_data <- adata$layers[[RNA_layer]]
  } else if(RNA_layer %in% names(adata$obsm)) {
    RNA_data <- adata$obsm[[RNA_layer]]
  } else {
    stop(paste("Could not find",RNA_layer,"in layers nor obsm"))
  }

  if(is.null(sce)){
    barcodes_use <- rownames(adata)
  } else {
    barcodes_use <- barcodes(colnames(sce), colnames(rownames(adata)))

    if(length(setdiff(barcodes_use, colnames(sce))) > 0){
      stop("Barcodes included in adata object that cannot be found in SCE. This is currently not supported.")
    }
    sce <- sce[, barcodes_use]
  }

  sce_VI <- SingleCellExperiment::SingleCellExperiment(
    list(
      counts=scutility::subset_matrix(t(RNA_data), barcodes=barcodes_use),
      logcounts=subset_matrix(log1p(t(RNA_data)), barcodes=barcodes_use)
    )
  )

  if(rownames_prefix != "") rownames(sce_VI) <- paste0(rownames_prefix, rownames(sce_VI))

  # If no existing sce is added, make the sce_VI the main object, else add as altExp
  if(is.null(sce)){
    sce <- sce_VI
  } else {
    altExp(sce, altexp_name) <- sce_VI[, barcodes_use]
  }

  reducedDim_add <- lapply(setNames(reducedDim_include,reducedDim_include), FUN=function(name){
    add <- as.data.frame(adata$obsm[[name]])
    rownames(add) <- rownames(adata)

    add <- add[barcodes_use, ]

    return(add)
  })

  names(reducedDim_add) <- gsub(reducedDimNames_replace, "", names(reducedDim_add))

  names(reducedDim_add) <- paste0(reducedDim_prefix, names(reducedDim_add))

  if(length(SingleCellExperiment::reducedDimNames(sce)) == 0){
    SingleCellExperiment::reducedDims(sce) <- reducedDim_add
  } else {
    SingleCellExperiment::reducedDims(sce) <- append(SingleCellExperiment::reducedDims(sce), reducedDim_add)
  }

  # Add colData (stored in adata$obs)
  if(obs_include[1] != "all"){
    adata$obs <- adata$obs[, intersect(colnames(adata$obs), obs_include), drop=FALSE]
  }

  ## prefix obs
  colnames(adata$obs) <- paste0(obs_prefix, colnames(adata$obs))

  if(obs_overwrite == TRUE){
    cols_add <- setdiff(colnames(adata$obs), obs_ignore)
    cols_keep <- setdiff(colnames(SingleCellExperiment::colData(sce)), cols_add)
  } else {
    cols_keep <- colnames(colData(sce))
    cols_add <- setdiff(setdiff(colnames(adata$obs), obs_ignore), cols_keep)
  }

  # Combine colData DataFrame
  colData(sce) <- cbind(colData(sce)[, cols_keep], adata$obs[barcodes_use, cols_add])

  return(sce)
}


#' Load totalVI output from adata object into SingleCellExperiment
#'
#' @param adata AnnData object containing scVI output
#' @param sce Existing sce object. If supplied, coldata is added to this object and scVI normalized expression is added as an altExp
#' @param reducedDims Which reducedDims should be transferred from adata$obsm
#' @param RNA_layer Name of adata$layers that contains the scVI normalized expression matrix
#' @param barcodes Function to handle if different barcodes are included in either SCE or adata objects? Intersect to include only barcodes present in both objects, union to include all barcodes - will fill missing values with NA and missing expression with 0 (default: intersect)
#' @param protein_obsm Name of adata$obsm that contains the totalVI protein denoised expression matrix
#' @param protein_fg_obsm Name of adata$obsm that contains the totalVI protein foreground probability matrix
#' @param protein_altexp_name Name of altExp to place totalVI protein data
#' @param protein_rownames_prefix Prefix for rownames in altExp containing totalVI protein data
#' @param ... Passed on to scVI_load
#'
#' @return Returns a SingleCellExperiment (SCE). If existing sce is not given the returned SCE contains totalVI normalized expression, else totalVI normalized expression and totalVI protein expression are added as altExps
#' @export
#'
totalVI_load <- function(adata, sce=NULL, reducedDim_include=c("totalVI", "X_umap"), RNA_layer="denoised_rna", barcodes=intersect, protein_obsm="denoised_protein", protein_fg_obsm="protein_fg_prob", protein_altexp_name="VI_ADT", protein_rownames_prefix="VI_", ...){

  # Load totalVI normalized RNA expression and obs data similar to scVI output
  sce <- scVI_load(adata, sce=sce, reducedDim_include=reducedDim_include, RNA_layer=RNA_layer, barcodes=barcodes, ...)

  if(is.null(sce)){
    barcodes_use <- rownames(adata)
  } else {
    barcodes_use <- barcodes(colnames(sce), colnames(rownames(adata)))

    if(length(setdiff(barcodes_use, colnames(sce))) > 0){
      stop("Barcodes included in adata object that cannot be found in SCE. This is currently not supported.")
    }
    sce <- sce[, barcodes_use]
  }

  # Add normalized protein expression as an altExp
  altExp(sce, protein_altexp_name) <- SingleCellExperiment(
    list(
      counts=scutility::subset_matrix(t(adata$obsm[[protein_obsm]]), barcodes=barcodes_use),
      logcounts=scutility::subset_matrix(log1p(t(adata$obsm[[protein_obsm]])), barcodes=barcodes_use),
      fg_prob=scutility::subset_matrix(t(adata$obsm[[protein_fg_obsm]]), barcodes=barcodes_use)
    )
  )

  if(protein_rownames_prefix != "") rownames(altExp(sce, protein_altexp_name)) <- paste0(protein_rownames_prefix, rownames(altExp(sce, protein_altexp_name)))

  return(sce)
}

