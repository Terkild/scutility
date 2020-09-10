#' Load Cell Ranger data and divide by modality
#'
#' @param path  Path to cellranger output (outs) folder
#' @param modalities  Vector of modalities to be extracted from Cell Ranger output
#' @param folder  Folder containing count matrix (raw_feature_bc_matrix or filtered_feature_bc_matrix)
#' @param hto.pattern Pattern in feature name distinguishing hastag (HTOs) from other antibody derived tags (ADTs)
#' @param gex.listname  Name of dataframe containing gene expression counts (after Seurat::Read10X)
#' @param adt.listname  Name of dataframe containing ADT counts (after Seurat::Read10X)
#'
#' @return list of modality count matrices
#' @importFrom Seurat Read10X
#' @export

modalities_load_cellranger_count <- function(path, modalities=c("RNA","ADT","HTO"), folder="raw_feature_bc_matrix", hto.pattern="^hto", gex.listname="Gene Expression", adt.listname="Antibody Capture"){
  data <- Seurat::Read10X(data.dir=file.path(path,folder))

  modality <- list()

  if("RNA" %in% modalities) modality[["RNA"]] <- ifelse(is.list(data) & gex.listname %in% names(data), data[[gex.listname]], data)

  if(is.list(data) & adt.listname %in% names(data)){

    hto.rows <- grep(hto.pattern,rownames(data[[adt.listname]]))

    if("ADT" %in% modalities) modality[["ADT"]] <- data[[adt.listname]][-hto.rows,]

    if("HTO" %in% modalities & length(hto.rows) > 0){
      modality[["HTO"]] <- data[[adt.listname]][hto.rows,]
    }

  }

  return(modality)
}

#' Load Kallisto data by modality
#'
#' @param paths Paths to kallisto output folders
#' @param modalities  Vector of modalities to be loaded from kallisto output
#' @param barcode_suffix  Suffix added to the end of each cell barcode name (to make it compatible with Cell Ranger output)
#' @param folder  Folder containing count matrix
#'
#' @return list of modality count matrices
#' @export

modalities_load_kallisto <- function(paths, modalities=c("ADT","HTO"), barcode_suffix="", folder="counts_unfiltered"){

  modality <- list()

  for(i in seq_along(paths)){
    modality[[modalities[i]]] <- read_kallisto_data(file.path(paths[[i]],folder), ...)
    colnames(modality[[modalities[i]]]) <- paste0(colnames(modality[[modalities[i]]]), barcode_suffix)
  }

  return(modality)
}

#' Load and reformat kallisto output for loading
#'
#' @param path Path to kallisto output folder
#' @param name  Name of count matrix files
#'
#' @return matrix containing kallisto counts
#' @importFrom Matrix t readMM
#' @export

read_kallisto_data <- function(path, name="cells_x_genes"){
  ## Load mtx and transpose it
  res_mat <- as(Matrix::t(Matrix::readMM(file.path(path,paste0(name,".mtx")))), 'CsparseMatrix')
  ## Attach genes
  rownames(res_mat) <- read.csv(file.path(path,paste0(name,".genes.txt")), sep = '\t', header = F)[,1]
  ## Attach barcodes
  colnames(res_mat) <- read.csv(file.path(path,paste0(name,".barcodes.txt")), header = F, sep = '\t')[,1]

  return(res_mat)
}
