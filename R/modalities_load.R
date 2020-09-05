#' Load Cell Ranger data and divide by modality
#'
#' @export
modalities_load_cellranger_count <- function(path, modalities=c("RNA","ADT","HTO"), hto.pattern="^hto", gex.listname="Gene Expression", adt.listname="Antibody Capture"){
  data <- Seurat::Read10X(data.dir=path)

  hto.rows <- grep(hto.pattern,rownames(data[[adt.listname]]))

  modality <- list()

  if("RNA" %in% modalities) modality[["RNA"]] <- data[[gex.listname]]

  if("ADT" %in% modalities) modality[["ADT"]] <- data[[adt.listname]][-hto.rows,]

  if("HTO" %in% modalities & length(hto.rows) > 0){
    modality[["HTO"]] <- data[[adt.listname]][hto.rows,]
  }

  return(modality)
}

#' Load Kallisto data by modality
#'
#' @export
modalities_load_kallisto <- function(paths, modalities=c("GEX","ADT","HTO")){

  modality <- list()

  for(i in seq_along(paths)){
    modality[[modalities[i]]] <- read_kallisto_data(paths[[i]])
  }

  return(modality)
}

#' Load and reformat kallisto output for loading
#'
read_kallisto_data <- function(path, name="cells_x_genes"){
  library("Matrix")
  ## Load mtx and transpose it
  res_mat <- as(Matrix::t(Matrix::readMM(file.path(path,paste0(name,".mtx")))), 'CsparseMatrix')
  ## Attach genes
  rownames(res_mat) <- read.csv(file.path(path,paste0(name,".genes.txt")), sep = '\t', header = F)[,1]
  ## Attach barcodes
  colnames(res_mat) <- read.csv(file.path(path,paste0(name,".barcodes.txt")), header = F, sep = '\t')[,1]

  return(res_mat)
}
