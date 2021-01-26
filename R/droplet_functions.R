#' Get debris barcodes based on rank cutoff and called cell barcodes
#'
#' @param barcodes_cells Barcodes for droplet called as cell-containing
#' @param debris_cutoff Cutoff for which droplets are likely to contain cells or debris (can be set to estimated number of cells captured)
#' @param ranking_matrix  Count matrix used for ranking droplets
#'
#' @return vector of barcodes
#' @export
#' @import dplyr

droplet_barcodes_debris <- function(barcodes_cells, debris_cutoff, ranking_matrix){
  barcodes_debris <- data.frame(DropletUtils::barcodeRanks(ranking_matrix)) %>%
    mutate(BC=rownames(.)) %>%
    arrange(rank) %>%
    filter(rank<=debris_cutoff & !(BC %in% barcodes_cells)) %>%
    select(BC) %>% .[[1]]

  return(barcodes_debris)
}

#' Droplet counts per compartment
#'
#' Calculate count sums for subsets of columns in a matrix
#'
#' @param matrix Count matrix
#' @param compartment A named list of compartments and their barcodes (i.e. list("cells"=c("BC1", "BC2, ...), "debris"=c("BC101", "BC222, ...)))
#' @param sum_function Function applied to each group of barcodes
#' @param remaining_include Should barcodes that are not included in compartment list be included as its own compartment (i.e. for empty droplets)
#' @param remaining_name If remaining barcodes are included as a compartment, what should it be called?
#'
#' @return Vector or matrix depending on the sum function
#' @export
droplet_counts_per_compartment <- function(matrix, compartment, sum_function=sum, remaining_include=TRUE, remaining_name="empty"){
  compartment_names <- names(compartment)
  if(is.null(compartment_names)) compartment_names <- seq_along(compartment)

  sums <- list()
  for(i in seq_along(compartment)){
    sums[[fraction_names[i]]] <- sum_function(matrix[,intersect(compartment[[i]],colnames(matrix))])
  }
  if(remaining_include == TRUE){
    sums[[remaining_name]] <- sum_function(matrix[,setdiff(colnames(matrix), unlist(compartment))])
  }
  return(sums)
}
