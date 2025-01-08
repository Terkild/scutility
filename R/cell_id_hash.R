#' Generate Truncated SHA-256 Hash for Each Row of Two Columns
#'
#' This function generates a truncated SHA-256 hash for each row by concatenating
#' a truncated version of the `barcode_column` (first 16 characters) with the `run_column`.
#' The resulting hash is truncated to the specified length.
#'
#' @param barcode_column A character vector containing barcodes.
#' @param run_column A character vector containing run identifiers.
#' @param length An integer specifying the length of the truncated hash. Default is 16.
#'
#' @return A character vector of truncated SHA-256 hashes, one for each row.
#' @examples
#' # Example usage
#' barcode_column <- c("ACGTGCTAGCTAGCTA-1", "GTCAGTCAGTCAGTCA-1", "TGCATGCATGCATGCA-1")
#' run_column <- c("Run1", "Run2", "Run3")
#' cell_id_hash(barcode_column, run_column)
#'
#' @import digest
#' @export
cell_id_hash <- function(barcode_column, run_column, length = 16) {
  # Ensure barcode_column and run_column are character vectors
  barcode_column <- as.character(barcode_column)
  run_column <- as.character(run_column)

  # Truncate barcode_column to the first 16 characters
  truncated_barcodes <- substr(barcode_column, 1, 16)

  # Concatenate the truncated barcodes and run_column
  concatenated <- paste0(truncated_barcodes, run_column)

  # Generate SHA-256 hash for each concatenated string
  hashed <- vapply(
    concatenated,
    function(x) substr(digest(x, algo = "sha256", serialize = FALSE), 1, length),
    FUN.VALUE = character(1) # Ensures the output is a character vector
  )

  return(hashed)
}
