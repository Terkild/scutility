% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtx_balance.R
\name{mtx_balance}
\alias{mtx_balance}
\title{Balance/Upsample mtx}
\usage{
mtx_balance(
  path,
  n_samples,
  path_save = NULL,
  file_matrix = "matrix.mtx",
  file_barcodes = "barcodes.tsv",
  file_features = "genes.tsv",
  matrix_sep = " ",
  matrix_colnames = c("feature", "barcode", "count"),
  expected_counts = c(),
  ...
)
}
\description{
Balance/Upsample mtx
}
