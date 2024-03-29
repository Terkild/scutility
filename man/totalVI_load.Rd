% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scVI.R
\name{totalVI_load}
\alias{totalVI_load}
\title{Load totalVI output from adata object into SingleCellExperiment}
\usage{
totalVI_load(
  adata,
  sce = NULL,
  reducedDim_include = c("totalVI", "X_umap"),
  RNA_layer = "denoised_rna",
  barcodes = intersect,
  protein_obsm = "denoised_protein",
  protein_fg_obsm = "protein_fg_prob",
  protein_altexp_name = "VI_ADT",
  protein_rownames_prefix = "VI_",
  ...
)
}
\arguments{
\item{adata}{AnnData object containing scVI output}

\item{sce}{Existing sce object. If supplied, coldata is added to this object and scVI normalized expression is added as an altExp}

\item{RNA_layer}{Name of adata$layers that contains the scVI normalized expression matrix}

\item{barcodes}{Function to handle if different barcodes are included in either SCE or adata objects? Intersect to include only barcodes present in both objects, union to include all barcodes - will fill missing values with NA and missing expression with 0 (default: intersect)}

\item{protein_obsm}{Name of adata$obsm that contains the totalVI protein denoised expression matrix}

\item{protein_fg_obsm}{Name of adata$obsm that contains the totalVI protein foreground probability matrix}

\item{protein_altexp_name}{Name of altExp to place totalVI protein data}

\item{protein_rownames_prefix}{Prefix for rownames in altExp containing totalVI protein data}

\item{...}{Passed on to scVI_load}

\item{reducedDims}{Which reducedDims should be transferred from adata$obsm}
}
\value{
Returns a SingleCellExperiment (SCE). If existing sce is not given the returned SCE contains totalVI normalized expression, else totalVI normalized expression and totalVI protein expression are added as altExps
}
\description{
Load totalVI output from adata object into SingleCellExperiment
}
