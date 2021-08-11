#' Helper function for calculating cutoff value
cutoff_set <- function(values, cutoff){
  if(is.numeric(cutoff)){
    cutoff <- cutoff
  } else if(substr(cutoff, 1, 1) == "q"){
    cutoff <- quantile(values, probs=(as.numeric(gsub("q","", cutoff))/100))
  } else {
    cutoff <- NULL
  }
}

#' Plot dimensional reduction from SCE
#'
#' @param object SingleCellExperiment object
#' @param colour_by Variable to use for coloring
#' @param by_exprs_values SCE Assay to use
#' @param max.cutoff  Cutoff value for upper limit (can be integer or 'q' value for quantile: 'q95' sets at 95% quantile)
#' @param min.cutoff  Cutoff value for lower limit (see max.cutoff)
#' @param order Should points plotted in order by their value
#' @param decreasing  Should ordering be decreasing
#' @param na.last Should NA values be last (on top)?
#' @param ... Passed on to scater::plotReducedDim
#'
#' @return ggplot2 object
#'
#' @importFrom scater plotReducedDim retrieveCellInfo
#' @import SingleCellExperiment
#' @export
plot_dimred <- function(object, colour_by, by_exprs_values="logcounts", max.cutoff=NA, min.cutoff=NA, order=FALSE, decreasing=FALSE, na.last=FALSE, ...){
  data <- retrieveCellInfo(object, colour_by, exprs_values=by_exprs_values)

  if(!is.na(max.cutoff)) max.cutoff <- cutoff_set(data$value, max.cutoff)
  if(!is.na(min.cutoff)) min.cutoff <- cutoff_set(data$value, min.cutoff)

  if(is.numeric(max.cutoff)){
    data$value[data$value > max.cutoff] <- max.cutoff
  }
  if(is.numeric(min.cutoff)){
    data$value[data$value < min.cutoff] <- min.cutoff
  }

  colData(object)[, colour_by] <- data$value

  if(order == TRUE){
    object <- object[, order(colData(object)[, colour_by], decreasing=decreasing, na.last=na.last)]
  }

  scater::plotReducedDim(object=object, colour_by=colour_by, by_exprs_values=by_exprs_values, ...)
}
