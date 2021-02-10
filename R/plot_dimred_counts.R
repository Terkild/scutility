#' DimensionReduction with counts plot
#'
#' Very crude and hardcoded function that calls Seurat::FeaturePlot. Needs a remake for more flexibility.
#'
#' @param object Seurat object
#' @param feature Seurat feature
#' @param split.by Meta.data column in seurat object for splitting
#' @param colors Color palette
#' @param limits Limits for color palette (for counts it makes sense that lower limit is 0; DEFAULT: c(0, NA))
#' @param order Should cells be ordered by expression
#' @param slot Seurat slot
#' @param max.cutoff max.cutoff passed on to FeaturePlot
#' @param pt.size pt.size passed on to FeaturePlot
#' @param legend.key.width Theme settings for width of legend
#' @param legend.key.height Theme settings for height of legend
#' @param legend.text Theme settings for
#' @param combine Should plots be combined or returned as a list of plots?
#' @param ... Passed on to Seurat::FeaturePlot
#'
#' @return ggplot2 object(s)
#' @export
#' @import ggplot2
#' @importFrom Seurat FeaturePlot FetchData
#' @importFrom patchwork wrap_plots

seurat_plot_dimred_counts <- function(object,
                               feature,
                               split.by=NULL,
                               colors=c("#000033","#3377FF","#33AAFF","#33CC33","orange","red"),
                               limits=c(0,NA),
                               order=FALSE,
                               slot="counts",
                               max.cutoff='q95',
                               pt.size=0.75,
                               legend.key.width=unit(5,"mm"),
                               legend.key.height=unit(2,"mm"),
                               legend.text=element_text(angle=45),
                               combine=TRUE,
                               ...){

  if(is.null(split.by)){
    split <- rep(feature, nrow(object))
  } else {
    split <- FetchData(object=object, vars=split.by)[,1]
  }

  plot_func <- function(x){
    Seurat::FeaturePlot(object=object,
                        features=feature,
                        slot=slot,
                        max.cutoff=max.cutoff,
                        order=order,
                        pt.size=pt.size,
                        cells=which(split == x),
                        combine=FALSE,
                        ...)[[1]] +
      ggtitle(x) +
      scale_color_gradientn(colours=colors, limits=limits) +
      theme_get() +
      theme(plot.background=element_blank(),
            panel.background=element_blank(),
            axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            legend.key.width=legend.key.width,
            legend.key.height=legend.key.height,
            legend.position=c(1,0),
            legend.justification=c(1,0),
            legend.background=element_blank(),
            legend.text=legend.text,
            legend.direction="horizontal")
  }

  plots <- lapply(levels(as.factor(split)), plot_func)

  if(combine == TRUE){
    return(patchwork::wrap_plots(plots))
  } else {
    return(plots)
  }
}
