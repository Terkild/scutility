
#' Make labels for groups in DimRed plots
#'
#' Assign labels to positions in reduced dimensional space. Can be used to label clusters by other meta.data than their coloring.
#'
#' @param object Seurat object
#' @param reduction Label of reduction object to use (uses the last created reduction by default)
#' @param dims  Which dimensions should be used?
#' @param group.by Meta data column to label clusters
#' @param FUN   Function called on reduced dimensions to calculate label position
#'
#' @return tibble of labels and their positions
#' @export
#' @importFrom Seurat FetchData Reductions Key
#' @import dplyr
seurat_dimred_cluster_labels <- function(object, reduction=NULL, dims=c(1,2), group.by="ident", FUN=median){
  if(is.null(reduction)) reduction <- last(Seurat::Reductions(object))

  reductions <- paste0(Seurat::Key(object[[reduction]]), dims)

  cluster_labels <- Seurat::FetchData(object, vars=append(reductions, group.by)) %>%
    group_by(!!as.name(group.by)) %>%
    summarize(dim1=FUN(!!as.name(reductions[1])),
              dim2=FUN(!!as.name(reductions[2])))

  return(cluster_labels)
}

#' Add labels to DimPlot (or similar ggplot2)
#'
#' Sets a new color scale for the added annotation and sets up default aesthetics
#'
#' @param labels data.frame of labels and their positions (i.e. created by scutility::seurat_dimred_cluster_labels)
#' @param colors named vector of colors for each label
#' @param label_by    name or index of the column in labels that should be used for labelling
#' @param color_by    name or index of the column in labels that should be used for coloring
#' @param dim1    name or index of the column in labels that should be used for x
#' @param dim2    name or index of the column in labels that should be used for y
#' @param ...   Passed on to ggrepel::geom_label_repel
#'
#' @return list of ggplot2 functions to be added to an existing plot
#' @export
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggnewscale new_scale_color
#' @importFrom ggplot2 scale_color_manual
dimred_labels_add <- function(labels, colors=NULL, label_by=1, color_by=1, dim1=2, dim2=3, seed=112, label.size=NA, fontface="bold", label.padding=unit(.01, units="npc"), na.rm=TRUE, fill=alpha(c("white"),0.75), ...){

  if(is.numeric(label_by)) label_by <- colnames(labels)[label_by]
  if(is.numeric(color_by)) color_by <- colnames(labels)[color_by]
  if(is.numeric(dim1)) dim1 <- colnames(labels)[dim1]
  if(is.numeric(dim2)) dim2 <- colnames(labels)[dim2]

  if(!is.null(colors)){
    colors <- ggplot2::scale_color_manual(values=colors)
  }

  list(ggnewscale::new_scale_color(),
       ggrepel::geom_label_repel(data=labels,
                                 aes(x=.data[[dim1]], y=.data[[dim2]], color=.data[[color_by]], label=.data[[label_by]]),
                                 seed=seed,
                                 label.size=label.size,
                                 fontface=fontface,
                                 label.padding=label.padding,
                                 na.rm=na.rm,
                                 fill=fill, ...),
       colors)
}

#' Add labels to DimPlot with annotation from seurat object
#'
#' Wrapper function for calling seurat_dimred_cluster_labels and dimred_labels_add
#'
#' @param object Seurat object
#' @param reduction Label of reduction object to use (uses the last created reduction by default)
#' @param dims  Which dimensions should be used?
#' @param group.by Meta data column to label clusters
#' @param FUN   Function called on reduced dimensions to calculate label position
#' @param ...   Passed on to scutility::dimred_labels_add
#'
#' @return list of ggplot2 functions to be added to an existing plot
#' @export
seurat_dimred_labels_add <- function(object, reduction=NULL, dims=c(1,2), group.by="ident", FUN=median, ...){
  labels <- seurat_dimred_cluster_labels(object=object, reduction=reduction, dims=dims, group.by=group.by)

  return(dimred_labels_add(labels=labels, ...))
}
