#' Encircle using KDE contours
#'
#' Calculate kernel density estimates and extract contour lines at a given level to make flexible encircling of points.
#'
#' Modified from: https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in
#'
#' @param data  matrix or data frame with x and y values as the first two columns
#' @param contour_level Character string for which contour level to use (Default: "3%")
#' @param contour_factor contour_levels only available for whole percentages. If an intermediate level or a level below "1%" is needed, the contour level value is divided by this factor.
#' @param bandwidth Bandwidth passeds on the ks::kde
#' @param ... Passed on toe ks:kde
#'
#' @return data.frame of contour points
#' @importFrom ks kde
#' @importFrom grDevices contourLines
#' @export

dimred_encircle <- function(x, y, ..., contour_level="3%", contour_factor=1, bandwidth=c(10,10)){
  kd <- ks::kde(data.frame(x=x, y=y), compute.cont=TRUE, h=bandwidth, ...)

  contour <- with(kd, contourLines(x=eval.points[[1]],
                                   y=eval.points[[2]],
                                   z=estimate,
                                   levels=cont[contour_level]/contour_factor)[[1]])

  contour <- data.frame(contour)

  return(contour)
}

#' Encircle groups using KDE contours
#'
#' For each groupings, calculate kernel density estimates and extract contour lines at a given level to make flexible encircling of points.
#'
#' @param data matrix or data.frame containing dimensions and grouping parameters. All columns that are not given as dim1 or dim2 (see below), will be used in grouping the resulting KDE contour.
#' @param dim1 column position or column name of the first dimension in `data`
#' @param dim2 column position or column name of the second dimension in `data`
#' @param ... Passed on to scutility::dimred_encircle() - CURRENTLY NOT IMPLEMENTED
#'
#' @return data.frame of contour points and grouping variables
#' @import dplyr
#' @export
#'
dimred_encircle_groups <- function(data, dim1=1, dim2=2, ...){

  if(!is.numeric(dim1)) dim1 <- which(colnames(data) == dim1)
  if(!is.numeric(dim2)) dim2 <- which(colnames(data) == dim2)

  colnames(data)[c(dim1, dim2)] <- c("dim1", "dim2")

  # A bit of a hack to allow "..." to be passed on the dimred_encircle function
  encircle <- function(.x, .y, ...){
    cbind(dimred_encircle(x=.x$dim1, y=.x$dim2, ...), .y)
  }

  contours <- data %>% group_by(across(c(-dim1, -dim2))) %>% group_map(encircle, ...)

  contours <- bind_rows(contours)

  return(contours)
}


#' Add encircling to ggplot
#'
#' Encircle by Kernel Density Estimate contours
#'
#' @param data matrix or data.frame containing dimensions and grouping parameters. All columns that are not given as dim1 or dim2 (see below), will be used in grouping the resulting KDE contour.
#' @param dim1 column position or column name of the first dimension in `data`
#' @param dim2 column position or column name of the second dimension in `data`
#' @param color_by column position or column name for coloring the geom_path
#' @param size Line size for geom_polygon
#' @param alpha Alpha for geom_polygon
#' @param show_legend Should legend be drawn?
#' @param rescale Should scales be adjusted to make sure contours are complete? This will override previous scale_x_ and scale_y_.
#' @param ... Passed on to scutility::dimred_encircle_groups()
#'
#' @return list of ggplot elements
#' @export
#' @import ggplot2
#' @importFrom ggnewscale new_scale

dimred_encircle_add <- function(data, dim1=1, dim2=2, color_by=NULL, colors=NULL, size=1, alpha=0, line_alpha=0.5, show_legend=FALSE, rescale=TRUE, ...){
  contours <- dimred_encircle_groups(data=data, dim1=dim1, dim2=dim2, ...)

  if(is.null(color_by)){
    # a bit of a hack to allow having no "color_by".
    contours$color_by <- "ROI"
    color_by <- "color_by"
  }

  aesthetics <- aes(x=x, y=y, color=.data[[color_by]], group=.data[[color_by]], fill=.data[[color_by]])

  if(!is.null(colors) & !is.null(color_by)){
    fill <- scale_fill_manual(values=colors)
    colors <- scale_color_manual(values=alpha(colors, line_alpha))
  } else {
    fill <- NULL
  }

  if(rescale==TRUE){
    rescale_x <- scale_x_continuous(limits=function(x){range(append(x, contours$x))})
    rescale_y <- scale_y_continuous(limits=function(y){range(append(y, contours$y))})
  } else {
    rescale_x <- NULL
    rescale_y <- NULL
  }

  p <- list(ggnewscale::new_scale(c("color", "fill")),
            geom_polygon(data=contours,
                      aesthetics,
                      size=size, alpha=alpha, show.legend=show_legend),
            colors,
            fill,
            expand_limits(x=range(contours$x), y=range(contours$y)),
            rescale_x,
            rescale_y,
            NULL)

  return(p)
}

#' Encircle with Seurat meta.data
#'
#' Encircle by Kernel Density Estimate contours
#'
#' @param object Seurat object
#' @param reduction Name of reduction in Seurat object. Uses last in list by default.
#' @param group.by Name of meta data column used for grouping contours
#' @param split.by Name of meta data column included in the ggplot data.frame (for usage in splitting)
#' @param cells Subset of cells to be used for calculating KDE contours. Cell barcodes or indices.
#'
#' @return list of ggplot elements
#' @export
#' @importFrom Seurat Reductions FetchData
seurat_dimred_encircle_add <- function(object, reduction=NULL, group.by=NULL, split.by=NULL, cells=NULL, ...){
  if(is.null(reduction)) reduction <- last(Reductions(object))

  dims <- colnames(object[[reduction]])[1:2]
  data <- Seurat::FetchData(object, vars=c(dims, group.by, split.by), cells=cells)

  #return(data)
  return(dimred_encircle_add(data, color_by=group.by, ...))
}


