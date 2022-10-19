#' Helper function for calculating cutoff value
cutoff_set <- function(values, cutoff){
  if(is.numeric(cutoff)){
    cutoff <- cutoff
  } else if(substr(cutoff, 1, 1) == "q"){
    cutoff <- quantile(values, probs=(as.numeric(gsub("q","", cutoff))/100), na.rm=TRUE)
  } else {
    cutoff <- NULL
  }
}

#' Plot dimensional reduction from SCE
#'
#' @param object SingleCellExperiment object
#' @param dimred Which dimensional reductions should be used?
#' @param columns Vector of columns in SCE to plot
#' @param colour_by Variable to use for coloring
#' @param text_by Variable to use for labelling
#' @param features_add Fetch additional features to data.frame given to ggplot (for specialized plotting)
#' @param by_exprs_values SCE Assay to use
#' @param max.cutoff  Cutoff value for upper limit (can be integer or 'q' value for quantile: 'q95' sets at 95% quantile)
#' @param min.cutoff  Cutoff value for lower limit (see max.cutoff)
#' @param order Should points plotted in order by their value
#' @param order_function What function should be used for ordering?
#' @param shuffle Should points be plotted in random order
#' @param point_size Size of individual points
#' @param point_alpha Alpha level of each point
#' @param scale_color Color scale to use (viridis by default)
#' @param text_by_dimred_summary_fun Which function should be used to place text_by (default median)
#' @param text_by_size Text size if text_by is set
#' @param text_by_fontface Fontface is text_by is set (default "bold")
#' @param text_by_fill Label fill if text_by is set (alpha("white", 0.75) default)
#' @param text_by_padding Label padding if text_by is set.
#' @param text_by_colour Label color if text_by is set (black by default)
#' @param text_repel Should text labels be repelled from each other?
#' @param text_repel_max_overlaps Max overlaps for ggrepel
#' @param seed Seed for random operations
#' @param coldata_exclude_class By default loads all colData except the columns of the classes included in this vector
#' @param rasterise Should points be rasterised (ggrastr)?
#' @param rasterise_dev What device should be used for rasterisation?
#' @param rasterise_dpi What DPI should be used for rasterisation?
#' @param rasterise_scale What scale should be used for rasterisation?
#' @param ... Passed on to makePerCellDF
#'
#' @return ggplot2 object
#'
#' @importFrom scater makePerCellDF
#' @importFrom ggplot2 guides guide_legend
#' @importFrom ggrastr rasterise
#' @import SingleCellExperiment
#' @export
plot_dimred <- function(object, colour_by, features_add=c(), dimred="UMAP", columns=NULL, by_exprs_values="logcounts", seed=12232,
                        order=FALSE, order_function=function(x){order(x, decreasing=FALSE, na.last=FALSE)}, shuffle=FALSE,
                        point_size=0.5, point_alpha=1, max.cutoff=NA, min.cutoff=NA, scale_color=NULL,
                        text_by=NULL, text_by_colour="black", text_by_dimred_summary_fun=median, text_by_size=NA, text_by_fontface="bold", text_by_padding=unit(.01, units="npc"), text_by_fill=alpha(c("white"),0.75), text_repel=TRUE, text_repel_max_overlaps=Inf,
                        coldata_exclude_class=c("CompressedSplitDFrameList"),
                        rasterise=FALSE, rasterise_dev="cairo", rasterise_dpi=300, rasterise_scale=1, ...){

  features <- unique(append(c(colour_by, text_by),features_add))

  # including all colData unless the column class is tagged to not be included
  colData_remove <- which(unlist(lapply(colData(object), class)) %in% coldata_exclude_class)

  if(length(colData_remove) > 0){
    data <- colData(object)[,-colData_remove, drop=FALSE] %>% as.data.frame()
  } else {
    data <- colData(object) %>% as.data.frame()
  }

  ## Alternatively, only return metadata needed
  #data <- colData(object) %>% .[,intersect(features, colnames(.)), drop=FALSE] %>% as.data.frame()

  find_features <- setdiff(features,colnames(data))
  if(length(find_features) < 1) find_features <- NULL
  # if features are not included in colData, fetch them together with dim reduction
  data %<>% cbind(makePerCellDF(x=object, features=find_features, use.dimred=dimred, use.coldata=FALSE, assay.type=by_exprs_values, ...))

  if(!is.null(columns)) data %<>% .[columns,]

  #print(dim(data))
  # if cutoffs are set, calculate cutoffs
  if(!is.na(max.cutoff)) max.cutoff <- cutoff_set(data[[colour_by]], max.cutoff)
  if(!is.na(min.cutoff)) min.cutoff <- cutoff_set(data[[colour_by]], min.cutoff)

  # set values above or below cutoffs to the cutoff values
  if(is.numeric(max.cutoff)){
    data[[colour_by]][data[[colour_by]] > max.cutoff] <- max.cutoff
  }
  if(is.numeric(min.cutoff)){
    data[[colour_by]][data[[colour_by]] < min.cutoff] <- min.cutoff
  }

  if(order == TRUE){
    data <- data[order_function(data[[colour_by]]),]
  }

  if(shuffle == TRUE){
    set.seed(seed)
    data <- data[sample(nrow(data), nrow(data)),]
  }

  if(!is.null(text_by)){
    data_text <- data %>% group_by(.[[text_by]]) %>% summarize(across(starts_with(dimred), text_by_dimred_summary_fun))
    colnames(data_text)[1] <- text_by

    if(text_repel == FALSE){
      text_func <- geom_label
    } else {
      text_func <- ggrepel::geom_label_repel
    }

    text_add <- text_func(data=data_text,
                                          aes(label=.data[[text_by]]),
                                          color=text_by_colour,
                                          seed=seed,
                                          label.size=text_by_size,
                                          fontface=text_by_fontface,
                                          label.padding=text_by_padding,
                                          fill=text_by_fill,
                                          max.overlaps=text_repel_max_overlaps,
                                          na.rm=TRUE)
  } else {
    text_add <- NULL
  }

  if(is.null(scale_color)){
    if(class(data[[colour_by]]) == "numeric"){
      scale_color <- scale_color_viridis_c(na.value="lightgrey")
    } else {
      if(length(unique(data[[colour_by]])) > 8){
        scale_color <- scale_color_viridis_d(na.value="lightgrey")
      } else {
        scale_color <- scale_color_brewer(palette="Set1", na.value="lightgrey")
      }
    }
  }

  geom <- geom_point(alpha=point_alpha, size=point_size)

  if(rasterise == TRUE){
    geom <- rasterise(geom, dev=rasterise_dev, dpi=rasterise_dpi, scale=rasterise_scale)
  }

  plot <- ggplot(data, aes(x=.data[[paste0(dimred, ".1")]],
                           y=.data[[paste0(dimred, ".2")]],
                           col=.data[[colour_by]])) +
      geom +
      scale_color +
      text_add +
      theme_get()

  return(plot)
}
