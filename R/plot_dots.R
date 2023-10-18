#' Plot dots
#'
#' Plots population/cluster statistics as percentage positive (dot size) and mean expression within positive population (color)
#'
#' @param sce SingleCellExperiment object containing expression matrices
#' @param coldata DataFrame containing cluster information for each column in the provided matrices. Defaults to colData(sce)
#' @param rows colData column defining clustering to be used as rows
#' @param rows_colors rows can be color annotated by providing a named color palette here
#' @param rows_colors_width Relative (to 1) width of rows color annotation
#' @param rows_distribution_by Plot how cluster is distributed among another colData column
#' @param rows_distribution_fill Colors palette (color vector) for cluster distribution
#' @param rows_distribution_width Relative (to 1) width of rows distribution plot
#' @param colnames_fun Function to clean up rownames (such as trimming prefix)
#' @param assay_color Name of SCE assay to use for coloring dots
#' @param assay_positive Name of SCE assay to use for gating positive cells
#' @param threshold_positive Threshold for gating positive cells
#' @param threshold_color Only positive populations having a mean expression above this threshold will be displayed
#' @param threshold_fraction Minimum fraction of cells to be included in the plot
#' @param feature_annotation Annotation data.frame for features (columns) to be shown as a heatmap
#' @param feature_annotation_height Relative (to 1) height of feature annotation heatmap
#' @param max.cutoff  Cutoff value for upper limit (can be integer or 'q' value for quantile: 'q95' sets at 95\% percentile)
#' @param gaps_y Not currently used
#' @param scale_fill ggplot scale_fill_continuous function for coloring dots
#'
#' @return ggplot2 combined with patchwork
#' @export
plot_dots <- function(sce,
                      rows=SingleCellExperiment::colLabels(sce),
                      rows_colors=NA,
                      rows_colors_width=0.02,
                      rows_distribution_by=NA,
                      rows_distribution_fill=NULL,
                      rows_distribution_width=0.05,
                      colnames_fun=function(x){return(as.factor(x))},
                      coldata=SingleCellExperiment::colData(sce),
                      assay_color="logcounts",
                      assay_positive="fg_prob",
                      threshold_positive=0.95,
                      threshold_color=2,
                      threshold_fraction=0.05,
                      feature_annotation=NULL,
                      feature_annotation_height=0.1,
                      max.cutoff=NA,
                      gaps_y=NA,
                      scale_fill=scale_fill_viridis_c(option="turbo")
){

  # Get matrices and convert to long format
  mtx_color <- assay(sce, assay_color) %>% as.matrix %>% as.data.frame %>%  tibble::rownames_to_column("feature") %>% tidyr::pivot_longer(-feature, names_to="BC", values_to="color_by")
  mtx_positive <- assay(sce, assay_positive) %>% as.matrix %>% as.data.frame %>% tibble::rownames_to_column("feature") %>% tidyr::pivot_longer(-feature, names_to="BC", values_to="positive")


  # Get row information
  if(length(rows) == ncol(sce)){
    data_meta <- data.frame("cluster"=rows)
    rows <- "cluster"
  } else if(rows %in% colnames(coldata)){
    data_meta <- coldata[, rows, drop=FALSE] %>% as.data.frame()
  } else {
    stop("Rows has to be a vector of length equal to ncol(sce) or a column in the coldata")
  }

  data_meta[[rows]] <- as.factor(data_meta[[rows]])

  # Join matrices and metadata
  mtx <- cbind(mtx_color, positive=mtx_positive$positive, data_meta)

  features <- levels(colnames_fun(mtx_color$feature))

  if(!is.na(gaps_y)){
    gaps_y_plot <- geom_hline(data=gaps_y, aes(yintercept=cluster_count_cumsum+0.5), color=alpha("black", 0.5), size=0.5)
  } else {
    gaps_y_plot <- NULL
  }

  scale_y <- NULL#scale_y_discrete(expand=c(0,0,0,0))

  # Calculate metrics
  data_plot <- mtx %>%
    group_by(feature, .data[[rows]], .drop=FALSE) %>% mutate(cluster_count=n()) %>%
    # only include cells considered positive
    filter(positive > threshold_positive) %>%
    # calculate mean expression and percent positive
    group_by(feature, .data[[rows]], cluster_count, .drop=FALSE) %>%
    summarize(positive_count=n(), positive_mean=mean(color_by)) %>%
    mutate(positive_frac=positive_count/cluster_count) %>%
    filter(positive_frac > threshold_fraction & positive_mean > threshold_color)

  # Allow colnames to be changed (ordering, trimming etc.)
  data_plot$feature <- colnames_fun(data_plot$feature)
  data_plot$feature <- factor(data_plot$feature, levels=features)


  # set values above or below cutoffs to the cutoff values
  if(!is.na(max.cutoff)) max.cutoff <- cutoff_set(data[[colour_by]], max.cutoff)

  if(is.numeric(max.cutoff)){
    data_plot[["positive_mean"]][data_plot[["positive_mean"]] > max.cutoff] <- max.cutoff
  }

  # Make plot
  plots <- list()
  plots[["dots"]] <- data_plot %>%
    ggplot(aes(x=feature, y=.data[[rows]])) +
    geom_point(aes(size=positive_frac, fill=positive_mean), shape=21, color=alpha("black", 0.25)) +
    scale_x_discrete(drop=FALSE) +
    gaps_y_plot +
    scale_y +
    scale_fill +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1),
          plot.margin=margin(0,0,0,0, "pt"),
          panel.grid.major=element_line(color=alpha("black", 0.05)))

  # Draw distribution plot for rows
  if(!is.na(rows_distribution_by)){
    if(length(rows_distribution_by) == ncol(sce)){

      data_distribution <- cbind(data_meta[[rows]], rows_distribution_by)
      rows_distribution_by <- "distribution"

    } else if(rows %in% colnames(coldata)){
      data_distribution <- cbind(data_meta[[rows]],
                                 coldata[,c(rows_distribution_by), drop=FALSE] %>% as.data.frame())
    } else {
      stop("Data_distribution_by has to be: NA, a vector of length equal to ncol(sce) or a column in the coldata")
    }

    colnames(data_distribution) <- c(rows, rows_distribution_by)

    if(!is.null(rows_distribution_fill)){
      rows_distribution_scale_fill <- scale_fill_manual(values=rows_distribution_fill)
    } else {
      rows_distribution_scale_fill <- NULL
    }

    plots[["rows_distribution"]] <- data_distribution %>%
      group_by(.data[[rows]]) %>% mutate(cluster_count=n()) %>%
      group_by(.data[[rows]], cluster_count, .data[[rows_distribution_by]]) %>% summarize(count=n()) %>%
      mutate(freq=count/cluster_count) %>%
      ggplot(aes(x=.data[[rows]], y=freq, fill=.data[[rows_distribution_by]])) +
      geom_col(color=alpha("black", 0.25), width=1) +
      scale_y +
      gaps_y_plot +
      rows_distribution_scale_fill +
      coord_flip() +
      #scale_fill_manual(values=metadata(sce)$colors$tissue) +
      theme(axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line=element_blank(),
            panel.border=element_blank(),
            plot.margin=margin(0,0,0,0, "pt"),
            legend.position="right")

  }

  # Draw color legend for rows
  if(!is.na(rows_colors[1]) & all(names(rows_colors) %in% levels(data_plot[[rows]]))){

    data_rows <- data.frame(cluster=levels(as.factor(data_meta[[1]])) %>% forcats::fct_inorder())

    plots[["rows_colors"]] <- data_rows %>%
      ggplot(aes(x=1, y=cluster, fill=cluster)) +
      geom_tile(color=alpha("black", 0.5)) +
      scale_y +
      gaps_y_plot +
      scale_fill_manual(values=rows_colors) +
      theme(axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(),
            axis.ticks.x=element_blank(),
            axis.line=element_blank(),
            panel.border=element_blank(),
            plot.margin=margin(0,0,0,0, "pt"),
            legend.position="none")
  }



  # Draw feature annotation heatmap
  if(!is.null(feature_annotation)){
    feature_annotation$feature <- colnames_fun(feature_annotation$feature)

    # Make sure the same columns are displayed
    feature_annotation <- filter(feature_annotation, feature %in% features)
    feature_annotation$feature <- factor(feature_annotation$feature, levels=features)

    plots[["feature_annotation"]] <- ggplot(feature_annotation, aes(x=feature, y=name, fill=value)) +
      geom_tile(color=alpha("black", 0.5)) +
      scale_y_discrete(position = "right", expand=c(0,0,0,0)) +
      scale_x_discrete(drop=FALSE) +
      scale_fill_viridis_c(option="turbo") +
      theme(axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(),
            axis.ticks.x=element_blank(),
            axis.line=element_blank(),
            panel.border=element_blank(),
            plot.margin=margin(0,0,0,0, "pt"),
            legend.position="none")
  }


  spacer <- plot_spacer() + theme(plot.margin=margin(0,0,0,0, "pt"))

  if(is.null(plots[["feature_annotation"]])){
    plots[["feature_annotation"]] <- spacer
    feature_annotation_height <- 0
  }

  if(is.null(plots[["rows_distribution"]])){
    plots[["rows_distribution"]] <- spacer
    rows_distribution_width <- 0
  } else {
    # remove axis labels if "annotation" axes are used
    plots[["dots"]] <- plots[["dots"]] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(is.null(plots[["rows_colors"]])){
    plots[["rows_colors"]] <- spacer
    rows_colors_width <- 0
  } else {
    # remove axis labels if "annotation" axes are used
    plots[["dots"]] <- plots[["dots"]] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    plots[["rows_distribution"]] <- plots[["rows_distribution"]] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  plot <- patchwork::wrap_plots(spacer,spacer,plots[["feature_annotation"]],
                                plots[["rows_colors"]],plots[["rows_distribution"]],plots[["dots"]],
                                widths=c(rows_colors_width,rows_distribution_width,1),
                                heights=c(feature_annotation_height, 1),
                                guides="collect")

  return(plot)
}
