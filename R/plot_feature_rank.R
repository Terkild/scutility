#' Feature Rank plot
#'
#' @param value             Value used for ranking and x-axis (i.e. UMI count).
#' @param group             Grouping variable for comparison (i.e. titration).
#' @param group_color       Named vector for coloring groups. Used in density plot.
#' @param group_names       Rename group levels (i.e. can be used to append concentration or other info that is marker specific)
#' @param celltype_group    Celltype for grouping cells (i.e. lineage).
#' @param celltype_group_color Named vector for coloring celltype_group. Used in violin plot fill.
#' @param cell_color_by     Variable for coloring cells (i.e. celltype). Used in "barcode" plot.
#' @param cell_color        Named vector for coloring cells. "Barcode" colors.
#' @param barcode_subsample Subsample up to this number of cells from each celltype_group to display in barcode plot
#' @param barplot_group     Grouping variable for barplots (usually the same as group or celltype_group)
#' @param barplot_stack     Subgrouping variable for "stacks" within barplot (usually the same as cell_color_by or celltype_group)
#' @param barplot_stack_color Named vector for coloring barplot_stack levels
#' @param barplot_show_total Should total sums be displayed within barplot? (most meaningful for non-normalized counts)
#' @param barplot_total_fontface Font face for total in barplot. Default "bold".
#' @param barplot_show_ymax Should maximum y-value tick be shown for barplot? (most meaningful for non-normalized counts)
#' @param barplot_ymax_size Font size for ymax text
#' @param show_density      Should density plot be drawn?
#' @param show_barplot      Should barplot be drawn?
#' @param show_barcode      Should barcode plot be drawn?
#' @param show_violin       Should violin plot be drawn?
#' @param show_squares      Should group squares above barplot be drawn?
#' @param combine           Should plot components be combined (default). If FALSE, a list of plot components will be returned.
#' @param ties_method       ties.method passed on to rank().
#' @param rank_normalized   Should ranks be normalized by group (TRUE) or show as absolute counts (FALSE). Default TRUE.
#' @param threshold_y       If not NULL, a threshold line will be drawn at a the set rank
#' @param threshold_color   Threshold line color
#' @param threshold_size    Threshold line size
#' @param threshold_linetype Threshold line type
#' @param threshold_show_celltype_pct Show percent above threshold in celltype component
#' @param threshold_celltype_pct_size Label size for percent above threshold in celltype component
#' @param ...               Passed on to plot_feature_rank_combine(). Useful for setting relative component sizes.
#'
#' @return ggplot2 object or list of ggplot2 objects (if combine=FALSE)
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %<>%

plot_feature_rank_single <- function(value,
                                     group=NULL,
                                     group_color=NULL,
                                     group_names=NA,
                                     celltype_group=NULL,
                                     celltype_group_color=NULL,
                                     cell_color_by=NULL,
                                     cell_color=NULL,
                                     barcode_subsample=200,
                                     barplot_group=NULL,
                                     barplot_stack=NULL,
                                     barplot_stack_color=NULL,
                                     barplot_show_total=TRUE,
                                     barplot_total_fontface="bold",
                                     barplot_show_ymax=TRUE,
                                     show_density=TRUE,
                                     show_barplot=TRUE,
                                     show_barcode=TRUE,
                                     show_violin=TRUE,
                                     show_squares=TRUE,
                                     combine=TRUE,
                                     ties_method="first",
                                     rank_normalize=TRUE,
                                     threshold_y=NA,
                                     threshold_color=alpha("black",0.5),
                                     threshold_size=1,
                                     threshold_linetype="dashed",
                                     threshold_show_celltype_pct=TRUE,
                                     threshold_celltype_pct_size=(theme_get()$text$size/(1/0.352777778)),
                                     labels=c(),
                                     ...){

  # A bit of a hack to allow only some labels to be changed.
  labels_default <- c("x.rank"="Value", "y.rank"="Rank", "y.density"="Density", "y.barplot"="Sum")
  if(length(labels) > 0) labels_default[names(labels)] <- labels
  labels <- labels_default

  if(is_null(group)) group <- 1
  if(is_null(celltype_group)) celltype_group <- 1
  if(is_null(cell_color_by)){
    cell_color_by <- celltype_group
    cell_color <- celltype_group_color
  }
  if(is_null(barplot_group)){
    barplot_group <- celltype_group
  }
  if(is_null(barplot_stack)){
    barplot_stack <- cell_color_by
    barplot_stack_color <- cell_color
  }


  ## Calculate rank
  plotData <- data.frame(value=value,
                         group=factor(group),
                         celltype_group=celltype_group,
                         cell_color_by=cell_color_by,
                         barplot_group=barplot_group,
                         barplot_stack=barplot_stack) %>%
    group_by(group) %>% mutate(rank=rank(value, ties.method=ties_method))

  if(rank_normalize == TRUE){
    plotData %<>% group_by(group) %>% mutate(rank=rank/n())
  }

  scale_y <- scale_y_continuous(expand=c(0.001, 0, 0.001,0))
  scale_x <- scale_x_continuous(trans="biexp", expand=c(0.001, 0, 0.01,0))

  components <- list()

  group_sum <- plotData %>% group_by(group) %>% summarize(total=sum(value)) %>% arrange(desc(total))

  if(is.na(group_names)) group_names <- waiver()

  if(!is.null(group_color)){
    scale_fill <- scale_fill_manual(values=group_color, labels=group_names)
  } else {
    scale_fill <- NULL
  }

  ### DENSITY PLOT
  if(show_density == TRUE){
    aes_density <- aes(x=value, y=..density.., group=group, fill=group)

    if(is.null(group_color)) aes_density$fill <- NULL

    components[['density']] <- plotData %>% mutate(group=factor(group, levels=group_sum[['group']])) %>% ggplot(aes_density) +
      geom_density(alpha=0.75) +
      scale_x +
      scale_y_continuous(expand=c(0,0,0,0)) +
      scale_fill +
      labs(y=labels["y.density"]) +
      theme_get() +
      theme(legend.position="none", #c(1,1)
            legend.justification=c(1,1),
            legend.background=element_blank(),
            legend.title=element_blank(),
            panel.grid=element_blank(),
            plot.background=element_blank(),
            plot.margin=unit(c(0,0,0,0), "points"),
            panel.border=element_blank(),
            axis.line=element_line(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank())

    if(is.na(labels["y.density"])){
      components[['density']] <- components[['density']] +
        theme(axis.title.y=element_blank())
    }
  }


  ### RANK PLOT
  guide_reverse <- FALSE
  if(group_sum[['group']][1] != levels(plotData$group)[1]) guide_reverse <- TRUE

  components[['rank']] <- plotData %>% mutate(group=factor(group, levels=group_sum[['group']])) %>% ggplot(aes(ymin=rank, ymax=rank, xmin=-Inf, xmax=value, x=value, y=rank, fill=group, group=group)) +
    #geom_line() +
    geom_ribbon(orientation="y", color="black", alpha=0.9) +
    scale_fill +
    scale_x +
    scale_y +
    guides(fill=guide_legend(reverse=guide_reverse)) +
    labs(x=labels["x.rank"], y=labels["y.rank"]) +
    theme_get() +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background=element_rect(fill=alpha("white",0.75)),
          legend.title=element_blank(),
          legend.margin=margin(1,1,1,1,"mm"),
          legend.spacing.y=unit(1, "mm"),
          plot.background=element_blank(),
          plot.margin=unit(c(1,1,1,1), "points"),
          panel.grid.minor=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank())

  ## Remove x-label if set to NULL
  if(is.na(labels["x.rank"])){
    components[['rank']] <- components[['rank']] +
      theme(axis.title.x=element_blank())
  }
  if(is.na(labels["y.rank"])){
    components[['rank']] <- components[['rank']] +
      theme(axis.title.y=element_blank())
  }

  ### GROUP SQUARES
  if(show_squares == TRUE){
    components[['squares']] <- plotData %>% group_by(group) %>% summarize(count=n()) %>%
      ggplot() +
      geom_col(aes(x=1, y=1, fill=group), color="black") +
      facet_grid(cols=vars(group)) +
      scale_fill +
      scale_x_discrete(expand=c(0,0,0,0)) +
      scale_y_continuous(expand=c(0,0,0,0)) +
      theme_void() +
      theme(plot.margin=unit(c(0,0,1,0), "points"),
            panel.spacing=unit(1, "points"),
            legend.position="none",
            strip.text=element_blank(),
            strip.background=element_blank())

  }

  ### BAR PLOT
  if(show_barplot == TRUE){
    barplotData <- plotData %>%
      group_by(group, barplot_group, barplot_stack) %>%
      summarize(value=sum(value))


    if(!is.null(barplot_stack_color)){
      scale_fill <- scale_fill_manual(values=barplot_stack_color)
    } else {
      scale_fill <- NULL
    }

    components[['barplot']] <- ggplot(barplotData) +
      geom_col(aes(y=value, x=barplot_group, fill=barplot_stack)) +
      geom_vline(xintercept=Inf) +
      scale_fill +
      scale_x_discrete(expand=c(0,0,0,0)) +
      scale_y_continuous(expand=c(0,0,0.01,0)) +
      theme_get() +
      theme(legend.position="none",
            strip.text=element_blank(),
            strip.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.spacing=unit(1, "points"),
            panel.border=element_blank(),
            axis.line=element_line(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin=unit(c(1,0,0,0), "points")) +
      facet_grid(cols=vars(group))

    ## Should total counts be displayed within the barplot?
    if(barplot_show_total == TRUE){
      components[['barplot']] <- components[['barplot']] +
        geom_text(data=group_sum, aes(y=Inf, x=((length(unique(barplot_group))+1)/2) ,label=total),
                  hjust=1, angle=90, vjust=0.35, fontface=barplot_total_fontface, size=(theme_get()$text$size/(1/0.352777778)))
    }

    if(!is.na(labels["y.barplot"])){
      # A bit of a hack to get the y-axis label to overlap with the neighbouring plot (to minimize whitespace)
      components[['barplot']] <- components[['barplot']] +
        inset_element(p=ggplot() + annotate("text", x=1, y=1, label=labels["y.barplot"],
                                            angle=90, vjust=-0.3, hjust=0.5, size=(theme_get()$axis.text.y$size/(1/0.352777778))) +
                        theme_void() +
                        coord_fixed(clip='off'),
                      clip=FALSE, left=0, right=0, top=1, bottom=0, align_to="plot", ignore_tag=TRUE)

    }

    ## Should the y-axis maximum value be show (makes it much easier to interpret the plot)
    if(barplot_show_ymax == TRUE){
      barplot_ymax <- barplotData %>% group_by(group, barplot_group) %>% summarize(sum=sum(value)) %>% .[["sum"]] %>% max()

      # A bit of a hack to get the y-axis label to overlap with the neighbouring plot (to minimize whitespace)
      components[['barplot']] <- components[['barplot']] +
        inset_element(p=ggplot() + annotate("text", x=1, y=1, label=barplot_ymax,
                                            hjust=1, vjust=0.5, size=(theme_get()$text$size/(1/0.352777778))) +
                        # Attempt to add "0" at y-axis. But gets too cramped.
                        #annotate("text", x=1, y=1, label="0",
                        #         hjust=0, vjust=0, size=(theme_get()$text$size/(1/0.352777778))) +
                        theme_void() +
                        coord_fixed(clip='off'),
                      clip=FALSE, left=0, right=0, top=1, bottom=1, align_to="plot", ignore_tag=TRUE)
    }

  }


  ### CELLTYPE PLOT
  if(show_violin == TRUE | show_barcode == TRUE){

    components[['celltype']] <- ggplot(plotData, aes(y=rank, x=celltype_group)) +
      scale_y +
      theme_get() +
      theme(legend.position="none",
            plot.background=element_blank(),
            strip.text=element_blank(),
            strip.background=element_blank(),
            panel.grid.major.y=element_line(color="lightgrey"),
            panel.grid.major.x=element_blank(),
            panel.grid.minor=element_blank(),
            panel.spacing=unit(1, "points"),
            axis.text.y=element_blank(),
            axis.title=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin=unit(c(1,1,1,0), "points")) +
      facet_grid(cols=vars(group))

    ## VIOLIN PLOT
    if(show_violin == TRUE){
      if(!is.null(celltype_group_color)){
        scale_fill <- scale_fill_manual(values=celltype_group_color)
      } else {
        scale_fill <- NULL
      }

      components[['celltype']] <- components[['celltype']] +
        geom_violin(aes(fill=celltype_group), alpha=0.25, color=alpha("black", 0.1)) +
        scale_fill
    }

    ## BARCODE PLOT
    if(show_barcode == TRUE){
      if(!is.null(cell_color)){
        scale_color <- scale_color_manual(values=cell_color)
      } else {
        scale_color <- NULL
      }

      if(barcode_subsample > 0) barcodeData <- plotData %>%
          group_by(group, celltype_group) %>% sample_n(min(n(),barcode_subsample))

      components[['celltype']] <- components[['celltype']] +
        geom_point(data=barcodeData, aes(color=cell_color_by), pch="_", size=2, alpha=0.75) +
        scale_color
    }
  }



  if(!is.na(threshold_y)){
    ## DOES NOT MAKE SENSE WITH MULTIPLE GROUPS!
    #threshold_rank <- min(plotData$rank[plotData$value >=threshold_value])
    threshold_rank <- threshold_y

    ## Add threshold line to rank plot component
    components[['rank']] <- components[['rank']] +
      geom_hline(yintercept=threshold_rank, color=threshold_color, linetype=threshold_linetype, size=threshold_size)

    ## Add threshold line to celltype component
    components[['celltype']] <- components[['celltype']] +
      geom_hline(yintercept=threshold_rank, color=threshold_color, linetype=threshold_linetype, size=threshold_size)

    if(threshold_show_celltype_pct == TRUE){
      ## Calculate percent above threshold for each celltype_group
      threshold_celltype_pct <- plotData %>%
        group_by(group, celltype_group) %>%
        mutate(celltype_count=n()) %>%
        filter(rank >= threshold_rank) %>%
        summarize(above_threshold=n()/max(celltype_count))

      ## Add percent above threshold to plot
      components[['celltype']] <- components[['celltype']] +
        geom_text(data=threshold_celltype_pct, aes(y=(max(plotData$rank)*0.99), label=sprintf("%.1f",above_threshold*100)),
                  angle=90, hjust=1, vjust=0.3, size=threshold_celltype_pct_size)
    }
  }

  if(combine == FALSE){
    return(components)
  } else {
    return(plot_feature_rank_combine(components, ...))
  }
}

#' Combine Feature Rank plot
#'
#' Function for combining components of the feature rank plot.
#'
#' @param components    List of components (returned from feature_rank_plot() when combine=FALSE)
#' @param height_density Relative height of density plot
#' @param height_rank   Relative height of rank plot
#' @param width_rank    Relative width of rank plot
#' @param width_violin  Relative width of violin plot
#' @param hjust         Horizontal justification of title. Passed on to cowplot::plot_grid()
#' @param vjust         Vertical justification of title. Passed on to cowplot::plot_grid()
#' @param label_x       Label x-position. Passed on to cowplot::plot_grid()
#' @param ...           Passed on to cowplot::plot_grid()
#'
#' @export
#' @importFrom patchwork wrap_plots plot_annotation

plot_feature_rank_combine <- function(components,
                                      title=NA,
                                      subtitle=NULL,
                                      height_rank=8,
                                      height_density=2,
                                      height_barplot=1.6,
                                      width_rank=7.8,
                                      width_celltype=2.2,
                                      hjust=0.5,
                                      vjust=1,
                                      label_x=0.55,
                                      ...){

  labels <- NULL
  if(is.na(title)) title <- NULL

  patchwork::wrap_plots(components[["density"]],
                        patchwork::wrap_plots(components[["squares"]], components[["barplot"]],
                                              heights=c((height_density-height_barplot),height_barplot), ncol=1),
                        components[["rank"]],
                        components[["celltype"]],
                        ncol=2,
                        widths=c(width_rank,width_celltype),
                        heights=c(height_density, height_rank)) +
    patchwork::plot_annotation(title=title, subtitle=subtitle, theme=theme(plot.subtitle=element_text(hjust=0.5)))
}

#' Feature Rank plot
#'
#' @param value             Value used for ranking and x-axis (i.e. UMI count).
#' @param group             Grouping variable for comparison (i.e. titration).
#' @param group_names       Rename group levels (i.e. can be used to append concentration or other info that is marker specific). Named vector (or list of such if split is set).
#' @param split             Splitting variable (generates individual plot for each group)
#' @param celltype_group    Celltype for grouping cells (i.e. lineage).
#' @param cell_color_by     Variable for coloring cells (i.e. celltype). Used for "barcode" plot.
#' @param threshold_y       If not NULL, a threshold line will be drawn at a the set rank. If a single value, same threshold is set for all splits. When splitting, multiple values can be given in a vector (named by split)
#' @param combine           Should plot components be combined (default). If FALSE, a list of plot components will be returned.
#' @param nrow              Number of rows to use when wrapping
#' @param title             If not NA, title will be set. If a single value title only first panel will get title. When splitting, multiple values can be given in a vector (named by split)
#' @param remove_ylabel     Should y-labels for all but the first panel be removed?
#' @param remove_ylabel_width Decrease relative width of non-first plots by this factor to account for removed y-labels in others
#' @param ...               Passed on to plot_feature_rank_single().
#'
#' @return ggplot2 object or list of ggplot2 objects (if combine=FALSE)
#' @export

plot_feature_rank <- function(value, group=NULL, group_names=NA, split=NULL, celltype_group=NULL, cell_color_by=NULL, combine=TRUE, nrow=1, threshold_y=NA, title=NA, remove_ylabel=FALSE, remove_ylabel_width=0.92, labels=c(), widths=NULL, ...){
  plotData <- data.frame(value=value)

  if(!is.null(group)) plotData$group=group
  if(!is.null(celltype_group)) plotData$celltype_group=celltype_group
  if(!is.null(cell_color_by)) plotData$cell_color_by=cell_color_by

  if(length(split) == length(value)){
    data_list <- plotData %>% mutate(split=as.factor(split)) %>%
      group_by(split) %>% group_split() %>%
      setNames(levels(split))

    labels_list <- lapply(data_list, function(x) labels)
    widths <- rep(1, length(data_list))

    if(remove_ylabel == TRUE){
      labels_list[-1] <- lapply(labels_list[-1], function(x) c(c("y.rank"=NA, "y.density"=NA),labels))

      widths[-1] <- sapply(widths[1], function(x) remove_ylabel_width)
    }

    ## Check if different thresholds are needed for each split
    if(length(threshold_y) != length(data_list)){
      threshold_y <- rep(threshold_y, length(data_list))
    }
    if(is.null(names(threshold_y))) names(threshold_y) <- names(data_list)

    ## If only a single title is set, only use it for first panel
    if(length(title) != length(data_list)){
      title_empty <- rep("", length(data_list))
      title_empty[1] <- title
      title <- title_empty
    }
    if(is.null(names(title))) names(title) <- names(data_list)

    if(!is.list(group_names)){
      group_names <- rep(group_names, length(data_list))
    }
    if(is.null(names(group_names))) names(group_names) <- names(data_list)

    plot_list <- lapply(names(data_list),
                        function(x) with(data_list[[x]],
                                         plot_feature_rank_single(value=value,
                                                                  group=group,
                                                                  group_names=group_names[[x]],
                                                                  celltype_group=celltype_group,
                                                                  cell_color_by=cell_color_by,
                                                                  subtitle=x,
                                                                  threshold_y=threshold_y[x],
                                                                  combine=combine,
                                                                  title=title[x],
                                                                  labels=labels_list[[x]],
                                                                  ...)))

    if(combine == TRUE){
      plot <- cowplot::plot_grid(plotlist=plot_list, nrow=nrow, rel_widths=widths)
    } else {
      plot <- plot_list
    }

  } else {
    plot <- plot_feature_rank_single(value=value,
                                     group=group,
                                     group_names=group_names,
                                     celltype_group=celltype_group,
                                     cell_color_by=cell_color_by,
                                     threshold_y=threshold_y,
                                     combine=combine,
                                     title=title,
                                     labels=labels,
                                     ...)
  }

  return(plot)
}

#' Feature Rank plot Seurat
#'
#' Create feature rank plot from Seurat object
#'
#' @param object            Seurat object
#' @param feature           Meta.data column for ranking and x-axis (i.e. UMI count).
#' @param group.by          Meta.data column for grouping variable for comparison (i.e. titration).
#' @param split.by          Meta.data column for splitting plots
#' @param celltype_group    Meta.data column celltype for grouping cells (i.e. lineage).
#' @param cell_color_by     Meta.data column for coloring cells (i.e. celltype). Used for "barcode" plot.
#' @param slot              Slot that the value should be extracted from (default "counts")
#'
#' @return ggplot2 object or list of ggplot2 objects (if combine=FALSE)
#' @export
#' @importFrom Seurat FetchData

seurat_plot_feature_rank <- function(object,
                                     feature,
                                     group.by=NULL,
                                     split.by=NULL,
                                     celltype_group=NULL,
                                     cell_color_by=NULL,
                                     slot="counts",
                                     ...){
  if(!is.null(split.by)){
    split <-  Seurat::FetchData(object, split.by)[[1]]
  } else {
    split <- NULL
  }

  if(!is.null(cell_color_by)) cell_color_by <- Seurat::FetchData(object, cell_color_by)[[1]]
  if(!is.null(celltype_group)) celltype_group <- Seurat::FetchData(object, celltype_group)[[1]]
  if(!is.null(group.by)) group.by <- Seurat::FetchData(object, group.by)[[1]]

  plot_feature_rank(value=Seurat::FetchData(object, feature, slot=slot)[[1]],
                    group=group.by,
                    split=split,
                    celltype_group=celltype_group,
                    cell_color_by=cell_color_by,
                    ...)
}
