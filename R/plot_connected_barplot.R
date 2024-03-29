#' Plot connected barplot from seurat object
#'
#' @param object Seurat object
#' @param group.by Meta data column containing population/cluster assignment (y-axis)
#' @param split.by Meta data column containing group or sample assignment (x-axis)
#' @param wrap.by Meta data column to split into multiple independent plots (i.e. patient or tissue information)
#' @param normalize_to  Data frame containing normalization factor (i.e. total cell count from sample) as first column, split.to (and wrap.by) in following columns. Important that colnames match split.to and wrap.by names
#' @param cells Vector of rows or cell barcodes to be included in plot.
#' @param combine Boolean of whether multiple plots (when wrap.by is set) should be combined into a single plot (using cowplot)
#' @param wrap_add  Integer of how much width should be allocated to y-axis text relative to the with of a bar (for combining plots)
#'
#' @returns ggplot object or list of ggplot objects
#' @importFrom Seurat FetchData
#' @export

seurat_plot_connected_barplot <- function(object, group.by="ident", split.by, wrap.by=NULL, normalize_to=NULL, cells=c(), combine=TRUE, wrap_add=0.5, ...){
  getData <- Seurat::FetchData(object, vars=c(group.by, split.by, wrap.by))
  colnames(getData)[1:2] <- c("group.by","split.by")

  if(length(cells)>0){
    getData <- getData[cells,]
  }

  if(!is.null(wrap.by)){
    colnames(getData)[3] <- "wrap.by"

    dataList <- split(getData, getData$wrap.by)
    plots <- lapply(dataList, function(x){
        ## If normalize_to is set, extract data frame related to current "wrap"
        if(!is.null(normalize_to)) normalize_to <- normalize_to[normalize_to[[wrap.by]] == x$wrap.by[1], c(1,2)]
        plot_connected_barplot(population=x$group.by, group=x$split.by, normalize_to=normalize_to, ...) + ggtitle(x$wrap.by[1])
      })

    if(combine == TRUE){
      wrapSize <- getData %>%
        group_by(wrap.by, split.by) %>%
        summarize(splitSize=n()) %>%
        group_by(wrap.by) %>%
        summarize(wrapSize=n())

      plots <- cowplot::plot_grid(plotlist=plots[wrapSize$wrap.by],
                         align="hv", axis="tblr", nrow=1,
                         rel_widths=(wrapSize$wrapSize+wrap_add))
    }
  } else {
    plots <- plot_connected_barplot(population=getData$group.by, group=getData$split.by, normalize_to=normalize_to, ...)
  }

  return(plots)
}

#' Plot connected barplot
#'
#' @param population  vector containing population/cluster assignment (y-axis)
#' @param group vector containing group or sample assignment (x-axis)
#' @param connected Should bars from same group be connected?
#' @param bar_width Relative width of bars (integer between 0 and 1)
#' @param normalize_to Data frame containing normalization factor (i.e. total cell count from sample) as first column and group in the second column.
#' @param y_value what statistic to plot ("percent" or "count")
#' @param order should populations be ordered by total count (across groups)
#' @param colors (named) vector of colors for populations
#'
#' @returns ggplot object
#' @import ggplot2
#' @export

plot_connected_barplot <- function(population, group, connected=TRUE, bar_width=0.4, normalize_to=NULL, y_value="percent", order=FALSE, colors=c(), label=FALSE){
  getData <- data.frame(group=group, population=population)

  if(class(group) != "factor") getData$group <- as.factor(getData$group)

  plotData <- getData %>%
    group_by(group) %>%
    mutate(groupCount=n()) %>%
    group_by(group, groupCount, population) %>%
    summarize(populationCount=n()) %>%
    ## Arrange by total counts within a population across groups
    group_by(population) %>%
    arrange(desc(populationCount))

  ## Vector of population names ordered by total counts within a population across groups
  populations <- unique(plotData$population)
  population_n <- length(populations)

  ## If colors are insufficient get a distinct set of colors
  if(length(colors) < population_n){
    colors <- scutility::colors_get_distinct()
  }

  ## If not already assigned, assign colors to population
  if(length(setdiff(populations,names(colors)))>0){
    colors <- colors[1:population_n]
    names(colors) <- populations
  }
  population_order <- names(colors)

  ## If set, order populations based on total counts within a population across groups
  if(order == TRUE){
    population_order <- populations
  }

  plotData <- plotData %>%
    mutate(population=factor(population, levels=population_order),
           group=group)


  if(y_value == "percent"){
    if(!is.null(normalize_to)){
      colnames(normalize_to) <- c("norm", "group")

      plotData <- plotData %>%
        left_join(normalize_to, by="group") %>%
        mutate(pct=populationCount/groupCount*norm) %>%
        mutate(value=pct)
    } else {
      plotData <- plotData %>%
        mutate(pct=populationCount/groupCount*100) %>%
        mutate(value=pct)
    }

    y_label <- "%"
  } else {
    plotData <- plotData %>%
      mutate(value=populationCount)

    y_label <- "Cells"
  }


  plot <- ggplot(plotData, aes(x=group, y=value, fill=population, alluvium=population, stratum=population)) +
    ggalluvial::geom_stratum(color=alpha("black",0.5), width=bar_width) +
    scale_fill_manual(values=colors) +
    labs(y=y_label) +
    scale_y_continuous(expand=c(0,0,0.0,0)) +
    scale_x_discrete(expand=c(0,0,0,0)) +
    ggplot2::theme(axis.title.x=element_blank())

  if(connected == TRUE){
    plot <- plot + ggalluvial::geom_alluvium(alpha=0.6, color=alpha("grey",0.5), width=bar_width)
  }

  if(label == "last"){
    labelData <- plotData[plotData$group == last(levels(plotData$group)),]
    plot <- plot + ggrepel::geom_text_repel(data=labelData, stat = "stratum", aes(label=population))
  }

  return(plot)
}
