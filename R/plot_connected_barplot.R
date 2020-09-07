#' Plot connected barplot from seurat object
#'
#' @param object Seurat object
#' @param group.by Meta data column containing population/cluster assignment (y-axis)
#' @param split.by Meta data column containing group or sample assignment (x-axis)
#'
#' @returns ggplot object
#' @importFrom Seurat FetchData
#' @export

seurat_plot_connected_barplot <- function(object, group.by="ident", split.by, ...){
  getData <- Seurat::FetchData(object, vars=c(group.by, split.by))
  colnames(getData)[1:2] <- c("group.by","split.by")

  plot_connected_barplot(population=getData$group.by, group=getData$split.by, ...)
}

#' Plot connected barplot
#'
#' @param population  vector containing population/cluster assignment (y-axis)
#' @param group vector containing group or sample assignment (x-axis)
#' @param value what statistic to plot ("percent" or "count")
#' @param order should populations be ordered by total count (across groups)
#' @param color (named) vector of colors for populations
#'
#' @returns ggplot object
#' @import ggplot2
#' @export

plot_connected_barplot <- function(population, group, value="percent", order=FALSE, colors=c(), label=FALSE){
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


  if(value == "percent"){
    plotData <- plotData %>%
      mutate(pct=populationCount/groupCount*100) %>%
      mutate(value=pct)
  } else {
    plotData <- plotData %>%
      mutate(value=populationCount)
  }


  plot <- ggplot(plotData, aes(x=group, y=value, fill=population, alluvium=population, stratum=population)) +
    ggalluvial::geom_alluvium(alpha=0.6, color=alpha("grey",0.5), width=0.4) +
    ggalluvial::geom_stratum(color=alpha("black",0.5), width=0.4) +
    scale_fill_manual(values=colors) +
    #labs(y="% of PBMC", x="Days post onset") +
    guides(fill=F) +
    scale_y_continuous(expand=c(0,0,0.05,0))

  if(label == "last"){
    labelData <- plotData[plotData$group == last(levels(plotData$group)),]
    plot <- plot + ggrepel::geom_text_repel(data=labelData, stat = "stratum", aes(label=population))
  }

  return(plot)
}
