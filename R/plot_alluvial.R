#' Plot alluvial from seurat data
#'
#' @param object	Seurat object
#' @param group1 	columns to be included as group1 (x-axis bar 1)
#' @param group2 	columns to be included as group2 (x-axis bar 2)
#'
#' @importFrom Seurat FetchData
#' @export

seurat_plot_alluvial <- function(object, group1, group2){
  getData <- Seurat::FetchData(object, vars=c(group1,group2))

  plot_alluvial(group1=getData[,1], group2=getData[,2])
}

#' Plot alluvial from seurat data
#'
#' @param Seurat object
#' @param group1 	vector to be included as group1 (x-axis bar 1)
#' @param group2 	vector to be included as group2 (x-axis bar 2)
#'
#' @import ggalluvial
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @export
plot_alluvial <- function(group1, group2){

  plotData <- data.frame(group1=group1, group2=group2) %>%
    group_by(group1, group2) %>%
    summarize(count=n())

  ggplot(plotData, aes(y=count, axis1=group1, axis2=group2, fill=group1)) +
    ggalluvial::geom_alluvium(width = 1/8) +
    ggalluvial::geom_stratum(width = 1/8) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(values=scutility::colors_get_distinct())
}
