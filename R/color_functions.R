##color_functions.R


clusters <- FetchData(object, vars=c("celltype", "ident"))

#' Expand color vector for subclusters
#'
#' Make a vector of colors for a vector of subclusters grouped by their
#' assignment to a cluster
#'
#' @param subcluster  vector of subcluster assignments (can be duplicated across cell assignments)
#' @param cluster vector of cluster assignment for each subcluster (vector of same length as subcluster)
#' @param cluster_colors (named) vector of colors for each subcluster

color_subcluster <- function(subcluster, cluster, cluster_colors=c()){
  clusters <- data.frame(subcluster=as.character(subcluster), cluster=as.character(cluster)) %>%
    group_by(cluster, subcluster) %>%
    # get number of cells in each subcluster
    summarize(count=n()) %>%
    arrange(cluster, desc(count)) %>%
    ungroup() %>%
    mutate(clustertype=make.unique(subcluster)) %>%
    group_by(cluster) %>%
    mutate(rank=rank(count)) %>%
    mutate(color=colorspace::lighten(cluster_colors[cluster],(rank-(max(rank)/2))*0.15))

  colors.clustertype <- clusters[['color']]
  names(colors.clustertype) <- clusters[['subcluster']]

  return(colors.clustertype)
}

## Make subclusters within each celltype same color with different brightness

