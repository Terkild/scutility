#' Expand color vector for subclusters
#'
#' Make a vector of colors for a vector of subclusters grouped by their
#' assignment to a cluster
#'
#' @param subcluster  vector of subcluster assignments (can be duplicated across cell assignments)
#' @param cluster vector of cluster assignment for each subcluster (vector of same length as subcluster)
#' @param cluster_colors (named) vector of colors for each cluster
#' @param brightness_change how big a change in color brightness should be applied at each step (between 0 and 1)
#'
#' @return vector of colors
#' @export

color_subcluster <- function(subcluster, cluster, cluster_colors=c(), brightness_change=0.15){
  clusters <- data.frame(subcluster=as.character(subcluster), cluster=as.character(cluster)) %>%
    group_by(cluster, subcluster) %>%
    # get number of cells in each subcluster
    summarize(count=n()) %>%
    arrange(cluster, desc(count)) %>%
    ungroup() %>%
    mutate(clustertype=make.unique(subcluster)) %>%
    group_by(cluster) %>%
    mutate(rank=rank(count)) %>%
    mutate(color=colorspace::lighten(cluster_colors[cluster],(rank-(max(rank)/2))*brightness_change))

  colors.clustertype <- clusters[['color']]
  names(colors.clustertype) <- clusters[['subcluster']]

  return(colors.clustertype)
}

#' Color parent clusters by children
#'
#' Merges subcluster colors to a mixed color for parent
#'
#' @param clusters  vector of cluster annotations with length equal to subclusters
#' @param subclusters vector of subclusters
#' @param subcluster_colors vector of colors for each subcluster
#'
#' @import DescTools
#' @import magrittr
#' @export


color_parentcluster <- function(clusters, subclusters, subcluster_colors){
  newcolors <- data.frame(colors=subcluster_colors, subclusters=subclusters, clusters=clusters) %>%
    group_by(clusters) %>%
    summarize(newcolor=Reduce("MixColor", colors))

  newcolor <- newcolors$newcolor
  names(newcolor) <- newcolors$cluster

  return(newcolor)
}
