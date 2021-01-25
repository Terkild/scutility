#' Expand color vector for subclusters
#'
#' Make a vector of colors for a vector of subclusters grouped by their
#' assignment to a cluster
#'
#' @param subcluster  vector of subcluster assignments (can be duplicated across cell assignments)
#' @param cluster vector of cluster assignment for each subcluster (vector of same length as subcluster)
#' @param cluster_colors (named) vector of colors for each cluster
#' @param order should colors be ordered by subcluster size? Set to FALSE to use factor levels
#' @param brightness_change how big a change in color brightness should be applied at each step (between 0 and 1)
#'
#' @return vector of colors
#' @export

color_subcluster <- function(subcluster, cluster, cluster_colors=c(), order=TRUE, brightness_change=NA, max_change=1){

  clusters <- data.frame(subcluster=subcluster, cluster=cluster) %>%
    group_by(cluster, subcluster) %>%
    summarize(count=n())

  if(order == TRUE){
    clusters <- clusters %>%
      arrange(cluster, desc(count)) %>%
      group_by(cluster) %>%
      mutate(rank=rank(count))
  } else {
    clusters <- clusters %>%
      group_by(cluster) %>%
      mutate(rank=rank(subcluster))
  }


  if(!is.na(brightness_change)){
    ## Deprecated, but kept for consistency
    clusters <- clusters %>% mutate(color=colorspace::lighten(cluster_colors[cluster], amount=-(median(rank)-rank)*brightness_change))
  } else {
    clusters <- clusters %>%
      mutate(brightness_change=max_change/max(rank)) %>%
      mutate(color=colorspace::lighten(cluster_colors[cluster], amount=(median(rank)-rank)*brightness_change))
  }

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
#' @param subcluster_colors Named vector of colors for each subcluster or vector of colors of equal length as subclusters
#' @param mix Should colors of the two largest subclusters be mixed (by DescTools::MixColor)? If FALSE, the color from the largest subcluster will be used
#'
#' @import DescTools
#' @import magrittr
#' @export


color_parentcluster <- function(clusters, subclusters, subcluster_colors, mix=TRUE){

  if(length(subcluster_colors) != length(subclusters)){
    subcluster_colors <- subcluster_colors[subclusters]
  }

  newcolors <- data.frame(colors=subcluster_colors, subclusters=subclusters, clusters=clusters) %>%
    group_by(subclusters) %>% mutate(subcluster_count=n()) %>%
    arrange(subcluster_count) %>%
    group_by(clusters)

  if(mix == TRUE){
    newcolors %<>%
      summarize(newcolor=Reduce("MixColor", unique(colors)))
  } else {
    newcolors %<>%
      summarize(newcolor=colors[which.max(subcluster_count)])
  }

  newcolor <- newcolors$newcolor
  names(newcolor) <- newcolors$cluster

  return(newcolor)
}
