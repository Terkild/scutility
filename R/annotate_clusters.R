#' Merge annotation by cluster
#'
#' Merges individual cell annotations by cluster using majority rule (the entire cluster is annotated as the individual annotation represented by most cells within the cluster)
#'
#'

annotate_merge_by_cluster <- function(annotation, cluster, threshold="max", collapse="/"){

  group_merge <- data.frame(ann=annotation, cluster=cluster) %>%
    filter((ann  %in% exclude) == FALSE) %>%
    group_by(cluster) %>% mutate(cluster_count=n()) %>%
    group_by(ann, cluster, cluster_count) %>%
    summarize(count=n()) %>% mutate(freq=count/cluster_count) %>%
    group_by(cluster)

  if(threshold == "max"){
    group_merge <- group_merge %>% summarize(celltype=ann[which.max(count)])
  } else {
    group_merge <- group_merge %>% summarize(celltype=paste(ann[which(freq >= threshold)], collapse=collapse))
  }


  newcluster <- data.frame(annotation=annotation) %>%
    left_join(group_merge) %>%
    select(celltype) %>% .[[1]]

  return(newcluster)
}

#' Number subclusters by parent cluster

subcluster_number <- function(subcluster, cluster){
  data_subcluster <- data.frame(subcluster=subcluster, cluster=cluster) %>%
    group_by(cluster, subcluster) %>% summarize(count=n()) %>%
    group_by(cluster) %>% mutate(rank = rank(-count, ties.method='first'), num_rank=n()) %>%
    mutate(subcluster_name = ifelse(num_rank>1, paste0(cluster,".",rank), cluster))

  newcluster <- data.frame(subcluster=subcluster) %>%
    left_join(data_subcluster) %>%
    select(subcluster_name) %>% .[[1]]

  return(newcluster)
}



