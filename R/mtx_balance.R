#' Load MTX file as long format data.frame
#'
#' @export
mtx_load_long <- function(path, file_matrix="matrix.mtx", file_barcodes="barcodes.tsv", file_features="genes.tsv", matrix_sep=" ", matrix_skip=3, matrix_colnames=c("feature", "barcode", "count")){
  mtx <- read.table(file.path(path,file_matrix), skip=matrix_skip, header=FALSE, sep=matrix_sep)
  barcodes <- read.table(file.path(path,file_barcodes), header=FALSE)[[1]]
  features <- read.table(file.path(path,file_features), header=FALSE)[[1]]
  colnames(mtx) <- matrix_colnames

  mtx[,1] <- factor(features[mtx[,1]])
  mtx[,2] <- factor(barcodes[mtx[,2]])

  return(mtx)
}

#' Save long format mtx data.frame as mtx
#'
#' @export
mtx_long_save <- function(mtx, path, file_matrix="matrix.mtx", file_barcodes="barcodes.tsv", file_features="genes.tsv", matrix_header=c("%%MatrixMarket matrix coordinate integer general","%"), matrix_sep=" "){
  dir.create(path=path,showWarnings=FALSE, recursive=TRUE)

  write.table(matrix_header, file.path(path, file_matrix), col.names=FALSE, row.names=FALSE, quote=FALSE)

  mtx[,1] <- droplevels(mtx[,1])
  mtx[,2] <- droplevels(mtx[,2])

  write.table(data.frame(levels(mtx[,1]),levels(mtx[,1])), file.path(path, file_features), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  write.table(levels(mtx[,2]), file.path(path, file_barcodes), col.names=FALSE, row.names=FALSE, quote=FALSE)


  mtx[,1] <- as.integer(mtx[,1])
  mtx[,2] <- as.integer(mtx[,2])

  write.table(paste(c(max(mtx[,1]), max(mtx[,2]), nrow(mtx)), collapse=matrix_sep), file.path(path, file_matrix), col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
  write.table(mtx, file.path(path, file_matrix), col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep=matrix_sep)
}

#' Balance/Upsample mtx
#'
#' @import dplyr
#' @importFrom purrr map2_dfr
#' @export
mtx_balance <- function(path,
                        n_samples,
                        path_save=NULL,
                        file_matrix="matrix.mtx",
                        file_barcodes="barcodes.tsv",
                        file_features="genes.tsv",
                        matrix_sep=" ",
                        matrix_colnames=c("feature", "barcode", "count"),
                        expected_counts=c(),
                        ...){

  # Load matrix from files
  mtx <- mtx_load_long(path, file_matrix=file_matrix, file_barcodes=file_barcodes, file_features=file_features, matrix_sep=matrix_sep, matrix_colnames=matrix_colnames, ...)

  ## Determine how many UMIs to sample from each tag
  feature_sum <- mtx %>%
    group_by(feature) %>%
    summarize(sum=sum(count)) %>%
    # Only include the top N tags where N is the number of samples given
    top_n(n=n_samples, wt=sum) %>%
    ungroup() %>%
    mutate(sample_size=max(sum)-sum) %>%
    as.data.frame()

  # if expected cell counts from each sample is not equal correct for this
  if(length(expected_counts) == nrow(feature_sum)){
    if(length(intersect(names(expected_counts), feature_sum$feature)) != nrow(feature_sum)) stop("expected_counts needs to have same names as tags")

    feature_sum <- feature_sum %>%
      mutate(sum_per_cell=sum/expected_counts[as.character(feature)]) %>%
      mutate(sample_size_per_cell=max(sum_per_cell)-sum_per_cell,
             sample_size=sample_size_per_cell*expected_counts[as.character(feature)]) %>%
      as.data.frame()
  }

  # Upsampling
  mtx_up <- mtx %>%
    filter(feature %in% feature_sum$feature) %>%
    group_split(feature) %>%
    map2_dfr(feature_sum$sample_size, ~ if(.y > 0){slice_sample(.x, n = .y, weight_by=.x$count, replace=TRUE)}) %>%
    select(feature, barcode) %>%
    group_by(feature, barcode) %>%
    summarize(count_add=n())

  # Make new long format matrix
  mtx_new <- mtx %>%
    filter(feature %in% feature_sum$feature) %>%
    left_join(mtx_up) %>%
    mutate(sum=count+ifelse(is.na(count_add),0,count_add))

  # Save balanced matrix to files
  if(!is.null(path_save)){
    mtx_long_save(mtx_new[,c(matrix_colnames[1:2],"sum")], path=path_save, file_matrix=file_matrix, file_barcodes=file_barcodes, file_features=file_features, matrix_sep=matrix_sep)
    return(mtx_new)#return(TRUE)
  } else {

    # or return balanced matrix in long format data.frame
    return(mtx_new)
  }
}


#' Balance/Upsample mtx
#'
#' @import dplyr
#' @importFrom purrr map2_dfr
#' @export
mtx_balance_matrix <- function(matrix,
                                n_samples=NULL,
                                expected_counts=c(),
                                ...){

  # Load matrix
  mtx <- data.frame(feature=rownames(matrix), barcode=rep(colnames(matrix), each=nrow(matrix)), count=matrix[seq(1:length(matrix))])

  ## Determine how many UMIs to sample from each tag
  feature_sum <- mtx %>%
    group_by(feature) %>%
    summarize(sum=sum(count))


  if(length(setdiff(names(expected_counts), feature_sum$feature)) < 1){

    # only include tags expected to have cell counts
    feature_sum <- feature_sum %>% filter(feature %in% names(expected_counts))

  } else if(is.null(n_sample)){

    stop("n_sample or expected_counts has to be set")

  } else {
    # Only include the top N tags where N is the number of samples given
    # top_n(n=n_samples, wt=sum) %>%
    feature_sum <- feature_sum %>% top_n(n=n_samples, wt=sum)
  }

  feature_sum <- feature_sum %>% ungroup() %>%
    mutate(sample_size=max(sum)-sum) %>%
    as.data.frame()

  # if expected cell counts from each sample is not equal correct for this
  if(length(expected_counts) == nrow(feature_sum)){
    if(length(intersect(names(expected_counts), feature_sum$feature)) != nrow(feature_sum)) stop("expected_counts needs to have same names as tags")

    feature_sum <- feature_sum %>%
      mutate(sum_per_cell=sum/expected_counts[as.character(feature)]) %>%
      mutate(sample_size_per_cell=max(sum_per_cell)-sum_per_cell,
             sample_size=sample_size_per_cell*expected_counts[as.character(feature)]) %>%
      as.data.frame()
  }

  #print(feature_sum)

  # Upsampling
  mtx_up <- mtx %>%
    filter(feature %in% feature_sum$feature) %>%
    group_split(feature) %>%
    map2_dfr(feature_sum$sample_size, ~ if(.y > 0){slice_sample(.x, n = .y, weight_by=.x$count, replace=TRUE)}) %>%
    select(feature, barcode) %>%
    group_by(feature, barcode) %>%
    summarize(count_add=n())

  #print(mtx_up %>% group_by(feature) %>% summarize(count=sum(count_add)))

  # Make new long format matrix
  mtx_new <- mtx %>%
    filter(feature %in% feature_sum$feature) %>%
    left_join(mtx_up) %>%
     mutate(count_add=ifelse(is.na(count_add),0,count_add)) %>%
     mutate(sum=count+count_add)

  #print(mtx_new %>% group_by(feature) %>% summarize(count=sum(count),count_add=sum(count_add),sum=sum(sum)))

  return(mtx_new)
}
