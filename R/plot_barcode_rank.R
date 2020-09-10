#' Barcode-Rank plot
#'
#' @param data  Count matrix
#' @param lower Lower threshold passed on to DropletUtils::barcodeRanks()
#' @param show.inflection Plot line at inflection point?
#' @param show.knee Plot line at knee point?
#'
#' @return ggplot object
#' @import ggplot2
#' @importFrom DropletUtils barcodeRanks
#' @importFrom tibble tibble
#' @export

plot_barcode_rank <- function(data, lower=10, show.inflection=TRUE, show.knee=FALSE){
  require(ggplot2)

  bc_rank <- DropletUtils::barcodeRanks(data, lower=lower)

  knee_plt <- tibble::tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>%
    distinct() %>%
    dplyr::filter(total > 0)

  meta <- S4Vectors::metadata(bc_rank)
  annot <- tibble::tibble(inflection = meta[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > meta[["inflection"]]]),
                  knee = meta[["knee"]],
                  knee_cutoff = max(bc_rank$rank[bc_rank$total > meta[["knee"]]]))

  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")

  if(show.inflection == TRUE){
    p <- p + geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
      geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
      geom_label(aes(y=rank_cutoff,x=Inf,label=rank_cutoff), data = annot, hjust=1, vjust=0) +
      geom_label(aes(y=0,x=inflection,label=inflection), data = annot, hjust=1, vjust=0)
  }

  if(show.knee == TRUE){
    p <- p + geom_hline(aes(yintercept = knee_cutoff), data = annot, linetype = "dotted", col="red") +
      geom_vline(aes(xintercept = knee), data = annot, linetype = "dotted", col="red") +
      geom_label(aes(y=knee_cutoff,x=Inf,label=knee_cutoff), data = annot, hjust=1, vjust=1, col="red") +
      geom_label(aes(y=0,x=knee,label=knee), data = annot, hjust=0, vjust=0, col="red")
  }

  return(p)
}
