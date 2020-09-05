#' Change default ggplot theme
ggplot_set_theme <- function(){
  theme_set(theme_bw() +
              theme(
                panel.grid.minor=element_blank(),
                strip.background=element_blank(),
                plot.background=element_blank(),
                legend.position = "right",
                axis.text.x = element_text(angle=45,hjust=1),
                plot.margin = unit(c(1,1,1,1),"mm")))
}

#' Get a vector of more or less distinct colors
#'
#' Could be improved my manual curation.
colors_get_distinct <- function(){
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  return(col_vector)
}
