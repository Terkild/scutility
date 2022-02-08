#' Plot quadrants
#'
#' @param plotData data frame containing data to be plotted
#' @param x column name for x data
#' @param y column name for y data
#' @param color column name for colouring
#' @param wrap  column name for wrapping plot
#' @param trans transformation function for x and y scales (log1p default)
#' @param breaks breaks for x and y scales
#' @param threshold_x x gating threshold
#' @param threshold_y y gating threshold
#' @param text_pos_lower text position
#' @param jitter jitter function for geom_points
#'
#' @return ggplot
#' @export
#'
plot_quadrants <- function(plotData, x, y, color, wrap=NULL,
                           trans="log1p", breaks=c(0,1,2,3,5,10,25,50,100,200,500,1000),
                           threshold_x=0.5, threshold_y=0.5,
                           text_pos_lower=log1p(-0.25), jitter=position_jitter(width=0.18, height=0.18, seed=124)){
  if(is.null(wrap)){
    plotData$wrap <- "all"
    wrap <- "wrap"
    add_wrap <- NULL
  } else {
    add_wrap <- facet_wrap(~ .data[[wrap]])
  }

  quadrant_stats <- plotData %>%
    mutate(xpos=as.integer(.data[[x]] > threshold_x),
           ypos=as.integer(.data[[y]] > threshold_y)) %>%
    group_by(.data[[wrap]], xpos, ypos) %>%
    summarize(quadrant_count=n()) %>%
    group_by(.data[[wrap]]) %>%
    mutate(group_sum=sum(quadrant_count),
           quadrant_pct=quadrant_count/group_sum)

  plotData %>% .[sample(nrow(plotData),nrow(plotData)),] %>%
    ggplot(aes(x=.data[[x]], y=.data[[y]], color=.data[[color]])) +
    geom_point(position=jitter, size=0.5, alpha=0.5) +
    geom_text(data=quadrant_stats, aes(x=ifelse(xpos==1, Inf, text_pos_lower), y=ifelse(ypos==1, Inf, text_pos_lower), hjust=xpos*1.5, vjust=ypos*1.5, label=sprintf("%.2f", quadrant_pct*100)), color="black") +
    scale_y_continuous(trans=trans, breaks=breaks[which(breaks <= max(plotData[[y]]))]) +
    scale_x_continuous(trans=trans, breaks=breaks[which(breaks <= max(plotData[[x]]))]) +
    #geom_vline(xintercept=eval(parse(text=paste0(trans,"(threshold_x)")))) + geom_hline(yintercept=eval(parse(text=paste0(trans,"(threshold_y)")))) +
    geom_vline(xintercept=threshold_x) + geom_hline(yintercept=threshold_y) +
    add_wrap
}
