

#' Render Stacked Plot
#'
#' This function renders a stacked plot based on data and change points
#'
#' @param trueData Data.frame of data, the first column is date and the second
#'     column is the value
#' @param parCPData Data.frame of change points. First column is start point,
#'     second column is end point, third is beta value, fourth is variance1, and
#'     fifth is variance2.
#' @param title String to be title of the plot
#' @param subtitle String to be subtitle of the plot
#' @param dataName String to give name for coloring
#' @param varPlots (Optional) Boolean to indicate if variance plots should also
#'     be built. Default is TRUE.
#'
#' @return GGplot indicating data wiht change points, beta values, and
#'     potentially variance values
#' @export
#'
#' @examples
#'  data <- UKcovid[UKcovid$nation=='Wales',c(2,3)]
#'  data <- data[order(data$date),]
#'
#'  result_Vost <- detectAllChangePoints(fullData=data,
#'                         method='Vostrikova',
#'                         lower=c(-Inf, 0, 10^-8, -Inf),
#'                         upper=c(Inf, Inf, Inf, Inf),
#'                         alpha=0.05, nStart=NA, nEnd=NA)
#'  data_plot <- stackedPlot(trueData = data,
#'                         parCPData = result_Vost[,c(1:2,7:10)],
#'                         title = NULL, subtitle = NULL,
#'                         dataName = nation, varPlots = FALSE)
renderStackedPlot <- function(trueData, parCPData, title, subtitle,
                        dataName, varPlots = TRUE){
  ## Setup Data
  colnames(trueData) <-  c('Date','Value')
  colnames(parCPData) <- c('st','en','b1','v1','v2')


  colors <- RColorBrewer::brewer.pal(5, 'Set1')

  if(varPlots){
    dataPlot <-
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x=Date, y=Value, col='Data',
                                      linetype='Data'), trueData,
                key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = 'blank') +
      ggplot2::geom_blank(ggplot2::aes(col='beta', linetype='beta'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var1', linetype='var1'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var2', linetype='var2'),
                 key_glyph = "blank") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank()) +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",'data', 'var1', 'var2'),
                          values = c(colors[c(1,5)],'black',colors[3:4])) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', "data", 'var1', 'var2'),
                            values = c('solid',"dashed",'solid','solid','solid'))

    bPlot <-
      ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = st, y = b1, xend = en, yend = b1,
                       colour = "beta",  linetype='beta'), data = parCPData,
                   key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='Data', linetype='Data'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var1', linetype='var1'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var2', linetype='var2'),
                 key_glyph = "blank") +
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank()) +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",'data', 'var1', 'var2'),
                          values = c(colors[c(1,5)],'black',colors[3:4])) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', "data", 'var1', 'var2'),
                            values = c('solid',"dashed", 'solid', 'solid', 'solid'))

    v1Plot <-
      ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = st, y = v1, xend = en, yend = v1,
                       colour = "var1",  linetype='var1'), data = parCPData,
                   key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='Data', linetype='Data'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='beta', linetype='beta'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var2', linetype='var2'),
                 key_glyph = "blank") +
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank()) +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",'data', 'var1', 'var2'),
                          values = c(colors[c(1,5)],'black',colors[3:4])) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', "data", 'var1', 'var2'),
                            values = c('solid',"dashed", 'solid', 'solid', 'solid'))

    v2Plot <-
      ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = st, y = v2, xend = en, yend = v2,
                       colour = "var2",  linetype='var2'), data = parCPData,
                   key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='Data', linetype='Data'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='beta', linetype='beta'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col='var1', linetype='var1'),
                 key_glyph = "blank") +
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank()) +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",'data', 'var1', 'var2'),
                          values = c(colors[c(1,5)],'black',colors[3:4])) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', "data", 'var1', 'var2'),
                            values = c('solid',"dashed", 'solid', 'solid', 'solid'))
    # Combine
    plot <-
      dataPlot + bPlot + v1Plot + v2Plot +
      patchwork::plot_layout(nrow = 5, guides = "collect",
                  heights = c(3,1,1,1)) + #unit(c(1,2),c('cm','null'))) +
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                      plot.subtitle = ggplot2::element_text(hjust = 0.5))
      ) &
      ggplot2::theme(legend.position='bottom')
  } else{
    dataPlot <-
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x=Date, y=Value, col=dataName,
                                      linetype=dataName), trueData,
                key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = 'blank') +
      ggplot2::geom_blank(ggplot2::aes(col='beta', linetype='beta'),
                 key_glyph = "blank") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(size=14)) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",dataName),#c('beta',dataName,"CPs"),
                          values = c(colors[c(1,5)],'black')) +#c('black','black','red')) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', dataName),#c('beta',dataName,"CPs"),
                            values = c('solid', 'dashed',"solid"))#xxx

    bPlot <-
      ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = st, y = b1, xend = en, yend = b1,
                       colour = "beta",  linetype='beta'), data = parCPData,
                   key_glyph = "path") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=parCPData$st[-1], col='CPs',
                                       linetype='CPs'),
                 key_glyph = "blank") +
      ggplot2::geom_blank(ggplot2::aes(col=dataName, linetype=dataName),
                 key_glyph = "blank") +
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(size=14)) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::scale_colour_manual(name = "Method:",
                          labels = c('beta',"CPs",dataName),#c('beta',dataName,"CPs"),
                          values = c(colors[c(1,5)],'black')) +#c('black','black','red')) +
      ggplot2::scale_linetype_manual(name = "Method:",
                            labels = c('beta','CPs', dataName),#c('beta',dataName,"CPs"),
                            values = c('solid', 'dashed',"solid"))#xxx

    # Combine
    plot <-
      dataPlot + bPlot +
      patchwork::plot_layout(nrow = 2, guides = "collect",
                  heights = c(2,1)) + #unit(c(1,2),c('cm','null'))) +
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=22),
                      plot.subtitle = ggplot2::element_text(hjust = 0.5))
      ) &
      #theme(legend.position='bottom')
      ggplot2::theme(legend.position='none')
  }

  plot
}
