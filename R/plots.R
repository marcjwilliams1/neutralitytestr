#' Plot cumulative distribution
#' \code{lsq_plot} Plots the cumulative distribution of the data as well as the best fit linear model line.
#'
#' @param object neutrality test object
#' @return ggplot object.
#' @examples
#' lsq_plot(neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25))
#' @export
lsq_plot <- function(object) {

  # Set values for u and v based on the input data:
  u <- max( object$cumulativefrequency$f )
  v <- min( object$cumulativefrequency$f )

  # Return nothing if less than 5 sample points
  if( nrow( object$cumulativefrequency ) < 5 ) return(NULL)

  # Get good positions for breaks and set label variables:
  breaks    <- c(u,round((u-v) / 4,2),v)
  breakPos  <- 1 / breaks - 1 / u
  breakLab  <- paste("1/", breaks,sep="")

  formula <- object$cumulativefrequency$M_f ~ object$cumulativefrequency$inv_f + 0
  # Create the plot:
  p <- ggplot2::ggplot( object$cumulativefrequency, ggplot2::aes( x=inv_f, y=M_f, col = "1") ) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x + 0, se=FALSE)   +
    ggplot2::geom_point(ggplot2::aes(colour="2")) +
    ggplot2::scale_colour_manual(values = c("firebrick","black"),
                        labels = c("Best fit line", "Data"),
                        name = "") +
    ggplot2::xlab( "Inverse allelic frequency 1/f" ) +
    ggplot2::ylab( "Cumulative number \nof mutations M(f)" ) +
    ggplot2::ggtitle("Linear model best fit") +
    ggplot2::scale_x_continuous( trans=scales::identity_trans(), breaks=breakPos,
                            labels=breakLab  ) +
    ggpmisc::stat_poly_eq(ggplot2::aes(label =  paste(..eq.label..)),
                 formula = formula, parse = TRUE, eq.with.lhs = "italic(mu/beta)~`=`~",
                 label.y.npc = 0.8, col = "Black") +
    ggpmisc::stat_poly_eq(formula = formula,
                 parse = TRUE,
                 label.y.npc = 0.9, col = "Black") +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy", minor = "none") +
    ggplot2::theme(legend.position = c(0.5, 0.15))


  return(p)
}

#' Plot normalized cumulative distribution
#' \code{normalized_plot} Plots the (normalized) cumulative distribution of the data as well as the theoretical expectation from
#' a neutral evolutionary model.
#'
#' @param object neutrality test object
#' @return ggplot object.
#' @examples
#' normalized_plot(neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25))
#' @export
normalized_plot <- function(object){

  # Set values for u and v based on the input data:
  u <- max( object$cumulativefrequency$f )
  v <- min( object$cumulativefrequency$f )

  # Return nothing if less than 5 sample points
  if( nrow( object$cumulativefrequency ) < 5 ) return(NULL)

  # Format data into dataframe and group via empirical vs Theroetical for ggplot
  df <- data.frame(M_f = c(object$cumulativefrequency$nM_f, object$cumulativefrequency$tM_f))
  df$data <- c(rep("Empirical",length(object$cumulativefrequency$nM_f)), rep("Theoretical",length(object$cumulativefrequency$nM_f)))
  df$inv_f <- rep(object$cumulativefrequency$inv_f,2)

  # Get good positions for breaks and set label variables:
  breaks    <- c(u,round((u-v) / 4,2),v)
  breakPos  <- 1 / breaks - 1 / u
  breakLab  <- paste("1/", breaks,sep="")
  lab_metric  <- paste( "Area = ", sprintf( "%0.3f", object$area$metric ),
                       "\nDk = ", sprintf("%0.3f", object$Dk$metric ) )
  xpos_metric<- 1 / (u - 0.1)
  ypos_metric <- 0.9 * max( df$M_f )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = inv_f, y = M_f, col = data)) +
    ggplot2::geom_line(size = 2, alpha = 0.5) +
    ggplot2::xlab("Time")+
    ggplot2::ylab("Population size") +
    ggplot2::scale_colour_manual(values=c("firebrick","black"), name = "") +
    ggplot2::scale_fill_manual(values=c("firebrick","black"), name = "") +
    ggplot2::xlab( "Inverse allelic frequency 1/f" ) +
    ggplot2::ylab( "Normalized M(f)" ) +
    ggplot2::ggtitle("Normalized cumulative distribution" ) +
    ggplot2::scale_x_continuous(trans=scales::identity_trans(),
                       breaks=breakPos,
                       labels=breakLab) +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "xy", minor = "none") +
    ggplot2::theme(legend.position = c(0.5, 0.15))

  return(p)

}

#' Plot VAF histogram
#' \code{vaf_histogram} Plots a histogram of the variant allele frequencies.
#'
#' @param object neutrality test object
#' @return ggplot object.
#' @examples
#' vaf_histogram(neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25))
#' @export
vaf_histogram <- function(object){
  p <- ggplot2::ggplot( data.frame(x=object$VAF), ggplot2::aes(x=x) ) +
    ggplot2::geom_histogram(binwidth=0.01) +
    ggplot2::xlab( "Allelic frequency f" ) +
    ggplot2::ylab("Number of \nmutations") +
    ggplot2::xlim( -0.01, min(round(max(object$VAF), 1) + 0.1, 1.0)) +
    ggplot2::ggtitle("VAF histogram") +
    cowplot::theme_cowplot() + cowplot::background_grid(major = "xy", minor = "none")

  return(p)
}

#' Plot all plots in the package and make composite figure.
#' \code{plot_all} Plots histogram, linear model best fit plot
#' and normalized plot and plot and makes composite figure.
#'
#' @param object neutrality test object
#' @return ggplot object.
#' @examples
#' plot_all(neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25))
#' @export
#' @export
plot_all <- function(object){

  p1 <- vaf_histogram(object)
  p2 <- lsq_plot(object)
  p3 <- normalized_plot(object)

  p <- cowplot::plot_grid(p1, p2, p3,
            labels = c("A", "B", "C"), ncol = 3)

  return(p)

}

