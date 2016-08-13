lsq_cumulative_plot <- function(object ) {

  # Set values for u and v based on the input data:
  u <- max( object$cumulativefrequency$f )
  v <- min( object$cumulativefrequency$f )

  # Return nothing if less than 5 sample points
  if( nrow( object$cumulativefrequency ) < 5 ) return(NULL)

  # Get good positions for breaks and set label variables:
  breaks    <- c(u,round((u-v) / 4,2),v)
  breakPos  <- 1 / breaks - 1 / u
  breakLab  <- paste("1/", breaks,sep="")
  labR2_mu  <- paste0( "atop( R^{2}==", sprintf( "%0.3f", object$rsq$metric ),
                       ", mu/beta == ", sprintf("%0.3f", object$mu ), ")" )
  xposR2_mu <- 1 / (u - 0.2)
  yposR2_mu <- 0.9 * max( object$cumulativefrequency$M_f )

  # Create the plot:
  p <- ggplot( object$cumulativefrequency, aes( x=inv_f ) ) +
        geom_abline( aes(slope=object$mu, intercept=0, color="1"), linetype=1,size=1 ) +
        geom_point(aes(y=M_f,colour="2")) +
        scale_colour_manual(values = c("firebrick","black"),
                        labels = c("Best fit line", "Data"),
                        name = "") +
        xlab( "Inverse allelic frequency 1/f" ) +
        ylab( "Cumulative number of mutations M(f)" ) +
        ggtitle("Linear model best fit") +
        scale_x_continuous( trans=identity_trans(), breaks=breakPos,
                            labels=breakLab  ) +
        annotate( "text", x=xposR2_mu, y=yposR2_mu, label=labR2_mu, parse=TRUE, hjust=1 ) +
        theme_bw()



  return(p)
}

norm_cumulative_plot <- function(object){

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

  p <- ggplot(df, aes(x = inv_f, y = M_f, col = data)) +
    geom_line(size = 2, alpha = 0.5) +
    theme_bw() +
    xlab("Time")+
    ylab("Population size") +
    scale_colour_manual(values=c("dodgerblue","firebrick")) +
    scale_fill_manual(values=c("dodgerblue","firebrick")) +
    xlab( "Inverse allelic frequency 1/f" ) +
    ylab( "Normalized M(f)" ) +
    ggtitle("Normalized cumulative distribution" ) +
    scale_x_continuous( trans=identity_trans(), breaks=breakPos,
                        labels=breakLab  ) +
    annotate( "text", x=xpos_metric, y=ypos_metric, label=lab_metric, hjust=1 ) +
    theme(legend.title=element_blank())

  return(p)

}

vaf_histogram <- function(object){
  p <- ggplot( data.frame(x=object$VAF), aes(x=x) ) +
          geom_histogram(binwidth=0.01) +
          xlab( "Allelic frequency f" ) +
          ylab("Number of mutations") +
          xlim( -0.01, 1.01) +
          ggtitle("Variant allele frequency histogram") +
          theme_bw()

  return(p)
  }
