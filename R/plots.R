lsq_plot <- function(object ) {

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
  p <- ggplot( object$cumulativefrequency, aes( x=inv_f, y=M_f, col = "1") ) +
    geom_smooth(method = "lm", formula = y ~ x + 0, se=FALSE)   + 
    geom_point(aes(colour="2")) +
    scale_colour_manual(values = c("firebrick","black"),
                        labels = c("Best fit line", "Data"),
                        name = "") +
    xlab( "Inverse allelic frequency 1/f" ) +
    ylab( "Cumulative number \nof mutations M(f)" ) +
    ggtitle("Linear model best fit") +
    scale_x_continuous( trans=identity_trans(), breaks=breakPos,
                            labels=breakLab  ) +
    stat_poly_eq(aes(label =  paste(..eq.label..)), 
                 formula = formula, parse = TRUE,
                 label.y.npc = 0.8, col = "Black") +
    stat_poly_eq(formula = formula, 
                 parse = TRUE,
                 label.y.npc = 0.9, col = "Black") +
    theme(legend.position = c(0.8, 0.15))


  return(p)
}

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

  p <- ggplot(df, aes(x = inv_f, y = M_f, col = data)) +
    geom_line(size = 2, alpha = 0.5) +
    xlab("Time")+
    ylab("Population size") +
    scale_colour_manual(values=c("firebrick","black"), name = "") +
    scale_fill_manual(values=c("firebrick","black"), name = "") +
    xlab( "Inverse allelic frequency 1/f" ) +
    ylab( "Normalized M(f)" ) +
    ggtitle("Normalized cumulative distribution" ) +
    scale_x_continuous(trans=identity_trans(), 
                       breaks=breakPos,
                       labels=breakLab) +
    theme(legend.position = c(0.8, 0.15))

  return(p)

}

vaf_histogram <- function(object){
  p <- ggplot( data.frame(x=object$VAF), aes(x=x) ) +
          geom_histogram(binwidth=0.01) +
          xlab( "Allelic frequency f" ) +
          ylab("Number of \nmutations") +
          xlim( -0.01, 1.01) +
          ggtitle("VAF histogram")

  return(p)
}


plot_all <- function(object){
  
  p1 <- vaf_histogram(object)
  p2 <- lsq_plot(object)
  p3 <- normalized_plot(object)
  
  p <- plot_grid(p1, p2, p3,
            labels = c("A", "B", "C"), ncol = 3)
  
  return(p)
  
}

