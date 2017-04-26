#
# Collection of useful functions
#

oddsratioWald.proc <- function(n00, n01, n10, n11, alpha = 0.05){
  #
  #  Compute the odds ratio between two binary variables, x and y,
  #  as defined by the four numbers nij:
  #
  #    n00 = number of cases where x = 0 and y = 0
  #    n01 = number of cases where x = 0 and y = 1
  #    n10 = number of cases where x = 1 and y = 0
  #    n11 = number of cases where x = 1 and y = 1
  #
  OR <- (n00 * n11)/(n01 * n10)
  #
  #  Compute the Wald confidence intervals:
  #
  siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
  zalph <- qnorm(1 - alpha/2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  #
  ORlo <- exp(loglo)
  ORhi <- exp(loghi)
  #
  oframe <- data.frame(LowerCI = ORlo, OR = OR, UpperCI = ORhi, alpha = alpha)
  oframe
}

# Given a set of p-values, returns them into segments
discretePvals <- function( pVals_n_iP,
                           lims_iL = c( 0.001, 0.0001, 0.00001 ) ){
  pVals_n_iP[ is.nan( pVals_n_iP ) ] <- 1;
  pVals_n_iP[ is.na( pVals_n_iP ) ] <- 1;
  limsS_iL <- c( paste( ">", lims_iL[1], sep = "" ), 
                 paste( "<", lims_iL, sep = "" ) )
  pVals_s_iP <- sapply( pVals_n_iP,
                        function( x ){
                          r <- p( ">", lims_iL[1] );
                          for( iL in 1:length( lims_iL ) ){
                            if( x < lims_iL[iL] ){
                              r <- p( "<", lims_iL[iL] );
                            }
                          }
                          return( r );
                        } );
  pVals_s_iP <- factor( pVals_s_iP,
                        levels = limsS_iL );
  return( pVals_s_iP );
}

# Transform a sensible matrix into the dataframe oriented crap that ggplot needs
unfold_v2 <- function( data_iSiF,
                       iFs_unfold = 1:dim( data_iSiF )[2],
                       iFs_keep = NULL ){
  
  # Add dim names if they don't exist
  if( is.null( dimnames(data_iSiF)[[1]] ) ){
    dimnames(data_iSiF)[[1]] <- as.character( 1:dim( data_iSiF )[1] )
  }
  if( is.null( dimnames(data_iSiF)[[2]] ) ){
    dimnames(data_iSiF)[[2]] <- as.character( 1:dim( data_iSiF )[2] )
  }
  
  # Five dim names if they don't exist
  if( length( dimnames( data_iSiF )[[1]] ) == 0 &&
    is.null( dimnames( data_iSiF )[[1]] ) ){
    dimnames( data_iSiF )[[1]] <- as.character( 1:dim( data_iSiF )[1] );
  }
  if( length( dimnames( data_iSiF )[[2]] ) == 0 &&
    is.null( dimnames( data_iSiF )[[2]] ) ){
    dimnames( data_iSiF )[[2]] <- as.character( 1:dim( data_iSiF )[2] );
  }
  
  # Transform into indexes if names are provided
  if( is.character( iFs_unfold ) ){
    iFs_unfold <- match( iFs_unfold, dimnames( data_iSiF )[[2]] );
  }
  if( is.character( iFs_keep ) ){
    iFs_keep <- match( iFs_keep, dimnames( data_iSiF )[[2]] );
  }
  
  # Unfold matrix
  data_iSiP <- data.frame( );
  for( iF in 1:length( iFs_unfold ) ){
    data_iS <- data_iSiF[ ,iFs_unfold[iF]];
    data_iSiP <- rbind( data_iSiP , 
                        cbind( data_iS, 
                               rep( dimnames( data_iSiF )[[2]][iFs_unfold[iF]], 
                                    length( data_iS ) ),
                               data_iSiF[ , iFs_keep ],
                               dimnames( data_iSiF )[[1]] ) );
  }
  dimnames( data_iSiP )[[2]] <- c( "unfoldedValue", "oldColNames", dimnames( data_iSiF )[[2]][iFs_keep], "oldRowNames" );
  
  # Short dimentions
  data_iSiP[ ,"unfoldedValue"] <- as.numeric( as.character( data_iSiP[ ,"unfoldedValue"] ) );
  data_iSiP[ ,"oldRowNames"] <- factor( data_iSiP[ ,"oldRowNames"],
                                        levels = dimnames( data_iSiF )[[1]] );
  data_iSiP[ ,"oldColNames"] <- factor( data_iSiP[ ,"oldColNames"],
                                        levels = dimnames( data_iSiF )[[2]] );
  
  #Return stuff
  return( data_iSiP );
}


# Transform a sensible matrix into the dataframe oriented crap that ggplot needs
unfold <- function( data_iSiF,
                    iFs_unfold = 1:dim( data_iSiF )[2],
                    iFs_keep = NULL ){
  # Transform into indexes if names are provided
  if( is.character( iFs_unfold ) ){
    iFs_unfold <- match( iFs_unfold, dimnames( data_iSiF )[[2]] );
  }
  if( is.character( iFs_keep ) ){
    iFs_keep <- match( iFs_keep, dimnames( data_iSiF )[[2]] );
  }
  # Unfold matrix
  data_iSiP <- data.frame( );
  for( iF in 1:length( iFs_unfold ) ){
    data_iS <- data_iSiF[ ,iFs_unfold[iF]];
    data_iSiP <- rbind( data_iSiP , 
                        cbind( data_iS, 
                               rep( dimnames( data_iSiF )[[2]][iFs_unfold[iF]], 
                                    length( data_iS ) ),
                               data_iSiF[ , iFs_keep ]) );
  }
  dimnames( data_iSiP )[[2]] <- c( "unfoldedValue", "unfoldedName", dimnames( data_iSiF )[[2]][iFs_keep] );
  return( data_iSiP );
}

# Sensible paste
p = function(...){
  paste( ..., sep = "", collapse = "" );
}

# Converts a data frame into a 2D array
frame2array = function(dataFrame_iRiC) {
  array_iRiC = array(dim=dim(dataFrame_iRiC));
  for( iC in 1:dim(dataFrame_iRiC)[2] ){
    array_iRiC[,iC] <- dataFrame_iRiC[,iC];
  }
  dimnames(array_iRiC) <- dimnames(dataFrame_iRiC);
  return(array_iRiC);
}

# # Matlab repmat equivalent
# repmat = function(X,m,n){
#   mx = dim(X)[1]
#   nx = dim(X)[2]
#   matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
# }

# Matlab imagesc equivalent
# lims = Limits in x for the color scale
# main = title
# colorbar = TRUE if you want the colorbar to be shown
library(gplots)
image2F = function(x_iCiR,lims=NA,main="",Rowv=FALSE,Colv=FALSE,dendrogram="none",density.info=c("histogram","density","none")){
  #Check lims
  if( any(is.na(lims)) ){
    lims <- range(x_iCiR,na.rm=TRUE);
  }
  if( lims[1]>=lims[2] ){
    stop("lims[1] must be lower than lims[2]");
  }
  if( lims[2]<min(x_iCiR,na.rm=TRUE) |  lims[1]>max(x_iCiR,na.rm=TRUE) ){
    stop("range(x_iCiR) is out of the lims[c(1,2)]");
  }
  #Organise labels and some other fields
  iCs <- seq(1,dim(x_iCiR)[1],max(1,round(dim(x_iCiR)[1]/20)));
  iRs <- seq(1,dim(x_iCiR)[2],max(1,round(dim(x_iCiR)[2]/20)));
  labelsY <- NA*(1:dim(x_iCiR)[1]);
  if( is.null(rownames(x_iCiR)) ){
    labelsY[iCs] <- iCs;
  }else{
    labelsY[iCs] <- rownames(x_iCiR)[iCs];}
  labelsX <- NA*(1:dim(x_iCiR)[2]);
  if( is.null(colnames(x_iCiR)) ){
    labelsX[iRs] <- iRs;
  }else{
    labelsX[iRs] <- colnames(x_iCiR)[iRs];}
  names <- names(dimnames(x_iCiR));
  if( is.null(names) ){
    names <- c("colums","rows");
  }
  #Check NAs
  anynans <- function(x) any(is.na(x));
  if( dendrogram=="both" || dendrogram=="row" ){
    x_iCiR <- x_iCiR[,!apply(x_iCiR,2,anynans)];
  }
  if( dendrogram=="both" || dendrogram=="column" ){
    x_iCiR <- x_iCiR[!apply(x_iCiR,1,anynans),];
  }
  #Plot image
  heatmap.2(x_iCiR,
            Rowv=Rowv,
            Colv=Colv,
            dendrogram=dendrogram,
            trace="none",
            col=rainbow(300)[200:1],
            na.color="#000000FF",
            denscol="#000000FF",
            density.info=density.info,
            labCol=labelsX,
            labRow=labelsY,
            breaks=seq(lims[1],lims[2],(lims[2]-lims[1])/200),
            xlab=names[2],
            ylab=names[1]);
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Shows in a sevillan plot 
sevillanplot <- function( pVals_iSiA ,
                          power_iSiA = NA,
                          iAs_plot = NA,
                          pBreaks_iP = c( 0.05 , 0.005 , 0.0005 ),
                          extraFactors_iSiF = NA ){
  
  # Required libraries
  require( pracma );
  require( ggplot2 );
  require( grid );
  
  # Check inputs
  if( length( dim( pVals_iSiA ) ) != 2 ){
    stop( "The input parameter 'pVals_iSiA' must have 2 dimentions.");
  }
  #   if( length( pBreaks_iP ) != 3 ){
  #     stop( "Three breaks required" );
  #   }
  if( is.na( iAs_plot )[1] ){
    iAs_plot <- 1:dim( pVals_iSiA )[2];
  }
  if( is.na( power_iSiA ) ){
    power_iSiA <- ones( dim( pVals_iSiA )[1], dim( pVals_iSiA )[2] );
  }
  
  # Create names for breaks
  pBreaksNames_iP <- format( pBreaks_iP , scientific = FALSE );
  pBreaksNames_iP <- paste( ">" , 
                            format( pBreaksNames_iP , scientific = FALSE ) , 
                            sep = "" );
  pBreaksNames_iP[ length(pBreaksNames_iP) + 1 ] <- 
    paste( "<=" , 
           format( pBreaks_iP[ length(pBreaks_iP) ] , scientific = FALSE ) , 
           sep = "" );
  #pBreaksNames_iP[5] <- "too few samples";
  #pBreaksNames_iP[6] <- "NaN";
  pBreaksNames_iP[5] <- "Not sampled";
  
  # Create dot locations and colours
  mg_F <- meshgrid( x = 1:size( pVals_iSiA[ , iAs_plot ] )[2],
                    y = 1:size( pVals_iSiA[ , iAs_plot ] )[1] );
  factorise <- function(x){
    if( is.na(x) ){ return( pBreaksNames_iP[5] ); };
    if( is.nan(x) ){ return( pBreaksNames_iP[5] ); };
    if( x>pBreaks_iP[1] ){ return( pBreaksNames_iP[1] ); };
    if( x>pBreaks_iP[2] ){ return( pBreaksNames_iP[2] ); };
    if( x>pBreaks_iP[3] ){ return( pBreaksNames_iP[3] ); };
    if( x<=pBreaks_iP[3] ){ return( pBreaksNames_iP[4] ); };
  }
  df <- data.frame( 
    x = as.vector( mg_F$X ),
    y = as.vector( mg_F$Y ),
    colour = sapply( as.numeric( pVals_iSiA[ ,iAs_plot] ), factorise ),
    size = as.vector( as.numeric( power_iSiA[ ,iAs_plot] ) ) );
  if( any( !is.na( extraFactors_iSiF ) ) ){
    df <- cbind( df, extraFactors_iSiF );
  }
  limsSize <- c( 1.0, 3.0 );
  df$size[ df$size < limsSize[1] ] <- limsSize[1]; #cut limits to avoid problems
  df$size[ df$size > limsSize[2] ] <- limsSize[2]; #cut limits to avoid problems
  for( iA in 1:length( dimnames( pVals_iSiA )[[2]] ) ){
    df[ , dimnames( pVals_iSiA )[[2]][iA] ] <- pVals_iSiA[ df$y , iA ];
  }
  
  # Add optionals
  if( unique( power_iSiA[1:numel(power_iSiA)] ) > 1 ){
    gSize <- scale_size( name = "e-size with 80% power",
                         range = c( 7, 3 ),
                         limits = limsSize );
  }else{
    gSize <- NULL;
  }
  
  # Plot
  g <- 
    ggplot( df ) + 
    aes( x = x , 
         y = y , 
         colour = colour ) + 
    geom_point( size = 5 ) +
    scale_colour_manual( name = "p-value",
                         values = c( hsv(0,0,1.0), hsv(0,0.2,0.95), hsv(0,0.5,0.9), hsv(0,1,0.8), hsv(1,0,0.8), hsv(1,0,0.8) ), 
                         breaks = pBreaksNames_iP,
                         limits = pBreaksNames_iP ) +
    gSize +
    theme( axis.text.x = element_text( angle=0, vjust=0.5, hjust=0.5, size=12 ),
           axis.text.y = element_text( angle=0, vjust=0.5, hjust=1, size=12 ),
           panel.margin = unit( 1, "lines" ) ) +
    scale_x_continuous("", 
                       breaks= 1:dim( pVals_iSiA[ , iAs_plot] )[2],
                       labels= dimnames( pVals_iSiA[ , iAs_plot] )[[2]] ) +
    scale_y_continuous("",
                       breaks= 1:dim( pVals_iSiA )[1],
                       labels= dimnames( pVals_iSiA )[[1]] ); g
  #   scale_x_continuous("", 
  #                      breaks= 1:dim( pVals_iSiA[ , iAs_plot] )[2],
  #                      labels= dimnames( pVals_iSiA[ , iAs_plot] )[[2]], 
  #                      limits= c( 0.5, dim( pVals_iSiA[ , iAs_plot] )[2] + 0.5 ) ) +
  #     scale_y_continuous("",
  #                        breaks= 1:dim( pVals_iSiA )[1],
  #                        labels= dimnames( pVals_iSiA )[[1]], 
  #                        limits= c( 0.5, dim( pVals_iSiA )[1] + 0.5 ) );
  
  # Return results
  return( g );
}
# 
# # Shows in a sevillan plot 
# sevillanplot <- function( pVals_iSiA ,
#                           power_iSiA = NA,
#                           iAs_plot = NA,
#                           pBreaks_iP = c( 0.05 , 0.005 , 0.0005 ) ){
#   
#   # Required libraries
#   require( pracma );
#   require( ggplot2 );
#   require( grid );
#   
#   # Check inputs
#   if( length( dim( pVals_iSiA ) ) != 2 ){
#     stop( "The input parameter 'pVals_iSiA' must have 2 dimentions.");
#   }
#   if( length( pBreaks_iP ) != 3 ){
#     stop( "Three breaks required" );
#   }
#   if( is.na( iAs_plot )[1] ){
#     iAs_plot <- 1:dim( pVals_iSiA )[2];
#   }
#   if( is.na( power_iSiA ) ){
#     power_iSiA <- ones( dim( pVals_iSiA )[1], dim( pVals_iSiA )[2] );
#   }
#   
#   # Create names for breaks
#   pBreaksNames_iP <- format( pBreaks_iP , scientific = FALSE );
#   pBreaksNames_iP <- paste( ">" , 
#                             format( pBreaksNames_iP , scientific = FALSE ) , 
#                             sep = "" );
#   pBreaksNames_iP[ length(pBreaksNames_iP) + 1 ] <- 
#     paste( "<=" , 
#            format( pBreaks_iP[ length(pBreaks_iP) ] , scientific = FALSE ) , 
#            sep = "" );
#   pBreaksNames_iP[5] <- "too few samples";
#   #pBreaksNames_iP[6] <- "NaN";
#   
#   # Create dot locations and colours
#   mg_F <- meshgrid( x = 1:size( pVals_iSiA[ , iAs_plot ] )[2],
#                     y = 1:size( pVals_iSiA[ , iAs_plot ] )[1] );
#   factorise <- function(x){
#     if( is.na(x) ){ return( pBreaksNames_iP[5] ); };
#     if( is.nan(x) ){ return( pBreaksNames_iP[6] ); };
#     if( x>pBreaks_iP[1] ){ return( pBreaksNames_iP[1] ); };
#     if( x>pBreaks_iP[2] ){ return( pBreaksNames_iP[2] ); };
#     if( x>pBreaks_iP[3] ){ return( pBreaksNames_iP[3] ); };
#     if( x<pBreaks_iP[3] ){ return( pBreaksNames_iP[4] ); };
#   }
#   df <- data.frame( 
#     x = as.vector( mg_F$X ),
#     y = as.vector( mg_F$Y ),
#     colour = sapply( as.numeric( pVals_iSiA[ ,iAs_plot] ), factorise ),
#     size = as.vector( as.numeric( power_iSiA[ ,iAs_plot] ) ) );
#   limsSize <- c( 1.0, 3.0 );
#   df$size[ df$size < limsSize[1] ] <- limsSize[1]; #cut limits to avoid problems
#   df$size[ df$size > limsSize[2] ] <- limsSize[2]; #cut limits to avoid problems
#   for( iA in 1:length( dimnames( pVals_iSiA )[[2]] ) ){
#     df[ , dimnames( pVals_iSiA )[[2]][iA] ] <- pVals_iSiA[ df$y , iA ];
#   }
#   
#   # Plot
#   g <- 
#     ggplot( df ) + 
#     aes( x = x , 
#          y = y , 
#          colour = colour,
#          size = size ) + 
#     geom_point( ) +
#     scale_colour_manual( name = "p-value",
#                          values = c( hsv(0,0,1.0), hsv(0,0.2,0.95), hsv(0,0.5,0.9), hsv(0,1,0.8), hsv(1,0,0.8), hsv(1,0,0.8) ), 
#                          breaks = pBreaksNames_iP,
#                          limits = pBreaksNames_iP ) +
#     scale_size( name = "e-size with 80% power",
#                 range = c( 7, 3 ),
#                 limits = limsSize ) +
#     theme( axis.text.x = element_text( angle=0, vjust=0.5, hjust=0.5, size=12 ),
#            axis.text.y = element_text( angle=0, vjust=0.5, hjust=1, size=12 ),
#            panel.margin = unit( 1, "lines" ) ) +
#     scale_x_continuous("", 
#                        breaks= 1:dim( pVals_iSiA[ , iAs_plot] )[2],
#                        labels= dimnames( pVals_iSiA[ , iAs_plot] )[[2]], 
#                        limits= c( 0.5, dim( pVals_iSiA[ , iAs_plot] )[2] + 0.5 ) ) +
#     scale_y_continuous("",
#                        breaks= 1:dim( pVals_iSiA )[1],
#                        labels= dimnames( pVals_iSiA )[[1]], 
#                        limits= c( 0.5, dim( pVals_iSiA )[1] + 0.5 ) );
#   
#   # Return results
#   return( g );
# }