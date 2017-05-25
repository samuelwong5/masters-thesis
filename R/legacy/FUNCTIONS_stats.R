 

#
# Calculates statistical significance of correlations
#
testCorr = function( x_iSiX,
                     y_iSiY,
                     corrs_iXiY ){
  
  # Get dimentions
  iSnum <- dim(x_iSiX)[1]
  iXnum <- length(x_iSiX)
  iYnum <- length(y_iSiY)
  min_iXYnum <- min(iXnum, iYnum)
  m <- iSnum - 3/2 - (iXnum + iYnum)/2
  
  # Integrate
  ev_iXiY <- (1 - corrs_iXiY^2)
  w_iS <- rev(cumprod(rev(ev_iXiY)))
  
  # initialize
  d1 <- d2 <- f <- vector("numeric", min_iXYnum)
  
  # Compute stuff
  for (i in 1:min_iXYnum) {
    s <- sqrt((iXnum^2 * iYnum^2 - 4)/(iXnum^2 + iYnum^2 - 5))
    si <- 1/s
    d1[i] <- iXnum * iYnum
    d2[i] <- m * s - iXnum * iYnum/2 + 1
    r <- (1 - w_iS[i]^si)/w_iS[i]^si
    f[i] <- r * d2[i]/d1[i]
    iXnum <- iXnum - 1
    iYnum <- iYnum - 1
  }
  
  # Computer more stuff
  pv <- pf(f, d1, d2, lower.tail = FALSE)
  (dmat <- cbind(WilksL = w_iS, F = f, df1 = d1, df2 = d2, iXnum = pv));
  return( dmat );
  
}