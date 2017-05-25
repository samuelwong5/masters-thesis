#
# Functions to do some more advanced stats
#
# v1 = From scratch
# V2 = The single-GLM step is now paralellised
# v3 = Using RLM
# v4 = Attempting speed GLM
# v5 = Adding a simplified version of the systematic GLM
# v6 = systematicGLM_simple generalised to accept NAs and to extract also levels of iX when iX is a factor
# v7 = Small bug corrected in "systematicGLM_simple" when no "nKs" were used
# v8 = Compatibility improvements in "systematicGLM_simple"
# v9 = Paralelissing "systematicGLM_simple"
#
#For installing further libraries: 
#source("http://bioconductor.org/biocLite.R")
#biocLite("wathever")
#library(wathever)


# library( pracma )
# x_iSiX <- rand( 1000, 100 ) - 0.5;
# covs_iSiC <- rand( 1000, 3 ) - 0.5;
# y_iS <- ( x_iSiX[ ,1:10] %*% t(t( 0.1 * (1:10) )) 
#           + covs_iSiC %*% t(t( 0.1 * (1:3) )) );
# dimnames( y_iS )[[2]] <- list( "y" );
# dimnames( covs_iSiC )[[2]] <- paste( "cov",
#                                      1:3,
#                                      sep = "");
# dimnames( x_iSiX )[[2]] <- paste( "x",
#                                   1:100,
#                                   sep = "");
# data_iSiX <- cbind( y_iS,
#                     covs_iSiC,
#                     x_iSiX );
# yName <- dimnames( y_iS )[[2]];
# covNames_iC <- dimnames( covs_iSiC )[[2]];
# varNames_iV <- dimnames( x_iSiX )[[2]];
# data_iSiX <- as.data.frame( data_iSiX );

# Import stuff
source( "FUNCTIONS_utils.R" );

cleanGLM = function( data_iSiC, 
                     nYs, 
                     nCs_eliminate = NULL ){
  # Cleans data by using a GLM model
  # WARNING = The function to keep the correlates of "nCs_keep" is still not implemented
  
  # Check input
  stopifnot( all( nYs %in% dimnames( data_iSiC )[[2]] ) );
  stopifnot( all( nCs_eliminate %in% dimnames( data_iSiC )[[2]] ) );
  
  # Make var names compatible
  dimnames( data_iSiC )[[2]] <- paste( "mekagoenlaputa",
                                       dimnames( data_iSiC )[[2]],
                                       sep = "" ); 
  nYs <- paste( "mekagoenlaputa",
                nYs,
                sep = "" );
  nCs_eliminate <- paste( "mekagoenlaputa",
                          nCs_eliminate,
                          sep = "" );
  
  # Clean all the variables
  dataClean_iSiY <- array( data = NA,
                           dim = c( dim( data_iSiC )[1], 
                                    length( nYs ) ),
                           dimnames = list( subject = dimnames( data_iSiC )[[1]],
                                            nYs = nYs ) );
  for( iY in 1:length( nYs ) ){
    
    # Hale the user
    if( mod( iY, 100 ) == 1 ){
      print( p( "Analysing variable iY ", iY, " of ", length( nYs ) ) );
    }
    
    # Select family
    if( class( data_iSiC[ ,nYs[iY]] ) == "factor" |
          class( data_iSiC[ ,nYs[iY]] ) == "character" ){
      family = "binomial";
    }else{
      family = "gaussian";
    }
    
    # Run model
    #     glmAll_F <- glm( formula = as.formula( paste( c( p( nYs[iY], " ~ " ), 
    #                                                      nCs_keep,
    #                                                      nCs_eliminate ),
    #                                                   collapse = " + " ) ),
    #                      family = family,
    #                      data = data_iSiC );
    #     glmKee_F <- glm( formula = as.formula( paste( c( p( nYs[iY], " ~ " ), 
    #                                                      nCs_keep ),
    #                                                   collapse = " + " ) ),
    #                      family = family,
    #                      data = data_iSiC );
    glmEli_F <- glm( formula = as.formula( paste( c( p( nYs[iY], " ~ " ), 
                                                     nCs_eliminate ),
                                                  collapse = " + " ) ),
                     family = family,
                     data = data_iSiC );
    dataClean_iSiY[ names( glmEli_F$residuals ), iY ] <- glmEli_F$residuals;
  }
  
  # Return
  dimnames( dataClean_iSiY )[[2]] <- gsub( "mekagoenlaputa",
                                           "", 
                                           dimnames( dataClean_iSiY )[[2]] );
  return( dataClean_iSiY );
}


systematicGLM_oneX = function( data_iSiC = data_iSiC, 
                               nYs, 
                               nX, 
                               nKs,
                               verbose = FALSE ){
  # A version of "systematicGLM_simple" for only one "nX"
  #
  
  # Eliminate samples that are NA or Inf in either iY or iX
  bIsNa_iS <- zeros( dim( data_iSiC )[1], 1 ) != 0
  if( sum( bIsNa_iS ) > 0 ){
    warning( p( "Eliminating ", sum( bIsNa_iS ), 
                " samples out of ", length( bIsNa_iS ),
                " that are NA/null/inf for nY ", nYs[iY], 
                " and nX ", nXs[nX] ) )
  }          
  
  # Skip if outcome is equal to one of the variables
  if( nYs[iY] == nX ){
    r <- list();
    r$p_iX <- c( );
    r$e_iX <- c( );
    r$t_iX <- c( );
    r$names_iX <- c( );
    return( r );
  }
  
  # Run GLM
  if( is.null( nKs ) ){
    glm_F <- glm( formula = as.formula( paste( c( p( nYs[iY], " ~ ", nX ) ),
                                               collapse = " + " ) ),
                  family = family,
                  data = data_iSiC[!bIsNa_iS, ] );
  }else{
    nKs_c <- setdiff( nKs,
                      c( nYs[iY],
                         nX ) );
    glm_F <- glm( formula = as.formula( paste( c( p( nYs[iY], " ~ ", nX ),
                                                  nKs_c ),
                                               collapse = " + " ) ),
                  family = family,
                  data = data_iSiC[!bIsNa_iS, ] );
  }
  
  # Extract results using vectors, just in case X had levels
  iCs_isX <- which( grepl( nX,
                           dimnames( coef( summary(glm_F) ) )[[1]] ) );
  stopifnot( length( iCs_isX ) > 0 );
  p_iX <- coef( summary(glm_F) )[ iCs_isX, 4 ];
  e_iX <- coef( summary(glm_F) )[ iCs_isX, 1 ];
  t_iX <- coef( summary(glm_F) )[ iCs_isX, 3 ];
  names_iX <- dimnames( coef( summary(glm_F) ) )[[1]][iCs_isX];
  
  # Return result
  r <- list();
  r$p_iX <- p_iX;
  r$e_iX <- e_iX;
  r$t_iX <- t_iX;
  r$names_iX <- names_iX;
  return( r );
  
}

systematicGLM_simple = function( data_iSiC, 
                                 nYs, 
                                 nXs, 
                                 nKs = NULL,
                                 verbose = FALSE ){
  # A simpler version of "systematicGLM"
  # This version simply systematically test for any effects of the "nXs" variables into the "nYs" variables, while controlling
  # for the variables listed in "nKs". These vectors ( nYs, nXs and nKs ) whould contain the names of teh variables.
  # The colums of "data_iSiC" should contain all these names.
  #
  
  # Check input
  stopifnot( all( nYs %in% dimnames( data_iSiC )[[2]] ) );
  stopifnot( all( nXs %in% dimnames( data_iSiC )[[2]] ) );
  stopifnot( all( nKs %in% dimnames( data_iSiC )[[2]] ) );
  
  # Prepare parallelism
  library( parallel );
  cl <- makeCluster( getOption( "cl.cores", 
                                3 ) );
  glmOneX = function( iX ){
    r <- systematicGLM_oneX( data_iSiC = data_iSiC, 
                             nYs = nYs, 
                             nX = nXs[iX], 
                             nKs = nKs,
                             verbose = FALSE );
    return( r );
  }
  library( MASS );
  library( e1071 );
  clusterExport( cl, 
                 c( "data_iSiC",
                    "nYs",
                    "nXs",
                    "nKs",
                    "systematicGLM_oneX",
                    "mod",
                    "p",
                    "zeros" ),
                 envir = environment() );  
  
  # Run GLMs
  p_iYiX <- array( data = NA,
                   dim = c( length( nYs ), 
                            0 ),
                   dimnames = list( nYs = nYs,
                                    nXs = c() ) );
  e_iYiX <- p_iYiX;
  t_iYiX <- p_iYiX; 
  for( iY in 1:length( nYs ) ){
    
    # Hale the user
    print( p( "Running iY ", iY, " of ", length( nYs ) ) )
    
    # Select family
    if( class( data_iSiC[ ,nYs[iY]] ) == "factor" |
        class( data_iSiC[ ,nYs[iY]] ) == "character" |
          class( data_iSiC[ ,nYs[iY]] ) == "logical" |
          length( unique( data_iSiC[ ,nYs[iY]] ) ) == 2 ){
      family = "binomial";
    }else{
      family = "gaussian";
    }
    
    # Send data to cores
    clusterExport( cl, 
                   c( "iY",
                      "family" ),
                   envir = environment() );  
    
    # Run glm in parallel
    cross_jF <- parLapply( cl,
                           X = 1:length( nXs ),
                           fun = glmOneX );
    extractField = function( z_jZ,
                             field ){
      field_iZ <- sapply( z_jZ,
                          function( z ){
                            return( z[field] );
                          } );
      field_iZ <- unlist( field_iZ );
      return( field_iZ );
    }
    p_iX <- extractField( cross_jF, "p_iX" );
    e_iX <- extractField( cross_jF, "e_iX" );
    t_iX <- extractField( cross_jF, "t_iX" );
    names_iX <- extractField( cross_jF, "names_iX" );
                    
    # Store results
    expandMatrix = function( p_iYiX,
                             names_iX ){
      bAdd_iX <- !sapply( names_iX,
                          function( nX ){
                            return( nX %in% dimnames( p_iYiX )[[2]] );
                          } );
      if( all( !bAdd_iX ) ){
        return( p_iYiX );
      }
      newNames_iX <- c( dimnames( p_iYiX )[[2]],
                        names_iX[ bAdd_iX ] );
      p_iYiX <- cbind( p_iYiX,
                       repmat( NaN,
                               dim( p_iYiX )[1],
                               sum( bAdd_iX ) ) )
      dimnames( p_iYiX )[[2]] <- newNames_iX;
      return( p_iYiX );
    }  
    p_iYiX <- expandMatrix( p_iYiX,
                            names_iX );
    e_iYiX <- expandMatrix( e_iYiX,
                            names_iX );
    t_iYiX <- expandMatrix( t_iYiX,
                            names_iX );
    p_iYiX[iY,names_iX] <- p_iX;
    e_iYiX[iY,names_iX] <- e_iX;
    t_iYiX[iY,names_iX] <- t_iX;
  }
  
  # Return results
  r <- list( p_iYiX = p_iYiX,
             e_iYiX = e_iYiX,
             t_iYiX = t_iYiX )
  return( r );
}

systematicGLM = function( data_iSiX,
                          yName,
                          covNames_iC,
                          varNames_iV ){
  # Builds a systematic GLM in 2 steps
  # 1) With many GLMs, finds all the variables in "x" that have significant effects on "y" when accounting for the covariants on "covs"
  # 2) It builds a single GLM with all the variables in "x" that were significant, plus "covs"
  # Then it returns the pvalues of the second GLM for each variable in "x"
  #
  
  # Check input
  stopifnot( yName %in% dimnames( data_iSiX )[[2]] );
  stopifnot( all( covNames_iC %in% dimnames( data_iSiX )[[2]] ) );
  stopifnot( all( varNames_iV %in% dimnames( data_iSiX )[[2]] ) );
  
  # Eliminate inapropriate chars
  eliminateUncompatibleChars <- function( str_iS ){
    str_iS <- gsub( "(\\'|\\[|\\]|\\:)",
                    "",
                    str_iS,
                    fixed = FALSE )
    return( str_iS );
  }
  dimnames( data_iSiX )[[2]] <- eliminateUncompatibleChars( dimnames( data_iSiX )[[2]] );
  yName <- eliminateUncompatibleChars( yName );
  covNames_r_iC <- eliminateUncompatibleChars( covNames_iC );
  varNames_r_iV <- eliminateUncompatibleChars( varNames_iV );
  
  # Prepare paralellisation
  library( MASS )
  library( stats )
  library( parallel )
  library( speedglm );
  cl <- makeCluster( getOption( "cl.cores", 
                                3 ) );
  clusterExport( cl, c( "yName",
                        "covNames_r_iC",
                        "varNames_r_iV",
                        "rlm",
                        "dt",
                        "speedglm" ),
                 envir = environment() );
  
  # Parallel function
  runGLM = function( iV ){
    
    # Run model
    glm_F <- speedglm( formula = as.formula( paste( yName, 
                                               " ~ ",
                                               paste( covNames_r_iC,
                                                      collapse = " + " ), 
                                               " + ",
                                               varNames_r_iV[iV] ) ),
                  data = data_p_iSiX );
    
    # Collect results
    coef_iRiC <- coef( summary( glm_F ) );
    iR_var <- grep( varNames_r_iV[iV],
                    dimnames( coef_iRiC )[[1]] );
    if( length( iR_var ) == 0 ){
      p = 1;
      e = NA;
    }else{
      iR_min <- which.min( coef( summary( glm_F ) )[iR_var,4] );
      e <-coef( summary( glm_F ) )[ iR_var[iR_min], 1 ];
      p <-coef( summary( glm_F ) )[ iR_var[iR_min], 4 ];
      p <- as.numeric( as.character( p ) );
#       iR_min <- which.min( coef_iRiC[ iR_var, 3 ] );
#       e <- coef_iRiC[ iR_var[iR_min], 1 ];
#       p <- coef_iRiC[ iR_var[iR_min], 3 ];
#       p <- dt( x = p,
#                df = dim( data_iSiX )[1]-2 );
    }
    
    # Return result
    ep <- c( e, p );
    names( ep ) <- c( "e", "p" );
    return( ep )
    
  }
  
  # Run function in paralell
  print( paste( "Running single GLM for iVnum ", length( varNames_iV ) ,"variables " ) );
  iVend_iV2 <- seq( 1, 
                    length( varNames_r_iV ),
                    50 );
  p_iV <- vector( mode = "numeric",
                  length = length( varNames_r_iV ) );
  e_iV <- p_iV;
  for( iV2 in 2:length( iVend_iV2 ) ){
    
    # Send needed data
    print( p( "Running iVs from ", iVend_iV2[iV2-1],
              " to ", iVend_iV2[iV2],
              " of ", length( varNames_r_iV ) ) );
    iVs <- iVend_iV2[iV2-1] : iVend_iV2[iV2];
    data_p_iSiX <- data_iSiX[ , c( yName,
                                   covNames_r_iC,
                                   varNames_r_iV[ iVs ] ) ];
    clusterExport( cl, c( "yName",
                          "covNames_r_iC",
                          "varNames_r_iV",
                          "data_p_iSiX" ),
                   envir = environment() );
    
    # Execute in paralell
    pS_iRiV <- parSapply( cl = cl, 
                          iVs, 
                          runGLM );
    p_iV[iVs] <- pS_iRiV["p", ];
    e_iV[iVs] <- pS_iRiV["e", ];
    
  }
  names( p_iV ) <- varNames_iV;
  names( e_iV ) <- varNames_iV;
  stopifnot( length( e_iV ) == length( varNames_iV ) );
  
  # Close clusters
  stopCluster( cl );
  
  # Now build big model with all significant vars
  pTh <- 0.01;
  bSignificant_iV <- p_iV < pTh;
  bSignificant_iV[ is.na( p_iV ) ] <- FALSE;
  print( paste( "Running full GLM for ",
                sum( bSignificant_iV ),
                " of ",
                length( bSignificant_iV ),
                " variables.",
                sep = "" ) );
  glm_F <- speedglm( formula = as.formula( paste( yName, 
                                             " ~ ",
                                             paste( covNames_r_iC,
                                                    collapse = " + " ), 
                                             " + ",
                                             paste( varNames_r_iV[ bSignificant_iV ],
                                                    collapse = " + " ) ) ),
                data = data_iSiX );
  
  # Collect the p-value of the full model
  pFull_iV <- p_iV;
  eFull_iV <- e_iV;
  for( iV in 1:length( varNames_iV ) ){
    if( bSignificant_iV[iV] ){
      coef_iRiC <- coef( summary( glm_F ) );
      iR_var <- grep( varNames_iV[iV],
                      dimnames( coef_iRiC )[[1]] );
      if( length( iR_var ) == 0 ){
        pFull_iV[iV] <- 1.0;
        eFull_iV[iV] <- NA;
      }else{
        stopifnot( length( iR_var ) > 0 )
        iR_min <- which.min( coef( summary( glm_F ) )[iR_var,4] );
        eFull_iV[iV] <-coef( summary( glm_F ) )[ iR_var[iR_min], 1 ];
        pFull_iV[iV] <-coef( summary( glm_F ) )[ iR_var[iR_min], 4 ];
        pFull_iV[iV] <- as.numeric( as.character( pFull_iV[iV] ) );
#         iR_min <- which.min( coef_iRiC[ iR_var, 3 ] );
#         eFull_iV[iV] <- coef_iRiC[ iR_var[iR_min], 1 ];
#         pFull_iV[iV] <- coef_iRiC[ iR_var[iR_min], 3 ];
#         pFull_iV[iV] <- dt( x = pFull_iV[iV],
#                             df = dim( data_iSiX )[1]-2 );
      }
    }
  }
  
  # Return
  pe_iViP <- cbind( pValue_singleGLM = p_iV,
                    pValue_fullGLM = pFull_iV,
                    effect_singleGLM = e_iV,
                    effect_fullGLM = eFull_iV );
  return( pe_iViP );
  
}





















