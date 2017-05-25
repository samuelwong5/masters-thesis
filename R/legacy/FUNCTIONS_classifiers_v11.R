#
# A variety of classifications techniques
#
# v2 = Getting preselection step out of the stepForwards loop
# v3 = Using absolute testing set
# v4 = Trying to solve the potentia bug
# v5 = Trying to make deep step forwards
# v6 = Tiding up the classifier and feature selector
# v7 = Trying to fix it (last version gives change accuracy on testing set)
# v8 = Adding step forwards
# v9 = Also returning features
# v10 = In "classifier", using the best list for each number of features. Samples can be balanced
# v11 = In "systematicClassifiers", if only 1 crossvalidation fold is requested, the classifiers themselves are returned
#

library( pracma )

systematicClassifiers = function( dData_iSiG,
                                  bClas_iS,
                                  sClassifiers_iC = c( "lda" ),
                                  sPresselecto_iP = c( "LASSO" ),
                                  iNumFolds = 10,
                                  iNumFeatures = 10,
                                  iNumRepetitions = 10 ){
  #
  # Calculate accuracy for each combination of classifier and preselector.
  # It uses crossvalidation with "iNumFolds" folds.
  # iNumFeatures = The maximun number of features selected by the preselector
  # If "iNumFolds" is 1, all the data will be used as training and testing (i.e., no crossvalidation).
  # In addition, selected features will be returned
  #
  
  # Check input
  source( "FUNCTIONS_utils.R" );
  stopifnot( dim( dData_iSiG )[1] == length( bClas_iS ) );
  stopifnot( length( unique( bClas_iS ) ) == 2 );
  stopifnot( class( bClas_iS ) == "logical" );
  stopifnot( iNumFolds > 0 )
  
  # Shuffle data
  iSs <- randperm( 1:dim( dData_iSiG )[1] );
  dData_iSiG <- dData_iSiG[iSs, ]
  bClas_iS <- bClas_iS[iSs]
  
  # Divide population for crossvalidation, balancing classes
  iFold_iS <- array( data = NA,
                     dim = c( dim( dData_iSiG )[1] ) );
  iFold_iS[ bClas_iS ] <- ( 1 : sum( bClas_iS ) - 1 ) /  sum( bClas_iS );
  iFold_iS[ !bClas_iS ] <- ( 1 : sum( !bClas_iS ) - 1 ) /  sum( !bClas_iS );
  iFold_iS <- floor( iFold_iS * iNumFolds );
  stopifnot( diff( range( table( iFold_iS[ bClas_iS ] ) ) ) < 2 );
  stopifnot( diff( range( table( iFold_iS[ !bClas_iS ] ) ) ) < 2 );
  stopifnot( sort( unique( iFold_iS ) ) == 1 : iNumFolds - 1 )
  bTrain_iSiF <- !sapply( unique( iFold_iS ),
                          function( iF ){
                            return( iF == iFold_iS );
                          } );
  stopifnot( dim( bTrain_iSiF )[2] == iNumFolds )
  stopifnot( all( rowSums( bTrain_iSiF ) == iNumFolds - 1 ) )
  
  # If only 1 fold, all sample
  if( iNumFolds == 1 ){
    print( p( "iNumFolds is 1. No crossvalidation will be applied. " ) );
    bTrain_iSiF[ , ] <- TRUE;
    bTests_iSiF <- bTrain_iSiF;
  }else{
    bTests_iSiF <- !bTrain_iSiF;
    stopifnot( all( colSums( bTrain_iSiF ) > 0 ) );
  }
  
  # Prepare parallelism
  library( parallel );
  cl <- makeCluster( getOption( "cl.cores", 
                                4 ) );
  crossValidation = function( iF ){
    r <- classifier( dData_iSiG = dData_iSiG,
                     bClas_iS = bClas_iS,
                     bTrain_iS = bTrain_b_iSiF[ ,iF],
                     bTests_iS = bTests_b_iSiF[ ,iF],
                     sClassifier = sClassifiers_iC[iC],
                     sPresselect = sPresselecto_iP[iP],
                     iNumFeatures = iNumFeatures );
    return( r );
  }
  library( MASS );
  library( e1071 );
  clusterExport( cl, 
                 c( "classifier",
                    "bTrain_iSiF",
                    "bTests_iSiF",
                    "dData_iSiG",
                    "bClas_iS",
                    "sClassifiers_iC",
                    "sPresselecto_iP",
                    "iNumFeatures",
                    "lda",
                    "qda",
                    "svm",
                    "p",
                    "stepforwards",
                    "glmnetRank_2F",
                    "randperm" ),
                 envir = environment() );
  
  # Run crossvalidation in parallel for all combinations "preselector X classifier"
  accTr_iPiCiG <- array( data = NA,
                         dim = c( length( sPresselecto_iP ),
                                  length( sClassifiers_iC ),
                                  iNumFeatures ),
                         dimnames = list( preselector = sPresselecto_iP,
                                          classifier = sClassifiers_iC,
                                          numFeatures = 1:iNumFeatures ) );
  accTe_iPiCiG <- accTr_iPiCiG;
  features_jFiG <- NULL
  classifi_jF <- NULL
  for( iP in 1:length( sPresselecto_iP ) ){
    for( iC in 1:length( sClassifiers_iC ) ){
      
      # Hale the user
      print( p( "Running crossvalidation with preselector iP ", iP, " and classifier iC ", iC ) );
      
      # Run several times randomly balancing the samples
      accTr_iPiCiG[iP,iC, ] <- 0.0;
      accTe_iPiCiG[iP,iC, ] <- 0.0;
      for( iB in 1:iNumRepetitions ){
        
        # Hale the user
        print( p( "Repetition iB ", iB, " of ", iNumRepetitions ) );
        
        # Balance classes in a random fashion
        bTrain_b_iSiF <- bTrain_iSiF;
        bTests_b_iSiF <- bTests_iSiF;
        if( iNumRepetitions > 1 ){
          # Function to claculte who we have to use to balance the classes in "bClass0_iS"
          balance = function( bClass0_iS ){
            iSnum <- min( table( bClass0_iS ) );
            iSs_true <- randperm( which( bClass0_iS ) );
            iSs_false <- randperm( which( !bClass0_iS ) );
            bUse_iS <- rep( FALSE, 
                            length( bClass0_iS ) );
            bUse_iS[ iSs_true[ 1:iSnum ] ] <- TRUE;
            bUse_iS[ iSs_false[ 1:iSnum ] ] <- TRUE;
            return( bUse_iS )
          }
          # Balance all divsiions of the crossvalidation
          for( iF in 1:dim( bTrain_iSiF )[2] ){
            bTrain_b_iSiF[ bTrain_iSiF[ ,iF], iF ] <- balance( bClas_iS[ which( bTrain_iSiF[ bTrain_iSiF[ ,iF], iF ] ) ] )
            bTests_b_iSiF[ bTests_iSiF[ ,iF], iF ] <- balance( bClas_iS[ which( bTests_iSiF[ bTests_iSiF[ ,iF], iF ] ) ] )
            stopifnot( all( which( bTrain_b_iSiF[ ,iF] ) %in% which( bTrain_iSiF[ ,iF] ) ) );
            stopifnot( all( which( bTests_b_iSiF[ ,iF] ) %in% which( bTests_iSiF[ ,iF] ) ) );
            stopifnot( numel( intersect( which( bTrain_b_iSiF[ ,iF] ), which( bTests_iSiF[ ,iF] ) ) ) == 0 );
            stopifnot( numel( intersect( which( bTests_b_iSiF[ ,iF] ), which( bTrain_iSiF[ ,iF] ) ) ) == 0 );
            stopifnot( numel( intersect( which( bTests_b_iSiF[ ,iF] ), which( bTrain_b_iSiF[ ,iF] ) ) ) == 0 );
            stopifnot( 0 == diff( table( bClas_iS[ which( bTrain_b_iSiF[ bTrain_iSiF[ ,iF], iF ] ) ] ) ) );
            stopifnot( 0 == diff( table( bClas_iS[ which( bTests_b_iSiF[ bTests_iSiF[ ,iF], iF ] ) ] ) ) );
          } 
        }                     
        
        # Send classifier to the cores
        clusterExport( cl, 
                       c( "iP",
                          "iC",
                          "bTrain_b_iSiF",
                          "bTests_b_iSiF" ),
                       envir = environment() );      
        
        # Run crossvalidation in parallel
        cross_iF <- parLapply( cl,
                               X = 1:dim( bTrain_iSiF )[2],
                               fun = crossValidation );
        #cross_iF <- vector( mode = "list",
        #                    length  = dim( bTrain_iSiF )[2] )
        #for( iF in 1:dim( bTrain_iSiF )[2] ){
        #  cross_iF[iF] <- crossValidation(iF)
        #}
        
        # Collect output
        acc_iGiF  <- sapply( cross_iF, 
                             function( s ) {
                               s$accTr_iG 
                             } );
        accTr_iPiCiG[iP,iC, ] <- accTr_iPiCiG[iP,iC, ] + rowMeans( acc_iGiF ) / iNumRepetitions;
        acc_iGiF  <- sapply( cross_iF, 
                             function( s ) {
                               s$accTe_iG 
                             } );
        accTe_iPiCiG[iP,iC, ] <- accTe_iPiCiG[iP,iC, ] + rowMeans( acc_iGiF ) / iNumRepetitions;
        
        # Collect features if no crossvalidation was applied
        if( iNumFolds == 1 ){
          stopifnot( length( sClassifiers_iC ) == 1 ); #Only 1 classifier can be used if you want to extract the selected features
          stopifnot( length( sPresselecto_iP ) == 1 ); #Only 1 preselector can be used if you want to extract the selected features
          stopifnot( length( cross_iF ) == 1 )
          features_jFiG <- cross_iF[[1]]$iGpreselect_jFiG;
          classifi_jF <- cross_iF[[1]]$classifier_jF;
        }
        
      }
    }
    
  }
  
  # Return results
  stopCluster( cl )
  r <- list( accTr_iPiCiG = accTr_iPiCiG,
             accTe_iPiCiG = accTe_iPiCiG ,
             features_jFiG = features_jFiG  ,
             classifi_jF = classifi_jF );
  return( r )
}


stepforwards = function( x_iSiG,# = as.matrix( dData_iSiG[bTrain_iS, ] ),
                         y_iS,# = bClas_iS[bTrain_iS] * 1,
                         sClassifier,
                         iNumFeatures ){
  #
  # Ranks the features of "x_iSiG" using the step forwards algorithm for predicting "y_iS"
  # The used classifier will be "sClassifier", while the total number of features selected will be "iNumFeatures"
  #
  
  # Check input
  stopifnot( dim( x_iSiG )[1] == length( y_iS ) );
  stopifnot( iNumFeatures <= dim( x_iSiG )[2] )
  
  # Select best feature one by one
  iGs_best_iF <- array( data = NA,
                        dim = c( iNumFeatures ) );
  for( iF in 1:iNumFeatures ){
    
    # Measrue accuracy of all features
    acc_iG <- -1 * ones( dim( x_iSiG )[2],
                         1 );
    for( iG in 1:dim( x_iSiG )[2] ){
      
      # Skip if this feature was already selected
      if( iG %in% iGs_best_iF ){
        next;
      }
      
      # Train
      iGs_use <- c( iGs_best_iF[ !is.na( iGs_best_iF ) ],
                    iG )
      if( sClassifier == "lda" ){
        classif_F <- lda( x = as.matrix( x_iSiG[ , iGs_use ] ),
                          grouping = y_iS );
      }else if( sClassifier == "qda" ){
        classif_F <- qda( x = as.matrix( x_iSiG[ , iGs_use ] ),
                          grouping = y_iS );
      }else if( sClassifier == "svm" ){
        classif_F <- svm(  x = as.matrix( x_iSiG[ , iGs_use ] ),
                           y = as.numeric( y_iS ) );
      }else{
        stop( p( "The classifier method ", methodClassifier, " is not recognised." ) );
      }
      
      # Test
      pre_F <- predict( object = classif_F,
                        newdata = as.matrix( x_iSiG[ , iGs_use ] ) );
      
      # Calculate accuracy
      if( sClassifier == "svm" ){
        pre2_F <- round( pre_F );
      }else{
        pre2_F <- pre_F$class;
      }
      stopifnot( length( pre2_F ) == length( y_iS ) );
      acc_iG[iG] <- mean( pre2_F == y_iS,
                          na.rm = TRUE );
      
    }
    
    # Select the feature with highest accuracy
    r <- sort( acc_iG,
               index.return = TRUE,
               decreasing = TRUE );
    iGs_best_iF[iF] <- r$ix[1];
    
  }
  
  # Return results
  r <- list()
  if( length( dimnames( x_iSiG )[[2]] ) > 1 ){
    r$nGs_ranked = iGs_best_iF
  }else{
    r$nGs_ranked = dimnames( x_iSiG )[[2]][ iGs_best_iF ];
  }
  return( r )
  
}

glmnetRank_2F = function( glmnet, first=T, names=T,...) {
  # Redefinition of glmnetRank (it has a bug in the repository!)
  beta <- glmnet$beta[, order(glmnet$lambda, decreasing=T)]
  if (first == T) {
    # Return variables that pops-up first
    coef <- unique(unlist(apply(beta, 2, function(x) {
      o <- order(abs(x), decreasing=T)
      return (o[abs(x[o]) > 0])
    })))
    nrcoef = length(coef)
    coef <- c(coef, setdiff(1:nrow(glmnet$beta), coef))
  } else {
    coef <- order(abs(beta[, ncol(beta)]), decreasing=T)
    nrcoef = length(coef)
  }
  if (names) {
    # Return variable names that are non-zero at the end
    coef <- rownames(glmnet$beta)[coef]
  }
  glmnetR=list()
  glmnetR$coef = coef
  glmnetR$nrcoef = nrcoef
  return (glmnetR)
}

classifier = function( dData_iSiG,
                       bClas_iS,
                       bTrain_iS,
                       bTests_iS,
                       sClassifier = "lda",
                       sPresselect = "LASSO",
                       iNumFeatures = iNumFeatures,
                       balance = FALSE ){
  #
  # Classifies data according to the specified classifier and the specified preselector
  #
  
  # Check intputs
  stopifnot( iNumFeatures < dim( dData_iSiG )[2] );
  stopifnot( length( bClas_iS ) == length( bTrain_iS ) )
  stopifnot( length( bClas_iS ) == length( bTests_iS ) )
  stopifnot( length( bClas_iS ) == dim( dData_iSiG )[1] )
  #stopifnot( all( xor( bTrain_iS, 
  #                     bTests_iS ) ) )
  
  # Balance classes
  if( balance ){
    # Auxiliary function
    balance = function( bUse_iS,
                        bClas_iS ){
      # Limit "bUse_iS" bot "bClas_iS[bUse_iS]" to be balanced
      iSnum_use <- min( table( bClas_iS[ bUse_iS ] ) )
      iSs_all <- randperm( which( bUse_iS ) )
      iSs_true <- which( bUse_iS & bClas_iS )
      iSs_fals <- which( bUse_iS & !bClas_iS )
      if( length( iSs_true ) > 1 ){
        iSs_true <- randperm( iSs_true )
      }
      if( length( iSs_fals ) > 1 ){
        iSs_fals <- randperm( iSs_fals )
      }
      stopifnot( length( intersect( iSs_true, iSs_fals ) ) == 0 )
      stopifnot( all( intersect( iSs_true, iSs_fals ) %in% iSs_all ) )
      iSs_true <- iSs_true[ 1:iSnum_use ]
      iSs_fals <- iSs_fals[ 1:iSnum_use ]
      bUse_iS[ ] <- FALSE
      bUse_iS[ c( iSs_true, iSs_fals ) ] <- TRUE
      return( bUse_iS )
    }
    # Balance classes
    iSnum_train <- min( table( bClas_iS[ bTrain_iS ] ) )
    stopifnot( iSnum_train > 5 )
    bTrain_iS <- balance( bUse_iS = bTrain_iS, 
                          bClas_iS )
    bTests_iS <- balance( bUse_iS = bTests_iS, 
                          bClas_iS )
    # Check balance is OK
    stopifnot( sum( bTests_iS ) >= 2 )
    stopifnot( length( unique( table( bClas_iS[ bTrain_iS ] ) ) ) == 1 )
    stopifnot( length( unique( table( bClas_iS[ bTests_iS ] ) ) ) == 1 )
  }
  
  # Preselect features
  library( SurvRank )
  library( glmnet )
  library( pracma )
  iGpreselect_jFiG <- vector( mode = "list")
  for( iF in 1:iNumFeatures ){
    if( sPresselect == "LASSO" ){
      # Preselect with LASSO
      l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                     y = bClas_iS[bTrain_iS] * 1,
                     family = "binomial",
                     alpha = 1,
                     dfmax = iF );
      r <- glmnetRank_2F( l_F,
                          first = FALSE );
      iGpreselect_jFiG[[iF]] = r$coef[1:iF];
    }else if( sPresselect == "RIDGE" ){
      # Preselect with RIDGE
      l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                     y = bClas_iS[bTrain_iS] * 1,
                     family = "binomial",
                     alpha = 0.001,
                     dfmax = iF );
      r <- glmnetRank_2F( l_F,
                          first = FALSE );
      iGpreselect_jFiG[[iF]] = r$coef[1:iF];
    }else if( toupper( sPresselect ) == "FORWARDS" ){  
      r <- stepforwards( x_iSiG = as.matrix( dData_iSiG[bTrain_iS, ] ),
                         y_iS = bClas_iS[bTrain_iS] * 1,
                         sClassifier = sClassifier,
                         iNumFeatures = iF );
      iGpreselect_jFiG[[iF]] = dimnames( dData_iSiG )[[2]][ r$nGs_ranked ];
    }else if( toupper( sPresselect ) == "NONE" ){  
      l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                     y = bClas_iS[bTrain_iS] * 1,
                     family = "binomial",
                     alpha = 1,
                     dfmax = iF );
      r <- glmnetRank_2F( l_F,
                          first = FALSE );
      iGpreselect_jFiG[[iF]] = randperm( r$coef )[1:iF];
      #iGpreselect_jFiG[[iF]] = r$coef[1:iF];
      #     if( is.null( dimnames( dData_iSiG )[[2]] ) ){
      #       iGpreselect_iG = randperm( 1:dim( dData_iSiG )[2] );
      #     }else{
      #       iGpreselect_iG = randperm( dimnames( dData_iSiG )[[2]] );
      #     }
    }else{
      stop( p( "Preselector ", sPresselect, " does not exist." ) );
      break;
    }
  }
  #iGpreselect_jFiG[[1]] <- iGpreselect_jFiG[[2]][1]
  
  # Measure accuracy for each feature set
  library( MASS );
  library( e1071 );
  accTr_iG <- array( data = NA,
                     dim = c( length( iGpreselect_jFiG ) ) );
  accTe_iG <- accTr_iG;
  classif_jF <- list( mode = "list",
                      length = length( iGpreselect_jFiG ) )
  for( iF in 1:length( iGpreselect_jFiG ) ){
    
    # Train
    if( sClassifier == "lda" ){
      classif_jF[[iF]] <- lda( x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_jFiG[[iF]] ] ),
                             grouping = bClas_iS[ bTrain_iS ] );
    }else if( sClassifier == "qda" ){
      classif_jF[[iF]] <- qda( x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_jFiG[[iF]] ] ),
                             grouping = bClas_iS[ bTrain_iS ] );
    }else if( sClassifier == "svm" ){
      classif_jF[[iF]] <- svm(  x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_jFiG[[iF]] ] ),
                              y = as.numeric( bClas_iS[ bTrain_iS ] ) );
    }else{
      stop( p( "The classifier method ", methodClassifier, " is not recognised." ) );
    }
    
    # Test
    preTr_F <- predict( object = classif_jF[[iF]],
                        newdata = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_jFiG[[iF]] ] ) );
    preTe_F <- predict( object = classif_jF[[iF]],
                        newdata = as.matrix( dData_iSiG[ bTests_iS, iGpreselect_jFiG[[iF]] ] ) );
    
    # Calculate accuracy
    if( sClassifier == "svm" ){
      preTe2_F <- round( preTe_F );
      preTr2_F <- round( preTr_F );
    }else{
      preTe2_F <- preTe_F$class;
      preTr2_F <- preTr_F$class;
    }
    stopifnot( length( preTr2_F ) == length( bClas_iS[ bTrain_iS ] ) );
    stopifnot( length( preTe2_F ) == length( bClas_iS[ bTests_iS ] ) );
    accTr_iG[iF] <- mean( preTr2_F == bClas_iS[ bTrain_iS ],
                          na.rm = TRUE );
    accTe_iG[iF] <- mean( preTe2_F == bClas_iS[ bTests_iS ],
                          na.rm = TRUE );
    
  }
  
  # Return result
  r$accTr_iG <- accTr_iG;
  r$accTe_iG <- accTe_iG;
  r$iGpreselect_jFiG <- iGpreselect_jFiG;
  r$classifier_jF <- classif_jF;
  return( r );
  
}
