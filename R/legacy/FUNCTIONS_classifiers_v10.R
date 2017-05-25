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
  # If "" is 1, all the data will be used as training and testing (i.e., no crossvalidation).
  # In addition, selected features will be returned
  #
  
  # Check input
  source( "FUNCTIONS_utils.R" );
  stopifnot( dim( dData_iSiG )[1] == length( bClas_iS ) );
  stopifnot( length( unique( bClas_iS ) ) == 2 );
  stopifnot( class( bClas_iS ) == "logical" );
  stopifnot( iNumFolds > 0 )
  
  # Divide population for crossvalidation
  iFold_iS <- array( data = NA,
                     dim = c( dim( dData_iSiG )[1] ) );
  iFold_iS[ bClas_iS ] <- ( 1 : sum( bClas_iS ) - 1 ) /  sum( bClas_iS );
  iFold_iS[ !bClas_iS ] <- ( 1 : sum( !bClas_iS ) - 1 ) /  sum( !bClas_iS );
  iFold_iS <- floor( iFold_iS * iNumFolds );
  stopifnot( diff( range( table( iFold_iS ) ) ) < 5 );
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
                                3 ) );
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
                    "stepforwards" ),
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
  features_iPiCiG <- accTr_iPiCiG;
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
          for( iF in 1:dim(bTrain_iSiF )[2] ){
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
          stopifnot( length( cross_iF ) == 1 )
          features_iPiCiG[iP,iC, ] <- cross_iF[[1]]$iGpreselect_iG
        }
        
      }
    }
    
  }
  
  # Return results
  stopCluster( cl )
  r <- list( accTr_iPiCiG = accTr_iPiCiG,
             accTe_iPiCiG = accTe_iPiCiG ,
             features_iPiCiG = features_iPiCiG );
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


classifier = function( dData_iSiG,
                       bClas_iS,
                       bTrain_iS,
                       bTests_iS,
                       sClassifier = "lda",
                       sPresselect = "LASSO",
                       iNumFeatures = iNumFeatures ){
  #
  # Classifies data according to the specified classifier and the specified preselector
  #
  
  # Check intputs
  stopifnot( iNumFeatures < dim( dData_iSiG )[2] );
  
  # Preselect features
  library( SurvRank )
  library( glmnet )
  library( pracma )
  if( sPresselect == "LASSO" ){
    # Preselect with LASSO
    l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                   y = bClas_iS[bTrain_iS] * 1,
                   family = "binomial",
                   alpha = 1 );
    r <- glmnetRank( l_F );
    iGpreselect_iG = r$coef;
  }else if( sPresselect == "RIDGE" ){
    # Preselect with RIDGE
    l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                   y = bClas_iS[bTrain_iS] * 1,
                   family = "binomial",
                   alpha = 0 );
    r <- glmnetRank( l_F );
    iGpreselect_iG = r$coef;
  }else if( toupper( sPresselect ) == "FORWARDS" ){  
    r <- stepforwards( x_iSiG = as.matrix( dData_iSiG[bTrain_iS, ] ),
                       y_iS = bClas_iS[bTrain_iS] * 1,
                       sClassifier = sClassifier,
                       iNumFeatures = iNumFeatures );
    iGpreselect_iG = dimnames( dData_iSiG )[[2]][ r$nGs_ranked ];
  }else if( toupper( sPresselect ) == "NONE" ){  
    l_F <- glmnet( x = as.matrix( dData_iSiG[bTrain_iS, ] ),
                   y = bClas_iS[bTrain_iS] * 1,
                   family = "binomial",
                   alpha = 0 );
    r <- glmnetRank( l_F );
    iGpreselect_iG = randperm( r$coef );
    #     if( is.null( dimnames( dData_iSiG )[[2]] ) ){
    #       iGpreselect_iG = randperm( 1:dim( dData_iSiG )[2] );
    #     }else{
    #       iGpreselect_iG = randperm( dimnames( dData_iSiG )[[2]] );
    #     }
  }else{
    stop( p( "Preselector ", sPresselect, " does not exist." ) );
  }
  
  # Measure accuracy for each feature set
  library( MASS );
  library( e1071 );
  iGpreselect_iG <- iGpreselect_iG[ 1:iNumFeatures ]
  accTr_iG <- array( data = NA,
                     dim = c( length( iGpreselect_iG ) ) );
  accTe_iG <- accTr_iG;
  for( iG in 1:length( iGpreselect_iG ) ){
    
    # Train
    if( sClassifier == "lda" ){
      classif_F <- lda( x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_iG[1:iG] ] ),
                        grouping = bClas_iS[ bTrain_iS ] );
    }else if( sClassifier == "qda" ){
      classif_F <- qda( x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_iG[1:iG] ] ),
                        grouping = bClas_iS[ bTrain_iS ] );
    }else if( sClassifier == "svm" ){
      classif_F <- svm(  x = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_iG[1:iG] ] ),
                         y = as.numeric( bClas_iS[ bTrain_iS ] ) );
    }else{
      stop( p( "The classifier method ", methodClassifier, " is not recognised." ) );
    }
    
    # Test
    preTr_F <- predict( object = classif_F,
                        newdata = as.matrix( dData_iSiG[ bTrain_iS, iGpreselect_iG[1:iG] ] ) );
    preTe_F <- predict( object = classif_F,
                        newdata = as.matrix( dData_iSiG[ bTests_iS, iGpreselect_iG[1:iG] ] ) );
    
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
    accTr_iG[iG] <- mean( preTr2_F == bClas_iS[ bTrain_iS ],
                          na.rm = TRUE );
    accTe_iG[iG] <- mean( preTe2_F == bClas_iS[ bTests_iS ],
                          na.rm = TRUE );
    
  }
  
  # Return result
  r$accTr_iG <- accTr_iG;
  r$accTe_iG <- accTe_iG;
  r$iGpreselect_iG <- iGpreselect_iG;
  return( r );
  
}
