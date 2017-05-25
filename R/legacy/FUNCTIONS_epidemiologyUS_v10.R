#
# Calculates epidemiology relationships between diseases
#
# v1 = Copied from .../EMIF/iteration2/epiPath_diseaseEpi_v3.R
# v2 = dont know
# v3 = comorbidities() added for MBI
# v4 = comorbidities() divided into 2 functions
# v5 = underreporting parametrised per disease
# v6 = Also corrects for tail effect
# v7 = Making functions more general
# v9 = Nothing changed for now
# v10 = Correcting oddsRatio()
#
# INPUT
# - names_iD = Names of the diseases, as listed in GWAS cataloge

source( "FUNCTIONS_utils.R" );
library( pracma );
library( parallel );
try( strRep <- strrep,
     silent = TRUE );


translateMultnum = function( drugs_iD ){
  #
  # Translate Multnum IDs into sensible drug names
  #
  require( plyr );
  
  # Load drug names
  # For earlier tha 2006 we should use another one http://www.cdc.gov/nchs/ahcd/ahcd_database.htm
  drugs_iL = readLines( "G:/Data/National Ambulatory Medical Care Survey/ICPSR_28521/DS0001/drugs.txt" );
  drugs <- paste( drugs_iL,
                  collapse = " " );
  re_F <- gregexpr( "[0-9]{4,}",
                    drugs,
                    perl = TRUE );
  re_F[[1]] <- c( re_F[[1]],
                  nchar( drugs ) );
  med_iM <- c();
  for( iM in 1:( length( re_F[[1]] ) - 1 ) ){
    med_iM <- c( med_iM,
                 substr( drugs,
                         re_F[[1]][iM],
                         re_F[[1]][iM+1]-1 ) );
  }
  med_iMiR <- data.frame( code = sapply( med_iM,
                                         function( s ){
                                           s <- substr( s,
                                                        1,
                                                        5 );
                                           return( as.numeric( s ) );
                                         } ),
                          medication = sapply( med_iM,
                                               function( s ){
                                                 substr( s,
                                                         6,
                                                         nchar( s ) );
                                               } ) );
  
  # Find the names of the drugs
  drugs_iD <- mapvalues( drugs_iD,
                         as.character( med_iMiR$code ),
                         as.character( med_iMiR$medication ),
                         warn_missing = FALSE );
  
  # Return result
  return( drugs_iD );
  
}

expandCombinations = function( expansion_iSiE ){
  #
  # Expand all the combinations of values present in the columns of "expansion_iSiE"
  #
  
  # Calculate all combinations of expanded variables
  if( all( 0 == expansion_iSiE ) ){
    # Get everything if there are no expansion variables
    expandedCombiantions_iEEiE <- as.matrix( "everything" );
    bHasCombination_iSiEE <- array( data = TRUE,
                                    dim = c( dim( expansion_iSiE )[1],
                                             1 ) );
  }else{
    # Calculate all combinations of variables
    print( "Calculating the patients that have the different combinations of expansion variables" );
    expansion_jEiV <- apply( expansion_iSiE,
                             c( 2 ),
                             function( e_iS ){
                               unique( e_iS );
                             } );    
    expandedCombiantions_iEEiE <- expand.grid( expansion_jEiV );
    # Calculate which subjects comply with each combination of variables
    bHasCombination_iSiEE <- array( data = TRUE,
                                    dim = c( dim( expansion_iSiE )[1],
                                             dim( expandedCombiantions_iEEiE )[1] ),
                                    dimnames = list( subject = dimnames( expansion_iSiE )[[1]],
                                                     expandedCombination = apply( expandedCombiantions_iEEiE,
                                                                                  c( 1 ),
                                                                                  function( e_iE ){
                                                                                    paste( e_iE,
                                                                                           collapse = " ~ " );
                                                                                  } ) ) );
    for( iE in 1:dim( expandedCombiantions_iEEiE )[2] ){
      for( iV in 1:length( expansion_jEiV[[iE]] ) ){
        bIsValue_iEE <- expandedCombiantions_iEEiE[ ,iE] == expansion_jEiV[[iE]][iV];
        bIsValue_iS <- expansion_iSiE[ ,iE] == expansion_jEiV[[iE]][iV];
        bHasCombination_iSiEE[ !bIsValue_iS, bIsValue_iEE ] <- FALSE;
      }
    }
    stopifnot( all( 1 == apply( bHasCombination_iSiEE, 
                                c(1), 
                                sum ) ) );
  }
  
  # Return result
  r <- list( bHasCombination_iSiEE = bHasCombination_iSiEE,
             expandedCombiantions_iEEiE = expandedCombiantions_iEEiE );
  return( r );
  
}

probabilities = function( codes_jD = NA,
                          epi_iSiC,
                          expansion_iSiE = pracma::repmat( 0,
                                                           dim( epi_iSiC )[1],
                                                           1 ) ){
  #
  # Calcualtes the probabilities of each disease present in "codes_jD"
  #
  
  # Check input
  stopifnot( dim( epi_iSiC )[1] == dim( epi_iSiC )[1] );
  
  # If no codes inputed, they will be all existing ones
  if( all( is.na( codes_jD ) ) ){
    codes_iD <- as.character( unique( unlist( epi_iSiC ) ) );
    codes_jD <- as.list( codes_iD );
    names( codes_jD ) <- codes_iD;
  }
  
  # Get expanded variables
  r <- expandCombinations( expansion_iSiE );
  bHasCombination_iSiEE <- r$bHasCombination_iSiEE;
  expandedCombiantions_iEEiE <- r$expandedCombiantions_iEEiE;
  
  
  # Find the cases for each disease
  bCase_iSiD <- findCases( codes_jD = codes_jD,
                           epi_iSiC = epi_iSiC );
  stopifnot( dim( bCase_iSiD )[1] == dim( epi_iSiC )[1] );
  
  # Calculate odds ratio for each combination of variables
  print( "Calculating the probabilities for each combination of expanded variables" );
  probs_iDiE <- data.frame();
  for( iD in 1:length( codes_jD ) ){
    # Calculate numbers of people
    probs_iEE <- apply( bHasCombination_iSiEE,
                        c( 2 ),
                        function( b_iS ){
                          mean( bCase_iSiD[ b_iS, iD ] );
                        } );
    num_iEE <- apply( bHasCombination_iSiEE,
                      c( 2 ),
                      function( b_iS ){
                        sum( bCase_iSiD[ b_iS, iD ] );
                      } );
    numSamples_iEE <- colSums( bHasCombination_iSiEE );
    stopifnot( all( !is.na( numSamples_iEE ) ) );
    stopifnot( length( numSamples_iEE ) == length( probs_iEE ) );
    stopifnot( length( numSamples_iEE ) == length( num_iEE ) );
    # Save
    probs_iDiE <- rbind( probs_iDiE,
                         cbind( disease = names( codes_jD )[iD],
                                probDis = probs_iEE,
                                numDis = num_iEE,
                                expandedCombiantions_iEEiE,
                                numSamples = numSamples_iEE ) );
  }
  
  # Eliminate annoying factors
  probs_iDiE$probDis <- as.numeric( as.character( probs_iDiE$probDis ) );
  probs_iDiE$numDis <- as.numeric( as.character( probs_iDiE$numDis ) );
  
  # Return the results
  return( probs_iDiE );
  
}

oddsRatio = function( codes_jD,
                      epi_iSiC,
                      expansion_iSiE = pracma::repmat( 0,
                                                       dim( epi_iSiC )[1],
                                                       1 ) ){
  #
  # Calculates the odds ratio between the diseases in "codes_jD1" and in "codes_jD2". 
  # These probabilities are calculated separatedly for each combination of variables in "expansion_iSiC"
  #
  
  # Check input
  stopifnot( dim( epi_iSiC )[1] == dim( expansion_iSiE )[1] );
  
  # Get expanded variables
  r <- expandCombinations( expansion_iSiE );
  bHasCombination_iSiEE <- r$bHasCombination_iSiEE;
  expandedCombiantions_iEEiE <- r$expandedCombiantions_iEEiE;
  
  
  # Find the cases for each disease
  bCase_iSiD <- findCases( codes_jD = codes_jD,
                           epi_iSiC = epi_iSiC );
  stopifnot( dim( bCase_iSiD )[1] == dim( epi_iSiC )[1] );
  
  # Calculate odds ratio for each combination of variables
  print( "Calculating the OR for each combination of expanded variables" );
  or_iDDiE <- data.frame();
  for( iD1 in 1:length( codes_jD ) ){
    print( p( "Calculating the OR with disease iD1 ",
              iD1, 
              " of ",
              length( codes_jD ) ) );
    for( iD2 in iD1:length( codes_jD ) ){
      if( iD1 != iD2 ){
        
        # Calculate numbers of people
        yesD1yesD2_iS <- bCase_iSiD[ , iD1 ] & bCase_iSiD[ , iD2 ];
        yesD1notD2_iS <- bCase_iSiD[ , iD1 ] & !bCase_iSiD[ , iD2 ];
        notD1yesD2_iS <- !bCase_iSiD[ , iD1 ] & bCase_iSiD[ , iD2 ];
        notD1notD2_iS <- !bCase_iSiD[ , iD1 ] & !bCase_iSiD[ , iD2 ];
        yesD1yesD2_iEE <- apply( bHasCombination_iSiEE,
                                 c( 2 ),
                                 function( b_iS ){
                                   sum( yesD1yesD2_iS[b_iS] );
                                 } );
        yesD1notD2_iEE <- apply( bHasCombination_iSiEE,
                                 c( 2 ),
                                 function( b_iS ){
                                   sum( yesD1notD2_iS[b_iS] );
                                 } );
        notD1yesD2_iEE <- apply( bHasCombination_iSiEE,
                                 c( 2 ),
                                 function( b_iS ){
                                   sum( notD1yesD2_iS[b_iS] );
                                 } );
        notD1notD2_iEE <- apply( bHasCombination_iSiEE,
                                 c( 2 ),
                                 function( b_iS ){
                                   sum( notD1notD2_iS[b_iS] );
                                 } );
        
        # Calculate OR
        stopifnot( length( yesD1yesD2_iEE ) == length( yesD1notD2_iEE ) );
        stopifnot( length( yesD1yesD2_iEE ) == length( notD1yesD2_iEE ) );
        stopifnot( length( yesD1yesD2_iEE ) == length( notD1notD2_iEE ) );
        or_iEE <- ( ( yesD1yesD2_iEE / yesD1notD2_iEE ) 
                    / ( notD1yesD2_iEE / notD1notD2_iEE ) );
        
        # Calculate numbers of people
        numD1_iEE <- apply( bHasCombination_iSiEE,
                            c( 2 ),
                            function( b_iS ){
                              sum( bCase_iSiD[ b_iS, iD1 ] );
                            } );
        numD2_iEE <- apply( bHasCombination_iSiEE,
                            c( 2 ),
                            function( b_iS ){
                              sum( bCase_iSiD[ b_iS, iD2 ] );
                            } );
        numSamples_iEE <- colSums( bHasCombination_iSiEE );
        stopifnot( all( !is.na( numSamples_iEE ) ) );
        stopifnot( length( numSamples_iEE ) == length( numD1_iEE ) );
        stopifnot( length( numSamples_iEE ) == length( numD2_iEE ) );
        stopifnot( length( numSamples_iEE ) == length( or_iEE ) );
        stopifnot( sum( numSamples_iEE ) == dim( bHasCombination_iSiEE )[1] );
        stopifnot( all( numD1_iEE == ( yesD1yesD2_iEE + yesD1notD2_iEE ) ) );
        stopifnot( all( numD2_iEE == ( notD1yesD2_iEE + yesD1yesD2_iEE ) ) );
        
        # Save
        or_iDDiE <- rbind( or_iDDiE,
                           cbind( disease1 = names( codes_jD )[iD1],
                                  disease2 = names( codes_jD )[iD2],
                                  yesD1yesD2 = yesD1yesD2_iEE,
                                  yesD1notD2 = yesD1notD2_iEE,
                                  notD1yesD2 = notD1yesD2_iEE,
                                  notD1notD2 = notD1notD2_iEE,
                                  numDis1 = numD1_iEE,
                                  numDis2 = numD2_iEE,
                                  or_iEE,
                                  expandedCombiantions_iEEiE,
                                  numSamples = numSamples_iEE ) )
      }
    }
  }
  
  # Return the results
  return( or_iDDiE );
  
}

loadNationalAmbulatoryMedicalCareSurvey <- function( ){
  
  # Get all the subdirectories
  mainDir <- "G:/Data/National Ambulatory Medical Care Survey";
  dirs_iD <- list.files( mainDir );
  bIsZip_iD <- sapply( dirs_iD,
                       function( s ){
                         stopifnot( "ICPSR_" == substr( s,
                                                        1,
                                                        6 ) )
                         return( substr( s,
                                         nchar( s ) - 3,
                                         nchar( s ) ) == ".zip" );
                       } );
  dirs_iD <- dirs_iD[ !bIsZip_iD ]
  
  # Load files directory by directory
  epi_iRiC <- data.frame();
  for( iD in 1:length( dirs_iD ) ){
    
    # Hale the user
    print( p( "Loading dir ",
              iD,
              " of ",
              length( dirs_iD ) ) );
    
    # Load data for this directory
    studyID <- strsplit( dirs_iD[iD],
                         "_" )[[1]][2];
    epi0_iRiC <- read.delim( p( mainDir, 
                                "/",
                                dirs_iD[iD],
                                "/DS0001/",
                                studyID,
                                "-0001-Data.tsv" ),
                             sep = "\t",
                             header = TRUE,
                             quote = "" );
    
    # Get the interesting part, replacing unexisting columns with NA dummies
    columns_iC <- c( "VYEAR",
                     "VMONTH",
                     "AGE",
                     "SEX",
                     "RACE",
                     "USETOBAC",
                     sapply( 1:3,
                             function( n ){
                               return( p( "DIAG", 
                                          n ) );
                             } ),
                     sapply( 1:8,
                             function( n ){
                               return( p( "MED", 
                                          n ) );
                             } ) );
    dimnames( epi0_iRiC )[[2]] <- sapply( dimnames( epi0_iRiC )[[2]],
                                          toupper );
    bExistingColums_iC <- columns_iC %in% dimnames( epi0_iRiC )[[2]];
    stopifnot( sum( bExistingColums_iC ) > 10 );
    columns2_iC <- columns_iC;
    columns2_iC[ !bExistingColums_iC ] <- "naDummy";
    epi0_iRiC$naDummy <- NA;
    epi0_c_iRiC <- epi0_iRiC[ , columns2_iC ];
    dimnames( epi0_c_iRiC )[[2]] <- columns_iC
    
    # Accumulate result
    epi_iRiC <- rbind( epi_iRiC,
                       epi0_c_iRiC );
    
  };
  
  # Return loaded data
  return( epi_iRiC );
  
}


prevalenceFullTail <- function( bCase_iSiD,
                                epi_iSiC,
                                codesA_jD ){
  #
  # Calculates the prevalences of the diseases listed in "codesA_jD" when the list of diseases are full
  #
  
  # Calculate probabilities of having full tail
  fileName <- p( "comorbidities_v8_numDiseases_iSnum",
                 dim( bCase_iSiD )[1],
                 "_iDnum",
                 dim( bCase_iSiD )[2],
                 ".Rdata" );
  if( !file.exists( fileName ) ){
    print( "File with statistical curves does not exist. Calculating it will take a bit of time. " );
    numDiseases_iS <- apply( epi_iSiC[ , c("DIAGNOS1", "DIAGNOS2", "DIAGNOS3", "DIAGNOS4", "DIAGNOS5", "DIAGNOS6", "DIAGNOS7" ) ],
                             c( 1 ),
                             function( x_iX ){
                               sum( x_iX != " " );
                             } );
    tailProb_iDiT <- array( data = NA,
                            dim = c( dim( bCase_iSiD )[2],
                                     7 ),
                            dimnames = list( disease = dimnames( bCase_iSiD )[[2]],
                                             tail = c("DIAGNOS1", "DIAGNOS2", "DIAGNOS3", "DIAGNOS4", "DIAGNOS5", "DIAGNOS6", "DIAGNOS7" ) ) );
    save( file = fileName,
          list = c( "numDiseases_iS",
                    "tailProb_iDiT" ) );
  }else{
    load( fileName );
  }
  
  # Calculate prevalences when the tail is full
  tailFit_iDiT  <- NA * tailProb_iDiT;
  a_iD <- rep( NA, dim( bCase_iSiD )[2] );
  prevWhenTailFull_iD <- rep( NA, dim( bCase_iSiD )[2] );
  for( iD in 1:dim( bCase_iSiD )[2] ){
    # Calculate tail
    for( iT in 1:dim( tailProb_iDiT )[2] ){
      #tailProb_iDiT[iD,iT] <- sum( epi_iSiC[ , dimnames( tailProb_iDiT )[[2]][iT] ] %in% codesA_jD[[iD]] ) / sum( numDiseases_iS >= iT );
      #tailProb_iDiT[iD,iT] <- sum( epi_iSiC[ numDiseases_iS == 7, dimnames( tailProb_iDiT )[[2]][iT] ] %in% codesA_jD[[iD]] ) / sum( numDiseases_iS == 7 );
      tailProb_iDiT[iD,iT] <- mean( epi_iSiC[ numDiseases_iS == 7, dimnames( tailProb_iDiT )[[2]][iT] ] %in% codesA_jD[[iD]] );
    }
    # Fit tail with poisson
    kKs_fit <- 3:7;
    poisson = function( a,
                        kKs = 1:7 ){
      ps <- ( ( a ^ ( kKs ) ) / ( fact( kKs ) ) ) * exp( -a )
      return( ps / sum( ps ) );
    }
    a_iA <- 2 ^ seq( 0.1, 10, 0.1 );
    e_iA <- NA * a_iA;
    s <- sum( tailProb_iDiT[iD,  kKs_fit] );
    for( iA in 1:length( a_iA ) ){
      e_iA[iA] <- sum( ( tailProb_iDiT[iD,  kKs_fit] - poisson( a_iA[iA],  kKs_fit ) * s ) ^ 2 );
    }
    kA_ini <- 0.01 * a_iA[ which.min( e_iA ) ];
    kA <- kA_ini;
    e <- e_iA[ which.min( e_iA ) ];
    a <- a_iA[ which.min( e_iA ) ];
    iI <- 0;
    while( iI < 100 ){
      iI <- iI + 1;
      if( mod( iI, 10 ) == 1 ){
        print( p( " e ", e, "; a ",  a, "; kA ", kA ) );
      }
      eNew <- sum( ( tailProb_iDiT[iD,  kKs_fit] - poisson( a + kA,  kKs_fit ) * s ) ^ 2 );
      eNew_opposite <- sum( ( tailProb_iDiT[iD,  kKs_fit] - poisson( a - kA,  kKs_fit ) * s ) ^ 2 );
      if( eNew_opposite < eNew ){
        kA <- -sign( kA ) * min( kA, kA_ini );
      }
      if( e > min( eNew, eNew_opposite ) ){
        kA <- kA * 1.1;
        a <- a + kA;
        e <- min( eNew, eNew_opposite );
      }else{
        kA <- kA / 1.1;
      }
    }
    tailFit_iDiT[iD, ] <- poisson( a ) * s;
    # Use the fit to calculate the probability of having the disease when the tail is saturated
    prevWhenTailFull_iD[ iD ] <- s * sum( poisson( a,
                                                   8:50 ) );
    a_iD[iD] <- a;
  }
  
  # Return result
  names( prevWhenTailFull_iD ) <- dimnames( bCase_iSiD )[[2]];
  return( prevWhenTailFull_iD );
  
}

comorbidities <- function( bCase0_iSiD,
                           iSs0_control_jDiSiK,
                           reportingRate_iD = rep( 1,
                                                dim(  bCase0_iSiD )[1] ),
                           bTailFull_iS = rep( FALSE, dim(  bCase0_iSiD )[1] ),
                           tailProb_iD = rep( 0, dim(  bCase0_iSiD )[2] ) ){
  #
  # Calculate comorbidities scores
  #
  
  # Repair "iSs_control_jDiSiK" when iKnum is 1 (it has to be a matrix for "comorbidities()" to work, but stupid-R tends to simplify it into a "numeric" type)
  iSs_control_jDiSiK <- iSs0_control_jDiSiK;
  for( iDa in 1:length( iSs0_control_jDiSiK ) ) {
    if( length( dim( iSs0_control_jDiSiK[[iDa]] ) ) == 1 ){
      iSs_control_jDiSiK[[iDa]] <- array( data = iSs0_control_jDiSiK[[iDa]],
                                          dim = c( length( iSs0_control_jDiSiK[[iDa]] ),
                                                   1 ),
                                          dimnames = list( subject = names( iSs0_control_jDiSiK[[iDa]] ),
                                                           iK = c( 1 ) ) );
    }
  }
  
  # Check input is OK
  print( "Checking input is OK" );
  iSs_controls <- unlist( iSs_control_jDiSiK );
  iSs_controls <- iSs_controls[ !is.na( iSs_controls ) ];
  iSs_controls <- iSs_controls[ randi( length( iSs_controls ),
                                       100,
                                       1 ) ];
  stopifnot( all( iSs_controls %in% dimnames(  bCase0_iSiD )[[1]] ) );
  stopifnot( all( ( reportingRate_iD <= 1 ) & ( reportingRate_iD > 0 ) ) );
  
  # Calculate comorbidity table for all pairs of diseases
  com_iDiC <- data.frame( );
  for( iDa in 1:length( iSs_control_jDiSiK ) ){
    # Calculate indexes of controls
    print( p( "Calcualting comorbidities of iDa ", iDa ) );
    bBad_iS <- apply( iSs_control_jDiSiK[[iDa]],
                      c( 1 ),
                      function( x_iX ){
                        any( is.na( x_iX ) );
                      } );
    iSs_Acontrols <- unlist( iSs_control_jDiSiK[[iDa]][!bBad_iS, ] )
    iSs_Acases <- dimnames( iSs_control_jDiSiK[[iDa]] )[[1]][!bBad_iS];
    stopifnot( ( length( iSs_Acontrols ) / length( iSs_Acases ) ) == round( length( iSs_Acontrols ) / length( iSs_Acases ) ) ); 
    # Calcualte epidemiology tats
    for( iDb in 1:length( iSs_control_jDiSiK ) ){
      if( iDa != iDb ){
        
        # Calculate comorbidity table, without using controls
        com_i1iC <- data.frame( nameA = dimnames(  bCase0_iSiD )[[2]][iDa],
                                nameB = dimnames(  bCase0_iSiD )[[2]][iDb],
                                disA = sum(  bCase0_iSiD[ ,iDa] ),
                                disB = sum(  bCase0_iSiD[ ,iDb] ),
                                AB = sum(  bCase0_iSiD[ ,iDa] &  bCase0_iSiD[ ,iDb] ),
                                AnotB = sum(  bCase0_iSiD[ ,iDa] & ! bCase0_iSiD[ ,iDb] ),
                                BnotA = sum( ! bCase0_iSiD[ ,iDa] &  bCase0_iSiD[ ,iDb] ),
                                notAnotB = sum( ! bCase0_iSiD[ ,iDa] & ! bCase0_iSiD[ ,iDb] ) );
        com_i1iC$oddsRatio_raw <- ( ( com_i1iC$AB / com_i1iC$BnotA ) 
                                    / (com_i1iC$AnotB / com_i1iC$notAnotB ) );
        
        # Calculate comorbidity table, using controls
        bHasB_Acontrols <-  bCase0_iSiD[ iSs_Acontrols ,iDb];
        bHasB_Acases <-  bCase0_iSiD[ iSs_Acases ,iDb];
        com_i1iC$B_Acontrols <- sum( bHasB_Acontrols );
        com_i1iC$B_Acases <- sum( bHasB_Acases );
        com_i1iC$notB_Acontrols <- sum( !bHasB_Acontrols );
        com_i1iC$notB_Acases <- sum( !bHasB_Acases );
        com_i1iC$oddsRatio_CC_raw <- ( ( com_i1iC$B_Acases / com_i1iC$B_Acontrols ) 
                                       / (com_i1iC$notB_Acases / com_i1iC$notB_Acontrols ) );
        
        # Correct for tail if required
        if( any( tailProb_iD[ c( iDa, iDb ) ] != 0 ) ){
          
          # Correct the population based OR and A-B counts
          stopifnot( all( names( tailProb_iD ) == dimnames(  bCase0_iSiD )[[2]] ) );
          iSnum <- dim(  bCase0_iSiD )[1];
          com_i1iC$AB <- ( com_i1iC$AB 
                           + sum(  bCase0_iSiD[ ,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb]
                           + sum(  bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
                           + sum( ! bCase0_iSiD[ ,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb] * tailProb_iD[iDa] );
          com_i1iC$AnotB <- ( com_i1iC$AnotB
                              - sum(  bCase0_iSiD[ ,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb]
                              + sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa] );# * ( 1 - tailProb_iD[iDb] ) );
          com_i1iC$BnotA <- ( com_i1iC$BnotA 
                            - sum(  bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
                            + sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDb] );# * ( 1 - tailProb_iD[iDa] ) );
          com_i1iC$notAnotB <- ( com_i1iC$notAnotB 
                                 - sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
                                 - sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDb] );
          com_i1iC$oddsRatio_raw <- ( ( com_i1iC$AB / com_i1iC$BnotA ) 
                                      / (com_i1iC$AnotB / com_i1iC$notAnotB ) );
          
          # Correct the case control OR and A-B counts
          stopifnot( all( names( tailProb_iD ) == names( iSs0_control_jDiSiK[[iDa]] ) ) );
          stopifnot( all( names( tailProb_iD ) == names( iSs0_control_jDiSiK[[iDb]] ) ) );
          iSnum <- dim(  bCase0_iSiD )[1];
          iSs <- unique( c( iSs_Acases,
                            iSs_Acontrols ) );
          com_i1iC$B_Acases <- ( com_i1iC$B_Acases 
                                 + sum(  bCase0_iSiD[iSs_Acases,iDa] & ! bCase0_iSiD[iSs_Acases,iDb] & bTailFull_iS[iSs_Acases] ) * tailProb_iD[iDb]
                                 + sum(  bCase0_iSiD[iSs_Acases,iDb] & ! bCase0_iSiD[iSs_Acases,iDa] & bTailFull_iS[iSs_Acases] ) * tailProb_iD[iDa]
                                 + sum( ! bCase0_iSiD[iSs_Acases,iDa] & ! bCase0_iSiD[iSs_Acases,iDb] & bTailFull_iS[iSs_Acases] ) * tailProb_iD[iDb] * tailProb_iD[iDa] );
          com_i1iC$notB_Acases <- ( com_i1iC$notB_Acases
                                    - sum(  bCase0_iSiD[iSs_Acases,iDa] & ! bCase0_iSiD[iSs_Acases,iDb] & bTailFull_iS[iSs_Acases] ) * tailProb_iD[iDb]
                                    + sum( ! bCase0_iSiD[iSs_Acases,iDb] & ! bCase0_iSiD[iSs_Acases,iDa] & bTailFull_iS[iSs_Acases] ) * tailProb_iD[iDa] );# * ( 1 - tailProb_iD[iDb] ) );
          com_i1iC$B_Acontrols <- ( com_i1iC$B_Acontrols 
                                    - sum(  bCase0_iSiD[iSs_Acontrols,iDb] & ! bCase0_iSiD[iSs_Acontrols,iDa] & bTailFull_iS[iSs_Acontrols] ) * tailProb_iD[iDa]
                                    + sum( ! bCase0_iSiD[iSs_Acontrols,iDb] & ! bCase0_iSiD[iSs_Acontrols,iDa] & bTailFull_iS[iSs_Acontrols] ) * tailProb_iD[iDb] );# * ( 1 - tailProb_iD[iDa] ) );
          com_i1iC$notB_Acontrols <- ( com_i1iC$notB_Acontrols 
                                       - sum( ! bCase0_iSiD[iSs_Acontrols,iDb] & ! bCase0_iSiD[iSs_Acontrols,iDa] & bTailFull_iS[iSs_Acontrols] ) * tailProb_iD[iDa]
                                       - sum( ! bCase0_iSiD[iSs_Acontrols,iDb] & ! bCase0_iSiD[iSs_Acontrols,iDa] & bTailFull_iS[iSs_Acontrols] ) * tailProb_iD[iDb] );
#           com_i1iC$B_Acases <- ( com_i1iC$B_Acases 
#                                  + sum(  bCase0_iSiD[iSs,iDa] & ! bCase0_iSiD[iSs,iDb] & bTailFull_iS[iSs] ) * tailProb_iD[iDb]
#                                  + sum(  bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDa]
#                                  + sum( ! bCase0_iSiD[iSs,iDa] & ! bCase0_iSiD[iSs,iDb] & bTailFull_iS[iSs] ) * tailProb_iD[iDb] * tailProb_iD[iDa] );
#           com_i1iC$notB_Acases <- ( com_i1iC$notB_Acases
#                                     - sum(  bCase0_iSiD[iSs,iDa] & ! bCase0_iSiD[iSs,iDb] & bTailFull_iS[iSs] ) * tailProb_iD[iDb]
#                                     + sum( ! bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDa] );# * ( 1 - tailProb_iD[iDb] ) );
#           com_i1iC$B_Acontrols <- ( com_i1iC$B_Acontrols 
#                                     - sum(  bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDa]
#                                     + sum( ! bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDb] );# * ( 1 - tailProb_iD[iDa] ) );
#           com_i1iC$notB_Acontrols <- ( com_i1iC$notB_Acontrols 
#                                        - sum( ! bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDa]
#                                        - sum( ! bCase0_iSiD[iSs,iDb] & ! bCase0_iSiD[iSs,iDa] & bTailFull_iS[iSs] ) * tailProb_iD[iDb] );
#           com_i1iC$B_Acases <- ( com_i1iC$B_Acases 
#                                  + sum(  bCase0_iSiD[,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb]
#                                  + sum(  bCase0_iSiD[,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
#                                  + sum( ! bCase0_iSiD[,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb] * tailProb_iD[iDa] );
#           com_i1iC$notB_Acases <- ( com_i1iC$notB_Acases
#                                     - sum(  bCase0_iSiD[,iDa] & ! bCase0_iSiD[ ,iDb] & bTailFull_iS ) * tailProb_iD[iDb]
#                                     + sum( ! bCase0_iSiD[,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa] * ( 1 - tailProb_iD[iDb] ) );
#           com_i1iC$B_Acontrols <- ( com_i1iC$B_Acontrols 
#                                     - sum(  bCase0_iSiD[,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
#                                     + sum( ! bCase0_iSiD[,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDb] * ( 1 - tailProb_iD[iDa] ) );
#           com_i1iC$notB_Acontrols <- ( com_i1iC$notB_Acontrols 
#                                        - sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDa]
#                                        - sum( ! bCase0_iSiD[ ,iDb] & ! bCase0_iSiD[ ,iDa] & bTailFull_iS ) * tailProb_iD[iDb] );
          com_i1iC$oddsRatio_CC_raw <- ( ( com_i1iC$B_Acases / com_i1iC$B_Acontrols ) 
                                         / (com_i1iC$notB_Acases / com_i1iC$notB_Acontrols ) );
          
        }

#         # Correct for low reporting rate if required. WARNING
#         if( any( reportingRate_iD[ c( iDa, iDb ) ] != 1 ) ){
#           iSnum <- dim(  bCase0_iSiD )[1];
#           comR_iC <- correctTable( AB = com_i1iC$AB,
#                                    AnotB = com_i1iC$AnotB,
#                                    BnotA = com_i1iC$BnotA,
#                                    notAnotB = com_i1iC$notAnotB,
#                                    iSnum = dim(  bCase0_iSiD )[1],
#                                    reportingRate_iAB = reportingRate_iD[ c( iDa, iDb ) ] );
#           com_i1iC[ ,c( "AB", "AnotB", "BnotA", "notAnotB" )] <- comR_iC[c( "AB", "AnotB", "BnotA", "notAnotB" )];
#           com_i1iC$oddsRatio_raw <- ( ( com_i1iC$AB / com_i1iC$BnotA ) 
#                                       / (com_i1iC$AnotB / com_i1iC$notAnotB ) );
#           comR_iC <- correctTable( AB = com_i1iC$B_Acases,
#                                    AnotB = com_i1iC$notB_Acases,
#                                    BnotA = com_i1iC$B_Acontrols,
#                                    notAnotB = com_i1iC$notB_Acontrols,
#                                    iSnum = length( iSs_Acases ) + length( iSs_Acontrols ),
#                                    reportingRate_iAB = reportingRate_iD[ c( iDa, iDb ) ],
#                                    checkOutput = FALSE );
#           com_i1iC$oddsRatio_CC_raw <- ( ( com_i1iC$B_Acases / com_i1iC$B_Acontrols ) 
#                                          / (com_i1iC$notB_Acases / com_i1iC$notB_Acontrols ) );
#           com_i1iC[ ,c( "B_Acases", "notB_Acases", "B_Acontrols", "notB_Acontrols" )] <- comR_iC[c( "AB", "AnotB", "BnotA", "notAnotB" )];
#           com_i1iC[ ,c( "B_Acases", "notB_Acases", "B_Acontrols", "notB_Acontrols" )] <- sapply( com_i1iC[ ,c( "B_Acases", "notB_Acases", "B_Acontrols", "notB_Acontrols" )],
#                                                                                                  function( x ){
#                                                                                                    if( x < 0 ){
#                                                                                                      x = 0;
#                                                                                                    };
#                                                                                                    return( x );
#                                                                                                  } );
#         }

        # Save
        com_iDiC <- rbind( com_iDiC,
                           com_i1iC );
        
      }
    }
  }
  
  # Calculate statistics
  com_iDiC <- testComorbidities( bCase_iSiD = bCase0_iSiD,
                                 com_iDiC );
  
  # Return result
  return( com_iDiC );
  
}


correctTable <- function( AB,
                          AnotB,
                          BnotA,
                          notAnotB,
                          iSnum,
                          reportingRate_iAB,
                          checkOutput = TRUE ){
  #
  # Correct comorbidities table for low diagnosis reporting rate
  #
  
  # Check input
  stopifnot( AB <= iSnum );
  stopifnot( AnotB <= iSnum );
  stopifnot( BnotA <= iSnum );
  stopifnot( notAnotB <= iSnum );
  stopifnot( AB >= 0 );
  stopifnot( AnotB >= 0 );
  stopifnot( BnotA >= 0 );
  stopifnot( notAnotB >= 0 );
  stopifnot( sum( c( AB, AnotB, BnotA, notAnotB ) ) == iSnum );
  stopifnot( all( ( reportingRate_iAB <= 1 ) & ( reportingRate_iAB > 0 ) ) );
  
  # Transform into probabilities
  AB <- AB / iSnum;
  AnotB <- AnotB / iSnum;
  BnotA <- BnotA / iSnum;
  notAnotB <- notAnotB / iSnum;
  # Find real probabilities from probabilities of reported
  AB <- AB / ( reportingRate_iAB[1] * reportingRate_iAB[2] );
  AnotB <- ( AnotB - AB * ( 1 - reportingRate_iAB[2] ) * reportingRate_iAB[1] ) / reportingRate_iAB[1];
  BnotA <- ( BnotA - AB * ( 1 - reportingRate_iAB[1] ) * reportingRate_iAB[2] ) / reportingRate_iAB[2];
  notAnotB <- notAnotB - AB * ( 1 - reportingRate_iAB[1] ) * ( 1 - reportingRate_iAB[2] ) - AnotB * ( 1 - reportingRate_iAB[1] ) - BnotA * ( 1 - reportingRate_iAB[2] );
  # Transform again into numbers
  AB <- AB * iSnum;
  AnotB <- AnotB * iSnum;
  BnotA <- BnotA * iSnum;
  notAnotB <- notAnotB * iSnum;
  
  # Check output
  if( checkOutput ){
    stopifnot( AB <= iSnum );
    stopifnot( AnotB <= iSnum );
    stopifnot( BnotA <= iSnum );
    stopifnot( notAnotB <= iSnum );
    stopifnot( AB >= 0 );
    stopifnot( AnotB >= 0 );
    stopifnot( BnotA >= 0 );
    stopifnot( notAnotB >= 0 );
    stopifnot( round( sum( c( AB, AnotB, BnotA, notAnotB ) ) ) == iSnum );
  }
  
  # Return
  com_i1iC <- c( AB, AnotB, BnotA, notAnotB );
  names( com_i1iC ) <- c( "AB", "AnotB", "BnotA", "notAnotB" );
  return( com_i1iC );
    
}


testComorbidities <- function( bCase_iSiD,
                               com_iDiC,
                               doCaseControl = TRUE ){
  #
  # Run statistical tests for the given comorbidities
  #
  
  com2_iDiC <- data.frame( );
  for( iD in 1:dim( com_iDiC )[1] ){
    # Pupulation oriented stats
    com_i1iC <- com_iDiC[ iD, ];
    or_F <- epitools::oddsratio( x = round( rbind( c( com_i1iC$AB,  com_i1iC$AnotB ),
                                                   c( com_i1iC$BnotA, com_i1iC$notAnotB ) ) ),
                                 method = "fisher",
                                 conf.level = 0.99 );
    com_i1iC$oddsRatio <- or_F$measure["Exposed2","estimate"];
    com_i1iC$vorMin <- or_F$measure["Exposed2","lower"];
    com_i1iC$orMax <- or_F$measure["Exposed2","upper"];
    com_i1iC$pValue <- or_F$p.value["Exposed2","fisher.exact"];
    com_i1iC$pValueFDR <- p.adjust( com_i1iC$pValue,
                                    method = "fdr",
                                    n = floor( 0.5 * ( dim( bCase_iSiD )[2] ^ 2 ) - 0.5 * dim( bCase_iSiD )[2] ) );
    # Case control oriented stats
    if( doCaseControl ){
      or_F <- epitools::oddsratio( x = round( rbind( c( com_i1iC$B_Acases,  com_i1iC$notB_Acases ),
                                                     c( com_i1iC$B_Acontrols, com_i1iC$notB_Acontrols ) ) ),
                                   method = "fisher",
                                   conf.level = 0.99 );
      com_i1iC$oddsRatio_CC <- or_F$measure["Exposed2","estimate"];
      com_i1iC$vorMin_CC <- or_F$measure["Exposed2","lower"];
      com_i1iC$orMax_CC <- or_F$measure["Exposed2","upper"];
      com_i1iC$pValue_CC <- or_F$p.value["Exposed2","fisher.exact"];
      com_i1iC$pValueFDR_CC <- p.adjust( com_i1iC$pValue_CC,
                                         method = "fdr",
                                         n = floor( 0.5 * ( dim( bCase_iSiD )[2] ^ 2 ) - 0.5 * dim( bCase_iSiD )[2] ) );
    }
    # Save progress
    com2_iDiC <- rbind( com2_iDiC,
                        com_i1iC );
  }
  
  # Return results
  return( com2_iDiC );
  
}


translateICD9 = function( codes_jD ){
  #
  # Translates standard ICD9 codes into the US nomenclature
  #
  
  library( "pracma" );
  codes_r_jD <- codes_jD;
  for( iD in 1:length( codes_jD ) ){
    codes_r_jD[[iD]] <- sapply( codes_r_jD[[iD]],
                                function( s ){
                                  if( grepl( ".", s , fixed = TRUE ) == FALSE ){
                                    s <- p( s, "." );
                                  }
                                  dot <- regexec( ".", 
                                                  s, 
                                                  fixed = TRUE )[[1]][1];
                                  while( nchar( s ) - dot < 2 ){
                                    s <- p( s, 
                                            "-" );
                                    dot <- regexec( ".", 
                                                    s, 
                                                    fixed = TRUE )[[1]][1];
                                  }
                                  return( s );
                                } );
    codes_r_jD[[iD]] <- strRep( codes_r_jD[[iD]],
                                ".",
                                "" );
#     stopifnot( all( sapply( codes_r_jD[[iD]],
#                             nchar ) <= 5 ) );
  }
  
  # Reutrn result
  return( codes_r_jD )
}


prevalences = function( codes_jD,
                        epi_iSiC ){
  #
  # Calculates prevalences
  #
  
  library( epitools );
  p_iDi3 <- array( data = NA,
                   dim = c( length( codes_jD ),
                            3 ),
                   dimnames = list( disease = names( codes_jD ),
                                    result = c( "prevalence",
                                                "ciMin",
                                                "ciMax" ) ) );
  for( iD in 1:dim( p_iDi3 )[1] ){
    # Hale the user
    #     if( mod( iD, 100 ) == 1 ){
    print( p( "Calculating stuff for disease iD ", iD, " of ", dim( p_iDi3 )[1] ) );
    #     }
    # Get people
    bHasD_iS <- apply( epi_iSiC,
                       c( 1 ),
                       function( d_iD ){
                         any( d_iD %in% codes_jD[[iD]] );
                       } );
    # Calculate prevalence
    if( sum( bHasD_iS ) > 0 ){
      p_F <- binom.approx( x = as.double( sum( bHasD_iS ) ),
                           n = as.double( length( bHasD_iS ) ) );
      p_iDi3[iD,"prevalence"] <- p_F$proportion;
      p_iDi3[iD,"ciMin"] <- p_F$lower;
      p_iDi3[iD,"ciMax"] <- p_F$upper;
    }
  }
  
  # Return results
  return( p_iDi3 );
  
}


findCases = function( codes_jD,
                      epi_iSiC ){
  #
  # Gets the populations of people with the diseases listed in codes_jD
  #
  
  # Build populations
  print( p( "Building populations for disease list A" ) );
  library( parallel )
  bCase_iSiD <- array( data = NA,
                       dim = c( dim( epi_iSiC )[1], 
                                length( codes_jD ) ),
                       dimnames = list( subject = dimnames( epi_iSiC )[[1]],
                                        disease = names( codes_jD ) ) );
  cl <- makeCluster( getOption( "cl.cores", 
                                4 ) );
  clusterExport( cl = cl,
                 varlist = c( "codes_jD" ),
                 envir = environment() );
  for( iD in 1:length( codes_jD ) ){
    # Select population
    print( p( "Extracting population with ", names( codes_jD )[iD] ) );
    clusterExport( cl = cl,
                   varlist = c( "iD" ),
                   envir = environment() );
    bCase_iS <- parApply( cl = cl,
                          X = epi_iSiC,
                          MARGIN = 1,
                          FUN = function( e_iC ){
                            any( !is.na( match( codes_jD[[iD]],
                                                e_iC ) ) );
                          } );
    bCase_iSiD[ ,iD] <- bCase_iS;
  }
  stopCluster( cl );
  
  # Return results
  return( bCase_iSiD );
  
}
  

findControls = function( bCase_iSiD,
                         epi_iSiC,
                         iKnum = 1,
                         numCores = 4 ){
  #
  # For each patient (dim 1) and each disease (dim 2) represented in "bCase_iSiD", it finds as many controls as indicated in "iKnum"
  #
  
  # Check
  stopifnot( all( dimnames( bCase_iSiD )[[1]] == dimnames( epi_iSiC )[[1]] ) );
  
  # Set up cluster
  stopifnot( dim( bCase_iSiD )[1] == dim( epi_iSiC )[1] );
  iSnum <- dim( bCase_iSiD )[1];
  cl <- makeCluster( getOption( "cl.cores", 
                                numCores ) );
  clusterExport( cl = cl,
                 varlist = c( "epi_iSiC",
                              "iKnum",
                              "mod",
                              "p",
                              "randi",
                              "iKnum",
                              "bCase_iSiD" ),
                 envir = environment() );
  
  # Find the controls for the diseases in A
  print( p( "Finding controls for for disease list A" ) );
  library("pracma");
  iSs_control_jDiSiK <- vector( mode = "list",
                                length = dim( bCase_iSiD )[2] );
  for( iD in 1:dim( bCase_iSiD )[2] ){
    # Initialise
    print( p( "Calculating controls for disease iD ", iD ) );
    iSs_cases <- which( bCase_iSiD[ ,iD] );
    iSs_control_jDiSiK[[iD]] <- array( data = NA,
                                       dim = c( length( iSs_cases ),
                                                iKnum ),
                                       dimnames = list( iS_case = names( iSs_cases ),
                                                        iK_control = 1:iKnum ) );
    iSstart_iC <- 1 + floor( ( length( iSs_cases ) / numCores ) * ( 0:( numCores - 1 ) ) );
    iSend_iC <- 1 + ceil( ( length( iSs_cases ) / numCores ) * ( 1:numCores ) );
    iSend_iC[ length( iSend_iC ) ] <- length( iSs_cases );
    clusterExport( cl = cl,
                   varlist = c( "iSs_cases",
                                "iD",
                                "iSstart_iC",
                                "iSend_iC" ),
                   envir = environment() );
    # Find all codes
    # for( iS in 1:length( iSs_cases ) )
    findControl <- function( iC ){
      iSs <- iSstart_iC[iC]:iSend_iC[iC];
      stopifnot( all( bCase_iSiD[ iSs_cases[iSs] ,iD] ) ); #Check that all these guys are cases
      iSs_control_iSiK <- array( data = NA,
                                 dim = c( length( iSs ),
                                          iKnum ),
                                 dimnames = list( iS_case = names( iSs_cases[iSs] ),
                                                  iK_control = 1:iKnum ) );
      iS0 <- 1;
      for( iS in iSstart_iC[iC]:iSend_iC[iC] ){
        bGoodCases_iS <- ( epi_iSiC[ , "AGE_YRS" ] == epi_iSiC[ iSs_cases[iS], "AGE_YRS" ] ) & 
          ( epi_iSiC[ , "SEX" ] == epi_iSiC[ iSs_cases[iS], "SEX" ] ) &
          ( epi_iSiC[ , "RACE" ] == epi_iSiC[ iSs_cases[iS], "RACE" ] ) &
          ( epi_iSiC[ , "SVYEAR" ] == epi_iSiC[ iSs_cases[iS], "SVYEAR" ] ) &
          ( !bCase_iSiD[ ,iD] );
        iSs_goodCases <- NA;
        if( sum( bGoodCases_iS ) > 0 ){
          iSs_goodCases <- dimnames( epi_iSiC )[[1]][ which( bGoodCases_iS ) ];
          iSs_goodCases <- iSs_goodCases[ randi( length( iSs_goodCases ), iKnum, 1 ) ];
          if( sum( bGoodCases_iS ) <= iKnum ){
            warning( p( "Only ", sum( bGoodCases_iS ), " controls found for iS ", iS, " iD ", iD ) );
          }
        }
        iSs_control_iSiK[iS0, ] <- iSs_goodCases;
        iS0 <- iS0 + 1;
      }
      return( iSs_control_iSiK );
    }
    iSs_control_iSiK_iC <- parSapply( cl = cl,
                                      X = 1:numCores,
                                      FUN = findControl );
    
    # Merge results
    for( iC in 1:numCores ){
      stopifnot( all( ( dimnames( iSs_control_iSiK_iC[[iC]] )[[1]] ) %in% ( dimnames( iSs_control_jDiSiK[[iD]] )[[1]] ) ) );
      iSs_control_jDiSiK[[iD]][ dimnames( iSs_control_iSiK_iC[[iC]] )[[1]], ] <- iSs_control_iSiK_iC[[iC]];
    }
    dimnames( iSs_control_jDiSiK[[iD]] )[[1]] <- names( iSs_cases );
    stopifnot( all( names( iSs_cases ) == dimnames( bCase_iSiD )[[1]][iSs_cases] ) );
    stopifnot( all( names( iSs_cases ) == dimnames( epi_iSiC )[[1]][iSs_cases] ) );
    # Check all is OK
    iSs_allControls <- unlist( iSs_control_jDiSiK[[iD]] );
    iSs_allControls <- iSs_allControls[ !is.na( iSs_allControls ) ];
    stopifnot( all( iSs_allControls %in% ( dimnames( epi_iSiC )[[1]] ) ) );
  }
  stopCluster( cl );
  
  # Check things are OK
  print( p( "Checking stuff" ) );
  for( iD in 1:length( iSs_control_jDiSiK )  ){
    # Check cases are real cases
    iSs <- randi( dim( iSs_control_jDiSiK[[iD]] )[1],
                  10,
                  1 );
    iSs_cases <- dimnames( iSs_control_jDiSiK[[iD]] )[[1]][iSs];
    stopifnot( all( bCase_iSiD[ iSs_cases, iD ] ) );
    # Check subjects have matched variables
    for( iS in iSs ){
      iS_controls <- iSs_control_jDiSiK[[iD]][iS, ];
      iS_case <- dimnames(iSs_control_jDiSiK[[iD]])[[1]][iS];
      if( all( !is.na( iS_controls ) ) & length( iS_controls ) > 0 & all( !is.na( epi_iSiC[ iS_case, ] ) ) ){
        stopifnot( all( epi_iSiC[ iS_controls, "AGE" ] == epi_iSiC[ iS_case, "AGE"] ) );
        stopifnot( all( epi_iSiC[ iS_controls, "SEX" ] == epi_iSiC[ iS_case, "SEX"] ) );
        stopifnot( all( epi_iSiC[ iS_controls, "AGEUNITS" ] == epi_iSiC[ iS_case, "AGEUNITS"] ) );
        stopifnot( all( epi_iSiC[ iS_controls, "RACE" ] == epi_iSiC[ iS_case, "RACE"] ) );
        stopifnot( all( epi_iSiC[ iS_controls, "cohort" ] == epi_iSiC[ iS_case, "cohort"] ) );
      }
    }
  }
  
  # Return results
  return( iSs_control_jDiSiK );
  
}
  

runDWAS = function( iSs_control_jDiSiK,
                    bCase_iSiD,
                    codesB_jB,
                    epi_iSiC,
                    iKnum,
                    conf.level = 1 - ( 0.05 / ( length( codesB_jB ) * length( codesB_jB ) ) ) ){
  #
  # Run a DWAS analysis 
  #
  
  # Run DWAS with disease iD in 2:end for bCase_iSiD)[:,iD]
  print( p( "Calculating odds ratio between lists A and  B" ) );
  library("pwr");
  p_iAiBi4 <- array( data = NA,
                     dim = c( length( iSs_control_jDiSiK ),
                              length( codesB_jB ),
                              6 ),
                     dimnames = list( A = names( iSs_control_jDiSiK ),
                                      B = names( codesB_jB ),
                                      result = c( "oddsRatio_raw",
                                                  "oddsRatio",
                                                  "orMin",
                                                  "orMax",
                                                  "numD_InCases",
                                                  "numD_InControls" ) ) );
  for( iA in 1:dim( p_iAiBi4 )[1] ){
    # Get populations for disease A
    all( dimnames( bCase_iSiD )[[1]][ bCase_iSiD[ ,iA] ] == dimnames( iSs_control_jDiSiK[[iA]] )[[1]] );
    iSs_cases <- which( bCase_iSiD[ , iA ] );
    bIncomplete <- is.na( iSs_cases ) | apply( iSs_control_jDiSiK[[iA]],
                                                      c( 1 ),
                                                      function( x_iX ){
                                                        any( is.na( x_iX ) );
                                                      } );
    iSs_cases <- iSs_cases[ !bIncomplete ];
    iSs_controls <- as.numeric( iSs_control_jDiSiK[[iA]][ !bIncomplete, ] );
    stopifnot( iKnum * length( iSs_cases ) == length( iSs_controls ) );
    stopifnot( all( !is.na( iSs_controls ) ) );
    stopifnot( all( !is.na( iSs_cases ) ) );
    # Calculate odds ratio with disease B
    for( iB in 1:dim( p_iAiBi4 )[2] ){
      # Hale the user
      if( mod( iB, 100 ) == 1 ){
        print( p( "Calculating stuff for disease iA ", iA, " and iB ", iB ) );
      }
      # Get people
      bHasD_cases_iS <- apply( epi_iSiC[iSs_cases, c("DIAGNOS1", "DIAGNOS2", "DIAGNOS3", "DIAGNOS4", "DIAGNOS5", "DIAGNOS6", "DIAGNOS7" ) ],
                               c( 1 ),
                               function( d_iD ){
                                 any( d_iD %in% codesB_jB[[iB]] );
                               } );
      bHasD_controls_iS <- apply( epi_iSiC[iSs_controls, c("DIAGNOS1", "DIAGNOS2", "DIAGNOS3", "DIAGNOS4", "DIAGNOS5", "DIAGNOS6", "DIAGNOS7" ) ],
                                  c( 1 ),
                                  function( d_iD ){
                                    any( d_iD %in% codesB_jB[[iB]] );
                                  } );
      # Run statistics if we have people
      if( ( sum( bHasD_cases_iS ) + sum( bHasD_controls_iS ) ) > 0 ){
        or_F <- epitools::oddsratio( x = rbind( c( sum( bHasD_cases_iS ),  sum( !bHasD_cases_iS ) ),
                                                c( sum( bHasD_controls_iS ), sum( !bHasD_controls_iS ) ) ),
                                     method = "fisher",
                                     conf.level = conf.level );
        p_iAiBi4[iA,iB,"oddsRatio_raw"] <- mean( bHasD_cases_iS ) / mean( bHasD_controls_iS );
        p_iAiBi4[iA,iB,"oddsRatio"] <- or_F$measure["Exposed2","estimate"];
        p_iAiBi4[iA,iB,"orMin"] <- or_F$measure["Exposed2","lower"];
        p_iAiBi4[iA,iB,"orMax"] <- or_F$measure["Exposed2","upper"];
        p_iAiBi4[iA,iB,"numD_InCases"] <- sum( bHasD_cases_iS );
        p_iAiBi4[iA,iB,"numD_InControls"] <- sum( bHasD_controls_iS );
      }
    }
  }
  
  # Return results
  return( p_iAiBi4 );
  
}

