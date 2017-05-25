#
# Doing stuff with UK Biobank
#
# v1 = From scratch
# v2 = Adding functions for boolean matrixes
#


mergeBooleanMatrix = function( matrix_iSiD = bDrug_jViSiD[[1]],
                               oldNames_iD = drugInfo_iDiC[ ,"UKBdrugID"],
                               newNames_iD = drugInfo_iDiC[ ,drugLevel],
                               ignoreCase = TRUE ){
  # Calculates a new boolean matrix "matrix_m_iSiDu", where each column "iDu" represents the conjunction (boolean OR) of all new names
  # For instance, if:
  # oldNames_iD = c( "A", "B", "B", "C", "D" );
  # newNames_iD = c( "a", "b", "c", "d", "d" );
  # The columns of the new matrix will be:
  # a = patient was taking drug A
  # b = patient was taking drug B
  # c = patient was taking drug B
  # d = patient was taking either drug C OR drug D
  #
  
  # Put everything into lower case
  if( ignoreCase ){
    dimnames( matrix_iSiD )[[2]] <- tolower( dimnames( matrix_iSiD )[[2]] );
    oldNames_iD <- tolower( oldNames_iD );
    newNames_iD <- tolower( newNames_iD );
  }
  
  # Check inputs
  oldNames_iD <- as.character( oldNames_iD );
  newNames_iD <- as.character( newNames_iD );
  stopifnot( length( oldNames_iD ) == length( newNames_iD ) );
  stopifnot( sum( !( oldNames_iD %in% dimnames( matrix_iSiD )[[2]] ) ) < 10 );
  
  # Get only available names
  bAvailable_iD <- oldNames_iD %in% dimnames( matrix_iSiD )[[2]];
  oldNames_iD <- oldNames_iD[ bAvailable_iD ];
  newNames_iD <- newNames_iD[ bAvailable_iD ];
  
  # Build new boolean matrix with new names
  newNames_u_iDu <- unique( newNames_iD );
  matrix_m_iSiDu <- array( data = FALSE,
                           dim = c( dim( matrix_iSiD )[1],
                                    length( newNames_u_iDu ) ),
                           dimnames = list( patient = dimnames( matrix_iSiD )[[1]],
                                            drug = newNames_u_iDu ) );
  for( iDu in 1:length( newNames_u_iDu ) ){
    iDs_isDrug <- oldNames_iD[ newNames_iD == newNames_u_iDu[iDu] ];
    matrix_m_iSiDu[ ,iDu] <- rowSums( as.matrix( matrix_iSiD[ ,iDs_isDrug ] ) ) > 0;
  }
  
  # Return result
  return( matrix_m_iSiDu );
}


cutBooleanMatrix = function( matrix_iSiD = bDrug_iSiD,
                             minPower = 0.8,
                             d = 0.25, 
                             sig.level = 0.05 ){
  # Calculates a new boolean matrix "bDrug_c_iSiD", where all the columns with power lower than "minPower" have been eliminated
  
  # Check input
  stopifnot( minPower > 0 );
  stopifnot( minPower < 1 );
  
  # Libraries
  source( "FUNCTIONS_utils.R" )
  library( pwr )
  
  # Calculate power
  count_iD <- colSums( matrix_iSiD );
  pwr_iD <- sapply( count_iD,
                    function( c ){
                      if( c < 2 ){
                        return( 0 );
                      }
                      p <- pwr.t2n.test( n1 = c, 
                                         n2= dim( matrix_iSiD )[1] - c, 
                                         d = d, 
                                         sig.level = sig.level, 
                                         alternative = "two.sided" );
                      return( p$power );
                    } );
  
  # Get only the columns with enough power
  iDs_mostFrequent <- names( pwr_iD )[ pwr_iD >= minPower ];
  print( p( "We have found ", length( iDs_mostFrequent ), " drugs with power > ", 100 * minPower, "%." ) );
  matrix_c_iSiD <- matrix_iSiD[ ,iDs_mostFrequent];
  
  # Return result
  return( matrix_c_iSiD );
}


covariantsBooleanMatrix = function( matrix1_iSiD = bDrug_iSiD,
                                    matrix2_iSiX = bDrug_iSiD,
                                    numCovs = 10,
                                    backupFile = "covariantsBooleanMatrix",
                                    useBackup = TRUE ){
  # Calculates the ORs between all columns of "matrix1_iSiD" with all columns of "matrix2_iSiX".
  # Then, it returns per column in "matrix1_iSiD" the names of the top "numCovs" columns of "matrix2_iSiX" that had highest statistical power for a binomial test
  
  # Calculate fast coincidences matrixes
  backupFile <- p( backupFile,
                   ".Rdata" );
  if( file.exists( backupFile ) & useBackup ){
    print( "Backup file found. Loading councidences matrixes from there." );    
    load( backupFile );
  }else{
    library( pracma );
    print( "Starting to calculate fast matrixes. This is going to take a while." );
    print( "1/4 matrix" );
    drug1dis1_iDiX <- t(  matrix1_iSiD ) %*%   matrix2_iSiX; 
    print( "2/4 matrix" );
    drug1dis0_iDiX <- t(  matrix1_iSiD ) %*%  !matrix2_iSiX; 
    print( "3/4 matrix" );
    drug0dis1_iDiX <- t( !matrix1_iSiD ) %*%   matrix2_iSiX; 
    print( "4/4 matrix" );
    drug0dis0_iDiX <- t( !matrix1_iSiD ) %*%  !matrix2_iSiX; 
    print( "done!" );
    save( file = backupFile,
          list = c( "drug1dis1_iDiX",
                    "drug1dis0_iDiX",
                    "drug0dis1_iDiX",
                    "drug0dis0_iDiX" ) );
  }
  
  # Get odds ratios from fast matrixes
  or_iDiX <- ( drug1dis1_iDiX / drug1dis0_iDiX ) / ( drug0dis1_iDiX / drug0dis0_iDiX );
  or_iDiX[ drug1dis1_iDiX < 100 ] <- NA; #With too few samples the ORs becomes too unstable and noisy
  or_iDiX[ drug1dis0_iDiX < 100 ] <- NA;
  or_iDiX[ drug0dis1_iDiX < 100 ] <- NA;
  or_iDiX[ drug0dis0_iDiX < 100 ] <- NA;
  or_iDiX[ or_iDiX > 100 ] <- 100;
  #stop( "Check that ORs are OK ")
  #image2F(or_iDiX)
  
  
  # Get most coincident events
  covDix_iDiN <- apply( or_iDiX,
                        c( 1 ),
                        function( or_iX ){
                          or_iX <- sort( or_iX,
                                         decreasing = TRUE );
                          return( names( or_iX )[1:numCovs] );
                        } );
  covDix_iDiN <- t( covDix_iDiN );
  
  return( covDix_iDiN );
}

drugsNamesTable = function( drugs_iD ){
  #
  # Gives the table with the IDs, names and sensible names of UK Biobank drugs
  #
  
  # Load the talbe provided by UK Biobank
  library( stringr )
  library( plyr )
  drugInfo_iDiC <- read.csv( "G:/Data/02UKB/Doc/ukb5792x_coding.csv",
                             header = TRUE,
                             stringsAsFactors = FALSE,
                             skip = 3 );
  
  # Produce initial table
  drugInfo_iDiC <- drugInfo_iDiC[ !is.na( drugInfo_iDiC[ ,1] ), ]; #The last row is NA
  dimnames( drugInfo_iDiC )[[1]] <- drugInfo_iDiC[ ,"Code"];
  drugInfo_iDiC[ ,"SensibleName"] <- drugInfo_iDiC[ ,"Meaning"];
  
  # Calculate sensible names
  eliminate_iE <- c( "\\bs\\/f\\b",
                     "\\bm\\/r\\b",
                     "((\\d)+(\\.))?(\\d)+micrograms\\b", # METRIC UNITS---------
                     "((\\d)+(\\.))?(\\d)+microgram\\b",
                     "((\\d)+(\\.))?(\\d)+g(\\/ml|\\/l)?\\b",
                     "((\\d)+(\\.))?(\\d)+mg(\\/ml|\\/l)?\\b",
                     "((\\d)+(\\.))?(\\d)+mgs(\\/ml|\\/l)?\\b",
                     "((\\d)+(\\.))?(\\d)+mcg(\\/ml|\\/l)?\\b",
                     "((\\d)+(\\.))?(\\d)+ml\\b",
                     "((\\d)+(\\.))?(\\d)+ku\\b",
                     "((\\d)+(\\.|\\,))?(\\d)+units(\\/ml|\\/g|\\/l)?\\b",
                     "((\\d)+(\\.|\\,))?(\\d)+u(\\/ml|\\/g|\\/l)?\\b",
                     "\\bmicrograms\\b",
                     "\\bmicrogram\\b",
                     "\\bnanograms\\b",
                     "\\bgm\\b",
                     "\\bmgs\\b",
                     "\\bml\\b",
                     "\\sl\\b",
                     "\\bunits\\b",
                     "\\bunit\\b",
                     "\\bsachet\\b",
                     "\\bmr\\b",
                     "((\\d)+(\\.))?(\\d)%",
                     "(\\d)+(\\.|\\,)(\\d)+",
                     "(\\b|\\-)(\\d){2,}\\b",
                     "\\b(\\d)*million\\b",
                     "(\\d){2,}",
                     "\\(", # GRAMMAR --------------------------------------------
                     "\\)",
                     "\\[",
                     "\\]",
                     "\\bin\b",
                     "\\bof\b",
                     "\\bfrom\b",
                     "\\bear\\b", # BODY PARS ------------------------------------
                     "\\bproduct\\b", # NAMES ------------------------------------
                     "\\bi(\\-)?v\\b",
                     "\\binhaler\\b",
                     "\\beffervescent\\b",
                     "\\bsuppository\\b",
                     "\\bdome\\b",
                     "\\binj\\b",
                     "\\bsyrup\\b",
                     "(\\b|\\+|\\-)buffer\\b",
                     "(\\b|\\+|\\-)diluent\\b",
                     "(\\b|\\+|\\-)kit\\b",
                     "(\\b|\\+|\\-)saline\\b",
                     "(\\b|\\+|\\-)diluent\\b",
                     "(\\b|\\+|\\-)organon\\b",
                     "\\bsyrup\\b",
                     "\\bdispersible\\b",
                     "\\badult\\b",
                     "\\bpaediatric\\b",
                     "\\bforte\\b",
                     "\\bsyrup\\b",
                     "\\bpdr\\b",
                     "\\bfor\\b",
                     "\\brecon\\b",
                     "\\bsuspension\\b",
                     "\\bextra\\b",
                     "\\bsoluble\\b",
                     "\\bintravenous\\b",
                     "\\bretard\\b",
                     "\\bshampoo\\b",
                     "\\btablet\\b",
                     "\\bsalts\\b",
                     "\\bspray\\b",
                     "\\bnasal\\b",
                     "\\bcompound\\b",
                     "\\bspacehaler\\b",
                     "\\beye\\b",
                     "\\be/c\\b",
                     "\\bpreparation\\b",
                     "\\baqueous\\b",
                     "\\bextract\\b",
                     "\\bliquid\\b",
                     "\\bevohaler\\b",
                     "\\boral\\b",
                     "\\bgel\\b",
                     "\\bcream\\b",
                     "\\bointment\\b",
                     "\\binfusion\\b",
                     "\\bconcentratet\\b",
                     "\\bpowder\\b",
                     "\\bcapsule\\b",
                     "\\bturbohaler\\b",
                     "\\bmineral\\b",
                     "\\bsupplement\\b",
                     "\\bdrops\\b",
                     "\\bgranules\\b",
                     "\\bchewable\\b",
                     "\\binjection\\b",
                     "\\bneonatal\\b",
                     "\\binjection\\b",
                     "\\bset\\b",
                     "\\bintramuscular\\b",
                     "\\bim\\b",
                     "\\belixir\\b",
                     "\\bpessary\\b",
                     "\\bsoft\\b",
                     "\\biu\\b",
                     "\\bmiu\\b",
                     "\\bcomplex\\b",
                     "\\bpack\\b",
                     "\\bsublingual\\b",
                     "\\bdose\\b",
                     "\\bvial\\b",
                     "\\bsolution\\b",
                     "\\bcombination\\b" );
  for( iE in 1:length( eliminate_iE ) ){
    drugInfo_iDiC[ ,"SensibleName"] <- gsub( eliminate_iE[iE], 
                                             " ",
                                             drugInfo_iDiC[ ,"SensibleName"],
                                             perl = TRUE,
                                             ignore.case = TRUE );
  }
  drugInfo_iDiC[ ,"SensibleName"] <- gsub( "\\/", 
                                           " ",
                                           drugInfo_iDiC[ ,"SensibleName"],
                                           perl = TRUE,
                                           ignore.case = TRUE );
  drugInfo_iDiC[ ,"SensibleName"] <- gsub( "\\s\\w\\s", 
                                           " ",
                                           drugInfo_iDiC[ ,"SensibleName"],
                                           perl = TRUE,
                                           ignore.case = TRUE );
  drugInfo_iDiC[ ,"SensibleName"] <- str_trim( drugInfo_iDiC[ ,"SensibleName"] );
  drugInfo_iDiC[ ,"SensibleName"] <- gsub( "\\s+", 
                                           " ",
                                           drugInfo_iDiC[ ,"SensibleName"],
                                           perl = TRUE,
                                           ignore.case = TRUE );
  
  # Return table
  return( drugInfo_iDiC )
  
}