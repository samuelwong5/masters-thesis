#
# Functions for doing stuff with the Gene ontology
#
# v1 = from scratch
# v2 = speeding stuff up
#

allDOmappings = function( with = "ICD9CM" ){
  # Finds all the mappings with class "ICD9CM" of all DO terms
  #
  # Args:
  #   with: the name of the IDs to map to
  #
  # Returns:
  #   A data frame with 2 colums
  #   with: The original id
  #   DO: The new DO id
  
  # Get all DO terms
  source( "FUNCTIONS_utils.R" );
  fileName <- p( "allDOmappings_v2_with", 
                 with,
                 ".Rdata" );
  if( file.exists( fileName ) ){
    load( fileName );
  }else{
    # Load all DO terms
    print( p( "Warning: Translationg file not found for ", with, ". Building it for the first time will take some time. " ) );
    require( pracma );
    require( DO.db );
    require( jsonlite );
    doTerm_iTiC <- dbGetQuery( DO_dbconn(), 
                               "SELECT * FROM do_term" );
    ids_iIi2 <- data.frame();
    # Function for finding the matches of a DO term
    cl <- makeCluster( getOption( "cl.cores", 
                                  4 ) );
    clusterExport( cl, 
                   c( "with",
                      "fromJSON",
                      "grepl",
                      "cbind",
                      "p" ) )
    getMatching = function( doTerm ){
      # Get matching IDs
      ids_i1i2 <- NA;
      metadata_iM <- NA;
      try( 
        metadata_iM <- fromJSON( p( "http://www.disease-ontology.org/api/metadata/",
                                    doTerm ) ),
        silent = TRUE );
      if( any( !is.na( metadata_iM ) ) & any( "xrefs" == names( metadata_iM ) ) ){
        bMatch_iX <- grepl( with,
                            metadata_iM[["xrefs"]],
                            fixed = TRUE );
        if( any( bMatch_iX ) ){
          ids_i1i2 <-  cbind( with = metadata_iM[["xrefs"]][ bMatch_iX ],
                              DO = doTerm );
        }
      }
      return( ids_i1i2 );
    }
    # Find all terms matching any DO
    require( parallel );
    iTs_iT2 <- seq( 1,
                    dim( doTerm_iTiC )[1],
                    100 );
    if( iTs_iT2[ length( iTs_iT2 ) ] != dim( doTerm_iTiC )[1] ){
      iTs_iT2 <- c( iTs_iT2,
                    dim( doTerm_iTiC )[1] + 1 );
    }
    ids_iIi2 <- data.frame();
    for( iT2 in 1:( length( iTs_iT2 ) - 1 ) ){
      # Hale the user
      iTs <- iTs_iT2[iT2]:( iTs_iT2[iT2+1] - 1 );
      print( p( "Mapping term iT ", iTs[1] ) );
      # Match ids in parallel
      ids_iTi1i2 <- parSapply( cl = cl,
                               X = doTerm_iTiC[ iTs,
                                                "do_id" ],
                               FUN = getMatching );
      # Merge all
      for( iT in 1:length( ids_iTi1i2 ) ){
        if( any( !is.na( ids_iTi1i2[[iT]] ) ) ){
          ids_iIi2 <- rbind( ids_iIi2,
                             ids_iTi1i2[[iT]] );
        }
      }
      rm( ids_iTi1i2 );
    }
    stopCluster( cl );
    # Now save results
    save( file = fileName,
          list = c( "ids_iIi2" ) );
  }
  
  # Make sure we return characters
  ids_iIi2[ ,1] <- as.character( ids_iIi2[ ,1] );
  ids_iIi2[ ,2] <- as.character( ids_iIi2[ ,2] );
  
  # Return result
  return( ids_iIi2 );
}


translateFromDo = function( id_iI,
                            to = "ICD9CM" ){
  # For each DO id included in "id_iI", returns all the equivalent ids of type "to"
  #
  # Args:
  #   id_iI: The original list of ids
  #   to: the name of the type of ID to translate to
  #
  # Returns:
  #   do2something_jI: A list of the new ids per ID in "id_iI"
  
  # Get all DO terms if this translation has never been attempted before
  ids_iIi2 <- allDOmappings( to );
  
  # Now translate all terms
  do2something_jI <- vector( mode = "list",
                             length = length( id_iI ) );
  names( do2something_jI ) <- id_iI;
  for( iI in 1:length( do2something_jI ) ){
    iIs <-match( id_iI[iI],
                 ids_iIi2[ ,"DO"] );
    iIs <- iIs[ !is.na( iIs ) ];
    if( length( ids_iIi2[iIs,"with"] ) == 0 ){
      do2something_jI[[iI]] <- NA;
    }else{
      do2something_jI[[iI]] <- c( ids_iIi2[iIs,"with"] );
    }
  }
  
  # Return result
  return( do2something_jI );
  
}

translate2do = function( id_iI0,
                         from = "EFO" ){
  # Translates a list of IDs to Disease Ontology ids
  #
  # Args:
  #   id_iI0: The original list of ids
  #   from: the name of the type of ID used in "id_iI0". The name has to correspond to any of the ones listed in http://disease-ontology.org/, in the Xrefs section
  #
  # Returns:
  #   A list with 3 elements
  #   iI0_iD: the indexes of the ids in "id_iI0" that was possible to translate into DO
  #   idFrom_iD: The original id
  #   idDO_iD: The new DO id
  
  # Get all DO terms if this translation has never been attempted before
  ids_iIi2 <- allDOmappings( from );
  
  # Now translate all terms
  iI_iI0 <- match( id_iI0,
                   ids_iIi2[ ,1] );
  r <- list();
  r$iI0_iD <- which( !is.na( iI_iI0 ) );
  r$idFrom_iD <- as.character( ids_iIi2[ iI_iI0[ !is.na( iI_iI0 ) ], 1 ] );
  r$idDO_iD <- as.character( ids_iIi2[ iI_iI0[ !is.na( iI_iI0 ) ], 2 ] );
  stopifnot( length( r$iI0_iD ) == length( r$idFrom_iD ) );
  stopifnot( length( r$iI0_iD ) == length( r$idDO_iD ) );
  
  # Return result
  return( r );
  
}