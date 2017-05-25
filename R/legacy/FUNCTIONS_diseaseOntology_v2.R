#
# Functions for doing stuff with the Gene ontology
#
# v1 = from scratch
# v2 = speeding stuff up
#
#For installing further libraries: 
#source("http://bioconductor.org/biocLite.R")
#biocLite("wathever")
#library(wathever)

require(stringr)


diseaseGenes_DO = function( do_iD ){
  #
  # Given a list of diseases, returns the genes associated in DO for each one of those diseases
  #
  
  # Get the OMINs
  omins_jDiG <- translateFromDo( do_iD,
                                 to = "OMIM",
                                 useChildren = TRUE,
                                 useParentUpToLevel = 1 );
  
  # Load OMIM mappings to genes
  omim_iRiC <- read.delim( file = "G:/Data/omim/morbidmap.txt",
                           sep = "|",
                           header  = FALSE );
#   omim_iRiC <- read.delim( file = "G:/Data/omim/geneMap2.txt",
#                            sep = "|",
#                            header  = FALSE );
#   omim_iRiC <- omim_iRiC[ grepl( "gene",
#                                  omim_iRiC[ ,"Type" ],
#                                  fixed = TRUE ), ];
#   dimnames( omim_iRiC )[[1]] <- as.character( omim_iRiC[ ,"X..MIM.Number"] );
  
  # Find the genes of each disease
  entrez_jDiG <- vector( mode = "list",
                       length = length( do_iD ) );
  for( iD in 1:length( do_iD ) ){
    
    # Get OMIM number
    omimnNum_iG <- gsub( "OMIM:",
                         "",
                         omins_jDiG[[iD]],
                         fixed = TRUE );
    omimRegExp <- paste( omimnNum_iG, 
                         collapse = "|" );
  
    # Transform each into a gen
#     bMatch12_iR <- grep( omimRegExp,
#                          omim_iRiC[ ,"V12"] );
#     bMatch11_iR <- grep( omimRegExp,
#                          omim_iRiC[ ,"V11"] );
#     genes_iR <- omim_iRiC[ c( bMatch12_iR,
#                               bMatch11_iR ), "V6" ];
    bMatch1_iR <- grep( omimRegExp,
                         omim_iRiC[ ,"V1"] );
    genes_iR <- omim_iRiC[ bMatch1_iR, "V2" ];
    genes_iR <- as.character( genes_iR );
    genes_iR <- sapply( genes_iR,
                        function( ss ){
                          ss <- gsub( " ",
                                      "",
                                      ss,
                                      fixed = TRUE );
                          ss <- strsplit( ss,
                                          "," );
                          return( ss );
                        } )
    genes_iR <- unlist( genes_iR );
    
    # Transform into entrez
    genes_iR <- mget( as.character( genes_iR ), 
                      org.Hs.egALIAS2EG,
                      ifnotfound = NA );
    entrez_jDiG[[ iD ]] <- unlist( genes_iR );
    
    # Replace if NA if there is nothing
    if( length( entrez_jDiG[[ iD ]] ) == 0 ){
      entrez_jDiG[[ iD ]] <- NA;
    }
    
  }

  names( entrez_jDiG ) <- do_iD;
  
  # Return result
  return( entrez_jDiG );
  
}

findParentsAndChilds = function( fileName ){
  #
  # Introduces into a single list per DO term, all the children and parents of that DO term
  #
  # Args:
  #   fileName: File where the results will be stored
  #
  
  # Build the whole tree of DO if it has not been done in the past
  print( p( "Warning: The hiearchy tree for DO has never been built before. Building will take some time " ) );
  require( pracma );
  require( DO.db );
  require( jsonlite );
  doTerm_iTiC <- dbGetQuery( DO_dbconn(), 
                             "SELECT * FROM do_term" );
  families_iIi3 <- data.frame();
  
  # Function for finding the parents of a DO term
  require( parallel );
  cl <- makeCluster( getOption( "cl.cores", 
                                4 ) );
  clusterExport( cl, 
                 c( "with",
                    "fromJSON",
                    "grepl",
                    "cbind",
                    "p" ) )
  getParents = function( doTerm ){
    # Get matching IDs
    ids_i1i2 <- NA;
    metadata_iM <- NA;
    try( 
      metadata_iM <- fromJSON( p( "http://www.disease-ontology.org/api/metadata/",
                                  doTerm ) ),
      silent = TRUE );
    if( any( !is.na( metadata_iM ) ) & any( "parents" == names( metadata_iM ) ) ){
      ids_i1i2 <-  cbind( childID = doTerm,
                          parentsNames = metadata_iM[["parents"]][ ,2],
                          parentsIDs = metadata_iM[["parents"]][ ,3] );
    }
    return( ids_i1i2 );
  }
  # Find all terms matching any DO
  iTs_iT2 <- seq( 1,
                  dim( doTerm_iTiC )[1],
                  100 );
  if( iTs_iT2[ length( iTs_iT2 ) ] != dim( doTerm_iTiC )[1] ){
    iTs_iT2 <- c( iTs_iT2,
                  dim( doTerm_iTiC )[1] + 1 );
  }
  families_iIi3 <- data.frame();
  for( iT2 in 1:( length( iTs_iT2 ) - 1 ) ){
    # Hale the user
    iTs <- iTs_iT2[iT2]:( iTs_iT2[iT2+1] - 1 );
    print( p( "Finding parents of term iT ", iTs[1] ) );
    # Match ids in parallel
    ids_iTi1i2 <- parSapply( cl = cl,
                             X = doTerm_iTiC[ iTs,
                                              "do_id" ],
                             FUN = getParents );
    # Merge all
    for( iT in 1:length( ids_iTi1i2 ) ){
      if( any( !is.na( ids_iTi1i2[[iT]] ) ) ){
        families_iIi3 <- rbind( families_iIi3,
                                ids_iTi1i2[[iT]] );
      }
    }
    rm( ids_iTi1i2 );
  }
  stopCluster( cl );
  
  # Cluster diseases
  families_iIi3 <- apply( families_iIi3,
                          c( 2 ),
                          as.character )
  disIDs_u_iD <- unique( families_iIi3[ ,"childID"] )
  parentsIDs_u_iD <- vector( mode = "list",
                             length = length( disIDs_u_iD ) );
  parentsLevels_u_iD <- vector( mode = "list",
                                length = length( disIDs_u_iD ) );
  done_iD <- rep( FALSE,
                  length( parentsIDs_u_iD ) );
  iK <- 0;
  while( any( !done_iD ) ){
    iK <- iK + 1;
    for( iD in 1:length( disIDs_u_iD ) ){
      # Hale the user
      if( mod( iD, 1000 ) == 1 ){
        print( p( "Running iteration iK ", iK, ", disease iD ", iD, " of ", length( parentsIDs_u_iD ) ) );
      }
      # Expand childs that still have parents
      if( ! done_iD[iD] ){
        # Get direct parents, if this is the first iteration
        if( is.null( parentsIDs_u_iD[[iD]] ) ){
          iIs <- which( families_iIi3[ ,"childID"] == disIDs_u_iD[iD] );
          parentsIDs_u_iD[[iD]] <- unlist( families_iIi3[iIs,"parentsIDs"] );
          parentsLevels_u_iD[[iD]] <- rep( 1,
                                           length( parentsIDs_u_iD[[iD]] ) );
        }
        # Get parents parents
        stopifnot(  max( parentsLevels_u_iD[[iD]] ) < 100 );
        iPs <- which( parentsLevels_u_iD[[iD]] == max( parentsLevels_u_iD[[iD]] ) );
        iIs <- match( parentsIDs_u_iD[[iD]][iPs],
                      families_iIi3[ ,"childID"] );
        iIs <- iIs[ !is.na( iIs ) ];
        newParents_iI <- families_iIi3[iIs,"parentsIDs"];
        newParents_iI <- newParents_iI[ !is.na( newParents_iI ) ];
        bAlreadyExists_iI <- !is.na( match( newParents_iI,
                                            parentsIDs_u_iD[[iD]] ) );
        if( length( newParents_iI ) == 0 
            || sum( !bAlreadyExists_iI ) == 0 ){
          # There are no more parents to expand
          done_iD[iD] = TRUE;
        }else{
          # Expand the parents that still have parents themselves
          parentsIDs_u_iD[[iD]] <- c( parentsIDs_u_iD[[iD]],
                                      newParents_iI[ !bAlreadyExists_iI ] );
          parentsLevels_u_iD[[iD]] <- c( parentsLevels_u_iD[[iD]],
                                         rep( max( parentsLevels_u_iD[[iD]] ) + 1,
                                              sum( !bAlreadyExists_iI ) ) );
        }
        
      }
    }
  }
  
  # Check all is OK
  stopifnot( length( disIDs_u_iD ) == length( parentsIDs_u_iD ) );
  stopifnot( length( disIDs_u_iD ) == length( parentsLevels_u_iD ) );
  
  # Build matrix with all relationships
  # TO IMPROVE: Not all DO terms are actually in "disIDs_u_iD". I have seen some parents that are not there
  familyLevel_iPiC <- array( data = NA,
                             dim = c( length( disIDs_u_iD ),
                                      length( disIDs_u_iD ) ),
                             dimnames = list( parent = disIDs_u_iD,
                                              child = disIDs_u_iD ) );
  for( iC in 1:dim( familyLevel_iPiC )[2] ){
    iPs <- match(  parentsIDs_u_iD[[iC]], dimnames( familyLevel_iPiC )[[1]] );
    familyLevel_iPiC[ iPs[ !is.na( iPs ) ], iC ] <- parentsLevels_u_iD[[iC]][ !is.na( iPs ) ];
    stopifnot( dimnames( parentsIDs_u_iD )[[2]][iC] == disIDs_u_iD[iC] )
  }
  
  # Now build lists of children
  childrenIDs_u_iD <- vector( mode = "list",
                              length = length( disIDs_u_iD ) );
  childrenLevels_u_iD <- vector( mode = "list",
                                 length = length( disIDs_u_iD ) );
  for( iP in 1:dim( familyLevel_iPiC )[1] ){
    iDs_children <- which( !is.na( familyLevel_iPiC[iP, ] ) );
    childrenIDs_u_iD[[iP]] <- dimnames( familyLevel_iPiC )[[2]][ iDs_children ];
    childrenLevels_u_iD[[iP]] <- familyLevel_iPiC[ iP, iDs_children ];
  }
  
  # Save the results
  names( parentsIDs_u_iD ) <- disIDs_u_iD;
  names( parentsLevels_u_iD ) <- disIDs_u_iD;
  names( childrenIDs_u_iD ) <- disIDs_u_iD;
  names( childrenLevels_u_iD ) <- disIDs_u_iD;
  save( file = fileName,
        list = c( "families_iIi3",
                  "disIDs_u_iD",
                  "parentsIDs_u_iD",
                  "parentsLevels_u_iD",
                  "familyLevel_iPiC",
                  "childrenIDs_u_iD",
                  "childrenLevels_u_iD" ) );
  
  
}

classifyDisease_DO = function( doID_iI,
                               class_jCiD ){
  #
  # Given a list of DO codes for diseases in "doID_iI", finds which of these diseases bears any of the calssifications indicated in "class_jCiD"
  #
  # Args:
  #   doID_iI: List of DO codes for diseases, for instance c( "DOID:8454", "DOID:13909", "DOID:13198" )
  #   class_jCiD: The classifications you want to recognise. For instance 
  #   
  #        class_jCiD = list( cancer = c( "cancer",
  #                                       "neoplasm" ),
  #                           nervoud = c( "nervous system disease" ),
  #                           immune = c( "immune system disease" ) );
  #
  
  # Load the DO tree
  source( "FUNCTIONS_utils.R" );
  fileName <- "classifyDisease_v2.Rdata";
  if( !file.exists( fileName ) ){
    findParentsAndChilds( fileName );
  }
  load( fileName );
  
  # Once we have the DO tree, find which diseases can be classified
  classification_iIiC <- data.frame( stringsAsFactors  = FALSE );
  for( iI in 1:length( doID_iI ) ){
    # See if any parent belongs is any of the categories
    parentsiD_iP <- parentsIDs_u_iD[ doID_iI[iI] ];
    stopifnot( length( parentsiD_iP ) == 1 );
    bClass_iC <- !is.na( match( class_jCiD,
                                parentsiD_iP[[1]] ) );
    # If we found a parent, add it
    if( any( bClass_iC ) ){
      classification_iIiC <- rbind( classification_iIiC,
                                    cbind( ID = doID_iI[iI],
                                           class = names( class_jCiD )[ bClass_iC ] ) );        
    }
  }
  
  # Get rid of the stupid factors
  classification_iIiC[ ,"ID"] <- as.character( classification_iIiC[ ,"ID"] );
  classification_iIiC[ ,"class"] <- as.character( classification_iIiC[ ,"class"] );
  
  # Return result
  return( classification_iIiC );
  
  
}


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
                            to = "ICD9CM",
                            useChildren = FALSE,
                            useParentUpToLevel = 0 ){
  #
  # For each DO id included in "id_iI", returns all the equivalent ids of type "to"
  #
  # Args:
  #   id_iI: The original list of ids
  #   to: the name of the type of ID to translate to. For instance  to = "ICD9CM"
  #   useChildren: If TRUE, the children of each term will also be used in the translation
  #
  # Returns:
  #   do2something_jI: A list of the new ids per ID in "id_iI"
  #
  
  # Get all DO terms if this translation has never been attempted before
  ids_iIi2 <- allDOmappings( to );
  
  # Get also the children of each term, if it was requested
  if( useChildren || useParentUpToLevel > 0 ){
    fileName <- "classifyDisease_v2.Rdata";
    if( !file.exists( fileName ) ){
      findParentsAndChilds( fileName );
    }
    load( fileName );
  }
  
  # Now translate all terms
  do2something_jI <- vector( mode = "list",
                             length = length( id_iI ) );
  names( do2something_jI ) <- id_iI;
  for( iI in 1:length( do2something_jI ) ){
    
    # Get direct mappings
    iIs <- which( id_iI[iI] == ids_iIi2[ ,"DO"] );
    iIs <- iIs[ !is.na( iIs ) ];
    # Add children if requested
    if( useChildren ){
      iIs <- c( iIs,
                which( ids_iIi2[ ,"DO"] %in% unlist( childrenIDs_u_iD[ id_iI[iI] ] ) ) );
    }
    
    # Add 1st parent if requested
    if( length( iIs ) & useParentUpToLevel > 0 ){
      parents_iP <- unlist( parentsIDs_u_iD[ id_iI[iI] ] );
      parents_iP <- parents_iP[ unlist( parentsLevels_u_iD[ id_iI[iI] ] ) <= useParentUpToLevel ];
      iIs_parents <- which( ids_iIi2[ ,"DO"] %in% parents_iP );
      iIs <- c( iIs,
                iIs_parents );      
    }
    
    # Save
    if( length( ids_iIi2[iIs,"with"] ) == 0 ){
      do2something_jI[[iI]] <- NA;
    }else{
      do2something_jI[[iI]] <- c( ids_iIi2[iIs,"with"] );
    }
    
  }
  
  # Chekc
  stopifnot( length( id_iI ) == length( do2something_jI ) );
  
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