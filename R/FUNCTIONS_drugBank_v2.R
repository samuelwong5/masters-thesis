#
# Doing stuff with drugbank
#
# v1 = From scratch
# v2 = Adding "find synonyms" capability
#

getDrugsFrame = function(){
  #
  # Gets the table with all the information of DrugBank per drug
  #
  
  # Dependencies
  require( pracma );
  source( "FUNCTIONS_utils.R" );
  
  # Get the translation table
  file <- "getDrugsFrame_v1.Rdata";
  if( !file.exists( file ) ){
    
    # Load data
    warning( "you will need 12GB RAM for this!!!" );
    require( "XML" )
    doc <- xmlTreeParse( "G:/Data/DrugBank/drugbank.xml",
                         getDTD = F );
    r <- xmlRoot( doc );
    
    # Parse the data
    names_jDiN <- vector( mode = "list",
                        length = xmlSize( r ) );
    syns_jDiS <- names_jDiN;
    producst_jDiP <- names_jDiN;
    categories_jDiC <- names_jDiN;
    for( iD in 1:xmlSize( r ) ){
      
      # Hale the user
      if( mod( iD, 100 ) == 1 ){
        print( p( "Parsing node iD ",
                  iD,
                  " of ",
                  xmlSize( r ) ) );
      }
      
      # Get the interesting data
      names_jDiN[[iD]] <- xmlSApply( r[[iD]][[ "name" ]], 
                                     xmlValue );
      syns_jDiS[[iD]] <- xmlSApply( r[[iD]][[ "synonyms" ]], 
                                    xmlValue );
      producst_jDiP[[iD]] <- xmlSApply( r[[iD]][[ "products" ]], 
                                        function( x ){
                                          return( xmlValue( x[["name"]] ) )
                                        } );
      producst_jDiP[[iD]] <- unique( producst_jDiP[[iD]] );
      categories_jDiC[[iD]] <- xmlSApply( r[[iD]][[ "categories" ]], 
                                          xmlValue );
      
    }
    
    # Merge into data frame
    stopifnot( length( names_jDiN ) == length( unlist( names_jDiN ) ) );
    frame_iDiC <- data.frame( row.names = unlist( names_jDiN ) );
    frame_iDiC$synonims <- syns_jDiS;
    frame_iDiC$products <- producst_jDiP;
    frame_iDiC$categories <- categories_jDiC;
    
    # Save the result
    save( file = file,
          list = c( "frame_iDiC" ) );
    
  }else{
    
    load( file );
    
  }
  
  # Return table
  return( frame_iDiC );
  
}

getDrugsTable = function(){
  #
  # Gets the table with all the information of DrugBank per drug
  #
  
  # Dependencies
  require( pracma );
  source( "FUNCTIONS_utils.R" );
  
  # Get the translation table
  file <- "classifyDrug_v1.Rdata";
  if( !file.exists( file ) ){
    
    # Load data
    warning( "you will need 12GB RAM for this!!!" );
    require( "XML" )
    doc <- xmlTreeParse( "G:/Data/DrugBank/drugbank.xml",
                         getDTD = F );
    r <- xmlRoot( doc );
    
    # Parse the data
    table_iDiC <- data.frame();
    for( iN in 1:xmlSize( r ) ){
      
      # Hale the user
      if( mod( iN, 100 ) == 1 ){
        print( p( "Parsing node iN ",
                  iN,
                  " of ",
                  xmlSize( r ) ) );
      }
      
      # Get the interesting data
      names_iN <- xmlSApply( r[[iN]][[ "name" ]], 
                             xmlValue );
      syns_iS <- xmlSApply( r[[iN]][[ "synonyms" ]], 
                            xmlValue );
      categories_iC <- xmlSApply( r[[iN]][[ "categories" ]], 
                                  xmlValue );
      
      # Save into a table
      if( length( categories_iC ) > 0 ){
        for( iC in 1:length( categories_iC ) ){
          table_iDiC <- rbind( table_iDiC,
                               cbind( name = c( names_iN,
                                                syns_iS ),
                                      category = categories_iC[[iC]] ) );
        }
      }
    }
    
    # Save the result
    save( file = file,
          list = c( "table_iDiC" ) );
    
  }else{
    
    load( file );
    
  }
  
  # Return table
  return( table_iDiC );
  
}

findSynonims = function( medication_iM ){
  #
  # Find which drugs among the ones listed in "medication_iM" are synonims or ingredients of each other
  #
  
  # Get the translation table and extract the synonims from there
  frame_iDiC <- getDrugsFrame();
  
  # Extract the synonims from the drug frame
  synmonims_jDiS <- apply( frame_iDiC,
                           c( 1 ),
                           function( f_iC ){
                             return( c( unlist( f_iC[[ "synonims" ]] ),
                                        unlist( f_iC[[ "products" ]] ) ) );
                           } );
  synmonims_jDiS <- lapply( synmonims_jDiS,
                            function( s_iS ){
                              names( s_iS ) <- NULL;
                              s_iS <- unique( s_iS );
                              return( s_iS );
                            } )
  stopifnot( all( names( synmonims_jDiS ) == dimnames( frame_iDiC )[[1]] ) );
  
  # Unify synonims
  syns_iS <- unique( unlist( synmonims_jDiS ) );
  parent_iS <- vector( mode = "character",
                       length = length( syns_iS ) );
  names( parent_iS ) <- syns_iS;
  for( iD in 1:length( synmonims_jDiS ) ){
    parent_iS[ synmonims_jDiS[[iD]] ] <- names( synmonims_jDiS )[iD];
  }
  
  # Add parents as parents of themselves
  parents_iP <- names( synmonims_jDiS );
  names( parents_iP ) <- parents_iP;
  parent_iS <- c( parents_iP,
                  parent_iS );
  
  # Translate the medications given in "medication_iM" into DrugBank names
  library( stringr )
  medciations_t_iM <- medication_iM;
  for( iM in 1:length( medication_iM ) ){
    
    # Hale the user
    if( mod( iM, 100 ) == 1 ){
      print( p( "Translating medication iM ", iM,
                " of ", length( medication_iM ),
                " into a DrugBank name." ) );
    }    
    
    # Find parent synonim for this drug
    medication_iM[iM] <- str_trim( medication_iM[iM] );
    iSs <- grep( p( "\\b", medication_iM[iM], "\\b" ),
                 names( parent_iS ),
                 ignore.case = TRUE );
    #stopifnot( length( unique( parent_iS[iSs] ) ) <= 1 );
    parent <- parent_iS[ iSs[1] ];
    if( length( parent ) == 0 ){
      medciations_t_iM[iM] <- "NA";
    }else{
      medciations_t_iM[iM] <- parent;
    }
    
  }
  
  # Return results
  parentSynonim_iMi2 <- data.frame( medication = medication_iM,
                                    parentSynonim = medciations_t_iM );
  return( parentSynonim_iMi2 );
  
}

#cbind( t(t( medciations_t_iM[1:10])), t(t( medication_iM )) )

classifyDrug = function( medication_iM ){
  #
  # Classify the drugs given in "medication_iM"
  #
  
  # Dependencies
  require( pracma );
  source( "FUNCTIONS_utils.R" );
  
  # Get the translation table
  table_iDiC <- getDrugsTable();
  
  # Translate the drugs
  require( stringr )
  classes_jM <- vector( mode = "list",
                        length = length( medication_iM ) );
  for( iM in 1:length( medication_iM ) ){
    
    # Hale the user
    if( mod( iM, 100 ) == 1 ){
      print( p( "Classifying drug iM ",
                iM,
                " of ",
                length( medication_iM ) ) );
    }    
    
    # Classify this drug
    medication_iM[iM] <- str_trim( medication_iM[iM] );
    iDs <- grep( p( "\\b", medication_iM[iM], "\\b" ),
                 table_iDiC[ ,"name"],
                 ignore.case = TRUE );
    classes_jM[[iM]] <- as.character( unique( table_iDiC[ iDs, "category" ] ) );
    names( classes_jM[[iM]] ) <- NULL;
    
  }
  
  # Return results
  names( classes_jM ) <- medication_iM;
  return( classes_jM );
  
}
  