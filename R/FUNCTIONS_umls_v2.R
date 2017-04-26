#
# Functions for doing stuff with the UMLS metathesaurus
#
# Responsabilities
# - Calculate pathway levels
#
# v1 = from scratch
# v2 = "icdCodes" added
#

createMapping = function( ){
  #
  # Loads the important mappings from the UMLS RRF file, and save the results in "createMapping_v1.Rdata"
  #
  
  # Load file with mappings
  print( "About to load the mapping file. This will take a while (200MBs or so)")
  map_iTiC <- read.delim( file = "G:/Project MBI/umls/2015AB/META/MRMAP.RRF",
                          sep = "|",
                          quote = "",
                          header = FALSE,
                          stringsAsFactors = FALSE );
  colNames_iCi2 <- data.frame( colNumber = c( 1,
                                              2,
                                              7,
                                              10,
                                              13,
                                              15,
                                              18, 
                                              24 ),
                               colName = c( "MAPSETCUI", 
                                            "MAPSETSAB",
                                            "FROMID",
                                            "FROMTYPE",
                                            "REL",
                                            "TOID",
                                            "TOTYPE",
                                            "MAPTYPE"),
                               stringsAsFactors = FALSE );
  colnames( map_iTiC )[ colNames_iCi2[ ,"colNumber"] ] <- colNames_iCi2[ ,"colName"];
  
  # Save the mappings of ICD10
  bSNOMEDmappings_iT <- map_iTiC$MAPSETCUI == "C4048021";
  bICD10mappings_iT <- map_iTiC$MAPSETCUI == "C3860204" ;
  map_SNOMED_ICD9_iTiC <- map_iTiC[ bSNOMEDmappings_iT,  ];
  map_ICD10_ICD9_iTiC <- map_iTiC[ bICD10mappings_iT,  ];
  save( file = "createMapping_v1.Rdata",
        list = c( "map_SNOMED_ICD9_iTiC",
                  "map_ICD10_ICD9_iTiC" ) )
  
}

icdCodes = function( umls_iD ){
  
  # This is the damn file I was looking for!
  source( "FUNCTIONS_utils.R" );
  print( "Loading UMLS database. This is going to take a couple of minutes" );
  mrconso_iTiC <- read.delim( file = "G:/Project MBI/umls/2015AB/META/MRCONSO.RRF.aa",
                              sep = "|",
                              quote = "",
                              header = FALSE,
                              stringsAsFactors = FALSE );
  dimnames( mrconso_iTiC )[[2]] <- c( "CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE", "STR", "SRL", "SUPPRESS", "CVF", "wtf"  );
  
  # Keep only the ICD rows
  print( "Cutting down database. This will take ~ 0.5 minutes" );
  bICD_iT = mrconso_iTiC[ ,"SAB"] == "ICD9CM";
  mrconso_iTiC <- mrconso_iTiC[bICD_iT, ];
  
  # Get the ICD codes of each requested UMLS
  library( pracma );
  icds_jD <- vector( mode = "list",
                     length = length( umls_iD ) );
  names( icds_jD ) <- umls_iD
  for( iD in 1:length( umls_iD ) ){
    if( mod( iD, 1000 ) == 1 ){
      print( p( "Extracting ICD9 codes for iD ", iD, " of ", length( umls_iD ) ) );
    }
    bUMLS_iT <- mrconso_iTiC[ ,"CUI"] == umls_iD[iD];
    if( sum( bUMLS_iT ) > 0 ){
      icds_jD[[iD]] <- unique( mrconso_iTiC[bUMLS_iT,"SDUI"] );
    }
  }
  
  # Return results
  return( icds_jD );
  
}



translateSNOME2ICD9 = function( snomed_jDiC ){
  #
  # Using UMLS, translates any list of SNOMED codes intot he equivalent ICD9 codes
  #
  # Args:
  #  snomed_jDiC: List of lists of snomed codes. For instance, the input:
  #
  #   snomed_jDiC <- list( t1d = c( "154673001",
  #                                 "190322003",
  #                                 "190362004",
  #                                 "267469001",
  #                                 "46635009" ),
  #                        ad = c( "154998003",
  #                                "267688001",
  #                                "26929004",
  #                                "73768007" ),
  #                        t2d = c( "54672006",
  #                                 "190323008",
  #                                 "190384004",
  #                                 "267468009",
  #                                 "44054006" ) );
  #
  #   will create the result:
  #
  #   icd9_jDiC <- list( t1d = c( "250.01" ),
  #                      ad = c( "331.0" ),
  #                      t2d = c( "250.00" ) );
  #
  library( pracma );
  source( "FUNCTIONS_utils.R" );
  
  # Load mappings
  load( "createMapping_v1.Rdata" );
  
  # Translate the mappinigs for each disease
  icd9_jDiC <- vector( mode = "list",
                       length = length( snomed_jDiC ) );
  names( icd9_jDiC ) <- names( snomed_jDiC );
  for( iD in 1:length( snomed_jDiC ) ){
    
    # Hale the user
    if( mod( iD, 100 ) == 1 ){
      print( p( "Executing iD ", iD, " of ", length( snomed_jDiC ) ) );
    }
    
    # Find direct mappings between SNOMED and ICD9
    icd9_jC <- sapply( snomed_jDiC[[iD]],
                       function( c ){
                         iTs <- which( c == map_SNOMED_ICD9_iTiC$FROMID );
                         return( map_SNOMED_ICD9_iTiC$TOID[iTs] );
                       } )
    icd9_jDiC[[iD]] <- unlist( icd9_jC );
    
    # Try to enrich the list by mapping back and forth to ICD10
    icd10_jC <- sapply( icd9_jDiC[[iD]],
                        function( c ){
                          iTs <- which( c == map_ICD10_ICD9_iTiC$TOID );
                          return( map_ICD10_ICD9_iTiC$FROMID[iTs] );
                        } );
    icd10_jC <- unlist( icd10_jC );
    icd9_jC <- sapply( icd10_jC,
                       function( c ){
                         iTs <- which( c == map_ICD10_ICD9_iTiC$FROMID );
                         return( map_ICD10_ICD9_iTiC$TOID[iTs] );
                         
                       } );
  }
  
  # Return results
  return( icd9_jDiC );
  
}










