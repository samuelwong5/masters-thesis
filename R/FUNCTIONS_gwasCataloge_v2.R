#
# Functions for doing stuff with the GWAS cataloge
#
# Responsabilities
# - Calculate pathway levels
#
# v1 = first function copied from EMIF/it2/epiPath_diseasePathLoad_v3.R
# v2 = diseasePathLoad now searchs diseases by EFO code
# v3 = Split disease path load function into two
#

classifyDisease = function( efoID_iI,
                            class_jCiD ){
  #
  # Given a list of EFO codes for diseases in "efoID_iI", finds which of these diseases bears any of the calssifications indicated in "class_jCiD"
  #
  # Args:
  #   efoID_iI: List of EFO codes for diseases, for instance c( "EFO_0003888", "EFO_0004197, "EFO_0004895" )
  #   class_jCiD: The classifications you want to recognise. For instance 
  #   
  #        class_jCiD = list( cancer = c( "cancer",
  #                                       "neoplasm" ),
  #                           nervoud = c( "nervous system disease" ),
  #                           immune = c( "immune system disease" ) );
  #
  
  # Check into
  stopifnot( class( efoID_iI ) == "character" );
  
  # Load EFO ontology
  library( ontoCAT );
  print( "Loading EFO ontology" );
  efo <- getEFO(); #documentation = http://rpackages.ianhowson.com/bioc/ontoCAT/
  
  # Classify all the IDs
  classification_iIiC <- data.frame( );
  for( iI in 1:length( efoID_iI ) ){
    # get term
    term <- NA;
    try(
      term <- getTermById( efo,
                           as.character( efoID_iI[iI] ) ),
      silent = TRUE
    )
    # Explore term if it exists in the ontology
    if( !is.na( term ) ){
      # Find the classes of this  disease
      parents_iP <- getAllTermParents( efo,
                                       term );
      parents_iP <- sapply( parents_iP,
                            getLabel );
      bClass_iC <- sapply( class_jCiD,
                           function( class_iD ){
                             any( !is.na( match( class_iD,
                                                 parents_iP ) ) );
                           } );
      # Store if it belongs to any of the classes
      bDisease <- "disease" %in% parents_iP;
      if( bDisease & any( bClass_iC ) ){
        classification_iIiC <- rbind( classification_iIiC,
                                      cbind( ID = efoID_iI[iI],
                                             class = names( class_jCiD )[ bClass_iC ] ) );
      }
    }
  }
  
  # Return result
  return( classification_iIiC );
  
}

getAllTraits = function( ){
  #
  # Gets all the traits existing in GWAS cataloge, including their EFO id
  # This function needs the tsv file from GWAS cataloge website, whose content loads into "gwas_iRiC"
  #
  # Returns:
  #   diseases_iDiC: A data frame with 2 colums, one (ID) representing the EFO ids, and another (trait) representing the name of that trait
  #
  
  # Load GWAS cataloge
  gwas_iRiC <- read.table( file = "H:/Data/Gwas catalog/gwas_catalog_v1.0.1-downloaded_2015-11-26.tsv",
                           sep = "\t",
                           header = TRUE,
                           comment.char = "",
                           quote = "" );
  
  # Find terms 
  iRs <- 1:dim( gwas_iRiC )[1];
  efoID_0_iT <- as.character( gwas_iRiC[ iRs, "MAPPED_TRAIT_URI" ] );
  efoName_0_iT <- as.character( gwas_iRiC[ iRs, "MAPPED_TRAIT" ] );
  gwasName_0_iT <- as.character( gwas_iRiC[ iRs, "DISEASE.TRAIT" ] );
  
  # Eliminate the complex terms ( we only want pure diagnosis )
  bComplex_iT <- grepl( ",",
                        efoName_0_iT,
                        fixed = TRUE );
  efoID_0_iT <- efoID_0_iT[ !bComplex_iT ];
  efoName_0_iT <- efoName_0_iT[ !bComplex_iT ];
  gwasName_0_iT <- gwasName_0_iT[ !bComplex_iT ];
  stopifnot( length( efoID_0_iT ) == length( efoName_0_iT ) );
  stopifnot( length( efoID_0_iT ) == length( gwasName_0_iT ) );
  
  # Extract ID number
  efoID_0_iT <- sapply( efoID_0_iT,
                        function( s ){
                          ss <- strsplit( s,
                                          "/")[[1]];
                          return( ss[ length( ss ) ] );
                        } );  
  
  # Return results
  diseases_iDiC <- data.frame( ID = efoID_0_iT,
                               efoTrait = efoName_0_iT,
                               gwasTrait = gwasName_0_iT ) ;
  diseases_iDiC <- unique( diseases_iDiC );
  return( diseases_iDiC );
  
}


loadDiseaseGenes = function( ){
  #
  # Loads all the diseases that there exist in GWAS cataloge, groups them into EFO codes, and then gets the genes of each disease.
  # All results stored in "loadDiseaseGenes_v2.Rdata"
  #
  
  # Get diseases from database
  gwas_iRiC <- read.table( file = "H:/Data/Gwas catalog/gwas_catalog_v1.0.1-downloaded_2015-11-26.tsv",
                           sep = "\t",
                           header = TRUE,
                           comment.char = "",
                           quote = "" );
  efo_iR <- sapply( as.character( gwas_iRiC[ ,"MAPPED_TRAIT_URI"] ),
                    function( s_iS ){
                      ss <- strsplit( s_iS,
                                      "/" )[[1]];
                      return( ss[ length( ss ) ] );
                    } );
  stopifnot( length( efo_iR ) == dim( gwas_iRiC )[1] );
  gwas_iRiC$efo <- as.character( efo_iR );
  efo_iI <- unique( efo_iR );
  
  # Get the best genes of each disease
  source( "FUNCTIONS_utils.R" );
  library( pracma );
  disGenes_iI <- vector( mode = "list",
                         length = length( efo_iI ) );
  disPvals_iI <- vector( mode = "list",
                         length = length( efo_iI ) );
  disNumSt_iI <- vector( mode = "list",
                         length = length( efo_iI ) );
  names( disGenes_iI ) <- efo_iI;
  names( disPvals_iI ) <- efo_iI;
  names( disNumSt_iI ) <- efo_iI;
  for( iI in 1:length( efo_iI ) ){
    
    # Hale the user
    if( mod( iI, 100 ) == 1 ){
      print( p( "Loading disease iI ", iI, " of ", length( efo_iI ) ) );
    }
    
    # Get genes
    iRs_dis <- which( efo_iI[iI] == gwas_iRiC[ ,"efo"] );
    genes_s <- paste( gwas_iRiC[iRs_dis,"REPORTED.GENE.S."],
                      sep = ",",
                      collapse = "," );
    
    # Slip if there are no genes
    if( genes_s == "" ){
      disGenes_iI[[iI]] <- NA;
      disPvals_iI[[iI]] <- NA;
      disNumSt_iI[[iI]] <- NA;  
      next;
    }
    
    # Clean
    genes_iG <- strsplit( genes_s,
                          split = "," )[[1]];
    genes_iG <- strTrim( genes_iG );
    iGs_good1 <- !grepl( "*",
                         genes_iG,
                         fixed = TRUE );
    iGs_good2 <- ( !grepl( "Intergenic",
                         genes_iG )
                   & !grepl( "intergenic",
                           genes_iG ) );
    iGs_good3 <- !genes_iG == "";
    iGs_good4 <- genes_iG != "NR";
    iGs_good5 <- sapply( genes_iG,
                         function( x_iX ){
                           length( x_iX ) > 0;
                         } );
    genes_iG <- genes_iG[ iGs_good1 & iGs_good2 & iGs_good3 & iGs_good4 & iGs_good5 ];
    if( length( genes_iG ) > 0 ){
      genes_iG <- unique( genes_iG );
      
      # Get corresponding p-val
      pVal_iG <- array( data = NA,
                        dim = c( length( genes_iG ), 1 ) );
      numStudies_iG <- array( data = NA,
                              dim = c( length( genes_iG ), 1 ) );
      for( iG in 1:length( genes_iG ) ){
        iRs_gene <- grep( genes_iG[iG],
                          gwas_iRiC[iRs_dis,"REPORTED.GENE.S."] );
        pVal_iR <- gwas_iRiC[ iRs_dis[iRs_gene], "P.VALUE" ];
        pVal_iR <- as.numeric( as.character( pVal_iR ) );
        pVal_iG[iG] <- mean( log( pVal_iR ) / log( 10 ), 
                             na.rm = TRUE );
        pVal_iG[iG] <- 10 ^ pVal_iG[iG];
        numStudies_iG[iG] <- length(iRs_gene);
      }
      
      # Save
      disGenes_iI[[iI]] <- genes_iG;
      disPvals_iI[[iI]] <- pVal_iG;
      disNumSt_iI[[iI]] <- numStudies_iG;
    }
    
  }
  
  # Save results
  save( file = "loadDiseaseGenes_v2.Rdata",
        list = c( "disGenes_iI", 
                  "disPvals_iI", 
                  "disNumSt_iI",
                  "efo_iI" ) );
}


diseaseGenes_DGN = function( do_iD ){
  #
  # Given a list of diseases, returns the genes for those diseases in disGeNet
  #
  # Args:
  #   do_iD: The names of the diseases to be used
  #
  # Returns:
  #   genes_jD: genes in Entrex
  #
  
  # Load DGN
  source( "FUNCTIONS_diseaseOntology_v2.R" );
  dgn_iRiC <- read.table( file = "./curated_gene_disease_associations.txt",
                          sep = "\t",
                          header = TRUE,
                          comment.char = "",
                          quote = "" );
  dgn_iRiC[ ,"diseaseId"] <- sapply( dgn_iRiC[ ,"diseaseId"],
                                     function( s ){
                                       s <- gsub( "umls:",
                                                  "UMLS_CUI:",
                                                  s,
                                                  fixed = TRUE );
                                       return( s );
                                     } );
  umlsID_iD <- translateFromDo( id_iI = do_iD,
                                to = "UMLS_CUI" );
  stopifnot( length( umlsID_iD ) == length( do_iD ) );
  
  # Get the genes of each disease
  disGenes_iI <- vector( mode = "list",
                         length = length( umlsID_iD ) );
  names( disGenes_iI ) <- do_iD;
  for( iD in 1:length( umlsID_iD ) ){
    iRs <- which( dgn_iRiC[ ,"diseaseId"] %in% umlsID_iD[[iD]] );
    disGenes_iI[[iD]] <- dgn_iRiC[ iRs, "geneId" ];
    stopifnot( all( !is.na( disGenes_iI[[iD]] ) ) );
    disGenes_iI[[iD]] <- as.character( disGenes_iI[[iD]] );
  }
  stopifnot( all( sapply( disGenes_iI,
                          length ) > 0 ) );
  
  # Return result
  return( disGenes_iI );
  
}

diseaseGenes_GWAS = function( efo_iD ){
  #
  # Given a list of diseases, returns the genes for those diseases
  #
  # Args:
  #   efo_iD: The names of the diseases to be used. An example c( "EFO_0004540", "EFO_0004621", "EFO_0004763" )
  #
  # Returns:
  #   genes_jD: genes in Entrex
  #
  
  # Load GWAS genes from GWAS cataloge
  if( !file.exists( "loadDiseaseGenes_v2.Rdata" ) ){
    loadDiseaseGenes();
  }
  load( file = "loadDiseaseGenes_v2.Rdata" );
  
  # Get only the diseases we want
  bUse_iD <- match( efo_iD,
                    efo_iI );
  stopifnot( all( !is.na( bUse_iD ) ) );
  disGenes_iI <- disGenes_iI[ bUse_iD ];
  efo_iI <- efo_iI[ bUse_iD ];
  
  # Translate stupid genes into Entrez
  library( org.Hs.eg.db )
  for( iI in 1:length( disGenes_iI ) ){
    
    # Skip if no genes
    if( is.na( disGenes_iI[[iI]] ) ){
      next;
    }
    
    # Translate into entrez
    stopifnot( all( !is.na( disGenes_iI[[iI]] ) ) );
    genesEntrz_iG <- mget( disGenes_iI[[iI]], 
                           org.Hs.egSYMBOL2EG,
                           ifnotfound = NA );
    
    # Keep duplicates, but eliminate NAs
    genesEntrz_iG <- unlist( genesEntrz_iG );
    genesEntrz_iG <- genesEntrz_iG[ !is.na( genesEntrz_iG ) ];
    
    # Save
    disGenes_iI[[iI]] <- genesEntrz_iG;
    
  }
  
  # Return result
  return( disGenes_iI );
  
}


diseasePathLoad = function( efo_iD ){
  #
  # Given a list of diseases, finds the load of their GWAS genes on each KEGpathway
  # This function needs data from "diseaseSignature_v2.Rdata"
  #
  # Args:
  #   efo_iD: The names of the diseases to be used. An example c( "EFO_0004540", "EFO_0004621", "EFO_0004763" )
  #
  # Returns:
  #   pathLoad_iDiP: A 2D array, where rows are diseases and colums and paths. The contents are the proportion of GWAS genes of each disease on each pathway
  #
  
  stop( "Obsolete!" );
  
  # Load GWAS genes from GWAS cataloge
  if( !file.exists( "loadDiseaseGenes_v2.Rdata" ) ){
    loadDiseaseGenes();
  }
  load( file = "loadDiseaseGenes_v2.Rdata" );
  
  # Get only the diseases we want
  bUse_iD <- match( efo_iD,
                    efo_iI );
  stopifnot( all( !is.na( bUse_iD ) ) );
  disGenes_iI <- disGenes_iI[ bUse_iD ];
  disPvals_iI <- disPvals_iI[ bUse_iD ];
  disNumSt_iI <- disNumSt_iI[ bUse_iD ];
  efo_iI <- efo_iI[ bUse_iD ];
  
  # Translate stupid genes into Entrez
  library( org.Hs.eg.db )
  disGenes2_iI <- vector( mode = "list", length = length( disGenes_iI ) );
  disPvals2_iI <- vector( mode = "list", length = length( disGenes_iI ) );
  disNumSt2_iI <- vector( mode = "list", length = length( disGenes_iI ) );
  names( disGenes2_iI ) <- names( disGenes_iI );
  names( disPvals2_iI ) <- names( disGenes_iI );
  names( disNumSt2_iI ) <- names( disGenes_iI );
  for( iI in 1:length( disGenes_iI ) ){
    
    # Eliminate stupid NaNs
    if( length( disGenes_iI[[iI]] ) == 0 || all( is.na( disGenes_iI[[iI]]  ) ) ){
      disGenes2_iI[[iI]] <- c( "dummyGen" );
      disPvals2_iI[[iI]] <- c( 1 );
      disNumSt2_iI[[iI]] <- c( 0 );
      next;
    }
    iGs_bad <- is.na( disPvals_iI[[iI]] ) | is.nan( disPvals_iI[[iI]] );
    disGenes2_iI[iI] <- list( disGenes_iI[[iI]][!iGs_bad] );
    disPvals2_iI[iI] <- list( disPvals_iI[[iI]][!iGs_bad] );
    disNumSt2_iI[iI] <- list( disNumSt_iI[[iI]][!iGs_bad] );
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disPvals2_iI[[iI]] ) );
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disNumSt2_iI[[iI]] ) );
    
    # Give names for control
    names( disGenes2_iI[[iI]] ) <- disGenes2_iI[[iI]];
    names( disPvals2_iI[[iI]] ) <- disGenes2_iI[[iI]];
    names( disNumSt2_iI[[iI]] ) <- disGenes2_iI[[iI]];
    
    # Translate into entrez
    stopifnot( names( disGenes_iI )[iI] == names( disGenes2_iI )[iI] );
    stopifnot( all( !is.na( disGenes_iI[[iI]] ) ) );
    genesEntrz_iG <- mget( disGenes2_iI[[iI]], 
                           org.Hs.egSYMBOL2EG,
                           ifnotfound = NA );
    genesEntrz2_iG <- as.vector( genesEntrz_iG, mode = "character" );
    names( genesEntrz2_iG ) <- names( genesEntrz_iG );
    stopifnot( names( genesEntrz2_iG ) == names( disGenes2_iI[[iI]] ) )
    disGenes2_iI[[iI]] <- genesEntrz2_iG;
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disPvals2_iI[[iI]] ) );
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disNumSt2_iI[[iI]] ) );
    stopifnot( all( names( disGenes2_iI[[iI]] ) == names( disPvals2_iI[[iI]] ) ) );
    stopifnot( all( names( disGenes2_iI[[iI]] ) == names( disNumSt2_iI[[iI]] ) ) );
    
    # Eliminate NAs again
    iGs_bad <- disGenes2_iI[[iI]] == "NA";
    disGenes2_iI[iI] <- list( disGenes2_iI[[iI]][!iGs_bad] );
    disPvals2_iI[iI] <- list( disPvals2_iI[[iI]][!iGs_bad] );
    disNumSt2_iI[iI] <- list( disNumSt2_iI[[iI]][!iGs_bad] );
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disPvals2_iI[[iI]] ) );
    stopifnot( length( disGenes2_iI[[iI]] ) == length( disNumSt2_iI[[iI]] ) );
    stopifnot( all( names( disGenes2_iI[[iI]] ) == names( disPvals2_iI[[iI]] ) ) );
    stopifnot( all( names( disGenes2_iI[[iI]] ) == names( disNumSt2_iI[[iI]] ) ) );
    
  }
  
  # Eliminate genes that are not in KEGG
  paths_jI <- mget( unlist( disGenes2_iI ), 
                    org.Hs.egPATH,
                    ifnotfound=NA );
  bHasPath_iI <- is.na( paths_jI );
  disGenes2b_iI <-  disGenes2_iI;
  disPvals2b_iI <-  disPvals2_iI;
  disNumSt2b_iI <-  disNumSt2_iI;
  for( iI in 1:length( disGenes2_iI ) ){
    bUse <- bHasPath_iI[ disGenes2_iI[[iI]] ];
    if( length( bUse ) == 0 ){
      disGenes2b_iI[[iI]] <- NA;
      disPvals2b_iI[[iI]] <- NA;
      disNumSt2b_iI[[iI]] <- NA;
    }else{
      disGenes2b_iI[[iI]] <- disGenes2_iI[[iI]][ bHasPath_iI[ disGenes2_iI[[iI]] ] ];
      disPvals2b_iI[[iI]] <- disPvals2_iI[[iI]][ bHasPath_iI[ disGenes2_iI[[iI]] ] ];
      disNumSt2b_iI[[iI]] <- disNumSt2_iI[[iI]][ bHasPath_iI[ disGenes2_iI[[iI]] ] ]; stopifnot( length( disGenes2_iI[[iI]] ) == length( disPvals2_iI[[iI]] ) );
    }
    stopifnot( length( disGenes2b_iI[[iI]] ) == length( disPvals2b_iI[[iI]] ) );
    stopifnot( length( disGenes2b_iI[[iI]] ) == length( disNumSt2b_iI[[iI]] ) );
    stopifnot( all( names( disGenes2b_iI[[iI]] ) == names( disPvals2b_iI[[iI]] ) ) );
    stopifnot( all( names( disGenes2b_iI[[iI]] ) == names( disNumSt2b_iI[[iI]] ) ) );
  }
  
  # Eliminate diseases without genes
  iGnum = 0;
  iIs_noGenes <- sapply( disGenes2_iI,
                         function( x_iX ){
                           length( x_iX ) <= iGnum;
                         } );
  disGenes3_iI <-  disGenes2_iI[ !iIs_noGenes ];
  disPvals3_iI <-  disPvals2_iI[ !iIs_noGenes ];
  disNumSt3_iI <-  disNumSt2_iI[ !iIs_noGenes ];
  
  # Get only the best genes
  for( iI in 1:length( disGenes3_iI ) ){
    disPvals3_s <- sort( disPvals3_iI[[iI]] );
    goodNames_iI <- names( disPvals3_s );#[ 1:iGnum ];
    disGenes3_iI[[iI]] <- disGenes3_iI[[iI]][ goodNames_iI ];
    disNumSt3_iI[[iI]] <- disNumSt3_iI[[iI]][ goodNames_iI ];
    disPvals3_iI[[iI]] <- disPvals3_iI[[iI]][ goodNames_iI ];
    #stopifnot( length( disGenes3_iI[[iI]] ) == iGnum );
    stopifnot( length( disGenes3_iI[[iI]] ) == length( disNumSt3_iI[[iI]] ) );
    stopifnot( length( disGenes3_iI[[iI]] ) == length( disPvals3_iI[[iI]] ) );
    stopifnot( all( names( disGenes3_iI[[iI]] ) == names( disNumSt3_iI[[iI]] ) ) );
    stopifnot( all( names( disGenes3_iI[[iI]] ) == names( disPvals3_iI[[iI]] ) ) );
  }
  
  # Calculate connectivity
  source( "FUNCTIONS_utils.R" );
  source( "FUNCTIONS_dataPathways_v3.R" );
  allGenes_iG <- unique( unlist( disGenes3_iI ) );
  bConnected_iGiP <- keggGeneConnectivity( allGenes_iG );
  #bConnected_iGiP <- reactomeGeneConnectivity( allGenes_iG );
  
  # Transform into distance matrix
  library( pracma );
  iConnected_iGiP <- NA * bConnected_iGiP;
  iConnected_iGiP[ bConnected_iGiP ] = 1;
  iConnected_iGiP[ !bConnected_iGiP ] = 0;
  connections_iGiG <- iConnected_iGiP %*% t( iConnected_iGiP );
  stopifnot( all( connections_iGiG == t(connections_iGiG) ) );
  #   iConnected_2_iGiP <- iConnected_iGiP %*% ( t( iConnected_iGiP ) %*% iConnected_iGiP );
  #   iConnected_2_iGiP <- iConnected_2_iGiP %*% ( t( iConnected_2_iGiP ) %*% iConnected_2_iGiP );
  #   iConnected_2_iGiP <- iConnected_2_iGiP %*% ( t( iConnected_2_iGiP ) %*% iConnected_2_iGiP );
  #   iConnected_2_iGiP <- iConnected_2_iGiP %*% ( t( iConnected_2_iGiP ) %*% iConnected_2_iGiP );
  #   iConnected_iGiP <- iConnected_2_iGiP;
  
  #   # Normalise
  #   for( iP in 1:dim( iConnected_iGiP )[2] ){
  #     iConnected_iGiP[ ,iP] <- iConnected_iGiP[ ,iP] / sum( iConnected_iGiP[ ,iP] );
  #   }
  
  # Get proportion of shared hits
  pathLoad_iDiP <- array( data = NA,
                          dim = c( length( disGenes3_iI ),
                                   dim( iConnected_iGiP )[2] ),
                          dimnames = list( diseases = names( disGenes3_iI ),
                                           paths = dimnames( iConnected_iGiP )[[2]] ) );
  for( iI in 1:length( disGenes3_iI ) ){
    pathLoad_iDiP[iI, ] <- colSums( as.matrix( iConnected_iGiP[ disGenes3_iI[[iI]], ] ) );
    pathLoad_iDiP[iI, ] <- pathLoad_iDiP[iI, ];# / sum( pathLoad_iDiP[iI, ]  );
    #pathLoad_iDiP[iI, ] <- pathLoad_iDiP[iI, ] / length( disGenes3_iI[[iI]] );
    #pathLoad_iDiP[iI, ] <- pathLoad_iDiP[iI, ] / mean( pathLoad_iDiP[iI, ] );
  }
  
  # Return results
  return( pathLoad_iDiP );
  
}






