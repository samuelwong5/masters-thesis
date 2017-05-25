#
# Collection of functions related to pathway analysis
#
# v1 = From scratch
#
#For installing further libraries: 
#source("http://bioconductor.org/biocLite.R")
#biocLite("wathever")
#library(wathever)

#DEPENDENCIES
library (org.Hs.eg.db)
library(KEGGgraph)
library(KEGG.db)
library(reactome.db)
library(data.table)
source( "FUNCTIONS_utils.R" )

loadString_v1 = function( scoreTh = 400 ){
  #
  # Returns the protein-protein connectivity according to STRING. 
  #
  # v1 = For now, it only loads the protein interactions that have direct experimental support
  #
  
  # Return pre-calculated data if it already exists
  fileName <- p( "loadString_v1_scoreTh_", scoreTh, ".Rdata" );
  if( file.exists( fileName ) ){
    print( p( "Preprocessed file ", fileName, " found. Loading..." ) ); 
    load( fileName );
    return( connected_iEiE );
  }
  
  # Load data
  print( "Loading STING database" )
  string_iRiC <- fread( input = "G:/Data/STRING/9606.protein.links.detailed.v10.txt",
                        sep = " ",
                        header = TRUE );
  string_iRiC <- as.data.frame(string_iRiC);
  
  # Get only best known interactions
  string_c_iRiC <- string_iRiC[ string_iRiC[ ,"experimental"] > scoreTh, ]
  if( dim( string_c_iRiC )[1] < 100 ){
    stop( "With threshold ", scoreTh, " we can only find ", dim( string_c_iRiC )[1], " protein-protein interactions.")
  }
  
  # Get the esemble names
  string_c_iRiC[ ,"protein1"] <- gsub( "9606\\.",
                                       "",
                                       string_c_iRiC[ ,"protein1"] )
  string_c_iRiC[ ,"protein2"] <- gsub( "9606\\.",
                                       "",
                                       string_c_iRiC[ ,"protein2"] )
  names_iP <- unique( c( unique(  string_c_iRiC[ ,"protein1"] ), 
                         unique(  string_c_iRiC[ ,"protein2"] ) ) )
  
  # Find the gene coding each protein
  print( "Loading biomart mapping database" )
  library( biomaRt )
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  esProt2esGene_iSpiM <- getBM( attributes = c( "ensembl_gene_id", 
                                                "ensembl_peptide_id" ),
                                filters = "ensembl_peptide_id", 
                                values = names_iP,
                                mart = mart )
  esProts_iP <- unique( esProt2esGene_iSpiM[ ,"ensembl_peptide_id"] )
  esProt2esGene_jSp <- lapply( esProts_iP,
                               function( esProt ){
                                 bProt_iSp <- esProt2esGene_iSpiM[ ,"ensembl_peptide_id"] == esProt
                                 return( esProt2esGene_iSpiM[ bProt_iSp, "ensembl_gene_id" ] )
                               } );
  stopifnot( length( esProt2esGene_jSp ) == length( esProts_iP ) );
  names( esProt2esGene_jSp ) <- esProts_iP; 
  
  # Translate gene names into entrez
  library( plyr )
  esGene2enGene_jSg <- mget( esProt2esGene_iSpiM[ ,"ensembl_gene_id"], 
                             org.Hs.egENSEMBL2EG,
                             ifnotfound = NA );
  enGene_iNg <- unlist( esGene2enGene_jSg );
  
  # Create connectivity matrix
  print( "Transforming protein interactions into a connectivity matrix" )
  library( pracma )
  connected_iEiE <- array( data = FALSE,
                           dim = c( length( enGene_iNg ),
                                    length( enGene_iNg ) ),
                           dimnames = list( proteinA = enGene_iNg,
                                            proteinB = enGene_iNg ) );
  for( iP in 1:length( esProts_iP ) ){
    # Hale the user
    if( mod( iP, 100 ) == 0 ){
      print( p( "Processing protein iP ", iP, " of ", length( esProts_iP ) ) );
    }
    # Find interacting proteins
    main_esProt <- esProts_iP[iP]
    bInteracting_iR <- ( string_c_iRiC[ ,"protein1"] == main_esProt
                         | string_c_iRiC[ ,"protein2"] == main_esProt );
    interacting_esProt_iI <- unique( c( string_c_iRiC[ bInteracting_iR, "protein1" ],
                                        string_c_iRiC[ bInteracting_iR, "protein2" ] ) );
    # Translate proteins to genes
    interacting_esGene_iI <- unlist( esProt2esGene_jSp[ interacting_esProt_iI ] )
    interacting_esGene_iI <- interacting_esGene_iI[ !is.na( interacting_esGene_iI ) ]
    main_esGene <- unlist( esProt2esGene_jSp[ main_esProt ] )
    main_esGene <- main_esGene[ !is.na( main_esGene ) ]
    stopifnot( length( main_esGene ) == 1 )
    # Translate genes to entrez
    interacting_enGene_iI <- unlist( esGene2enGene_jSg[ interacting_esGene_iI ] );
    interacting_enGene_iI <- interacting_enGene_iI[ !is.na( interacting_enGene_iI ) ]
    main_enGene_iM <- unlist( esGene2enGene_jSg[ main_esGene ] );
    main_enGene_iM <- main_enGene_iM[ !is.na( main_enGene_iM ) ]
    if( length( main_enGene_iM ) > 0 ){
      # Mark as connection if the entrez gene exists
      stopifnot( all( interacting_enGene_iI %in% dimnames( connected_iEiE )[[1]] ) )
      stopifnot( all( main_enGene_iM %in% dimnames( connected_iEiE )[[1]] ) )
      connected_iEiE[ main_enGene_iM, interacting_enGene_iI ] <- TRUE;
      connected_iEiE[ interacting_enGene_iI, main_enGene_iM ] <- TRUE;
    }
  }
  
  # Save results
  fileName <- p( "loadString_v1_scoreTh_", scoreTh, ".Rdata" );
  save( file = fileName,
        list = c( "connected_iEiE" ) );
  
  # Return results
  return( connected_iEiE );
  
}





