#
# Collection of functions related to pathway analysis
#
# v3 = Reactome introduced
# v4 = Function "calculatePathHits" added
# v5 = Now "enrichKolg" can also work with reactome
# v6 = Connectivity for GWAS catalogue incorporated
# v7 = getHumanMouseHortologs added
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
source( "FUNCTIONS_utils.R" )

getHumanMouseHortologs( ){
  # Return a ta ble with the orthologs between mouse and humans
  
  library( biomaRt );
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
  hort_iGi2 <- getLDS( attributes = c( "entrezgene" ), 
                       mart = human, 
                       attributesL = c( 
  return( hort_iGi2 )
}

enrichKolg = function( stat_iG,
                       bConnected_iGiP = NA,
                       alternative = "less" ){
  # This function runs enrichment analysis on the Entrez genes listed in "names( stat_iG )", using the values
  # stored in "stat_iG". A boolean matrix describing which genes are part of which pathways can be introduced
  # in "bConnected_iGiP". The return "p_iPiC" will list the p-values and effect sizes corresponding to each pathway.
  # The core test used in "kolmogorov-smirnov"
  
  # Check input is OK
  stopifnot( is.na( bConnected_iGiP ) || ( length( stat_iG ) == dim( bConnected_iGiP )[1] ) )
  
  # Get the connectivity matrix "bConnected_iGiP", which indicates whichg enes are part to which pathways
  if( all( is.na( bConnected_iGiP ) ) ){
#     allEntrez_iB <- names( stat_iG );
#     allEntrez_iB <- allEntrez_iB[ !is.na( allEntrez_iB ) ];
#     pathNames_jGiP <- mget( as.character( allEntrez_iB ), 
#                             org.Hs.egPATH,
#                             ifnotfound = NA );
#     bConnected_iGiP <- calculateConnectivity( allEntrez_iB,
#                                               pathNames_jGiP );
    bConnected_iGiP = keggGeneConnectivity( names( stat_iG ) )
  }
  stopifnot( numel( bConnected_iGiP ) > 10 );
  
  # Test whether any pathway is enriched
  stopifnot( all( names( stat_iG ) == dimnames( bConnected_iGiP )[[1]] ) );
  p_iPiC <- array( data = NA,
                   dim = c( dim( bConnected_iGiP )[2],
                            2 ),
                   dimnames = list( dimnames( bConnected_iGiP )[[2]],
                                    c( "p",
                                       "PathName" ) ) );
  for( iP in 1:dim( bConnected_iGiP )[2] ){
    # Test wether pathway iP is enriched
    ks_F <- ks.test( x = abs( stat_iG[ bConnected_iGiP[ ,iP] ] ),
                     y = abs( stat_iG[ !bConnected_iGiP[ ,iP] ] ),
                     alternative = alternative );
    p_iPiC[iP,"p"] <- ks_F$p.value;
  }
  
  # Get the names of the pathways
  p_iPiC[ ,"PathName"] <- dimnames( bConnected_iGiP )[[2]]
  
  #Return stuff
  return( p_iPiC );
  
}

enrich = function( genes_iG,
                   background_iB,
                   bConnected_iGiP = NA ){
  # This function runs enrichment analysis on the Entrez genes listed in "genes_iG", while considering the background
  # as "background_iB". A boolean matrix describing which genes are part of which pathways can be introduced
  # in "bConnected_iGiP". The return "p_iPiC" will list the p-values and effect sizes corresponding to each pathway.
  # The core test used is "binomial"
  
  # Check input is OK
  stopifnot( all( genes_iG %in% background_iB ) );
  background_iB <- unique( background_iB );
  genes_iG <- unique( genes_iG );
  
  # Get the connectivity matrix "bConnected_iGiP", which indicates whichg enes are part to which pathways
  if( all( is.na( bConnected_iGiP ) ) ){
    allEntrez_iB <- background_iB;
    allEntrez_iB <- allEntrez_iB[ !is.na( allEntrez_iB ) ];
    pathNames_jGiP <- mget( as.character( allEntrez_iB ), 
                            org.Hs.egPATH,
                            ifnotfound = NA );
    bConnected_iGiP <- calculateConnectivity( allEntrez_iB,
                                              pathNames_jGiP );
  }
  
  # Calculate the number of genes per pathway
  stopifnot( all( !is.na( bConnected_iGiP ) ) );
  geneCount_iP <- apply( bConnected_iGiP,
                         c( 2 ),
                         function( b_iG ){
                           sum( 1 * b_iG[ genes_iG ] );
                         } );
  backgroundCount_iP <- apply( bConnected_iGiP,
                         c( 2 ),
                         function( b_iG ){
                           sum( 1 * b_iG[ background_iB ] );
                         } );
  
  # Test whether any pathway is enriched
  stopifnot( length( geneCount_iP ) == length( backgroundCount_iP ) );
  stopifnot( all( !is.na( geneCount_iP ) ) );
  stopifnot( all( !is.na( backgroundCount_iP ) ) );
  p_iPiC <- array( data = NA,
                 dim = c( length( geneCount_iP ),
                          2 ),
                 dimnames = list( dimnames( bConnected_iGiP )[[2]],
                                  c( "p",
                                     "PathName" ) ) );
  for( iP2 in 1:length( geneCount_iP ) ){
    # Test wether pathway iP2 is enriched
    bt_F <- binom.test( x = geneCount_iP[iP2],
                        n = length( background_iB ),
                        p = backgroundCount_iP[iP2] / length( background_iB ),
                        alternative = "greater" );
    p_iPiC[iP2,"p"] <- bt_F$p.value;
  }
  
  # Get the names of the pathways
  library( KEGG.db );
  keggid2keggname <- as.list(KEGGPATHID2NAME);
  p_iPiC[ ,"PathName"] <- unlist( keggid2keggname[ dimnames( p_iPiC )[[1]] ] );
   
  # Return
  return( p_iPiC );
  
}


countPathHits_GE = function( entrez_jDiG,
                             normalise = TRUE ){
  #
  # As in "countPathHits", but it uses gene expression to deduce pathways
  #
  
  # Load gen expression
  print( "countPathHits_GE = Loading LUK dataset" );
  load( "loadAE_luk_security.Rdata" );
  
  # Translate names into entrez
  print( "countPathHits_GE = Translating genes into entrez" );
  geneNames_iG <- dimnames( expression_iGiS )[[1]];
  geneNames_iG <- sapply( geneNames_iG,
                          function( s ){
                            ss <- strsplit( s, 
                                            ".", 
                                            fixed = TRUE );
                            ss <- unlist( ss );
                            return( ss[1] )
                          } );
  geneNames_iG[ geneNames_iG == "" ] <- NA;
  geneNames_iG[ !is.na( geneNames_iG )] <- mget( geneNames_iG[ !is.na( geneNames_iG )], 
                                                 org.Hs.egREFSEQ2EG,
                                                 ifnotfound = NA );
  expression_iGiS <- expression_iGiS[ !is.na( geneNames_iG ), ];
  geneNames_iG <- geneNames_iG[ !is.na( geneNames_iG ) ];
  
  # Select only the genes that are present in the input
  stopifnot( dim( expression_iGiS )[1] == length( geneNames_iG ) );
  geneNames_in_iG <- unique( unlist( entrez_jDiG ) );
  bBoth_iG <- geneNames_iG %in% geneNames_in_iG;
  expression_iGiS <- expression_iGiS[ bBoth_iG, ];
  geneNames_iG <- geneNames_iG[ bBoth_iG ];
  
  # Average the expression of duplicated genes
  print( "countPathHits_GE = Averaging duplicated genes" );
  require( pracma )
  bDuplicateed_iG <- duplicated( geneNames_iG );
  geneNames_dup_iG <- unique( unlist( geneNames_iG[ bDuplicateed_iG ] ) );
  for( iG in 1:length( geneNames_dup_iG ) ){
    iGs <- which( geneNames_iG == geneNames_dup_iG[iG] );
    expression_iGiS[ iGs, ] <- repmat( t( as.matrix( colMeans( expression_iGiS[ iGs, ],
                                                               na.rm = TRUE ) ) ),
                                       n = length( iGs ),
                                       m = 1 );
  }
  expression_iGiS <- expression_iGiS[ !bDuplicateed_iG, ];
  geneNames_iG <- geneNames_iG[ !bDuplicateed_iG ]
  stopifnot( dim( expression_iGiS )[1] == length( unique( geneNames_iG ) ) );
  
  if( !file.exists( "pca.Rdata" ) ){
    # Zscore
    print( "countPathHits_GE = Zscoring data" );
    expression_n_iGiS <- apply( expression_iGiS,
                                c( 1 ),
                                function( x_iX  ){
                                  x_iX <- ( x_iX - mean( x_iX,
                                                         na.rm = TRUE ) ) / sd( x_iX,
                                                                                na.rm = TRUE );
                                  return( x_iX );
                                } );
    expression_n_iGiS <- t( expression_n_iGiS );
    
    #   # Calculate gen gen correlation LONG
    #   cor_iGiG <- cor( x = t( expression_iGiS ),
    #                    use = "everything" );
    #   save( file = "cor_iGiG.Rdata",
    #         list = c( "cor_iGiG" ) );
    
    # Calculate PCA, and get the 100 first components as pathways
    library( fastICA );
    print( "countPathHits_GE = Calculating PCA. This is going to take a while..." );
#     pca_F <- princomp( x = t( expression_n_iGiS ) );
#     save( file = "pca.Rdata",
#           list = c( "pca_F",
#                     "expression_n_iGiS" ) );
    ica_F <- fastICA( X = t( expression_n_iGiS ),
                      n.comp = 20 );
    save( file = "ica.Rdata",
          list = c( "ica_F",
                    "expression_n_iGiS" ) );
    
  }else{
    #load( "pca.Rdata" );
    load( "ica.Rdata" );
  }

pca_F <- ica_F;
pca_F$scores <- ica_F$S;
  
  # Transform PCAs into connectivity matrixes
  numPaths <- 20;
  l_iGiP <- t( cor( pca_F$scores[,1:numPaths],
                    t( expression_n_iGiS ) ) );
  #bConnected2_iGiP <- apply( pca_F$loadings[ ,1:numPaths],
  bConnected_iGiP <- apply( l_iGiP,
                            c( 2 ),
                            function( x_iG ){
                              #x_iG <- rand( length( x_iG ), 1 );
                              abs( x_iG ) >  quantile( abs( x_iG ),
                                                       0.9 );
                            } );
#    bConnected_iGiP <- abs( l_iGiP ) > quantile( as.numeric( l_iGiP ),
#                                                          0.99 );
  
  # Apply the entrez names to the connectivity matrix
  stopifnot( all( sapply( geneNames_iG, 
                          length ) == 1 ) );
  dimnames( bConnected_iGiP )[[1]] <- unlist( geneNames_iG );
  
  # Now calculate the path counts
  print( "countPathHits_GE = Calculating pathway load, using PCs as pathways" );
  count_iDiP <- countPathHits( entrez_jDiG,
                               bConnected_iGiP );
  
  # Normalise counts to the number of genes
  if( normalise ){
    numGenes_iP <- colSums( bConnected_iGiP );
    print( p( "Average number of genes per path ", mean( numGenes_iP ) ) );
    stopifnot( all( !is.na( numGenes_iP ) ) );
    stopifnot( all( numGenes_iP != 0 ) );
    count_iDiP <- apply( count_iDiP,
                         c( 1 ),
                         function( c_iP ){
                           return( c_iP / numGenes_iP );
                         } );
    count_iDiP <- t( count_iDiP );
  }
  
  
  # Return results
  r <- list( count_iDiP = count_iDiP,
             bConnected_iGiP = bConnected_iGiP );
  return( r )
  
}

countPathHits = function( entrez_jDiG,
                          bConnected_iGiP = NULL ){
  #
  # Given a list of groups of genes (e.g. each group of gene corresponds to a disease), this function calculates the number of hits per pathway for that list of genes
  #
  
  # Get connectivity matrix
  if( all( is.na( bConnected_iGiP ) ) ){
    allEntrez_iG <- unique( unlist( entrez_jDiG ) );
    allEntrez_iG <- allEntrez_iG[ !is.na( allEntrez_iG ) ];
    pathNames_jGiP <- mget( as.character( allEntrez_iG ), 
                            org.Hs.egPATH,
                            ifnotfound = NA );
    bConnected_iGiP <- calculateConnectivity( allEntrez_iG,
                                              pathNames_jGiP );
  }
  
  # Count hits per disease in each path
  count_iDiP <- array( data = NA,
                       dim = c( length( entrez_jDiG ),
                                dim( bConnected_iGiP )[2] ),
                       dimnames = list( disease = names( entrez_jDiG ),
                                        path = dimnames( bConnected_iGiP )[[2]] ) ); 
  for( iD in 1:length( entrez_jDiG ) ){
    count_iP <- apply( bConnected_iGiP,
                       c( 2 ),
                       function( connected_iG ){
                         sum( connected_iG[ unique( entrez_jDiG[[iD]] ) ],
                              na.rm = TRUE );
                       } );
    count_iDiP[iD, ] <- count_iP;
  }
  
  # Return results
  return( count_iDiP );
  
}

# === Claculates connectivity between genes and pathways
# CAUTION: Genes not assigned to any pathway in KEGG will be ignored
# geneLevels = Gene expression levels that have to be collapsed into pathway levels
# geneLevels_ref = Refference (i.e. control) gene expression levels
# returning$pathLevels = Empty matrix for the pathway levels
# returning$bConnect_iP2iG2 = Connectivity between pathways and genes
keggGeneConnectivity = function( geneNames_iG ) {
  
  # Get paths
  pathNames_jGiP <- mget(as.character(geneNames_iG), org.Hs.egPATH,ifnotfound=NA);
  
  #Calculate connectivity
  bConnected_iGiP <- calculateConnectivity( geneNames_iG,
                                            pathNames_jGiP)
  
  # Give names
  library( KEGG.db );
  keggid2keggname <- as.list(KEGGPATHID2NAME);
  dimnames( bConnected_iGiP )[[2]] <- unlist( keggid2keggname[ dimnames( bConnected_iGiP )[[2]] ] );
  
  #Return
  return( bConnected_iGiP );
  
}

reactomeGeneConnectivity = function( geneNames_iG ) {
  
  # Get paths
  pathNames_jGiP <- mget(as.character(geneNames_iG), reactomeEXTID2PATHID,ifnotfound=NA);
  
  #Calculate connectivity
  bConnected_iGiP <- calculateConnectivity( geneNames_iG,
                                            pathNames_jGiP)
  
  library( reactome.db );
  keggid2keggname <- as.list( reactomePATHID2NAME );
  dimnames( bConnected_iGiP )[[2]] <- unlist( keggid2keggname[ dimnames( bConnected_iGiP )[[2]] ] );
  
  #Return
  return( bConnected_iGiP );
  
}

dgnGeneConnectivity = function( entrez_iG ) {
  
  # Get diseases from database
  dgn_iRiC <- read.table( file = "G:/Data/Disgenet/all_gene_disease_associations.tsv",
                          sep = "\t",
                          header = TRUE,
                          comment.char = "#",
                          quote = "" );
  
  # Cut DGN into the genes that exist in the input list "entrez_iG"
  bExist_iR <- dgn_iRiC[ ,"geneId"] %in% entrez_iG
  dgn_c_iRiC <- dgn_iRiC[bExist_iR, ];
  
  # Get only the strong assocaitions
  source( "FUNCTIONS_utils.R" )
  sTh <- 0.1;
  bStrong_iR <- dgn_c_iRiC[ ,"score"] > sTh
  dgn_cc_iRiC <- dgn_c_iRiC[ bStrong_iR, ]
  diseases_iD <- as.character( unique( dgn_cc_iRiC[ ,"diseaseName"] ) )
  print( p( "At score threshold ", sTh, " we find ", length( diseases_iD ), " DGN diseases" ) )
  
  # Calculate the connectivity
  bConnected_iGiD <- array( data = FALSE,
                            dim = c( length( entrez_iG ),
                                     length( diseases_iD ) ),
                            dimnames = list( genes = entrez_iG ,
                                             diseses = diseases_iD ) );
  for( iG in 1:dim( bConnected_iGiD ) ){
    bHasGene_iR <- dgn_cc_iRiC[ ,"geneId"] %in% entrez_iG[iG];
    uDs_hasGene_iD2 <- as.character( dgn_cc_iRiC[ bHasGene_iR, "diseaseName" ] );
    if( length( uDs_hasGene_iD2 ) > 0 ){
      bConnected_iGiD[iG,uDs_hasGene_iD2] <- TRUE;
    }
  }
  
  #Return
  return( bConnected_iGiD );
  
}

gwasGeneConnectivity = function( entrez_iG ) {

  # Get diseases from database
  gwas_iRiC <- read.table( file = "G:/Data/Gwas catalog/gwas_catalog_v1.0.1-downloaded_2015-11-26.tsv",
                           sep = "\t",
                           header = TRUE,
                           comment.char = "",
                           quote = "" );
  
  # Translate gene names into Symbolic
  symbol_jGiS <- mget( entrez_iG, 
                       org.Hs.egSYMBOL,
                       ifnotfound = NA );
  
  # Get only the GWAS that contain any of these genes
  symbol_jGs <- unlist( symbol_jGiS )
  bExists_iR <- sapply( as.character( gwas_iRiC[ ,"MAPPED_GENE"] ),
                        function( s ){
                          ss <- strsplit( s,
                                          " - " )[[1]];
                          return( any( ss %in% symbol_jGs ) )
                        } );
  gwas_c_iRiC <- gwas_iRiC[bExists_iR, ]
  
  # Find unfolded traits
  traits_iT <- unique( unlist( sapply( as.character( gwas_c_iRiC[ ,"MAPPED_TRAIT"] ),
                                       function( s ){
                                         ss <- strsplit( s,
                                                         ", ")
                                         return( ss )
                                       } ) ) );

  # Translate into connectivity matrix
  print( "Calculating GWAS cataloge conectivity to genes")
  disNames_iD <- as.character( unique( gwas_c_iRiC[ ,"MAPPED_TRAIT"] ) );
  bConnected_iGiT <- array( data = FALSE,
                            dim = c( length( entrez_iG ),
                                     length( traits_iT ) ),
                            dimnames = list( gene = entrez_iG,
                                             disease = traits_iT ) );
  stopifnot( all( names( symbol_jGiS ) == entrez_iG ) )
  for( iT in 1:dim( bConnected_iGiT )[2] ){
    # Find which studies have this trait
    bTrait_iR <- sapply( as.character( gwas_c_iRiC[ ,"MAPPED_TRAIT"] ),
                         function( s ){
                           ss <- strsplit( s,
                                           ", ")[[1]];
                           return( any( ss == traits_iT[iT] ) )
                         } );
    # Get all the genes of those stadies as connected to trait iT
    traitGenes_iG2 <- as.character( gwas_c_iRiC[ bTrait_iR, "MAPPED_GENE" ] );
    traitGenes_iG2 <- unlist( strsplit( traitGenes_iG2, " - ") )
    bGene_iG <- sapply( 1 : dim( bConnected_iGiT )[1],
                        function( iG ){
                          any( symbol_jGiS[[iG]] %in% traitGenes_iG2 );
                        } );
    stopifnot( any( bGene_iG ) )
    bConnected_iGiT[ bGene_iG, iT ] <- TRUE
  }
  
  # Return connectivity
  return( bConnected_iGiT )
  
}

# === Calculates the connectivity between genes according to where they share GO temrs

ontologyGeneConnectivity = function( geneNames_iG ){
  
  # Types of evidence (http://geneontology.org/page/guide-go-evidence-codes#ida)
  # NAS Non-traceable Author Statement
  # IEA Inferred from Electronic Annotation
  # IDA Inferred from Direct Assay
  # IBA Inferred from Biological aspect of Ancestor
  # TAS Traceable Author Statement
  # IGI Inferred from Genetic Interaction
  # IMP Inferred from Mutant Phenotype *
  # IC Inferred by Curator
  # IEP Inferred from Expression Pattern
  # IPI Inferred from Physical Interaction
  # ND No biological Data available
  # ISS Inferred from Sequence or structural Similarity **
  # EXP Inferred from Experiment
  #allowedEvidence_iE <- c( "IEA", "IDA", "IMP", "ISS" );
  #allowedEvidence_iE <- c( "NAS", "IEA", "IDA", "IBA", "TAS", "IGI", "IMP", "IC", "IEP", "IPI", "ND", "ISS", "EXP" );
  experimentalEvidence_iE <- c( "EXP", "IDA", "IPI", "IMP", "IGI", "IEP" ); #Experimental codes
  analysisEvidence_iE <- c( "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA" ); #***Analysis evidence codes
  statementEvidence_iE <- c( "TAS", "NAS" ); #**Statement evidence codes
  curatorEvidence_iE <- c( "IC", "ND" ); #Curator statement codes
  electronicEvidence_iE <- c( "IEA" ); #Electronic
  allowedEvidence_iE <- c( statementEvidence_iE );#analysisEvidence_iE, statementEvidence_iE );
  # Types of ontologies
  # MF - molecular function, 
  # BP - biological process, or
  # CC - cellular component
  allowedOntology_iO <- c( "MF", "BP", "CC" );
  #allowedOntology_iO <- c( "CC" );
  #geneNames_iG <- allGenes_iG[1:10]
  
  
  #ls("package:GO.db")
  
  # Get paths, evidence and ontology
  go_jGjPi3 <- mget(as.character(geneNames_iG), org.Hs.egGO,ifnotfound=NA);
  evidence_jGiP <- lapply( go_jGjPi3,
                           function( go_iT ){
                             if( ( length( go_iT ) == 1 ) && all( is.na( go_iT ) ) ){
                               return( NA );
                             }
                             sapply( go_iT,
                                     function( go ){
                                       go$Evidence;
                                     } )
                           } );
  ontology_jGiP <- lapply( go_jGjPi3,
                           function( go_iT ){
                             if( ( length( go_iT ) == 1 ) && all( is.na( go_iT ) ) ){
                               return( NA );
                             }
                             sapply( go_iT,
                                     function( go ){
                                       go$Ontology;
                                     } )
                           } );
  stopifnot( length( geneNames_iG ) == length( evidence_jGiP ) );
  stopifnot( length( geneNames_iG ) == length( ontology_jGiP ) );
  
  # Get only allowed evidence and ontology
  goID_jGiP <- lapply( 1:length( evidence_jGiP ),
                       function( iG ){
                         iE_iS <- match( evidence_jGiP[[iG]], allowedEvidence_iE );
                         iO_iS <- match( ontology_jGiP[[iG]], allowedOntology_iO );
                         return( unique( names( ontology_jGiP[[iG]][ !is.na( iE_iS ) & !is.na( iO_iS ) ] ) ) );
                       } );
  stopifnot( length( goID_jGiP ) == length( evidence_jGiP ) );
  stopifnot( length( goID_jGiP ) == length( ontology_jGiP ) );
  
  #   # Calculate higest categories
  #   aMF <- as.list( GOMFANCESTOR );
  #   aCC <- as.list( GOCCANCESTOR );
  #   aBP <- as.list( GOBPANCESTOR );
  #   offMF <- as.list( GOMFCHILDREN );
  #   offCC <- as.list( GOCCCHILDREN );
  #   offBP <- as.list( GOBPCHILDREN );
  #   topCategories1_iA <- c(   unlist( offMF["GO:0003674"] ),              
  #                             unlist( offCC["GO:0005575"] ), 
  #                             unlist( offBP["GO:0008150"] ) );
  #   topCategories2_iA <- c( unlist( offMF[ unlist( offMF["GO:0003674"] ) ] ), 
  #                           unlist( offCC[ unlist( offCC["GO:0005575"] ) ] ), 
  #                           unlist( offBP[ unlist( offBP["GO:0008150"] ) ] ) );
  #   topCategories2_iA <- setdiff( topCategories2_iA, 
  #                                 topCategories1_iA );
  #   topCategories2_iA <- topCategories1_iA;
  #   topCategories2_iA <- topCategories2_iA[ !is.na( topCategories2_iA ) ];
  #   topCategories2_iA <- topCategories2_iA[ topCategories2_iA != "all" ];
  #   allGO_iO <- unique( unlist(goID_jGiP) );
  #   bAncestor_iAiO <- array( data = FALSE,
  #                            dim = c( length( topCategories2_iA ),
  #                                     length( allGO_iO ) ),
  #                            dimnames = list( topCategory = topCategories2_iA,
  #                                             goTerm = allGO_iO ) );
  #   for( iO in 1:dim( bAncestor_iAiO )[2] ){
  #     ancestors_iA2 <- c( unlist( aMF[ allGO_iO[iO] ] ), 
  #                         unlist( aCC[ allGO_iO[iO] ] ),
  #                         unlist( aBP[ allGO_iO[iO] ] ) );
  #     bIsAncestro_iA <- !is.na( match( topCategories2_iA, ancestors_iA2 ) );
  #     bAncestor_iAiO[ , iO ] <- bIsAncestro_iA;
  #   }
  #   bAncestor_iAiO[ bAncestor_iAiO ] = 1;
  #   
  #   # Transform into highest categories
  #   goID2_jGiP <- lapply( goID_jGiP,
  #                         function( s_iS ){
  #                           if( length( s_iS ) == 0 ){
  #                             return( NA );
  #                           }
  #                           if( length( s_iS ) == 1 ){
  #                             numHist_iA <- bAncestor_iAiO[ ,s_iS];   
  #                           }else{
  #                             numHist_iA <- rowSums( bAncestor_iAiO[ ,s_iS] );                            
  #                           }
  #                           return( dimnames( bAncestor_iAiO )[[1]][ numHist_iA > 0 ] );
  #                         } );
  #   stopifnot( length( goID2_jGiP ) == length( goID_jGiP ) );
  
  goID2_jGiP <- goID_jGiP;
  
  # Clean unlinked terms
  goID2_jGiP <- lapply( goID2_jGiP,
                        function( s_iS ){
                          if( length( s_iS ) == 0 ){
                            return( NA );
                          }
                          return( s_iS );
                        } );
  stopifnot( length( goID2_jGiP ) == length( evidence_jGiP ) );
  stopifnot( length( goID2_jGiP ) == length( ontology_jGiP ) );
  
  #Calculate connectivity
  bConnected_iGiP <- calculateConnectivity( geneNames_iG,
                                            goID2_jGiP)
  
  # Set names properly
  library(GO.db)
  goterms <- Term(GOTERM)
  dimnames( bConnected_iGiP )[[2]] <- sapply( 1:dim( bConnected_iGiP )[2],
                                              function( iG ){
                                                paste( dimnames( bConnected_iGiP )[[2]][iG],
                                                       goterms[ dimnames( bConnected_iGiP )[[2]][iG] ] )
                                              } );
  
  #Return
  return( bConnected_iGiP );
}

# === Internal function to calculate the connectivity between genes and paths, as described in "pathNames_jGiP"
calculateConnectivity = function( geneNames_iG,
                                  pathNames_jGiP ){
  
  # Check that things are OK
  if( length(pathNames_jGiP) != length(geneNames_iG) ){
    stop("The numero of genes has changed. There is something going wrong.")
  }
  pathNames_iP <- vector("character",length=0);
  for( iG in 1:length(geneNames_iG) ){
    pathNames_iP <- c(pathNames_iP,pathNames_jGiP[[iG]]);
  }
  pathNames_iP <- unique(pathNames_iP);
  pathNames_iP <- pathNames_iP[!is.na(pathNames_iP)];
  
  # Find connectivity
  bConnected_iGiP <- array(FALSE,dim=c(length(geneNames_iG),length(pathNames_iP)));
  for( iG in 1:length(geneNames_iG) ){
    iPs <- match(pathNames_jGiP[[iG]],pathNames_iP);
    bConnected_iGiP[iG,iPs] <- TRUE;
  }
  
  #Return results
  dimnames(bConnected_iGiP) <- list( gene = geneNames_iG , path = pathNames_iP );
  return(bConnected_iGiP);
  
}