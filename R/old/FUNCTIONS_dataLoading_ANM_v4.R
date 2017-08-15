#
# Actor class to load and save DataVST objects
#
# Responsabilities
# - Load data from ANM different sources (batch 1, batch 2...)
# - Save data
#
# v2 = Demographics unmerged
# v3 = Adding SNP file split
# v4 = Added function to fix screwed up sample names
#
#For installing further libraries: 
#source("http://bioconductor.org/biocLite.R")
#biocLite("wathever")
#library(wathever)

#DEPENDENCIES
source("FUNCTIONS_utils.R")
library(gdata);
library(pracma);

fixSampleNames = function( nSs ){
  # This function tries to fix the very awfully formated sample names of the different versions of ANM
  
  # Homogenize stupid iSs names
  homogenizeNames = function( s ){
    status <- "";
    if( substr( s, 4, 4 ) == "C" ){
      status <- "CTL";
    }else if( substr( s, 4, 4 ) == "A" ){
      status <- "ADC";
    }else if( substr( s, 4, 4 ) == "M" ){
      status <- "MCI";
    }   
    s <- gsub( "CTL|ADC|AD|MCI|\\.",
               "",
               s );
    r_F <- regexpr( "[A-z]+",
                    s,
                    perl = TRUE );
    site <- substr( s,
                    r_F[1],
                    r_F[1] + attr( r_F, "match.length" ) - 1 );
    site <- substr( site,
                    1,
                    min( 3,
                         nchar( site ) ) );
    r_F <- regexpr( "[0-9]+",
                    s,
                    perl = TRUE );
    number <- substr( s,
                      r_F[1],
                      r_F[1] + attr( r_F, "match.length" ) - 1 );
    number <- as.numeric( number );   
    toupper( p( site, status, number ) );
  }
  nSs_fixed <- sapply( nSs,
                       homogenizeNames );
  stopifnot( all( !is.na( nSs_fixed ) ) )
  stopifnot( all( nSs_fixed != "" ) );
  stopifnot( all( nSs_fixed != "NA" ) );
  stopifnot( all( nSs_fixed != "<NA>" ) );
  return( nSs_fixed )
  
}

loadANM_snp = function( ){
  
  # Load file (THIS TAKES A LONG WHILE)
  library( snpMatrix );
  sample <- read.plink( bed = "G:/Data/ANM/anm_imputed/anm_batch1_batch2_merged.bed", 
                        bim = "G:/Data/ANM/anm_imputed/anm_batch1_batch2_merged.bim", 
                        fam = "G:/Data/ANM/anm_imputed/anm_batch1_batch2_merged.fam" );
  
  # Save into sentible pieces
  map <- sample$map;
  save( file = "causalNetwork_v1_map.Rdata",
        list = c( "map" ) );
  fam <- sample$fam;
  save( file = "causalNetwork_v1_fam.Rdata",
        list = c( "fam" ) );
  chunks_iK <- round( dim( sample$genotypes )[2] * ( 0 : 100 ) / 100 );
  chunks_iK[1] <- 1;
  for( iK in 1:( length( chunks_iK ) - 1 ) ){
    print( p( "Saving chung ", iK ) );
    snp_iSiN <- sample$genotypes[ , chunks_iK[iK] : chunks_iK[iK+1] ];
    save( file = p( "causalNetwork_v1_snp_chunk", iK, ".Rdata" ),
          list = c( "snp_iSiN" ) );    
  }
  
}

# === Load data from different ANM batches
loadANM_expression = function( batch = 1 ) {
  
  # Load data batch 1
  print( "Loading Files" );
  if( batch == 1 ){
    array()
    expression_iGiS    <- read.table("G:/Project WT-Inf/Data/ANM batch 1/Pre_processed_gene_expression_ANM_BATCH_1_UPDATED_BEFORE_CUTOFF.csv",
                                     row.names = 1, 
                                     header = TRUE, 
                                     colClasses = c(NA,"numeric"),
                                     sep = ",");
    subjects_iSiC      <- read.csv( "G:/Project WT-Inf/Data/ANM batch 1/Subject_identifier_gene_expression_ANM_BATCH_1.csv",
                                    header = TRUE,
                                    stringsAsFactors = FALSE );
  }else if( batch  == 2 ){
    expression_iGiS    <- read.csv( "G:/Project WT-Inf/Data/ANM batch 2/Pre_processed_gene_expression_ANM_BATCH_2_UPDATED_BEFORE_CUTOFF.csv",
                                    row.names = 1, 
                                    header = TRUE, 
                                    colClasses = c(NA,"numeric"),
                                    sep = ",")
    subjects_iSiC      <- read.csv( "G:/Project WT-Inf/Data/ANM batch 2/Subject_Identifier_BATCH_2.csv",
                                    header = TRUE,
                                    stringsAsFactors = FALSE);
    dimnames( subjects_iSiC )[[2]] <- mapvalues( dimnames( subjects_iSiC )[[2]],
                                                 c( "Subject.ID",
                                                    "Chip.ID"),
                                                 c( "EU_Code",
                                                    "ID") );
  }else{
    stop( p( "Batch number ",
             batch,
             " does not exist" ) );
  }
  genes1_iGiC         <- read.csv( "G:/Project WT-Inf/Data/ANM batch 1/Gene_probe_identifier_BATCH_1.csv",
                                   header = TRUE,
                                   stringsAsFactors = FALSE );
  genes2_iGiC         <- read.csv( "G:/Project WT-Inf/Data/ANM batch 2/Gene_probe_identifier_UPDATED_BATCH_2.csv",
                                   header = TRUE,
                                   stringsAsFactors = FALSE );
  
  # Merge the ID-name tables for the genes
  print( "Translating names" );
  stopifnot( all( !duplicated( genes1_iGiC$nuID ) ) );
  stopifnot( all( !duplicated( genes2_iGiC$nuID ) ) );
  dimnames( genes1_iGiC )[[1]] <- genes1_iGiC$nuID;
  dimnames( genes2_iGiC )[[1]] <- genes2_iGiC$nuID;
  stopifnot( all( genes1_iGiC[ genes2_iGiC$nuID, "ENTREZ_GENE_ID" ] == genes2_iGiC$ENTREZ_GENE_ID,
                  na.rm = TRUE ) );
  stopifnot( all( genes2_iGiC[ genes1_iGiC$nuID, "ENTREZ_GENE_ID" ] == genes1_iGiC$ENTREZ_GENE_ID,
                  na.rm = TRUE ) )
  genes_iGiC <- rbind( genes1_iGiC,
                       genes2_iGiC );
  genes_iGiC <- unique( genes_iGiC );
  stopifnot( all( !duplicated( genes_iGiC$nuID ) ) );
  dimnames( genes_iGiC )[[1]] <- genes_iGiC$nuID;
  
  # Transform genes names into entrez, eliminating genes without entrez equivalent
  library( org.Hs.eg.db )
  namesRaw_iG <- dimnames( expression_iGiS )[[1]];
  namesEnt_iG <- genes_iGiC[ namesRaw_iG, "ENTREZ_GENE_ID" ];
  expression_c_iGiS <- expression_iGiS[ !is.na( namesEnt_iG ), ];
  expression_c_iGiS$entrezName <- namesEnt_iG[ !is.na( namesEnt_iG ) ];
  stopifnot( sum( !is.na( namesEnt_iG ) ) > 6000 );
  
  # Translate subject names into standard IDs
  stopifnot( all( !duplicated( subjects_iSiC$ID ) ) );
  dimnames( subjects_iSiC )[[1]] <- subjects_iSiC$ID;
  namesRaw_iS <- dimnames( expression_c_iGiS )[[2]];
  namesRaw_iS <- sapply( namesRaw_iS,
                         function( s ){
                           substr( s,
                                   2,
                                   nchar( s ) );
                         } );
  namesIDs_iS <- subjects_iSiC[ namesRaw_iS, "EU_Code" ];
  namesIDs_iS[ namesRaw_iS == "ntrezName" ] <- "entrezName";
  stopifnot( all( !is.na( namesIDs_iS ) ) );
  dimnames( expression_c_iGiS )[[2]] <- namesIDs_iS;
  
  # Load demographics
  print( "Adding demographics" );
  diagnoses_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/CaseHistory/data.csv",header=FALSE,stringsAsFactors=FALSE);
  conclussi_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/Conclusion/data.csv",header=FALSE,stringsAsFactors=FALSE);
  demograph_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/Demographics/data.csv",header=FALSE,stringsAsFactors=FALSE);
  apoe_iSiC         <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/APOE/data.csv",header=FALSE,stringsAsFactors=FALSE);  
  
  # Get demographics
  demographics_iSiC           <- array( data="na" ,
                                        dim = c( dim( expression_c_iGiS )[2] , 
                                                 13 ),
                                        dimnames = list( subjectID = dimnames( expression_c_iGiS )[[2]],
                                                         feature = c("subjectID",  #1
                                                                     "PatientState",
                                                                     "PatientFinalState",
                                                                     "gender",
                                                                     "age",  #5
                                                                     "APOE",
                                                                     "centre",
                                                                     "ageOnset",
                                                                     "education",
                                                                     "profession", #10
                                                                     "accommodation",
                                                                     "PatientState_detail",
                                                                     "PatientFinalState_detail") ) );
  bNoDemographics_iSiD            <- array( data = TRUE, 
                                            dim = c( dim( expression_c_iGiS )[2] , 5 ) );
  subjectID_iS <- dimnames( expression_c_iGiS )[[2]];
  for( iS in 1:length( subjectID_iS ) ){
    
    # Ignore if we are in the row with the gene names
    if( subjectID_iS[iS] == "entrezName" ){
      next;
    }
    
    # Get demographics from diagnoses file
    bSs <- diagnoses_iSiC[,1] == subjectID_iS[iS];
    if( any(bSs) ){
      demographics_iSiC[iS,1] <- subjectID_iS[iS];
      demographics_iSiC[iS,8] <- as.character( min( as.numeric( diagnoses_iSiC[bSs,10] ),
                                                    na.rm=TRUE ) );
      bNoDemographics_iSiD[iS,1] <- FALSE;
    }
    
    # Get data from conclussions file
    bSs <- conclussi_iSiC[,1] == subjectID_iS[iS];
    if( any(bSs) ){
      visitNum_iV <- as.numeric(conclussi_iSiC[bSs,2]);
      diagnosis_iV <- as.character(conclussi_iSiC[bSs,7]);
      demographics_iSiC[iS,2] <- diagnosis_iV[visitNum_iV==min(visitNum_iV)];
      demographics_iSiC[iS,3] <- diagnosis_iV[visitNum_iV==max(visitNum_iV)];
      diagnosisDetail_iV <- as.character(conclussi_iSiC[bSs,12]);
      demographics_iSiC[iS,12] <- diagnosisDetail_iV[visitNum_iV==min(visitNum_iV)];
      demographics_iSiC[iS,13] <- diagnosisDetail_iV[visitNum_iV==max(visitNum_iV)];
      bNoDemographics_iSiD[iS,2] <- FALSE;
    }
    
    # Get demographics from demographics file
    bSs <- demograph_iSiC[,1] == subjectID_iS[iS];
    if( any(bSs) ){
      birthDate <- as.numeric(demograph_iSiC[bSs,21]);
      questionsDate <- as.numeric(substring(demograph_iSiC[bSs,2],1,4));
      demographics_iSiC[iS,4] <- as.character(demograph_iSiC[bSs,22]);
      demographics_iSiC[iS,5] <- questionsDate-birthDate;
      demographics_iSiC[iS,9] <- as.character(demograph_iSiC[bSs,7]);
      demographics_iSiC[iS,10] <- as.character(demograph_iSiC[bSs,9]);
      demographics_iSiC[iS,11] <- as.character(demograph_iSiC[bSs,11]);
      bNoDemographics_iSiD[iS,3] <- FALSE;
    }
    
    # Get demographics from apoe file
    bSs <- apoe_iSiC[,1] == subjectID_iS[iS];
    if( any(bSs) ){
      demographics_iSiC[iS,6] <- as.character(apoe_iSiC[bSs,3]);
      bNoDemographics_iSiD[iS,4] <- FALSE;
    }
    
    # Get demographics from other sources
    demographics_iSiC[iS,7] <- substr(subjectID_iS[iS],1,3);
  }
  
  # Merge demographics with gene expression
  expSubject_iS <- dimnames( expression_c_iGiS )[[2]];
  expAndDem_c_iGDiS <- rbind( expression_c_iGiS,
                              t( demographics_iSiC[ expSubject_iS, ] ) );
  
  # PATCH: Unmerge demographics
  exp_iSiG <- t( expAndDem_c_iGDiS[ !is.na( expAndDem_c_iGDiS$entrezName ), ] );
  dimnames( exp_iSiG )[[2]] <- exp_iSiG[ "entrezName", ];
  exp_iSiG <- exp_iSiG[ -dim( exp_iSiG )[1], ];
  dem_iSiD <- t( expAndDem_c_iGDiS[ c( "subjectID",
                                       "PatientState",            
                                       "PatientFinalState",
                                       "gender",
                                       "age",
                                       "APOE",                    
                                       "centre",
                                       "ageOnset",
                                       "education",
                                       "profession",              
                                       "accommodation",
                                       "PatientState_detail",
                                       "PatientFinalState_detail" ), -dim( expAndDem_c_iGDiS )[2] ] );
  stopifnot( all( dimnames( dem_iSiD )[[1]] == dimnames( exp_iSiG )[[1]] ) );
  
  # Return result
  r <- list( dem_iSiD = dem_iSiD,
             exp_iSiG = apply( exp_iSiG,
                               c( 1, 2 ),
                               as.numeric ) );
  return( r );
}


# === Load data from different ANM batches
loadANM_protein = function( ) {
  
  # Load data 
  print( "Loading Files" );
  data_iSiC    <- read.table( "G:/Data/ANM/SomaLogic.AD.Ben/Somalogic.AD.Ben.csv",
                              row.names = 1, 
                              header = TRUE, 
                              sep = "," );
  meta_iRiP    <- read.csv( "G:/Data/ANM/SomaLogic.AD.Ben/Protein_names.csv",
                            header = FALSE,
                            row.names = 1 );
  
  # Separate demographics from potrs
  prot_iSiG <- data_iSiC[ ,-c(1:6)];
  demo_iSiD <- data_iSiC[ , c(1:6)];
  
  # Eliminate stupid factors
  protID_iP <- lapply( meta_iRiP["Target", ],
                       as.character );
  homogenizeNames = function( s_iS ){
                      s_iS <- gsub( "\\.| |\\-|\\:|\\,|\\/|0",
                            "",
                            as.character( s_iS ),
                            ignore.case = TRUE );
                      s_iS <- tolower( s_iS );
                      return( s_iS );
                    };
  protID_iP <- homogenizeNames( protID_iP );
  protID_iG <- homogenizeNames( dimnames( prot_iSiG )[[2]] );
  bFound_iG <- protID_iG %in% protID_iP;
  table_F <- table( bFound_iG )
  print( p( "Identified ", table_F["TRUE"], " of ", sum( table_F ), " proteins" ) );
  
  # Get only identified proteins
  prot_c_iSiG <- prot_iSiG[ , bFound_iG ];
  protID_c_iG <- protID_iG[ bFound_iG ];
  
  # Translate protein names
  library( plyr );
  gen_iP <- sapply( meta_iRiP["EntrezGeneID", ],
                    function( s ){
                      s <- gsub( "\\,",
                                 " ",
                                 s );
                      strsplit( s,
                                " " )[[1]][[1]]
                    } );
  gen_c_iG <- mapvalues( protID_c_iG,
                         protID_iP,
                         gen_iP,
                         warn_missing = FALSE );
  stopifnot( length( gen_c_iG ) == dim( prot_c_iSiG )[2] );
  dimnames( prot_c_iSiG )[[2]] <- gen_c_iG;
  
  # Return result
  stopifnot( all( "numeric" == apply( prot_c_iSiG, c(2), class ) ) );
  stopifnot( all( dimnames( prot_c_iSiG )[[1]] == dimnames( demo_iSiD )[[1]] ) );
  r <- list( dem_iSiD = demo_iSiD,
             exp_iSiG = prot_c_iSiG );
  return( r );
}
