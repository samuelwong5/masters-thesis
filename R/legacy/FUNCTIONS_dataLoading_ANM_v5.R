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

synthetiseData = function( target_iSiG,
                           source_iSiP ,
                           iReps = 10 ){
  #
  # Creates synthetic copies of data "target_iSiG" by predicting each
  # variable in "target_iSiG" with multiple imputations from "source_iSiP"
  #
  
  # Calculate neccessary data structures 
  target_iSiG <- as.matrix( target_iSiG )
  source2_iSiP <- rbind( source_iSiP,
                         source_iSiP )
  source2_iSiP <- cbind( imp = NA,
                         source2_iSiP )
  target_si_iSiG <- repmat( target_iSiG, 
                            n = iReps, 
                            m = 1 )
  
  # Create the synthetic data
  for( iG in 1:dim( target_iSiG )[2] ){
    # Add to be imputed column
    print( paste( "Imputing variable iG ", iG, " of ", dim( target_iSiG )[2] ) )
    source2_iSiP[ 1:dim( source_iSiP )[1], "imp" ] <- target_iSiG[ ,iG]
    # Impute column gen_ckd_iSiG source2_iSiP
    rM <- mice( source2_iSiP,
                m = iReps,
                maxit = 10,
                meth = 'norm.boot',
                seed = 500,
                printFlag = FALSE )
    # Recover imputed column
    target_i_iSiG <- complete( rM,
                               action = "long" )
    bImputed_iS <- 1:dim( source2_iSiP )[1] > dim( source_iSiP )[1]
    bImputed_iS <- rep( bImputed_iS,
                        iReps )
    stopifnot( length( bImputed_iS ) == dim( target_i_iSiG )[1] )
    stopifnot( sum( bImputed_iS ) == dim( target_si_iSiG )[1] )
    target_si_iSiG[ ,iG] <- target_i_iSiG[ bImputed_iS, "imp" ]
  }
  
  # Return results
  return( target_si_iSiG )
}


mergeCohorts = function( dat1_iSiG,
                         dat2_iSiG ){
  #
  # Merge cohorts with intersecting population and columns
  #
  
  # Convert data into matrix
  if( class( dat1_iSiG ) != "matrix" ){
    warning( "Data needs to be a matrix. Converting automatically into matrix now.")
    dat1_iSiG <- as.matrix( dat1_iSiG )
  }
  if( class( dat2_iSiG ) != "matrix" ){
    warning( "Data needs to be a matrix. Converting automatically into matrix now.")
    dat2_iSiG <- as.matrix( dat2_iSiG )
  }
  
  # Calculate dimentions of the common array
  print( "Calculate merging dimentions" )
  stopifnot( !any( duplicated( dimnames( dat1_iSiG )[[1]] ) ) )
  stopifnot( !any( duplicated( dimnames( dat2_iSiG )[[1]] ) ) )
  nSs <- unique( c( dimnames( dat2_iSiG )[[1]],
                    dimnames( dat1_iSiG )[[1]] ) )
  stopifnot( !any( duplicated( dimnames( dat1_iSiG )[[2]] ) ) )
  stopifnot( !any( duplicated( dimnames( dat2_iSiG )[[2]] ) ) )
  nGs <- unique( c( dimnames( dat2_iSiG )[[2]],
                    dimnames( dat1_iSiG )[[2]] ) )
  
  # Copy data into a common array
  print( "Copying data into merging array" )
  datM_iSiG <- array( data = NA,
                      dim = c( length( nSs ),
                               length( nGs ) ),
                      dimnames = list( nSs = nSs,
                                       nGs = nGs ) )
  datM_iSiG[ dimnames( dat1_iSiG )[[1]],
             dimnames( dat1_iSiG )[[2]] ] <- dat1_iSiG
  datM_iSiG[ dimnames( dat2_iSiG )[[1]],
             dimnames( dat2_iSiG )[[2]] ] <- dat2_iSiG
  
  # Return result
  return( datM_iSiG )
}

fixSampleNames = function( nSs,
                           bStatus = TRUE ){
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
    if( bStatus ){
      return( toupper( p( site, status, number ) ) );
    }else{
      return( toupper( p( site, number ) ) );
    }
  }
  nSs_fixed <- sapply( nSs,
                       homogenizeNames );
  stopifnot( all( !is.na( nSs_fixed ) ) )
  stopifnot( all( nSs_fixed != "" ) );
  stopifnot( all( nSs_fixed != "NA" ) );
  stopifnot( all( nSs_fixed != "<NA>" ) );
  return( nSs_fixed )
  
}

fixSampleNames2 = function( nSs ){
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
    s <- gsub( "CTL|ADC|MCI|\\.",
               "",
               s );
    r_F <- regexpr( "[A-z]+",
                    s,
                    perl = TRUE );
    site <- substr( s,
                    r_F[1],
                    r_F[1] + attr( r_F, "match.length" ) - 1 );
    #     site <- substr( site,
    #                     1,
    #                     min( 3,
    #                          nchar( site ) ) );
    r_F <- regexpr( "[0-9]+",
                    s,
                    perl = TRUE );
    number <- substr( s,
                      r_F[1],
                      r_F[1] + attr( r_F, "match.length" ) - 1 );
    number <- as.numeric( number );   
    return( toupper( p( site, number ) ) );
  }
  nSs_fixed <- sapply( nSs,
                       homogenizeNames );
  stopifnot( all( !is.na( nSs_fixed ) ) )
  stopifnot( all( nSs_fixed != "" ) );
  stopifnot( all( nSs_fixed != "NA" ) );
  stopifnot( all( nSs_fixed != "<NA>" ) );
  return( nSs_fixed )
  
}



loadANM_demo = function( subjectID_iS,
                         correctNames = TRUE ){
  
  # Load demographics
  print( "Loading demographics" );
  diagnoses_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/CaseHistory/data.csv",header=FALSE,stringsAsFactors=FALSE);
  conclussi_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/Conclusion/data.csv",header=FALSE,stringsAsFactors=FALSE);
  demograph_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/Demographics/data.csv",header=FALSE,stringsAsFactors=FALSE);
  apoe_iSiC         <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/APOE/data.csv",header=FALSE,stringsAsFactors=FALSE);  
  diagnoses_iSiC    <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/CaseHistory/data.csv",header=FALSE,stringsAsFactors=FALSE);
  adas_iSiC         <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/ADAS-COG/data.csv",header=TRUE,stringsAsFactors=FALSE);
  cdr_iSiC         <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/CDR-Subject/data.csv",header=TRUE,stringsAsFactors=FALSE);
  mmse_iSiC         <- read.csv("G:/Data/ANM/Data from Richard/CohortExplorer-141277251712804/MMSE/data.csv",header=TRUE,stringsAsFactors=FALSE);
  imaging_iSiC      <- read.csv("G:/Project - Me Elena 1 GSK/Data/All data for analysis.csv",header=TRUE,stringsAsFactors=FALSE);  
  
  # Correct names
  print( "Correcting names" );
  subjectID_iS <- fixSampleNames( nSs = subjectID_iS );
  diagnoses_iSiC[,1] <- fixSampleNames( nSs = diagnoses_iSiC[,1] );
  conclussi_iSiC[,1] <- fixSampleNames( nSs = conclussi_iSiC[,1] );
  demograph_iSiC[,1] <- fixSampleNames( nSs = demograph_iSiC[,1] );
  apoe_iSiC[,1] <- fixSampleNames( nSs = apoe_iSiC[,1] );
  adas_iSiC[,1] <- fixSampleNames( nSs = adas_iSiC[,1] );
  cdr_iSiC[,1] <- fixSampleNames( nSs = cdr_iSiC[,1] );
  mmse_iSiC[,1] <- fixSampleNames( nSs = mmse_iSiC[,1] );
  imaging_iSiC[,"Code"] <- fixSampleNames( nSs = imaging_iSiC[,"Code"] );
  
  # Get demographics
  print( "Curating data" );
  features_iF <- c("subjectID",  #1
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
                   "PatientFinalState_detail",
                   "WholeBrain",
                   "entorhinal",
                   "entorhinal_normalised",
                   "HippLeft",
                   "HippRight",
                   "WholeBrain",
                   "HippLeft_normalised",
                   "HippRight_normalised",
                   "Hipp",
                   "Hipp_normalised",
                   "mmse",
                   "mmse_change",
                   "adas",
                   "adas_change" )
  demographics_iSiC           <- array( data="na" ,
                                        dim = c( length( subjectID_iS ), 
                                                 length( features_iF ) ),
                                        dimnames = list( subjectID = subjectID_iS,
                                                         feature = features_iF ) );
  bNoDemographics_iSiD            <- array( data = TRUE, 
                                            dim = c( length( subjectID_iS ), 5 ) );
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
      questionsDate <- as.numeric(substring(demograph_iSiC[bSs,2],1,4))[1];
      demographics_iSiC[iS,4] <- as.character(demograph_iSiC[bSs,22])[1];
      demographics_iSiC[iS,5] <- questionsDate-birthDate[1];
      demographics_iSiC[iS,9] <- as.character(demograph_iSiC[bSs,7])[1];
      demographics_iSiC[iS,10] <- as.character(demograph_iSiC[bSs,9])[1];
      demographics_iSiC[iS,11] <- as.character(demograph_iSiC[bSs,11])[1];
      bNoDemographics_iSiD[iS,3] <- FALSE;
    }
    
    # Get demographics from apoe file
    bSs <- apoe_iSiC[,1] == subjectID_iS[iS];
    if( any(bSs) ){
      demographics_iSiC[iS,6] <- as.character(apoe_iSiC[bSs,3])[1];
      bNoDemographics_iSiD[iS,4] <- FALSE;
    }
    
    # Get demographics from imaging file
    bSs <- imaging_iSiC[ ,"Code"] == subjectID_iS[iS];
    stopifnot( sum( bSs ) <= 1 )
    if( any(bSs) ){
      demographics_iSiC[iS,"WholeBrain"] <- imaging_iSiC[ bSs, "WholeBrain" ];
      demographics_iSiC[iS,"entorhinal"] <- imaging_iSiC[ bSs, "entorhinal" ];
      demographics_iSiC[iS,"entorhinal_normalised"] <- imaging_iSiC[ bSs, "entorhinal" ] / imaging_iSiC[ bSs, "WholeBrain" ];
      demographics_iSiC[iS,"HippLeft"] <- imaging_iSiC[ bSs, "HippLeft" ];
      demographics_iSiC[iS,"HippRight"] <- imaging_iSiC[ bSs, "HippRight" ];
      demographics_iSiC[iS,"Hipp"] <- mean( c( imaging_iSiC[ bSs, "HippLeft" ],
                                               imaging_iSiC[ bSs, "HippRight" ] ) );
      demographics_iSiC[iS,"HippLeft_normalised"] <- imaging_iSiC[ bSs, "HippLeft_normalised" ];
      demographics_iSiC[iS,"HippRight_normalised"] <- imaging_iSiC[ bSs, "HippRight_normalised" ];
      demographics_iSiC[iS,"Hipp_normalised"] <- mean( c( imaging_iSiC[ bSs, "HippLeft_normalised" ],
                                                          imaging_iSiC[ bSs, "HippRight_normalised" ] ) );
    }    
    
    # Get demographics from imaging file
    bSs <- mmse_iSiC[ ,"entity_id"] == subjectID_iS[iS];
    if( any(bSs) ){
      # Get mean MMSE
      date_iV <- as.numeric( as.Date(mmse_iSiC[bSs,"QuestionnaireRun.timeStart"]) );
      visitNum_iV <- as.numeric(mmse_iSiC[bSs,"visit"]);
      mmse_iV <- as.character(mmse_iSiC[bSs,"MMSE_Total"]);
      demographics_iSiC[iS,"mmse"] <- mmse_iV[visitNum_iV==min(visitNum_iV)];
      # Get change
      if( sum( bSs ) > 2 ){
        r <- lm( mmse_iV ~ date_iV )
        demographics_iSiC[iS,"mmse_change"] <- coef( summary(r) )["date_iV","Estimate"]
      }
    }     
    
    # Get demographics from imaging file
    bSs <- adas_iSiC[ ,"entity_id"] == subjectID_iS[iS];
    if( any(bSs) ){
      # Get mean MMSE
      visitNum_iV <- as.numeric(adas_iSiC[bSs,"visit"]);
      mmse_iV <- as.character(adas_iSiC[bSs,"ADAS_COG_Total"]);
      demographics_iSiC[iS,"adas"] <- mmse_iV[visitNum_iV==min(visitNum_iV)];
      # Get change
      date_iV <- adas_iSiC[bSs,"QuestionnaireRun.timeStart"]
      mmse_iV <- mmse_iV[ date_iV != "" ]
      date_iV <- date_iV[ date_iV != "" ]
      if( length( date_iV ) > 2 ){
        date_iV <- as.numeric( as.Date( date_iV ) );
        r <- lm( mmse_iV ~ date_iV )
        demographics_iSiC[iS,"adas_change"] <- coef( summary(r) )["date_iV","Estimate"]
      }
    }       
    
    # Get demographics from other sources
    demographics_iSiC[iS,7] <- substr(subjectID_iS[iS],1,3);
  }
  
  # Return data
  return( demographics_iSiC );
  
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
