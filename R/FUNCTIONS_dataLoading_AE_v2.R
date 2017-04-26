#
# Actor class to load and save DataVST objects
#
# v2 = Now getDemographics_SDRF also takes into account that sometimes each subject appears several times in the SDRF file
#
# Responsabilities
# - Load data from Array Express sources
# - Save data
#
#For installing further libraries: 
#source("http://bioconductor.org/biocLite.R")
#biocLite("wathever")
#library(wathever

#source("CLASS_DataVST.R")

loadAE_48350 = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress GEOD48350's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "C:/Users/alejon/Documents/samuel/masters/data/E-GEOD-48350/samples",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "C:/Users/alejon/Documents/samuel/masters/data/E-GEOD-48350/A-AFFY-44.adf.txt",
                                         columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Change gene names from RefSeq into Entrez
  library( org.Hs.eg.db );
  entrez_iG <- as.list( org.Hs.egREFSEQ2EG );
  iGs_exist <- names( entrez_iG ) %in% dimnames( expression_e_iGiS )[[1]];
  entrez_iG <- entrez_iG[ iGs_exist ]; 
  entrez_iG <- unlist( entrez_iG );
  expression_te_iGiS <- expression_e_iGiS[ names( entrez_iG ), ];
  dimnames( expression_te_iGiS )[[1]] <- entrez_iG;
  
  # Save stuff
  save( file = "loadAE_E-GEOD-48350_security.Rdata", list = c( "expression_te_iGiS" ) );
  rm( expression_iGiS );
  rm( expression_e_iGiS );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_te_iGiS, 
                                             factorName_sdrf_iF = c( "Term.Source.REF.", 
                                                                     "Characteristics.age.yrs.", 
                                                                     "Characteristics.sex.",
                                                                     "Characteristics.braak.stage."), 
                                             factorName_iF = c( "brainRegion", "age", "gender", "ad" ), 
                                             fileSDRF = "C:/Users/alejon/Documents/samuel/masters/data/E-GEOD-48350/E-GEOD-48350.sdrf.txt" );
  demographics_iSiC[ ,"age"] <- sapply( demographics_iSiC[ ,"age"],
                                        function( s ){
                                          return( strsplit( s, " ")[[1]][1] );
                                        } );
  demographics_iSiC[ demographics_iSiC[ ,"age"] == ">90", "age" ] <- 90;
  
  
  # Put all data together
  stopifnot( all( dimnames( demographics_iSiC )[[1]] %in% dimnames( expression_te_iGiS )[[2]] ) );
  all_iSiC <- as.data.frame( t( expression_te_iGiS ) );
  all_iSiC <- cbind( demographics_iSiC,
                     all_iSiC[ dimnames( demographics_iSiC )[[1]], ] );
  all_iSiC$age <- as.numeric( as.character( all_iSiC$age ) );
  
  # Return
  return( all_iSiC );
  
}

loadAE_hamill = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Hamill's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "G:/Project EMIF/AEdata/Hamill5281/E-GEOD-5281.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "G:/Project EMIF/AEdata/Hamill5281/A-AFFY-44.adf.txt",
                                         columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Change gene names from RefSeq into Entrez
  library( org.Hs.eg.db );
  entrez_iG <- as.list( org.Hs.egREFSEQ2EG );
  iGs_exist <- names( entrez_iG ) %in% dimnames( expression_e_iGiS )[[1]];
  entrez_iG <- entrez_iG[ iGs_exist ]; 
  entrez_iG <- unlist( entrez_iG );
  expression_te_iGiS <- expression_e_iGiS[ names( entrez_iG ), ];
  dimnames( expression_te_iGiS )[[1]] <- entrez_iG;
  
  # Save stuff
  save( file = "loadAE_hamill_security.Rdata", list = c( "expression_te_iGiS" ) );
  rm( expression_iGiS );
  rm( expression_e_iGiS );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_te_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics.organ.region.", 
                                                                     "Characteristics.age.", 
                                                                     "Characteristics.sex.",
                                                                     "Characteristics.disease.state."), 
                                             factorName_iF = c( "brainRegion", "age", "gender", "ad" ), 
                                             fileSDRF = "G:/Project EMIF/AEdata/Hamill5281/E-GEOD-5281.sdrf.txt" );
  demographics_iSiC[ ,"age"] <- sapply( demographics_iSiC[ ,"age"],
                                        function( s ){
                                          return( strsplit( s, " ")[[1]][1] );
                                        } );
  demographics_iSiC[ demographics_iSiC[ ,"age"] == ">90", "age" ] <- 90;
  
  
  # Put all data together
  stopifnot( all( dimnames( demographics_iSiC )[[1]] %in% dimnames( expression_te_iGiS )[[2]] ) );
  all_iSiC <- as.data.frame( t( expression_te_iGiS ) );
  all_iSiC <- cbind( demographics_iSiC,
                     all_iSiC[ dimnames( demographics_iSiC )[[1]], ] );
  all_iSiC$age <- as.numeric( as.character( all_iSiC$age ) );
  
  # Return
  return( all_iSiC );
  
}

loadAE_zhu = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Zhu's. ",sep="") );
  
  # Get gene expression files
  expression_1_iGiS <- getPreprocessedFiles( folder = "G:/Project EMIF/AEdata/Zhu4372/E-GEOD-44772.processed.1",
                                             nameColum = "Reporter.Identifier",
                                             expColum = "VALUE" );
  expression_2_iGiS <- getPreprocessedFiles( folder = "G:/Project EMIF/AEdata/Zhu4372/E-GEOD-44772.processed.2",
                                             nameColum = "Reporter.Identifier",
                                             expColum = "VALUE" );
  stopifnot( all( dimnames( expression_1_iGiS )[[1]] == dimnames( expression_2_iGiS )[[1]] ) );
  expression_iGiS <- cbind( expression_1_iGiS,
                            expression_2_iGiS );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "G:/Project EMIF/AEdata/Zhu4372/A-GEOD-4372.adf_cut.txt",
                                         columnName = "Reporter.Database.Entry..genbank." );
  
  # Change gene names from RefSeq into Entrez
  library( org.Hs.eg.db );
  entrez_iG <- as.list( org.Hs.egREFSEQ2EG );
  iGs_exist <- names( entrez_iG ) %in% dimnames( expression_e_iGiS )[[1]];
  entrez_iG <- entrez_iG[ iGs_exist ]; 
  entrez_iG <- unlist( entrez_iG );
  expression_te_iGiS <- expression_e_iGiS[ names( entrez_iG ), ];
  dimnames( expression_te_iGiS )[[1]] <- entrez_iG;
  
  # Save stuff
  save( file = "loadAE_zhu_security.Rdata", list = c( "expression_te_iGiS" ) );
  rm( expression_iGiS );
  rm( expression_e_iGiS );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_te_iGiS, 
                                             factorName_sdrf_iF = c( "Comment..Sample_description.", 
                                                                     "Characteristics..age.", 
                                                                     "Characteristics..disease.status.",
                                                                     "Characteristics..sex." ), 
                                             factorName_iF = c( "brainReion", 
                                                                "age", 
                                                                "ad", 
                                                                "gender" ), 
                                             fileSDRF = "G:/Project EMIF/AEdata/Zhu4372/E-GEOD-44772.sdrf.txt",
                                             nameExtension = "2" );
  
  # Put all data together
  stopifnot( all( dimnames( demographics_iSiC )[[1]] %in% dimnames( expression_te_iGiS )[[2]] ) );
  all_iSiC <- as.data.frame( t( expression_te_iGiS ) );
  all_iSiC <- cbind( demographics_iSiC,
                     all_iSiC[ dimnames( demographics_iSiC )[[1]], ] );
  all_iSiC$age <- as.numeric( as.character( all_iSiC$age ) );
  
  # Return
  return( all_iSiC );
  
}

loadAE_berchtold = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Berchtold's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "G:/Project EMIF/AEdata/Berchtold48350/E-GEOD-48350.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "G:/Project EMIF/AEdata/Berchtold48350/A-AFFY-44.adf.csv",
                                         columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Change gene names from RefSeq into Entrez
  library( org.Hs.eg.db );
  entrez_iG <- as.list( org.Hs.egREFSEQ2EG );
  iGs_exist <- names( entrez_iG ) %in% dimnames( expression_e_iGiS )[[1]];
  entrez_iG <- entrez_iG[ iGs_exist ]; 
  entrez_iG <- unlist( entrez_iG );
  expression_te_iGiS <- expression_e_iGiS[ names( entrez_iG ), ];
  dimnames( expression_te_iGiS )[[1]] <- entrez_iG;
    
  # Save stuff
  save( file = "loadAE_berchtold_security.Rdata", list = c( "expression_te_iGiS" ) );
  rm( expression_iGiS );
  rm( expression_e_iGiS );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_te_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.yrs.", 
                                                                     "Characteristics..apoe.genotype.", 
                                                                     "Characteristics..braak.stage.",
                                                                     "Characteristics..brain.region.",
                                                                     "Comment..Sample_title."), 
                                             factorName_iF = c( "age", "apoe", "braak", "brainRegion", "gender" ), 
                                             fileSDRF = "G:/Project EMIF/AEdata/Berchtold48350/E-GEOD-48350.sdrf.txt" );
  demographics_iSiC <- cbind( demographics_iSiC,
                              ad = demographics_iSiC[ ,"braak"] == "  " );
  demographics_iSiC[ ,"gender"] <- sapply( demographics_iSiC[ ,"gender"],
                                           function( s ){
                                             if( grepl( "female",
                                                       s,
                                                       fixed = TRUE ) ){
                                               return( "female" );
                                             }else{
                                               return( "male" );
                                             }
                                           } );
  stopifnot( sum( demographics_iSiC[ ,"ad"] == "TRUE" ) > 100 );
  stopifnot( sum( demographics_iSiC[ ,"gender"] == "female" ) > 100 );
  stopifnot( sum( demographics_iSiC[ ,"gender"] == "male" ) > 100 );
  
  # Put all data together
  stopifnot( all( dimnames( demographics_iSiC )[[1]] %in% dimnames( expression_te_iGiS )[[2]] ) );
  all_iSiC <- as.data.frame( t( expression_te_iGiS ) );
  all_iSiC <- cbind( demographics_iSiC,
                     all_iSiC[ dimnames( demographics_iSiC )[[1]], ] );
  all_iSiC$age <- as.numeric( as.character( all_iSiC$age ) );
  
  # Return
  return( all_iSiC );
  
}

loadAE_luk = function( ){
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Lek's. ",sep="") );
  
  # Get gene expression files
  #expression_iGiS <- read.csv( file = "G:/Project WT-Inf/Data/AE Luk 185/E-TABM-185.processed.csv",
  names_i1iS <- read.csv( file = "G:/Project WT-Inf/Data/AE Luk 185/E-TABM-185.processed.2/E-TABM-185.processed.txt",
                          row.names = 1,
                          nrows = 1,
                          header = FALSE,
                          sep = "\t",
                          stringsAsFactors = FALSE );
  expression_iGiS <- read.csv( file = "G:/Project WT-Inf/Data/AE Luk 185/E-TABM-185.processed.2/E-TABM-185.processed.txt",
                               row.names = 1,
                               skip = 1,
                               sep = "\t" ); #nrows = 10000 );
  stopifnot( dim( names_i1iS )[2] == dim( expression_iGiS )[2] );
  dimnames( expression_iGiS )[[2]] <- as.character( names_i1iS );
  expression_iGiS <- as.matrix( expression_iGiS );
  expression_iGiS <- as.array( expression_iGiS );
  dimnames( expression_iGiS )[[2]] <- sapply( dimnames( expression_iGiS )[[2]], 
                                              function( s ){ 
                                                substr( s, 1, nchar(s)-4 ); 
                                              } );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "G:/Project WT-Inf/Data/AE Luk 185/A-AFFY-33.adf_cut.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_luk_security.Rdata", 
        list = c( "expression_iGiS" ) );
  
  # Return results
  return( expression_iGiS );
  
  
}

loadAE_johannes = function(){
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Ivan's. ",sep="") );
  
  # Get gene expression files
  expression23_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Johannes 13691/E-GEOD-13691.processed.1",
                                            nameColum = "ID_REF",
                                            expColum = "VALUE" );
  expression24_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Johannes 13691/E-GEOD-13691.processed.2",
                                            nameColum = "ID_REF",
                                            expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression23_iGiS <- getGeneNames_ADF( expression_iGiS = expression23_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Johannes 13691/A-AFFY-23.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  expression24_iGiS <- getGeneNames_ADF( expression_iGiS = expression24_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Johannes 13691/A-AFFY-24.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Merge
  iGs <- unique( c( dimnames( expression23_iGiS )[[1]], dimnames( expression24_iGiS )[[1]] ) );
  iG_iG23 <- match( dimnames( expression23_iGiS )[[1]], iGs );
  iG_iG24 <- match( dimnames( expression24_iGiS )[[1]], iGs );
  expression_iGiS <- array( data = NA, dim = c( length( iGs ), dim( expression23_iGiS )[2] + dim( expression24_iGiS )[2] ) );
  expression_iGiS[ iG_iG23, 1:dim( expression23_iGiS )[2] ] <- expression23_iGiS;
  expression_iGiS[ iG_iG24, dim( expression23_iGiS )[2] + ( 1:dim( expression24_iGiS )[2] ) ] <- expression24_iGiS;
  dimnames( expression_iGiS )[[1]] <- iGs;
  dimnames( expression_iGiS )[[2]] <- c( dimnames( expression23_iGiS )[[2]], dimnames( expression24_iGiS )[[2]] );
  
  # Save stuff
  save( file = "loadAE_johannes_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "status" ), 
                                             factorName_iF = c( "status" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Johannes 13691/E-GEOD-13691.sdrf.txt" );
  
  # Tidy demographics up a bit
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "transgenic for Ubb+1" ,"status"] <- "Ubb+1";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "wildtype" ,"status"] <- "wild";
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_ivan = function(){
  # AD11 mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Ivan's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Ivan 63617/E-GEOD-63617.processed.1",
                                           nameColum = "Reporter.Identifier",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Ivan 63617/A-MEXP-724.adf.txt",
                                       columnName = "Reporter.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_hyemyung_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.", 
                                                                     "Characteristics..genotype.", 
                                                                     "Characteristics..organism.part." ), 
                                             factorName_iF = c( "age", "status", "tissue" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Ivan 63617/E-GEOD-63617.sdrf.txt" );
  
  # Tidy demographics up a bit
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "AD11 anti-NGF transgenic mouse" ,"status"] <- "AD11";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "wildtype" ,"status"] <- "wild";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "VH transgenic mouse" ,"status"] <- "VH";
  demographics_iSiC[ demographics_iSiC[ ,"tissue"] == "pool of whole brain from 2 mice" ,"tissue"] <- "brain"; 
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_phillip = function(){
  # Tg2576 mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Phillip's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Phillip 36237/E-GEOD-36237.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Phillip 36237/A-AFFY-45.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_hyemyung_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics.age.", 
                                                                     "Characteristics.variation.", 
                                                                     "Characteristics.organism.part.",
                                                                     "Characteristics.treated.with."), 
                                             factorName_iF = c( "age", "status", "tissue", "extraTreatment" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Phillip 36237/E-GEOD-36237.sdrf.txt" );
  
  # Tidy demographics up a bit
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "Tg+ (over-expressing human mutant amyloid precursor protein)" ,"status"] <- "Tg2576";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "Tg- (non-transgenic WT littermates)" ,"status"] <- "wild";
  cutName <- function(x) { iX_start <- 1;
                           iX_end <- regexpr( " ", x, fixed=TRUE )[1];
                           return( substr( x, iX_start, iX_end-1 ) ); };
  demographics_iSiC[ ,"extraTreatment"] <- sapply( demographics_iSiC[ ,"extraTreatment"], cutName );  
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_giovanni = function(){
  # APP mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Giovanni's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- read.table( file = "F:/Project DKK1/programsR_it1/Data/Giovanni 14499/E-GEOD-14499.processed.1/E-GEOD-14499-processed-data-1721514524.txt",
                                 header = TRUE, 
                                 sep = "\t");
  expression_iGiS <- as.matrix( expression_iGiS );
  expression_iGiS <- as.array( expression_iGiS );
  dimnames( expression_iGiS )[[1]] <- expression_iGiS[ ,1];
  dimnames( expression_iGiS )[[2]] <- sapply( dimnames( expression_iGiS )[[2]], function( s ){ substr( s, 1, nchar(s)-4 ); });
  expression_iGiS <- expression_iGiS[-1,-1];
  expression2_iGiS <- as.numeric( expression_iGiS );
  dim( expression2_iGiS ) <- dim( expression_iGiS );
  dimnames( expression2_iGiS ) <- dimnames( expression_iGiS );
  expression_iGiS <- expression2_iGiS;
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Giovanni 14499/A-AFFY-45.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_giovanni_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  dimnames( expression_iGiS )[[2]] <- paste( "GSE14499", dimnames( expression_iGiS )[[2]], sep = "" );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Comment..Sample_characteristics.", "Comment..Sample_characteristics.", "Comment..Sample_characteristics." ), 
                                             factorName_iF = c( "status", "tissue", "extraTreatment" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Giovanni 14499/E-GEOD-14499.sdrf.healed.txt" );
  
  # Polish the demographics of these guys
  cutName <- function(x) { iX_start <- regexpr( "Genotype: ", x, fixed=TRUE )[1] + 10;
                           iX_end <- regexpr( "; Treatment:", x, fixed=TRUE )[1];
                           return( substr( x, iX_start, iX_end-1 ) ); };
  demographics_iSiC[ ,"status"] <- sapply( demographics_iSiC[ ,"status"], cutName );
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "APP" ,"status"] <- "J20";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "non-transgenic" ,"status"] <- "wild";
  cutName <- function(x) { iX_start <- regexpr( "Region: ", x, fixed=TRUE )[1] + 8;
                           iX_end <- regexpr( "; Genotype:", x, fixed=TRUE )[1];
                           return( substr( x, iX_start, iX_end ) ); };
  demographics_iSiC[ ,"tissue"] <- sapply( demographics_iSiC[ ,"tissue"], cutName );
  cutName <- function(x) { iX_start <- regexpr( "Treatment: ", x, fixed=TRUE )[1] + 11;
                           iX_end <- nchar( x );
                           return( substr( x, iX_start, iX_end ) ); };
  demographics_iSiC[ ,"extraTreatment"] <- sapply( demographics_iSiC[ ,"extraTreatment"], cutName );
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_hyemyung2 = function(){
  # Tg6799 mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Hyemyung's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Hyemung2 52022/E-GEOD-52022.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Hyemung2 52022/A-AFFY-45.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_hyemyung2_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.", "Characteristics..genotype.", "Characteristics..organism.part." ), 
                                             factorName_iF = c( "age", "status", "tissue" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Hyemung2 52022/E-GEOD-52022.sdrf.txt" );
  
  # Fix demographics
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "AD mutant" ,"status"] <- "5xFAD";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "littermate" ,"status"] <- "wild";
  demographics_iSiC[ demographics_iSiC[ ,"tissue"] == "Brain Hippocampus" ,"tissue"] <- "hippocampus";
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_hyemyung = function(){
  # Tg6799 mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Hyemyung's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Hyemyung 52024/E-GEOD-52024.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Hyemyung 52024/A-AFFY-45.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_hyemyung_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.", "Characteristics..genotype.", "Characteristics..organism.part." ), 
                                             factorName_iF = c( "age", "status", "tissue" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Hyemyung 52024/E-GEOD-52024.hyb.sdrf.txt" );
  
  # Fix demographics
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "AD mutant" ,"status"] <- "5xFAD";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "littermate" ,"status"] <- "wild";
  demographics_iSiC[ demographics_iSiC[ ,"tissue"] == "Brain Hippocampus" ,"tissue"] <- "hippocampus";
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyemyung" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_hyslop = function(){
  # TgCRND8  mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Becker's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Hyslop 31372/E-GEOD-31372.processed.1",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Hyslop 31372/A-AFFY-45.adf.txt",
                                       columnName = "Composite.Element.Database.Entry.refseq." );
  
  # Save stuff
  save( file = "loadAE_hyslop_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.", "Characteristics..genetic.background.", "Characteristics..organism.part." ), 
                                             factorName_iF = c( "age", "status", "tissue" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Hyslop 31372/E-GEOD-31372.sdrf_healed.txt" );
  
  # Fix demographics
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "TgCRND8 transgenic mouse" ,"status"] <- "TgCRND8";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "non-transgenic littermate mouse" ,"status"] <- "wild";
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_hyslop" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

loadAE_becker = function(){
  # 3xTgAD mice model
  
  # Some libs
  require( R.utils );
  require( data.table );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress Becker's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "F:/Project DKK1/programsR_it1/Data/Becker 60911/E-GEOD-60911.processed.1",
                                           nameColum = "Reporter.Identifier",
                                           expColum = "VALUE" );
  
  # Change the dimention names into Entrez names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                       fileADF = "F:/Project DKK1/programsR_it1/Data/Becker 60911/A-MEXP-1174.adf_cut.txt",
                                       columnName = "Reporter.Database.Entry.refseq." );
                                       #columnName = "Reporter.Database.Entry.hugo." );
  
  # Save stuff
  save( file = "loadAE_becker_security.Rdata", list = c( "expression_iGiS" ) );
  
  # Get demographics
  print( paste(match.call()[1]," = Loading demographics. ",sep="") );
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..genotype.", "Characteristics..organism.part.", "Characteristics..genotype."), 
                                             factorName_iF = c( "status", "tissue", "extraTreatment" ), 
                                             fileSDRF = "F:/Project DKK1/programsR_it1/Data/Becker 60911/E-GEOD-60911.sdrf.txt" );
  # Fix demographics
  demographics_iSiC[ ,"status"] <- sapply( demographics_iSiC[ ,"status"], function( s ) { substr( s, 1, 4 );} );
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "3XTG" ,"status"] <- "3xTg";
  demographics_iSiC[ demographics_iSiC[ ,"status"] == "POLB" ,"status"] <- "wild";
  demographics_iSiC[ ,"extraTreatment"] <- sapply( demographics_iSiC[ ,"extraTreatment"], function( s ) { substr( s, nchar(s)-6, nchar(s) );} );
  demographics_iSiC[ demographics_iSiC[ ,"extraTreatment"] == "ld type" ,"extraTreatment"] <- "none";  
  
  # Put everything into a DataVST structure
  data <- DataVST$new();
  data$meaningDimV_c <- meaningDimV$GENE_EXPRESSION;
  data$history_jHiC <- list( "loadAE_becker" );
  data$setValue( expression_iGiS );
  data$demographics_iSiC <- demographics_iSiC;
  
  # Return
  return( data );
  
}

getRawFiles = function( folder ){
  
  require(affy);
  require(mogene11stv1cdf)
  
  # Load file names
  files_iF <- dir(folder);
  files_iF <- sapply( files_iF , function (x) { paste( folder, "/", x , sep="" ) } );
  
  a <- ReadAffy( filenames = files_iF[[1]] );
  b <- rma( a );
  
  # Load gene expression from files
  print( paste( match.call()[1], " = Loading expression file by file. This takes a while.", sep="") );
  for( iF in 1:length(files_iF) ){
    
    # Hale the user
    if( mod( iF, 100 ) == 1 ){
      print( paste( match.call()[1], " = Loading file ", iF, " of ", length(files_iF), sep="") );
    }
    
    # Load this subject
    expression_iGiC <- read.table( file = files_iF[iF],
                                   header = TRUE, 
                                   sep = "\t");
    stopifnot( dimnames( expression_iGiC )[[2]][1] == nameColum );
    stopifnot( dimnames( expression_iGiC )[[2]][2] == expColum );
    
    # Add this subject
    if( iF == 1 ){
      expression_iGiS <- array( data = NA, dim = c( dim( expression_iGiC )[1], length(files_iF) ) );
      dimnames( expression_iGiS )[[1]] <- expression_iGiC[ ,nameColum]
      dimnames( expression_iGiS )[[2]] <- sapply( files_iF, function(s) { 
        iS_start <- regexpi( s, "/GSM")$start[1];
        iS_end <- regexpi( s, "_sample")$start[1];
        return( substr( s, iS_start+1, iS_end-1 ) );
      });
    }
    expression_iGiS[ ,iF] <- expression_iGiC[ ,expColum];
    stopifnot( all( dimnames( expression_iGiS )[[1]] == expression_iGiC[ ,nameColum] ) );
    
  }
  
  # Final check
  print( paste( match.call()[1], " = Running final data check.", sep="") );
  stopifnot( sum( is.na( c( expression_iGiS ) ) ) == 0 );
  
  # Return
  return( expression_iGiS );
  
}

# Given a folder with all the subject files (usual method used in ArrayExpress), this function loads the expression
# of all those files and store them in a single matrix
# - folder = Folder with all the files
# - nameColum = Name of the colum with the probe names
# - expColum = name of the columns with expression values
getPreprocessedFiles = function( folder, nameColum, expColum ){
  
  # Load file names
  files_iF <- dir(folder);
  files_iF <- sapply( files_iF , function (x) { paste( folder, "/", x , sep="" ) } );
  
  # Load gene expression from files
  print( paste( match.call()[1], " = Loading expression file by file. This takes a while.", sep="") );
  for( iF in 1:length(files_iF) ){
    
    # Hale the user
    if( mod( iF, 100 ) == 1 ){
      print( paste( match.call()[1], " = Loading file ", iF, " of ", length(files_iF), sep="") );
    }
    
    # Load this subject
    expression_iGiC <- read.table( file = files_iF[iF],
                                   header = TRUE, 
                                   sep = "\t");
    stopifnot( dimnames( expression_iGiC )[[2]][1] == nameColum );
    stopifnot( dimnames( expression_iGiC )[[2]][2] == expColum );
    
    # Add this subject
    if( iF == 1 ){
      expression_iGiS <- array( data = NA, dim = c( dim( expression_iGiC )[1], length(files_iF) ) );
      dimnames( expression_iGiS )[[1]] <- expression_iGiC[ ,nameColum]
      dimnames( expression_iGiS )[[2]] <- sapply( files_iF, function(s) { 
        iS_start <- regexpi( s, "/GSM")$start[1];
        iS_end <- regexpi( s, "_sample")$start[1];
        return( substr( s, iS_start+1, iS_end-1 ) );
      });
    }
    expression_iGiS[ ,iF] <- expression_iGiC[ ,expColum];
    stopifnot( all( dimnames( expression_iGiS )[[1]] == expression_iGiC[ ,nameColum] ) );
    
  }
  
#   # Final check
#   print( paste( match.call()[1], " = Running final data check.", sep="") );
#   stopifnot( sum( is.na( c( expression_iGiS ) ) ) == 0 );
  
  # Return
  return( expression_iGiS );
  
}

# Given a matrix of gene expression "expression_iGiS", whose dimnames()[[2]] are names of subjects, 
# this function returns a number of demographic factors for those subjects as recorded in the corresponding
# SDRF file
# - expression_iGiS = matrix of expression
# - fileSDRF = name of the SDRF file with demographics informacion
# - factorName_sdrf_iF = names of the demographics that you want to collect, as recorded in "fileSDRF"
# - factorName_iF = more friendly names for the same demographics, which will be used as dimnames for the returned matrix
getDemographics_SDRF = function( expression_iGiS, 
                                 factorName_sdrf_iF, 
                                 factorName_iF, 
                                 fileSDRF,
                                 nameExtension = "1" ){
  
  # Eliminate useless lines
  # CUE = If there are weird characters, read.table stops reading there
  demographics_iRiC <- read.table( file = fileSDRF,
                                   header = TRUE,
                                   fill = TRUE,
                                   sep = ",",
                                   comment.char = "",
                                   quote = "");
  
  # Find index correspondence
  
  print(dimnames(demographics_iRiC)[[2]])
  subjectID_iR <- as.character( demographics_iRiC[ ,"Source.Name"] );
  subjectID_iS <- dimnames( expression_iGiS )[[2]];
  subjectID_iS <- paste( subjectID_iS,
                         nameExtension,
                         sep = " " )
  iR_iS <- match( subjectID_iS, subjectID_iR );

  stopifnot( length( iR_iS ) == length( unique( iR_iS ) ) ); # stop("At least one subject is represented more than 1");
  
  # Get interesting demographics
  demographics_iSiC <- array( data = NA, dim = c( dim( expression_iGiS )[2], length( factorName_sdrf_iF ) ) );
  for( iF in 1:length( factorName_sdrf_iF ) ){
    stopifnot( any( factorName_sdrf_iF[iF] == dimnames( demographics_iRiC )[[2]] ) ); #This colum does not exist
    factor_iS <- as.character( demographics_iRiC[iR_iS,factorName_sdrf_iF[iF]] );
    demographics_iSiC[ ,iF] <- factor_iS;
  }
  dimnames( demographics_iSiC )[[1]] <- dimnames( expression_iGiS )[[2]];
  dimnames( demographics_iSiC )[[2]] <- factorName_iF;
  
  #Return
  stopifnot( dim( demographics_iSiC )[1] == dim( expression_iGiS )[2] );
  return( demographics_iSiC );
}

# Given a matrix of gene expression "expression_iGiS", whose dimnames()[[1]] are names of affymetrix probes, 
# this function transform these probe names into ensembl gene names
# - expression_iGiS = matrix of expression
# - fileADF = name of the ADF file with information about the affymetrix probes
getGeneNames_ADF = function( expression_iGiS, fileADF, columnName = "Composite.Element.Database.Entry.entrez." ){
  
  # Load the affymetrix probe IDs
  reporterID_iG <- dimnames( expression_iGiS )[[1]];
  geneIDs_iRiC <- read.table( file = fileADF,
                              header = TRUE,
                              fill = TRUE,
                              sep = "\t");
  iR_iG <- match( reporterID_iG, geneIDs_iRiC[ ,1] );
  entrezID_iG <- geneIDs_iRiC[iR_iG,columnName];
  # Eliminate genes with no entrez name
  stopifnot( sum( is.na(entrezID_iG) ) < round( 0.5 * length( entrezID_iG ) ) ); #More than half the probes have no entrez name
  expression_iGiS <- expression_iGiS[!is.na( entrezID_iG ),];
  entrezID_iG <- entrezID_iG[!is.na( entrezID_iG )];
  
  # Eliminate repeated names
  entrezID_u_iG <- unique( entrezID_iG );
  iGs_u <- match( entrezID_u_iG, entrezID_iG );
  expression_iGiS <- expression_iGiS[iGs_u, ];
  # Set new dimnames
  dimnames( expression_iGiS )[[1]] <- entrezID_u_iG;
  
#   # Duplicate reporters with 2 Entrez names. VERY IMPORTANT! Do so in the other cohorts
#   iSnum <- dim( expression_iGiS )[2];
#   entrezID_u_iG <- unlist(entrezID_iG);
#   entrezID_u2_iG <- entrezID_u_iG;
#   expression2_iGiS <- matrix( data= NA ,
#                               nrow= length(entrezID_u_iG) ,
#                               ncol= iSnum );
#   dimnames(expression2_iGiS)[[1]] <- entrezID_u_iG;
#   dimnames(expression2_iGiS)[[2]] <- dimnames(expression_iGiS)[[2]];
#   bEmpty_iG <- array(TRUE,c(length(entrezID_u_iG),1));
#   for( iG in 1:dim(expression_iGiS)[1] ){
#     
#     # Duplicate data for all existing entrez names, getting care not to erase 2 entrezs that belong to different ensembls
#     iGs <- match(entrezID_iG[[iG]],entrezID_u2_iG);
#     entrezID_u2_iG[iGs] <- "assigned";
#     duplicated_iSiG <- replicate(length(iGs),expression_iGiS[iG,]);
#     expression2_iGiS[iGs,] <- t( duplicated_iSiG );
#     if( dim(duplicated_iSiG)[1] != dim(expression2_iGiS)[2] ||
#           dim(duplicated_iSiG)[2] != length(iGs)  ){
#       stop("Wrong dims");
#     }
#     
#   }
#   expression_iGiS <- expression2_iGiS;
#   rm( expression2_iGiS );
  
  # Return result
  return( expression_iGiS );
  
}






