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

loadAE_63063 = function() {
  # Import libraries
  
  require( R.utils );
  require( data.table );
  require( pracma );
  
  expression_iGiS <- getPreprocessedFiles( folder = "C:/Users/cwong/masters-thesis/data/E-GEOD-63063/samples/1/",
                                           nameColum = "Reporter.Identifier",
                                           expColum = "VALUE" );
  
  expression_e_iGiS <- getGeneNames_ADF(expression_iGiS = expression_iGiS,
                                        fileADF = "C:/Users/cwong/masters-thesis/data/E-GEOD-63063/A-GEOD-10558.adf.txt",
                                        columnName = "Reporter.Database.Entry..genbank.");
  
  library( org.Hs.eg.db );
  entrez_iG <- as.list(org.Hs.egACCNUM2EG);
  print(entrez_iG[1:10])
  iGs_exist <- names( entrez_iG ) %in% dimnames( expression_e_iGiS )[[1]];
  entrez_iG <- entrez_iG[ iGs_exist ]; 
  entrez_iG <- unlist( entrez_iG );


  expression_te_iGiS <- expression_e_iGiS[ names( entrez_iG ), ];

  dimnames( expression_te_iGiS )[[1]] <- entrez_iG;
  print(dim(expression_te_iG))
  print(expression_te_iGiS[1:10, 1:10])
  save( file = "loadAE_E-GEOD-48350_security.Rdata", list = c( "expression_te_iGiS" ) );
  rm( expression_iGiS );
  rm( expression_e_iGiS );
  
  demographics_iSiC <- getDemographics_SDRF( expression_iGiS = expression_te_iGiS, 
                                             factorName_sdrf_iF = c( "Characteristics..age.", 
                                                                     "Characteristics..sex."),
                                             factorName_iF = c( "age", "gender" ), 
                                             fileSDRF = "C:/Users/cwong/masters-thesis/data/E-GEOD-63063/E-GEOD-63063.sdrf.txt" );
  demographics_iSiC[ ,"age"] <- sapply( demographics_iSiC[ ,"age"],
                                        function( s ){
                                          return( strsplit( s, " ")[[1]][1] );
                                        } );
  demographics_iSiC[ demographics_iSiC[ ,"age"] == ">90", "age" ] <- 90;
  
  stopifnot( all( dimnames( demographics_iSiC )[[1]] %in% dimnames( expression_te_iGiS )[[2]] ) );
  all_iSiC <- as.data.frame( t( expression_te_iGiS ) );
  all_iSiC <- cbind( demographics_iSiC,
                     all_iSiC[ dimnames( demographics_iSiC )[[1]], ] );
  all_iSiC$age <- as.numeric( as.character( all_iSiC$age ) );
  
  write.csv(all_iSiC, file="C:/Users/cwong/masters-thesis/data/E-GEOD-63063/E-GEOD-63063-combined.csv")
}


pretty_print_progress = function(text, curr, total) {
  cat(paste("\r", text, ": ", curr, "/", total, sep=""));
  if (curr == total) {
    cat("\n");
  }
}

# Loads data from ArrayExpress
#   Preconditions
#     - Folder structure:
#         folder/
#           |----samples/
#           |      |----source_1_sample_table.txt
#           |      |----source_2_sample_table.txt
#           |     ...
#           |      |----source_n_sample_table.txt
#           |----dataset.sdrf.txt
#           |----chip.adf.txt
#
#   Parameters
#     folder     - top level folder (e.g. "C:/Users/user/data/E-GEOD-884422")
#     dataset    - name of dataset (e.g. "E-GEOD-884422")
#     funcs      - vector of functions to parse metadata from sdrf file
#     adf        - filename of ADF file
#     adf_sep    - delimiter in adf file
#     sample_sep - delimiter in sample files
#     sdrf_sep   - delimiter in sdrf file

load_data = function(folder,
                     dataset,
                     funcs,
                     adf,
                     adf_sep=",",
                     sample_sep=",",
                     sdrf_sep=",") {
  require(R.utils);
  require(data.table);
  require(pracma);
  
  print(paste("Combining data from ", folder, sep=""));
  samples_dir <- paste(folder, "samples", sep="/");
  samples_files <- list.files(samples_dir);
  print(paste("Number of samples detected: ", length(samples_files), sep=""))
  
  print("Reading ADF file")
  adf <- paste(folder, adf, sep="/");
  adf_data <- read.table(file = adf,
                         header = TRUE,
                         sep = adf_sep,
                         fill = TRUE,
                         na.strings = "");

  library( org.Hs.eg.db );
  entrez_table <- as.list( org.Hs.egREFSEQ2EG );
  
  combined <- NULL;
  samples_count <- length(samples_files);
  
  for (i in 1:samples_count) {
    pretty_print_progress("Processing sample files", i, samples_count);
    data <- read.table(file = paste(samples_dir, samples_files[i], sep="/"),
                       header = TRUE, 
                       sep = sample_sep,
                       na.strings = "");
    probe_ids <- data[, 1];is
    
    # Match sample probes with Refseq IDs
    probe_matches <- match(probe_ids, adf_data[, 1]);

    
    # Strip all non-existing refseq rows
    probe_entrez <- adf_data[probe_matches, "Composite.Element.Database.Entry.refseq."];
    data <- data[!is.na(probe_entrez), ];
    probe_entrez <- probe_entrez[!is.na(probe_entrez)];
  
    
    # Unique Refseq IDs
    probe_entrez_u <- unique(probe_entrez);
    u <- match(probe_entrez_u, probe_entrez);
    data <- data[u, 2, drop=FALSE];
    
    dimnames(data)[[1]] <- probe_entrez_u;
    dimnames(data)[[2]] <- c(samples_files[i]);
    
    if (dim(data)[[1]] > 0) {
      if (is.null(combined)) {
        combined <- data;
      } else {
        d_a <- dimnames(combined)[[1]];
        d_b <- dimnames(data)[[1]];
        m_a <- match(d_a, d_b);
        m_b <- match(d_b, d_a);
        combined <- combined[!is.na(m_a), , drop=FALSE];
        data <- data[!is.na(m_b), , drop=FALSE];
        combined[,samples_files[i]] <- data[, 1];
      }
    }
  }
  
  
  
  print("Converting gene identifiers into Entrez IDs");
  library( org.Hs.eg.db );
  entrez <- as.list( org.Hs.egREFSEQ2EG );
  entrez_exist <- names( entrez ) %in% dimnames( combined )[[1]];
  entrez <- entrez[ entrez_exist ]; 
  entrez <- unlist( entrez );
  combined <- combined[names(entrez), ];
  entrez_u <- unique(entrez);
  u <- match(entrez_u, entrez);
  combined <- combined[u, ];

  dimnames(combined)[[1]] <- entrez_u;

  
  # Get patient data
  sdrf_fp <- paste(folder, "/", dataset, ".sdrf.txt", sep="");
  sdrf_data <- read.csv(file = sdrf_fp,
                          header = TRUE,
                          sep = sdrf_sep,
                          fill = TRUE,
                          na.strings = "");
  dimnames(sdrf_data)[[1]] = sdrf_data[, 1];

  c_dim <- dim(combined)[[1]];
  
  for (i in 1:dim(combined)[2]) {
    
    pretty_print_progress("Getting metadata from SDRF file", i, samples_count);
    metadata <- sdrf_data[paste(strsplit(dimnames(combined)[[2]][i], "_")[[1]][1], " 1", sep=""), ];
    for (j in 1:length(funcs)) {
      combined[c_dim + j, i] <- funcs[[j]](metadata);
    }
  }
  data <- t(combined[,1:ncol(combined)]);
  rm(combined);
  return(data);
  
}

test = function() {
  if (!file.exists("C:/Users/cwong/masters-thesis/data/E-GEOD-84422/E-GEOD-84422-combined-data.csv")) {
    #84422
    funcs_84422 = c(
      function(x) {   # Gender
        if (as.character(x[["Term.Accession.Number"]]) == "female") {
          return(0);
        } 
        return(1);
      },
      function(x) {   # AGE
        return(strsplit(as.character(x[["Comment..Sample_title."]]), " ")[[1]][2]);
      },
      function(x) {   # AD
        return(switch(as.character(x[["ad"]]), normal=0, 1));
      }
    )
    data_84422 = load_data("C:/Users/cwong/masters-thesis/data/E-GEOD-84422",
                     "E-GEOD-84422",
                     funcs_84422,
                     "A-AFFY-33.adf.txt",
                     adf_sep=",",
                     sample_sep="\t",
                     sdrf_sep=",")
    write.csv(data_84422, file="C:/Users/cwong/masters-thesis/data/E-GEOD-84422/E-GEOD-84422-combined-data.csv");
  }
  
  
  if (!file.exists("C:/Users/cwong/masters-thesis/data/E-GEOD-48350/E-GEOD-48350-combined-data.csv")) {
    # 48350
    funcs_48350 = c(
      function(x) {   # Gender
        if (as.character(x[["Characteristics..sex."]]) == "female") {
          return(0);
        } 
        return(1);
      },
      function(x) {   # AGE
        return(strsplit(as.character(x[["Characteristics..age.yrs."]]), " ")[[1]][2]);
      },
      function(x) {   # AD
        if (x[["Characteristics..braak.stage."]] > 0) {
          return(1);
        }
        return(0);
      }
    )
    data_48350= load_data("C:/Users/cwong/masters-thesis/data/E-GEOD-48350",
                           "E-GEOD-48350",
                           funcs_48350,
                           "A-AFFY-44.adf.txt",
                           adf_sep="\t",
                           sample_sep="\t",
                           sdrf_sep=",")
    write.csv(data_48350, file="C:/Users/cwong/masters-thesis/data/E-GEOD-48350/E-GEOD-48350-combined-data.csv");
  }
}

loadAE_84422 = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress GEOD84422's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "C:/Users/cwong/masters-thesis/data/E-GEOD-84422/samples",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "C:/Users/cwong/masters-thesis/data/E-GEOD-84422/A-AFFY-33.adf.txt",
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
                                                                     "Comment..Sample.Table.", 
                                                                     "ad"), 
                                             factorName_iF = c( "brainRegion", "age", "gender", "ad" ), 
                                             fileSDRF = "C:/Users/cwong/masters-thesis/data/E-GEOD-48350/E-GEOD-48350.sdrf.txt" );
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


loadAE_48350 = function(){
  # Postmortem brain in AD. 250 people
  
  # Some libs
  require( R.utils );
  require( data.table );
  require( pracma );
  
  # Hale the user
  print( paste(match.call()[1]," = Loading data from ArrayExpress GEOD48350's. ",sep="") );
  
  # Get gene expression files
  expression_iGiS <- getPreprocessedFiles( folder = "C:/Users/alejon/Documents/samuel/masters/data/E-GEOD-84422/samples",
                                           nameColum = "ID_REF",
                                           expColum = "VALUE" );
  
  # Change the dimention names into RefSeq names
  print( paste(match.call()[1]," = Loading gene names.  ",sep="") );
  expression_e_iGiS <- getGeneNames_ADF( expression_iGiS = expression_iGiS, 
                                         fileADF = "C:/Users/alejon/Documents/samuel/masters/data/E-GEOD-84422/A-AFFY-44.adf.txt",
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
                                   sep = "\t",
                                   comment.char = "",
                                   quote = "");
  # Find index correspondence
  
  subjectID_iR <- as.character( demographics_iRiC[ ,"Source.Name"] );

  print(dimnames(demographics_iRiC)[[2]]);
  subjectID_iS <- dimnames( expression_iGiS )[[2]];
  subjectID_iS <- paste( subjectID_iS,
                         nameExtension,
                         sep = " " )
  iR_iS <- match( subjectID_iS, subjectID_iR );

  #stopifnot( length( iR_iS ) == length( unique( iR_iS ) ) ); # stop("At least one subject is represented more than 1");
  
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
                              sep = ",");
  iR_iG <- match( reporterID_iG, geneIDs_iRiC[ , 1] );

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






