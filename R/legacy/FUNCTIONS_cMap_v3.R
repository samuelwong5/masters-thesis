#
# This program builds the drug signatures by using C-map data
# INPUT = Drug names (drugName_iD) and cell lines (drugCell_iS). Use the names indicated in "cmap_instances.csv".
# OUTPUT = The signatures, saved in "drugSignatures_v3.Rdata"
#
# v2 = This is actually the first version. Modified from "FUNCTIONS_cMap.R" from "WT/it6"
# v3 = Fixing some stupind bug in "loadExpression()"
#

mainDir <- "G:\\Project - Imelda WT-Inf\\";

# Find files
findFiles = function( drugName_iD ){
  # Given a list of drug names, finds all the Cmap files related to each drug
  #
  # Args:
  #   drugName_iD: Names of the drugs
  #
  # Returns:
  #   file_iEiC: A data frame with the names of the files and some other features of each file
  #
  
  # Find the names of the files corresponding to each experiment
  experimentInfo_iEiC    <- read.csv( paste( mainDir, "Data/CmapFull/cmap_instances.csv", sep = "" ),
                                      header = TRUE,
                                      stringsAsFactors = FALSE);
  # bDrug_iE <- !grepl("^[[:digit:]]",experimentInfo_iEiC[ ,"cmap_name"]) & experimentInfo_iEiC[ ,"cmap_name"] != "";
  # drugName_iD <- unique( experimentInfo_iEiC[bDrug_iE,"cmap_name"] );
  bDrugWithoutFile_iD <- rep( FALSE, length( drugName_iD ) );
  file_iE <- c();
  drug_iE <- c();
  vehicle_iE <- c();
  concentration_iE <- c();
  scanner_iE <- c();
  vendor_iE <- c();
  batch_iE <- c();
  cell_iE <- c();
  for( iD in 1:length(drugName_iD) ){
    
    # Get experiments
    print( paste( "Loading C-map data for ", drugName_iD[iD], sep="" ) );
    iEs <- which( ( grepl( tolower( drugName_iD[iD] ), 
                           experimentInfo_iEiC[,"cmap_name"], 
                           fixed = TRUE ) ) &
                    ( experimentInfo_iEiC[,"array3"] == "HT_HG-U133A" ) );
    print( paste( "Found ", as.character( length( iEs ) ), " experiments." ) );
    if( length(iEs)==0 ){
      #This drug hasn't got experiments that we can use. Skip and go to the next one
      bDrugWithoutFile_iD[ iD ] <- TRUE;
      next;
    }
    
    # Add drug files
    file_iE <- c( file_iE, experimentInfo_iEiC[iEs,"perturbation_scan_id"] );
    drug_iE <- c( drug_iE, rep( drugName_iD[iD], length(iEs) ) );
    vehicle_iE <- c( vehicle_iE, rep( FALSE, length(iEs) ) );
    concentration_iE <- c( concentration_iE, experimentInfo_iEiC[iEs,"concentration..M."] );
    scanner_iE <- c( scanner_iE, experimentInfo_iEiC[iEs,"scanner"] );
    vendor_iE <- c( vendor_iE, experimentInfo_iEiC[iEs,"vendor"] );
    batch_iE <- c( batch_iE, experimentInfo_iEiC[iEs,"batch_id"] );
    cell_iE <- c( cell_iE, experimentInfo_iEiC[iEs,"cell2"] );
    
    # Add vehicle files
    for( iE in iEs ){
      
      # Get prefix
      dots <- gregexpr( ".", experimentInfo_iEiC[iE,"perturbation_scan_id"], fixed=TRUE );
      if( length(dots) != 1 ){
        stop("A file name has more than 1 dot");
      }
      filePrefix <- experimentInfo_iEiC[iE,"perturbation_scan_id"];
      filePrefix <- strsplit( filePrefix, ".", fixed = TRUE )[[1]][1];
      
      # Get sufix
      fileSufix_iE <- strsplit( experimentInfo_iEiC[iE,"vehicle_scan_id4"], ".", fixed = TRUE )[[1]];
      fileSufix_iE <- fileSufix_iE[ 2:length(fileSufix_iE) ];
      if( length( fileSufix_iE ) == 0 | any( sapply( fileSufix_iE, length ) == 0 ) ){
        stop("Sothing is going wrong with the names of the vehicle files");
      }
      
      # Get all file atributes
      file_iE <- c( file_iE, sapply( fileSufix_iE, function(suffix) { paste( filePrefix, ".", suffix , sep = ""); } ) );
      drug_iE <- c( drug_iE, rep( drugName_iD[iD], length( fileSufix_iE ) ) );
      vehicle_iE <- c( vehicle_iE, rep( TRUE, length( fileSufix_iE ) ) );
      concentration_iE <- c( concentration_iE, rep( experimentInfo_iEiC[iE,"concentration..M."], length( fileSufix_iE ) ) );
      scanner_iE <- c( scanner_iE, rep( experimentInfo_iEiC[iE,"scanner"], length( fileSufix_iE ) ) );
      vendor_iE <- c( vendor_iE, rep( experimentInfo_iEiC[iE,"vendor"], length( fileSufix_iE ) ) );
      batch_iE <- c( batch_iE, rep( experimentInfo_iEiC[iE,"batch_id"], length( fileSufix_iE ) ) );
      cell_iE <- c( cell_iE, rep( experimentInfo_iEiC[iE,"cell2"], length( fileSufix_iE ) ) );
      
    }
  }
  
  # Correct file names
  file_iE <- sapply( file_iE, function(x) { sub( "'", "", x, fixed=TRUE ) } );
  
  # Check file names are OK
  if( length( file_iE ) != length( drug_iE ) ||
        length( file_iE ) != length( vehicle_iE ) ||
        length( file_iE ) != length( concentration_iE ) ||
        length( file_iE ) != length( scanner_iE ) ||
        length( file_iE ) != length( vendor_iE ) ||
        length( file_iE ) != length( batch_iE ) ||
        length( file_iE ) != length( cell_iE ) ){
    stop("The number of elements in some of the file variables do not match the rest.")
  }
  
  # Merge all
  file_iEiC <- data.frame( file = file_iE,
                           drug = drug_iE,
                           vehicle = vehicle_iE,
                           concentration = concentration_iE,
                           scanner = scanner_iE,
                           vendor = vendor_iE,
                           batch = batch_iE,
                           cell = cell_iE,
                           stringsAsFactors = FALSE );
  return( file_iEiC )
  
}

loadExpression = function( file_iE ){
  # Loads all the drug files given in "file_iE"
  #
  # Args:
  #   file_iE: Names of the drug files, such as "5500024030760072207028.H09", "5500024030760072207028.B1"...
  #
  # Returns:
  #   exp_iEiG: A data frame with the expression level of each file "iE" for each gene "iG"
  #
  
  # Dependencies
  library( affy );
  
  # Find the folder of each file
  folders_iF <- c( "cmap_build02.volume1of7", "cmap_build02.volume2of7",
                   "cmap_build02.volume3of7", "cmap_build02.volume4of7",
                   "cmap_build02.volume5of7", "cmap_build02.volume6of7",
                   "cmap_build02.volume7of7" );
  fileFolder_iE <- rep( NA, length( file_iE ) );
  for( iF in 1:length( folders_iF ) ){
    # See which files exist in folder iF
    cmapFiles_iF <- dir( paste( mainDir, "Data/CmapFull/", folders_iF[iF], sep = "" ) );
    exists_iE <- sapply( file_iE, 
                         function( s ) any( grepl( s, cmapFiles_iF, fixed = TRUE ) ) );
    stopifnot( all( is.na( fileFolder_iE[ exists_iE ] ) ) ); #some file has been found in more than one folder
    fileFolder_iE[ exists_iE ] <- iF;
  } 
  stopifnot( all( !is.na( fileFolder_iE ) ) ); #some file wasn't found in any folder
  
  # Load each file
  exp_iEiG <- array( data=NA, 
                     dim = c( length( file_iE ), 22277 ),
                     dimnames = list( file = file_iE,
                                      gen = NULL ) );
  geneNames_iG <- NA;
  for( iE in 1:length( file_iE ) ){
    
    # Hale the user
    if( iE %% 100 == 1 ){
      print( paste( iE, " of ", length( file_iE ) , " loaded. ", sep = "" ) );
    }
    
    # Load file
    print( paste( "------------- Loading experiment iE ", as.character(iE)," of ", length(file_iE), sep=""));
    affydata <- ReadAffy( filenames = paste( mainDir, 
                                             "Data/CmapFull/", 
                                             folders_iF[ fileFolder_iE[iE] ], 
                                             "/",
                                             file_iE[iE], 
                                             ".CEL", 
                                             sep="" ) );
    nomalisedVals <- affy::rma( affydata );
    exp_iEiG[iE, ] <- exprs( nomalisedVals );
    
    # Check
    if( all( is.na( geneNames_iG ) ) ){
      geneNames_iG <- featureNames( affydata );
    }else if( any( geneNames_iG != featureNames(affydata) ) ){
      stop("Some of the expreiments have different features than the others.");
    }
    
  }
  
  # Load Affimetrix annotations
  annotations_iRiC    <- read.csv( paste( mainDir, 
                                          "Data/CmapFull/HT_HG-U133A.na34.annot.csv",
                                          sep = "" ),
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   comment.char = "#" );
  dimnames( annotations_iRiC )[[1]] <- annotations_iRiC[,1];
  
  # Get gene names
  #geneNames_iG <- annotations_iRiC[ geneNames_iG, "Gene.Symbol" ];
  geneNames_1_iG <- annotations_iRiC[ , "Entrez.Gene" ];
  geneNames_2_iG <- sapply( geneNames_1_iG, 
                            function(x) {
                              a <- strsplit( sub( " ",
                                                  "",
                                                  x ),
                                             "///" ); 
                              return( a[[1]][1] );
                            } );
  names( geneNames_2_iG ) <- NULL;
  dimnames( exp_iEiG )[[2]] <- geneNames_2_iG;
  exp_iEiG <- exp_iEiG[ , geneNames_2_iG != "---" ];
  
  #Return results
  return( exp_iEiG );
  
}

testDrugs = function( exp_iEiG,
                      file_iEiC ){
  # For all genes "iG" and all drugs "iD", it tests whether the drugs changes the expression of the gen
  #
  # Args:
  #   exp_iEiG: Expression of each gen "iG" from each file "iE". This matrix is built by function "loadExpression"
  #   file_iEiC: Matrix with information of each one of the files. This matrix is built by function "findFiles"
  #
  # Returns:
  #   pVals_iDiG: P values for each drug "iD" and each gen "iG"
  #
  
  # Check input
  stopifnot( all( file_iEiC[ ,"file"] %in% dimnames( exp_iEiG )[[1]] ) );
  stopifnot( all( file_iEiC[ ,"file"] == dimnames( exp_iEiG )[[1]] ) );
  
  # Test all drugs available in Cmap
  mainDir <- "G:/Project WT-Inf/";
  drugName_iD <- unique( file_iEiC[ ,"drug"] );
  cell_iL <- unique( file_iEiC[ ,"cell"] );
  pVals_iDiG <- array( data = NA, 
                       dim = c( length( drugName_iD ), 
                                dim( exp_iEiG )[2] ),
                       dimnames = list( drug = drugName_iD,
                                        gen = dimnames( exp_iEiG )[[2]] ) );
  effec_iDiG <- pVals_iDiG;
  numSamples_iD <- array( data = NA, 
                          dim = c( length( drugName_iD ), 
                                   1 ) );
  numscanners_iD <- array( data = NA, 
                           dim = c( length( drugName_iD ), 
                                    1 ) );
  cellLines_iDiC <- array( data = FALSE, 
                           dim = c( length( drugName_iD ), 
                                    length( cell_iL ) ),
                           dimnames = list( drug = drugName_iD,
                                            cell = cell_iL ) );
  for( iD in 1:length( drugName_iD ) ){
    
    # Create formula
    iEs <- which( file_iEiC[ ,"drug"] == drugName_iD[iD] );
    formulat <- 'exp_iE2 ~ file_iEiC[iEs,"vehicle"]';
    if( length( unique( file_iEiC[iEs,"concentration"] ) ) > 1 ){
      formulat <- paste( formulat, ' + file_iEiC[iEs,"concentration"]', sep = "" );
    }
    if( length( unique( file_iEiC[iEs,"scanner"] ) ) > 1 ){
      formulat <- paste( formulat, ' + file_iEiC[iEs,"scanner"]', sep = "" );
    }
    if( length( unique( file_iEiC[iEs,"cell"] ) ) > 1 ){
      formulat <- paste( formulat, ' + file_iEiC[iEs,"cell"]', sep = "" );
    }
    vehicle_iE2 <- file_iEiC[iEs,"vehicle"];
    cell_iE2 <- file_iEiC[iEs,"cell"];
    numSamples_drug_iC <- sapply( unique( cell_iE2 ), function(drug) { sum( !vehicle_iE2[ drug == cell_iE2 ] ) } );
    numSamples_veh_iC <- sapply( unique( cell_iE2 ), function(drug) { sum( vehicle_iE2[ drug == cell_iE2 ] ) } );
    if( any( numSamples_drug_iC ) == 0 || any( numSamples_veh_iC ) == 0 ){
      stop("One of the cell lines is not balanced. i.e. all its samples are vehicles or drugs.")
    }
    numSamples_iD[iD] <- length( iEs );
    numscanners_iD[iD] <- length( unique( file_iEiC[iEs,"scanner"] ) );
    cellLines_iDiC[ iD, unique( file_iEiC[iEs,"cell"] ) ] <- TRUE;
    
    # Run GLM for all experiments available with this drug
    print( paste( "Testing drug ", drugName_iD[iD], sep = "" ) );
    a_iPEiG <- apply( exp_iEiG[iEs, ], 
                      c(2), 
                      function( exp_iE2 ){
                        model <- glm( formula = as.formula( formulat ),
                                      family = gaussian );  
                        p_iG <- coef( summary( model ) )[ 'file_iEiC[iEs, "vehicle"]TRUE', "Pr(>|t|)" ];
                        e_iG <- coef( summary( model ) )[ 'file_iEiC[iEs, "vehicle"]TRUE', "Estimate" ];
                        return(  c( p_iG, 
                                    e_iG ) );
                      }  );
    pVals_iDiG[iD, ] <- a_iPEiG[1, ];
    effec_iDiG[iD, ] <- a_iPEiG[2, ];
    
  }
  
  # Return results
  drugFeatures_iDiC <- cbind( numSamples = numSamples_iD,
                         numScanners = numscanners_iD,
                         cellLines_iDiC );
  dimnames( drugFeatures_iDiC )[[2]][1] <- "numSamples";
  dimnames( drugFeatures_iDiC )[[2]][2] <- "numScanners";
  r <- list( drugFeatures_iDiC = drugFeatures_iDiC,
             pVals_iDiG = pVals_iDiG,
             effec_iDiG = effec_iDiG );
  return( r );
  
}

findDrugs = function( results_iDiC,
                      targetEffects_iG ){
  
  
}



# 
# 
# 
# a_iPEiG <- apply( exp_iEiG[iEs,1:100], 
#                   c(2), 
#                   function( exp_iE2 ){
#                     model <- glm( formula = as.formula( formulat ),
#                                   family = gaussian );  
#                     p_iG <- coef( summary( model ) )[ 'file_iEiC[iEs, "vehicle"]TRUE', "Pr(>|t|)" ];
#                     e_iG <- coef( summary( model ) )[ 'file_iEiC[iEs, "vehicle"]TRUE', "Estimate" ];
#                     return(  p_iG, 
#                              e_iG );
#                   } );
# 
# 
# 
# 
# 
# file_iEiC <- findFiles( c("paracetamol") )
# exp_iEiG <- loadExpression( file_iEiC[ ,1] );
# results_iDiC <- testDrugs( exp_iEiG,
#                            file_iEiC );
