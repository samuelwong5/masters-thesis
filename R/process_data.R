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


# =====================E-GEOD-84422 dataset=====================
if (!file.exists("C:/Users/cwong/masters-thesis/data/E-GEOD-84422/E-GEOD-84422-combined-data.csv")) {
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

# =====================E-GEOD-48350 dataset=====================
if (!file.exists("C:/Users/cwong/masters-thesis/data/E-GEOD-48350/E-GEOD-48350-combined-data.csv")) {
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

