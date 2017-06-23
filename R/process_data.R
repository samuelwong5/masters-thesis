# ======================UTILITY FUNCTIONS=========================
pretty_print_progress = function(text, curr, total) {
  cat(paste("\r", text, ": ", curr, "/", total, sep=""));
  if (curr == total) {
    cat("\n");
  }
}


vert_stack = function(a, b) {
  d_a <- colnames(a);
  d_b <- colnames(b);
  m_a <- match(d_a, d_b);
  m_b <- match(d_b, d_a);
  a <- a[,!is.na(m_a)];
  b <- b[,!is.na(m_b)];
  return(rbind(a,b));
}
# =====================END UTILITY FUNCTIONS======================




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
#     folder       - top level folder (e.g. "C:/Users/user/data/E-GEOD-884422")
#     dataset      - name of dataset (e.g. "E-GEOD-884422")
#     funcs        - vector of functions to parse metadata from sdrf file
#     adf_id_col   - column name of intermediate identifiers
#     entrez_table - table to convert from intermediate identifiers to entrez
#     adf          - filename of ADF file
#     adf_sep      - delimiter in adf file
#     sample_sep   - delimiter in sample files
#     sdrf_sep     - delimiter in sdrf file

load_data = function(folder,
                     dataset,
                     adf_id_col,
                     entrez_table,
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
  combined <- NULL;
  samples_count <- length(samples_files);
  
  for (i in 1:samples_count) {
    pretty_print_progress("Processing sample files", i, samples_count);
    data <- read.table(file = paste(samples_dir, samples_files[i], sep="/"),
                       header = TRUE, 
                       sep = sample_sep,
                       na.strings = "");
    probe_ids <- data[, 1];
    
    # Match sample probes IDs
    probe_matches <- match(probe_ids, adf_data[, 1]);
    
    # Strip all non-existing rows
    probes <- adf_data[probe_matches, adf_id_col];
    data <- data[!is.na(probes), ];
    probes <- probes[!is.na(probes)];

    # Unique probe IDs
    probes_u <- unique(probes);
    u <- match(probes_u, probes);
    data <- data[u, 2, drop=FALSE];
    
    dimnames(data)[[1]] <- probes_u;
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
  
  #print("Converting gene identifiers into Entrez IDs");
  #entrez <- entrez_table;
  #entrez_exist <- names( entrez ) %in% dimnames( combined )[[1]];
  #entrez <- entrez[ entrez_exist ]; 
  #entrez <- unlist( entrez );
  #combined <- combined[names(entrez), ];
  #entrez_u <- unique(entrez);
  #u <- match(entrez_u, entrez);
  #combined <- combined[u, ];
  
  #dimnames(combined)[[1]] <- entrez_u;
  
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









library( org.Hs.eg.db );
refseq_entrez <- as.list( org.Hs.egREFSEQ2EG );
genbank_entrez <- as.list( org.Hs.egACCNUM2EG );

data_dir <- paste(dirname(dirname(sys.frame(1)$ofile)), 'data', sep='/')
data_fp = function(relative_path) {
  return(paste(data_dir, relative_path, sep='/'));
}

# =====================E-GEOD-84422 dataset=====================
if (!file.exists(data_fp("E-GEOD-84422/E-GEOD-84422-combined-data.csv"))) {
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
      return(grepl("normal", x[["Term.Accession.Number"]], fixed=TRUE));
    }
  )
  data_84422 = load_data(data_fp("E-GEOD-84422"),
                         "E-GEOD-84422",
                         "Composite.Element.Database.Entry.refseq.",
                         refseq_entrez,
                         funcs_84422,
                         "A-AFFY-33.adf.txt",
                         adf_sep="\t",
                         sample_sep="\t",
                         sdrf_sep="\t")
  
  colnames(data_84422)[length(colnames(data_84422))-2] = "gender";
  colnames(data_84422)[length(colnames(data_84422))-1] = "age";
  colnames(data_84422)[length(colnames(data_84422))] = "ad";
  write.csv(data_84422, file=data_fp("E-GEOD-84422/E-GEOD-84422-combined-data.csv"));
} else {
  data_84422 = read.csv(file=data_fp("E-GEOD-84422/E-GEOD-84422-combined-data.csv"));
  colnames(data_84422) <- substring(colnames(data_84422), 2);
}
# =====================END OF E-GEOD-84422======================


# =====================E-GEOD-48350 dataset=====================
if (!file.exists(data_fp("E-GEOD-48350/E-GEOD-48350-combined-data.csv"))) {
  funcs_48350 = c(
    function(x) {   # Gender
      if (grepl("female", x[["Comment..Sample_source_name."]], fixed=TRUE)) {
        return(0);
      } else if (grepl("male", x[["Comment..Sample_source_name."]], fixed=TRUE)) {
        return(1);
      } else if (grepl("female", x[["Characteristics..age.yrs."]], fixed=TRUE)) {
        return(0);
      } 
      return(1);
    },
    function(x) {   # AGE
      if (!is.na(as.numeric(x[["Characteristics..age.yrs."]]))) {
        return(as.numeric(x[["Characteristics..age.yrs."]]));
      } else {
        return(as.numeric(strsplit(x[["Characteristics..apoe.genotype."]], ' ')[[1]][2]));
      }
    },
    function(x) {   # AD
      if (!is.na(as.numeric(x[["Characteristics..braak.stage."]]))) {
        return(1);
      }
      return(0);
    }
  )
  data_48350 = load_data(data_fp("E-GEOD-48350"),
                         "E-GEOD-48350",
                         "Composite.Element.Database.Entry.refseq.",
                         refseq_entrez,
                         funcs_48350,
                         "A-AFFY-44.adf.txt",
                         adf_sep="\t",
                         sample_sep="\t",
                         sdrf_sep="\t")
  colnames(data_48350)[length(colnames(data_48350))-2] = "gender";
  colnames(data_48350)[length(colnames(data_48350))-1] = "age";
  colnames(data_48350)[length(colnames(data_48350))] = "ad";
  write.csv(data_48350, file=data_fp("E-GEOD-48350/E-GEOD-48350-combined-data.csv"));
} else {
  data_48350 = read.csv(file=data_fp("E-GEOD-48350/E-GEOD-48350-combined-data.csv"));
  colnames(data_48350) <- substring(colnames(data_48350), 2);
}
# =====================END OF E-GEOD-48350======================


# =====================E-GEOD-63063 dataset=====================
if (!file.exists(data_fp("E-GEOD-63063/E-GEOD-63063-combined-data.csv"))) {
  funcs_63063 = c(
    function(x) {   # Gender
      if (as.character(x[["Characteristics..sex."]]) == "female") {
        return(0);
      } 
      return(1);
    },
    function(x) {   # AGE
      return(x[["Characteristics..age."]]);
    },
    function(x) {   # AD
      if (as.character(x[["Characteristics..status."]]) == "CTL") { return(0); }
      return(1);
    }
  )
  data_63063 = load_data(data_fp("E-GEOD-63063"),
                         "E-GEOD-63063",
                         "Reporter.Database.Entry..genbank.",
                         genbank_entrez,
                         funcs_63063,
                         "A-GEOD-10558.adf.txt",
                         adf_sep="\t",
                         sample_sep="\t",
                         sdrf_sep="\t")
  colnames(data_63063)[length(colnames(data_63063))-2] = "gender";
  colnames(data_63063)[length(colnames(data_63063))-1] = "age";
  colnames(data_63063)[length(colnames(data_63063))] = "ad";
  write.csv(data_63063, file=data_fp("E-GEOD-63063/E-GEOD-63063-combined-data.csv"));
} else {
  data_63063 = read.csv(file=data_fp("E-GEOD-63063/E-GEOD-63063-combined-data.csv"));
  colnames(data_63063) <- substring(colnames(data_63063), 2);
}
# =====================END OF E-GEOD-63063======================


# =====================E-GEOD-63063-2 dataset=====================
if (!file.exists(data_fp("E-GEOD-63063-2/E-GEOD-63063-2-combined-data.csv"))) {
  funcs_63063_2 = c(
    function(x) {   # Gender
      if (as.character(x[["Characteristics..sex."]]) == "female") {
        return(0);
      } 
      return(1);
    },
    function(x) {   # AGE
      return(x[["Characteristics..age."]]);
    },
    function(x) {   # AD
      if (as.character(x[["Characteristics..status."]]) == "CTL") { return(0); }
      return(1);
    }
  )
  data_63063_2 = load_data(data_fp("E-GEOD-63063-2"),
                         "E-GEOD-63063",
                         "Reporter.Database.Entry..genbank.",
                         genbank_entrez,
                         funcs_63063_2,
                         "A-GEOD-10558.adf.txt",
                         adf_sep="\t",
                         sample_sep="\t",
                         sdrf_sep="\t")
  colnames(data_63063_2)[length(colnames(data_63063_2))-2] = "gender";
  colnames(data_63063_2)[length(colnames(data_63063_2))-1] = "age";
  colnames(data_63063_2)[length(colnames(data_63063_2))] = "ad";
  write.csv(data_63063_2, file=data_fp("E-GEOD-63063-2/E-GEOD-63063-2-combined-data.csv"));
} else {
  data_63063 = read.csv(file=data_fp("E-GEOD-63063-2/E-GEOD-63063-2-combined-data.csv"));
  colnames(data_63063_2) <- substring(colnames(data_63063_2), 2);
}
# =====================END OF E-GEOD-63063======================
combined <- vert_stack(data_84422, data_48350);
combined <- vert_stack(combined, data_63063);
write.csv(combined, file=data_fp("combined.csv"));
