#!/usr/bin/env Rscript

# Load optparse package to read in command-line arguments.
library("optparse")

version <- "1.0"

option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Location of input file (required).", metavar="PATH"),
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Location of output file (required).", metavar="PATH"),
  
  make_option(c("-f", "--feature_col"), type="integer", default=1,
              help=paste("Last column containing feature name information. The script assumes",
                         "that the feature name and description will be in columns ranging from 1:feature_col.",
                         sep=" "), metavar="INT"),
  
  make_option(c("-p", "--prop"), type="numeric", default=0.2,
              help="Proportion of samples that must have the feature for it to be retained (default: 0.2).",
              metavar="NUMERIC"),
  
  make_option(c("-e", "--exclude"), type="character", default=NULL,
              help="Any features matching the strings in this comma-delimited list will be excluded.",
              metavar="STRING"),
  
  make_option(c("-r", "--rel"), action = "store_true", type="logical", default=FALSE,
              help="Flag to indicate that output table should be re-normalized in relative abundance (default: FALSE).", 
              metavar="BOOL")
  )

  opt_parser <- OptionParser(
    option_list=option_list, 
    usage = "%prog [options] -i FILE",
    description = paste(
      "\nScript to filter a tab-delimited file with features as rows and columns as samples.\n",
      "The first column is assumed to be the feature names unless the feature_col argument is changed.\n",
      "Note that by default the output will not be re-normalized into relative abundance",
      sep=" ")
  )
  
  opt <- parse_args(opt_parser)
  
# Check that input and output files set:
if(is.null(opt$input)) { stop("path to input file needs to be set.") }
if(is.null(opt$output)) { stop("path to output file needs to be set.") }
  
# Read in input file. Note that rownames are not assigned.
intable <- read.table(opt$input, header=T, sep="\t", stringsAsFactors = FALSE, 
                      quote="", comment.char="")

# Get feature names / descriptions only in new df.
feature_info <- intable[ , 1:opt$feature_col, drop=FALSE]

# Get feature abundances only in separate df.
intable_counts <- intable[, -c(1:opt$feature_col), drop=FALSE]

# Identify all features that are observed in fewer than the specified proportion of samples.
features2remove <- which(rowSums(intable_counts > 0) < opt$prop*ncol(intable_counts))

# If features to exclude by name were also set then also identify those rows.
if(! is.null(opt$exclude)) {
  str2exclude <- strsplit(x = opt$exclude, split = ",")[[1]]

  # Loop over all feature columns and strings to exclude and find all matches.  
  for (feature_col in colnames(feature_info)){
      for(exclude in str2exclude) {
        features2remove <- c(features2remove, grep(exclude, feature_info[,feature_col]))
      }
  }

  # Get only unique rows to discard.
  features2remove <- unique(features2remove)
  
}

if(length(features2remove) > 0) {
  # Remove features from table.
  intable_filt <- intable[-features2remove,,drop=FALSE]
} else {
  intable_filt <- intable
}

# Re-normalize data into relative abundance per-sample if option set.
if(opt$rel){
  value_col <- c((opt$feature_col+1):ncol(intable_filt))
  
  intable_filt[,value_col] <- sweep(intable_filt[,value_col], 
                                    2, 
                                    colSums(intable_filt[,value_col]), 
                                    '/') * 100
}

# Write out filtered table to file.
write.table(x = intable_filt, 
            file = opt$output, 
            quote = FALSE, 
            sep = "\t",
            col.names = TRUE, 
            row.names = FALSE)
