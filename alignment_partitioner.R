if(!require("optparse")) {install.packages("optparse", repos = "http://cran.us.r-project.org")}
library(optparse)

option_list = list(
  make_option(c("-a", "--alignment"), type = "character", default = NULL, 
              help = "alignment file name", metavar = "character"),
  make_option(c("-p", "--partition_file"), type = "character", default = NULL, 
              help = "partition file name", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL, 
              help = "output file name", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$alignment)) {
  print_help(opt_parser)
  stop("Alignment file name must be supplied.n", call. = FALSE)
}
if(is.null(opt$partition_file)) {
  print_help(opt_parser)
  stop("Partition file name must be supplied.n", call. = FALSE)
}
if(is.null(opt$out) & !is.null(opt$alignment)) {
  opt$out <- paste(strsplit(opt$alignment, ".", fixed = T)[[1]][1], "partitioned.phy", sep = "_")
}

non_vect_every_third <- function(x) {
  rawToChar(charToRaw(x)[c(TRUE, FALSE, FALSE)])
}
every_third <- Vectorize(non_vect_every_third)

regroup_parts <- function(alignment, partition_file, output_name) {
  # Read in the alignment:
  if(!file.exists(paste("~/", alignment, sep = ""))) {
    alignment <- readline(prompt = "File not found. Please tell me the full path to file: ") 
  }
  # Read in the partition file:
  if(!file.exists(paste("~/", partition_file, sep = ""))) {
    partition_file <- readline(prompt = "File not found. Please tell me the full path to file: ") 
  }
  aln <- read.table(paste("~/", alignment, sep = ""), header = T, stringsAsFactors = F)
  parts <- read.table(paste("~/", partition_file, sep = ""), sep = "\t", stringsAsFactors = F)
  
  # Make a list whose number of elements is equal to the total number of partitions, with
  # each element being a data frame whose row(s) correspond to one of the blocks that make
  # up each individual partition. Each row of each data frame is associated with a boolean
  # that indicates whether the range specified in the other two columns is to be split up
  # by codon or not.
  splitlist <- list()
  for(i in 1:nrow(parts)) {
    temp_vect <- strsplit(strsplit(parts[i,1], " = ", fixed = T)[[1]][2], ", ", fixed = T)[[1]]
    block_table <- data.frame(matrix(NA, nrow = length(temp_vect), ncol = 3))
    for(j in 1:length(temp_vect)) {
      if(grepl("\\", temp_vect[j], fixed = T) == TRUE) {
        block_table[j,1] <- strsplit(strsplit(temp_vect[j], "\\", fixed = T)[[1]][1], "-", fixed = T)[[1]][1]
        block_table[j,2] <- strsplit(strsplit(temp_vect[j], "\\", fixed = T)[[1]][1], "-", fixed = T)[[1]][2]
        block_table[j,3] <- TRUE
      } else {
        block_table[j,1] <- strsplit(temp_vect[j], "-", fixed = T)[[1]][1]
        block_table[j,2] <- strsplit(temp_vect[j], "-", fixed = T)[[1]][2]
        block_table[j,3] <- FALSE
      }
    }
    splitlist[[i]] <- block_table
  }
  
  # Create a list whose elements are the already concatenated, contiguous partitions:
  part_list <- list()
  for(i in 1:length(splitlist)) {
    # Initialize an empty data frame of the right size. The column(s) of this data frame
    # will correspond to (one or more) contiguous blocks making up a given partition.
    partition_frame <- data.frame(matrix(NA, nrow = nrow(aln), ncol = nrow(splitlist[[i]])), stringsAsFactors = FALSE)
    for(j in 1:nrow(splitlist[[i]])) {
      if(splitlist[[i]][j,3] == TRUE) {
        partition_frame[,j] <- every_third(substr(aln[,2], splitlist[[i]][j,1], splitlist[[i]][j,2]))
      } else {
        partition_frame[,j] <- substr(aln[,2], splitlist[[i]][j,1], splitlist[[i]][j,2])
      }      
    }
    # Create a new data frame whose first column is populated with the taxon names.
    # We then concatenate all columns of partition_frame row-wise so that we end up with
    # a single vector of strings rather than a multi-column data frame, and assign this
    # vector to the new data frame as its second column using cbind.
    part_list[[i]] <- cbind(aln[,1], apply(partition_frame, 1, paste, collapse = ""))
  }
  
  # Calculate the lengths of individual partitions and store them in a vector
  partition_lengths <- vector()
  for(i in 1:length(part_list)) {
    partition_lengths <- c(partition_lengths, nchar(part_list[[i]][1]))
  }
  
  # In PAML-compatible alignments, the individual partitions/datasets are separated by a
  # blank line followed by another line which reports the number of taxa in the partition
  # and the partition length. To insert these two lines between data frames when printing
  # to file, we open a connection and fill in the file using a loop:
  
  fileprefix <- strsplit(partition_file, "_", fixed = T)[[1]][1]
  
  header <- file(output_name, "w")
  for(i in 1:length(part_list)) {
    writeLines(paste("", "", sep = ""), header) # blank line
    writeLines(paste("       ", (nrow(part_list[[i]])), "    ", nchar(part_list[[i]][1,2]), sep = ""), header)
    write.table(part_list[[i]], header, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
  }
  close(header)
}

regroup_parts(opt$alignment, opt$partition_file, opt$out)
