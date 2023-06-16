#!/usr/bin/env Rscript

#############
## MESSAGE ##
#############

cat(' ____  _____ _____ _____ ____ ___ ____  _____
|  _ \\| ____| ____|  ___|  _ |_ _|  _ \\| ____|
| |_) |  _| |  _| | |_  | |_) | || |_) |  _|   
|  _ <| |___| |___|  _| |  __/| ||  __/| |___  
|_| \\_|_____|_____|_|   |_|  |___|_|   |_____|\n\n')


###############
## FUNCTIONS ##
###############

# Print header
print_header <- function(number){
  headers <- c("\n _ _ _     \n| (_) |__  \n| | | \'_ \\ \n| | | |_) |\n|_|_|_.__/ \n\nInstalling and loading all packages\n",
               "\n ___ _   _ _ __  _ __  \n/ __| | | | '_ \\| '_ \\ \n\\__ \\ |_| | |_) | |_) |\n|___/\\__,_| .__/| .__/ \n          | |   | |    \n          |_|   |_|    \n\n")
  
  cat(headers[number])}

# Download and load R packages
install_and_load_packages <- function(packages) {
  installed_packages <- installed.packages()[, 'Package']                   # Summarize all installed packages
  
  for (package in packages) {
    if (!(package %in% installed_packages)) {                               # If package not installed
      if (package %in% c('Biostrings')) {     
        options("BioC_mirror" = "https://bioconductor.org")                 # Set the Bioconductor image
        BiocManager::install(package)                                       # Install package with BiocManager
      } else {                                                      
        options(repos = "https://cran.rstudio.com/")                        # Set the cran mirror
        install.packages(package)                                           # Install package
      }
    }
    
    suppressPackageStartupMessages(library(package, character.only = TRUE)) # Silently load package
  }
}


####################################
## PARSING COMMAND LINE ARGUMENTS ##
####################################

install_and_load_packages(c('argparse', 'stringr'))

parser <- ArgumentParser(description = 'Reefpipe2.R command line arguments')

# Mandatory command line arguments
parser$add_argument('-b', '--base_dir', metavar = 'BASEDIR', type = 'character', required = TRUE, help = 'The base directory path for the analysis.')
parser$add_argument('-t', '--taxtable', type = 'character', required = TRUE, help = 'The general name of the taxonomic table you want to use.')
parser$add_argument('-e', '--environment', type = 'character', choices = c('linux', 'windows', 'mac'), required = TRUE, help = 'The operating system you are currently using. The choices are linux, windows or mac.')

# Optional command line arguments
parser$add_argument('-T', '--taxlevels', type = 'character', default = 'Phylum,Class,Order,Family,Genus,Species', help = 'The taxonomic levels present in the taxonomic tables. Default levels are Phylum,Class,Order,Family,Genus,Species.')
parser$add_argument('-s', '--similarity', type = 'numeric', required = FALSE, default = 97, help = 'The percentage similarity to cluster sequences.Default is 97.')

# Parse the arguments
args <- parser$parse_args()

# Access the argument values
mainpath <- normalizePath(args$base_dir)
similarity <- args$similarity
taxtable <- args$taxtable
taxlevels <- gsub(' ', '', str_to_title(strsplit(args$taxlevels, ",")[[1]]))
envir <- args$environment


################
## R PACKAGES ##
################

print_header(1)

# Install and load
install_and_load_packages(c('BiocManager','Biostrings','bioseq','openxlsx','dplyr','progress','dplyr','purrr'))


#########################
## MAIN PIPELINE PATH ##
########################

args <- commandArgs()
pipeline_path <- NULL
for (arg in args) {
  if (startsWith(arg, "--file=")) {
    pipeline_path <- sub("^--file=", "", arg)
    break
  }
}

if (!is.null(pipeline_path)) {
  pipeline_path <- normalizePath(pipeline_path)
} else {
  stop("Pipeline path not found.")
}


#############################
## FIND THE RESULTS FOLDER ##
#############################

paths <- normalizePath(list.dirs(mainpath, full.names = TRUE, recursive = TRUE), winslash = '/')

paths.result <- paths[basename(paths) == 'Results']

if(length(paths.result) == 0){
  stop('No Results folder could be found')
}

iter <- 1

print_header(2)

for(path.result in paths.result){
  
  # Iteration message
  cat(paste0('\nIteration ', iter, ' out of ', length(paths.result), ': ', basename(dirname(path.result)), '\n'))
  
  # Iteration label
  label <- paste0('\n[', iter, '/', length(paths.result), ' - Step ')
  
  # Add 1 to iter
  iter <- iter + 1
  
  # Step counter
  step <- 0
  
  
  #############################
  ## CONSTRUCT OUTPUT FOLDER ##
  #############################
  
  path.output <- file.path(path.result, '04.Additional')
  
  if(!dir.exists(path.output)){
    cat(paste('Creating output directory:', path.output, '\n'))
    dir.create(path.output)
  } else{
    
    # Check if there are files in the output folder
    output.contents <- list.files(path.output)
    
    if(length(output.contents) != 0){
      cat(paste(path.output, 'already contains information.\nProceed to next iteration.\n'))
      next
    }
  }
  
  
  ##############################
  ## FIND THE MULTIFASTA FILE ##
  ##############################
  
  # Check if the multifasta file exists
  path.multifasta <- normalizePath(list.files(file.path(path.result, '01.ASV'), pattern = '.fasta$', full.names = T))
  
  if(length(path.multifasta) == 0){
    cat('No ASV multifasta file could be found.\nProceed to next iteration.\n')
    next
  } else if (length(path.multifasta) > 1){
    cat('Multiple ASV multifasta files could be found.\nProceed to next iteration.\n')
    next
  }
  
  
  ##############################
  ## FIND THE TAXONOMIC TABLE ##
  ##############################
  
  # Check if the taxonomic table exist
  taxtable_path <- list.files(path.result, pattern = taxtable, full.names = TRUE, recursive = TRUE)
  
  if(length(taxtable_path) == 0){
    cat(paste('No file called', taxtable, 'could be found.\nProceed to next iteration.\n'))
    next
  } else if(length(taxtable_path) > 1){
    cat(paste('More than 1 file called', taxtable, 'found.\nProceed to next iteration.\n'))
    next
  }
  
  
  ################################################################################
  ## MERGING AND SUBSEQUENT MULTIPLE SEQUENCE ALIGNMENT OF ALL MULTIFASTA FILES ##
  ################################################################################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Performing MSA\n'))
  
  # Define the clustalo executable path
  if(envir == 'linux'){exec <- 'clustalo'
  } else if(envir == 'windows'){exec <- file.path(dirname(pipeline_path), 'dependencies/clustal-omega-1.2.2-win64/clustalo.exe')
  } else if(envir == 'mac'){exec <- file.path(dirname(pipeline_path), 'dependencies/clustal-omega-1.2.3-macosx')}
  
  result <- system2(command = 'python3', args = c(paste0("\"", file.path(dirname(pipeline_path), 'dependencies/MSA.py'), "\""), 
                                                 '-m', paste0("\"", path.multifasta, "\""), 
                                                 '-o', paste0("\"", path.output, "\""), 
                                                 '-e', envir, 
                                                 '-c', paste0("\"", exec, "\"")))
  
  # If the exit status of MSA.py is not 0, halt the script
  if(result != 0){
    cat('MSA.py script execution failed.\nProceed to next iteration.\n')
    next
  }
  
  
  ##########################
  ## CLUSTERING SEQUENCES ##
  ##########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Clustering sequences\n'))
  
  # Import aligned sequence data
  dna_seq <- readDNAStringSet(file.path(path.output, 'Alignment.fasta'))
  
  # Transform aligned sequence data
  dna_seq <- dna(as.vector(dna_seq))
  
  # Calculate the correct threshold value for clustering
  threshold <- 1 - (similarity/100)
  
  # Cluster sequences based on similarity
  clusters <- seq_cluster(x = dna_seq, threshold = threshold, method = "single")
  cat(paste('Clustered', length(clusters), 'ASVs into', max(clusters), 'clusters\n'))
  
  
  #############################
  ## MAIN GENERAL EXCEL FILE ##
  #############################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Supplementing and grouping taxonomic tables\n'))
  
  cat(paste('Generating', file.path(path.output, 'Comprehensive.xlsx\n')))
  
  # Constructing data frame with ID, sequence and cluster comprehensive
  comprehensive <- data.frame(Sequence = gsub('-', '', as.character(dna_seq)), Cluster = clusters)
  rownames(comprehensive) <- names(dna_seq)
  
  # Remove the > character from the ASV ID
  rownames(comprehensive) <- gsub('^(&gt;|>)', '', rownames(comprehensive))
  
  # Reading in all of the data
  taxonomic_table <- read.xlsx(xlsxFile = taxtable_path, sheet = 1)
  
  # Check if the taxonomic levels are present in the taxonomic table
  not_present <- FALSE
  
  for(taxlevel in taxlevels){
    if(!taxlevel %in% colnames(taxonomic_table)){
      
      not_present <- TRUE
    } 
  }
  
  if(not_present){
    cat(paste(taxlevel), "can not be found in", taxtable, '\nProceed to next iteration.\n')
    next
  }
  
  # Remove the > character from the ASV ID
  taxonomic_table$ID <- gsub('^(>|&gt;)', '', taxonomic_table$ID)
  
  # Get the column index of the ID column in taxonomic_table
  taxonomic_index <- which(colnames(taxonomic_table) %in% c('ID', 'NA.'))
  
  # Get all of the column names from taxonomic_table
  all_names <- colnames(taxonomic_table)
  all_names <- all_names[-which(all_names %in% c('ID', 'NA.'))]
  
  # Add these column names to comprehensive data frame
  comprehensive[, all_names] <- NA
  
  for(name in rownames(comprehensive)){
    row <- taxonomic_table[taxonomic_table$ID == name,]
    row <- row[-taxonomic_index]
    row <- as.character(row)
    
    comprehensive[name, all_names] <- row
  }
  
  # Get paths to seqtab.rds files
  seqtab_file <- list.files(path.result, pattern = "seqtab.rds$", recursive = TRUE, full.names = TRUE)
  
  if(length(seqtab_file) == 0){
    cat('No seqtab.rds file could be found.\nProceed to next iteration.\n')
    next
  } else if(length(seqtab_file) > 1){
    cat('Multiple seqtab.rds files could be found.\nProceed to next iteration.\n')
    next
  }
  
  # Read all of the seqtab.rds files
  seqtab <- readRDS(seqtab_file)
  
  # Loop over all the comprehensive rows
  for(sequence in comprehensive$Sequence){
    samples <- c()
    
    # Loop over all the samples of the sequence table
    for(sample in rownames(seqtab)){
      
      if(seqtab[sample, sequence] > 0){
        
        samples <- c(samples, sample)
        
        comprehensive[comprehensive$Sequence == sequence, sample] <- seqtab[sample, sequence]
      } else{
        comprehensive[comprehensive$Sequence == sequence, sample] <- 0
      }
    }
    
    # Sort the samples
    samples <- sort(samples)
    
    # Concatenate the samples
    sample_string <- paste(samples, collapse = ';')
    
    # Add the concatenated samples to the comprehensive data frame
    comprehensive[comprehensive$Sequence == sequence, 'Samples'] <- sample_string 
  }
  
  # Get the sample names
  sample_names <- rownames(seqtab)
  
  # Rearrange the sample columns
  sample_names <- sort(sample_names)
  
  # Put the sample columns in a separate data frame
  sample_columns <- comprehensive[, sample_names, drop = FALSE]
  
  # Remove the sample columns from the comprehensive data frame
  comprehensive <- comprehensive[,!colnames(comprehensive) %in% sample_names]
  
  # Add the sample columns (sorted) to the comprehensive data frame
  comprehensive <- data.frame(comprehensive, sample_columns)
  
  # Writing this to an Excel file
  write.xlsx(x = comprehensive, file = file.path(path.output, 'Comprehensive.xlsx'), rowNames = TRUE)
  
  
  ######################################
  ## GROUP ALL ASVs BASED ON TAXONOMY ##
  ######################################
  
  cat(paste('Generating', file.path(path.output, 'GroupedTaxa.xlsx\n')))
  
  # Make an comprehensive backup
  backup_comprehensive <- comprehensive
  
  # Replace NA values with a placeholder value (e.g., "NA")
  comprehensive[is.na(comprehensive)]<- "#N/A"
  
  # Specify the formula for grouping
  formula <- as.formula(paste(". ~", paste(taxlevels, collapse = " + ")))
  
  # Join rows based on identical taxlevels
  joined <- aggregate(formula, data = comprehensive, FUN = paste, collapse = ', ')
  
  # Replace 'NA' values back to NA values
  joined[joined == '#N/A'] <- NA
  
  # Sort the rows by taxonomic level
  joined <- arrange(joined, !!!syms(taxlevels))
  
  # Write the taxlevel grouped data frame to an Excel file
  write.xlsx(x = joined, file = file.path(path.output, 'GroupedTaxa.xlsx'), rowNames = F)
  
  
  ############################
  ## GET SAMPLE INFORMATION ##
  ############################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Getting sample information from ENA\n'))
  
  system2(command = 'python3', args = c(paste0("\"", file.path(dirname(pipeline_path), 'dependencies/location.py'), "\""), 
                                       '-o', paste0("\"", file.path(path.output, 'samples.xlsx'), "\""), 
                                       '-e', paste(sample_names, collapse = ',')))
  
  cat('\n')
}