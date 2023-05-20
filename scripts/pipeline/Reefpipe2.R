############################
## DOWNLOADING R PACKAGES ##
############################

cat('\n _ _ _     
| (_) |__  
| | | \'_ \\ 
| | | |_) |
|_|_|_.__/ \n')

cat('\nInstalling and loading all packages\n')

pkg <- installed.packages()[,'Package']

ToInstall <- c(
  'argparse',
  'Biostrings',
  'bioseq',
  'xlsx',
  'dplyr',
  'progress',
  'stringr',
  'dplyr', 
  'purrr'
)

for (item in ToInstall){
  if (!item %in% pkg) {
    install.packages(item)
  }
  
  suppressPackageStartupMessages(library(item, character.only = T))
}


####################################
## PARSING COMMAND LINE ARGUMENTS ##
####################################

parser <- ArgumentParser(description = 'Reefpipe2 command line arguments')

# Mandatory command line arguments
parser$add_argument('-b', '--base_dir', metavar = 'BASEDIR', type = 'character', required = TRUE, help = 'The base directory path for the analysis.')
parser$add_argument('-t', '--taxtable', type = 'character', required = TRUE, help = 'The general name of the taxnomic tables you want to use.')
parser$add_argument('-T', '--taxlevels', type = 'character', default = 'Phylum,Class,Order,Family,Genus,Species', help = 'The taxonomic levels extracted from the taxonomic tables. Default levels are Phylum,Class,Order,Family,Genus,Species.')

# Optional command line arguments
parser$add_argument('-s', '--similarity', type = 'numeric', required = FALSE, default = 97, help = 'The percentage similarity to cluster sequences.')

# Parse the arguments
args <- parser$parse_args()

# Access the argument values
mainpath <- normalizePath(args$base_dir)
similarity <- args$similarity
taxtable <- args$taxtable
taxlevels <- gsub(' ', '', str_to_title(strsplit(args$taxlevels, ",")[[1]]))


####################################
## ABSOLUTE PATH OF MAIN WORKFLOW ##
####################################

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


##################################
## CONSTRUCT MAIN OUTPUT FOLDER ##
##################################

path.unified <- file.path(mainpath, 'Unified')

if(!dir.exists(path.unified)){
  dir.create(path.unified)
}


###############################################################################
## MERGING AND SUBSEQUENT MULTIPLE SEQUENCE ALIGNMENT OF ALL MULTIFASTA FILES##
###############################################################################

cat("\n _ __ ___  ___  __ _ 
| '_ ` _ \\/ __|/ _` |
| | | | | \\__ \\ (_| |
|_| |_| |_|___/\\__,_|\n\n")

result <- system2(command = 'python', args = c(file.path(dirname(pipeline_path), 'dependencies/MSA.py'), '-b', mainpath, '-o', path.unified))

# If the exit status of MSA.py is not 0, throw halt the script
if(result != 0){
  stop('MSA.py script execution failed.')
}


##########################
## CLUSTERING SEQUENCES ##
##########################

cat("\n      _       
     | |      
  ___| |_   _ 
 / __| | | | |
| (__| | |_| |
 \\___|_|\\__,_|\n\n")

# Message
cat('Clustering the sequences\n')

# Specify taxonomic levels
levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Import aligned sequence data
dna_seq <- readDNAStringSet(file.path(path.unified, 'AlignedCOIASVS.fasta'))

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

cat("\n      _          
     | |         
__  _| |_____  __
\\ \\/ / / __\\ \\/ /
 >  <| \\__ \\>  < 
/_/\\_\\_|___/_/\\_\\\n\n")

cat(paste('Generating', file.path(path.unified, 'Union.xlsx\n')))

# make new ASV IDs
ASVID_new <- paste0("ASV", seq(1, length(as.character(dna_seq))))

# Constructing data frame with ID, sequence and cluster information
information <- data.frame(ASVID = ASVID_new, Sequence = gsub('-', '', as.character(dna_seq)), Cluster = clusters)
rownames(information) <- names(dna_seq)

# Getting separate taxonomic classification tables
taxtable_paths <- list.files(mainpath, pattern = taxtable, full.names = TRUE, recursive = TRUE)

# Remove the > character from the ASV ID
rownames(information) <- gsub('>', '', rownames(information))

# Constructing empty data frame
all_taxonomy <- data.frame()

# Adding all taxonomic classification tables together
for(taxtable_path in taxtable_paths){
  temporary_dataframe <- read.xlsx(file = taxtable_path, sheetIndex = 1)
  all_taxonomy <- rbind(all_taxonomy, temporary_dataframe)
}

# Remove the > character from the ASV ID
all_taxonomy$ID <- gsub('>', '', all_taxonomy$ID)

# Get the column index of the ID column in all_taxonomy
all_IDIndex <- which(colnames(all_taxonomy) %in% c('ID', 'NA.'))

# Get all of the column names from all_taxonomy
all_names <- colnames(all_taxonomy)
all_names <- all_names[-which(all_names %in% c('ID', 'NA.'))]

# Add these column names to information data frame
information[, all_names] <- NA

for(name in rownames(information)){
  
  row <- all_taxonomy[all_taxonomy$ID == name,]
  
  if(nrow(row) > 1){
    stop('Multiple rows for the same ASV detected. This is not possible.')
  }
  
  if(nrow(row) == 1){
    row <- row[-all_IDIndex]
    row <- as.character(row)
    
    information[name, all_names] <- row
    
  }
  
  else{
    warning(paste('No matching row for ASV', name))
  }
}

# Get paths to seqtab.rds files
seqtab_files <- list.files(mainpath, pattern = "seqtab.rds$", recursive = TRUE, full.names = TRUE)

# Read all of the seqtab.rds files
seqtab_list <- map(seqtab_files, function(x){readRDS(x)})

cat('Determining samples\n')

sample_names <- c()

# Loop over all the information rows
for(row in 1:nrow(information)){
  sequence <- information[row, 'Sequence']
  samples <- c()
  
  # Loop over all the sequence tables
  for(seqtab in seqtab_list){
    
    # Check if the ASV sequence is present in the sequence table
    if(sequence %in% colnames(seqtab)){
      
      # Loop over all the samples of the sequence table
      for(seqtab_row in 1:nrow(seqtab)){
        
        # Check if the abundance of the sequence in the sample is bigger than 0
        if(seqtab[seqtab_row, sequence] > 0){
          samples <- c(samples, rownames(seqtab)[seqtab_row])
          information[information$Sequence == sequence, rownames(seqtab)[seqtab_row]] <- seqtab[seqtab_row, sequence]
          
          if(!rownames(seqtab)[seqtab_row] %in% sample_names){
            sample_names <- c(sample_names, rownames(seqtab)[seqtab_row])
          }
        }
      }
    }
  }
  
  # Sort the samples
  samples <- sort(samples)
  
  # Concatenate the samples
  sample_string <- paste(samples, collapse = ';')
  
  # Add the concatenated samples to the information data frame
  information[information$Sequence == sequence, 'Samples'] <- sample_string 
}

# Rearrange the sample columns
sample_names <- sort(sample_names)

# Put the sample columns in a separate data frame
sample_columns <- information[, sample_names]

# Remove the sample columns from the information data frame
information <- information[,!colnames(information) %in% sample_names]

# Add the sample columns (sorted) to the information data frame
information <- data.frame(information, sample_columns)

# Writing this to an Excel file
write.xlsx(x = information, file = file.path(path.unified, 'Union.xlsx'))


######################################
## GROUP ALL ASVs BASED ON TAXONOMY ##
######################################

cat(paste('Generating', file.path(path.unified, 'GroupedTaxa.xlsx\n')))

# Make an information backup
backup_information <- information

# Replace NA values with a placeholder value (e.g., "NA")
information[is.na(information)]<- "#N/A"

# Specify the formula for grouping
formula <- as.formula(paste(". ~", paste(taxlevels, collapse = " + ")))

# Join rows based on identical taxlevels
joined <- aggregate(formula, data = information, FUN = paste, collapse = ', ')

# Replace 'NA' values back to NA values
joined[joined == '#N/A'] <- NA

# Sort the rows by taxonomic level
joined <- arrange(joined, !!!syms(taxlevels))

# Add up the numbers in the sample columns
joined[, sample_names] <- apply(joined[, sample_names], 1, function(x){
  as.character(x)
  
})

# Write the taxlevel grouped data frame to an Excel file
write.xlsx(x = joined, file = file.path(path.unified, 'GroupedTaxa.xlsx'), row.names = F)


##########################################
## GROUP ALL ASVs BASED ON ASV SEQUENCE ##
##########################################

cat(paste('Generating', file.path(path.unified, 'GroupedSequences.xlsx\n')))

# Use backup_information to override information
information <- backup_information

# Make a column that contains the concatenated string of all taxonomy columns
information$Taxonomy <- apply(information[, taxlevels], 1, function(x){
  paste(as.character(x), collapse = ' ')
})

# Drop columns specified in taxlevels from the data frame information
information <- information[, !(colnames(information) %in% taxlevels)]

# Replace NA values with a placeholder value (e.g., "NA")
information[is.na(information)]<- "#N/A"

# Join rows based on identical ASV sequence
joined <- aggregate(. ~ Sequence, data = information, FUN = paste, collapse = ', ')

# Replace 'NA' values back to NA values
joined[joined == '#N/A'] <- NA

# Sort the rows by ASV sequence
joined <- arrange(joined, Sequence)

# Write the sequence grouped data frame to an Excel file
write.xlsx(x = joined, file = file.path(path.unified, 'GroupedSequences.xlsx'), row.names = F)