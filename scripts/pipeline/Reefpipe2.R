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
  'stringr'
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

cat('\n __  __ ____    _    
|  \\/  / ___|  / \\   
| |\\/| \\___ \\ / _ \\  
| |  | |___) / ___ \\ 
|_|  |_|____/_/   \\_\\\n\n')

result <- system2(command = 'python', args = c(file.path(dirname(pipeline_path), 'dependencies/MSA.py'), '-b', mainpath, '-o', path.unified))

# If the exit status of MSA.py is not 0, throw halt the script
if(result != 0){
  stop('MSA.py script execution failed.')
}

##########################
## CLUSTERING SEQUENCES ##
##########################

cat('\n  ____ _    _   _ 
 / ___| |  | | | |
| |   | |  | | | |
| |___| |__| |_| |
 \\____|_____\\___/ \n\n')

# Message
cat('Clustering the sequences\n')

# Specify taxonomic levels
levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Import aligned sequence data
dna_seq <- readDNAStringSet(file.path(path.unified, 'AlignedCOIASVS.fasta'))


###############################
## CLUSTER THE ASV SEQUENCES ##
###############################

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

# Constructing data frame with ID, sequence and cluster information
information <- data.frame(Sequence = gsub('-', '', as.character(dna_seq)), Cluster = clusters)
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

# Writing this to an Excel file
write.xlsx(x = information, file = file.path(path.unified, 'Union.xlsx'))

warnings()