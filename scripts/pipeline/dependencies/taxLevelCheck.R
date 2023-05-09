###################################
## CHECK IF --taxlevels IS VALID ##
###################################

# Extract the individual levels from --fuseLevels
individual_levels <- gsub(' ', '', str_to_title(strsplit(fuseLevels, ",")[[1]]))

# Get the configuration file
config_file <- normalizePath(file.path(dirname(pipeline_path), '../../data/reference/config.txt'))

# Read in the csv file
config_df <- read.csv(config_file)

for(i in 1:nrow(config_df)){
  
  # Extract the taxonomic levels per reference database from config_df and convert them to 1 string
  ref_levels <- gsub(' ', '', str_to_title(strsplit(config_df[i, "TaxonomicLevel"], ";")[[1]]))
  
  for(individual_level in individual_levels){
    if(!individual_level %in% ref_levels){
      stop(paste('Invalid taxonomic levels for --f argument: The level of', 
                 individual_level, 
                 'cannot be used for merging taxonomic tables as it is not present in all reference databases.'))
    }
  }
}

##################################################################
## CHECK IF REFERENCE DATABASES ARE COMPATIBLE WITH CONFIG FILE ##
##################################################################

# Get a list of all reference database files in the reference database directory
ref_files <- list.files(path = normalizePath(file.path(dirname(pipeline_path), '../../data/reference/')), full.names = T)

# Remove config file from ref_files
ref_files <- ref_files[!basename(ref_files) == 'config.txt']

# Read in the config csv file
config_df <- read.csv(config_file)

# Get the database names
config_refs <- config_df[,'Database']

# Loop over every database in the reference directory
for(ref_file in ref_files){
  
  # Check if the names of the reference databases match with those in the config file
  if(!basename(ref_file) %in% config_refs){
    stop(paste('The reference database config file is missing information or contains a typo related to the file', basename(ref_file), '. Please review the file and try again.'))
  } 
  
  # Check if the names of the reference databases occur only onece in the config file
  else if(sum(grepl(pattern = basename(ref_file), x = config_refs)) > 1){
    stop(paste('The reference database config file has multiple entries for the', basename(ref_file), 'database. Please review the file and try again.'))
  }
}