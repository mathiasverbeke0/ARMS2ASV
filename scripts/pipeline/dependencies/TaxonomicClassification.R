#library(dada2)
#library(Biostrings)
#library(xlsx)

# Get a list of all reference database files in the reference database directory
ref_files <- list.files(path = normalizePath(file.path(dirname(pipeline_path), '../../data/reference/')), full.names = T)

# Get the configuration file
config_file <- normalizePath(file.path(dirname(pipeline_path), '../../data/reference/config.txt'))

# Remove config file from ref_files
ref_files <- ref_files[!basename(ref_files) == 'config.txt']

# Read in the csv file
config_df <- read.csv(config_file)

# Check if the names of the reference databases match with those in the config file
config_ref <- config_df[,'Database']
print(config_ref)

# Read in the input multifasta file
seqs <- readDNAStringSet(path.ASV)
IDs <- names(seqs)

# Loop through each reference database file
for (ref_file in ref_files) {
  
  # Get the reference database base and name
  ref_base <- basename(ref_file)
  ref_name <- (basename(fs::path_ext_remove(ref_file)))
  
  # Extract the taxonomic levels
  tax_levels <- gsub(' ', '', str_to_title(strsplit(config_df[config_df$Database == ref_base, "TaxonomicLevel"], ";")[[1]]))
  
  cat(paste('Reference database:', ref_name, '\n'))
  
  if(!dir.exists(normalizePath(file.path(path.taxon, ref_name)))){
    cat(paste('Creating output directories:', normalizePath(file.path(path.taxon, ref_name), '\n')))
    dir.create(normalizePath(file.path(path.taxon, ref_name)))
  }
  
  # Create the taxonomy table using the DADA2 assignTaxonomy function
  tax_table <- assignTaxonomy(seqs, ref_file, taxLevels = tax_levels, minBoot = minBoot)
  tax_table <- cbind(ID = IDs, tax_table)
  
  # Print the taxonomy table to an Excel file
  write.xlsx(tax_table, file = normalizePath(file.path(path.taxon, ref_name, paste0(ref_name,'_tax_classification', '.xlsx'))))
}
