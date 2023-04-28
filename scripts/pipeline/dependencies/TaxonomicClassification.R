library(dada2)
library(Biostrings)
library(xlsx)

# Get a list of all reference database files in the reference database directory
ref_files <- list.files(path = '~/Traineeship/ARMSProject/data/databases/', full.names = T)

# Read in the input multifasta file
seqs <- readDNAStringSet("/home/guest/Traineeship/ARMSProject/data/output/Gothenburg/Gothenburg/06.Seq_Table/COI_ASVS.fasta")
IDs <- names(seqs)

# Loop through each reference database file
for (ref_file in ref_files) {
  ref_name <- (basename(fs::path_ext_remove(ref_file)))
  cat(paste('Reference database', ref_name, '\n'))
  
  # Create the taxonomy table using the DADA2 assignTaxonomy function
  tax_table <- assignTaxonomy(seqs, ref_file, minBoot = 90)
  tax_table <- cbind(ID = IDs, tax_table[,1:7])
  
  # Print the taxonomy table to an Excel file
  write.xlsx(tax_table, file = file.path('/home/guest/Traineeship/ARMSProject/data/output/Gothenburg/Gothenburg/07.Taxonomy', paste0(ref_name,'tax_classification', '.xlsx')))
  break
}

for(ref_file in ref_files){
  motiv <- readLines(ref_file)
  motiv <- length(strsplit(x = motiv, split = ';'))
  print(motiv[1])
}
