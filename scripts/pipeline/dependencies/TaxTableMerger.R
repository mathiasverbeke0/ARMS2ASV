##########################
## Specifying variables ##
##########################
fuseLevels2 <- gsub(' ', '', str_to_title(strsplit(fuseLevels, ",")[[1]]))

# Check if genus and species are in the fuseLevels2
genus_species <- FALSE

if('Genus' %in% fuseLevels2 & 'Species' %in% fuseLevels2){
  genus_species <- TRUE
}

#######################################
## Fetching all taxonomy table paths ##
#######################################

# Get a list of all the taxonomic files (excel) within the 07.Taxonomy subdirectories
tax_files <- list.files(path.taxon, pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)

# Ignore previously generated excel files upon rerunning script in the same parent folder
tax_files <- tax_files[!basename(tax_files) %in% c('Ambiguous.xlsx', 'Consensus.xlsx')]


###########################################
## Reading in BOLDSYSTEMS classification ##
###########################################

# If the BOLDSYSTEMS directory exists, load the taxonomic file
if (boldigger ==T){
  
  bold_dir <- file.path(path.taxon, 'BOLDSYSTEMS')
  
  # Get the BOLDSYSTEMS taxonomic file path
  bold_file <- list.files(bold_dir, pattern = "\\.xlsx$", full.names = TRUE)
  
  # Read in the BOLDSYSTEMS taxonomic file
  bold_taxa <- read_xlsx(bold_file, sheet = 'Sheet1')[,c("ID", fuseLevels2, "Similarity")]
  
  # Add the source
  bold_taxa$source <- 'BOLDSYSTEMS'
  
  # Replace all NA values in similarity column with 0
  bold_taxa[is.na(bold_taxa$Similarity), 8] <- as.character(0)
  
  # Nullify all taxonomies where similarity is less than 97% (= species level)
  bold_taxa[as.numeric(bold_taxa$Similarity) < 97, 2:8] <- NA
  
  # Remove the BOLDSYSTEMS directory path
  bold_index <- grepl(pattern = bold_dir, x = tax_files)
  tax_files <- tax_files[-bold_index]
}


######################################################
## Reading and processing the other taxonomic files ## 
######################################################
if(reference == T){
  
  # Iterate over the taxonomic files and store the resulting data frames in a list
  tax_tables <- map(tax_files, function(x){read_xlsx(x, sheet = 'Sheet1')[,c("ID", fuseLevels2)]})
  
  # Name the taxonomic tables using the base names of their directories
  tax_names <- tools::file_path_sans_ext(basename(tax_files))
  names(tax_tables) <- tax_names
  
  # Add the source to taxonomic tables
  tax_tables <- map(tax_names, function(x){
    tax_tables[[x]]$source <- x
    return(tax_tables[[x]])
  })
  
  
  ################################################
  ## Data manipulation of other taxonomic files ##
  ################################################
  
  # Replace all #N/A values with NA values
  tax_tables <- map(tax_tables, function(x){
    x <- replace(x, x == '#N/A', NA)
    return(x)
  })
  
  # Remove the genus from the species name
  if(genus_species == T){
    tax_tables <- map(tax_tables, function(x){
      x['Species'] <- apply(x, 1, function(row){
        row['Species'] <- trimws(gsub(pattern = row['Genus'], replacement = '', x = row['Species']))
      })
      return(x)
    })
  }
  
  # Remove duplicate rows in the taxonomic tables (these should not be there in the first place)
  tax_tables <- map(tax_tables, function(x){
    x[!duplicated(x), ]
  })
  
  # Define a function to merge two data frames
  merge_tax <- function(df1, df2) {
    # Merge rows that have the same classification in both data frames
    merge_cols <- c('ID', fuseLevels2)
    merged_raw <- full_join(df1, df2, by = merge_cols)
    
    # Add all sources to the merged_raw data frame
    merged_raw$source <- ifelse(
      is.na(merged_raw$source.x), 
      merged_raw$source.y, 
      ifelse(
        is.na(merged_raw$source.y),
        merged_raw$source.x,
        paste(merged_raw$source.x, merged_raw$source.y, sep = ";")
      )
    )
    
    # Remove individual sources
    merged_raw <- merged_raw[, -which(names(merged_raw) %in% c("source.x", "source.y"))]
    
    # Sort the rows by taxonomic level
    merged_sorted <- arrange(merged_raw, ID, !!!syms(fuseLevels2))
    
    # Create an empty data frame to store the merged data
    merged_final <- data.frame()
    
    # Loop through each unique ID in the merged data frame
    for(ID in unique(merged_sorted$ID)){
      
      # Select only the rows with the current ID
      merged_ID <- merged_sorted[merged_sorted$ID == ID,]
      
      # Select only the taxonomic columns
      sub_merged_ID <- merged_ID[,fuseLevels2]
      
      # Create a single taxonomic line for each row (without NA values)
      tax_lines <- apply(sub_merged_ID, 1, function(row) {
        non_na_values <- row[!is.na(row)]
        tax_line <- paste(non_na_values, collapse = ' ')
        return(tax_line)
      })
      
      # Check if there are both non-empty and empty strings in tax_lines
      if(any(tax_lines != "") & any(tax_lines == "")) {
        # If so, filter out empty strings from tax_lines
        tax_lines <- tax_lines[tax_lines != ""]
      }
      
      # Create empty vectors to store the taxonomic lines and their indexes
      tax_dump <- c()
      tax_index <- c()
      
      # Loop through each taxonomic line
      for(i in 1:length(tax_lines)){
        if(length(tax_dump) == 0){
          tax_dump <- tax_lines[i]
          tax_index <- c(i)
          next
        }
        
        # Check if the current line is not a substring of any item already in tax_dump, 
        # or if it is a substring that appears more than once; if not, add it to tax_dump
        if(!sum(grepl(pattern = tax_lines[i], tax_dump)) >= 1){
          tax_dump <- c(tax_dump, tax_lines[i])
          tax_index <- c(tax_index, i)
        }
      }
      
      # Add the selected rows to the final merged data frame
      merged_final <- rbind(merged_final, merged_ID[tax_index,])
    }
    
    # Return the final merged data frame
    return(merged_final)
  }
  
  # Merge the data frames
  merged <- reduce(tax_tables, merge_tax)
  
  
  #########################################
  ## Look for ambiguous classifications ###
  #########################################
  
  # Ignore previously generated log file upon rerunning script in the same folder
  if (file.exists(file.path(path.taxon, "TaxLog.txt"))){
    file.remove(file.path(path.taxon, "TaxLog.txt"))
    cat(paste('Ambiguous taxonomies in reference database tables:\n',
              '--------------------------------------------------\n'), 
        file = file.path(path.taxon, 'TaxLog.txt'))
  } else {
    cat(paste0('Ambiguous taxonomies in reference database tables:\n',
              '--------------------------------------------------\n'), 
        file = file.path(path.taxon, 'TaxLog.txt'))
  }
  
  # Make data frame to store ambiguous taxonomies
  ambiguous <- data.frame()
  
  # Loop over every unique ID
  for(ID in unique(merged$ID)){
    
    # Check if there is more than 1 row with the same ID
    if(nrow(merged[merged$ID == ID,])>1){
      
      # Add the ambiguous columns to the ambiguous data frame
      ambiguous <- rbind(ambiguous, merged[merged$ID == ID,])
      
      # Remove the ambiguous classifications from the merged data frame
      merged <- merged[merged$ID != ID,]
      
      # Write the IDs for ambiguous taxonomies to a log file
      cat(ID, 
          file = file.path(path.taxon, 'TaxLog.txt'), 
          append = T)
    }
  }
  
  if(nrow(ambiguous) != 0){
    # Storing ambiguous taxonomy table in excel file if it is not empty
    write.xlsx(ambiguous, file = file.path(path.taxon, 'Ambiguous.xlsx'))
  } else{
    cat('None\n', 
        file = file.path(path.taxon, 'TaxLog.txt'))
  }
}


######################################################################
## Supplement BOLDSYSTEMS taxonomy table with other taxonomy tables ## 
######################################################################

if(boldigger == T & reference == T){
  
  # Iterate over rows of BOLDSYSTEMS taxonomy table
  for(row in 1:nrow(bold_taxa)){
    
    # ASV ID
    ASV_ID <- sub("^>", "", bold_taxa$ID[row])
    
    ####################################################
    ## If all tax levels for BOLDSYSTEMS ASV are full ##
    ####################################################
    
    # Continue to the next row
    if(all(!is.na(bold_taxa[row, fuseLevels2]))){
      next
    }
    
    #####################################################
    ## If all tax levels for BOLDSYSTEMS ASV are empty ##
    #####################################################
    
    # If an ASV in the bold taxonomy has all taxonomic levels as NA and is unambiguous in the other taxonomy data frame...
    else if(all(is.na(bold_taxa[row, fuseLevels2])) & ASV_ID %in% unique(merged$ID)){
      
      # ... check if the ASV has any levels in the other taxonomy data frame
      if(!all(is.na(merged[merged$ID == ASV_ID, fuseLevels2]))){
        
        # Merge the other taxonomy row into the BOLDSYSTEMS taxonomy row, replacing fuseLevels2 and source columns
        bold_taxa[row, c(fuseLevels2, 'source')] <- merged[merged$ID == ASV_ID, c(fuseLevels2, 'source')]
      }
    }
    
    else if(all(is.na(bold_taxa[row, fuseLevels2]))){
      
      cat(paste(ASV_ID, '\n'), file = file.path(path.taxon, 'TaxLog.txt'), append = T)
    }
    
    ############################################################################################
    ## If tax levels of BOLDSYSTEMS ASV is less specific than merged reference database table ##
    ############################################################################################
    
    # If an ASV in the bold taxonomy has taxonomic levels and is unambiguous in the other taxonomy data frame...
    else if(!all(is.na(bold_taxa[row, fuseLevels2])) & ASV_ID %in% unique(merged$ID)){
      
      # Check if the reference databases have more taxonomic levels
      if(sum(!is.na(bold_taxa[row, fuseLevels2])) < sum(!is.na(merged[merged$ID == ASV_ID, fuseLevels2]))){
        
        # Extract the BOLDSYSTEMS and reference database taxonomies
        bold_taxonomy_na <- bold_taxa[row, fuseLevels2]
        reference_taxonomy_na <- merged[merged$ID == ASV_ID, fuseLevels2]
        
        # Remove all missing values
        bold_taxonomy <- bold_taxonomy_na[!is.na(bold_taxonomy_na)]
        reference_taxonomy <- reference_taxonomy_na[!is.na(reference_taxonomy_na)]
        
        # Some species levels are written as genus and species in the reference database tables
        # This needs to updated to only the species
        if(length(reference_taxonomy) == length(fuseLevels2) & genus_species == T){
          reference_taxonomy[length(fuseLevels2)] <- sub(pattern = reference_taxonomy[length(fuseLevels2) - 1], replacement = '', x = reference_taxonomy[length(fuseLevels2)])
          reference_taxonomy[length(fuseLevels2)] <- sub(pattern = ' ', replacement = '', x = reference_taxonomy[length(fuseLevels2)])
        }
        
        # Construct a single string for the bold taxonomy and merged reference taxonomy
        bold_string <- paste(bold_taxonomy, collapse = ' ')
        reference_string <- paste(reference_taxonomy, collapse = ' ')
        
        # Check if BOLDSYSTEMS taxonomy matches with merged reference database table's extended taxonomy
        if(grepl(pattern = bold_string, x = reference_string)){
          
          # Merge BOLDSYSTEMS and reference table sources
          source <- paste(c(bold_taxa[row, 'source'], merged[merged$ID == ASV_ID, 'source']), collapse = ';')
          
          # Update the sources in the BOLDSYSTEMS taxonomic table
          bold_taxa[row, c('source')] <- source
          
          # Update the taxonomic levels in the BOLDSYSTEMS taxonomic table
          bold_taxa[row, c(fuseLevels2)] <- as.list(reference_taxonomy_na)
          
          # Determine which taxonomic levels where supplemented to the BOLDSYSTEMS taxonomy
          bold_length <- length(bold_taxonomy)
          reference_length <- length(reference_taxonomy)
          extra_levels <- fuseLevels2[(bold_length + 1): reference_length]
          extra_string <- paste(extra_levels, collapse = ';')
          
          # Update the supplemented taxonomic levels in the BOLDSYSTEMS table
          bold_taxa[row, 'Supplemented'] <- extra_string
        }
      }
    }
    
    #######################################################################################################
    ## If tax levels of BOLDSYSTEMS ASV is less specific than ambiguous merged reference databases table ##
    #######################################################################################################
    
    else if(!all(is.na(bold_taxa[row, fuseLevels2])) & ASV_ID %in% unique(ambiguous$ID)){
      
      # Extract the BOLDSYSTEMS and ambiguous reference database taxonomies
      bold_taxonomy_na <- bold_taxa[row, fuseLevels2]
      reference_taxonomy_na <- ambiguous[ambiguous$ID == ASV_ID, fuseLevels2]
      
      # Remove all missing values and create single taxonomic line
      bold_taxonomy <- bold_taxonomy_na[!is.na(bold_taxonomy_na)]
      bold_string <- paste(bold_taxonomy, collapse = ' ')
      
      # Create a single taxonomic line for each row for ambiguous table(without NA values)
      reference_string <- apply(reference_taxonomy_na, 1, function(row){
        non_na_values <- row[!is.na(row)]
        
        # Some species levels are written as genus and species in the reference database tables
        # This needs to updated to only the species
        if(length(row) == length(fuseLevels2) & genus_species == T){
          row[length(fuseLevels2)] <- sub(pattern = row[length(fuseLevels2) - 1], replacement = '', x = row[length(fuseLevels2)])
          row[length(fuseLevels2)] <- sub(pattern = ' ', replacement = '', x = row[length(fuseLevels2)])
        }
        
        # Construct a single string for the bold taxonomy and merged reference taxonomy
        row_string <- paste(row, collapse = ' ')
      
        return(row_string)
      })
      
      if(sum(grepl(pattern = bold_string, x = reference_string)) == 1){
        index <- grep(pattern = bold_string, x = reference_string)
        
        # Merge BOLDSYSTEMS and reference table sources
        source <- paste(c(bold_taxa[row, 'source'], ambiguous[ambiguous$ID == ASV_ID, 'source'][index,]), collapse = ';')
        
        # Update the sources in the BOLDSYSTEMS taxonomic table
        bold_taxa[row, c('source')] <- source
        
        # Update the taxonomic levels in the BOLDSYSTEMS taxonomic table
        bold_taxa[row, c(fuseLevels2)] <- as.list(reference_taxonomy_na[index,])
        
        # Determine which taxonomic levels where supplemented to the BOLDSYSTEMS taxonomy
        bold_length <- length(bold_taxonomy)
        reference_length <- sum(!is.na(reference_taxonomy_na[index,]))
        extra_levels <- fuseLevels2[(bold_length + 1): reference_length]
        extra_string <- paste(extra_levels, collapse = ';')
        
        # Update the supplemented taxonomic levels in the BOLDSYSTEMS table
        bold_taxa[row, 'Supplemented'] <- extra_string
        
        # Add an ambiguity warning
        bold_taxa[row, 'Ambiguity'] <- 'Added'
      } else{
        bold_taxa[row, 'Ambiguity'] <- 'Available'
      }
    }
  }
  
  # Storing supplemented BOLDSYSTEMS taxonomic table in excel file
  write.xlsx(bold_taxa, file = file.path(path.taxon, 'Consensus.xlsx'))
} else if(boldigger == F & reference == T){
  write.xlsx(merged, file = file.path(path.taxon, 'Consensus.xlsx'))
}