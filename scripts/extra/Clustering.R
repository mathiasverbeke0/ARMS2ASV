library(Biostrings)
library(bioseq)
library(xlsx)
library(dplyr)
library(progress)

# Specify taxonomic levels
levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Import aligned sequence data
dna_seq <- readDNAStringSet("/home/guest/Traineeship/ARMSProject/scripts/extra/aligned_sequences.fasta")


###############################
## CLUSTER THE ASV SEQUENCES ##
###############################

# Transform aligned sequence data
dna_seq <- dna(as.vector(dna_seq))

# Cluster sequences based on similarity (>= 97%)
clusters <- seq_cluster(x = dna_seq, threshold = 0.03, method = "single")


#############################################
## READ AND ALTER THE CONSENSUS EXCEL FILE ##
#############################################

# Read in the consensus Excel file
consensus <- read.xlsx(file = "/home/guest/Traineeship/ARMSProject/data/temp2/ARMS_SWC_Gbg3_20200206_20200529/07.Taxonomy/Consensus.xlsx", sheetIndex = 1)

# Add a cluster column
consensus$Cluster <- clusters

# Extract all ASV IDs
ASV_IDs <- gsub(pattern = '>', replacement = '', x = consensus[, 'ID'])


###############################################
## COMBINE CONSENSUS AND CLUSTER INFORMATION ##
###############################################

# Initialize a progress bar
pb <- progress_bar$new(total = length(unique(clusters)))

# Loop over all clusters
for(cluster in unique(clusters)){
  
  # Update the progress bar
  pb$tick()
  
  # Check if there is more than one sequence in a cluster
  if(length(clusters[clusters == cluster]) > 1){
    
    
    #################################################
    ## EXTRACT INFORMATION FROM MULTI ASV CLUSTERS ##
    #################################################
    
    # Extract the ASV_IDs from all sequences in the cluster
    cluster_members <- names(clusters[clusters == cluster])
    
    # Make a cluster data frame
    cluster_df <- data.frame()
    
    # Add the cluster members to the cluster data frame
    for(cluster_member in cluster_members){
      cluster_df <- rbind(cluster_df, consensus[which(cluster_member == ASV_IDs),])
    }
    
    # Sort the cluster members based on the taxonomic levels
    clusterSorted_df <- arrange(cluster_df, !!!syms(levels))
    
    # Construct taxonomic strings for each row without NA values
    cluster_strings <- apply(clusterSorted_df[,levels], 1, function(row){
      non_na_values <- row[!is.na(row)]
      cluster_line <- paste(non_na_values, collapse = ' ')
      return(cluster_line)
    })
    
    # Remove all empty strings
    cluster_strings <- cluster_strings[cluster_strings != '']
    
    
    #########################################################
    ## SUPPLEMENT TAXONOMY TABLE USING CLUSTER INFORMATION ##
    #########################################################
    
    # Skip to the next cluster if there are no non-empty strings
    if(length(cluster_strings) == 0){
      next
    } 
    
    # Check if there is only 1 non-empty string
    else if(length(cluster_strings) == 1){
      
      # Assign the taxonomy of the only identified ASV to all taxonomically unidentified ASVs within the same cluster
      clusterSorted_df[,levels] <- clusterSorted_df[1, levels]
      clusterSorted_df[2:nrow(cluster_df),'source'] <- 'Cluster representative sequence'
    } 
    
    # Check if there is only 1 unique nom-empty string
    else if(length(unique(cluster_strings)) == 1){
      
      # If so, check if the length of the strings equals the number of ASVs in the cluster
      if(length(cluster_strings) == nrow(clusterSorted_df)){
        
        # If so, Skip to the next cluster
        next
      }
      
      else{
        
        # Define the all.na function
        all.na <- function(x) {
          all(is.na(x))
        }
        
        # Identify rows with all empty levels
        all_nas <- apply(clusterSorted_df[, levels], 1, all.na)
        
        # Assign taxonomic info to unidentified ASVs in incomplete rows
        clusterSorted_df[all_nas, levels] <- clusterSorted_df[1, levels]
        clusterSorted_df[all_nas, 'source'] <- 'Cluster root sequence'
      }
    }
    
    ################################
    ## UPDATE THE CONSENSUS TABLE ##
    ################################
    
    for(row in 1:nrow(clusterSorted_df)){
      consensus[consensus$ID == clusterSorted_df[row, 'ID'],] <- clusterSorted_df[row,]
    }
  }
}

consensus[consensus$source == 'Cluster root sequences',]

write.xlsx(consensus, file = '/home/guest/Traineeship/ARMSProject/data/temp2/ARMS_SWC_Gbg3_20200206_20200529/07.Taxonomy/Clustered.xlsx')
