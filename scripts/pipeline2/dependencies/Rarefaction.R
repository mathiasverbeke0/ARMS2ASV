suppressWarnings({
  # Remove rows where the sum of the observed count data across entire row is smaller than 2 
  seqtab.nochim <- seqtab.nochim[!rowSums(seqtab.nochim) <= 1,]
   
  S <- specnumber(seqtab.nochim) # observed number of ASV sequences, same as apply(seqtab.nochim, 1, function(x){sum(x != 0)})
  
  (raremax <- min(rowSums(seqtab.nochim))) # rarefaction threshold
  
  # The rarefy() function is used to estimate the expected number of species 
  # (or ASV sequences) in random subsamples of a specified size from a community 
  # or dataset. It is a common method for comparing species richness or diversity 
  # across different samples or communities. In the line below, rarefy(seqtab.nochim, raremax) 
  # calculates the rarefied species richness by simulating the expected number of 
  # species in random subsamples of size raremax from the seqtab.nochim dataset. 
  # The result is assigned to the variable Srare. The rarefaction process involves 
  # randomly selecting a subset of individuals (or ASV sequences) from the dataset 
  # without replacement, and then calculating the number of unique species (or ASV 
  # sequences) present in each subsample. This process is repeated multiple times 
  # to estimate the expected species richness for each subsample size.
  
  Srare <- rarefy(seqtab.nochim, raremax) # rarefy the data
  
  pdf(file = file.path(path.asv_plot, 'ASVRarefaction.pdf')) ## Write to pdf
    
    plot(S, Srare, xlab = "Observed No. of ASVs", ylab = "Rarefied No. of ASVs")
    abline(0, 1) # reference line
    
    rarecurve(seqtab.nochim, step = 20, sample = raremax, col = "blue", cex = 0.6, ylab = 'ASVs') # rarefaction curve
  
  dev.off()
})