library(vegan)

# Read in the data (this is not yet species data)
data <- readRDS(file = '/home/guest/Traineeship/ARMSProject/data/temp2/ARMS_SWC_Gbg2_20200206_20200520/06.Seq_Table/seqtab.rds')

# Remove rows where the sum of the observed count data across entire row is smaller than 2 
data <- data[!rowSums(data) <= 1,]

S <- specnumber(data) # observed number of ASV sequences, same as apply(data, 1, function(x){sum(x != 0)})

(raremax <- min(rowSums(data))) # rarefaction threshold


# The rarefy() function is used to estimate the expected number of species 
# (or ASV sequences) in random subsamples of a specified size from a community 
# or dataset. It is a common method for comparing species richness or diversity 
# across different samples or communities. In the line below, rarefy(data, raremax) 
# calculates the rarefied species richness by simulating the expected number of 
# species in random subsamples of size raremax from the data dataset. 
# The result is assigned to the variable Srare. The rarefaction process involves 
# randomly selecting a subset of individuals (or ASV sequences) from the dataset 
# without replacement, and then calculating the number of unique species (or ASV 
# sequences) present in each subsample. This process is repeated multiple times 
# to estimate the expected species richness for each subsample size.

Srare <- rarefy(data, raremax) # rarefy the data

plot(S, Srare, xlab = "Observed No. of ASVs", ylab = "Rarefied No. of ASVs")
abline(0, 1) # reference line

rarecurve(data, step = 20, sample = raremax, col = "blue", cex = 0.6, ylab = 'ASVs') # rarefaction curve
