library(vegan)

data <- readRDS(file = '/home/guest/Internship/data/output/RUN2/06.Seq_Table/seqtab.rds')

rowSums(data)

S <- specnumber(data[c(1:2, 4:5),]) # observed number of species
(raremax <- min(rowSums(data[c(1:2, 4:5),])))
Srare <- rarefy(data[c(1:2, 4:5),], raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(data[c(1:2, 4:5),], step = 20, sample = raremax, col = "blue", cex = 0.6)

rowSums(data[])


# Examle code
install.packages('vegan')

library(vegan)

data(BCI)

BCI['Abarema.macradenia']


S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)