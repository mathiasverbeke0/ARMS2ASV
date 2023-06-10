######################
## Loading packages ##
######################

library(scatterpie)
library(ggplot2)
library(openxlsx)
library(ape)
#install.packages("remotes")
#remotes::install_github("YuLab-SMU/ggtree")
library(ggtree)
library(RColorBrewer)
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


df1 <- read.xlsx(xlsxFile = '/home/guest/Traineeship/ARMSProject/data/base2/Unified/samples.xlsx', sheet = 1, rowNames = F)
df1$long <- as.numeric(df1$long)
df1$lat <- as.numeric(df1$lat)

df1[c(4,5,6), 'lat'] <- 59
df2 <- read.xlsx(xlsxFile = '/home/guest/Traineeship/ARMSProject/data/base2/Unified/Union.xlsx', sheet = 1, rowNames = F)

if(!'sample' %in% colnames(df1)){
  stop('sample names not in df1.')
} else{
  sample_names <- df1$sample
}

if('Species' %in% colnames(df2) & 'Genus' %in% colnames(df2) & 'Sequence' %in% colnames(df2)){
  df2 <- df2[!is.na(df2$Species), c('Sequence', 'Genus', 'Species', sample_names)]
  df2$Genus_Species <- apply(df2, 1, function(x){paste(x['Genus'], x['Species'])})
} else{stop('Species and/or Genus and/or Sequence not present in df2.')}

for(Genus_Species in unique(df2$Genus_Species)){
  
  HT <- unique(df2[df2$Genus_Species == Genus_Species, 'Sequence'])
  HT_names <- paste0('HT', seq(1:length(HT)))
  
  if(length(HT_names) < 2){next}
  
  HT <- setNames(HT, HT_names)
  
  location <- df1
  rownames(location) <- location$sample
  location <- location[,-grepl(pattern = '^sample', x = colnames(location))]
  location[, HT_names] <- 0
  
  for(haplotype in names(HT)){
    
    HT_abundance <- as.character(df2[df2$Genus_Species == Genus_Species & df2$Sequence == HT[haplotype], sample_names][1,])
    
    HT_abundance <- replace(HT_abundance, is.na(HT_abundance) | HT_abundance == "NA", 0)
    
    location[sample_names, haplotype] <- as.numeric(HT_abundance)
  }

  
  ##############################
  ## GROUPING CLOSE LOCATIONS ##
  ##############################
  
  distance <- function(lat1, lat2, lon1, lon2) {
    # Convert degrees to radians
    lon1 <- lon1 * (pi/180)
    lon2 <- lon2 * (pi/180)
    lat1 <- lat1 * (pi/180)
    lat2 <- lat2 * (pi/180)
    
    # Haversine formula
    dlon <- lon2 - lon1
    dlat <- lat2 - lat1
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    
    c <- 2 * asin(sqrt(a))
    
    # Radius of the Earth in kilometers
    r <- 6371
    
    # Calculate the result
    return(c * r)
  }
  
  # Specify where the HT columns can be found
  HT_columns <- grepl(pattern = '^HT', colnames(location))
  
  # Row to remove
  discard <- c()
  
  for(i in 1:(nrow(location) - 1)){
    if(i %in% discard){next}
    
    for(j in (i+1):nrow(location)){
      if(distance(location$lat[i], location$lat[j], location$long[i], location$long[j]) < 5) {
        location[i,HT_columns] <- location[i, HT_columns] + location[j, HT_columns]
        
        discard <- c(discard, j)
      }
    }
  }
  
  # Remove rows in discard variable
  if(length(discard) > 0){
    location <- location[-discard,]
  }
  
  ################################
  ## PLOTTING PHYLOGENETIC TREE ##
  ################################
  
  ###########################################
  ## PLOTTING REGIONAL HAPLOTYPE ABUNDANCE ##
  ###########################################
  
  # Calculate x-axis limits
  x_min <- min(location$long) - 2
  x_max <- max(location$long) + 2
  
  # Calculate y-axis limits
  y_min <- min(location$lat) - 2
  y_max <- max(location$lat) + 2
  
  # Calculate size range based on the difference between minimum and maximum longitude or latitude
  size_range <- max(x_max - x_min, y_max - y_min)
  
  # Set a constant factor for determining the size of the pie charts
  size_factor <- 0.04
  
  # Calculate the size of the pie charts based on the zoom level
  pie_sizes <- size_range * size_factor
  
  # Retrieve map data for the world
  world <- map_data('world')
  
  p <- ggplot(world, aes(long, lat)) +
    geom_map(map=world, aes(map_id=region), fill='#ECEADE', color = '#B3B1A5') +
    coord_fixed(xlim = c(x_min, x_max), ylim = c(y_min, y_max), ratio = 1)
  
  p <- p + geom_scatterpie(aes(x=long, y=lat, group=region, r = pie_sizes),
                           data=location, cols= HT_names, color= 'black', size = 0.05) +
    scale_fill_manual(values = col_vector[1:length(HT_names)])
  
  # Remove the default background
  p <- p + theme(panel.background = element_blank(), 
                 plot.background = element_blank())
  
  # Add a white background
  p <- p + theme(plot.background = element_rect(fill = "white"))
  
  # Add centered title using theme modification
  p <- p + labs(title = paste("Regional Haplotype Abundance of", Genus_Species)) +
    theme(plot.title = element_text(hjust = 0.5)) 
    
  if(length(HT_names) > 5){
    p <- p + guides(fill = guide_legend(title = "Haplotype", ncol = 2))
  } else {p <- p + guides(fill = guide_legend(title = "Haplotype"))}
  
  ggsave(filename = paste0('~/Downloads/', gsub(pattern = ' ', replacement = '', x = Genus_Species), '.png'), plot = p, dpi = 1000)
  
}

