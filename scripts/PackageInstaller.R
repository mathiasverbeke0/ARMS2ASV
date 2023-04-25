pkg <- installed.packages()[,'Package']
ToInstall <- c(
  'argparse',
  'dada2',
  'ggplot2',
  'stats', 
  'Biostrings',
  'ShortRead'
)

for (item in ToInstall){
  if (!item %in% pkg) {
    install.packages(item)
  }
}



