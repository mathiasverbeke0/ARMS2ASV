#!/usr/bin/env Rscript

#############
## MESSAGE ##
#############

cat(' ____  _____ _____ _____ ____ ___ ____  _____
|  _ \\| ____| ____|  ___|  _ |_ _|  _ \\| ____|
| |_) |  _| |  _| | |_  | |_) | || |_) |  _|   
|  _ <| |___| |___|  _| |  __/| ||  __/| |___  
|_| \\_|_____|_____|_|   |_|  |___|_|   |_____|\n')

####################################
## Parsing command line arguments ##
####################################

suppressWarnings({
  suppressPackageStartupMessages(library(argparse))
})

parser <- ArgumentParser(description = 'Arguments for DADA2 pipeline')

parser$add_argument('-b', '--base_dir', metavar = 'dir', type = 'character', required = TRUE, help = 'Base directory path')
parser$add_argument('-d', '--download', metavar = 'FILENAME', required = FALSE, help = 'Download data from ENA for the given accessions listed in FILENAME')
parser$add_argument('-r', '--run_mode', choices = c('single', 'multi'), required = TRUE, help = 'Specify whether the script should be run once or for multiple runs. Use \'single\' to run the script once, and \'multi\' to run it multiple times.')
parser$add_argument('-u', '--user', type = 'character', required = FALSE, help = 'BOLDSYSTEMS user ID')
parser$add_argument('-P', '--password', type = 'character', required = FALSE, help = 'BOLDSYSTEMS password')
parser$add_argument('-t', '--trim_primers', action = 'store_true', help = 'Specify if you want to trim primers or not. Default is FALSE.')
parser$add_argument('-p', '--primers', nargs=2, type='character', metavar=c('Fwd_Primer', 'Rev_Primer'), help='Forward and reverse primer sequences for trimming.')
parser$add_argument('-s', '--rm_singleton', action = 'store_false', help = 'Turn of singleton removal in the sequence table.')
parser$add_argument('-l', '--trunclen', nargs=2, type='integer', metavar=c('Fwd', 'Rev'), default=c(200,140), help='Set the maximum length for trimmed reads. Default truncation lengths are 200 bases for forward reads and 140 bases for reverse reads.')
parser$add_argument('-m', '--minlen', type='integer', metavar='minlen', default=50, help='Minimum length threshold for the trimmed reads. Default is 50.')
parser$add_argument('-B', '--BOLDigger', action = 'store_true', help = 'Perform taxonomic classification using boldigger.')


args <- parser$parse_args()

# Access the argument values
mainpath <- args$base_dir
download <- args$download
run_mode <- args$run_mode
trim_primers <- args$trim_primers
trunclen <- args$trunclen
primers <- args$primers
minlen <- args$minlen
rm_singleton <- args$rm_singleton
boldigger <- args$BOLDigger
user <- args$user
password <- args$password


########################################
## COMMAND LINE CONDITIONS AND CHECKS ##
########################################

# Mainpath must be existing directory
file_info <- file.info(mainpath)

if(is.na(file_info$isdir)){
  stop(paste(mainpath, 'does not exist.'))
} else if (file_info$isdir == F){
  stop(paste(mainpath, 'is not a directory.'))
}

# Primers must be provided if they must be trimmed
if(trim_primers == T & is.null(primers)){
  cat(paste(parser$format_usage(), '\n'))
  stop('If --trim_primers is specified, --primers must also be specified.')
}

# Taxonomic classification with boldigger requires username and password
if(boldigger == T & (is.null(password) | is.null(user))){
  cat(paste(parser$format_usage(), '\n'))
  stop('If --BOLDigger is specified, --user and --password also need to be specified.')
}


#########################
## hardcoded variables ## 
#########################

dirnames = c('01.Prefiltered', 
             '02.Filtered_Trimmed', 
             '03.Error_Rates', 
             '04.Sample_Inference', 
             '05.Merged_Reads', 
             '06.Seq_Table')


#############################
## DOWNLOAD ENA ACCESSIONS ##
#############################

if(!is.null(download)){
  
  cat('\n[Step 0] Fetching fastq files from ENA\n')
  
  script_path <- commandArgs()[4]
  script_path <- gsub(pattern = '\\\\', replacement = '/', script_path)
  script_path <- gsub(pattern = '--file=', replacement = '', script_path)
  
  # Source the R-script
  source(file.path(dirname(script_path), 'dependencies/ENAFetcher.R'))
  
  # Set run_mode to multi
  run_mode = 'multi'
}


###################################
## SINGLE OR MULTI RUN EXECUTION ##
###################################

if (run_mode == 'single'){
  paths = mainpath
} else if(run_mode == 'multi'){
  paths = list.dirs(path = mainpath, full.names = TRUE, recursive = FALSE)
}


############################
## Downloading R packages ##
############################

cat('\n[Step 1] Installing and loading all packages\n')

pkg <- installed.packages()[,'Package']
ToInstall <- c(
  'argparse',
  'xlsx',
  'dada2',
  'ggplot2',
  'stats', 
  'Biostrings',
  'ShortRead',
  'vegan'
)

for (item in ToInstall){
  if (!item %in% pkg) {
    install.packages(item)
  }
}


###########################################
## LOADING PACKAGES AND SOURCING SCRIPTS ##
###########################################

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(xlsx))


########################################
## ITERATION FOR EVERY SEQUENCING RUN ##
########################################

for(iter in 1:length(paths)){
  
  # Iteration message
  cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter])))
  
  # Iteration label
  label <- paste0('\n[', iter, '/', length(paths), ' - Step')
  
  ##########################
  ## FETCHING FASTQ FILES ##
  ##########################
  
  cat(paste(label, '2] Fetching fastq file paths\n'))
  
  # Forward and reverse read paths
  FwdRead <- sort(list.files(paths[iter], pattern="_1.fastq", full.names = TRUE))
  RevRead <- sort(list.files(paths[iter], pattern="_2.fastq", full.names = TRUE))
  
  # Extract sample names 
  sample.names <- sapply(strsplit(basename(FwdRead), "_"), `[`, 1)
  sample.names.check <- sapply(strsplit(basename(RevRead), "_"), `[`, 1)
  
  # Check if forward and reverse sample names match
  if (length(sample.names)!= length(sample.names.check)){
    
    stop('Amount of forward and reverse files do not match!')
    
  } else {
    
    for (i in 1:length(sample.names)){
      
      if(sample.names[i] != sample.names.check[i]){
        
        stop(paste0('The forward read file of ', 
                    sample.names[i], 
                    ' matched with the reverse read file of ', 
                    sample.names.check[i], '. Something in the main sample folder is wrong and needs to be fixed first.'))
      }
    }
  }
  
  
  ##################################
  ## REMOVE PRIMERS WITH CUTADAPT ##
  ##################################
  
  if (trim_primers == T){
    
    cat(paste(label, '3] Removing primers with cutadapt and prefiltering the reads\n\n'))
    
    # Make a new directory to store prefiltered sequences
    path.cut <- file.path(paths[iter], '01.Prefiltered')
    
    if(!dir.exists(path.cut)){
      dir.create(path.cut)
    }
    
    # Make files to store the prefiltered sequences
    FwdRead.cut <- file.path(path.cut, basename(FwdRead))
    RevRead.cut <- file.path(path.cut, basename(RevRead))
    
    # Define forward and reverse primer AND construct the cutadapt arguments
    FWD <- primers[1] 
    REV <- primers[2]
    
    FWD.argument <- paste0('-g', ' ^', FWD)
    REV.argument <- paste0('-G', ' ^', REV)
    
    # Run cutadapt
    for(i in seq_along(FwdRead)){
      system2('cutadapt', args = c(FWD.argument, # Define the forward read
                                   REV.argument, # Define the reverse read
                                   '-m 1',   # Only keep reads with a minimal length of 1,
                                   # '--discard-untrimmed',
                                   '-o', FwdRead.cut[i], '-p', RevRead.cut[i], # output files
                                   FwdRead[i], RevRead[i])) # input files
    }
  }
  
  ##################################
  ## FILTER READS BASED ON LENGTH ##
  ##################################
  if (trim_primers == F){
    
    cat(paste(label, '3] Prefiltering the reads\n'))
    
    # Make directory path to store prefiltered sequences
    path.cut <- file.path(paths[iter], '01.Prefiltered')
    
    # Make files to store the prefiltered sequences
    FwdRead.cut <- file.path(path.cut, basename(FwdRead))
    RevRead.cut <- file.path(path.cut, basename(RevRead))
    
    length_filtered <- filterAndTrim(FwdRead, FwdRead.cut, RevRead, RevRead.cut, minLen = 1, multithread = TRUE)
    
    # Get the path of the forward and reverse trimmed fastq files
    FwdRead.cut.check <- sort(list.files(path.cut, pattern="_1.fastq", full.names = TRUE))
    RevRead.cut.check <- sort(list.files(path.cut, pattern="_2.fastq", full.names = TRUE))
    
    # Check if forward and reverse files match
    if(length(FwdRead.cut.check) != length(RevRead.cut.check)){
      stop('Something went wrong! Forward and reverse files do not match anymore!\n')
    }
  }
  
  ###################################
  ## INSPECT READ QUALITY PROFILES ##
  ###################################
  
  cat(paste(label, '4] Plotting quality of prefiltered reads\n'))
  
  suppressWarnings({
    
    pdf(file = file.path(path.cut, 'QualityProfilesPreFiltered.pdf'))
    
    # Forward reads
    if(length(FwdRead.cut) <= 10){
      
      plotQualityProfile(FwdRead.cut) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    } else {
      
      plotQualityProfile(FwdRead.cut[1:10]) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    }
    
    # Reverse reads
    if(length(RevRead.cut) <= 10){
      
      plotQualityProfile(RevRead.cut) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    } else {
      
      plotQualityProfile(RevRead.cut[1:10]) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    }
    
    dev.off()
  })
  
  #####################
  ## FILTER AND TRIM ##
  #####################
  
  cat(paste(label, '5] Trimming the prefiltered reads\n'))
  
  # Make directory path to store prefiltered sequences
  path.filt <- file.path(paths[iter], '02.Filtered_Trimmed')
  
  # Make files to store the trimmed sequences
  FwdRead.filt <- file.path(path.filt, paste0(sample.names, "_1_filt.fastq.gz"))
  RevRead.filt <- file.path(path.filt, paste0(sample.names, "_2_filt.fastq.gz"))
  
  names(FwdRead.filt) <- sample.names
  names(RevRead.filt) <- sample.names
  
  # Trim the prefiltered reads
  out <- filterAndTrim(FwdRead.cut,    # Input files forward reads
                       FwdRead.filt,   # Output files forward reads
                       RevRead.cut,    # Input files reverse reads
                       RevRead.filt,   # Output files reverse reads
                       truncLen=trunclen,   # Max length of the reads (F will be truncated to 200, reverse to 140)
                       maxN=0,              # Amount of ambiguous reads that are allowed
                       maxEE=c(2,4),        # Maximum error rates for F and R read
                       truncQ=2,            # Minimum quality score that each base should have
                       minLen = minlen,     # Minimum length of the reads after trimming
                       rm.phix=TRUE,        # Remove contaminant reads
                       compress=TRUE,       # Output files are compressed
                       multithread = T)     # On Windows set multithread = FALSE
  
  saveRDS(out, file.path(path.filt, 'Filtered_Trimmed_Logfile.rds'))
  
  
  ####################################################
  ## INSPECT READ QUALITY PROFILES OF TRIMMED READS ##
  ####################################################
  
  cat(paste(label, '6] Plotting quality of trimmed reads\n'))
  
  suppressWarnings({
    pdf(file = file.path(path.filt, 'QualityProfilesFilteredTrimmed.pdf'))
    
    # Forward reads
    if(length(FwdRead.filt) <= 10){
      
      plotQualityProfile(FwdRead.filt) +
        geom_hline(yintercept = 30)
      
    } else {
      
      plotQualityProfile(FwdRead.filt[1:10]) +
        geom_hline(yintercept = 30)
      
    }
    
    # Reverse reads
    if(length(RevRead.filt) <= 10){
      
      plotQualityProfile(RevRead.filt) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    } else {
      
      plotQualityProfile(RevRead.filt[1:10]) +
        geom_hline(yintercept = 30) +
        scale_color_manual(guide = "none")
      
    }
    
    dev.off()
  })
  
  
  ###########################
  ## LEARN THE ERROR RATES ##
  ###########################
  
  cat(paste(label, '7] Learning error rates\n'))
  
  # Make directory path to store the error plots
  path.error <- file.path(paths[iter], '03.Error_Rates')
  
  if(!dir.exists(path.error)){
    cat(paste('Creating output directory:', path.error,'\n'))
    dir.create(path.error)
  }
  
  # Learn the error rates
  set.seed(100)
  
  
  errF <- learnErrors(FwdRead.filt, multithread=TRUE)
  errR <- learnErrors(RevRead.filt, multithread=TRUE)
  
  
  # Construct and store the error plots
  pdf(file = file.path(path.error, paste0(sample.names[1], '.pdf')))
  plotErrors(errF, nominalQ=TRUE)
  plotErrors(errR, nominalQ=TRUE)
  dev.off()
  
  
  ######################
  ## SAMPLE INFERENCE ##
  ######################
  
  cat(paste(label, '8] Inferring sample\n'))
  
  # Make directory path to store the dada objects
  path.infer <- file.path(paths[iter], '04.Sample_Inference')
  
  if(!dir.exists(path.infer)){
    cat(paste('Creating output directory:', path.infer, '\n'))
    dir.create(path.infer)
  }
  
  # Infer the samples
  
  dadaFwd <- dada(FwdRead.filt, err=errF, multithread=TRUE, pool = "pseudo")
  cat('\n\n')
  dadaRev <- dada(RevRead.filt, err=errR, multithread=TRUE, pool = "pseudo")
  
  
  # Store the dada objects
  saveRDS(dadaFwd, file.path(path.infer, 'dadaFwd.rds'))
  saveRDS(dadaRev, file.path(path.infer, 'dadaRev.rds'))
  
  
  ########################
  ## MERGE PAIRED READS ##
  ########################
  
  cat(paste('\n', label, '9] Merging paired reads\n'))
  
  # Make directory path to store the merger object
  path.merge <- file.path(paths[iter], '05.Merged_Reads')
  
  if(!dir.exists(path.merge)){
    cat(paste('Creating output directory:', path.merge, '\n'))
    dir.create(path.merge)
  }
  
  # Merge the paired reads
  mergers <- mergePairs(dadaFwd, FwdRead.filt, dadaRev, RevRead.filt, minOverlap = 10, maxMismatch = 1, verbose=T)

  # Write mergers object to an rds file
  saveRDS(mergers, file.path(path.merge, 'mergers.rds'))
  
  
  ##############################
  ## CONSTRUCT SEQUENCE TABLE ##
  ##############################
  
  cat(paste(label, '10] Constructing sequence tables\n'))
  
  # Make directory path to store COI ASV sequences
  path.seq <- file.path(paths[iter], '06.Seq_Table')
  
  if(!dir.exists(path.seq)){
    cat(paste('Creating output directory:', path.seq, '\n'))
    dir.create(path.seq)
  }
  
  # Construct the sequence table
  
  seqtab <- makeSequenceTable(mergers)

  
  #####################
  ## REMOVE CHIMERAS ##
  #####################
  
  cat(paste(label, '11] Removing chimeras\n'))
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  # Remove singletons
  if (rm_singleton == T & dim(seqtab.nochim)[1] > 1){
    cat(paste(label, '12] Removing singletons\n'))
    
    mode(seqtab.nochim) = 'numeric'
    seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
  }
  
  # Get ASV sequences
  asv_seqs <- colnames(seqtab.nochim)
  
  # Produce ASV headers
  asv_headers <- vector(dim(seqtab.nochim)[2], mode = 'character')
  
  for (i in 1:dim(seqtab.nochim)[2]){
    asv_headers[i] <- paste0('>', basename(paths[iter]), ';ASV', i)
  }
  
  # Combine ASV sequences and headers
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  
  # write ASV sequences and headers to a text and rds file
  write(asv_fasta, file.path(path.seq, 'COI_ASVS.fasta'))
  saveRDS(seqtab.nochim, file.path(path.seq, 'seqtab.rds'))
}


##############################
## TAXONOMIC CLASSIFICATION ##
##############################

if(boldigger == T){
  cat(paste('\n[12] Assigning taxonomy\n'))
  
  for(iter in 1:length(paths)){
    
    # Getting path to ASV multifasta file
    ASV_paths <- file.path(paths[iter], '06.Seq_Table/COI_ASVS.fasta')
    
    # Construct path to directory where taxonomy is be stored
    path.taxon <- file.path(paths[iter], '07.Taxonomy')
    
    # Create directory where taxonomy is stored
    if(!dir.exists(path.taxon)){
      cat(paste('Creating output directory:', path.taxon, '\n'))
      dir.create(path.taxon)
    }
    
    # Taxonomic classification with BOLDigger
    if(boldigger == T){
      cat('[BOLDigger] ')
      system2(command = 'boldigger-cline', args = c('ie_coi', 
                                                    'mathiasverbeke', 
                                                    'BBD936vjl', 
                                                    ASV_paths[iter], 
                                                    path.taxon))
      
      system2(command = 'boldigger-cline', args = c('first_hit',
                                                    file.path(path.taxon, 'BOLDResults_COI_ASVS_part_1.xlsx')
                                                    )
      )
      
      bold_taxonomy <- read.xlsx(file = file.path(path.taxon, 'BOLDResults_COI_ASVS_part_1.xlsx'), sheetIndex = 'First hit')
      write.xlsx(x = bold_taxonomy, file = file.path(path.taxon, 'BOLD_first_hit.xlsx'))
      file.rename(file.path(path.taxon, 'BOLDResults_COI_ASVS_part_1.xlsx'), file.path(path.taxon, 'BOLD_all_hits.xlsx'))
    }
  }
}