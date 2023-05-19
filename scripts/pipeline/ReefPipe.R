#!/usr/bin/env Rscript

#############
## MESSAGE ##
#############

cat(' ____  _____ _____ _____ ____ ___ ____  _____
|  _ \\| ____| ____|  ___|  _ |_ _|  _ \\| ____|
| |_) |  _| |  _| | |_  | |_) | || |_) |  _|   
|  _ <| |___| |___|  _| |  __/| ||  __/| |___  
|_| \\_|_____|_____|_|   |_|  |___|_|   |_____|\n\n')


#####################################
## DOWNLOADING AND LOAD R PACKAGES ##
#####################################

cat('\n _ _ _     
| (_) |__  
| | | \'_ \\ 
| | | |_) |
|_|_|_.__/ \n')

cat('\nInstalling and loading all packages\n')

pkg <- installed.packages()[,'Package']

# Specify all packages
ToInstall <- c(
  'xlsx',
  'dada2',
  'ggplot2',
  'stats', 
  'Biostrings',
  'ShortRead',
  'vegan',
  'readxl',
  'stringr',
  'argparse',
  'purrr',
  'dplyr',
  'gtools'
)

for(item in ToInstall){
  
  # If not installed, install the package
  if(!item %in% pkg) {
    install.packages(item)
  }
  
  # Load the package (silently)
  suppressPackageStartupMessages(library(item, character.only = T))
}


####################################
## Parsing command line arguments ##
####################################

parser <- ArgumentParser(description = 'Reefpipe command line arguments')

# Mandatory command line arguments
parser$add_argument('-b', '--base_dir', metavar = 'BASEDIR', type = 'character', required = TRUE, help = 'The base directory path for the analysis.')
parser$add_argument('-d', '--download', metavar = 'FILENAME', required = FALSE, help = 'Download data from ENA for the given accessions listed in FILENAME.')
parser$add_argument('-r', '--run_mode', choices = c('single', 'multi'), required = TRUE, help = 'Specify whether to run the script for a single or multiple sequencing runs. Select \'single\' for analyzing a single run and \'multi\' for analyzing multiple runs.')

# Command line arguments for preprocessing
parser$add_argument('-t', '--trim_primers', action = 'store_true', help = 'Trim primers with Cutadapt. By default, primers are not trimmed.')
parser$add_argument('-p', '--primers', nargs=2, type='character', metavar=c('Fwd_Primer', 'Rev_Primer'), help='Forward and reverse primer sequences for trimming.')
parser$add_argument('-s', '--singleton', action = 'store_false', help = 'Keep singletons in the sequence table. By default, singletons are removed (except when only one sample is analyzed).')
parser$add_argument('-m', '--minlen', type='integer', metavar='minlen', default=50, help='The minimum length threshold for the trimmed reads. Default is 50.')
parser$add_argument('-l', '--trunclen', nargs=2, type='integer', metavar=c('Fwd', 'Rev'), default=c(200,140), help='The maximum length for trimmed reads. Default truncation lengths are 200 bases for forward reads and 140 bases for reverse reads.')

# Command line arguments for taxonomic classification
parser$add_argument('-B', '--BOLDigger', action = 'store_true', help = 'Perform taxonomic classification using BOLDigger.')
parser$add_argument('-U', '--user', type = 'character', required = FALSE, help = 'The BOLDSYSTEMS user ID.')
parser$add_argument('-P', '--password', type = 'character', required = FALSE, help = 'The BOLDSYSTEMS password.')
parser$add_argument('-R', '--reference', action = 'store_true', help = 'Perform taxonomic classification with DADA2 using own reference databases.')
parser$add_argument('-M', '--minBoot', type = 'numeric', default = 80, help = 'The minimal bootstrap value for taxonomic classification with DADA2. Default is 80.')

# Command line arguments for taxonomic table fusing
parser$add_argument('-F', '--fuse', action = 'store_true', help = 'Fuse the information of all taxonomic tables.')
parser$add_argument('-f', '--fuseLevels', type = 'character', default = 'Phylum,Class,Order,Family,Genus,Species', help = 'The taxonomic levels used for fusing all taxonomic tables. Default levels are Phylum,Class,Order,Family,Genus,Species.')

# Parse the arguments
args <- parser$parse_args()

# Access the argument values
mainpath <- normalizePath(args$base_dir)
download <- args$download
run_mode <- args$run_mode
trim_primers <- args$trim_primers
primers <- args$primers
minlen <- args$minlen
trunclen <- args$trunclen
singleton <- args$singleton
boldigger <- args$BOLDigger
user <- args$user
password <- args$password
reference <- args$reference
minBoot <- args$minBoot
fuse <- args$fuse
fuseLevels <- args$fuseLevels


####################################
## ABSOLUTE PATH OF MAIN PIPELINE ##
####################################

args <- commandArgs()
pipeline_path <- NULL
for (arg in args) {
  if (startsWith(arg, "--file=")) {
    pipeline_path <- sub("^--file=", "", arg)
    break
  }
}

if (!is.null(pipeline_path)) {
  pipeline_path <- normalizePath(pipeline_path)
  
} else {
  stop("Pipeline path not found.")
}


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

# Fusing taxonomic tables cannot be done if only the --BOLDigger option is selected
if(boldigger == T & fuse == T & reference == F){
  cat('Warning: you cannot merge taxonomic tables if only the BOLDSYSTEMS table is generated.')
  
  # Change fuse to false
  fuse = F
}

# If --reference is selected, it must contain reference databases
if(reference == T){
  # Get a list of all reference database files in the reference database directory
  ref_files <- list.files(path = file.path(dirname(pipeline_path), '../../data/reference/'), full.names = T)
  
  # Get the configuration file
  config_file <- file.path(dirname(pipeline_path), '../../data/reference/config.txt')
  
  # Check if the configuration file exists
  if(!file.exists(config_file)){
    stop(paste("The file", config_file,  "does not exist, which is necessary for the process of taxonomic classification using locally stored reference databases. It appears that the file may have been (re)moved or is missing from its expected location."))
  }
  
  # Remove config file from ref_files
  ref_files <- ref_files[!basename(ref_files) == 'config.txt']
  
  if(length(ref_files) == 0){
    stop('The reference database folder is empty. Please make sure there are reference databases available in the directory.')
  }
}

# Fusing taxonomic tables requires all reference databases to have the levels specified with --fuseLevels
if(fuse == TRUE){
  source(file.path(dirname(pipeline_path), 'dependencies/taxLevelCheck.R'))
}


#############################
## DOWNLOAD ENA ACCESSIONS ##
#############################

if(!is.null(download)){
  cat('\n  ___ _ __   __ _ 
 / _ \\ \'_ \\ / _` |
|  __/ | | | (_| |
 \\___|_| |_|\\__,_|\n')
  
  cat('\nFetching fastq files from ENA\n\n')
  
  # Source the R-script
  source(file.path(dirname(pipeline_path), 'dependencies/ENAFetcher.R'))
  
  # Set run_mode to multi
  run_mode = 'multi'
}


###################################
## SINGLE OR MULTI RUN EXECUTION ##
###################################

if (run_mode == 'single' & is.null(download)){
  paths = mainpath
} else if(run_mode == 'multi'){
  paths = list.dirs(path = mainpath, full.names = TRUE, recursive = FALSE)
}


##############################
## DIRECTORY CONTENTS CHECK ##
##############################

# Mainpath must have contents
if(length(list.files(mainpath)) == 0){
  cat('\n')
  stop('There are no contents in the base directory.')
}

# If run_mode is multi, the child directories should only contain _1.fastq and _2.fastq files
if(run_mode == 'multi'){
  for(path in paths){
    files <- list.files(path = path, full.names = TRUE)
    is_match <- grepl(pattern = "_1.fastq|_2.fastq", x = files) 
    
    if (any(!is_match)){
      cat('\n')
      stop(paste(path, 
                 'includes files other than being _1.fastq or _2.fastq files.', 
                 '\nPlease ensure that the base directory you are using only contains directories with exclusively fastq files if you are using the multi-option or downloading fastq files from ENA through this pipeline.',
                 '\n\nExample:',
                 '\nbase_directory/',
                 '\n├── ENA_accessions.txt',
                 '\n├── sample1/',
                 '\n│   ├── sample1_1.fastq',
                 '\n│   └── sample1_2.fastq',
                 '\n├── sample2/',
                 '\n│   ├── sample2_1.fastq',
                 '\n│   └── sample2_2.fastq',
                 '\n└── sample3/',
                 '\n    ├── sample3_1.fastq',
                 '\n    └── sample3_2.fastq\n\n'))
    }
  }
}


#############################################
## ASV GENERATION FOR EVERY SEQUENCING RUN ##
#############################################

cat('\n  __ _ _____   __
 / _` / __\\ \\ / /
| (_| \\__ \\\\ V / 
 \\__,_|___/ \\_/  \n')

for(iter in 1:length(paths)){
  
  ##################################
  ## ITERATION AND STEP VARIABLES ##
  ##################################
  
  # Iteration message
  cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter])))
  
  # Iteration label
  label <- paste0('\n[', iter, '/', length(paths), ' - Step ')
  
  # Step counter
  step <- 0
  
  ##########################
  ## FETCHING FASTQ FILES ##
  ##########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Fetching fastq file paths\n'))
  
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
    
    # Message
    step <- step + 1
    cat(paste0(label, step, '] Removing primers with cutadapt and prefiltering the reads\n'))
    
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
    
    # Message
    step <- step + 1
    cat(paste0(label, step, '] Prefiltering the reads\n'))
    
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Plotting quality of prefiltered reads\n'))
  
  # Make PDF with read quality profiles
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Trimming the prefiltered reads\n'))
  
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Plotting quality of trimmed reads\n'))
  
  # Make PDF with read quality profiles
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Learning error rates\n'))
  
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Inferring sample\n'))
  
  # Make directory path to store the dada objects
  path.infer <- file.path(paths[iter], '04.Sample_Inference')
  
  if(!dir.exists(path.infer)){
    cat(paste('Creating output directory:', path.infer, '\n'))
    dir.create(path.infer)
  }
  
  # Infer the samples
  dadaFwd <- dada(FwdRead.filt, err=errF, multithread=TRUE, pool = "pseudo")
  cat('\n')
  dadaRev <- dada(RevRead.filt, err=errR, multithread=TRUE, pool = "pseudo")
  cat('\n')
  
  # Store the dada objects
  saveRDS(dadaFwd, file.path(path.infer, 'dadaFwd.rds'))
  saveRDS(dadaRev, file.path(path.infer, 'dadaRev.rds'))
  
  
  ########################
  ## MERGE PAIRED READS ##
  ########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Merging paired reads\n'))
  
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
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Constructing sequence tables\n'))
  
  # Make directory path to store COI ASV sequences
  path.seq <- file.path(paths[iter], '06.Seq_Table')
  
  if(!dir.exists(path.seq)){
    cat(paste('Creating output directory:', path.seq, '\n'))
    dir.create(path.seq)
  }
  
  # Construct the sequence table
  
  seqtab <- makeSequenceTable(mergers)
  
  ##########################
  ## CHECK SEQUENCE TABLE ##
  ##########################
  
  # If there are no sequences in the sequence table, continue
  if(dim(seqtab)[2] == 0){
    
    cat(paste('\nASVS could not be generated for', basename(paths[iter]), 
              '\nExcluding', paths[iter], 'from paths to take into consideration.\n'))
    
    # Write failed ASV generation to output file
    if(!file.exists(file.path(mainpath, 'log.txt'))){
      file.create(file.path(mainpath, 'log.txt'))
    }
    
    # Find index of path to remove
    index_to_remove <- grep(pattern = paths[iter], paths)
    
    # Remove path from paths variable
    paths <- paths[-index_to_remove]
    
    # Continue
    next
  }
  
  # If there is only one row in the sequence table, add the sample name manually
  if(dim(seqtab)[1] == 1){
    rownames(seqtab) <- sample.names
    cat('test')
  }
  
  #####################
  ## REMOVE CHIMERAS ##
  #####################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Removing chimeras\n'))
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  # Remove singletons
  if (singleton == T & dim(seqtab.nochim)[1] > 1){
    
    # Message
    step <- step + 1
    cat(paste0(label, step, '] Removing singletons\n'))
    
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
  
  
  #####################
  ## ASV RAREFACTION ##
  #####################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] plotting ASV rarefaction\n'))
  
  if(nrow(seqtab.nochim) > 1){
    
    # Make directory path to store ecological anlysis results
    path.eco <- file.path(paths[iter], '08.EcologicalAnalysis')
    
    if(!dir.exists(path.eco)){
      cat(paste('Creating output directory:', path.eco, '\n'))
      dir.create(path.eco)
    }
    
    # Execute ASV rarefaction script
    source(file.path(dirname(pipeline_path), 'dependencies/Rarefaction.R'))
  } 
  
  else{
    cat('Rarefaction not possible due to the sequence table having only 1 row.\n')
  }
  
  cat('\n\n')
}


##############################
## TAXONOMIC CLASSIFICATION ##
##############################

# Check if taxonomic classification needs to be performed
if((reference == T | boldigger == T) & length(paths) > 0){
  
  # Taxonomy message
  cat(' _             
| |_ __ ___  __
| __/ _` \\ \\/ /
| || (_| |>  < 
 \\__\\__,_/_/\\_\\\n\n')
  
  # Get paths to ASV multifasta files
  paths.ASV <- file.path(paths, '06.Seq_Table/COI_ASVS.fasta')
  
  # Get paths to transformed ASV multifasta files (BOLDigger renames the used multifasta files to COI_ASVS_done.fasta)
    boldigger.paths.ASV <- file.path(paths, '06.Seq_Table/COI_ASVS_done.fasta')
  
  # Construct paths to directories where taxonomy files will be stored
  paths.taxon <- file.path(paths, '07.Taxonomy')
  
  # Create directory where taxonomy is stored
  for(path.taxon in paths.taxon){
    if(!dir.exists(path.taxon)){
      cat(paste('Creating output directories:', path.taxon, '\n'))
      dir.create(path.taxon)
    }
  }
  
  
  #############################################################
  ## TAXONOMIC CLASSIFICATION WITH LOCAL REFERENCE DATABASES ##
  #############################################################
  
  if(reference == T){
    cat('\n[Reference libraries]')
  

    ########################################
    ## ITERATION FOR EVERY SEQUENCING RUN ##
    ########################################
    
    for(iter in 1:length(paths)){
      
      # Print iteration message
      cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter]), '\n'))
      
      # Locations of the ASV multifasta file and taxonomy directory
      path.taxon <- paths.taxon[iter]
      path.ASV <- paths.ASV[iter]
      
      # Execute taxonomic classification with DADA2
      source(file.path(dirname(pipeline_path), 'dependencies/TaxonomicClassification.R'))
    }
  }
  
  
  ##############################################
  ## TAXONOMIC CLASSIFICATION WITH BOLDIGGER ##
  #############################################
  
  if(boldigger == T){
    cat('\n\n[BOLDigger]')
    
    #########################################
    ## ITERATION OVER EVERY SEQUENCING RUN ##
    #########################################
    
    for(iter in 1:length(paths)){
      
      # Print iteration message
      cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter]), '\n'))
      
      # Locations of the ASV multifasta file and taxonomy directory
      path.taxon <- paths.taxon[iter]
      path.ASV <- paths.ASV[iter]
      boldigger.path.ASV <- boldigger.paths.ASV[iter]
      
      # Execute the BOLDigger command line tool: find top 20 hits
      system2(command = 'boldigger-cline', args = c('ie_coi', 
                                                    user, 
                                                    password, 
                                                    path.ASV, 
                                                    path.taxon))
      
      # Get the file names in the Taxonomy directory
      files <- list.files(path = path.taxon, full.names = T)
      
      non_dirs <- c()
      
      for(file in files){
        if(dir.exists(paths = file)){next} 
        else{non_dirs <- c(non_dirs, file)}
      }
      
      # Search for the BOLDigger output Excel file using regular expressions
      excel_file <- files[grepl(".xlsx$", files)]
      
      # Execute the BOLDigger command line tool: find first hit (top-scoring match)
      system2(command = 'boldigger-cline', args = c('first_hit', excel_file))
      
      # Read in the second sheet of the BOLDigger output excel file
      bold_taxonomy <- read.xlsx(file = file.path(path.taxon, 'BOLDResults_COI_ASVS_part_1.xlsx'), sheetIndex = 'First hit')
      
      # Create directory to store only the first hits
      if(!dir.exists(file.path(path.taxon, 'BOLDSYSTEMS'))){
        cat(paste('Saving first hit outputs to:', file.path(path.taxon, 'BOLDSYSTEMS')), '\n')
        dir.create(file.path(path.taxon, 'BOLDSYSTEMS'))
      }
      
      # Write the contents of the second sheet to a separate excel file
      write.xlsx(x = bold_taxonomy, file = file.path(path.taxon, 'BOLDSYSTEMS', 'BOLD_first_hit.xlsx'))
      
      # Remove the original BOLDigger output files
      for(non_dir in non_dirs){unlink(x = non_dir)}
      
      file.rename(from = boldigger.path.ASV, to = path.ASV)
    }
  }
  
  if(fuse){
    cat('\n[Merging]')
    
    for(iter in 1:length(paths)){
      
      # Print iteration message
      cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter]), '\n'))
      
      # Locations of the taxonomy directory
      path.taxon <- paths.taxon[iter]
      
      # Source the taxonomic table merging script
      source(file.path(dirname(pipeline_path), 'dependencies/TaxTableMerger.R'))
    }
  }
}