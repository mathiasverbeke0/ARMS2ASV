#!/usr/bin/env Rscript

############
## MESAGE ##
############

if(basename(download) == download){
  download = normalizePath(file.path(mainpath, download))
}

data <- readLines(download, warn = F)

runs_samples <- list()
current_run <- NULL

for(run in data){
  
  # Remove whitespace characters
  run <- trimws(run)
  run <- gsub(pattern = ' ', replacement = '', run) 
  
  # Determine what the first character is
  firstCharacter = substr(run,1,1)
  
  # If the first character is a >, this is the current run
  if(firstCharacter == '>'){
    current_run <- sub(pattern = '>', replacement = '', run)
    runs_samples[[current_run]] <- c()
  }
  
  else if(run == ''){
    next
  }
  
  else if(is.null(current_run)){
    stop('The ENA input file appears to contain an error.')
  }
  
  else{
    runs_samples[[current_run]] = c(runs_samples[[current_run]], run)
  }
}

for(run in names(runs_samples)){
  
  path.download <- file.path(mainpath, run)
  
  if(!dir.exists(path.download)){
    dir.create(path.download)
  }
  
  for(sample in runs_samples[[run]]){
    
    six_letter_code <- substr(sample, start = 1, stop = 6)
    one_letter_code <- substr(sample, start = nchar(sample), stop = nchar(sample))
    
    fwd_file_name = paste0(sample, '_1.fastq.gz')
    rev_file_name = paste0(sample, '_2.fastq.gz')
    
    url_fwd = paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', six_letter_code, '/', '00', one_letter_code, '/', sample, '/', fwd_file_name)
    url_rev = paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', six_letter_code, '/', '00', one_letter_code, '/', sample, '/', rev_file_name)
    
    dest_fwd <- file.path(path.download, fwd_file_name)
    dest_rev <- file.path(path.download, rev_file_name)
    
    max_retries <- 3
    retry_count <- 0
    download_success <- FALSE
    
    while (!download_success && retry_count < max_retries) {
      retry_count <- retry_count + 1
      tryCatch({
        download.file(url_fwd, dest_fwd)
        download_success <- TRUE
      }, error = function(e) {
        message(paste("Error downloading file, retrying (", retry_count, "/", max_retries, ")"))
        Sys.sleep(10)
      })
    }
    
    if (!download_success) {
      stop("Failed to download file after", max_retries, "attempts")
    }
    
    max_retries <- 3
    retry_count <- 0
    download_success <- FALSE
    
    while (!download_success && retry_count < max_retries) {
      retry_count <- retry_count + 1
      tryCatch({
        download.file(url_rev, dest_rev)
        download_success <- TRUE
      }, error = function(e) {
        message(paste("Error downloading file, retrying (", retry_count, "/", max_retries, ")"))
        Sys.sleep(10)
      })
    }
    
    if (!download_success) {
      stop("Failed to download file after", max_retries, "attempts")
    }
  }
}