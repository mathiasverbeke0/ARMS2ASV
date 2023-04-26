#!/usr/bin/env Rscript

############
## MESAGE ##
############

if(basename(download) == download){
  download = file.path(mainpath, download)
}

data = readLines(download, warn = F)

runs_samples <- list()
current_run <- NULL

for(run in data){
  
  run = gsub(pattern = ' ', replacement = '', run) 
  
  firstCharacter = substr(run,1,1)
  
  if(firstCharacter == '>'){
    current_run <- sub(pattern = '>', replacement = '', run)
    runs_samples[[current_run]] <- c()
  }
  
  else if(run == ''){
    next
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
    
    download.file(url_fwd, dest_fwd)
    download.file(url_rev, dest_rev)
    
  }
}