#///////////////////////////////////////////////////////////////////////////////////
# function to create a file that is compatible for oligotyping by Meren
# using the output files of mothur

# input: src.folder = string; default is current working directory;
#                     path to where mothur output files are stored
#        taxon = vector of strings or list of string vectors; default NULL;
#                concatinate with c() if specifying more than one; 
#                 sets the taxon you want to select; 
#                 example: c('Bacteria','Firmicutes');
#                if want to process multiple taxonomic classes 
#                (and return multiple fasta files), provide taxon as list of vectors
#                example: list(c('Bacteria','Phylum1','Class1','Order1'),
#                              c('Bacteria','Phylum1','Class1','Order2'))
#        sample_name = vector of strings; default NULL;
#                     name of samples to be included in the fasta file
#        output = string or vector of strings; 
#                 default NULL (does not produce fasta file); 
#                 name of output file; .fasta extension required;
#                 when multiple taxon filters specified, 
#                 need to provide file name for each one 

# output: fasta file(s) formated to oligotyping program requirements

oligotype_input <- function(src.folder=getwd(), taxon=NULL, sample_name=NULL, 
                            output=NULL){
  
  # function requirements------------------------------------------------------------
  #checking for required packages; installing those not yet installed
  if(require(seqinr)==FALSE) install.packages('seqinr')
  if(require(stringr)==FALSE) install.packages('stringr')
  if(require(dplyr)==FALSE) install.packages('dplyr')
  
  #loading required packages
  library(seqinr)
  require(stringr)
  require(dplyr)
  
  # getting all file names in src.folder---------------------------------------------
  files <- list.files(path=src.folder)
  
  # getting file names of required files
  count.fname <- str_extract(files, '.*(\\.trim\\.contigs\\.good\\.unique\\.good\\.filter\\.unique\\.precluster\\.denovo\\.uchime\\.pick.count_table)')
  count.fname <- count.fname[!is.na(count.fname)]
  
  tax.fname <- str_extract(files, '.*(\\.trim\\.contigs\\.good\\.unique\\.good\\.filter\\.unique\\.precluster\\.pick\\.wang\\.wang\\.taxonomy)')
  tax.fname <- tax.fname[!is.na(tax.fname)]
  
  fasta.fname <- str_extract(files, '.*(\\.trim\\.contigs\\.good\\.unique\\.good\\.filter\\.unique\\.precluster\\.pick\\.fasta)')
  fasta.fname <- fasta.fname[!is.na(fasta.fname)]
  
  # reading count_table, taxonomy and fasta dna files--------------------------------
  count <- read.table(file.path(src.folder, count.fname), header=TRUE)
  tax <- read.table(file.path(src.folder, tax.fname))
  fasta <- read.table(file.path(src.folder, fasta.fname))
  
  # formating files into database-like dataframes------------------------------------
  ## count table
  count_df <- count
  colnames(count_df)[1] <- 'seqid'
  count_df$seqid <- as.character(count_df$seqid)
  
  # selecting count for specified sample_names
  if(!is.null(sample_name)) {
    count_df <- count_df[colnames(count_df) %in% c('seqid', 'total', sample_name)]
    msg <- sprintf("Only returning results for samples %s", 
                   str_c(sample_name, collapse=', '))
    message(msg)
  }
  
  ## taxonomy table
  # first, separating by ;
  tax_parse <- str_split_fixed(tax[,2], ';', n=8)

  # removing numbers in brackets
  taxonomy <- c()
  for(i in 1:nrow(tax_parse)) {
    entry <- str_extract_all(tax_parse[i,], '(.*(?=\\())|([a-z]+)', 
                               simplify=TRUE)
    entry <- t(entry[,1])
    
    taxonomy <- rbind(taxonomy, entry)
  }


  # dropping last column
  taxonomy <- taxonomy[,-8]
  
  # assigning column names
  tax_name <- c('kingdom','phylum','class','order','family',
                'genus','species')
  colnames(taxonomy) <- tax_name
  
  # putting parsed taxonomy with respective sequence id
  tax_df <- data.frame(seqid=as.character(tax[,1]))
  tax_df <- cbind(tax_df, taxonomy)

  ## fasta table
  id_index <- seq(1,nrow(fasta),2)
  seq_index <- seq(2, nrow(fasta),2)

  seqid <- str_extract_all(fasta[id_index,], '(?<=\\>).*', simplify=TRUE)

  fasta_df <- data.frame(seqid=seqid, seq=fasta[seq_index,])
  
  # filtering for specified taxon----------------------------------------------------
  # check if multiple taxons provided
  logic <- is.list(taxon)
  
  # when only one taxon filter provided, convert to list
  if(!logic) taxon <- list(taxon)
  
  # when multiple taxons provided, go through one taxon at a time
  for(k in 1:length(taxon)) {
    taxon_lev_num <- length(taxon[[k]])
    taxon_lev <- tax_name[taxon_lev_num]
    
    # returning comment on progress
    message(sprintf("Filtering for %s", taxon[[k]][taxon_lev_num]))
    
    tax_sub <- tax_df %>%
      filter(grepl(taxon[[k]][taxon_lev_num], tax_df[,taxon_lev]))
    
    # Retrieving sequence for taxon of interest and respective counts
    fasta_sub <- left_join(tax_sub, fasta_df, by='seqid')
    fasta_sub$seqid <- as.character(fasta_sub$seqid)
    
    count_sub <- left_join(fasta_sub, count_df, by='seqid')
  
    # repeating sequences according to the number reads and 
    # creating sequence headers for output fasta file----------------------------------
    sampleid <- colnames(count_df)[!colnames(count_df) %in% c('seqid','total')]
  
    samp_header <- c() # headers for each sequence in a sample
    header <- c() # headers for all samples
    
    samp_seq <- list() # sequences for each sequence in a sample
    sequences <- list() # sequences for all samples
    
    out <- c() # output for each sequence in sample
    out_all <- c() # output for all samples
    
    # going through one sample at a time
    for(i in sampleid) {
      
      # going through one row at a time
      for(j in 1:nrow(count_sub)) {
        
        # repeating sequences according to number of reads
        nrep <- count_sub[j,i]
        sequ <- count_sub[j,'seq']
        
        # repeating the sequence the number of times that sequence is found;
        # ignore sequence if sequence read 0 times
        if(nrep != 0) {
          
          seq_rep <- rep(sequ, count_sub[j,i])
          
          # creating sequence headers
          read <- str_c('Read', 1:nrep, sep='')
          header <- str_c(i, read, sep='_')
          
          # adding repeated sequences and header as entry into output
          entry_df <- data.frame(header=header, seq=seq_rep)
          
          out <- rbind(out, entry_df)    
        }
      }
      
      out_all <- rbind(out_all, out)
      
    }

    # returning sequences as list format to write as fasta
    sequences <- as.list(out_all$seq)
    header <- out_all$header
    
    # writing file to fasta
    if(!is.null(output)) {
      
      # check output file name supplied for every taxon filter
      check <- length(output) == length(taxon)
      if(!check) {
        stop('A filename must be provided for each fasta output', call.=FALSE)
      }
      else{
        write.fasta(sequences=sequences, names=header, file.out=output[k])
      }
    }
    
  }
  
#   taxon_lev_num <- length(taxon)
#   taxon_lev <- tax_name[taxon_lev_num]
#   
#   tax_sub <- tax_df %>%
#     filter(grepl(taxon[taxon_lev_num], tax_df[,taxon_lev]))
#   
#   # Retrieving sequence for taxon of interest and respective counts
#   fasta_sub <- left_join(tax_sub, fasta_df, by='seqid')
#   fasta_sub$seqid <- as.character(fasta_sub$seqid)
#   
#   count_sub <- left_join(fasta_sub, count_df, by='seqid')
#   
#   # repeating sequences according to the number reads and 
#   # creating sequence headers for output fasta file----------------------------------
#   sampleid <- colnames(count_df)[!colnames(count_df) %in% c('seqid','total')]
#   
#   samp_header <- c() # headers for each sequence in a sample
#   header <- c() # headers for all samples
#   
#   samp_seq <- list() # sequences for each sequence in a sample
#   sequences <- list() # sequences for all samples
# 
#   out <- c() # output for each sequence in sample
#   out_all <- c() # output for all samples
#   
#   # going through one sample at a time
#   for(i in sampleid) {
# 
#     # going through one row at a time
#     for(j in 1:nrow(count_sub)) {
# 
#       # repeating sequences according to number of reads
#       nrep <- count_sub[j,i]
#       sequ <- count_sub[j,'seq']
#       
#       # repeating the sequence the number of times that sequence is found;
#       # ignore sequence if sequence read 0 times
#       if(nrep != 0) {
#           
#         seq_rep <- rep(sequ, count_sub[j,i])
#         
#         # creating sequence headers
#         read <- str_c('Read', 1:nrep, sep='')
#         header <- str_c(i, read, sep='_')
#         
#         # adding repeated sequences and header as entry into output
#         entry_df <- data.frame(header=header, seq=seq_rep)
#     
#         out <- rbind(out, entry_df)    
#         }
#     }
# 
#     out_all <- rbind(out_all, out)
#     
#   }
#   
# # returning sequences as list format to write as fasta
# sequences <- as.list(out_all$seq)
# header <- out_all$header
# 
#   # writing file to fasta
#   if(!is.null(output)) {
#     write.fasta(sequences=sequences, names=header, file.out=output)
#   }
  
}