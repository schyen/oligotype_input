#///////////////////////////////////////////////////////////////////////////////////
# function to create a file that is compatible for oligotyping by Meren
# using the output files of mothur

# input: src.folder = string; path to where mothur output files are stored
#        taxon = string; concatinate with c() if specifying more than one; 
#                 sets the taxon you want to select; 
#                 example: c('Bacteria','Firmicutes')
#        output = string; default NULL (does not produce fasta file); 
#                 name of output file; .fasta extension required

# output: fasta file formated to oligotyping program requirements

oligotype_input <- function(src.folder, taxon=NULL, output=NULL){
  
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
  taxon_lev_num <- length(taxon)
  taxon_lev <- tax_name[taxon_lev_num]
  
  tax_sub <- tax_df %>%
    filter(grepl(taxon[taxon_lev_num], tax_df[,taxon_lev]))
  

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
  
  # going through one sample at a time
  for(i in sampleid) {
    
    # going through one row at a time
    for(j in 1:nrow(count_sub)) {
      nrep <- count_sub[j,i]
      
      # repeating sequences according to number of reads
      sequ <- count_sub[j,'seq']
      entry <- as.list(rep(sequ, count_sub[j,i]))
      
      samp_seq <- c(samp_seq, entry)

      # creating sequence headers
      read <- str_c('Read', 1:nrep, sep='')
      entry <- str_c(i,read, sep='_')
      
      samp_header <- c(samp_header, entry)
    }
    sequences <- c(sequences, samp_seq)
    header <- c(header, samp_header)
  }

  # writing file to fasta
  if(!is.null(output)) {
    write.fasta(sequences=sequences, names=header, file.out=output)
  }
  
  
}