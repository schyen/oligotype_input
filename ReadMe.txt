## ReadMe ##

This file contains information on the oligotype_input function.

There is no Execution Script associated with this file, but the function is executed like in any other function:
1) set working directory to where you want to your output file to go
2) source this function (function can be saved in a different working directory)
3) execute this function by using oligotype_input()

Please read below for argument details---------------------------------------------------

function to create a file that is compatible for oligotyping by Meren
using the output files of mothur

input: src.folder = string; default is current working directory;
                     path to where mothur output files are stored
        taxon = vector of strings or list of string vectors; default NULL;
                concatinate with c() if specifying more than one; 
                 sets the taxon you want to select; 
                 example: c('Bacteria','Firmicutes');
                if want to process multiple taxonomic classes 
                (and return multiple fasta files), provide taxon as list of vectors
                example: list(c('Bacteria','Phylum1','Class1','Order1'),
                              c('Bacteria','Phylum1','Class1','Order2'))
        sample_name = vector of strings; default NULL;
                     name of samples to be included in the fasta file
        output = string or vector of strings; 
                 default NULL (does not produce fasta file); 
                 name of output file; .fasta extension required;
                 when multiple taxon filters specified, 
                 need to provide file name for each one 

# output: fasta file(s) formated to oligotyping program requirements