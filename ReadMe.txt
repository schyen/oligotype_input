## ReadMe ##

This file contains information on the oligotype_input function.

There is no Execution Script associated with this file, but the function is executed like in any other function:
1) set working directory to where you want to your output file to go
2) source this function (function can be saved in a different working directory)
3) execute this function by using oligotype_input()

Please read below for argument details---------------------------------------------------

function to create a file that is compatible for oligotyping by Meren
using the output files of mothur

input: src.folder = string; path to where mothur output files are stored
        taxon = string; concatinate with c() if specifying more than one; 
                 sets the taxon you want to select; 
                 example: c('Bacteria','Firmicutes')
        output = string; default NULL (does not produce fasta file); 
                 name of output file; .fasta extension required

output: fasta file formated to oligotyping program requirements