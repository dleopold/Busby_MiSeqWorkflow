## Notes

This repository contains an example MiSeq data set and the associated scripts for creating a reproducable workflow from processing raw reads through data analysis. This example uses GNU Make for workflow management. Basically, this means that the Makefile in the main project directory runs all of the scripts in the code folder in the proper order ensuring that everything is up to date.  

To explore the scripts I recommend using RStudio. Once you have Rstudio on your computer double click the .Rproj file in this folder. This will open RStudio and set your working directory to the project folder. You can then explore the project folders and files using files tab in the lower right pane of RStudio (assuming you have not changed the defaults). 


Assuming you have all of the dependencies (see below), all of the output can be deleted and recreated by opening a terminal (in lower left pane of RStudio as long as you have a relatively recent version) and typing "make clean" (this will delete the output folder) and then "make analysis" (this will rebuild everything). You should have a minimum of 8 gb ram.

This example project has the following dependencies (If you are working on the Busby Lab workstation everything should be ready to go)

Programs:
pheniqs https://biosails.github.io/pheniqs/
cutadapt https://cutadapt.readthedocs.io/en/stable/
SeqPurge https://github.com/imgag/ngs-bits
ITSx https://microbiology.se/software/itsx/

R packages:
tidyverse
magrittr
dada2
ShortRead
phyloseq
Biostrings
foreach
vegan
DECHIPHER
speedyseq

