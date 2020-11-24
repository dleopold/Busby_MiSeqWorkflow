## A simple example MiSeq amplicon sequencing project

This repository contains an example MiSeq data set and the associated scripts for creating a reproducable workflow from processing raw reads through data analysis. In the project directory the `Makefile` controls the workflow by running scripts in the `code` folder. After cloning this repository, the complete example project can be recreated by running `make analysis` from the main project folder in a terminal (assuming all dependencies are installed, see below). Alternatively, the scripts and comands can be run interactively, using the `makefile` commands as a roadmap. 

The data included in this repository consists fungal ITS2 amplicons from *Populus trichocarpa* leaves from two different sites. The data and code provided are intended only as a simple example and not as a comprehensive tutorial. There is no one correct way to process an amplicon sequencing dataset and many good online tutorials are available. For example, [this tutorial](https://astrobiomike.github.io/amplicon/dada2_workflow_ex) by Mike Lee, [this tutorial](https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html) for microbiome data processing using R / Bioconductor, or the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html). 

To explore this example I recommend using RStudio. Once you have Rstudio on your computer and have cloned this repository double click the .Rproj file in this folder. This will open RStudio and set your working directory to the project folder. You can then explore the project folders and files using files tab in the lower right pane of RStudio (assuming you have not changed the defaults). 

This example project has the following dependencies (If you are working on the Busby Lab workstation everything should be ready to go):

* [pheniqs](https://biosails.github.io/pheniqs/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [SeqPurge](https://github.com/imgag/ngs-bits)
* [ITSx](https://microbiology.se/software/itsx/)
* [R](https://www.r-project.org/)
* R-packages:  
    * tidyverse  
    * magrittr  
    * dada2  
    * ShortRead  
    * phyloseq  
    * Biostrings  
    * vegan  
    * DECHIPHER  
    * speedyseq  

