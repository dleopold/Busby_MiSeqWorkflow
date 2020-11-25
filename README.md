## A simple MiSeq amplicon sequencing project

This repository contains an example MiSeq data set and the associated scripts for creating a reproducable workflow from processing raw reads through data analysis. This intended only as a simple example and not as a comprehensive tutorial. There is no one correct way to process an amplicon sequencing dataset and many good online tutorials are available. For example, [this tutorial](https://astrobiomike.github.io/amplicon/dada2_workflow_ex) by Mike Lee, [this tutorial](https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html) for microbiome data processing using R / Bioconductor, or the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html). 

To explore this example I recommend using RStudio. Once you have cloned this repository (see below), double click the .Rproj file in the project folder. This will open RStudio and set your working directory to the project folder. In the project directory, the `Makefile` controls the workflow by running command line programs or running scripts in the `code` folder. The `Makefile` contains the code to run each script and defines the dependencies among them. For example, if you enter `make denoise` at the terminal, the code to run the denoising script (`code/denoise.R`) will only be executed after ensuring that the earlier `make trim` and `make demux` steps are up to date. This document is not intended to be a GNU Make tutorial, but there are lots of online resources to learn about how Make works and why it is great for reproducable research (e.g., [here](https://the-turing-way.netlify.app/reproducible-research/make.html) & [here](https://kbroman.org/minimal_make/)). For an example of this in practice, check out [this repository](https://github.com/dleopold/SpeciesPoolAge), which contains all of the code and metadata for a recent publication of mine.

____

## Getting started

To clone this repository to your local computer, open a terminal and navigate to the location you want to put the project folder and run:

```
git clone https://github.com/dleopold/Busby_MiSeqWorkflow.git
```

After cloning this repository, the complete example project can be recreated by running `make analysis` from the main project folder (assuming all dependencies are installed, see below). Alternatively, the scripts and comands can be run interactively, using the `Makefile` commands (and this document) as a roadmap. 

____

### Data

The data included in this repository is a subset of a microbiome community sequncing study of *Populus trichocarpa* foliar fungi from 80 trees sampled across two different common garden sites. The Illumina sequencing library was dual-indexed and sequenced using MiSeq v2 2x250 bp chemistry. If you are working with your own data, you will need to download it from the sequencing center and put it in the `data/MiSeq_raw` folder. In the Busby Lab, this will most often require dowloading the data from the CGRB infrastructure, e.g.:

```
scp -P 732 leopoldd@files.cgrb.oregonstate.edu:/nfs2/hts/miseq/190807_M01498_0588_000000000-CJ9GL/Data/Intensities/BaseCalls/nodemultiplex/*.* ./data/MiSeq_raw
```

In the above example, you would need to replace `leopoldd` with your CGRB login and the path to the data `/nfs/hts/miseq...nodemultiplex/*.*` with the path provided by CGRB when the sequencing was completed.

____

### Demultiplexing

This project uses [pheniqs](https://biosails.github.io/pheniqs/) to assign reads to samples based on the index read files (i.e. demultiplexing). If you are working with data that is already demultiplexed you may not need this step. In order for [pheniqs](https://biosails.github.io/pheniqs/) to work you will need to prepare a json configuration file. You should take a look at the example in this project, `code/demux.config.json`, to get an idea how it should be formatted. The [pheniqs](https://biosails.github.io/pheniqs/) documentation is also very comprehensive. For a new project, you will need to modify the `barcode` lines in the configuration file to properly match the sample ids with the indexing barcodes (Note: the I1 reads correspond to the reverse complement of the index on the i7 sequencing adapter and the I2 reads match the multiplexing index adjacent to the i5 adapter in the library construct). In the configuration file, the last block of code specifies the demultiplexing algorithm (see the [pheniqs](https://biosails.github.io/pheniqs/) documenation for details), the expected "noise" (usually the proportion of phiX control added by the sequencing center), and the names of the files to place reads that could not be assigned to samples. The easiest way to create the configuration file is to start with the example provided in this example project and use a good text editor with multicursor editing (e.g., [sublime text](https://www.sublimetext.com/)), so you don't have to enter every line by hand. If you are having trouble getting your configuration file to work, you can use also use an online json linter (eg, [here](https://jsonlint.com/)) to look for json formatting errors. 

Once your configuration file is ready you can run pheniqs:

```
mkdir -p output/processing/demux
pheniqs mux -R output/processing/demux/demux.report.txt \
    -c code/demux.config.json \
    --base-input data/MiSeq_raw \
    --base-output output/processing/demux/
```
or, using the Makefile:
```
make demux
```
A summary of the demultiplexing will be produced, `output/processing/demux/demux.report.json`. The last block of the report contains the overall summary of how many sequences were successfully assigned to samples.

____

### Trimming

The `code/trim.sh` script loops over each pair of fastq files for each demultiplexed sample and first removes the gene primers from the 5' end of the samples using cutadapt. The library constructs in this study include 3-6 degenerate bases (N) that preceed the gene primers, so the cutadapt parameters are set to search for each possible variant anchored to the begining of each read. For example, the gene primer on R1 in this example is ITS4 (TCCTCCGCTTATTGATATGC). So, the cutadapt parameters used will trim any read that begins with one of these variants:

* NNNTCCTCCGCTTATTGATATGC
* NNNNTCCTCCGCTTATTGATATGC
* NNNNNTCCTCCGCTTATTGATATGC
* NNNNNNTCCTCCGCTTATTGATATGC

and will discard any other reads. Reads that pass the first trimming step are then trimmed on the 3' end to remove "read-through" contamination. This occurs when an amplicon is shorter than the number of sequenced bases (ie, < 250 bp) so that the end of the read includes the reverse complement of the opposite direction gene primer. This is done with `SeqPurge`, which uses both the alignment of the paired reads and reverse complemented primers to trim the 3' contamination. For the markers sequenced in this example (ITS2), read-through contamination is not a significant issue because the amplicon is almost always > 250 bp. However, this step will also remove short dimer contamination (< 100 bp). For libraries using different primers you will need to modify the script, replacing the primer sequences (and their reverse compliments) and take a look at the cutadapt and SeqPurge parameters to make sure they are reasonable for your data. 

To help understand what is happening in this script it is helpful to walk through it manually. RStudio has an integrated bash terminal, which makes this easy. Simply open the `code/trim.sh` script and use `ctr+enter` to send the current line to the terminal (or copy-paste). When you get to the loop (line 37, `for FWD in...`) you can manually set the `FWD` variable to one of the forward demultiplexed reads to walk through the loop line-by-line. For eaxmple:

```
FWD=`find ${muxIn} -name "*R1.fastq.gz" | grep -v "undetermined" | head -n1`
```
Then you can walk through the rest of the loop and inspect the output of each step as you go to ensure that things are working as expected. To run the entire script, you can simply run:
```
mkdir -p output/processing/trim
./code/trim.sh
```
or, with the Makefile:
```
make trim
```
____

### Denoising
The `code/denoising.R` script uses DADA2 to infer the "true" sequence variants present in the data set, also known as denoising. This script largely follows the workflow outlined in the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html), which provides a nice walkthrough that I will not reproduce here. One caveate is that markers that vary significantly in length, like fungal ITS sequences, present some unique challenges. Specifically, it is often desirable to truncate reads to remove low quality tails, particularly the reverse reads. The DADA2 function to truncate reads will discard any reads shorter than the truncation length, which is problematic is some real amplicons are shorter than the desired truncation cutoff. To address this issue, the `code/denoise.R` script uses cutadapt to perform read truncation. Setting reasonable truncation thresholds will be project specific and requires looking at the error profiles of some samples. By manually working through the script line-by-line in RStudio, profile plots of base quality will be produced and saved to `output/processing/dada/`. Other diagnostic plots and a summary of the sample processing will also be saved to this directory as the script runs. 

To update the output after modifying the script run:
```
mkdir -p output/processing/dada
Rscript code/denoise.R
```
or, with the Makefile:
```
make denoise
```

____

### Compile
This R script will read in the OTU table produced by the denoising script perform a nuber of post-processing steps:

* Filter non-fungal OTUs from the data set by matching the sequences against the [UNITE](https://unite.ut.ee/repository.php) "all eukaryotes" database using [vsearch](https://github.com/torognes/vsearch) to perform semi-global alignments.
* Trim the representative sequences to remove conserved regions flanking the ITS2 region using [ITSx](https://microbiology.se/software/itsx/). This helps with taxonomic assignment.
* Merge OTU table, trimmed sequences, and sample data into a [phyloseq](https://joey711.github.io/phyloseq/) object for downstream manipulation.
* Cluster similar sequences. In this example I cluster OTUs with > 99% sequence similarity.
* Assign taxonomy using the DADA2 implementation of the [RDP Classifier](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/) and the [UNITE](https://unite.ut.ee/repository.php) fungal species hypothesis database (excluding singletons and taxa only identified at the Kingdom level).
* Save the phyloseq object for downstream use, along the OTU and taxonomy tables as csv files, in `output/processing/compile`.

The specifics of this post-processing stage will be highly project dependent. If you are working with different organism, markers or primer, consult the literature or relevent experts to determine an apropreate post-processing and data curation workflow. 

To run the entire compile script either work through it directly in RStudio, or run:
```
mkdir -p output/processing/compiled/
Rscript code/compile.R
```
or, with the Makefile:
```
make denoise
```

____

### Analysis
The final step of this example project is a simple analysis contained in an RMarkdown file, which is rendered into a well formatted report that is saved in `output/html`. A simple figure is also produced and saved in `output/figs`. For a real project, I recommend organizing your analyses into multiple scripts, each accomplishing a single task or a set of related tasks. For an example of a project with multiple analysis scripts incorporated into a reproducable workflow with GNU Make, see [this repository](https://github.com/dleopold/Populus_priorityEffects). 

To run the analyses manually, open `code/analysis.Rmd` and run each code chunk. To generate the html report, run:

```
make analysis
```

____

## Dependencies
This example project has the following dependencies (If you are working on the Busby Lab workstation everything should already be installed):

* [pheniqs](https://biosails.github.io/pheniqs/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [SeqPurge](https://github.com/imgag/ngs-bits)
* [ITSx](https://microbiology.se/software/itsx/)
* [vsearch](https://github.com/torognes/vsearch)
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

____

If you find any bugs in this example code or have questions about processing your own data, feel free to message me, devin.leopold@gmail.com.

