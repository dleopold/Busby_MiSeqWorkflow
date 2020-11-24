#!/bin/bash

#This script will remove gene primers from both ends of paired-end read to prepare for denoising
#requires cutadapt https://cutadapt.readthedocs.io/en/stable/
#requires SeqPurge https://github.com/imgag/ngs-bits

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#Steps marked with alternating hash/dash may require adjustment for each data set

# This script can take 2 inputs, the path to the demux files and the output path for the final trimmed files.
# Default values are provided.
muxIn=${1:-output/processing/demux}
outDir=${2:-output/processing/trim}

#make directory for temporary files that will get deleted later
scratch="output/scratch/"
mkdir -p $scratch

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#set gene primers for trimming
#In this example the library prep merged the reverse gene primer (ITS4) with the forward sequencing primer.
#As a result, ITS4 is on the forward read (we will need to reverse complement later!)
fwd_primer="TCCTCCGCTTATTGATATGC" #ITS4
rc_fwd_primer="GCATATCAATAAGCGGAGGA"
rev_primer="CAHCGATGAAGAACRYAG" #ITS3_kyo1
rc_rev_primer="CTRYGTTCTTCATCGDTG"

#Make empty file for collecting summary data during trimming (e.g., how many reads pass each step).
log="${outDir}/summary.tab"
echo -e "Sample\tRaw\tPrimer.trim\tReadThrough.trim" > ${log}

#Find the sequencer id from a fastq header for counting the number of sequences with grep 
seqID=`pigz -dc $(find ${muxIn} -name "*R1.fastq.gz" | head -n1) | awk 'NR==1' | awk -F ":" '{print $1}'`

#Loop over the forward reads (R1) in the demux folder.
#Exclude the "undetermined" samples (ie, those that did not get assigned to a sample durring demultiplexing).
for FWD in $( find ${muxIn} -name "*R1.fastq.gz" | grep -v "undetermined" ); do
    REV=`echo $FWD | sed 's/R1.fastq.gz/R2.fastq.gz/'` # Identify path to the matching reverse read (R2)
    
    #Identify current sample ID from the file name.
    SAMP=`basename $FWD | sed 's/.R1.fastq.gz//'`  

    ### Trim 1 ###
    # Remove gene primers and preceeding 3-6 Ns using cutadapt
    
    #First make output file names
    FWDtrim1="${scratch}/${SAMP}.R1.trim1.fq.gz"
    REVtrim1="${scratch}/${SAMP}.R2.trim1.fq.gz"
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #important cutadapt setting: overlap = length of shorter primer
    #I have increased the default error rate (-e) to allow 2 mismatches when seaching for primers (but this depends on primer lengths!) 
    cutadapt --quiet --discard-untrimmed --error-rate=0.15 --cores=0\
        -g "^N{3}${fwd_primer}" -g "^N{4}${fwd_primer}" -g "^N{5}${fwd_primer}" -g "^N{6}${fwd_primer}" \
        -G "^N{3}${rev_primer}" -G "^N{4}${rev_primer}" -G "^N{5}${rev_primer}" -G "^N{6}${rev_primer}" \
         -o $FWDtrim1 -p $REVtrim1 $FWD $REV

    #Count seqs after first trim step (or return 0 if cutadapt did not produce any output)
    if [[ -f $FWDtrim1 && -f $REVtrim1 ]]; then
        trim1Ct=`gzip -cd $FWDtrim1 | grep -c $seqID`
    else
        trim1Ct=0
    fi

    ### Trim 2 ###
    # Remove 3' read-through contamination (and filter short dimers) with SeqPurge.
    
    #This step will be skipped if cutadapt did not produce any output.
    if [ $trim1Ct -gt 0 ]; then
        #Set output file names
        FWDtrim2="${outDir}/${SAMP}.R1.fq.gz"
        REVtrim2="${outDir}/${SAMP}.R2.fq.gz"
        #Run SeqPurge (min_len should be reduced if shorter amplicons are expected)
        SeqPurge -in1 $FWDtrim1 -in2 $REVtrim1 -out1 $FWDtrim2 -out2 $REVtrim2 \
            -a1 $rc_rev_primer -a2 $rc_fwd_primer \
            -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}/${SAMP}.SeqPurge.log.txt
    fi

    #Count seqs remaining after second trim 
    if [[ -f $FWDtrim2 && -f $REVtrim2 ]]; then
        trim2Ct=`gzip -cd $FWDtrim2 | grep -c $seqID`
    else
        trim2Ct=0
    fi

    #Compile summary and append to file
    RAW=`gzip -cd $FWD | grep -c $seqID` #Number of sequenced in untrimmed input
    echo -e $SAMP"\t"$RAW"\t"$trim1Ct"\t"$trim2Ct | cat >> ${log}

    #remove tmp files
    rm ${scratch}/${SAMP}.*
done

#delete scratch folder
rm -r ${scratch}


