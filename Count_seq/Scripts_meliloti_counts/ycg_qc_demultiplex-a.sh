#!/bin/bash

##############################
#########  SETTINGS  #########
##############################

ncpus=7

# flags

while getopts "p:" flag;
do
    case "${flag}" in
        p) pdir=${OPTARG};;
    esac
done

samp=`basename ${pdir}` #basename removes the absolute path of the directory 
echo ">>> Processing batch ${samp} <<<"

#################################
#######  DEMULTIPLEX AND  #######
#######  ADAPTER REMOVAL  #######
#################################

#remove 3' adapter. After running this script you should have several folders organised by the TAGs 
#and the sequences you should see is just the signature of the mutants 

echo -e "\n\n >>> REMOVE COMMON SEQUENCES <<< \n\n"
for fn in ${pdir}/sequence_cutted/M*.fastq.gz;
do 
name=`basename ${fn}`
echo -e "\n\n >> Processing sample ${name} << \n\n"
cutadapt \
    -a AAGCTT \
    -e 0 \
    -j $ncpus \
    -o ${pdir}/final/${name%.fastq.gz}_to_quantify.fastq.gz ${fn}
done

echo -e "\n\n\n >>> FINISHED DEMULTIPLEX <<< \n\n\n"



