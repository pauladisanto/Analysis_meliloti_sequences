#!/bin/bash

# usage:
#
# ./script.sh -p fastq-directory -i metadata file (equiv: ./script.sh -p [fastq directory])
#
# 
# assumptions:
# uses ncpus cores - use command 'lscpu' to check available cores on your computer - change script accordingly
#

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

# remove 5' adapter. 

echo -e "\n\n >>> REMOVE COMMON SEQUENCES <<< \n\n"

for fn in ${pdir}/demultiplexed_CONDITION/M*.fastq.gz;
do 
name=`basename ${fn}`
echo -e "\n\n >> Processing sample ${name} << \n\n"
cutadapt \
    -g TACTAGCTCTACGACGGTCCACCTAAGCTT \
    -e 0 \
    -j $ncpus \
    -o ${pdir}/sequence_cutted/${name%.fastq.gz}_adapter-trimmed.fastq.gz ${fn}
done

echo -e "\n\n\n >>> FINISHED DEMULTIPLEX <<< \n\n\n"




