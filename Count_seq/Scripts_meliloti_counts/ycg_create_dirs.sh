#!/bin/bash

# usage:
#
# ./script.sh -p fastq-directory -i metadata file (equiv: ./script.sh -p [fastq directory])
#

##############################
#########  SETTINGS  #########
##############################

# flags

while getopts "p:" flag;
do
    case "${flag}" in
        p) pdir=${OPTARG};;
    esac
done

samp=`basename ${pdir}` #basename removes the absolute path of the directory 
echo ">>> Processing batch ${samp} <<<"

# make directories 
echo -e "\n\n\n >>> Make some directories if they don't already exist <<< \n\n\n"
mkdir ${pdir}/fastqc
mkdir ${pdir}/trimmed
mkdir ${pdir}/merged
mkdir ${pdir}/demultiplexed_POOL
mkdir ${pdir}/demultiplexed_CONDITION
mkdir ${pdir}/sequence_cutted
mkdir ${pdir}/final
mkdir ${pdir}/quants

echo -e "\n >>> FINISHED CREATING DIRS <<< \n"

