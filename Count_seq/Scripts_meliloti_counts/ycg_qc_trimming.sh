#!/bin/bash

# usage:
#
# ./script.sh -p fastq-directory -i metadata file (equiv: ./script.sh -p [fastq directory]) 
#
# # assumptions:
# uses ncpus cores - use command 'lscpu' to check available cores on your computer - change script accordingly


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

###############################
####### QC AND TRIMMING #######
###############################

# fastqc reports
# trimming
# fastqc reports after trimming

echo -e "\n\n\n >>> FastQC reports <<< \n\n\n"

fastqc -t 8 ${pdir}/*_1.fastq.gz ${pdir}/*_2.fastq.gz -o ${pdir}/fastqc &&
    echo -e "\n\n\n >>> TRIMMING of bad quality bases <<< \n\n\n" &&
    trim_galore -q 30 -j $ncpus --fastqc --trim-n \
    --paired ${pdir}/*_1.fastq.gz ${pdir}/*_2.fastq.gz \
    -fastqc --fastqc_args "--outdir ${pdir}/fastqc" \
    -o ${pdir}/trimmed

# -q means quality
#-j nucleos utilizados
#trim-n removes NNNN form the sides of the sequence
#--paired ask for a minimal length
#fastqc-args for quality 
#-o for outputfile

echo -e "\n >>> FINISHED TRIMMING <<< \n"

