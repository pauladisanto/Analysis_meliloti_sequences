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

#the metadata2.fasta is provided in this repository. It has the 3'TAGs. Reverse complement. 


echo -e "\n\n >>> DEMULTIPLEX BATCH <<< \n\n"

files=${pdir}/demultiplexed_POOL/*Forward*fastq.gz

for f in $files
do
  name_seq=`basename ${f}`
  echo "Processing $name_seq file..."
  cutadapt \
    -a file:${pdir}/${samp}_metadata2.fasta \
    -e 0 \
    -O 26 \
    --rename='{id} {comment} sample={adapter_name}' --rc \
    -j ${ncpus} \
    -o ${pdir}/demultiplexed_CONDITION/{name}_${name_seq}_demux_a.fastq.gz $f
done

echo -e "\n\n\n >>> FINISHED DEMULTIPLEX <<< \n\n\n"


