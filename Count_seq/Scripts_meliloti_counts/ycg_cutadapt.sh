#!/bin/bash

# usage:
#
# ./script.sh -p fastq-directory -i metadata file (equiv: ./script.sh -p [fastq directory]) 
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

#After trimming (previous step) you will have a new file. You need to update that name where it says name_seq
name_seq=NG-32116_Meliloti_first_assay_lib658176_10164_3


#################################
#######  DEMULTIPLEX AND  #######
#######  ADAPTER REMOVAL  #######
#################################

# merge read1 and read2
echo -e "\n\n >>> MERGE READS <<< \n\n"

bbmerge \
in1=${pdir}/trimmed/${name_seq}_1_val_1.fq.gz \
in2=${pdir}/trimmed/${name_seq}_2_val_2.fq.gz \
out=${pdir}/merged/${name_seq}_trimmed_merged.fastq.gz \
outu1=${pdir}/merged/${name_seq}_trimmed_unmerged_1.fastq.gz \
outu2=${pdir}/merged/${name_seq}_trimmed_unmerged_2.fastq.gz


# demultiplex
#the metadata.fasta is provided in this repository. It has the 5'TAGs 
#After running this script you will have new folders according to your 5'TAGs. After --rc all your sequences will have the same orientation  
#IMP! it is important to set the parameters -e 0 and the parameter O 26 because the variations among primers is small (four nucleotides). 
#Default settings are very relaxed. See manual cutadapt: https://cutadapt.readthedocs.io/en/stable/guide.html
#It is important to change the name of metadata for  the name of the folder plus the metadata (Test_complete_sequences_metadata.fasta)

echo -e "\n\n >>> DEMULTIPLEX BATCH <<< \n\n"
cutadapt \
    -g file:${pdir}/${samp}_metadata.fasta \
    -e 0 \
    -O 26 \
    --rename='{id} {comment} sample={adapter_name}' --rc \
    -j ${ncpus} \
    -o ${pdir}/demultiplexed_POOL/{name}_trimmed_merged_demux.fastq.gz ${pdir}/merged/${name_seq}_trimmed_merged.fastq.gz

echo -e "\n\n\n >>> FINISHED DEMULTIPLEX <<< \n\n\n"

    

