#!/bin/bash



##############################
#########  SETTINGS  #########
##############################
#uptag_index= has the path for the indexed signatures. This folder is provided in this repository. 
#You have to change this path for your computer path!!!
#Important! I generated the salmon index like this :
#salmon index -t Firmas_h.fasta -i barcodes_firmas_h-meliloti -k 21 --keepDuplicates 
#It is important to know that the signatures H4K4_2A and H4K4_5C are the same so you will have the same counts for both. 
#You could generate a new index avoiding these two signatures or you can use this Salmon index taking into account that 
#signatures H4K4_2A and H4K4_5C have the same DNA sequence
#you can run this to create the salmon index without duplicates 
#salmon index -t Firmas_h_no_duplicates.fasta -i barcodes_firmas_h-meliloti_no_dup -k 21 
#(see manual https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon)
#https://combine-lab.github.io/salmon/getting_started/

ncpus=7

# flags

while getopts "p:" flag;
do
    case "${flag}" in
        p) pdir=${OPTARG};;
    esac
done

samp=`basename ${pdir}` 
echo ">>> Processing batch ${samp} <<<"

uptag_index=/home/paula/Back_up/pau/S_meliloti/Sequencias_meliloti_exp_1/Planillas_firmas/barcodes_firmas_h-meliloti
#uptag_index=/home/paula/Back_up/pau/S_meliloti/Sequencias_meliloti_exp_1/Planillas_firmas/barcodes_firmas_h-meliloti_no_dup


##############################
#######  COUNT MATRIX ########
##############################

# Quants

echo -e "\n\n\n >>> SALMON QUANTIFICATION <<< \n\n\n"
 
for fn in ${pdir}/final/*.fastq.gz; 
do
name=`basename ${fn}`
echo -e "\n\n # Processing sample ${name} # \n\n"
salmon quant -i $uptag_index -l A\
    -r $fn \
    -p $ncpus \
    --noLengthCorrection \
    --fullLengthAlignment \
    --numAuxModelSamples 100000 \
    --output ${pdir}/quants/${name%_trimmed_merged_demux_adapter-trimmed.fastq.gz}
done


# Merge quants - count matrix, one for TPM and one for raw count numbers

echo -e "\n\n\n >>> EXPORTING QUANTIFICATION MATRICES <<< \n\n\n"

salmon quantmerge \
    --column TPM \
    --quants ${pdir}/quants/* \
    -o ${pdir}/${samp}_tpm.quant

salmon quantmerge \
    --column NUMREADS \
    --quants ${pdir}/quants/* \
    -o ${pdir}/${samp}_numreads.quant


echo -e "\n\n\n >>> FINISHED SALMON <<< \n\n\n"



