
#Two things I installed from https://www.ncbi.nlm.nih.gov/books/NBK179288/
#sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

#sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

#To activate EDirect for this terminal session, please execute the following
#export PATH=${PATH}:${HOME}/edirect
#One installation is complete, run:

 # export PATH=${PATH}:${HOME}/edirect

#to set the PATH for the current terminal session.

mkdir ref_genomes
cd ref_genomes


# get the ftp links for all reference genomes. It seems to have all the genomes...
wget ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt


# take out links for Pseudomas aeruginosas species (the ~ symbol has to be changed to == if one wants only "Pseudomonas aeruginosas"), 
# E. Coli and Sino.; make table with species (col 8) and link (col 20) - only use complete genomes or chromosome
awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Pseudomonas aeruginosa PAO1") print $8 "\t"  $20}' > ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Escherichia coli str. K-12 substr. MG1655") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Sinorhizobium meliloti 1021") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Pseudomonas putida KT2440") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Pectobacterium carotovorum subsp. carotovorum PC1") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8=="Staphylococcus epidermidis RP62A") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8~"Staphylococcus aureus strain=ST20190863") print $8 "\t"  $20}' >> ncbi_links.txt

awk -F '\t' '{if($12=="Chromosome" || $12=="Complete Genome") print}' assembly_summary_refseq.txt | awk -F '\t' '{if($8~"Paraburkholderia caledonica NBRC 102488") print $8 "\t"  $20}' >> ncbi_links.txt

# make directories for files  
mkdir RefSeqProteins
mkdir RefSeqReports
mkdir db

# get protein fasta files and assembly reports. I checked the links and I eliminated some of them

for next in $(cut -f2 ncbi_links.txt); do \
wget -r -nd -np -l1 -e robots=off -P RefSeqProteins/ -A "*protein.faa.gz" "$next"; \
wget -r -nd -np -l1 -e robots=off -P RefSeqReports/ -A "*assembly_report.txt" "$next"; \
wget -r -nd -np -l1 -e robots=off -P RefSeqGene/ -A "*genomic.fna.gz" "$next"
done

gunzip RefSeqProteins/*.gz

# save all protein sequences in one file and
# change header to keep protein ID and species - removes rest of text and brackets

#gawk -F' ' '{if ($1~">") {match($0, /(\[.*\])/, a); print $1 "_" a[0]} else {print}}' RefSeqProteins/*.faa | awk '{gsub(/[][]/, ""); gsub(" ", "_"); print;}' > db/all_proteins.fasta 

gawk -F' ' '{if ($1~">") {match($0, /(\[.*\])/, a); print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 "_" a[0]} else {print}}' RefSeqProteins/*.faa | awk '{gsub(/[][]/, ""); gsub(" ", "_"); print;}' > db/all_proteins.fasta 


#check that line counts match
wc -l RefSeqProteins/*.faa
wc -l db/all_proteins.fasta

############################
## Make a blast database 
##Aca hace el blast con la base de datos YSS_protein_db que la hace usando todas las sequencias fasta que tiene all_proteins
##segun entiendo esto incluiria las proteinas de Pseudomonas, E. coli y S. meliloti 1021. 
#Esta base de datos queda guardada en carpeta db
#Agregué del NCBI el proteoma de S meliloti porque el NCBI tiene otra nomenclatura
cd db
makeblastdb -in all_proteins.fasta -dbtype 'prot' -out YSS_protein_db

# Parameters
# in: FASTA file including all sequences for you database
# out: Database name
# dbtype: Database type usually 'prot' for proteins or 'nucleotide' for nucleotides

############################

### prepare QUERY SEQUENCES with one fasta file per sequence
## Used PlasmoDB to search for interesting genes - EC number, phenotype, expression etc. 
## Downloaded list of genes and used Uniprot to translate gene ID to Uniprot ID

#YO CREO QUE NO NECESITO ESTA PARTE Y QUE PUEDO HACERLO CON LAS SECUENCIAS FASTA QUE TENGO DESCARGADAS
# Download list of genes from PlasmoDB- custom and clear all choices (gives you Refseq id $1 and Plasmo ID $2)
# Remove $2 and first row (header)
#sed 's/\t.*//g' 69_proteases.txt | tail +2 > RefSeqID_69proteases.txt

## Translate RefSeq ID to Uniprot ID https://www.uniprot.org/uploadlists/

mkdir seq
cd seq

## retrieve sequences

#subseccion de secuencias fasta solo de las proteinas que estan en el ranking obtenido por el test estadistico glmmSEQ
#sin cortar la secuencia fasta
#Este archivo "prot_ID_meliloti.txt" lo tengo que preparar yo. Le tuve que remover ademas los "." de los nombres
seqtk subseq Sinorhizobium_meliloti_1021_gca_000006965.ASM696v1.pep.all.faa prot_ID_meliloti.txt > Selected_proteins_from_assay

#cortando la secuencia fasta por algun motivo no me selecciona todas las que quiero
#grep -w -A 2 -f  prot_ID_meliloti.txt Sinorhizobium_meliloti_1021_gca_000006965.ASM696v1.pep.all.faa --no-group-separator > selected_proteins_from_assay.faa

############################

### BLAST

blastp -query Selected_proteins_from_assay.faa -db YSS_protein_db -out Results_aligments_meliloti_reference_strains.txt -evalue 10e-15 -word_size 3 -num_threads 7

blastp -query Selected_proteins_from_assay.faa -db YSS_protein_db -out Results_aligments_meliloti_reference_strains.tsv -outfmt 6 -evalue 10e-15 -word_size 3 -num_threads 7

# query: Input sequence
# db: BLAST database
# out: Name of the output file
# outfmt: Output format
# evalue: Arbitrary cut-off of sequence similarity. The e-value depends on your database size, so the large your database the smaller your e-value can be. I recommend something between 10-5 and 10-20.
# word_size: Number of nucleotides/amino acids, which resembles the smallest unit of your query
# num_alignments: Maximum number of alignment partner in the database for a single query sequence



