
library(dplyr)
library(readxl)
library(stringr)
library(tidyverse)


filename = c("stats_ordered_POOL1_Norm_DESeq2.csv", "stats_ordered_POOL2_Norm_DESeq2.csv", "stats_ordered_POOL3_Norm_DESeq2.csv")

data = list()

for(i in 1: length(filename)) {

data[[i]] = read.table(filename[i], header=TRUE, sep=",")

}

#############################################################################################


data1= rbind(data[[1]], data[[2]], data[[3]])

Tabla<-subset(data1, P_Time.Condition < 0.1, select=c(X, Time, Time:ConditionUntreated, P_Time,	P_Condition, P_Time.Condition))

Increased_mutants = subset(Tabla, Time > 0, select=c(X, Time, Time:ConditionUntreated, P_Time,	P_Condition, P_Time.Condition))

Decreased_mutants = subset(Tabla, Time < 0, select=c(X, Time, Time:ConditionUntreated, P_Time,	P_Condition, P_Time.Condition))



#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Seleccion_mutantes_candidatos/Increased_mutants_POOLs_1_2_3.csv"
#write.table(Increased_mutants, file = output_file_name, 
#           quote = FALSE, sep = "\t", row.names = FALSE)

#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Seleccion_mutantes_candidatos/Decreased_mutants_POOLs_1_2_3.csv"
#write.table(Decreased_mutants, file = output_file_name, 
#           quote = FALSE, sep = "\t", row.names = FALSE)

#####the next script is to transform the gene ID into Protein ID and retrieve the portein sequence#####



Increased_mutants[c('Number', 'gene_ID')] <- str_split_fixed(Increased_mutants$X, '_', 2)

Decreased_mutants[c('Number', 'gene_ID')] <- str_split_fixed(Decreased_mutants$X, '_', 2)

filename = "data_unificado_funtion_ID_gene_protein_plus_proteins_from_fasta_2.csv"

data_Gene_ID_Prot_ID= read.table(filename, header=TRUE, sep="\t")



data_unificado_increased = Increased_mutants %>% 
    left_join(data_Gene_ID_Prot_ID, by=c("gene_ID"))


data_unificado_decreased = Decreased_mutants %>% 
    left_join(data_Gene_ID_Prot_ID, by=c("gene_ID"))


Protein_names_increased = subset(data_unificado_increased, select=c(protein_id)) %>% 
separate(protein_id, c("protein_id", "number"))
Protein_names_increased =  subset(Protein_names_increased, select=c(protein_id))
Protein_names_increased = Protein_names_increased[!is.na(Protein_names_increased$protein_id),]


Protein_names_decreased = subset(data_unificado_decreased, select=c(protein_id)) %>% 
separate(protein_id, c("protein_id", "number"))
Protein_names_decreased =  subset(Protein_names_decreased, select=c(protein_id))
Protein_names_decreased = Protein_names_decreased[!is.na(Protein_names_decreased$protein_id),]


#####Retrieve the sequences using these lines in BASH#####
write.csv(Protein_names_increased,file="Protein_names_increased_POOL_1_2_3.csv",row.names=F)
write.csv(Protein_names_decreased,file="Protein_names_decreased_POOL_1_2_3.csv",row.names=F)


sort Protein_names_increased_POOL_1_2_3.csv | uniq | tr -d '"'>  Protein_names_increased_POOL_1_2_3.csv_unique.txt

sort Protein_names_decreased_POOL_1_2_3.csv | uniq | tr -d '"'>  Protein_names_decreased_POOL_1_2_3.csv_unique.txt

mkdir protein_sequences
cd protein_sequences


seqtk subseq Sinorhizobium_meliloti_1021_gca_000006965.ASM696v1.pep.all.faa Protein_names_increased_POOL_1_2_3.csv_unique.txt > Protein_sequence_Increased_mutants


seqtk subseq Sinorhizobium_meliloti_1021_gca_000006965.ASM696v1.pep.all.faa Protein_names_decreased_POOL_1_2_3.csv_unique.txt > Protein_sequence_Decreased_mutants



########################Kegg databse#########################################

#luego de mandar la lista de proteinas kegg me devuelve un txt que tiene el nombre de la proteina y el termino de kegg que necesito
#quiero seleccionar solo la segunda columna


#selecciona columna 2
awk '{print $2}' user_ko_decreased_mutants_POOL_1_2_3.txt > list_KEGG_decreased

#quita los espacios

sed -i '/^$/d' list_KEGG_decreased 
###########################################

awk '{print $2}' user_ko_increased_mutants_POOL_1_2_3.txt > list_KEGG_increased

#quita los espacios

sed -i '/^$/d' list_KEGG_increased 

