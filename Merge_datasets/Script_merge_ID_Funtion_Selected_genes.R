#Scritp to prepare the imput data for Gene analisys meliloti
library(dplyr)
library(readxl)
library(stringr)

filename1 = "Selected_genes_p_value_015.csv"
filename2 = "All_genes_ID_meliloti_Ensembl.xlsx"

data1 = read.table(filename1, header=TRUE, sep=",")
data2 = read_excel(filename2)


data1[c('Number', 'locus_tag')] <- str_split_fixed(data1$Gene_ID, '_', 2)
data1 <- subset(data1, select=c(locus_tag, Normalization,Pool))

################################unifico la informacion de la planilla 1 con la 2 a traves de tag_idfentifier######################

data_unificado = data1 %>% 
	left_join(data2, by=c("locus_tag"))

##############################################guardo la data_unificado en una tabla excel#############################

#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Meliloti_gene_barcodes/Planillas_filtradas_depuradas/data_unificado_firmas_ID.csv"
write.table(data_unificado, file = "Best_candidates_genes.xlsx", 
            quote = FALSE, sep = "\t", row.names = FALSE)