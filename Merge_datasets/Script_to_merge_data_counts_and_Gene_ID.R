library(dplyr)

filename1 = "Pobigaylo_y_Serratia_corrected.csv"
data1 = read.table(filename1, header=TRUE, sep=",")


filename2 = "Counts_data_POOL1_C1_DMSO.csv"
data2 = read.table(filename2, header=TRUE, sep=",")


data_unificado = data2 %>% 
	left_join(data1, by=c("transposon_set" , "tag_idfentifier"))


output_file_name <- "/home/paula/Back_up/pau/S_meliloti/ID_gene_meliloti_proteina/Merge_ID_gene_Signature_ESBL_table/Planilla_unificada.csv"
write.table(data_unificado, file = output_file_name, 
            quote = FALSE, sep = "\t", row.names = FALSE)

#=====================================excel=con=las=funciones==================================================
#Voy a usar esto en realidad cuando ya tenga los resulatdos de DESeq2
library("readxl")


df = read_xlsx("/home/paula/Back_up/pau/S_meliloti/ID_gene_meliloti_proteina/Merge_ID_gene_Signature_ESBL_table/All_genes_ID_meliloti_Ensembl.xlsx")

data3=as.data.frame(df)


data_unificado_funciones = data_unificado %>% 
	left_join(data3, by=c("gene_ID"))


	output_file_name <- "/home/paula/Back_up/pau/S_meliloti/ID_gene_meliloti_proteina/Merge_ID_gene_Signature_ESBL_table/Planilla_unificada_funciones.tsv"
write.table(data_unificado_funciones, file = output_file_name, 
            quote = FALSE, sep = "\t", row.names = FALSE)
