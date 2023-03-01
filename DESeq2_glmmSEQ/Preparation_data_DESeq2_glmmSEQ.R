#prepare the sf files to make the csv file for DESeq analysis

library(tximportData)
library(dplyr)
library(tidyverse)
library(purrr)
#merged the tables, after renaming the sf files from Salmon using a mini script in Python (also provided in this folder)

merged_df = list.files(pattern= "*") %>%
   set_names() %>%
   map_dfr(read.delim, .id="file_name")


setwd("/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/Planillas_intermedias")

#This file name only has the samples that I put in the misture to be sequenced
filename2 = "Experimento_1_DMSO_C1_NUEVA.csv"
data2 = read.table(filename2, header=TRUE, sep=",")

data_unificado = merged_df %>% 
   left_join(data2, by=c("file_name"))

#to remove data that does not match with our TAGs but have a higer number of reads (higer than 2 and they were not filtered in the Python script)
data_unificado <- data_unificado[!is.na(data_unificado$Replica),]

#to save the dataframe and take a detailed check
#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/Planillas_intermedias/Data_counts_Experimento1_C1_DMSO_2.csv"
#write.table(data_unificado, file = output_file_name, 
#           quote = FALSE, sep = "\t", row.names = FALSE)

#Second merge with the dataframe that contains the ID name. Merge by transposon_set(number of pool) and Name(signature)

filename3 = "Pobigaylo_y_Serratia_corrected.csv"
data3 = read.table(filename3, header=TRUE, sep=",")

#Rename of the column tag_idfentifier for name
data3 <- data3 %>%
       rename("Name" = "tag_idfentifier")

data_unificado2 = data_unificado %>% 
   left_join(data3, by=c("transposon_set","Name"))

#to remove data with no ID gene, since there is no chance of assigning a funtion
data_unificado2 <- data_unificado2[!is.na(data_unificado2$gene_ID),]

data_unificado2 = data_unificado2 %>% 
 mutate_if(is.numeric, round, digits=0)

#to save the dataframe and take a detailed check
#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/Data_counts_plus_gene_ID_round_values.csv"
#write.table(data_unificado2, file = output_file_name, 
 #          quote = FALSE, sep = "\t", row.names = FALSE)

#select the important columns for DESeq2: ID gene and counts

Tabla_Compound_DMSO_POOL1 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_DMSO_POOL1[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 1, select=c(gene_ID, NumReads))

}
}
#--------------------------------------------------------------------------------------------------------------------------------
Tabla_Compound_DMSO_POOL2 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_DMSO_POOL2[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 2, select=c(gene_ID, NumReads))


}
}
#--------------------------------------------------------------------------------------------------------------------------------


Tabla_Compound_DMSO_POOL3 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_DMSO_POOL3[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 3, select=c(gene_ID, NumReads))

}
}

#--------------------------------------------------------------------------------------------------------------------------------

Tabla_Compound_C1_POOL1 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_C1_POOL1[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 1, select=c(gene_ID, NumReads))

}
}
#--------------------------------------------------------------------------------------------------------------------------------

Tabla_Compound_C1_POOL2 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_C1_POOL2[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 2, select=c(gene_ID, NumReads))

}
}

#--------------------------------------------------------------------------------------------------------------------------------

Tabla_Compound_C1_POOL3 = list()
index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_C1_POOL3[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 3, select=c(gene_ID, NumReads))

}
}

#--------------------------------------------------------------------------------------------------------------------------------

#prepare the data sets for DESeq2
#====POOL1====DMSO=========================
DMSO_POOL1_R1_t1 = Tabla_Compound_DMSO_POOL1[[1]]
#DMSO_POOL1_R1_t2 = Tabla_Compound_DMSO_POOL1[[2]]
DMSO_POOL1_R1_t3 = Tabla_Compound_DMSO_POOL1[[3]]
DMSO_POOL1_R2_t1 = Tabla_Compound_DMSO_POOL1[[4]]
DMSO_POOL1_R2_t2 = Tabla_Compound_DMSO_POOL1[[5]]
DMSO_POOL1_R2_t3 = Tabla_Compound_DMSO_POOL1[[6]]
DMSO_POOL1_R3_t1 = Tabla_Compound_DMSO_POOL1[[7]]
DMSO_POOL1_R3_t2 = Tabla_Compound_DMSO_POOL1[[8]]
DMSO_POOL1_R3_t3 = Tabla_Compound_DMSO_POOL1[[9]]

#====POOL2====DMSO=========================
DMSO_POOL2_R1_t1 = Tabla_Compound_DMSO_POOL2[[1]]
#DMSO_POOL2_R1_t2 = Tabla_Compound_DMSO_POOL2[[2]]
DMSO_POOL2_R1_t3 = Tabla_Compound_DMSO_POOL2[[3]]
#DMSO_POOL2_R2_t1 = Tabla_Compound_DMSO_POOL2[[4]]
DMSO_POOL2_R2_t2 = Tabla_Compound_DMSO_POOL2[[5]]
DMSO_POOL2_R2_t3 = Tabla_Compound_DMSO_POOL2[[6]]
#DMSO_POOL2_R3_t1 = Tabla_Compound_DMSO_POOL2[[7]]
#DMSO_POOL2_R3_t2 = Tabla_Compound_DMSO_POOL2[[8]]
DMSO_POOL2_R3_t3 = Tabla_Compound_DMSO_POOL2[[9]]

#====POOL3====DMSO=========================
DMSO_POOL3_R1_t1 = Tabla_Compound_DMSO_POOL3[[1]]
#DMSO_POOL3_R1_t2 = Tabla_Compound_DMSO_POOL3[[2]]
DMSO_POOL3_R1_t3 = Tabla_Compound_DMSO_POOL3[[3]]
#DMSO_POOL3_R2_t1 = Tabla_Compound_DMSO_POOL3[[4]]
DMSO_POOL3_R2_t2 = Tabla_Compound_DMSO_POOL3[[5]]
DMSO_POOL3_R2_t3 = Tabla_Compound_DMSO_POOL3[[6]]
DMSO_POOL3_R3_t1 = Tabla_Compound_DMSO_POOL3[[7]]
DMSO_POOL3_R3_t2 = Tabla_Compound_DMSO_POOL3[[8]]
DMSO_POOL3_R3_t3 = Tabla_Compound_DMSO_POOL3[[9]]

#====POOL1====C1=========================
C1_POOL1_R1_t1 = Tabla_Compound_C1_POOL1[[1]]
#C1_POOL1_R1_t2 = Tabla_Compound_C1_POOL1[[2]]
C1_POOL1_R1_t3 = Tabla_Compound_C1_POOL1[[3]]
C1_POOL1_R2_t1 = Tabla_Compound_C1_POOL1[[4]]
C1_POOL1_R2_t2 = Tabla_Compound_C1_POOL1[[5]]
C1_POOL1_R2_t3 = Tabla_Compound_C1_POOL1[[6]]
C1_POOL1_R3_t1 = Tabla_Compound_C1_POOL1[[7]]
C1_POOL1_R3_t2 = Tabla_Compound_C1_POOL1[[8]]
C1_POOL1_R3_t3 = Tabla_Compound_C1_POOL1[[9]]

#====POOL2====C1=========================
C1_POOL2_R1_t1 = Tabla_Compound_C1_POOL2[[1]]
#C1_POOL2_R1_t2 = Tabla_Compound_C1_POOL2[[2]]
C1_POOL2_R1_t3 = Tabla_Compound_C1_POOL2[[3]]
C1_POOL2_R2_t1 = Tabla_Compound_C1_POOL2[[4]]
#C1_POOL2_R2_t2 = Tabla_Compound_C1_POOL2[[5]]
C1_POOL2_R2_t3 = Tabla_Compound_C1_POOL2[[6]]
C1_POOL2_R3_t1 = Tabla_Compound_C1_POOL2[[7]]
C1_POOL2_R3_t2 = Tabla_Compound_C1_POOL2[[8]]
C1_POOL2_R3_t3 = Tabla_Compound_C1_POOL2[[9]]

#====POOL3====C1=========================
C1_POOL3_R1_t1 = Tabla_Compound_C1_POOL3[[1]]
#C1_POOL3_R1_t2 = Tabla_Compound_C1_POOL3[[2]]
C1_POOL3_R1_t3 = Tabla_Compound_C1_POOL3[[3]]
C1_POOL3_R2_t1 = Tabla_Compound_C1_POOL3[[4]]
#C1_POOL3_R2_t2 = Tabla_Compound_C1_POOL3[[5]]
C1_POOL3_R2_t3 = Tabla_Compound_C1_POOL3[[6]]
C1_POOL3_R3_t1 = Tabla_Compound_C1_POOL3[[7]]
C1_POOL3_R3_t2 = Tabla_Compound_C1_POOL3[[8]]
C1_POOL3_R3_t3 = Tabla_Compound_C1_POOL3[[9]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
POOL_1_C1_DMSO_t1_t2_t3 = cbind(DMSO_POOL1_R1_t1,DMSO_POOL1_R1_t3,DMSO_POOL1_R2_t1,DMSO_POOL1_R2_t2,
   DMSO_POOL1_R2_t3,DMSO_POOL1_R3_t1,DMSO_POOL1_R3_t2,DMSO_POOL1_R3_t3,
   C1_POOL1_R1_t1, C1_POOL1_R1_t3,C1_POOL1_R2_t1,C1_POOL1_R2_t2, C1_POOL1_R2_t3,C1_POOL1_R3_t1,
   C1_POOL1_R3_t2,C1_POOL1_R3_t3)


POOL_1_C1_DMSO_t1_t2_t3 = POOL_1_C1_DMSO_t1_t2_t3[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)]
POOL_1_C1_DMSO_t1_t2_t3<-  setNames(POOL_1_C1_DMSO_t1_t2_t3, c("","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_2_C1_DMSO_t1_t2_t3 = cbind(DMSO_POOL2_R1_t1, DMSO_POOL2_R1_t3, DMSO_POOL2_R2_t2, DMSO_POOL2_R2_t3, DMSO_POOL2_R3_t3, 
   C1_POOL2_R1_t1, C1_POOL2_R1_t3, C1_POOL2_R2_t1, C1_POOL2_R2_t3, C1_POOL2_R3_t1, C1_POOL2_R3_t2, C1_POOL2_R3_t3)

POOL_2_C1_DMSO_t1_t2_t3 = POOL_2_C1_DMSO_t1_t2_t3[,c(1,2,4,6,8,10,12,14,16,18,20,22,24)]
POOL_2_C1_DMSO_t1_t2_t3<-  setNames(POOL_2_C1_DMSO_t1_t2_t3, c("","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_3_C1_DMSO_t1_t2_t3 = cbind(DMSO_POOL3_R1_t1, DMSO_POOL3_R1_t3, DMSO_POOL3_R2_t2, DMSO_POOL3_R2_t3, DMSO_POOL3_R3_t1, 
   DMSO_POOL3_R3_t2, DMSO_POOL3_R3_t3, C1_POOL3_R1_t1, C1_POOL3_R1_t3, C1_POOL3_R2_t1, C1_POOL3_R2_t3, C1_POOL3_R3_t1, 
   C1_POOL3_R3_t2, C1_POOL3_R3_t3)

POOL_3_C1_DMSO_t1_t2_t3 = POOL_3_C1_DMSO_t1_t2_t3[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
POOL_3_C1_DMSO_t1_t2_t3<-  setNames(POOL_3_C1_DMSO_t1_t2_t3, c("","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Save the cvs files it will be the countMatrix input for DESeq2 _Count_matrix_DESeq.csv
 
output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/0_POOL_1_C1_DMSO_t1_t2_t3_Count_matrix_DESeq.csv"
write.table(POOL_1_C1_DMSO_t1_t2_t3, file = output_file_name, 
           quote = FALSE, sep = "\t", row.names = FALSE)

output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/0_POOL_2_C1_DMSO_t1_t2_t3_Count_matrix_DESeq.csv"
write.table(POOL_2_C1_DMSO_t1_t2_t3, file = output_file_name, 
           quote = FALSE, sep = "\t", row.names = FALSE)

output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/0_POOL_3_C1_DMSO_t1_t2_t3_Count_matrix_DESeq.csv"
write.table(POOL_3_C1_DMSO_t1_t2_t3, file = output_file_name, 
           quote = FALSE, sep = "\t", row.names = FALSE)



