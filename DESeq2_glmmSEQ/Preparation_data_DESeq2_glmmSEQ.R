#prepare the sf files to make the csv file for DESeq analysis
#source("Preparation_data_DESeq2_glmmSEQ.R")
library(tximportData)
library(dplyr)
library(tidyverse)
library(purrr)
#merged the tables, after renaming the sf files from Salmon using a mini script in Python (also provided in this folder)

merged_df = list.files(pattern= "*") %>%
   set_names() %>%
   map_dfr(read.delim, .id="file_name")


#setwd("/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/Planillas_intermedias")
setwd("/home/paula/Back_up/pau/S_meliloti_pool_4_9/Data_preparation_merge_counts_information")
#This file name only has the samples that I put in the misture to be sequenced
filename2 = "Experimento_2_Y_3_DMSO_C1_NUEVA.csv"
data2 = read.table(filename2, header=TRUE, sep=",")
#generate a new coloumn with the informations of the coloumns of the primers

data2= data2 %>%
     unite('file_name',Int_label_1:Int_label_2, remove = FALSE)

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

Tabla_Compound_DMSO_POOL4 = list()
Tabla_Compound_DMSO_POOL5 = list()
Tabla_Compound_DMSO_POOL6 = list()
Tabla_Compound_DMSO_POOL7 = list()
Tabla_Compound_DMSO_POOL8 = list()
Tabla_Compound_DMSO_POOL9 = list()


index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_DMSO_POOL4[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 4, select=c(gene_ID, NumReads))
Tabla_Compound_DMSO_POOL5[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 5, select=c(gene_ID, NumReads))
Tabla_Compound_DMSO_POOL6[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 6, select=c(gene_ID, NumReads))
Tabla_Compound_DMSO_POOL7[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 7, select=c(gene_ID, NumReads))
Tabla_Compound_DMSO_POOL8[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 8, select=c(gene_ID, NumReads))
Tabla_Compound_DMSO_POOL9[[index]]<-subset(data_unificado2,Compuesto=="DMSO" & Tiempo== B[j] & Replica== A[i] & transposon_set == 9, select=c(gene_ID, NumReads))

}
}
#--------------------------------------------------------------------------------------------------------------------------------

Tabla_Compound_C1_POOL4 = list()
Tabla_Compound_C1_POOL5 = list()
Tabla_Compound_C1_POOL6 = list()
Tabla_Compound_C1_POOL7 = list()
Tabla_Compound_C1_POOL8 = list()
Tabla_Compound_C1_POOL9 = list()

index=0
B = c(1,2,3)
A = c("R1", "R2", "R3")
for(i in 1:3) {
   for(j in 1:3) {
     index=index + 1

Tabla_Compound_C1_POOL4[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 4, select=c(gene_ID, NumReads))
Tabla_Compound_C1_POOL5[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 5, select=c(gene_ID, NumReads))
Tabla_Compound_C1_POOL6[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 6, select=c(gene_ID, NumReads))
Tabla_Compound_C1_POOL7[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 7, select=c(gene_ID, NumReads))
Tabla_Compound_C1_POOL8[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 8, select=c(gene_ID, NumReads))
Tabla_Compound_C1_POOL9[[index]]<-subset(data_unificado2,Compuesto=="C1" & Tiempo== B[j] & Replica== A[i] & transposon_set == 9, select=c(gene_ID, NumReads))

}
}


#--------------------------------------------------------------------------------------------------------------------------------
#select the important columns for DESeq2: ID gene and counts for time 0
Tabla_time_0_POOL4<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 4, select=c(gene_ID, NumReads))
Tabla_time_0_POOL5<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 5, select=c(gene_ID, NumReads))
Tabla_time_0_POOL6<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 6, select=c(gene_ID, NumReads))
Tabla_time_0_POOL7<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 7, select=c(gene_ID, NumReads))
Tabla_time_0_POOL8<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 8, select=c(gene_ID, NumReads))
Tabla_time_0_POOL9<-subset(data_unificado2,Compuesto=="NONE" & Tiempo== 0 & Replica== "R1" & transposon_set == 9, select=c(gene_ID, NumReads))

#--------------------------------------------------------------------------------------------------------------------------------
# Create objects using a loop for collect dataframes DMSO
for (i in 1:length(Tabla_Compound_DMSO_POOL4)) {
  assign(paste0("my_dataframe_DMSO_POOL4_", i), Tabla_Compound_DMSO_POOL4[[i]])
  assign(paste0("my_dataframe_DMSO_POOL5_", i), Tabla_Compound_DMSO_POOL5[[i]])
  assign(paste0("my_dataframe_DMSO_POOL6_", i), Tabla_Compound_DMSO_POOL6[[i]])
  assign(paste0("my_dataframe_DMSO_POOL7_", i), Tabla_Compound_DMSO_POOL7[[i]])
  assign(paste0("my_dataframe_DMSO_POOL8_", i), Tabla_Compound_DMSO_POOL8[[i]])
  assign(paste0("my_dataframe_DMSO_POOL9_", i), Tabla_Compound_DMSO_POOL9[[i]])

}

#--------------------------------------------------------------------------------------------------------------------------------
# Create objects using a loop for collect dataframes C1

# Create objects using a loop for collect dataframes DMSO
for (i in 1:length(Tabla_Compound_C1_POOL4)) {
  assign(paste0("my_dataframe_C1_POOL4_", i), Tabla_Compound_C1_POOL4[[i]])
  assign(paste0("my_dataframe_C1_POOL5_", i), Tabla_Compound_C1_POOL5[[i]])
  assign(paste0("my_dataframe_C1_POOL6_", i), Tabla_Compound_C1_POOL6[[i]])
  assign(paste0("my_dataframe_C1_POOL7_", i), Tabla_Compound_C1_POOL7[[i]])
  assign(paste0("my_dataframe_C1_POOL8_", i), Tabla_Compound_C1_POOL8[[i]])
  assign(paste0("my_dataframe_C1_POOL9_", i), Tabla_Compound_C1_POOL9[[i]])

}

#--------------------------------------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


POOL_4_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL4, my_dataframe_DMSO_POOL4_1,my_dataframe_DMSO_POOL4_2,my_dataframe_DMSO_POOL4_3,my_dataframe_DMSO_POOL4_4,
   my_dataframe_DMSO_POOL4_5, my_dataframe_DMSO_POOL4_6, my_dataframe_DMSO_POOL4_7, my_dataframe_DMSO_POOL4_8, my_dataframe_DMSO_POOL4_9,
   my_dataframe_C1_POOL4_1, my_dataframe_C1_POOL4_2, my_dataframe_C1_POOL4_3, my_dataframe_C1_POOL4_4, my_dataframe_C1_POOL4_5, my_dataframe_C1_POOL4_6,
   my_dataframe_C1_POOL4_7, my_dataframe_C1_POOL4_8, my_dataframe_C1_POOL4_9, Tabla_time_0_POOL4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_5_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL5, my_dataframe_DMSO_POOL5_1,my_dataframe_DMSO_POOL5_2,my_dataframe_DMSO_POOL5_3,my_dataframe_DMSO_POOL5_4,
   my_dataframe_DMSO_POOL5_5, my_dataframe_DMSO_POOL5_6, my_dataframe_DMSO_POOL5_7, my_dataframe_DMSO_POOL5_8, my_dataframe_DMSO_POOL5_9,
   my_dataframe_C1_POOL5_1, my_dataframe_C1_POOL5_2, my_dataframe_C1_POOL5_3, my_dataframe_C1_POOL5_4, my_dataframe_C1_POOL5_5, my_dataframe_C1_POOL5_6,
   my_dataframe_C1_POOL5_7, my_dataframe_C1_POOL5_8, my_dataframe_C1_POOL5_9, Tabla_time_0_POOL5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL6~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_6_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL6, my_dataframe_DMSO_POOL6_1,my_dataframe_DMSO_POOL6_2,my_dataframe_DMSO_POOL6_3,my_dataframe_DMSO_POOL6_4,
   my_dataframe_DMSO_POOL6_5, my_dataframe_DMSO_POOL6_6, my_dataframe_DMSO_POOL6_7, my_dataframe_DMSO_POOL6_8, my_dataframe_DMSO_POOL6_9,
   my_dataframe_C1_POOL6_1, my_dataframe_C1_POOL6_2, my_dataframe_C1_POOL6_3, my_dataframe_C1_POOL6_4, my_dataframe_C1_POOL6_5, my_dataframe_C1_POOL6_6,
   my_dataframe_C1_POOL6_7, my_dataframe_C1_POOL6_8, my_dataframe_C1_POOL6_9, Tabla_time_0_POOL6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


POOL_7_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL7, my_dataframe_DMSO_POOL7_1,my_dataframe_DMSO_POOL7_2,my_dataframe_DMSO_POOL7_3,my_dataframe_DMSO_POOL7_4,
   my_dataframe_DMSO_POOL7_5, my_dataframe_DMSO_POOL7_6, my_dataframe_DMSO_POOL7_7, my_dataframe_DMSO_POOL7_8, my_dataframe_DMSO_POOL7_9, 
   my_dataframe_C1_POOL7_1, my_dataframe_C1_POOL7_2, my_dataframe_C1_POOL7_3, my_dataframe_C1_POOL7_4, my_dataframe_C1_POOL7_5, my_dataframe_C1_POOL7_6,
   my_dataframe_C1_POOL7_7, my_dataframe_C1_POOL7_8, my_dataframe_C1_POOL7_9,Tabla_time_0_POOL7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_8_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL8, my_dataframe_DMSO_POOL8_1,my_dataframe_DMSO_POOL8_2,my_dataframe_DMSO_POOL8_3,my_dataframe_DMSO_POOL8_4,
   my_dataframe_DMSO_POOL8_5, my_dataframe_DMSO_POOL8_6, my_dataframe_DMSO_POOL8_7, my_dataframe_DMSO_POOL8_8, my_dataframe_DMSO_POOL8_9,
   my_dataframe_C1_POOL8_1, my_dataframe_C1_POOL8_2, my_dataframe_C1_POOL8_3, my_dataframe_C1_POOL8_4, my_dataframe_C1_POOL8_5, my_dataframe_C1_POOL8_6,
   my_dataframe_C1_POOL8_7, my_dataframe_C1_POOL8_8, my_dataframe_C1_POOL8_9, Tabla_time_0_POOL8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_9_C1_DMSO_t0_t1_t2_t3 =cbind(Tabla_time_0_POOL9, my_dataframe_DMSO_POOL9_1,my_dataframe_DMSO_POOL9_2,my_dataframe_DMSO_POOL9_3,my_dataframe_DMSO_POOL9_4,
   my_dataframe_DMSO_POOL9_5, my_dataframe_DMSO_POOL9_6, my_dataframe_DMSO_POOL9_7, my_dataframe_DMSO_POOL9_8, my_dataframe_DMSO_POOL9_9, 
   my_dataframe_C1_POOL9_1, my_dataframe_C1_POOL9_2, my_dataframe_C1_POOL9_3, my_dataframe_C1_POOL9_4, my_dataframe_C1_POOL9_5, my_dataframe_C1_POOL9_6,
   my_dataframe_C1_POOL9_7, my_dataframe_C1_POOL9_8, my_dataframe_C1_POOL9_9,Tabla_time_0_POOL9)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set the names and eliminate the coloumns that are repeated (gene_ID)

POOL_4_C1_DMSO_t0_t1_t2_t3 = POOL_4_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_4_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_4_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09",  "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
POOL_5_C1_DMSO_t0_t1_t2_t3 = POOL_5_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_5_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_5_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL6~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_6_C1_DMSO_t0_t1_t2_t3 = POOL_6_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_6_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_6_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

POOL_7_C1_DMSO_t0_t1_t2_t3 = POOL_7_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_7_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_7_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
POOL_8_C1_DMSO_t0_t1_t2_t3 = POOL_8_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_8_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_8_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09",  "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
POOL_9_C1_DMSO_t0_t1_t2_t3 = POOL_9_C1_DMSO_t0_t1_t2_t3 [,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

POOL_9_C1_DMSO_t0_t1_t2_t3 <-  setNames(POOL_9_C1_DMSO_t0_t1_t2_t3, c("","Sample_0","Sample_01", "Sample_02", "Sample_03", "Sample_04", "Sample_05", 
   "Sample_06", "Sample_07", "Sample_08", "Sample_09", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_14", "Sample_15", "Sample_16", "Sample_17","Sample_18","Sample_19"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Since I have a complete dataset I will merge all the datasets into one

Data_set_DSeq2 =  rbind(POOL_4_C1_DMSO_t0_t1_t2_t3, POOL_5_C1_DMSO_t0_t1_t2_t3,POOL_6_C1_DMSO_t0_t1_t2_t3,POOL_7_C1_DMSO_t0_t1_t2_t3,POOL_8_C1_DMSO_t0_t1_t2_t3, POOL_9_C1_DMSO_t0_t1_t2_t3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Save the cvs files it will be the countMatrix input for DESeq2 _Count_matrix_DESeq.csv
 
output_file_name <- "/home/paula/Back_up/pau/S_meliloti_pool_4_9/Data_preparation_merge_counts_information/Data_set_DSeq2_Count_matrix_DESeq.csv"
write.table(Data_set_DSeq2, file = output_file_name, 
           quote = FALSE, sep = "\t", row.names = FALSE)


#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/0_POOL_2_C1_DMSO_t1_t2_t3_Count_matrix_DESeq.csv"
#write.table(POOL_2_C1_DMSO_t1_t2_t3, file = output_file_name, 
 #          quote = FALSE, sep = "\t", row.names = FALSE)

#output_file_name <- "/home/paula/Back_up/pau/S_meliloti/Normalizacion_Meliloti_primer_exp/otras_planillas/0_POOL_3_C1_DMSO_t1_t2_t3_Count_matrix_DESeq.csv"
#write.table(POOL_3_C1_DMSO_t1_t2_t3, file = output_file_name, 
 #          quote = FALSE, sep = "\t", row.names = FALSE)



