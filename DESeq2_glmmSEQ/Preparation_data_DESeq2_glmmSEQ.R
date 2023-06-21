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
#output_file_name <- "/home/paula/Back_up/pau/S_meliloti_pool_4_9/Data_preparation_merge_counts_information/data_unificado_NO_NA.csv"
#write.table(data_unificado, file = output_file_name, 
 #          quote = FALSE, sep = "\t", row.names = FALSE)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL4-9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create a function to generate the desired expression
generate_expression <- function(pool_num) {
  # Create an empty list to store the data frames
  data_frames <- list()
  
  # Loop through the desired data frames
  for (i in 1:9) {
    data_frames[[i]] <- get(paste0("my_dataframe_DMSO_POOL", pool_num, "_", i))
    data_frames[[i+9]] <- get(paste0("my_dataframe_C1_POOL", pool_num, "_", i))
  }
  
  # Add Tabla_time_0_POOLX at the beginning and end of the list
  data_frames <- c(get(paste0("Tabla_time_0_POOL", pool_num)), data_frames, get(paste0("Tabla_time_0_POOL", pool_num)))
  
  # Combine all data frames horizontally using cbind
  result <- do.call(cbind, data_frames)
  
  # Return the result
  return(result)
}

# Generate the expressions for the different pools
POOL4_C1_DMSO_t0_t1_t2_t3 <- generate_expression(4)
POOL5_C1_DMSO_t0_t1_t2_t3 <- generate_expression(5)
POOL6_C1_DMSO_t0_t1_t2_t3 <- generate_expression(6)
POOL7_C1_DMSO_t0_t1_t2_t3 <- generate_expression(7)
POOL8_C1_DMSO_t0_t1_t2_t3 <- generate_expression(8)
POOL9_C1_DMSO_t0_t1_t2_t3 <- generate_expression(9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POOL4-9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set the names and eliminate the coloumns that are repeated (gene_ID)

pools <- c("POOL4_C1_DMSO_t0_t1_t2_t3", "POOL5_C1_DMSO_t0_t1_t2_t3","POOL6_C1_DMSO_t0_t1_t2_t3", "POOL7_C1_DMSO_t0_t1_t2_t3", "POOL8_C1_DMSO_t0_t1_t2_t3","POOL9_C1_DMSO_t0_t1_t2_t3" )
data_frames <- list(POOL4_C1_DMSO_t0_t1_t2_t3, POOL5_C1_DMSO_t0_t1_t2_t3,POOL6_C1_DMSO_t0_t1_t2_t3, POOL7_C1_DMSO_t0_t1_t2_t3, POOL8_C1_DMSO_t0_t1_t2_t3,POOL9_C1_DMSO_t0_t1_t2_t3)

# Loop over the pools
for (i in seq_along(pools)) {
  pool <- pools[i]
  data_frame <- data_frames[[i]]
  
  # Eliminate repeated columns (gene_ID)
  data_frame <- data_frame[, c(1, seq(2, ncol(data_frame), by = 2))]
  
  # Set column names
  col_names <- c("", paste0("Sample_", seq(0, ncol(data_frame)-2)))
  names(data_frame) <- col_names
  
  # Update the data frame in the list
  data_frames[[i]] <- data_frame
  
  # Assign the updated data frame back to the original object name
  assign(pool, data_frame)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Since I have a complete dataset I will merge all the datasets into one
Data_set_DSeq2 =  rbind(data_frames[[1]],data_frames[[2]],data_frames[[3]],data_frames[[4]],data_frames[[5]],data_frames[[6]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Save the cvs files it will be the countMatrix input for DESeq2 _Count_matrix_DESeq.csv
 
output_file_name <- "/home/paula/Back_up/pau/S_meliloti_pool_4_9/Data_preparation_merge_counts_information/Data_set_DSeq2_Count_matrix_DESeq.csv"
write.table(Data_set_DSeq2, file = output_file_name, 
           quote = FALSE, sep = "\t", row.names = FALSE)





