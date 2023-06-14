#https://cran.r-project.org/web/packages/glmmSeq/vignettes/glmmSeq.html
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(glmmSeq)
library(DESeq2)
library(dplyr)
library(ggplot2)

filename1 = "Col_data_POOL_4_9.csv"
filename2 = "Data_set_DSeq2_Count_matrix_DESeq.csv"

coldata = read.table(filename1, header=TRUE, sep=",")
countdata = read.table(filename2, header=TRUE, sep=",")

rownames(coldata) <- coldata[,1]
coldata= coldata[,2:ncol(coldata)]

coldata$Replicate <- factor(coldata$Replicate)
coldata$Condition <- factor(coldata$Condition)

#Some GeneID are repeated in the data set. I add a number to the gene ID

countdata$NewID <- paste(countdata$X, countdata$X.1, sep="_")

countdata<- subset(countdata, select = -c(X, X.1))

new_order = sort(colnames(countdata))

#Data set oredered in a correct way
countdata <- countdata[, new_order]

rownames(countdata) <- countdata[,1]
countdata = countdata[,2:ncol(countdata)]

#if everything is TRUE is a good sign 
rownames(coldata)==colnames(countdata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Eliminate values cero from coundata
id <- which(apply(countdata, 1, function (x) all(abs(x) >= 1)))

countdata <- countdata[id, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Normalize the data with the total reads of each pool and condition CPM
sum_col_data = colSums(countdata)

countdata_div_total_reads  = sweep(countdata,2,sum_col_data,`/`)*1000000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Dispersion

dispersions <- apply(countdata_div_total_reads, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
  })
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fitting Models
#To fit a model for one gene over time we use a formula such as:
#gene expression ~ fixed effects + random effects

#It is very important that the countdata and dispersion have the same name: if (!all(rownames(countdata) %in% names(dispersion))) {
      #stop("Some dispersion values are missing")

#nrow(countdata)
nrow(countdata_div_total_reads)
length(dispersions)

results <- glmmSeq(~ Time * Condition + (1 | Condition),
                   countdata = countdata_div_total_reads,
                   metadata = coldata,
                   dispersion = dispersions,
                   progress = TRUE)


stats <- summary(results)

#ordered list according to p values 

stats_ordered = stats[order(stats[, 'P_Time:Condition']),]

write.csv(stats_ordered, file="stats_ordered_POOL4_9_glmm_CPM.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#For variables such as time, which are matched according to an ID (the random effect), 
#we can examine the fitted model using plots which show estimated means and confidence intervals 
#based on coefficients for the fitted regression model, overlaid upon the underlying data. 
#In this case the samples are matched longitudinally over time.
#vector with the names of the genes ordered by P_Time:Condition (p-value)
gene_name_ordered= rownames(stats_ordered)

for(i in 1:length(gene_name_ordered)) {
plotColours <- c("skyblue", "goldenrod1")
modColours <- c("dodgerblue3", "goldenrod3")
shapes <- c(17, 19)  
ggmodelPlot(results,
            geneName = gene_name_ordered[i],
            x1var = "Time",
            x2var="Condition",
            xlab="Time",
            colours = plotColours,
            shapes = shapes,
            lineColours = plotColours, 
            modelColours = modColours,
            modelSize = 10)
file_name = gene_name_ordered[i]
#file_name = paste(i,"_",file_name,".pdf", sep="")
file_name = paste(i,"_",file_name,"POOL3",".png", sep="")

ggsave(file_name)
}

