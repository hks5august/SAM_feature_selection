
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("impute")

install.packages('samr')

BiocManager::install("multtest")



library(multtest)
library(samr)

#setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Radiation_responsive_GBM/SF2_val_GSE7505/GSE7505")
setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Radiation_responsive_GBM/SF2_val_GSE7505/GSE7505/SAM_analysis/")
set.seed(7)
#output <- SAM(x,y,resp.type="Quantitative",nperms=50,fdr.output=.5)
#x: a matrix, a data frame, or an ExpressionSet object. Each row of data (or exprs(data), respectively) must correspond to a variable (e.g., a gene), and each column to a sample (i.e.\ an observation).
#y: outcome


#data <- read.table("CNS_cell_line_data_with_SF2_values_tpose.txt", header=T, sep="\t", row.names=1, check.names=F)
#data <- read.table("2446_genes_CNS_cell_line_data_with_gene_symb_data_without_SD_tpose_with_SF2_tpose", header=T, sep="\t", row.names=1, check.names=F)
#data <- read.table("2446_genes_CNS_cell_line_data_with_gene_symb_data_without_SD_tpose_with_SF2_tpose", header=T, sep="\t", row.names=1, check.names=F)
data <- read.table("2438_genes_CNS_cell_line_data_with_gene_symb_data_without_SD_tpose_with_SF2_tpose_copy.txt", header=T, sep="\t", row.names=1, check.names=F)


head(data[1:6],2)

out <-  read.table("SF2_val_data.txt", header=T, sep="\t", row.names=1, check.names=F)
head(out,2)

#library(dplyr)
#out1 <- out %>% mutate_at(vars(SF2), ~ as.integer(round(.x)))


out1 <- round(out$SF2,3)

head(out1)


############### using samr method #################


data2=list(x=as.matrix(data),y=out1, geneid=row.names(data), genenames=row.names(data), logged2=TRUE)

data2
#?samr
samr.obj2<-samr(data2,  resp.type="Quantitative", assay.type= "array",  nperms=1000, regression.method="standard",  random.seed=7)
samr.obj2

samr.obj2


delta=0.5
samr.plot(samr.obj2,delta)

jpeg(file="SAMR_plot.jpeg", units="in", width=10, height=10, res=300)
samr.plot(samr.obj2,delta)
dev.off()

delta.table2 <- samr.compute.delta.table(samr.obj2)

siggenes.table2<-samr.compute.siggenes.table(samr.obj2,delta, data2, delta.table2)
siggenes.table2

pos_corr_genes <- as.data.frame(siggenes.table2$genes.up)
pos_corr_genes
neg_corr_genes <- as.data.frame(siggenes.table2$genes.lo)
neg_corr_genes
dim(pos_corr_genes)
dim(neg_corr_genes)

#write into a file
write.table(pos_corr_genes,file="pos_corr_genes.txt", sep='\t',  quote = F,row.names = FALSE)
write.table(neg_corr_genes,file="neg_corr_genes.txt", sep='\t',  quote = F,row.names = FALSE)


################ using SAM method ######################

output <- SAM(as.matrix(data), out1, resp.type="Quantitative", nperms=50, fdr.output=.5)
output$siggenes.table
output1 <- SAM(as.matrix(data), out1, resp.type="Quantitative", nperms=50, 
               fdr.output=.5, regression.method="standard")

output1$siggenes.table
output1 <- SAM(as.matrix(data), out1, resp.type="Quantitative", nperms=100,
               genenames=row.names(data), geneid=row.names(data),logged2=T,
               fdr.output=.05, regression.method="standard" ,  random.seed=27)
output1$siggenes.table

data2



summary(output)
summary(output1)

# Generate the Delta plots for Delta = 0.2, 0.4, 0.6, ..., 2
plot(output, seq(0.2, 0.4, 2))

# Obtain the SAM plot for Delta = 2
plot(output, 2)


# Get information about the genes called significant using 
# Delta = 3.
sam.sum3 <- summary(output, 3, entrez=FALSE)

# Obtain the rows of golub containing the genes called
# differentially expressed
sam.sum3@row.sig.genes


delta=.4 
samr.plot(output,delta)
delta.table<-samr.compute.delta.table(output) 
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta,data,delta.table)


############ with full data #########

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Radiation_responsive_GBM/SF2_val_GSE7505/GSE7505/SAM_analysis/full_data")

set.seed(7)

Full_data <- read.table("CNS_cell_line_data_with_gene_symb_data_without_SD_tpose_with_all_SF_vals_tpose", header=T, sep="\t", row.names=1, check.names=F)


out <-  read.table("../SF2_val_data.txt", header=T, sep="\t", row.names=1, check.names=F)

out1 <- round(out$SF2,3)

head(out1)


############### using samr method #################


Full_data2 <- list(x=as.matrix(Full_data),y=out1, geneid=row.names(Full_data), genenames=row.names(Full_data), logged2=TRUE)

Full_data2
#?samr
Full_samr.obj <- samr(Full_data2,  resp.type="Quantitative", assay.type= "array",  nperms=1000, regression.method="standard",  random.seed=7)
Full_samr.obj


output_f <- SAM(as.matrix(Full_data), out1, resp.type="Quantitative", nperms=1000, 
               fdr.output=.5, regression.method="standard")

output_f$siggenes.table



delta=0.5
samr.plot(Full_samr.obj,delta)

jpeg(file="Full_SAMR_plot.jpeg", units="in", width=10, height=10, res=300)
samr.plot(Full_samr.obj,delta)
dev.off()

delta.table2 <- samr.compute.delta.table(Full_samr.obj)

siggenes.table2_f<-samr.compute.siggenes.table(Full_samr.obj,delta, Full_data2, delta.table2)
siggenes.table2_f

pos_corr_genes_f <- as.data.frame(siggenes.table2_f$genes.up)
pos_corr_genes_f
neg_corr_genes_f <- as.data.frame(siggenes.table2_f$genes.lo)
neg_corr_genes_f
dim(pos_corr_genes_f)
dim(neg_corr_genes_f)

#write into a file
write.table(pos_corr_genes_f,file="pos_corr_genes_f.txt", sep='\t',  quote = F,row.names = FALSE)
write.table(neg_corr_genes_f,file="neg_corr_genes_f.txt", sep='\t',  quote = F,row.names = FALSE)

