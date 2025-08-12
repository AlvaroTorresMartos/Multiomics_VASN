
#----- READ DATA

DOSAGE_GWAS <- read.table("/home/mireia.bustos/discos/sda/augusto/Descargas/GWAS_PUBMEP/gwas/PLINK/DOSAGE_GWAS_PUBMEP_MIND_0_2_GENO_0_05_MAF_0_01_HWE_0_0001.raw",header=T)

MAP_GWAS <- read.table("/home/mireia.bustos/discos/sda/augusto/Descargas/GWAS_PUBMEP/gwas/PLINK/GWAS_PUBMEP_MIND_0_2_GENO_0_05_MAF_0_01_HWE_0_0001.map",header=F)


MAP_GWAS[,2] <- gsub("-",".",MAP_GWAS[,2])
MAP_GWAS[,2] <- gsub(":",".",MAP_GWAS[,2]) 


#save minor alleles
SNPS_and_minor_alelle <- colnames(DOSAGE_GWAS)

colnames(DOSAGE_GWAS) <- gsub("_A","",colnames(DOSAGE_GWAS))

colnames(DOSAGE_GWAS) <- gsub("_C","",colnames(DOSAGE_GWAS))

colnames(DOSAGE_GWAS) <- gsub("_G","",colnames(DOSAGE_GWAS))

colnames(DOSAGE_GWAS) <- gsub("_T","",colnames(DOSAGE_GWAS))


colnames(DOSAGE_GWAS)[7:106] == MAP_GWAS[1:100,2]

table(colnames(DOSAGE_GWAS)[7:length(colnames(DOSAGE_GWAS))] == MAP_GWAS[,2])

colnames(DOSAGE_GWAS)[7:length(colnames(DOSAGE_GWAS))] <- MAP_GWAS[,2]


#chromosome (1-22, 23=X, 24=Y, 25=PAR region, 26=Mitochondrial or 0 if unplaced)


#ELIMINAR del dosageGWAS las ROWS DE NOSOTROS
#ELIMINAR del dosageGWAS y del MAP  los SNPS DE CHRs 0,23:26

codes_gwas <- DOSAGE_GWAS[1:138,1]
codes_gwas


rownames(DOSAGE_GWAS) <- DOSAGE_GWAS[,1]


DOSAGE_GWAS_ <- DOSAGE_GWAS[1:138,7:length(colnames(DOSAGE_GWAS))]
MAP_GWAS_ <- MAP_GWAS[-which(MAP_GWAS[,1] %in% c(0,23:26)),]
DOSAGE_GWAS_ <- DOSAGE_GWAS_[,MAP_GWAS_[,2]]

dim(DOSAGE_GWAS_)
#[1]   138 471192

dim(MAP_GWAS_)
#[1] 471192     4

############ ANALISIS MQTL on LIST1

#La generacion de las matrices counts se ha extraido de: https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#preparing-count-matrices-from-bam-files





## Import required libraries in R

library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("factoextra")





## Cargamos datos phenotype:

BdWide <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Escritorio/BBDD_PUBMEP/BBDD_PubmeP_RESULTADO/BASES_METILACION/BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_02_06_2020_EPIC.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) 
rownames(BdWide) = BdWide$Code_new_T2


#lista90inds <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/RESUMEN_RESULTS_EPIGENETICA_2021/2021/LISTAS_OMICAS/lista_LONGI_epicarray.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) #Datos input en formato wide
#lista90inds[,2]

#pheno_data  <- BdWide[which(BdWide$Code_new_T2 %in% lista90inds[,2]),]
#rownames(pheno_data) <- pheno_data$Code_new_T2
pheno_data  <- BdWide[which(BdWide$Code_new_T2 %in% codes_gwas),]
pheno_data$id <- rownames(pheno_data)
dim(pheno_data)
#pheno_data
#pheno_data$Cole_T1[which(pheno_data$Cole_T1 == "Obeso_T1")] <- "Sobrepeso_T1"
class(pheno_data$Cole_T1)
pheno_data$Cole_T1 <- as.character(pheno_data$Cole_T1)
pheno_data$Sex_T1 <- as.character(pheno_data$Sex_T1)
pheno_data$Origen_T1 <- as.character(pheno_data$Origen_T1)

dim(pheno_data)
pheno_data = pheno_data[sort(row.names(pheno_data)),]

names(DOSAGE_GWAS_)


# Carga de datos de Epigenetica:

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k


ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

mVals_BMIQ = read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/mVals_BMIQ_NEW.csv",row.names=1,header=T)
ann850k = ann850k[rownames(mVals_BMIQ),]

codespub <- sort(as.character(read.table("/home/mireia.bustos/discos/sdb/mQTLs_NEWLIST_2024/PUB_CODES.txt",row.names=1)[1,]))

#EPIGENETICA_Array <- mVals_BMIQ[rownames(ann850k),gsub("-",".",pheno_data$Code_new_T2)[gsub("-",".",pheno_data$Code_new_T2) %in% colnames(mVals_BMIQ)]]

EPIGENETICA_Array <- mVals_BMIQ[rownames(ann850k),codes_gwas]
EPIGENETICA_Array = EPIGENETICA_Array[,sort(colnames(EPIGENETICA_Array))]

table(rownames(ann850k) == rownames(EPIGENETICA_Array))

#table(rownames(ann850k) == rownames(EPIGENETICA_Array))
#EPIGENETICA_Array = EPIGENETICA_Array[,order(colnames(EPIGENETICA_Array))]

pheno_data$Code_new_T2 = gsub("-",".",pheno_data$Code_new_T2)
pheno_data = pheno_data[pheno_data$Code_new_T2 %in% colnames(EPIGENETICA_Array),]


table(colnames(EPIGENETICA_Array) == pheno_data$Code_new_T2)
#colnames(EPIGENETICA_Array) = rownames(pheno_data)


#lista_validacion <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/RESUMEN_RESULTS_EPIGENETICA_2021/2021/VALIDACION_LISTA_2021/LISTA1.csv", header=T, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))

lista_validacion <- read.csv2("/home/mireia.bustos/ANALYSIS/REPORTS/mQTLs_eQTLS_NEWLIST_2024_PUBMEP/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.csv", header=T, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))

#lista_validacion <- lista_validacion[,2:7]
#lista_validacion[,1]

table(unique(lista_validacion$Name) %in% rownames(EPIGENETICA_Array))


#iniCovPub = read.table("/home/mireia.bustos/discos/sda/augusto/Descargas/GWAS_PUBMEP/GWAS_mQTL/PUBERTAL/covariates.txt")
#write.table(iniCovPub,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/PUB_CODES.txt",row.names=F,col.names=F)


DOSAGE_GWAS_ <- DOSAGE_GWAS_[codespub,]

pheno_data <- pheno_data[codespub,]
#rownames(pheno_data) = pheno_data$Code_new_T2

# Load generated files: (esto se lanza si he cerrado sesion)

snv_dataframe_hg38 <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/ANOTACION_HG38_EPICARRAY/EpicArray_Annotation_HG38.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 





# Prepare annotation files to MatrixEQTL:

EpicArray_Annotation_HG38_MatrixEQTL <- snv_dataframe_hg38[,c(2,3,4)]
rownames(EpicArray_Annotation_HG38_MatrixEQTL) <- snv_dataframe_hg38[,2]
names(EpicArray_Annotation_HG38_MatrixEQTL) <- c("snp","chr","pos")
EpicArray_Annotation_HG38_MatrixEQTL[,2] <- as.character(EpicArray_Annotation_HG38_MatrixEQTL[,2])
EpicArray_Annotation_HG38_MatrixEQTL <- EpicArray_Annotation_HG38_MatrixEQTL[order(EpicArray_Annotation_HG38_MatrixEQTL$chr,EpicArray_Annotation_HG38_MatrixEQTL$pos),]
head(EpicArray_Annotation_HG38_MatrixEQTL)

sort(unique(EpicArray_Annotation_HG38_MatrixEQTL$chr) )

MAP_GWAS_[,1] <- paste("chr",MAP_GWAS_[,1],sep="")



# Prepare MatrixEQTL input files:

pheno_data$Sex_T2 <- as.numeric(as.factor(pheno_data$Sex_T1))
covariates <- as.data.frame(t(pheno_data))
covariates$id <- rownames(covariates)
covariates <- covariates[c("id","Origen_T1","Sex_T1","Age_T2"),c(103,1:102)]
colnames(covariates) <- gsub("sample_","",colnames(covariates))
head(covariates)
dim(covariates)

EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array[rownames(EpicArray_Annotation_HG38_MatrixEQTL),]


EPIGENETICA_Array_MatrixEQTL$id <- rownames(EPIGENETICA_Array_MatrixEQTL)
EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array_MatrixEQTL[,c(139,1:138)]


DOSAGE_GWAS <- as.data.frame(t(DOSAGE_GWAS_))
DOSAGE_GWAS$id <- rownames(DOSAGE_GWAS)
DOSAGE_GWAS <- DOSAGE_GWAS[,c(139,1:138)]

######ojo
#
#
#

table(colnames(DOSAGE_GWAS) == colnames(EPIGENETICA_Array_MatrixEQTL))
 
dim(EPIGENETICA_Array_MatrixEQTL)

dim(covariates)

dim(EpicArray_Annotation_HG38_MatrixEQTL)

MAP_GWAS_ <- MAP_GWAS[,c(2,1,4)]
names(MAP_GWAS_) <- c("snp","chr","pos")

EpicArray_Annotation_HG38_MatrixEQTL$s2 <- EpicArray_Annotation_HG38_MatrixEQTL$pos

names(EpicArray_Annotation_HG38_MatrixEQTL) <- c("geneid","chr","s1","s2")

head(EpicArray_Annotation_HG38_MatrixEQTL)


write.table(covariates,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_covariates.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE, quote=F)

write.table(EPIGENETICA_Array_MatrixEQTL,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)


write.table(DOSAGE_GWAS_,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_gwas.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)

write.table(MAP_GWAS_,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_gwas_location.txt", append = FALSE, sep = "\t", dec = ".",row.names = F, col.names = TRUE, quote=F)

write.table(EpicArray_Annotation_HG38_MatrixEQTL,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_location.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)


























######## ######## ######## ######## ######## ######## ######## ########
######## RUN MATRIXEQTL with only validation cpgs
######## ######## ######## ######## ######## ######## ######## ########
lista_validacion <- read.csv2("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.csv")

EPIGENETICA_Array_MatrixEQTL = read.table("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica.txt", header=T)
rownames(EPIGENETICA_Array_MatrixEQTL) = EPIGENETICA_Array_MatrixEQTL$id

EpicArray_Annotation_HG38_MatrixEQTL = read.table("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_location.txt",header=T)
rownames(EpicArray_Annotation_HG38_MatrixEQTL) = EpicArray_Annotation_HG38_MatrixEQTL$geneid

VALIDACION_EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array_MatrixEQTL[unique(lista_validacion$Name),]

VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL <- EpicArray_Annotation_HG38_MatrixEQTL[unique(lista_validacion$Name),]


write.table(VALIDACION_EPIGENETICA_Array_MatrixEQTL,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)

write.table(VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_location_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)


table(rownames(VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL) == rownames(VALIDACION_EPIGENETICA_Array_MatrixEQTL))

dim(VALIDACION_EPIGENETICA_Array_MatrixEQTL)
dim(EpicArray_Annotation_HG38_MatrixEQTL)



# 2. Test local and distand gene-SNP pairs separately and plot Q-Q plots of local and distant p-values
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_gwas.txt";
snps_location_file_name = "/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_gwas_location.txt";

# Gene expression file name
expression_file_name = "/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_VALIDACIONlist.txt";
gene_location_file_name = "/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_epigenetica_location_VALIDACIONlist.txt";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/PUB_covariates.txt";

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 0.5/(138*471192);

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 5e2;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

gene

#Outliers in expression. Quantile normalization. Outliers in expression data are usually harder to deal with. The accepted remedy by the GTEx consortium is the transformation of the measurements for each gene into normally distributed while preserving relative rankings. The target distribution may be the standard normal distribution or the normal distribution the mean and spread of the original measurements. Here is the code for such transformation:
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);
gene


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

#No significant associations were found.
#3
#Task finished in 4.383 seconds


unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#Analysis done in:  8.194  seconds 
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
#     snps       gene statistic     pvalue FDR      beta
#1 rs12282 cg21561989  2.194844 0.03048858   1 0.4441128

cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
#[1] snps      gene      beta      statistic pvalue    FDR      
#<0 rows> (or 0-length row.names)


## Plot the Q-Q plot of local and distant p-values


pdf("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_PUB_QQPlot_and_distant_pvals.pdf")
plot(me)
dev.off()

png("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_PUB_QQPlot_and_distant_pvals.png")
plot(me)
dev.off()


cpg_cis <- unique((me$cis$eqtls)[,2])
cpg_trans <- unique((me$trans$eqtls)[,2])




library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


write.csv2(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22,23)],"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_PUB_cpg_CIS_sig_annotated.csv")



write.csv2(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22,23)],"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_PUB_cpg_TRANS_sig_annotated.csv")



library("AnnotationDbi")

library("Homo.sapiens")


me$cis$eqtls[,2] <- paste(me$cis$eqtls[,2],gsub(";.*","",ann850k[me$cis$eqtls[,2],23]),sep="_")
meCiseqtls_genes = sub(".*?_","",me$cis$eqtls[,2])


me$trans$eqtls[,2] <- paste(me$trans$eqtls[,2],gsub(";.*","",ann850k[me$trans$eqtls[,2],23]),sep="_")
meTranseqtls_genes = sub(".*?_","",me$trans$eqtls[,2])



symbols_RNAseq_cis <- mapIds(Homo.sapiens,
                     keys=meCiseqtls_genes,
                     column="SYMBOL",
                     keytype="REFSEQ",
                     multiVals="first",
                     fuzzy=T)

symbols_RNAseq_trans <- mapIds(Homo.sapiens,
                     keys=meTranseqtls_genes,
                     column="SYMBOL",
                     keytype="REFSEQ",
                     multiVals="first")


me$trans$eqtls[,2] <- paste(me$trans$eqtls[,2],symbols_RNAseq_trans,sep="_")

me$cis$eqtls[,2] <- paste(me$cis$eqtls[,2],symbols_RNAseq_cis,sep="_")





write.csv2(me$cis$eqtls,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_mQTLS_cis_mainOUTPUT.csv")

#mQTLS_cis_mainOUTPUT_puber

write.csv2(me$trans$eqtls,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/mQTL/PUB/RESULT/1.3_NEW_LIST_p0.3_mQTLS_trans_mainOUTPUT.csv")



