#La generacion de las matrices counts se ha extraido de: https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#preparing-count-matrices-from-bam-files





## Import required libraries in R

library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("factoextra")
library(biomaRt)





## Cargamos datos phenotype:

#BdWide <- read.csv2("/home/augusto/Escritorio/BBDD_PUBMEP/BBDD_PubmeP_RESULTADO/BASES_METILACION/COPIASEGURIDAD_19_04_2021/BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_19_04_2021_EPIC.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) 
BdWide <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Escritorio/BBDD_PUBMEP/BBDD_PubmeP_RESULTADO/BASES_METILACION/COPIASEGURIDAD_19_04_2021/BASE_PUBMEP_LONGI_WIDEformat_NIÑAS_Y_NIÑOS_213inds_19_04_2021_EPIC.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) 


dim(BdWide)

CODESPAXGENE <- c("MR021","MR017","MR020","MR025","MR015","MR014","MR027","MR024","MR023","MR019","MR022","MR018","MR016","MR026","MS297","MS299","MS293","MS305","MS303","MS295","MS296","MS289","MS287","MS291","MS300","MS298","MS288","MS302","MS301","MS290","MS304","MS292","Z480","Z416","Z419","Z484","Z421","Z433","Z432","Z422","Z486","Z487","Z412","Z411")
CODESPAXGENE

pheno_data  <- BdWide[which(BdWide$Code_new_T2 %in% CODESPAXGENE),]
pheno_data  <- pheno_data[order(pheno_data$Code_new_T2),]
pheno_data$Code_new_T2
rownames(pheno_data) <- pheno_data$Code_new_T2
rownames(pheno_data) <- paste("sample_",rownames(pheno_data),sep="")
pheno_data$ids <- rownames(pheno_data)
dim(pheno_data)
pheno_data <- pheno_data[,c(3949,1:3948)]
pheno_data
#pheno_data$Cole_T2[which(pheno_data$Cole_T2 == "Obeso_T2")] <- "Sobrepeso_T2"
class(pheno_data$Cole_T2)
pheno_data$Cole_T2 <- as.character(pheno_data$Cole_T2)
pheno_data$Sex_T2 <- as.character(pheno_data$Sex_T2)
pheno_data$Origen_T1 <- as.character(pheno_data$Origen_T1)





## Cargamos datos de genoma de referencia:

#gtffile <- file.path("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/HUMAN_GENOME_REFERENCE_HISAT/Homo_sapiens.GRCh38.102.gtf")
gtffile <- file.path("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/HUMAN_GENOME_REFERENCE_HISAT/Homo_sapiens.GRCh38.102.gtf")


txdb <- makeTxDbFromGFF(gtffile, format="gtf")

ebg <- exonsBy(txdb, by="gene")
ebg





## Cargamos archivos BAM de carrera 201109:

#filenames_201109 <- file.path(paste0("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/my_rna_seq_ANALYSIS_1/201109/",pheno_data$Code_new_T2, ".bam"))
filenames_201109 <- file.path(paste0("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/my_rna_seq_ANALYSIS_1/201109/",pheno_data$Code_new_T2, ".bam"))


file.exists(filenames_201109)

bamfiles_201109 <- BamFileList(filenames_201109, yieldSize=2000000)

seqinfo(bamfiles_201109[1])





## Generar matriz de counts carrera 201109: Metodo SummarizeOverlaps 
#(Para hacerlo con este metodo tengo que saber si la librería que usó genyo es strand-specific, en cuyo caso, deberia marcar el argumento ignore.strand como FALSE).

register(MulticoreParam(30))

se_201109 <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles_201109,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE,
                        fragments=TRUE )

#save(se_201109, file = '/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201109/summOv_se_201109_rand.Rdata')
save(se_201109, file = '/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/1_COUNT_MATRIX/summOv_se_201109_rand_2024_12_13.Rdata')




## Generar matriz de counts carrera 201109: Metodo featureCounts
register(MulticoreParam(10))
fc_201109_complete <- featureCounts(files=filenames_201109, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

colnames(fc_201109_complete$counts) <- gsub(".bam","",colnames(fc_201109_complete$counts))
head(fc_201109_complete$counts)

fc_201109 <- fc_201109_complete$counts





## Cargamos archivos BAM de carrera 201116:

filenames_201116 <- file.path(paste0("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/my_rna_seq_ANALYSIS_1/201116/",pheno_data$Code_new_T2, ".bam"))

file.exists(filenames_201116)

bamfiles_201116 <- BamFileList(filenames_201116, yieldSize=2000000)

seqinfo(bamfiles_201116[1])





## Generar matriz de counts carrera 201116: Metodo SummarizeOverlaps 
#(Para hacerlo con este metodo tengo que saber si la librería que usó genyo es strand-specific, en cuyo caso, deberia marcar el argumento ignore.strand como FALSE).

register(MulticoreParam(6))

se_201116 <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles_201116,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE,
                        fragments=TRUE )

#save(se_201116, file = '/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201116/summOv_se_201116_rand.Rdata')





## Generar matriz de counts carrera 201116: Metodo featureCounts

fc_201116_complete <- featureCounts(files=filenames_201116, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

colnames(fc_201116_complete$counts) <- gsub(".bam","",colnames(fc_201116_complete$counts))

head(fc_201116_complete$counts)

fc_201116 <- fc_201116_complete$counts





## Generar PCA PLOT para ver si replicas entre carreras se agrupan: Matrices de count derivadas del metodo featureCounts
##https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html

colnames(fc_201116) == colnames(fc_201109)

rownames(fc_201116) == rownames(fc_201109)

fc_201116_v2 <- fc_201116

fc_201109_v2 <- fc_201109

colnames(fc_201116_v2) <- paste(colnames(fc_201116_v2),"rep1",sep="_")

colnames(fc_201109_v2) <- paste(colnames(fc_201109_v2),"rep2",sep="_")

datos <- rbind(pheno_data,pheno_data) 

countdata <- cbind(fc_201116_v2,fc_201109_v2)

#countdata <- rlog(countdata) esto tarda un rato y no es del todo necesario lanzarlo para comprobar el pca.

pca <- prcomp(as.data.frame(t(countdata)))

get_eigenvalue(pca)

pdf("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/PLOTS_QC/carreras_rnaseq_PUBMEP_2020_biplot_from_featurecount.pdf", width=40/2.54, height=28/2.54)

fviz_pca_ind(pca,col.ind=datos$Code_new_T2,repel=TRUE)

dev.off()





## Generar PCA PLOT para ver si replicas entre carreras se agrupan: Matrices de count derivadas del metodo SummarizeOverlaps 
##https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html

colnames(assay(se_201116)) == colnames(assay(se_201109))

rownames(assay(se_201116)) == rownames(assay(se_201109))

se_201116_v2 <- assay(se_201116)

se_201109_v2 <- assay(se_201109)

colnames(se_201116_v2) <- paste(colnames(se_201116_v2),"rep1",sep="_")

colnames(se_201109_v2) <- paste(colnames(se_201109_v2),"rep2",sep="_")

datos <- rbind(pheno_data,pheno_data) 

countdata_se <- cbind(se_201116_v2,se_201109_v2)

#countdata <- rlog(countdata) esto tarda un rato.

pca_se <- prcomp(as.data.frame(t(countdata_se)))

get_eigenvalue(pca_se)

pdf("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/PLOTS_QC/carreras_rnaseq_PUBMEP_2020_biplot_from_SE.pdf", width=40/2.54, height=28/2.54)

fviz_pca_ind(pca_se,col.ind=datos$Code_new_T2,repel=TRUE)

dev.off()





## Guardar counts matrices de metodo featurecounts:
fc_201116 <- fc_201116[order(rownames(fc_201116)),]
fc_201109 <- fc_201109[order(rownames(fc_201109)),]

#write.csv2(fc_201116,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201116/fc_201116_from_featurecount.csv")
#write.csv2(fc_201109,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201109/fc_201109_from_featurecount.csv")

head(fc_201116)

head(fc_201109)





## Guardar counts matrices de metodo SummarizeOverlaps:
fc_se_201109 <- assay(se_201109)
fc_se_201116 <- assay(se_201116)
fc_se_201116 <- fc_se_201116[order(rownames(fc_se_201116)),]
fc_se_201109 <- fc_se_201109[order(rownames(fc_se_201109)),]

#write.csv2(fc_se_201116,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201116/fc_201116_from_SummarizeOverlaps.csv")
#write.csv2(fc_se_201109,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/CARRERA_201109/fc_201109_from_SummarizeOverlaps.csv")

head(fc_se_201116)

head(fc_se_201109)





## Combinar carreras (media y suma): metodo featurecounts

fc_SUM <- fc_201116 + fc_201109 

X <- list(fc_201116, fc_201109)

Y <- do.call(cbind, X)

Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

fc_MEAN <- apply(Y, c(1, 2), mean, na.rm = TRUE)

fc_MEAN <- as.matrix(fc_MEAN)

colnames(fc_MEAN) <- colnames(fc_SUM)

rownames(fc_MEAN) <- rownames(fc_SUM)

head(fc_MEAN)

head(fc_SUM)

#write.csv2(fc_MEAN,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_MEAN_from_featurecount.csv")

#write.csv2(fc_SUM,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_SUM_from_featurecount.csv")




## Combinar carreras (media y suma): metodo SummarizeOverlaps

fc_se_SUM <- fc_se_201116 + fc_se_201109 

X <- list(fc_se_201116, fc_se_201109)

Y <- do.call(cbind, X)

Y <- array(Y, dim=c(dim(X[[1]]), length(X)))

fc_se_MEAN <- apply(Y, c(1, 2), mean, na.rm = TRUE)

fc_se_MEAN <- as.matrix(fc_se_MEAN)

colnames(fc_se_MEAN) <- colnames(fc_se_SUM)

rownames(fc_se_MEAN) <- rownames(fc_se_SUM)

head(fc_se_MEAN)

head(fc_se_SUM)

#write.csv2(fc_se_MEAN,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_se_MEAN_from_SummarizeOverlaps.csv")

#write.csv2(fc_se_SUM,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_se_SUM_from_SummarizeOverlaps.csv")





# Load generated files (esto solo es necesario si he cerrado la sesión).

fc_MEAN <- read.csv2("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_MEAN_from_featurecount.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 

fc_SUM <- read.csv2("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_SUM_from_featurecount.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 

fc_se_MEAN <- read.csv2("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_se_MEAN_from_SummarizeOverlaps.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 

fc_se_SUM <- read.csv2("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/COUNT_MATRIX/COMBINADAS/fc_se_SUM_from_SummarizeOverlaps.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 

geneids <- rownames(fc_se_MEAN)





#Para obtener los nombres de los genes:

library("AnnotationDbi")

library("Homo.sapiens")

symbols_RNAseq <- mapIds(Homo.sapiens,
                     keys=rownames(fc_se_MEAN),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")





# Annotate ensemble genes: TODA LA INFORMACION SOBRE COMO HACER QUERIES ESTA EN https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") #se pueden modificar la version y la especie con distintos comandos.

filters = listFilters(ensembl)

filters

filters[1:15,]

attributes = listAttributes(ensembl)

attributes

attributes[1:15,]

searchFilters(mart = ensembl, pattern = "ensembl.*id")

geneids <- rownames(fc_se_MEAN)

RNAseq_Annotation_HG38 <- getBM(attributes = attributes[c(1:2,9:12),1],
      filters = 'ensembl_gene_id',
      values = geneids, 
      mart = ensembl)

dim(RNAseq_Annotation_HG38)

length(geneids)

table(RNAseq_Annotation_HG38[,1] == geneids)

#write.csv2(RNAseq_Annotation_HG38,"/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/OUTPUTS/ANNOTATION_COUNT_MATRIX/RNAseq_Annotation_HG38.csv")

#write.csv2(RNAseq_Annotation_HG38,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/ANNOTATION_COUNT_MATRIX/RNAseq_Annotation_HG38.csv")




# Carga de datos de Epigenetica:

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k

mVals_BMIQ <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/PRIMER_ANALISIS/BASES_DATOS_METILACION_GENERADAS_PREPROCESADAS/mVals_BMIQ_NEW.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)

ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)

ann850k = ann850k[rownames(mVals_BMIQ),]

EPIGENETICA_Array <- mVals_BMIQ[rownames(ann850k),colnames(fc_MEAN)]

table(rownames(ann850k) == rownames(EPIGENETICA_Array))

#lista_validacion <- read.csv2("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.csv")
lista_validacion <- read.csv2("/home/mireia.bustos/ANALYSIS/REPORTS/mQTLs_eQTLS_NEWLIST_2024_PUBMEP/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.csv", header=T, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))

# lista_validacion <- lista_validacion[,2:7]
# lista_validacion[,1]

table(unique(lista_validacion$Name) %in% rownames(EPIGENETICA_Array))





# Liftover: Cambiar anotacion Epigenetica From hg19 to hg38.

library(rtracklayer)
library(liftOver)

from_hg19_to_hg38 = import.chain("/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/HUMAN_GENOME_REFERENCE_HISAT/hg19ToHg38.over.chain")

snv_dataframe_hg19 <- ann850k[,c(1,2,2,3,4,22)]
names(snv_dataframe_hg19) <- c("chrom","start","end","strand","class","genesymbol")

snv_dataframe_hg19 <- GRanges(snv_dataframe_hg19) # we need to convert our dataframe into a granges object

seqlevelsStyle(snv_dataframe_hg19) = "UCSC"  # required


snv_dataframe_liftover <- liftOver(snv_dataframe_hg19, from_hg19_to_hg38)

snv_dataframe_hg38 <- as.data.frame(snv_dataframe_liftover) # you need to convert it again to dataframe
snv_dataframe_hg38

table(le$Name %in% snv_dataframe_hg38[,2])

# To download more chain files, go to http://hgdownload.soe.ucsc.edu/downloads.html#liftover, find your original genome assembly version and click on "LiftOver files".

#write.csv2(snv_dataframe_hg38,"/home/mireia.bustos/discos/sda/augusto/Descargas/ALL_EPIC_ARRAY_DATA_PUBMEP/ANOTACION_HG38_EPICARRAY/EpicArray_Annotation_HG38.csv")
#write.csv2(snv_dataframe_hg38,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/ANOTACION_HG38_EPICARRAY/EpicArray_Annotation_HG38.csv")



# Load generated files: (esto se lanza si he cerrado sesion)

snv_dataframe_hg38 <- read.csv2("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/ANOTACION_HG38_EPICARRAY/EpicArray_Annotation_HG38.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 

RNAseq_Annotation_HG38 <- read.csv2("home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/ANNOTATION_COUNT_MATRIX/RNAseq_Annotation_HG38.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"), row.names = 1) 





# Prepare annotation files to MatrixEQTL:

EpicArray_Annotation_HG38_MatrixEQTL <- snv_dataframe_hg38[,c(2,3,4)]
rownames(EpicArray_Annotation_HG38_MatrixEQTL) <- snv_dataframe_hg38[,2]
names(EpicArray_Annotation_HG38_MatrixEQTL) <- c("snp","chr","pos")
EpicArray_Annotation_HG38_MatrixEQTL[,2] <- as.character(EpicArray_Annotation_HG38_MatrixEQTL[,2])
EpicArray_Annotation_HG38_MatrixEQTL <- EpicArray_Annotation_HG38_MatrixEQTL[order(EpicArray_Annotation_HG38_MatrixEQTL$chr,EpicArray_Annotation_HG38_MatrixEQTL$pos),]
head(EpicArray_Annotation_HG38_MatrixEQTL)

RNAseq_Annotation_HG38_MatrixEQTL <- RNAseq_Annotation_HG38[,c(1,3,4,5)]
RNAseq_Annotation_HG38_MatrixEQTL[,2] <- paste("chr",RNAseq_Annotation_HG38_MatrixEQTL[,2],sep="")
names(RNAseq_Annotation_HG38_MatrixEQTL) <- c("geneid","chr","s1","s2")
RNAseq_Annotation_HG38_MatrixEQTL <- RNAseq_Annotation_HG38_MatrixEQTL[order(RNAseq_Annotation_HG38_MatrixEQTL$chr,RNAseq_Annotation_HG38_MatrixEQTL$s1),]
RNAseq_Annotation_HG38_MatrixEQTL <- RNAseq_Annotation_HG38_MatrixEQTL[-which(!RNAseq_Annotation_HG38_MatrixEQTL[,2] %in% unique(EpicArray_Annotation_HG38_MatrixEQTL[,2])),]
table(RNAseq_Annotation_HG38_MatrixEQTL[,2])
head(RNAseq_Annotation_HG38_MatrixEQTL)

unique(EpicArray_Annotation_HG38_MatrixEQTL$chr) == unique(RNAseq_Annotation_HG38_MatrixEQTL[,2])





# Total de reads por archivo (expresado en millones):

((sum(colSums(fc_MEAN[,])))/(44))/1000000
#[1] 10.53729
((sum(colSums(fc_SUM[,])))/(44))/1000000
#[1] 21.07458
((sum(colSums(fc_se_MEAN[,])))/(44))/1000000
#[1] 8.918756
((sum(colSums(fc_se_SUM[,])))/(44))/1000000
#[1] 17.83751





# Prepare MatrixEQTL input files:

pheno_data$Sex_T2 <- as.numeric(as.factor(pheno_data$Sex_T2))
covariates <- as.data.frame(t(pheno_data))
covariates$id <- rownames(covariates)
covariates <- covariates[c("Origen_T1","Sex_T2","Age_T2"),c(45,1:44)]
colnames(covariates) <- gsub("sample_","",colnames(covariates))
head(covariates)



EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array[rownames(EpicArray_Annotation_HG38_MatrixEQTL),]
RNAseq_array_fc_MEAN_MatrixEQTL <- as.data.frame(fc_MEAN[RNAseq_Annotation_HG38_MatrixEQTL[,1],])
RNAseq_array_fc_SUM_MatrixEQTL <- as.data.frame(fc_SUM[RNAseq_Annotation_HG38_MatrixEQTL[,1],])
RNAseq_array_fc_se_MEAN_MatrixEQTL <- as.data.frame(fc_se_MEAN[RNAseq_Annotation_HG38_MatrixEQTL[,1],])
RNAseq_array_fc_se_SUM_MatrixEQTL <- as.data.frame(fc_se_SUM[RNAseq_Annotation_HG38_MatrixEQTL[,1],])
colnames(RNAseq_array_fc_se_MEAN_MatrixEQTL) <- gsub(".bam","",colnames(RNAseq_array_fc_se_MEAN_MatrixEQTL))
colnames(RNAseq_array_fc_se_SUM_MatrixEQTL) <- gsub(".bam","",colnames(RNAseq_array_fc_se_SUM_MatrixEQTL))

EPIGENETICA_Array_MatrixEQTL$id <- rownames(EPIGENETICA_Array_MatrixEQTL)
EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array_MatrixEQTL[,c(45,1:44)]

RNAseq_array_fc_MEAN_MatrixEQTL$id <- rownames(RNAseq_array_fc_MEAN_MatrixEQTL)
RNAseq_array_fc_MEAN_MatrixEQTL <- RNAseq_array_fc_MEAN_MatrixEQTL[,c(45,1:44)]
RNAseq_array_fc_SUM_MatrixEQTL$id <- rownames(RNAseq_array_fc_SUM_MatrixEQTL)
RNAseq_array_fc_SUM_MatrixEQTL <- RNAseq_array_fc_SUM_MatrixEQTL[,c(45,1:44)]
RNAseq_array_fc_se_MEAN_MatrixEQTL$id <- rownames(RNAseq_array_fc_se_MEAN_MatrixEQTL)
RNAseq_array_fc_se_MEAN_MatrixEQTL <- RNAseq_array_fc_se_MEAN_MatrixEQTL[,c(45,1:44)]
RNAseq_array_fc_se_SUM_MatrixEQTL$id <- rownames(RNAseq_array_fc_se_SUM_MatrixEQTL)
RNAseq_array_fc_se_SUM_MatrixEQTL <- RNAseq_array_fc_se_SUM_MatrixEQTL[,c(45,1:44)]

colnames(RNAseq_array_fc_se_SUM_MatrixEQTL) == colnames(EPIGENETICA_Array_MatrixEQTL)

rownames(RNAseq_array_fc_se_SUM_MatrixEQTL) == rownames(RNAseq_array_fc_SUM_MatrixEQTL)


dim(EPIGENETICA_Array_MatrixEQTL)
#[1] 809230     44
 
dim(RNAseq_array_fc_MEAN_MatrixEQTL)
#[1] 60058    44
dim(RNAseq_array_fc_SUM_MatrixEQTL)
#[1] 60058    44
dim(RNAseq_array_fc_se_MEAN_MatrixEQTL)
#[1] 60058    44
dim(RNAseq_array_fc_se_SUM_MatrixEQTL)
#[1] 60058    44




head(RNAseq_array_fc_MEAN_MatrixEQTL)

head(RNAseq_array_fc_SUM_MatrixEQTL)

head(RNAseq_array_fc_se_MEAN_MatrixEQTL)

head(RNAseq_array_fc_se_SUM_MatrixEQTL)

dim(EPIGENETICA_Array_MatrixEQTL)
dim(RNAseq_array_fc_MEAN_MatrixEQTL)
dim(RNAseq_array_fc_SUM_MatrixEQTL)
dim(RNAseq_array_fc_se_MEAN_MatrixEQTL)
dim(RNAseq_array_fc_se_SUM_MatrixEQTL)
dim(covariates)

dim(RNAseq_Annotation_HG38_MatrixEQTL)

dim(EpicArray_Annotation_HG38_MatrixEQTL)

# write.table(covariates,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/covariates.txt", append = FALSE, sep = "\t", dec = ".",row.names = F, col.names = TRUE, quote=F)
# write.table(EPIGENETICA_Array_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(RNAseq_array_fc_MEAN_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcMEAN.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(RNAseq_array_fc_SUM_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcSUM.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(RNAseq_array_fc_se_MEAN_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcse_MEAN.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(RNAseq_array_fc_se_SUM_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcse_SUM.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(RNAseq_Annotation_HG38_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_location.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)
# write.table(EpicArray_Annotation_HG38_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_location.txt", append = FALSE, sep = "\t", dec = ".",            row.names = F, col.names = TRUE, quote=F)




covariates = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/covariates.txt",header=T)

EPIGENETICA_Array_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica.txt", header=T)
  
RNAseq_array_fc_MEAN_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcMEAN.txt", header=T)
  
RNAseq_array_fc_SUM_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcSUM.txt", header=T)
  
RNAseq_array_fc_se_MEAN_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcse_MEAN.txt", header=T)
  
RNAseq_array_fc_se_SUM_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcse_SUM.txt", header=T)
  
RNAseq_Annotation_HG38_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_location.txt", header=T)

EpicArray_Annotation_HG38_MatrixEQTL = read.table("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_location.txt", header=T)

















######## ######## ######## ######## ######## ######## ######## ########
######### MATRIX EQTL CON TODAS LAS CPGS:
######## ######## ######## ######## ######## ######## ######## ########

# 
# 1. Test all gene-SNP pairs and plot a histogram of all p-values
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica.txt";

# Gene expression file name
expression_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcMEAN.txt";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/covariates.txt";

# Output file name
output_file_name = tempfile();
output_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/PUBMEP_eQTLs"

# Only associations significant at this level will be saved
pvOutputThreshold = 0.05/(809230*58824);

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


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

me = Matrix_eQTL_engine(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

pdf("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/0_ALL_NEW_LIST_PREPUB_QQPlot_and_distant_pvals.pdf")
plot(me)
dev.off()



me$all$eqtls[,2] <- paste(me$all$eqtls[,2],symbols_RNAseq_trans,sep="_")




write.csv2(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/0_ALL_eQTLs_cpg_ALL_sig_annotated_10kb.csv")



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## 2. Test local and distand gene-SNP pairs separately and plot Q-Q plots of local and distant p-values
#
## Matrix eQTL by Andrey A. Shabalin
## http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
##
## Be sure to use an up to date version of R and Matrix eQTL.
#
## source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica.txt";
snps_location_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica_location.txt";

# Gene expression file name
expression_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/rnaseq_fcMEAN.txt";
gene_location_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/rnaseq_location.txt";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/covariates.txt";

# Output file name
output_file_name_cis = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/0_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv"
#  tempfile();
output_file_name_tra = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/0_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv"
  #tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05/(809230*60058);
pvOutputThreshold_tra = 0.05/(809230*60058);

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 5e5;

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

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

#plot(me)

pdf("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/eQTLs_NEW_LIST_PREPUB_QQPlot_and_distant_pvals.pdf")
plot(me)
dev.off()


me$trans$eqtls[,2] <- paste(me$trans$eqtls[,2],symbols_RNAseq_trans,sep="_")
names(me$trans$eqtls)[1] = "Name"

me$cis$eqtls[,2] <- paste(me$cis$eqtls[,2],symbols_RNAseq_cis,sep="_")
names(me$cis$eqtls)[1] = "Name" 

cpg_cis <- unique((me$cis$eqtls)[,1])
cpg_trans <- unique((me$trans$eqtls)[,1])

res_cis = merge(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$cis$eqtls, by="Name")
res_trans = merge(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$trans$eqtls, by="Name")


write.csv2(res_cis,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv")
write.csv2(res_trans,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/ALL_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")










######## ######## ######## ######## ######## ######## ######## ########
######## RUN MATRIXEQTL with only validation cpgs
######## ######## ######## ######## ######## ######## ######## ########


library(MatrixEQTL)

rownames(EPIGENETICA_Array_MatrixEQTL) = EPIGENETICA_Array_MatrixEQTL$id
rownames(EpicArray_Annotation_HG38_MatrixEQTL) = EpicArray_Annotation_HG38_MatrixEQTL$snp
new_lista_validacion <- read.csv2("/home/mireia.bustos/ANALYSIS/REPORTS/mQTLs_eQTLS_NEWLIST_2024_PUBMEP/2024_03_13_lista_def_CpGs_augusto_ampliada_288filas_193CpGs_89genes.csv", header=T, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))


#VALIDACION_EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array_MatrixEQTL[unique(lista_validacion[,1]),]
#VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL <- EpicArray_Annotation_HG38_MatrixEQTL[unique(lista_validacion[,1]),]

VALIDACION_EPIGENETICA_Array_MatrixEQTL <- EPIGENETICA_Array_MatrixEQTL[unique(new_lista_validacion$Name),]
VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL <- EpicArray_Annotation_HG38_MatrixEQTL[unique(new_lista_validacion$Name),]

#write.table(VALIDACION_EPIGENETICA_Array_MatrixEQTL,"/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = TRUE, quote=F)
#write.table(EpicArray_Annotation_HG38_MatrixEQTL,"/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica_location_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = TRUE, quote=F)

VALIDACION_EPIGENETICA_Array_MatrixEQTL$id = rownames(VALIDACION_EPIGENETICA_Array_MatrixEQTL)
write.table(VALIDACION_EPIGENETICA_Array_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = TRUE, quote=F)
write.table(EpicArray_Annotation_HG38_MatrixEQTL,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_location_VALIDACIONlist.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = TRUE, quote=F)

table(rownames(VALIDACION_EpicArray_Annotation_HG38_MatrixEQTL) == rownames(VALIDACION_EPIGENETICA_Array_MatrixEQTL))

dim(VALIDACION_EPIGENETICA_Array_MatrixEQTL)

# 1. Test all gene-SNP pairs and plot a histogram of all p-values


# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# 
# # Genotype file name
SNP_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica_VALIDACIONlist.txt";

# # Gene expression file name
expression_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/rnaseq_fcMEAN.txt";
# 
# # Covariates file name
# # Set to character() for no covariates
covariates_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/covariates.txt";
# 
# # Output file name
output_file_name = tempfile();
# 
# # Only associations significant at this level will be saved
#pvOutputThreshold = 0.05/(267*60058);
pvOutputThreshold = 0;
# 
# # Error covariance matrix
# # Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# 
# 
# ## Load genotype data
# 
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# ## Load gene expression data
# 
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
# 
# #Outliers in expression. Quantile normalization. Outliers in expression data are usually harder to deal with. The accepted remedy by the GTEx consortium is the transformation of the measurements for each gene into normally distributed while preserving relative rankings. The target distribution may be the standard normal distribution or the normal distribution the mean and spread of the original measurements. Here is the code for such transformation:
for( sl in 1:length(gene) ) {
 mat = gene[[sl]];
 mat = t(apply(mat, 1, rank, ties.method = "average"));
 mat = qnorm(mat / (ncol(gene)+1));
 gene[[sl]] = mat;
}
rm(sl, mat);
gene
# 
# 
# ## Load covariates
# 
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

# 
# ## Run the analysis
# 
me = Matrix_eQTL_engine(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = output_file_name,
#pvOutputThreshold = pvOutputThreshold,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name);
# 
# ## Results:
# 
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
# 
# ## Plot the histogram of all p-values
# 
# plot(me)
# 



pdf("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/3.2_ALL_NEW_eQTLs_QQPlot_and_distant_pvals.pdf")
plot(me)
dev.off()

png("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/3.2_ALL_NEW_eQTLs_QQPlot_and_distant_pvals.png")
plot(me)
dev.off()

cpg_cis <- unique((me$cis$eqtls)[,1])
cpg_trans <- unique((me$trans$eqtls)[,1])



write.csv2(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.1_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv")
write.csv2(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.1_ALL_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")


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
#SNP_file_name = "/home/mireia.bustos/discos/sdb/MIREIA/RNASEQ_ADIPOSEQ/RNAseq_PUBMEP/MATRIXEQTL/INPUTS/epigenetica_VALIDACIONlist.txt";
SNP_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_VALIDACIONlist.txt"

snps_location_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/epigenetica_location_VALIDACIONlist.txt";

# Gene expression file name
expression_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_fcMEAN.txt";
gene_location_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/rnaseq_location.txt";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/INPUTS/covariates.txt";

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 0.05/(267*60058);

pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 1

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e4;

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

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values


pdf("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1_OnlyTargets_NEW_eQTLs_QQPlot_and_distant_pvals.pdf")
plot(me)
dev.off()

png("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1_OnlyTargets_NEW_eQTLs_QQPlot_and_distant_pvals.png")
plot(me)
dev.off()

# pdf("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_QQPlot_and_distant_pvals.pdf")
# plot(me)
# dev.off()
# 
# png("/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_QQPlot_and_distant_pvals.png")
# plot(me)
# dev.off()



me$trans$eqtls[,2] <- paste(me$trans$eqtls[,2],symbols_RNAseq_trans,sep="_")
names(me$trans$eqtls)[1] = "Name"

me$cis$eqtls[,2] <- paste(me$cis$eqtls[,2],symbols_RNAseq_cis,sep="_")
names(me$cis$eqtls)[1] = "Name" 

cpg_cis <- unique((me$cis$eqtls)[,1])
cpg_trans <- unique((me$trans$eqtls)[,1])

res_cis = merge(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$cis$eqtls, by="Name")
res_trans = merge(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$trans$eqtls, by="Name")


write.csv2(res_cis,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/2.1_TODOS_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv")
write.csv2(res_trans,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/2.1_TODOS_ALL_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")


cpg_cis <- unique((me$cis$eqtls)[,1])
cpg_trans <- unique((me$trans$eqtls)[,1])



write.csv2(res_cis,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.1_TODOS_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv")




#-------------------------
#-------------------------
combined_results = res_cis %>% as.data.frame() %>% arrange(pvalue)
a =  combined_results %>% mutate(UCSC_RefGene_Name = sub(";.*", "", UCSC_RefGene_Name),  # Eliminar todo después (y incluido) del ";"
                                 Name_gene = paste(Name, UCSC_RefGene_Name, sep = "_"))
ba = a[,c(1,15,17:21)] %>% distinct()

# Load packages -------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)

# cpgs = cpgs %>% 
#   gsub(pattern = "_prepubertal", replacement = "") %>% 
#   gsub(pattern = "_pubertal", replacement = "") %>% 
#   unique()

# Import eQTMs results -----
# eqtms = read.csv("./results/eQTM/2024_05_02_mQTLs_eQTLs_with_New_list.knit.csv")
eqtms = read.csv2("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.1_TODOS_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv",row.names = 1)
eqtms = as.data.frame(eqtms) %>% arrange(pvalue) %>% mutate(UCSC_RefGene_Name = sub(";.*", "", UCSC_RefGene_Name), Name_gene = paste(Name, UCSC_RefGene_Name, sep = "_"))
#colnames(eqtms)[2] = "gene"

cpgs = eqtms$Name

eqtms = eqtms %>% 
  dplyr::rename(transcript = gene) %>% 
  tidyr::separate(Name_gene, into = c("cpg", "gene"), sep = "_", remove = FALSE) %>% 
  tidyr::separate(transcript, into = c("rna", "gene2"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(new_cpg = stringr::str_c(gene, cpg, sep = "_"), 
                new_rna = stringr::str_c(gene2, rna, sep = "_"), 
                predictor_label = stringr::str_c(new_cpg, rna, sep = "_")
  ) %>% 
  dplyr::select(-c(FDR))

cols = c("Higher methylation - Highter gene expression" = "#FFD1A9", 
         "Higher methylation - Lower gene expression" = "#B8E6B6", 
         "ns" = "grey")
sizes = c("Higher methylation - Highter gene expression" = 3, 
          "Higher methylation - Lower gene expression" = 3, 
          "ns" = 1.5) 
alphas = c("Higher methylation - Highter gene expression" = 1, 
           "Higher methylation - Lower gene expression" = 1, 
           "ns" = 0.15)


eqtms %>%  
  dplyr::mutate(direction = ifelse(beta > 0 & pvalue < 0.05 & abs(beta) > 0.5, "Higher methylation - Highter gene expression", 
                                   ifelse(beta < 0 & pvalue < 0.05 & abs(beta) > 0.5, "Higher methylation - Lower gene expression", 
                                          "ns")), 
                predictor_label = ifelse(new_cpg %in% cpgs, predictor_label, NA
                )
  ) %>% 
  ggplot(aes(x = beta, y = -log10(pvalue), 
             fill = direction, 
             size = direction, 
             alpha = direction, 
             label = predictor_label)) + 
  geom_point(shape = 21, colour = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.7, color = "darkgrey") + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", size = 0.7, color = "darkgrey") + 
  scale_fill_manual(values = cols) + # point colour
  scale_size_manual(values = sizes) + # point size
  scale_alpha_manual(values = alphas) + # transparency 
  geom_label_repel(max.overlaps = Inf, size = 4, box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50') +
  theme_minimal(base_size = 14) +  # Tamaño base de 14 para mejor legibilidad
  theme(
    plot.caption = element_text(size = 10),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey90"),  # Líneas de cuadrícula principales suaves
    panel.grid.minor = element_blank(),  # Elimina líneas de cuadrícula menores
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Borde negro para el panel
    axis.title = element_text(face = "bold"),  # Ejes en negrita para destacar títulos
    axis.text = element_text(color = "black"),  # Texto del eje en negro para mejor contraste
    legend.position = "none",  # Elimina la leyenda
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Título centrado, negrita y tamaño grande
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12)  # Subtítulo centrado y en negrita para equilibrio visual
  ) +
  labs(title = "eQTMs", 
       # subtitle = "Fórmula: HOMA-IR z-score ~ Metabolito + Hospital + Punto Temporal + (1 | Id) (N = 278)",
       x = "Beta Coefficient",
       y = "-log10(p-value)", 
       fill = "Direction")  + 
  scale_y_continuous(limits = c(0, 2.5),
                     breaks = seq(0, 2.5, by = 0.5)
  ) + 
  scale_x_continuous(
    limits = c(-1.5, 1.5)
  )



#_---------------------------------------
#--------------------------


# 
# write.csv2(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.1_OnlyTargets_NEW_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")

# write.csv2(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_cpg_CIS_sig_annotated_10kb.csv")
# write.csv2(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")
# 



library("AnnotationDbi")
library("Homo.sapiens")


me$cis$eqtls[,1] <- paste(me$cis$eqtls[,1],gsub(";.*","",ann850k[me$cis$eqtls[,1],22]),sep="_")

me$trans$eqtls[,1] <- paste(me$trans$eqtls[,1],gsub(";.*","",ann850k[me$trans$eqtls[,1],22]),sep="_")



symbols_RNAseq_cis <- mapIds(Homo.sapiens,
                     keys=(me$cis$eqtls)[,2],
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

symbols_RNAseq_trans <- mapIds(Homo.sapiens,
                     keys=(me$trans$eqtls)[,2],
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


me$trans$eqtls[,2] <- paste(me$trans$eqtls[,2],symbols_RNAseq_trans,sep="_")
names(me$trans$eqtls)[1] = "Name"

me$cis$eqtls[,2] <- paste(me$cis$eqtls[,2],symbols_RNAseq_cis,sep="_")
names(me$cis$eqtls)[1] = "Name" 

cpg_cis <- unique((me$cis$eqtls)[,1])
cpg_trans <- unique((me$trans$eqtls)[,1])

res_cis = merge(ann850k[cpg_cis,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$cis$eqtls, by="Name")
res_trans = merge(ann850k[cpg_trans,c(1:4,9,18,19,25,26,27,28,29,38,42,22)],me$trans$eqtls, by="Name")


write.csv2(res_cis,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/2.1_TODOS_ALL_eQTLs_cpg_CIS_sig_annotated_10kb.csv")
write.csv2(res_trans,"/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/2.1_TODOS_ALL_eQTLs_cpg_TRANS_sig_annotated_10kb.csv")



# write.csv2(me$cis$eqtls,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_mQTLS_cis_mainOUTPUT_10kb.csv")
# write.csv2(me$trans$eqtls,"/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/mQTLs_NEWLIST_2024/NEW_RESULTS_TABLES/eQTL/2_OLD_eQTLs_mQTLS_trans_mainOUTPUT_10kb.csv")










#transold <- read.csv2("/home/mireia.bustos/discos/sda/augusto/Descargas/ACTIVE_BRAIN/AUGUSTO_METILACION_ACTIVEBRAINS_males/PRIMER_ANALISIS/RESULTADOS_DEF_05_03_2021/EQTMS_PROYECTO_RNA_AND_EPIGENETICS_intervencion/RESULTS/mQTLS_trans_mainOUTPUT_SEX_phv_GOOD_NEW.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"))








############### LANZAR ANALISIS ONTOLOGIA TRANS:






mQTLS_cis <- read.csv2("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.2_ALL_eQTLs_mQTLS_cis_mainOUTPUT_10kb.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) 
#mQTLS_cis <- mQTLS_cis[1:15,]
head(mQTLS_cis)




mQTLS_trans <- read.csv2("/home/mireia.bustos/mQTLs_eQTMS_VASN_multiOmics/PUBMEP/OUTPUTS/2_eQTLs/1.2_ALL_eQTLs_mQTLS_trans_mainOUTPUT_10kb.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA")) 
head(mQTLS_trans)



SYMBOLSTRANS <- matrix(unlist(strsplit(mQTLS_trans[,3], "_")), ncol=2, byrow=TRUE)[,2]



#SYMBOLSTRANS <- matrix(unlist(strsplit(transold[,3], "_")), ncol=2, byrow=TRUE)[,2]


SYMBOLScis <- matrix(unlist(strsplit(mQTLS_cis[,3], "_")), ncol=2, byrow=TRUE)[,2]



SYMBOLSTRANS <- c(SYMBOLSTRANS,SYMBOLScis)
SYMBOLSTRANS <- SYMBOLSTRANS[which(SYMBOLSTRANS != "NA")]


require("hgu133plus2.db")
require("pd.hg.u133.plus.2") # Note: These packages should be changed according to the affymetrix platform under study.
library("plyr")
library("readr")
library("DOSE")

	ae.annots <- AnnotationDbi::select(x= hgu133plus2.db, keys=as.vector(SYMBOLSTRANS),  columns = c("ENTREZID"),  keytype="SYMBOL") 
ae.annots[,2]


library(clusterProfiler)
library(org.Hs.eg.db)

# Crear lista de genes (o cargarla desde un fichero)
gene<-ae.annots[,2]
#Usaremos como referencia unos datos precargados en el paquete DOSE
data(geneList, package="DOSE")

# Ejecutamos análisis de enriquecimiento de términos de Gene Ontology, para Componente Celular
GO <- enrichGO(gene = gene, universe = names(geneList), OrgDb = org.Hs.eg.db,	ont= "CC", pAdjustMethod = "BH", pvalueCutoff = 0.99, qvalueCutoff = 0.99, readable = TRUE)

head(GO)



#"*********************
#"*********************
#"*********************
#"*********************
#write.csv2(head(GO), file="/home/mireia.bustos/discos/sda/augusto/Descargas/RESUMEN_RESULTS_EPIGENETICA_2021/2021/VALIDACION_LISTA_2021/eQTM/LISTA1/ENRICHMENT_TRANSCRIPT_GENES/enrichmentGO_FILTERED_10kb.csv", row.names=TRUE)
#"*********************CAMBIAR en cada caso"
#"*********************
#"*********************
#"*********************


enrichKEGG <- enrichKEGG(gene         = gene,
                 organism     = 'hsa', pvalueCutoff = 0.99, qvalueCutoff = 0.99)
head(enrichKEGG )


#"*********************
#"*********************
#"*********************
#"*********************
#write.csv2(head(enrichKEGG), file="/home/mireia.bustos/discos/sda/augusto/Descargas/RESUMEN_RESULTS_EPIGENETICA_2021/2021/VALIDACION_LISTA_2021/eQTM/LISTA1/ENRICHMENT_TRANSCRIPT_GENES/enrichmentKEGG_FILTERED_10fkb.csv", row.names=TRUE)
#"*********************CAMBIAR en cada caso"
#"*********************
#"*********************
#"*********************











