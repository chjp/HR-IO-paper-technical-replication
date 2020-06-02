
library(Seurat)
library(dplyr)
library(Matrix)
#install.packages("batchelor")

rds = readRDS(file = "/data/riazlab/projects/TCRseq/Input/sc_HRD_with_batch_updated.rds")
#str(rds)
rds3=UpdateSeuratObject(rds)
rds3[["percent.mt"]] = PercentageFeatureSet(rds3, pattern = "^mt-")
counts = rds3@assays$RNA@counts
#rm(rds,rds3)

##### change gene ID to gene name. #####
geneID = dimnames(counts)[[1]]
geneID = as.data.frame(geneID)
colnames(geneID)="ID"
library("biomaRt")
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ID2name <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
             filters = 'ensembl_gene_id',
             values = geneID,
             mart = ensembl)

#### failed solution ####
#ID2name = merge(geneID, ID2name, by.x = "ID", by.y = "ensembl_gene_id", all.x = T)
#ID2name[which(is.na(ID2name[,2])),2] = ID2name[which(is.na(ID2name[,2])),1] 
### it shuffled gene order
#### failed solution ####

#### failed solution ####
#library(hash)
#hID2name = hash(ID2name$ID, ID2name$external_gene_name)
#hID2name[[ "ENSMUSG00000095041"]]
#v = hID2name[[ dimnames(counts)[[1]] ]]
#t = as.data.frame(values(v))
#values(v, keys = "ENSMUSG00000095041")
# Error in get(k, x) : object 'ENSMUSG00000097347' not found
#### failed solution ####

colnames(ID2name) = c("ID","GeneName")
library(plyr)
ID2name = join(geneID, ID2name, by = "ID", type = "left", match = "all")
which(duplicated(ID2name$GeneName))
ID2name[which(duplicated(ID2name$GeneName)),]
ID2name[which(is.na(ID2name$GeneName)),2] = ID2name[which(is.na(ID2name$GeneName)),1]
which(duplicated(ID2name$GeneName))
dimnames(counts)[[1]] = ID2name$GeneName
##### change gene ID to gene name. #####




#rds3norml <- NormalizeData(rds3, normalization.method = "LogNormalize", scale.factor = 10000)
#counts1 = rds3norml@assays$RNA@counts
barcodes = dimnames(counts)[[2]]
unique(substr(barcodes,20,25)) # BRCA1 mut and BRCA2 mut: 17QQ and 21EE
C17QQ2 = counts[,grepl("17QQ2", barcodes)]
C17QQ3 = counts[,grepl("17QQ3", barcodes)]

C21E1 = counts[,grepl("21E1", barcodes)]
C21E2 = counts[,grepl("21E2", barcodes)]
C21E3 = counts[,grepl("21E3", barcodes)]

CPAR1 = counts[,grepl("PAR1", barcodes)]
CPAR2 = counts[,grepl("PAR2", barcodes)]
CPAR3 = counts[,grepl("PAR3", barcodes)]
str(CPAR3)

library(batchelor)
#out = fastMNN(C17QQ2,C17QQ3,C21E1,C21E2,C21E3,CPAR1,CPAR2,CPAR3)
#saveRDS(object=out, file="/data/riazlab/projects/TCRseq/output/MNNcorrect.rds")
#out=readRDS("/data/riazlab/projects/TCRseq/output/MNNcorrect.rds")
#RunPCA(object=out)
rm(rds3,rds,counts)
#summary(out)
#runPCA(out)
#class(out@assays)
#str(out@assays)
#dim(out@assays)
#colnames(out@assays@data@listData)
#str(out@assays@data@listData$reconstructed)
#getListElement(out@assays[3,3],2)
#names(out@assays)
#dimnames(out@assays@data@ListData)


Sdataobj <- CreateSeuratObject(counts = out@assays@data@listData$reconstructed)
saveRDS(Sdataobj, file="/data/riazlab/projects/TCRseq/output/tmp.rds")
#Sdataobj = readRDS(file="/data/riazlab/projects/TCRseq/output/tmp.rds")
Sdataobj <- FindVariableFeatures(Sdataobj, selection.method = "vst", nfeatures = 2000)

#https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/RunTSNE
all.genes <- rownames(Sdataobj)
Sdataobj <- ScaleData(Sdataobj, features = all.genes)
Sdataobj <- RunPCA(Sdataobj, features = VariableFeatures(object = Sdataobj))
#Sdataobj <- JackStraw(Sdataobj, num.replicate = 100)
#Sdataobj <- ScoreJackStraw(Sdataobj, dims = 1:20)
#pdf("/data/riazlab/projects/TCRseq/output/JackStraw.fastMNN.pdf")
#JackStrawPlot(Sdataobj, dims = 1:20)
#dev.off()


Sdataobj <- FindNeighbors(Sdataobj, dims = 1:20)
Sdataobj <- FindClusters(Sdataobj, resolution = 1)
Sdataobj = RunTSNE(Sdataobj, dims = 1:20)
pdf("/data/riazlab/projects/TCRseq/output/tsne.fastMNN.pdf")
DimPlot(Sdataobj, reduction = "tsne", label = T)
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/pca.fastMNN.pdf")
DimPlot(Sdataobj, reduction = "pca")
dev.off()


pdf("/data/riazlab/projects/TCRseq/output/CD8a.BC.pdf")
FeaturePlot(Sdataobj, features = "ENSMUSG00000053977", min.cutoff = "q9")
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/CD4.BC.pdf")
FeaturePlot(Sdataobj, features = "ENSMUSG00000023274", min.cutoff = "q9")
dev.off()

features = c("ENSMUSG00000041959","ENSMUSG00000050592")
pdf("/data/riazlab/projects/TCRseq/output/heatmap.pdf")
DoHeatmap(subset(Sdataobj, downsample = 10), features = features, size = 3)
dev.off()