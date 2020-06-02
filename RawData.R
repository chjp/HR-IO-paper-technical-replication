library(Seurat)
library(dplyr)
library(Matrix)

rds = readRDS(file = "/data/riazlab/projects/TCRseq/HRD_seu4_afterTSNE_P80_updated.rds")

class(rds@raw.data)
counts = rds@raw.data

geneID = dimnames(counts)[[1]]
geneID = as.data.frame(geneID)
which(duplicated(geneID))
colnames(geneID)="ID"
library("biomaRt")
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#saveRDS(ensembl, "/data/projects/TCRseq/output/MmusculusENSEMBL.rds")
ID2name <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
             filters = 'ensembl_gene_id',
             values = geneID,
             mart = ensembl)
which(duplicated(ID2name$external_gene_name))

#ID2name = merge(geneID, ID2name, by.x = "ID", by.y = "ensembl_gene_id", all.x = T)
#ID2name[which(is.na(ID2name[,2])),2] = ID2name[which(is.na(ID2name[,2])),1] 
colnames(ID2name) = c("ID","GeneName")
library(plyr)
ID2name = join(geneID, ID2name, by = "ID", type = "left", match = "all")
which(duplicated(ID2name$GeneName))
ID2name[which(duplicated(ID2name$GeneName)),]
ID2name[which(is.na(ID2name$GeneName)),2] = ID2name[which(is.na(ID2name$GeneName)),1]
which(duplicated(ID2name$GeneName))
dimnames(counts)[[1]] = ID2name$GeneName

Sdataobj = CreateSeuratObject(counts = counts)
object.size(Sdataobj)
Sdataobj[["percent.mt"]] <- PercentageFeatureSet(Sdataobj, pattern = "^mt-")
pdf("/data/projects/TCRseq/output/QC.pdf",
    width = 21,
    height = 7)
VlnPlot(Sdataobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Sdataobj <- subset(Sdataobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & percent.mt < 5)
Sdataobj <- NormalizeData(Sdataobj, normalization.method = "LogNormalize", scale.factor = 10000)


#We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset
#(i.e, they are highly expressed in some cells, and lowly expressed in others).
#We and others have found that focusing on these genes in downstream analysis
#helps to highlight biological signal in single-cell datasets.
Sdataobj <- FindVariableFeatures(Sdataobj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Sdataobj), 10)
# plot variable features with and without labels
pdf("/data/riazlab/projects/TCRseq/output/SigFeatures.pdf", width = 14)
plot1 <- VariableFeaturePlot(Sdataobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# Scaling the data
all.genes <- rownames(Sdataobj)
Sdataobj <- ScaleData(Sdataobj, features = all.genes)
Sdataobj <- RunPCA(Sdataobj, features = VariableFeatures(object = Sdataobj))
#saveRDS(Sdataobj, file="/data/riazlab/projects/TCRseq/output/tmp.rds")
Sdataobj = readRDS(file="/data/riazlab/projects/TCRseq/output/tmp.rds")


#pdf("/data/riazlab/projects/TCRseq/output/PCgenes.pdf", width = 9)
#VizDimLoadings(Sdataobj, dims = 1:2, reduction = "pca")
#dev.off()
Sdataobj <- RunUMAP(Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindNeighbors(object = Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindClusters(object = Sdataobj, resolution = 1)
#saveRDS(Sdataobj, file="/data/riazlab/projects/TCRseq/output/fastMNN.Sdataobj.rds")
Sdataobj = readRDS(file="/data/riazlab/projects/TCRseq/output/fastMNN.Sdataobj.rds")
which(Sdataobj@assays$RNA@counts@x < 0) # none
pdf("/data/riazlab/projects/TCRseq/output/umap.pdf")
DimPlot(Sdataobj, reduction = "umap", label = TRUE)
dev.off()

Sdataobj <- RunTSNE(object = Sdataobj, dims=1:20)
pdf("/data/riazlab/projects/TCRseq/output/tsne.pdf")
DimPlot(Sdataobj, reduction = "tsne", label = T)
dev.off()


#lable
pdf("/data/riazlab/projects/TCRseq/output/Cd8a.pdf")
FeaturePlot(Sdataobj, features = "Cd8a", min.cutoff = "q9")
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/Cd4.pdf")
FeaturePlot(Sdataobj, features = "Cd4", min.cutoff = "q9")
dev.off()
#pdf("/data/riazlab/projects/TCRseq/output/MarkerGenes.pdf")
#FeaturePlot(Sdataobj, features = c("Cd8a","Cd4","Pdcd1","Pdcd1","Il7r","Ccr7", 
#                                "Cd28",  "Mki67",  "Ly6c1",  "Ly6c2",  "Tcf7",  "Sell",
#                                "Nkg7",  "Gzmb",  "Isg15",  "Cd40lg",  "Lag3",  
#                                "Tnfrsf18",  "Icos",  "Tnfrsf9"),min.cutoff = "q9")
pdf("/data/riazlab/projects/TCRseq/output/tmp.pdf")
FeaturePlot(Sdataobj, features = "Cd8a", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Cd4", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Pdcd1", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Il7r", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Ccr7", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Cd28", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Cd8a", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Mki67", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Ly6c1", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Ly6c2", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Tcf7", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Sell", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Nkg7", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Gzmb", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Isg15", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Cd40lg", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Lag3", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Tnfrsf18", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Icos", min.cutoff = "q9")
FeaturePlot(Sdataobj, features = "Tnfrsf9", min.cutoff = "q9")
dev.off()

Idents(Sdataobj)

factor(Idents(Sdataobj))
markers = c("Cd8a", "Cd4", "Il7r", "Ccr7",
            "Tcf7", "Sell", "Nkg7","Cd44",
             "S100a4", "Cd14","Lyz2",
            "Ly6c1",  
             "Lag3", "Tnfrsf18",
             "Ifit1", 
           "Gzmb", "Ms4a1", "Ms4a7","Cst3")
            #"Il7r", "Klf6", "Lef1","Ass1",  "Ly6c2","Stat1", "Pdcd1","Ccr7", "Cd28", "Mki67",  "Gzmb", "Isg15", "Cd40lg", "Icos", "Tnfrsf9")

pdf("/data/riazlab/projects/TCRseq/output/Cluster_markers.pdf")
DotPlot(Sdataobj, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()

pdf("/data/riazlab/projects/TCRseq/output/Heatmap.pdf")
DoHeatmap(subset(Sdataobj, downsample = 100), 
    features = c("Cd4","Cd8a","S100a6","Acp5","Tubb5","Gramd3","St8sia6","Nkg7","Reep5","Crip1","B4galnt1","G0s2","Psen2","Gata3","H3f3b","Cmah","Igfbp4","Chd3"), 
    size = 3)
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/Violin.pdf")
VlnPlot(Sdataobj, features = c("Cd14"), pt.size = 0.1)
VlnPlot(Sdataobj, features = c("Cd33"), pt.size = 0.1) # Myeloid cell surface antigen CD33
VlnPlot(Sdataobj, features = c("Cd3d"), pt.size = 0.1) # arrested T cell differentiation
VlnPlot(Sdataobj, features = c("Cd79a"), pt.size = 0.1)
VlnPlot(Sdataobj, features = c("Ncr1"), pt.size = 0.1) # Cd335/ NKp46 NK cell
dev.off()

geneName = dimnames(Sdataobj@assays$RNA)[[1]]
ind = grep(dimnames(Sdataobj@assays$RNA)[[1]], pattern = "fce", ignore.case = T)
geneName[ind]