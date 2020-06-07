#  .libPaths(c("/opt/common/CentOS_7/R/R-3.5.1/lib64/R/library", .libPaths()))
library(Seurat)
library(dplyr)
library(Matrix)

rds = readRDS(file = "/data/riazlab/projects/TCRseq/Input/HRD_seu4_afterTSNE_P80_updated.rds")
counts = rds@raw.data
##### change ensemble ID to gene name #####
ID2name = readRDS("/data/riazlab/projects/TCRseq/Input/ID2name.rds")
dimnames(counts)[[1]] = ID2name$GeneName
rds = CreateSeuratObject(counts = counts)

rds <- NormalizeData(rds)

rds <- FindVariableFeatures(rds)
#### Batch correction #####
rds@meta.data$Experiment = substr(rownames(rds@meta.data),20,25)
library(SeuratWrappers)
rds <- RunFastMNN(object.list = SplitObject(rds, split.by = "Experiment"))
saveRDS(rds, file="/data/riazlab/projects/TCRseq/output/fastMNNcorrect.rds")
#### Batch correction #####

Sdataobj = rds
Sdataobj[["percent.mt"]] <- PercentageFeatureSet(Sdataobj, pattern = "^mt-")
pdf("/data/riazlab/projects/TCRseq/output/QC.v01.pdf",
    width = 21,
    height = 7)
VlnPlot(Sdataobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Sdataobj <- subset(Sdataobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & percent.mt < 5)

#We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset
#(i.e, they are highly expressed in some cells, and lowly expressed in others).
#We and others have found that focusing on these genes in downstream analysis
#helps to highlight biological signal in single-cell datasets.
Sdataobj <- FindVariableFeatures(Sdataobj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Sdataobj), 10)
# plot variable features with and without labels
pdf("/data/riazlab/projects/TCRseq/output/SigFeatures.v01.pdf", width = 14)
plot1 <- VariableFeaturePlot(Sdataobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot2
dev.off()

# Scaling the data
all.genes <- rownames(Sdataobj)
Sdataobj <- ScaleData(Sdataobj, features = all.genes)
Sdataobj <- RunPCA(Sdataobj, features = VariableFeatures(object = Sdataobj))


pdf("/data/riazlab/projects/TCRseq/output/PCgenes.v01.pdf", width = 9)
VizDimLoadings(Sdataobj, dims = 1:2, reduction = "pca")
dev.off()
Sdataobj <- RunUMAP(Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindNeighbors(object = Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindClusters(object = Sdataobj, resolution = 1)

pdf("/data/riazlab/projects/TCRseq/output/umap.v01.pdf")
DimPlot(Sdataobj, reduction = "umap", label = TRUE)
dev.off()

Sdataobj <- RunTSNE(object = Sdataobj, dims=1:20)
pdf("/data/riazlab/projects/TCRseq/output/tsne.v01.pdf")
DimPlot(Sdataobj, reduction = "tsne", label = T)
dev.off()


#lable
pdf("/data/riazlab/projects/TCRseq/output/MarkerTsne.v01.pdf")
    FeaturePlot(Sdataobj, features = "Cd8a", min.cutoff = "q9")
    FeaturePlot(Sdataobj, features = "Cd4", min.cutoff = "q9")
    FeaturePlot(Sdataobj, features = "Pdcd1", min.cutoff = "q9")
    FeaturePlot(Sdataobj, features = "Il7r", min.cutoff = "q9")
    FeaturePlot(Sdataobj, features = "Ccr7", min.cutoff = "q9")
    FeaturePlot(Sdataobj, features = "Cd28", min.cutoff = "q9")
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

markers = c("Cd8a", "Cd4", "Il7r", "Ccr7",
            "Tcf7", "Sell", "Nkg7", "Cd44",
            "Cd14", "Lyz2", "Fcgr3", "Itgax", # Fcgr3 is Cd16; Itgax is Cd11c
            "Ly6c1", "Cd40lg", "Cd79a", "Ncr1",
            "Lag3", "Ms4a1", "Cst3")
            #"Il7r", "Klf6", "Lef1","Ass1",  "Ly6c2","Stat1", "Pdcd1","Ccr7", "Cd28", "Mki67",  "Gzmb", "Isg15", "Cd40lg", "Icos", "Tnfrsf9")

pdf("/data/riazlab/projects/TCRseq/output/Cluster_markers.v02.pdf")
DotPlot(Sdataobj, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()
#
#pdf("/data/riazlab/projects/TCRseq/output/Heatmap.pdf")
#DoHeatmap(subset(Sdataobj, downsample = 100), 
#    features = c("Cd4","Cd8a","S100a6","Acp5","Tubb5","Gramd3","St8sia6","Nkg7","Reep5","Crip1","B4galnt1","G0s2","Psen2","Gata3","H3f3b","Cmah","Igfbp4","Chd3"), 
#    size = 3)
#dev.off()

pdf("/data/riazlab/projects/TCRseq/output/Violin.v01.pdf")
    VlnPlot(Sdataobj, features = c("Cd14"), pt.size = 0.1)
    VlnPlot(Sdataobj, features = c("Cd33"), pt.size = 0.1) # Myeloid cell surface antigen CD33
    VlnPlot(Sdataobj, features = c("Cd3d"), pt.size = 0.1) # arrested T cell differentiation
    VlnPlot(Sdataobj, features = c("Cd79a"), pt.size = 0.1)
    VlnPlot(Sdataobj, features = c("Ncr1"), pt.size = 0.1) # Cd335/ NKp46 NK cell
dev.off()
