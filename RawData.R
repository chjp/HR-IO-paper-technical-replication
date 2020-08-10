#  .libPaths(c("/opt/common/CentOS_7/R/R-3.5.1/lib64/R/library", .libPaths()))
library(Seurat)
library(dplyr)
library(Matrix)

rds = readRDS(file = "/data/lab/projects/TCRseq/Input/HRD_seu4_afterTSNE_P80_updated.rds")
counts = rds@raw.data
##### change ensemble ID to gene name #####
ID2name = readRDS("/data/lab/projects/TCRseq/Input/ID2name.rds")
dimnames(counts)[[1]] = ID2name$GeneName
rds = CreateSeuratObject(counts = counts)
rds <- NormalizeData(rds)
rds <- FindVariableFeatures(rds)
#### Batch correction #####
rds@meta.data$Experiment = substr(rownames(rds@meta.data),20,25)
library(SeuratWrappers)
rds <- RunFastMNN(object.list = SplitObject(rds, split.by = "Experiment"))
#saveRDS(rds, file="/data/lab/projects/TCRseq/output/fastMNNcorrect.rds")
#### Batch correction #####

Sdataobj = rds
Sdataobj[["percent.mt"]] <- PercentageFeatureSet(Sdataobj, pattern = "^mt-")
pdf("/data/lab/projects/TCRseq/output/QC.v01.pdf",
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
pdf("/data/lab/projects/TCRseq/output/SigFeatures.v01.pdf", width = 14)
plot1 <- VariableFeaturePlot(Sdataobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

# Scaling the data
all.genes <- rownames(Sdataobj)
Sdataobj <- ScaleData(Sdataobj, features = all.genes)
Sdataobj <- RunPCA(Sdataobj, features = VariableFeatures(object = Sdataobj))


pdf("/data/lab/projects/TCRseq/output/PCgenes.v01.pdf", width = 9)
VizDimLoadings(Sdataobj, dims = 1:2, reduction = "pca")
dev.off()
Sdataobj <- RunUMAP(Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindNeighbors(object = Sdataobj, reduction = "pca", dims = 1:20)
Sdataobj <- FindClusters(object = Sdataobj, resolution = 1)

pdf("/data/lab/projects/TCRseq/output/umap.v01.pdf")
DimPlot(Sdataobj, reduction = "umap", label = TRUE)
dev.off()

Sdataobj <- RunTSNE(object = Sdataobj, dims=1:20)
pdf("/data/lab/projects/TCRseq/output/tsne.v01.pdf")
DimPlot(Sdataobj, reduction = "tsne", label = T)
dev.off()


#lable
pdf("/data/lab/projects/TCRseq/output/MarkerTsne.v02.pdf")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Cd8a", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Cd4", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Pdcd1", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Il7r", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Ccr7", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Cd28", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Mki67", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Ly6c1", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Ly6c2", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Tcf7", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Sell", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Nkg7", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Gzmb", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Isg15", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Cd40lg", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Lag3", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Tnfrsf18", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Icos", min.cutoff = "q9")
    FeaturePlot(Sdataobj, reduction="tsne", features = "Tnfrsf9", min.cutoff = "q9")
dev.off()

markers = c("Cd8a", "Cd4", "Il7r", "Ccr7",
            "Tcf7", "Sell", "Nkg7", "Cd44",
            "Cd14", "Lyz2", "Fcgr3", "Itgax", # Fcgr3 is Cd16; Itgax is Cd11c
            "Ly6c1", "Cd40lg", "Cd79a", "Ncr1",
            "Lag3", "Ms4a1", "Cst3")
            #"Il7r", "Klf6", "Lef1","Ass1",  "Ly6c2","Stat1", "Pdcd1","Ccr7", "Cd28", "Mki67",  "Gzmb", "Isg15", "Cd40lg", "Icos", "Tnfrsf9")

pdf("/data/lab/projects/TCRseq/output/Cluster_markers.v02.pdf")
DotPlot(Sdataobj, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()
#
#pdf("/data/lab/projects/TCRseq/output/Heatmap.pdf")
#DoHeatmap(subset(Sdataobj, downsample = 100), 
#    features = c("Cd4","Cd8a","S100a6","Acp5","Tubb5","Gramd3","St8sia6","Nkg7","Reep5","Crip1","B4galnt1","G0s2","Psen2","Gata3","H3f3b","Cmah","Igfbp4","Chd3"), 
#    size = 3)
#dev.off()

pdf("/data/lab/projects/TCRseq/output/Violin.v01.pdf")
    VlnPlot(Sdataobj, features = c("Cd14"), pt.size = 0.1)
    VlnPlot(Sdataobj, features = c("Cd33"), pt.size = 0.1) # Myeloid cell surface antigen CD33
    VlnPlot(Sdataobj, features = c("Cd3d"), pt.size = 0.1) # arrested T cell differentiation
    VlnPlot(Sdataobj, features = c("Cd79a"), pt.size = 0.1)
    VlnPlot(Sdataobj, features = c("Ncr1"), pt.size = 0.1) # Cd335/ NKp46 NK cell
dev.off()

levels(x = Sdataobj)
Sdataobj2 <- RenameIdents(object = Sdataobj, 
    '0' = 'Cd4+ T cells 1', '1' = 'Cd8a+ T cells 1',
    '2' = 'Cd79a+ B cells 1', '3' = 'Cd4+ T cells 2', '4' = 'Cd4+ T cells 3', 
    '5' = 'Cd14 macrophages 1', '6' = 'Cd14 macrophages 2', '7' = 'Cd79a+ B cells 4', 
    '8' = 'Cd4+ T cells 4', '9' = 'Ncr1+ NK cells 1', '10' = 'Lyz2+ Dendritic', 
    '11' = 'Ncr1+ NK cells 2', '12' = 'Cd8a+ NK cells', '13' = '? Il7r+ cells', 
    '14' = 'Ambiguous1', '15' = 'Cd8a+ T cells 2', '16' = 'Cd79a+ B cells 3', 
    '17' = 'Ambiguous2', '18' = '? Cd8a Cd4 Ly6c1 Dendritic', 
    '19' = 'Ambiguous3', '20' = 'Ambiguous4', '21' = 'Ambiguous5', 
    '22' = '? Cd79a+ Ly6c1+', '23' = 'Cd79a+ B cells 2')
levels(x = Sdataobj2)
pdf("/data/lab/projects/TCRseq/output/tsne.lable.v00.pdf", width = 13)
DimPlot(Sdataobj2, reduction = "tsne", label = T)
dev.off()

#### Compare clusters Fig5 D #####
Cd4Tcell1 <- subset(Sdataobj2, idents = "Cd4+ T cells 1")
Cd4Tcell2 <- subset(Sdataobj2, idents = "Cd4+ T cells 2")
Cd4Tcell3 <- subset(Sdataobj2, idents = "Cd4+ T cells 3")
Cd4Tcell4 <- subset(Sdataobj2, idents = "Cd4+ T cells 4")

avg.c1 <- log1p(AverageExpression(Cd4Tcell1, verbose = FALSE)$RNA)
avg.c2 <- log1p(AverageExpression(Cd4Tcell2, verbose = FALSE)$RNA)
avg.c3 <- log1p(AverageExpression(Cd4Tcell3, verbose = FALSE)$RNA)
avg.c4 <- log1p(AverageExpression(Cd4Tcell4, verbose = FALSE)$RNA)
clusters = data.frame(gene=rownames(avg.c1), c1=avg.c1, c3=avg.c2, c4=avg.c3, c5=avg.c4)
colnames(clusters) = c("gene", "c1", "c2", "c3", "c4")

#SuppTable26 = read.csv("/data/lab/projects/TCRseq/Input/Supplementary_Table_26.csv", sep=",")
genes_focus = c("Macf1", "Cd55", "Fam78a", "Emp3", "S1pr1", "Jun", "Selenop", "Il6ra", "Ppp1r15a", 
                "Tspan32", "St8sia6", "Ccr7", "Foxp1", "Actn1", "Tcf7", "Lef1", "Igfbp4", "Nacc2", 
                "Tgfbr3", "Pdlim4", "Itm2a", "Ass1", "Cd40lg", "Cd4", "Zbtb7b", "Klf6", "Isg15", 
                "Ifit3", "Gbp7", "Usp18", "Irf7", "Ifit1", "Irf1", "Gbp4", "Stat1", "Gbp2") # genes from Fig5d
clusters = clusters[clusters$gene %in% genes_focus,]

library(pheatmap)
pdf("/data/lab/projects/TCRseq/output/Fig5d.v00.pdf")
pheatmap(clusters[,c(2,3,4,5)], scale = "column",
         cluster_rows = F, cluster_cols = T,
         show_rownames = T
)
dev.off()

### compare c1 and c4 ###
c1_c4 = FindMarkers(Sdataobj2, ident.1 = "Cd4+ T cells 1", ident.2 = "Cd4+ T cells 4")
head(c1_c4)
## scatter correlation plot ##
genes.to.label = c("Isg15","Ifit3","Gbp7","Usp18","Irf7","Irf1","Gbp4","Gbp2","Stat1")
library(ggplot2)
library(cowplot)
pdf("/data/lab/projects/TCRseq/output/scatter.v00.pdf")
p1 <- ggplot(clusters, aes(c1, c4)) + geom_point() + ggtitle("Comapre Cluster_1 and Cluster_4")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1 
dev.off()
#### Compare clusters Fig5 D #####
