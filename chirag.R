library(Seurat)
library(dplyr)
library(Matrix)


rds = readRDS(file = "/data/riazlab/projects/TCRseq/Input/sc_HRD_with_batch_updated.rds")
#str(rds)
rds3=UpdateSeuratObject(rds)
rm(rds)
rds3[["percent.mt"]] = PercentageFeatureSet(rds3, pattern = "^mt-")

changeID = function(geneID){
    tmp = geneID 
    library("biomaRt")
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    #saveRDS(ensembl, "/data/projects/TCRseq/output/MmusculusENSEMBL.rds")
    ID2name <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                filters = 'ensembl_gene_id',
                values = tmp,
                mart = ensembl)
    colnames(ID2name) = c("ID","GeneName")
    library(plyr)
    tmp = as.data.frame(tmp)
    colnames(tmp) = "ID"
    ID2name = join(tmp, ID2name, by = "ID", type = "left", match = "all")
    ID2name[which(duplicated(ID2name$GeneName)),]
    ID2name[which(is.na(ID2name$GeneName)),2] = ID2name[which(is.na(ID2name$GeneName)),1]
    as.vector(ID2name$GeneName)
}

rds3@assays$RNA@counts@Dimnames[[1]] = changeID(rds3@assays$RNA@counts@Dimnames[[1]])
rds3@assays$RNA@data@Dimnames[[1]] = changeID(rds3@assays$RNA@data@Dimnames[[1]])

id = grep(rds3@assays$RNA@data@Dimnames[[1]], pattern="Cd3")
genes = rds3@assays$RNA@data@Dimnames[[1]]
genes[id]


pdf("/data/riazlab/projects/TCRseq/output/chirag/tsne.v00.pdf")
DimPlot(rds3, reduction = "tsne", label = T)
dev.off()

#lable
pdf("/data/riazlab/projects/TCRseq/output/chirag/CD8a.pdf")
FeaturePlot(rds3, features = "Cd8a", min.cutoff = "q9")
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/chirag/CD4.pdf")
FeaturePlot(rds3, features = "Cd4", min.cutoff = "q9")
dev.off()

markers = c("Cd8a", "Cd4", "Il7r", "Ccr7",
            "Tcf7", "Sell", "Nkg7", "Cd44",
            "Cd14", "Cd3d", "Lyz2", "Fcgr3", "Itgax", # Fcgr3 is Cd16; Itgax is Cd11c
            "Ly6c1", "Cd40lg", "Cd79a", "Ncr1",
            "Lag3", "Ms4a1", "Cst3")
pdf("/data/riazlab/projects/TCRseq/output/chirag/DotPlot.v01.pdf")
DotPlot(rds3, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()

pdf("/data/riazlab/projects/TCRseq/output/chirag/Violin.pdf")
  VlnPlot(rds3, features = c("Cd14"), pt.size = 0.1)
  VlnPlot(rds3, features = c("Cd33"), pt.size = 0.1) # Myeloid cell surface antigen CD33
  VlnPlot(rds3, features = c("Cd3d"), pt.size = 0.1) # arrested T cell differentiation
  VlnPlot(rds3, features = c("Cd79a"), pt.size = 0.1)
  VlnPlot(rds3, features = c("Ncr1"), pt.size = 0.1) # Cd335/ NKp46 NK cell
  VlnPlot(rds3, features = c("Cst3"), pt.size = 0.1)
dev.off()

#### re-cluster myloid cells ####
myloid_cells <- SubsetData(object = rds3, ident.use = c(8,10,11), do.clean = TRUE, do.scale = TRUE)
myloid_cells = CreateSeuratObject(counts=myloid_cells@assays$RNA@counts)
saveRDS(myloid_cells, file="/data/riazlab/projects/TCRseq/output/chirag/MyloidCluster.rds")
myloid_cells=readRDS("/data/riazlab/projects/TCRseq/output/chirag/MyloidCluster.rds")
myloid_cells <- FindVariableFeatures(myloid_cells, selection.method = "vst", nfeatures = 2000)

top100 <- head(VariableFeatures(myloid_cells), 100)
pdf("/data/riazlab/projects/TCRseq/output/chirag/SigFeatures.v00.pdf", width = 14)
  plot1 <- VariableFeaturePlot(myloid_cells)
  plot2 <- LabelPoints(plot = plot1, points = top100, repel = TRUE)
  plot2
dev.off()

all.genes <- rownames(myloid_cells)

myloid_cells <- ScaleData(myloid_cells, features = all.genes)
myloid_cells <- RunPCA(myloid_cells, features = VariableFeatures(object = myloid_cells))
pdf("/data/riazlab/projects/TCRseq/output/chirag/PCgenes.v00.pdf", width = 9)
VizDimLoadings(myloid_cells, dims = 1:2, reduction = "pca")
dev.off()
myloid_cells <- FindNeighbors(object = myloid_cells, reduction = "pca", dims = 1:10)
myloid_cells <- FindClusters(object = myloid_cells, resolution = 0.8)


myloid_cells <- RunTSNE(object = myloid_cells, dims=1:10)
pdf("/data/riazlab/projects/TCRseq/output/chirag/tsne.v02.pdf")
DimPlot(myloid_cells, reduction = "tsne", label = T)
dev.off()

pdf("/data/riazlab/projects/TCRseq/output/chirag/Violin.v02.pdf")
  VlnPlot(myloid_cells, features = c("Cd14"), pt.size = 0.1)
  VlnPlot(myloid_cells, features = c("Cd160"), pt.size = 0.1) # see all.genes[grep(all.genes, pattern="Cd16")]
  VlnPlot(myloid_cells, features = c("Fcgr3"), pt.size = 0.1) # Synonymous of Cd16

  VlnPlot(myloid_cells, features = c("Csf1r"), pt.size = 0.1)
  VlnPlot(myloid_cells, features = c("Flt3"), pt.size = 0.1)
dev.off()

markers = c("Msr1", "Mrc1", "Arg1", "Nos2", # Cd206=Mrc1 Cd204=Msr1
            "Cd274", "Vegfa", "Apoe", "Trem2", #Cd274=Pdl1 
            "Mmp8", "Itgae", "Itga1", "Il1b",
            "Ifitm1", "Ifitm3", "Isg15", "Clec4e",
            "Stat1", "Txnip", "Irf8", "Ccl3")
pdf("/data/riazlab/projects/TCRseq/output/chirag/DotPlot.v03.pdf")
DotPlot(rds3, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()
#### re-cluster myloid cells ####