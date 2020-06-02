
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

#lable
pdf("/data/riazlab/projects/TCRseq/output/chirag/CD8a.pdf")
FeaturePlot(rds3, features = "Cd8a", min.cutoff = "q9")
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/chirag/CD4.pdf")
FeaturePlot(rds3, features = "Cd4", min.cutoff = "q9")
dev.off()

markers = c("Cd8a", "Cd4", "Il7r", "Ccr7",
            "Tcf7", "Sell", "Nkg7","Cd44",
             "S100a4", "Cd14","Lyz2",
            "Ly6c1",  
             "Lag3", "Tnfrsf18",
             "Ifit1", 
           "Gzmb", "Ms4a1", "Ms4a7","Cst3")
pdf("/data/riazlab/projects/TCRseq/output/chirag/DotPlot.pdf")
DotPlot(rds3, features = rev(markers), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()