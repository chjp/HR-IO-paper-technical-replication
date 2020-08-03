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


#### pheatmap ####
library(ggplot2)
library(cowplot)

c1 <- subset(rds3, idents = "1")
c3 <- subset(rds3, idents = "3")
c4 <- subset(rds3, idents = "4")
c5 <- subset(rds3, idents = "5")

avg.c1 <- log1p(AverageExpression(c1, verbose = FALSE)$RNA)
avg.c3 <- log1p(AverageExpression(c3, verbose = FALSE)$RNA)
avg.c4 <- log1p(AverageExpression(c4, verbose = FALSE)$RNA)
avg.c5 <- log1p(AverageExpression(c5, verbose = FALSE)$RNA)
clusters = data.frame(gene=rownames(avg.c1), c1=avg.c1, c3=avg.c3, c4=avg.c4, c5=avg.c5)
colnames(clusters) = c("gene", "c1", "c3", "c4", "c5")

SuppTable26 = read.csv("/data/riazlab/projects/TCRseq/Input/Supplementary_Table_26.csv", sep=",")
#genes_focus = unique(as.character(SuppTable26$V2))
genes_focus = c("Isg15","Ifit3","Gbp7","Usp18","Irf7","Irf1","Gbp4","Gbp2","Stat1",
                "Ccr9","Macf1"
)

clusters = clusters[clusters$gene %in% genes_focus,]
library(pheatmap)
pdf("/data/riazlab/projects/TCRseq/output/chirag/Fig5d.v01.pdf")
pheatmap(clusters[,c(2,3,4,5)], scale = "column",
         cluster_rows = T, cluster_cols = T,
         show_rownames = T
)
dev.off()
pdf("/data/riazlab/projects/TCRseq/output/chirag/Fig5d.v02.pdf")
pheatmap(clusters[,c(2,3,4,5)], scale = "row",
         cluster_rows = T, cluster_cols = T,
         show_rownames = T
)
dev.off()

### compare c1 and c4 ###
c1_c4 = FindMarkers(rds3, ident.1 = "1", ident.2 = "4")
head(c1_c4)
## scatter correlation plot ##
genes.to.label = c("Isg15","Ifit3","Gbp7","Usp18","Irf7","Irf1","Gbp4","Gbp2","Stat1")

install.packages("ggplot2",dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(ggplot2)
library(cowplot)
saveRDS(clusters,"~/tmp.rds")
clusters = readRDS("~/tmp.rds")
pdf("/data/riazlab/projects/TCRseq/output/chirag/scatter.v00.pdf")
p1 <- ggplot(clusters, aes(c1, c4)) + geom_point() + ggtitle("Comapre Cluster_1 and Cluster_4")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1 
dev.off()