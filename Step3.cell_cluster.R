#####3.single cell intergration and cell type annotation####
rm(list=ls())
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SingleR)
  library(scater)
  library(scran)
  library(dplyr)
  library(SingleCellExperiment)
  library(celldex)
  library(BiocParallel)
  library(scRNAseq)
})

#PAT=''#PATH of pipeline directory

PATH=paste0(PAT,'./data/')

dir = c(paste0(PATH, "GSE67835_count.csv"), 
        paste0(PATH, "GSE71315_count.csv"), 
        paste0(PATH, "GSE125288_count.csv"),
        paste0(PATH, "GSE75140_count.csv"))
names(dir) = c('GSE67835', 'GSE71315','GSE75140','GSE125288')

for(i in 1:length(dir)){
    counts <- read.table(dir[i],sep=',',header=TRUE,row.names = 1)
    tmp <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))
    tmp.data<-calculateQCMetrics(tmp)
    low_lib_tmp.data <- isOutlier(tmp.data$log10_total_counts, type="lower", nmad=3)
    low_genes_tmp.data <- isOutlier(tmp.data$log10_total_features_by_counts, type="lower", nmad=3)
    tmp.data <- tmp.data[, !(low_lib_tmp.data | low_genes_tmp.data)]
    tmp.data
    write.csv(tmp.data, file=dir[i])
}

scRNAlist <- list()
minGene=100
maxGene=20000
for(i in 1:length(dir)){
  counts <- read.table(dir[i],sep=',',header=TRUE,row.names = 1)
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=1, project  = names(dir[i]))
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-")
  p<-VlnPlot(scRNAlist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(scRNAlist[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scRNAlist[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p+plot1 + plot2)#set threshold
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(scRNAlist[[i]]@assays$RNA))
  HB.genes <- rownames(scRNAlist[[i]]@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes)
  col.num <- length(levels(as.factor(scRNAlist[[i]]@meta.data$orig.ident)))
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset=nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < 10)
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures=2000)
}
scRNAlist

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
scRNA3<- IntegrateData(anchorset = scRNA.anchors)
scRNA<-scRNA3 ####Intergrated data

inte <- SingleCellExperiment(assays = list(counts = as.matrix(scRNA@assays$integrated@data)))
factor<-computeSumFactors(inte)
summary(sizeFactors(factor))
factor1 <-as.matrix(factor@colData)
write.csv(factor1, file=paste0(PATH, "human_brain_factor.csv"))

scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA, ndims=50, reduction="pca")
plot1+plot2

pc.num=1:10#set threshold
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution =2)

metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
table(scRNA@meta.data$seurat_clusters)
metadata$Cluster_id<-metadata

scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
plot1 = DimPlot(scRNA, reduction = "tsne", label=T,pt.size = 1.5)
plot1

scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap') 
plot3 = DimPlot(scRNA, reduction = "umap", label=T)
plot3
plot3+plot1

cell<-list(
    Astrocyte<- c("GFAP","SLC1A2","AQP4","GJA1","GJB6","SLC4A4","SLC39A12","FABP7")
    Endothelial<- c("CD34","VWF","FTL1","RGS5","PTPRB","PALMD","APOLD1","NOSTRIN","ESAM","CLDN5","PECAM1")
    Microglia<- c("ITGAM","CX3CR1","CCL3","CSF1R","P2RY12","CCL4","IBA1","C1QA","LAPTM5","C1QB","C1QC","AIF1","HLA-DRA")
    Oligodendrocyte<- c("MOG","MBP","MAG","OPALIN","MOBP","CNP","CLDN11","PLP1","O4")
    OPC<- c("OLIG1","OLIG2","PDGFRA","VCAN","SOX10","GPR17")
    Excitatory<- c("NEUROD6","STMN2","SYT1")
    Inhibitory<- c("CALB2","GAD1","GAD2")
    Cajal_Retzuis<- c("RELN")
    Purkinje_cell<- c("PCP4")
    NPC<- c("PAX6","SOX2")
    VLMC<- c("OGN","LUM","DCN")
    Tancyte<- c("RAX","CX43")
    Schwann<- c("SOX10","FOXDD3")
    Ependymocytes<- c("CCDC153","AK7","PITO","DYNLRB2","CALML4","RBP1")
    Mural<- c("PDGFRB","RGS5")
    Neutrophils<- c("S100A9")
    Macrophages<- c("APOE","MS4A7","MS4A6C","LYZ2","TGFBI")
    NKT_cell<- c("PTPRC","CD3E","ITGAL")
)
for (i in 1:length(cell)) {
    ele<-cell[[i]]
    p3 <- VlnPlot(scRNA, features=ele, slot="counts", log=TRUE,pt.size=0, group.by="seurat_clusters", ncol=2)
    p3
}
ref <- HumanPrimaryCellAtlasData()
pred.scRNA1 <- SingleR(test = scRNA@assays$RNA@counts, ref = ref,labels = ref$label.main, clusters = scRNA@active.ident, fine.tune = TRUE, BPPARAM = MulticoreParam(40))
pred.scRNA1$pruned.labels
plotScoreHeatmap(pred.scRNA1)

cg <- LaMannoBrainData('human-es')
pred.scRNA <- SingleR(test = scRNA@assays$RNA@counts, ref = cg,labels = cg$label.main, clusters = scRNA@active.ident, fine.tune = TRUE, BPPARAM = MulticoreParam(40))
pred.scRNA$pruned.labels
plotScoreHeatmap(pred.scRNA)

celltype1 = data.frame(ClusterID=rownames(pred.scRNA),
                       celltype=pred.scRNA$pruned.labels,
                       celltype1=pred.scRNA1$pruned.labels,
                       stringsAsFactors = F)
celltype1

scRNA1<-scRNA
new.cluster.ids1<- c("Glutamatergic_excitatory", "GABAergicneuron","GABAergicneuron","NPC","Oligodendrocyte",
                     "Glutamatergic_excitatory","Glutamatergic_excitatory","Glutamatergic_excitatory","Astrocyte",
                     "OPC",
                     "proliferating","proliferating","Microglia","proliferating",
                     "GABAergicneuron","proliferating","intermediate progenitor","Glutamatergic_excitatory","proliferating",
                     "Astrocyte","Astrocyte","Endothelial")
names(new.cluster.ids1) <- levels(scRNA1)
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids1)

ploty1 <-DimPlot(scRNA1, reduction = "umap", label = FALSE, pt.size = 1.0,cols=c('#d9ef8b','#dfc27d','#f46d43',
                                                                                 '#92c5de','#c51b7d','#000000',
                                                                                 '#b2abd2','#ffeda0','#e0f3f8'))
ploty2 <-DimPlot(scRNA1, reduction = "tsne", label = FALSE, pt.size = 1.0,cols=c('#d9ef8b','#dfc27d','#f46d43',
                                                                                 '#92c5de','#c51b7d','#000000',
                                                                                 '#b2abd2','#ffeda0','#e0f3f8'))
ploty1+plot3+ploty2+plot1

metadata <-  FetchData(object = scRNA1, vars = c('ident',"tSNE_1", "tSNE_2", "UMAP_1","UMAP_2","nCount_RNA","seurat_clusters"))
write.csv(metadata, file=paste0(PATH, "brain_metadata.csv"))
combined =FindAllMarkers(scRNA, test="wilcox",min.pct = 0.25,only.pos = TRUE,  thresh.use = 0.25)
write.csv(combined.filter, file=paste0(PATH, "brain_DE.csv"))
combined.filter=combined %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC )
write.csv(combined.filter, file=paste0(PATH, "brain_marker.csv"))
saveRDS(scRNA, paste0(PATH, "brain_cluster.RDS"))




