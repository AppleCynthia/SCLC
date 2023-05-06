

library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(snow)

source('/src/utilities.R')

dir_path <- '/data/part3'
setwd(dir_path)


celltype <- 'Endo'
merged.se.obj <- readRDS(file.path('/data/part1', '1.merged.se.obj.RDS'))

obj <- subset(merged.se.obj, cells = 
                rownames(merged.se.obj@meta.data
                         [merged.se.obj@meta.data$all_CellType.1 %in% celltype &
                             merged.se.obj@meta.data$condition == 'Tumor', 
                         ]))

#run seurat
run_seurat <- function(se.obj){
	se.obj<- NormalizeData(se.obj, normalization.method = "LogNormalize", 
		scale.factor = 10000)
	se.obj <- FindVariableFeatures(se.obj, selection.method = "vst", nfeatures = 2000)
	se.obj <- ScaleData(se.obj)
	se.obj <- RunPCA(se.obj, features = VariableFeatures(object = se.obj))
	se.obj <- FindNeighbors(se.obj, dims = 1:15)
	se.obj <- FindClusters(se.obj, resolution = 0.5)
	se.obj <- RunUMAP(se.obj, dims = 1:15)
	se.obj <- RunTSNE(se.obj, dims = 1:15)
	return(se.obj)
}

obj <- run_seurat(obj)

pdf(file.path(dir_path,   'Endo',paste0(celltype, '.seurat_clusters.umap.pdf')),7.5,6)
DimPlot(obj, label =T, group.by ='seurat_clusters', label.size =5) 
dev.off()


use.markers <- c('KCNE3', 'DLL4', 'FCN3', 'HPGD', 'PGF',#tip 2,3,4
                 'MKI67','RRM2', #tumor 8
                 'ACKR1', 'SELP', 'C7', #STALK 1,5,6,7, (8,9))
                 'PROX1', 'CCL21', 'TFF3',#LYMP  0              
                 'CD48','CD7','CXCR4')#EPCs 9 

VlnPlot(obj, features = use.markers)
ggsave(file = file.path(dir_path,  'Endo', paste0(celltype, '.sub.markers.pdf')), 
  width = 12, height =10)
FeaturePlot(obj, features = use.markers)
ggsave(file = file.path(dir_path,  'Endo', paste0(celltype, '.sub.markers.featureplot.pdf')), 
  width = 12, height =10)


##############       annotate sub celltye
cellanno <- obj@meta.data 
cellanno$cell.name <- rownames(cellanno)

cellanno <- cellanno %>% mutate(subtypes = 
                                  ifelse(RNA_snn_res.0.5 %in% c(1,5,6,7), 'Stalk-like ECs', 
                                         ifelse(RNA_snn_res.0.5 %in% c(2,3,4), 'Tip-like ECs',
                                                ifelse(RNA_snn_res.0.5 ==0, 'Lymphatic ECs',
                                                       ifelse(RNA_snn_res.0.5 ==8, 'Tumor ECs', 
                                                              'EPCs')))))


identical(cellanno$cell.name, rownames(obj@meta.data))
obj@meta.data$subtypes <- cellanno$subtypes
saveRDS(obj,'Endo.rds')


endo_celltype_colors <- c('#DDA0DD','#D2B4BC','#FFC0CB','#B0C4DE','#EEE8AA')
names(endo_celltype_colors) <- c('EPCs',  'Lymphatic ECs', 'Stalk-like ECs', 'Tip-like ECs', 'Tumor ECs' )


identical(obj@meta.data$seurat_clusters, obj@meta.data$RNA_snn_res.0.5)



pdf(file.path(dir_path,   'Endo',paste0(celltype, '.subtypes.umap.pdf')),7.5,6)
DimPlot(obj, group.by = 'subtypes', label =T, label.size =5,cols = endo_celltype_colors) 
dev.off()

pdf(file.path(dir_path,   'Endo',paste0(celltype, '.P_ID.umap.pdf')),7.5,6)
DimPlot(obj, group.by ='Patient_ID', cols = patient_colors)
dev.off()

pdf(file.path(dir_path,   'Endo',paste0(celltype, '.Treatment.umap.pdf')),7.5,6)
DimPlot(obj, group.by = 'Treatment', cols =treatment_colors)
dev.off()


Idents(obj) <- 'subtypes'
avg.nec <- log1p(AverageExpression(obj, verbose = FALSE)$RNA)

pdf(file.path(dir_path,   'Endo',
 paste0('sub.marker.pdf')), width =4, height =6)
pheatmap(avg.nec[use.markers, ],
color = colorRampPalette(c("#004e66", "#e1eef6", "#ff5f2e"))(200),
scale = 'row', cluster_cols =F , cluster_rows =F)

dev.off()


##########   DGE in subtype

table(obj@meta.data$Treatment)

use.sub <- 'Stalk-like ECs'
use.sub <- 'Tip-like ECs'
use.sub <- 'Lymphatic ECs'



Idents(obj) <- 'subtypes'
DE_genes <- FindMarkers(obj, 
	ident.1 = "Chemotherapy", 
	ident.2 = "Treatment_naive",
	group = 'Treatment',
	subset.ident = use.sub,
	 verbose = FALSE,
	 logfc.threshold=0.25,
   min.diff.pct = 0.1)


write.csv(DE_genes, file = file.path(dir_path,   'Endo',
	paste0(use.sub, '.DE.genes.csv')))
dim(DE_genes)
dim(DE_genes[DE_genes$avg_log2FC > 0,])
dim(DE_genes[DE_genes$avg_log2FC < 0,])


subobj <- subset(obj,  cells = rownames(obj@meta.data[obj@meta.data$subtypes == use.sub, 
]))
Idents(subobj) <- 'Treatment'
avg.nec <- log1p(AverageExpression(subobj
	, verbose = FALSE)$RNA)

table(subobj@meta.data$Treatment)

pdf(file.path(dir_path,   'Endo',
 paste0(use.sub, '.DE.genes.heatmap.pdf')), width =4, height =8)
colors = colorRampPalette(c("#3b8686", "#79bd9a", "#cff09e", "#f9cdad","#fec9c9","#fe4365"))(200)
df <- avg.nec[head(rownames(DE_genes),100),]
df <- df[,c(2,1)]
pheatmap(df[,], color =colors,
fontsize =6, cluster_col =F)
dev.off()



#################################  GO and KEGG
#GO
GO_genes <- rownames(DE_genes)

ensembl.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
	keytype = "SYMBOL", column="ENSEMBL")

ensembl.genes <- na.omit(ensembl.genes)

GO<- enrichGO(gene = ensembl.genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = 'ALL', pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
GO <- setReadable(GO, 'org.Hs.eg.db', 'ENSEMBL')

write.csv(GO, file = file.path(dir_path,   'Endo',
	paste0(use.sub, '.GO.csv')))

pdf(file.path(dir_path,  'Endo',
 paste0(use.sub, '.dot.GO.pdf')),10,6)
dotplot(GO, showCategory = 25)
dev.off()

pdf(file.path( dir_path,  'Endo',
               paste0(use.sub, '.bar.GO.bar.pdf')),10,6)
barplot(GO, showCategory = 25)
dev.off()


#KEGG
entrez.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
	keytype = "SYMBOL", column="ENTREZID")
kk <- enrichKEGG(gene         = entrez.genes,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
write.csv(edox, file = file.path(dir_path,   'Endo',
 paste0(use.sub, '.kegg.csv')))

pdf(file = file.path(dir_path,   'Endo',
	paste0(use.sub, '.cnet.KEGG.pdf')))
cnetplot(edox, showCategory = 15)
dev.off()

pdf(file = file.path(dir_path,   'Endo',
                     paste0(use.sub, '.dot.KEGG.pdf')))
dotplot(edox, showCategory = 15)
dev.off()

pdf(file = file.path(dir_path,   'Endo',
                     paste0(use.sub, '.bar.KEGG.pdf')),10,6)
barplot(edox, showCategory = 15)
dev.off()






####  Stalk-like ECs
use.sub <- 'Stalk-like ECs'
use.sub <- 'Tip-like ECs'
use.sub <- 'Lymphatic ECs'

celltype <- "Stalk-like ECs"

stalk.obj <- subset(obj, cells = 
                rownames(obj@meta.data
                         [obj@meta.data$subtypes == "Stalk-like ECs" ,
                            
                         ]))

############ angiogenesis
gene_vln <- c("ID1","CTSH","CLDN5","ADAMTS1","KLF2","PTGIS","GATA2","TNFRSF12A")

gene_vln_data <- stalk.obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))


if (identical(rownames(gene_vln_data1),rownames(stalk.obj@meta.data))){
  gene_vln_data1$Treatment =stalk.obj@meta.data$Treatment
}else{
  print("not same rownames!")
}


dir.create(file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln')))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))

for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
     ggtheme+ggtitle(gene_vln[small_num])+
     theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 10)) +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  ggsave(file = file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln'),paste0('angiogenesis_',gene_vln[small_num],'.pdf')), 
         width = 12, height =8)
}


################# fluid  sheer stress
gene_vln <- c("DUSP1","CAV1","KLF2","BMPR2","CAV2","THBD")

gene_vln_data <- stalk.obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))

if (identical(rownames(gene_vln_data1),rownames(stalk.obj@meta.data))){
  gene_vln_data1$Treatment =stalk.obj@meta.data$Treatment
}else{
  print("not same rownames!")
}

dir.create(file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln')))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))

for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
     ggtheme+ggtitle(gene_vln[small_num])+
     theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 10)) +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  ggsave(file = file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln'),paste0('fluid_',gene_vln[small_num],'.pdf')), 
         width = 12, height =8)
}



#####################vascular endothelial growth factor receptor signaling pathway
use.sub <- 'Lymphatic ECs'
lym.obj <- subset(obj, cells = 
                      rownames(obj@meta.data
                               [obj@meta.data$subtypes == "Lymphatic ECs" ,
                                 
                               ]))


gene_vln <-  c("FLT4","DAB2IP","ROCK2","PAK2","HIF1A","PRKD2")
gene_vln_data <- lym.obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))

if (identical(rownames(gene_vln_data1),rownames(lym.obj@meta.data))){
  gene_vln_data1$Treatment =lym.obj@meta.data$Treatment
}else{
  print("not same rownames!")
}

dir.create(file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln')))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))

for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
     ggtheme+ggtitle(gene_vln[small_num])+
     theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 10)) +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  ggsave(file = file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln'),paste0('VGF_',gene_vln[small_num],'.pdf')), 
         width = 12, height =8)
}


############# Leukocyte transendothelial migration

gene_vln <- c("ROCK2","PLCG1","ARHGAP35","MAPK14","PXN","GNAI3")
gene_vln_data <- lym.obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))

if (identical(rownames(gene_vln_data1),rownames(lym.obj@meta.data))){
  gene_vln_data1$Treatment =lym.obj@meta.data$Treatment
}else{
  print("not same rownames!")
}

#dir.create(file.path(dir_path,   'Endo1', paste0(use.sub,'_Gene-baxplot')))
dir.create(file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln')))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))

for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
      ggtheme+ggtitle(gene_vln[small_num])+
      theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 10)) +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  ggsave(file = file.path(dir_path,   'Endo', paste0(use.sub,'_Gene-vln'),paste0('Leukocyte_',gene_vln[small_num],'.pdf')), 
         width = 12, height =8)
}

