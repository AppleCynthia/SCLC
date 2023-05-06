


library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(snow)




dir_path <- '/data/part2'
setwd(file.path(dir_path,'NE'))

merged.se.obj <- readRDS(file.path('/data/part1', '1.merged.se.obj.RDS'))


##  NE cells in tumor tissue
obj <- subset(merged.se.obj, cells = 
                rownames(merged.se.obj@meta.data
                         [merged.se.obj@meta.data$all_CellType.1 == 'NE' &
                             merged.se.obj@meta.data$condition == 'Tumor' ,
                              
                         ]))



## run seurat
run_seurat <- function(se.obj){
	se.obj<- NormalizeData(se.obj, normalization.method = "LogNormalize", 
	scale.factor = 10000)
	se.obj <- FindVariableFeatures(se.obj, selection.method = "vst", nfeatures = 1000)
	se.obj <- ScaleData(se.obj)
	se.obj <- RunPCA(se.obj, features = VariableFeatures(object = se.obj))
		return(se.obj)
}
obj <- run_seurat(obj)


ElbowPlot(obj)
ggsave(file=file.path(dir_path,  'NE','elbowplot1.pdf'))

obj <- FindNeighbors(obj, dims = 1:17)
obj <- FindClusters(obj, resolution = 0.7)
obj <- RunUMAP(obj, dims = 1:17, min.dist= 0.4)




pdf(file.path(dir_path,   'NE', paste0( 'cluster','.umap.pdf')),7,6)
DimPlot(obj, label =T, group.by ='seurat_clusters')
dev.off()



pdf(file.path(dir_path,   'NE', paste0('Epi', 'celltype','.umap.pdf')),7,6)
DimPlot(obj, label =T, group.by ='all_CellType.1', cols =CellType_colors[unique(obj@meta.data$all_CellType.1)])
dev.off()


pdf(file.path(dir_path,   'NE', paste0('Epi', 'P_ID','.umap.pdf')),7,6)
DimPlot(obj, group.by ='Patient_ID', cols = patient_colors)
dev.off()


pdf(file.path(dir_path,   'NE', paste0('Epi', 'treatment','.umap.pdf')),7,6)
DimPlot(obj, group.by = 'Treatment', cols =treatment_colors)
dev.off()

saveRDS(obj,file.path(dir_path, 'NE',  'NE.rds'))
celltype <- 'NE'



################# DE analysis between treated and untreated patients

Idents(obj) <- 'all_CellType.1'
use.sub <-'NE'

DE_genes <- FindMarkers(obj, 
	ident.1 = "Chemotherapy",
	ident.2 = "Treatment_naive", 
	group = 'Treatment',
  subset.ident = use.sub,
	 verbose = FALSE,
	 logfc.threshold=0.25,
                min.pct = 0.1,
                min.diff.pct = 0.2)

 write.csv(DE_genes, file = file.path(dir_path,   'NE',
	paste0(use.sub, '.DE.genes.csv')))



avg.nec <- log1p(AverageExpression(obj, verbose = FALSE,subset.ident = "NE",group.by="Treatment")$RNA)


length(which(DE_genes$avg_log2FC > 0))
length(which(DE_genes$avg_log2FC < 0))


pdf(file.path(dir_path,  'NE',
 paste0(use.sub, '.DE.genes.heatmap.pdf')), width =4, height =8)
colors = colorRampPalette(c("#3b8686", "#79bd9a", "#cff09e", "#f9cdad","#fec9c9","#fe4365"))(200)
df <- avg.nec[rownames(DE_genes),] 
df1 <- as.data.frame(df[,c(2,1)])
pheatmap(df1[,1:2], color =colors,
	fontsize =6, cluster_col =F)
dev.off()
write.csv(avg.nec, file = file.path(dir_path,   'NE',
                                     paste0('avg.nec.csv')))



#########GO

GO_genes <- rownames(DE_genes)

ensembl.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
	keytype = "SYMBOL", column="ENSEMBL")

ensembl.genes <- na.omit(ensembl.genes)

GO<- enrichGO(gene = ensembl.genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = 'ALL', pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
GO <- setReadable(GO, 'org.Hs.eg.db', 'ENSEMBL')

write.csv(GO, file = file.path(dir_path,  'NE',
	paste0('NE', '.GO.csv')))

pdf(file.path(dir_path,   'NE',
 paste0('NE', '.GO.pdf')),10,6)
dotplot(GO, showCategory = 20)
dev.off()

pdf(file.path( dir_path,  'NE',
               paste0(use.sub, '.GO.bar.pdf')),10,6)
barplot(GO, showCategory = 25)
dev.off()


#########KEGG 
entrez.genes <- mapIds(org.Hs.eg.db, keys =  rownames(DE_genes), 
	keytype = "SYMBOL", column="ENTREZID")
kk <- enrichKEGG(gene         = entrez.genes,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
write.csv(edox, file = file.path(dir_path,  'NE',
 paste0('NE', '.kegg.csv')))

pdf(file = file.path(dir_path,   'NE',
	paste0('NE', '.KEGG.pdf')))
cnetplot(edox, showCategory = 15)
dev.off()

pdf(file.path(dir_path,   'NE',
              paste0('NE', '.KEGG.dot.pdf')),10,6)
dotplot(edox, showCategory = 15)
dev.off()



gene_vln_0 <- c("NEAT1","NFKBIA","C2CD4B","GPRC5B","STAT3","HLA-E",
          "IRF1","HLA-A","HLA-E",
          "CD74","HLA-A","HLA-E","HLA-B")


gene_vln <- unique(gene_vln_0)


gene_vln_data <- obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))

if (identical(rownames(gene_vln_data1),rownames(obj@meta.data))){
  gene_vln_data1$Treatment =obj@meta.data$Treatment
}else{
  print("not same rownames!")
}

dir.create(file.path(dir_path, 'NE',  'Gene-vln'))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))

gene_vln[6] <- "HLA.E"
colnames(gene_vln_data1)[6] <- "HLA.E"
gene_vln[8] <- "HLA.A"
colnames(gene_vln_data1)[8] <- "HLA.A"
gene_vln[10] <- "HLA.B"
colnames(gene_vln_data1)[10] <- "HLA.B"


for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
       ggtheme+ggtitle(gene_vln[small_num])+
    theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 20)) +
    labs(y="Gene expression") +theme(axis.title.y = element_text(vjust = 2, size = 15,face = "bold"))  + 
    guides(fill = F)  +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
   
  ggsave(file = file.path(dir_path, 'NE',  'Gene-vln',paste0(gene_vln[small_num],'.pdf')), 
         width = 8, height =6)
}
 


###########   DNA repair  score 

gene.data <- as.data.frame(gene <- repair)
gene <-  as.list(gene.data)
gene
emt.obj <- AddModuleScore(object = obj,features = gene,crl=100,name = "repair_feature")
head(emt.obj@meta.data)

emt.obj@meta.data$Treatment1 <- factor(emt.obj@meta.data$Treatment,levels=c("Treatment_naive","Chemotherapy"))
pdf(file.path(dir_path,   'NE',
              paste0('score.repair.pdf')),8,8)
ggplot(emt.obj@meta.data,aes_string(x= "Treatment1", y = "repair_feature1", fill ="Treatment1")) +
  geom_violin() + scale_fill_manual(values = treatment_colors) +
  geom_boxplot(width=.02,col="black",fill="white")+
  ggtheme+ggtitle("DNA Repair")+
  theme(axis.text.x = element_text(vjust = 0.5),
        plot.title=element_text(hjust = 0.5,face = "bold",size = 20)) +
  labs(y="Score") +theme(axis.title.y = element_text(vjust = 2, size = 15,face = "bold"))  + 
  guides(fill = F)  +
  geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
dev.off()




######### target score
therapy_targets <- c('MYC','BCL2','TOP2B','KDM1A','DLL3','TOP1','PARP1','TOP2A','AURKB','AURKA','CHEK1','EZH2','VEGFA')

gene.data <- as.data.frame(gene <- therapy_targets)
gene <-  as.list(gene.data)
gene
emt.obj <- AddModuleScore(object = obj,features = gene,crl=100,name = "therapy_targets_feature")
head(emt.obj@meta.data)

emt.obj@meta.data$Treatment1 <- factor(emt.obj@meta.data$Treatment,levels=c("Treatment_naive","Chemotherapy"))
pdf(file.path(dir_path,   'NE',
              paste0('score.therapy_targets.pdf')),8,8)
ggplot(emt.obj@meta.data,aes_string(x= "Treatment1", y = "therapy_targets_feature1", fill ="Treatment1")) +
  geom_violin() + scale_fill_manual(values = treatment_colors) +
  geom_boxplot(width=.02,col="black",fill="white")+
  ggtheme+ggtitle("Therapy Targets")+
  theme(axis.text.x = element_text(vjust = 0.5),
        plot.title=element_text(hjust = 0.5,face = "bold",size = 20)) +
  labs(y="Score") +theme(axis.title.y = element_text(vjust = 2, size = 15,face = "bold"))  + 
  guides(fill = F)  +
  geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
dev.off()


######################### violin plot of cellular senescence associated genes 


Cellular_senescence <- c("GADD45A","GADD45B")
gene_vln <- Cellular_senescence
gene_vln_data <- obj@assays$RNA@data[gene_vln, ]
gene_vln_data1 <- as.data.frame(t( gene_vln_data ))

if (identical(rownames(gene_vln_data1),rownames(obj@meta.data))){
  gene_vln_data1$Treatment =obj@meta.data$Treatment
}else{
  print("not same rownames!")
}

dir.create(file.path(dir_path, 'NE',  'Gene-vln'))
gene_vln_data1$Treatment1 <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))


for (small_num in seq(1,length(gene_vln))){
  ggplot(gene_vln_data1,aes_string(x= "Treatment1", y = gene_vln[small_num], fill ="Treatment1")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
     ggtheme+ggtitle(gene_vln[small_num])+
     theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 20)) +
    labs(y="Gene expression") +theme(axis.title.y = element_text(vjust = 2, size = 15,face = "bold"))  + 
    guides(fill = F)  +
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  
  ggsave(file = file.path(dir_path, 'NE',  'Gene-vln',paste0(gene_vln[small_num],'.pdf')), 
         width = 8, height =6)
}



  