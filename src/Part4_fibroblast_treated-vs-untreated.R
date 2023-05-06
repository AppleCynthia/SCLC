



library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(snow)



source('/src/utilities.R')

dir_path <- '/data/part4'
setwd(dir_path)


celltype <- 'Fib'
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

  return(se.obj)
}

if(!dir.exists(file.path(dir_path,  'Fib'))){
  dir.create(file.path(dir_path,  'Fib'))
}else{
  print('directory exist')
}


obj <- run_seurat(obj)

ElbowPlot(obj)
ggsave(file=file.path(dir_path, 'Fib','elbowplot.pdf'))
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:15)


###===================================================================

use.markers <- c(  'ACTA2', 'MYH11', 'TAGLN', 'PLN',
                   'COL11A1','COL10A1','CNN1', 'MMP11', 'SYNPO2',# smoothe muscle cell
                   'RGS5','MCAM', 'CSPG4', 'HIGD1B','FAM162B',# pericyte
                   'CFD', 'C3', 'DCN', 'CLU', 'FBLN1', 'APOD','BMP5','TCF21', #Fibroblast
                   'MKI67', 'CCNB2','CDKN3',# prolife fib
                   
                   'CXCL9', 'CCL19', 'CXCL10','CXCL11'#active fib        
)



VlnPlot(obj, features = use.markers)
ggsave(file = file.path(dir_path,  'Fib',paste0(celltype, '.sub.markers.pdf')), 
  width = 12, height =18)

FeaturePlot(obj, features = use.markers)
ggsave(file = file.path(dir_path,  'Fib', paste0(celltype, '.sub.markers.featureplot.pdf')), 
       width = 12, height =18)


pdf(file.path(dir_path, 'Fib',paste0(celltype, '.umap.pdf')),7.5,6)
DimPlot(obj, group.by = 'seurat_clusters',label =T, label.size =5)
dev.off()


##annotate sub celltye========================
cellanno <- obj@meta.data 
cellanno$cell.name <- rownames(cellanno)
cellanno <- cellanno %>% mutate(subtypes = 
        ifelse(RNA_snn_res.0.5 %in% c(0,3,9,11), 'Fibroblast', 
          ifelse(RNA_snn_res.0.5 ==10, 'Proliferating Fib',
              ifelse(RNA_snn_res.0.5 %in% c(2,5,7), 'Pericytes', 
                ifelse(RNA_snn_res.0.5 %in%c(1,4,6), 'Myofibrolast',
                      ifelse(RNA_snn_res.0.5 == 8, 'Active Fib',
                      'Undetermined'))))))

if(identical(cellanno$cell.name, rownames(obj@meta.data))){
  obj@meta.data$subtypes <- cellanno$subtypes
}


saveRDS(obj,'Fib.rds')
##============================================
pdf(file.path(dir_path,  'Fib',paste0(celltype, '.subtypes.umap.pdf')),7.5,6)
DimPlot(obj, group.by ='subtypes', cols = )
dev.off()


pdf(file.path(dir_path,  'Fib',paste0(celltype, '.P_ID.umap.pdf')),7.5,6)
DimPlot(obj, group.by ='Patient_ID', cols = patient_colors)
dev.off()

pdf(file.path(dir_path, 'Fib',paste0(celltype, '.Treatment.umap.pdf')),7.5,6)
DimPlot(obj, group.by = 'Treatment', cols =treatment_colors)
dev.off()


Idents(obj) <- 'subtypes'

use.markers <- c(  'ACTA2', 'MYH11', 'TAGLN', 'PLN',
                'COL11A1','COL10A1','CNN1', 'MMP11', 'SYNPO2',#smoothe muscle cell
                'RGS5','MCAM', 'CSPG4', 'HIGD1B','FAM162B',#pericyte
                  'CFD', 'C3', 'DCN', 'CLU', 'FBLN1', 'APOD','BMP5','TCF21', #Fibroblast
                  'MKI67', 'CCNB2','CDKN3',#prolife fib
 
                   'CXCL9', 'CCL19', 'CXCL10','CXCL11'#active fib        
                )


avg.nec <- log1p(AverageExpression(obj, verbose = FALSE)$RNA)


pdf(file.path(dir_path,  'Fib',
 paste0('sub.marker.pdf')), width =4, height =6)
pheatmap(avg.nec[use.markers, ],
color = colorRampPalette(c("#004e66", "#e1eef6", "#ff5f2e"))(200),
scale = 'row', cluster_cols =F , cluster_rows =F)

dev.off()


#DE in subtypes++++++++

 use.sub <- "Fibroblast"
 use.sub <- "Myofibrolast"
 use.sub <- "Pericytes"

 use.sub <- "Fib"

DE_genes <- FindMarkers(obj, 
  ident.1 = "Chemotherapy",
  ident.2 = "Treatment_naive", 
  group = 'Treatment',
   verbose = FALSE,
   logfc.threshold=0.25,
                min.pct = 0.3,
                min.diff.pct = 0.2)


write.csv(DE_genes, file = file.path(dir_path,  'Fib',
  paste0(use.sub, '.DE.genes.csv')))


dim(DE_genes)
dim(DE_genes[DE_genes$avg_log2FC > 0,])
dim(DE_genes[DE_genes$avg_log2FC < 0,])



Idents(subobj) <- 'Treatment'
avg.nec <- log1p(AverageExpression(subobj
  , verbose = FALSE)$RNA)


pdf(file.path(dir_path,  'Fib',
 paste0(use.sub, '.DE.genes.heatmap.pdf')), width =4, height =8)
colors = colorRampPalette(c("#3b8686", "#79bd9a", "#cff09e", "#f9cdad","#fec9c9","#fe4365"))(200)
df <- avg.nec[c(head(rownames(DE_genes),50), tail(rownames(DE_genes),50)),] #Fib,B
pheatmap(df[,1:2], color =colors,
  fontsize =6, cluster_col =F)
dev.off()


#GO and KEGG=================================================
########### GO
GO_genes <- rownames(DE_genes)

ensembl.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
  keytype = "SYMBOL", column="ENSEMBL")

ensembl.genes <- na.omit(ensembl.genes)

GO<- enrichGO(gene = ensembl.genes,
                  OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont  = 'ALL', pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
GO <- setReadable(GO, 'org.Hs.eg.db', 'ENSEMBL')

write.csv(GO, file = file.path(dir_path,  'Fib',
  paste0(use.sub, '.GO.csv')))

pdf(file.path(dir_path, 'Fib',
 paste0(use.sub, '.dot.GO.pdf')),20,6)
dotplot(GO, showCategory = 25)
dev.off()

pdf(file.path( dir_path,  'Fib',
               paste0(use.sub, '.bar.GO.bar.pdf')),10,6)
barplot(GO, showCategory = 25)
dev.off()

#######KEGG
entrez.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
  keytype = "SYMBOL", column="ENTREZID")
kk <- enrichKEGG(gene         = entrez.genes,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
write.csv(edox, file = file.path(dir_path,  'Fib',
 paste0(use.sub, '.kegg.csv')))

pdf(file = file.path(dir_path,  'Fib',
  paste0(use.sub, '.KEGG.pdf')),15,15)
cnetplot(edox, showCategory = 5)
dev.off()

pdf(file = file.path(dir_path,   'Fib',
                     paste0(use.sub, '.dot.KEGG.pdf')))
dotplot(edox, showCategory = 15)
dev.off()

pdf(file = file.path(dir_path,   'Fib',
                     paste0(use.sub, '.bar.KEGG.pdf')),10,6)
barplot(edox, showCategory = 15)
dev.off()



#GSVA==========================================================
use.sub <- "Fibroblast"
use.sub <- "Myofibrolast"
use.sub <- "Pericytes"

cluster1 <- subset(obj, cells = rownames(
  obj@meta.data[obj@meta.data$subtypes == use.sub, ]))

cluster1 <- obj
expr <- as.data.frame(cluster1@assays$RNA@data)

meta <- cluster1@meta.data[, "Treatment"]

expr=as.matrix(expr)



kegggmt <- read.gmt("/data/h.all.v7.2.symbols.gmt")


colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)

 kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=15)

 saveRDS(kegg2, file = file.path(dir_path,  'Fib',
  paste0(use.sub,'.gsva.rds')))


exprSet <- kegg2

group <- factor(meta,levels = c("Treatment_naive",
  "Chemotherapy"),ordered = F)

design <- model.matrix(~group)

colnames(design) <- levels(group)

fit <- lmFit(exprSet,design)

fit2 <- eBayes(fit)

allDiff=topTable(fit2, adjust='fdr',coef=2,number=30,
  p.value=0.05, lfc = 0.1)



pathway_genes <- function(path, DE.gene){
  genes <- intersect(kegg_list[path][[1]],
    DE.gene)
  return(c(genes))
}

allDiff$pathway <- gsub('HALLMARK_', '', rownames(allDiff))
  
DE.gene <- rownames(DE_genes)
allDiff$overlappedDEgene <- lapply(rownames(allDiff), function(x)pathway_genes(x, DE.gene))

allDiff <- allDiff %>% mutate(group1 = ifelse(t <= 0, 'Down', 'Up'))
allDiff <- allDiff[order(allDiff$t, decreasing =T),]



y_nam <- allDiff[order(allDiff$t, decreasing =F),]$pathway

 ggplot(allDiff, aes(pathway, t, fill=factor(group1))) + 
  geom_bar(stat = 'identity') + 
    scale_fill_manual(values = as.character(treatment_colors)) +
    scale_x_discrete(limits = y_nam)  +

    theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
      legend.position = 'none', 
      legend.text=element_text(size=12),
      axis.title.x = element_blank(),
      axis.text.x =element_text(size=12),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size=8),
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))+
    coord_flip() 
ggsave(file =file.path(dir_path,  'Fib',
   paste0(use.sub, '.gsva.path.pdf')), width =4, height =2)

allDiff$overlappedDEgene <- as.character(allDiff$overlappedDEgene)
write.csv(allDiff, file = file.path(dir_path,  'Fib',
   paste0(use.sub, '.gsva.diff.csv')))


kegg3 <- as.data.frame(kegg2) 
kegg4 <- kegg3[paste0('HALLMARK_',allDiff$pathway),]

annotation_col <- cluster1@meta.data[, c('Treatment', 'Patient_ID')]
rownames(annotation_col) <- rownames(cluster1@meta.data)

rownames(kegg4) <- gsub('HALLMARK_', '', rownames(kegg4))

kk2 <- pheatmap(kegg4,
                border_color = NA,
                cluster_rows = F,
                cluster_cols = T,
                annotation_col =annotation_col,
                annotation_colors = list(Treatment= treatment_colors,
                  P_ID = patient_colors[unique(annotation_col$P_ID)]),
                annotation_legend=T, 
                show_rownames = T,
                show_colnames = F,
                color =colorRampPalette(c("blue", "white","red"))(100),
                cellwidth = 0.2, cellheight = 13,
                fontsize = 10)
pdf(file.path(dir_path,  'Fib',
  paste0(use.sub, ".gsva.heatmap.pdf")),width = 20,height = 6)
kk2
dev.off()


###########  Ridgeline plot
library(ggridges)
library(ggplot2)
library(reshape2)


gene_vln <- c('COL10A1','COL12A1','COL15A1','COL16A1','COL18A1','COL4A1','COL4A2','COL5A1','COL5A2','COL8A1')
gene_colors <-c("#3b818c","#8fb2c9","#887322","#6e8b74","#ab372f","#617172","#74759b","#5d3f51","#894e54","#681752")
names(gene_colors) <- c('COL10A1','COL12A1','COL15A1','COL16A1','COL18A1','COL4A1','COL4A2','COL5A1','COL5A2','COL8A1')

library(cowplot)

subobj <- subset(obj,  cells = rownames(obj@meta.data[obj@meta.data$Treatment == "Treatment_naive", 
]))
gene_vln_data  <- as.data.frame(subobj@assays$RNA@data[gene_vln, ])
gene_vln_data1 <- mutate(gene_vln_data, expression = rownames(gene_vln_data))
gene_vln_data1 <- cbind(rownames(gene_vln_data ),gene_vln_data)
gene_vln_data2 <- melt(gene_vln_data1, variable.name = "Expression")
colnames(gene_vln_data2)[1] <- "Gene" 

p1 <- ggplot(gene_vln_data2, aes(x = value, y = Gene , fill = Gene )) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")+scale_fill_manual(values = gene_colors)
subobj <- subset(obj,  cells = rownames(obj@meta.data[obj@meta.data$Treatment == "Chemotherapy" 
                                                        ,
]))
gene_vln_data <- as.data.frame(subobj@assays$RNA@data[gene_vln, ])
gene_vln_data1 <- mutate(gene_vln_data, expression = rownames(gene_vln_data))
gene_vln_data1 <- cbind(rownames(gene_vln_data ),gene_vln_data)
gene_vln_data2 <- melt(gene_vln_data1, variable.name = "Expression")
colnames(gene_vln_data2)[1] <- "Gene" 

p2 <- ggplot(gene_vln_data2, aes(x = value, y = Gene , fill = Gene )) +
  geom_density_ridges() +
  theme_ridges() +  
  theme(legend.position = "none") + scale_fill_manual(values = gene_colors)

plot_grid(p1,p2,nrow=2)
#p2
ggsave(file = file.path(dir_path,  'Fib',paste0('ShanJiTu.pdf')), 
       width = 6, height =8)
    


table(obj@meta.data$subtypes)

