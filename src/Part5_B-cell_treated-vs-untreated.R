




source('/src/utilities.R')

dir_path <- '/data/part5'
setwd(dir_path)


immune_colors <- c("#ec7696","#b2cf87","#a35c8f","#66a9c9","#e2c027")
names(immune_colors) <- c("B","Mast","Myeloid","NK","T")



Patient_ID_ID  <- c("P5","P6","P7","P8")
immune_obj <- subset(merged.se.obj, cells = rownames(merged.se.obj@meta.data[merged.se.obj@meta.data$all_CellType.1 == c("T","B" , "Myeloid", "Mast","NK")   &
                                                                          merged.se.obj@meta.data$condition == 'Tumor'&
                                                                          merged.se.obj@meta.data$Patient_ID %in% Patient_ID_ID, 
]))


### 
static_cell <- function(type,PaID){
  num_subtypes = length(which( immune_obj@meta.data$all_CellType.1 == type & immune_obj@meta.data$Patient_ID == PaID))
  num_Patient_ID= length(which(immune_obj@meta.data$Patient_ID == PaID ))

  percet_subtypes = num_subtypes/num_Patient_ID
  static_percet_subtypes = sprintf("%0.4f", percet_subtypes)
  return(static_percet_subtypes)
  
}

P_inf_PID <- c()
P_inf_subtypes <- c()
P_inf_percent <- c()

Patient_ID_ID <-  c("P6","P5","P7","P8") 

for (PID in c(Patient_ID_ID)){
  for (i in c("T","B" , "Myeloid", "Mast","NK"          )){
    P_inf_PID <- c(P_inf_PID,PID)
    P_inf_subtypes <- c(P_inf_subtypes,i)
    P_inf_percent <- c(P_inf_percent,static_cell(i,PID))
  }
  
}

perc_subtypes  <- data.frame(Patient = P_inf_PID,
                             Percent = as.numeric((P_inf_percent)),
                             subtypes = P_inf_subtypes
)

head(perc_subtypes)


perc_subtypes$Patient1 <- factor(perc_subtypes$Patient,levels=c("P8","P7","P5","P6") )
perc_subtypes$subtypes1 <- factor(perc_subtypes$subtypes,levels=c("T","B", "Myeloid", "Mast","NK" ) )

p<- ggplot(data = perc_subtypes,aes(Patient1,weight=Percent,fill=factor(subtypes1)))+geom_bar(width=0.6)+theme_gray()+scale_fill_manual(values = immune_colors)+
  theme(legend.title = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14))+ylab("Cellular fraction")+coord_cartesian(ylim = c(0,1.0)) +coord_flip()


pdf(
              paste0('Percent_Patient_tumor-tissue_immune-cell_naive-treat.pdf'),10,8)
print(p)
dev.off()


write.csv(perc_subtypes,file=file.path(dir_path, 
                                       paste0('Percent_Patient_tumor-tissue_immune-cell_naive-treat.csv')))



############### B cell 
Patient_ID_ID  <- c("P5","P6","P7","P8")
B_obj <- subset(merged.se.obj, cells = rownames(merged.se.obj@meta.data[merged.se.obj@meta.data$all_CellType.1 == 'B' &
                                                                          merged.se.obj@meta.data$condition == 'Tumor'&
                                                                          merged.se.obj@meta.data$Patient_ID %in% Patient_ID_ID, 
]))

CellType_Type <- 'B'
immuB_obj <- B_obj 

immuB_obj<- NormalizeData(immuB_obj, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
immuB_obj <- FindVariableFeatures(immuB_obj, 
                                  selection.method = "vst", nfeatures = 2000)
immuB_obj <- ScaleData(immuB_obj)
gc()
immuB_obj<- RunPCA(immuB_obj, 
                   features = VariableFeatures(object = immuB_obj))
gc()

ElbowPlot(immuB_obj)
ggsave(file=file.path(dir_path,   'elbowplot.pdf'))

immuB_obj <- FindNeighbors(immuB_obj, dims = 1:15)
immuB_obj <- FindClusters(immuB_obj, resolution = 0.7)    
immuB_obj <- RunUMAP(immuB_obj, dims = 1:15)


table(immuB_obj@meta.data$Patient_ID)
table(immuB_obj@meta.data$Treatment)
table(immuB_obj@meta.data$seurat_clusters)



use.markers <- c('XBP1', 'CD38',#plasma
                 'MS4A1','CXCR4','HLA-DRA','HLA-DRB1', #follicular
                 
                 'CD2','CD7'# pro-b  
                 
)

dir.create(file.path( dir_path,  
                     'markers'))


plot_markers <- function(x){
  FeaturePlot(immuB_obj, features = x,
              label=T) + NoLegend()
  ggsave(file =file.path( dir_path,  
                         'markers',paste0(x, '.png')))
}



for(i in c(use.markers)){
  plot_markers(i)
}



pdf(file.path(dir_path, 
              paste0('umap.cluster.1.pdf')),7.5,6)
DimPlot(immuB_obj, group.by = 'seurat_clusters',label =T, label.size =5) 
dev.off()


cellanno <- immuB_obj@meta.data 
cellanno$cell.name <- rownames(cellanno)
cellanno <- cellanno %>% mutate(subtypes = 
                                  ifelse(RNA_snn_res.0.7 %in% c(0,4,5,7), 'Plasma cells', 
                                         ifelse(RNA_snn_res.0.7  %in% c(6,8), 'Pro-B','Follicular B cells'
                                         )))



if(identical(cellanno$cell.name, rownames(immuB_obj@meta.data))){
  immuB_obj@meta.data$subtypes <- cellanno$subtypes
}


###===========================================

B_colors <- c("#fca106","#d2b116","#96c24e")
names( B_colors) <-c("Pro-B","Plasma cells","Follicular B cells")


pdf(file.path(dir_path,   paste0("subtypes", '.umap.pdf')),8,6)
DimPlot(immuB_obj, group.by = 'subtypes', label=T, label.size=5,cols=B_colors)
dev.off()

pdf(file.path(dir_path,   paste0("subtypes_no_label", '.umap.pdf')),8,6)
DimPlot(immuB_obj, group.by = 'subtypes', label=F, label.size=5,cols=B_colors)
dev.off()


patient_colors <- c("#91D1C2FF", "#8491B4FF", "#F39B7FFF" ,"#3C5488FF")
names(patient_colors) <- c("P5","P6","P7" ,"P8")

pdf(file.path(dir_path,   paste0("Patient_ID",'.umap.pdf')),8,6)
DimPlot(immuB_obj, group.by ='Patient_ID', cols = patient_colors)
dev.off()


pdf(file.path(dir_path,  paste0( "Treatment",'.umap.pdf')),8,6)
DimPlot(immuB_obj, group.by = 'Treatment', cols =treatment_colors)
dev.off()

saveRDS(immuB_obj,"immuB_obj.RDS")

########num
table(immuB_obj@meta.data$Treatment )
table(immuB_obj@meta.data$seurat_clusters )
table(immuB_obj@meta.data$subtypes )
length(which(immuB_obj@meta.data$Treatment == "Treatment_naive" & immuB_obj@meta.data$subtypes == "Plasma cells"))
length(which(immuB_obj@meta.data$Treatment == "Chemotherapy" & immuB_obj@meta.data$subtypes == "Plasma cells"))
table(immuB_obj@meta.data[immuB_obj@meta.data$subtypes == "Plasma cells",]$Patient_ID )



###################### heattmap

Idents(immuB_obj) <- 'seurat_clusters'
avg.nec <- log1p(AverageExpression(immuB_obj, verbose = FALSE)$RNA)


clust_subtypes <- data.frame(
  subtypes <- c(
    "Plasma cells","Plasma cells","Plasma cells","Plasma cells",
    "Pro-B","Pro-B",
    "Follicular B cells","Follicular B cells","Follicular B cells","Follicular B cells"
  )
)

colnames(clust_subtypes) <- "subtypes"         
rownames(clust_subtypes) <- c('0','4','5','7', 
                              '6','8',  
                              '1','2','3','9')

top_anno <- HeatmapAnnotation(
  df = clust_subtypes,

  col = list('subtypes' = B_colors))

pdf(file.path(dir_path,  
              paste0('sub.marker.pdf')),width =6, height =4)
pheatmap(avg.nec[use.markers, c('0','4','5','7', 
                                '6','8',  
                                '1','2','3','9')],
         
         top_annotation = top_anno, 
         color = colorRampPalette(c("#004e66", "#e1eef6", "#ff5f2e"))(200),
         scale = 'row', cluster_cols =F , cluster_rows =F,
         cellwidth = 20, cellheight = 20)

dev.off()


###################  DEG
Idents(immuB_obj) <- 'subtypes'
use.sub <- 'Plasma cells'
DE_genes <- FindMarkers(immuB_obj, 
                        ident.1 = "Chemotherapy", 
                        ident.2 = "Treatment_naive",
                        group = 'Treatment',
                        subset.ident = use.sub,
                        verbose = FALSE,
                        logfc.threshold=0.25,
                        min.pct = 0.1,
                        min.diff.pct = 0.1)


write.csv(DE_genes, file = file.path( dir_path,  
                                     paste0(use.sub, '.DE.genes.csv')))

head(DE_genes)
nrow(DE_genes)
nrow(DE_genes[DE_genes$avg_log2FC > 0,])

subobj <- subset(immuB_obj, idents = use.sub)

Idents(subobj) <- 'Treatment'
avg.nec <- log1p(AverageExpression(subobj
                                   , verbose = FALSE)$RNA)


pdf(file.path(dir_path,  
              paste0(use.sub, '.DE.genes.heatmap.pdf')), width =4, height =8)
colors = colorRampPalette(c("#3b8686", "#79bd9a", "#cff09e", "#f9cdad","#fec9c9","#fe4365"))(200)
df <- avg.nec[head(rownames(DE_genes),100),] 

df <- df[,c(2,1)]
pheatmap(as.matrix(df[,c(2,1)]), color =colors,
         fontsize =6, cluster_col =F)
dev.off()

dim(DE_genes)
length(which(DE_genes$avg_log2FC > 0))
length(which(DE_genes$avg_log2FC < 0))



################    heatmap
immuB_obj <-readRDS("immuB_obj.RDS")
DE_genes <- read.table(file.path( dir_path,  
           paste0(use.sub, '.DE.genes.csv')),row.names=1,header=T,sep=',')
head(DE_genes)
nrow(DE_genes)
nrow(DE_genes[DE_genes$avg_log2FC > 0,])

subobj <- subset(immuB_obj, idents = use.sub)

Idents(subobj) <- 'Treatment'
avg.nec <- log1p(AverageExpression(subobj
                                   , verbose = FALSE)$RNA)


pdf(file.path(dir_path,  
              paste0(use.sub, '.DE.genes.heatmap.pdf')), width =4, height =8)
colors = colorRampPalette(c("#3b8686", "#79bd9a", "#cff09e", "#f9cdad","#fec9c9","#fe4365"))(200)
df <- avg.nec[head(rownames(DE_genes),100),] 
df <- df[,c(1,2)]
pheatmap(as.matrix(df[,c(2,1)]), color =colors,
         fontsize =6, cluster_col =F,border=F)
dev.off()



#GO =================================================

GO_genes <- rownames(DE_genes)

ensembl.genes <- mapIds(org.Hs.eg.db, keys = GO_genes, 
                        keytype = "SYMBOL", column="ENSEMBL")

ensembl.genes <- na.omit(ensembl.genes)

GO<- enrichGO(gene = ensembl.genes,
              OrgDb  = org.Hs.eg.db, keyType = "ENSEMBL",
              ont  = 'ALL', pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
GO <- setReadable(GO, 'org.Hs.eg.db', 'ENSEMBL')

write.csv(GO, file = file.path( dir_path,  
                               paste0(use.sub, '.GO.csv')))

pdf(file.path( dir_path,  
              paste0(use.sub, '.GO.pdf')),10,6)
dotplot(GO, showCategory = 25)
dev.off()

pdf(file.path( dir_path,  
               paste0(use.sub, '.GO.pdf')),10,6)
barplot(GO, showCategory = 25)
dev.off()



######################### GSVA
use.sub <- "Plasma cells"
subobj <- subset(al.integrated, idents = use.sub)
cluster1 <- subobj
expr <- as.data.frame(cluster1@assays$RNA@data)

meta <- cluster1@meta.data[, "Treatment"]

expr=as.matrix(expr)


kegggmt <- read.gmt("/data/h.all.v7.2.symbols.gmt")


colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)

kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=15)

saveRDS(kegg2, file = file.path(dir_path,   
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
ggsave(file =file.path(dir_path,  
                       paste0(use.sub, '.gsva.path.pdf')), width =5, height =4)

allDiff$overlappedDEgene <- as.character(allDiff$overlappedDEgene)
write.csv(allDiff, file = file.path(dir_path,   
                                    paste0(use.sub, '.gsva.diff.csv')))






kegg3 <- scale(as.matrix(as.data.frame(kegg2) ))

kegg4 <- kegg3[paste0('HALLMARK_',allDiff$pathway),]

annotation_col <- cluster1@meta.data[, c('Treatment', 'Patient_ID')]
rownames(annotation_col) <- rownames(cluster1@meta.data)

rownames(kegg4) <- gsub('HALLMARK_', '', rownames(kegg4))

kk2 <- pheatmap(kegg4,
                border_color = NA,
                cluster_rows = T,
                cluster_cols = T,
                annotation_col =annotation_col,
                annotation_colors = list(Treatment= treatment_colors,
                                         Patient_ID = patient_colors[unique(annotation_col$Patient_ID)]),
                annotation_legend=T, 
                show_rownames = T,
                show_colnames = F,
                color =colorRampPalette(c("blue", "white","red"))(100),
                cellwidth = 0.4, cellheight = 13,
                fontsize = 10)
pdf(file.path(dir_path,  
              paste0(use.sub, ".gsva.heatmap.pdf")),width = 11,height = 6)
kk2
dev.off()


######## GSVA vln
kegg3 <- as.data.frame(kegg2) 
rownames(kegg3) <- gsub('HALLMARK_', '', rownames(kegg3))
kegg4 <- kegg3[paste0(allDiff$pathway),]
gene_vln <- rownames(kegg4)

gene_vln_data1 <- as.data.frame(t( kegg4[gene_vln,] ))

if (identical(rownames(gene_vln_data1),rownames(subobj@meta.data))){
  gene_vln_data1$Treatment = subobj$Treatment
}else{
  print("not same rownames!")
}
colnames(gene_vln_data1) <- c(gene_vln,"Treatment")

gene_vln_data1$Treatment <- factor(gene_vln_data1$Treatment,levels=c("Treatment_naive","Chemotherapy"))
compare <- 'ChemoVS.Naive'
for (small_num in seq(1,length(gene_vln)) ){  
  ggplot(gene_vln_data1,aes_string(x= "Treatment", y = gene_vln[small_num], fill ="Treatment")) +
    geom_violin() + scale_fill_manual(values = treatment_colors) +
    geom_boxplot(width=.02,col="black",fill="white")+
    ggtheme+ggtitle(gene_vln[small_num])+
     theme(axis.text.x = element_text(vjust = 0.5),
          plot.title=element_text(hjust = 0.5,face = "bold",size = 10)) +
    
    geom_signif(test = "wilcox.test",comparisons =  list(c("Treatment_naive", "Chemotherapy")))
  
  ggsave(file = file.path(dir_path,  ,'GSVA_DEG',paste0('vln_',compare,'.',gene_vln[small_num],'.pdf')), 
         width = 12, height =8)
}




