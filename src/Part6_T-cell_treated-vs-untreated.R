


library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(snow)


source('/src/utilities.R')

dir_path <- '/data/part6'
setwd(dir_path)

merged.se.obj <- readRDS(file.path('/Volumes/G/Prj_SCLC/219/preprocessing/', '1.merged.se.obj.RDS'))
Patient_ID_ID  <- c("P5","P6","P7","P8","P9","P10")
T_obj <- subset(merged.se.obj, cells = rownames(merged.se.obj@meta.data[merged.se.obj@meta.data$all_CellType.1 == 'T' &
                                                                        merged.se.obj@meta.data$condition == 'Tumor'&
                                                                          merged.se.obj@meta.data$Patient_ID %in% Patient_ID_ID, 
]))

table(T_obj@meta.data$Patient_ID)
table(T_obj@meta.data$Treatment )

Patient_ID_ID  <- c("P5","P6","P7","P8")
T_obj <- subset(T_obj, cells = rownames(T_obj@meta.data[T_obj@meta.data$Patient_ID %in% Patient_ID_ID, 
]))


if(!dir.exists(file.path(dir_path, 'markers'))){
  dir.create(file.path(dir_path, 'markers'))
}else{
  print('directory exist')
}


############################

T_obj <- CellCycleScoring(T_obj,
                          g2m.features = cc.genes$g2m.genes,
                          
                          s.features = cc.genes$s.genes)
"MLF1IP" %in%  cc.genes$s.genes
cc.genes$s.genes[which(cc.genes$s.genes == "MLF1IP")]
cc.genes$s.genes[which(cc.genes$s.genes == "MLF1IP")] <- "CENPU"

cc.genes$g2m.genes[which(cc.genes$g2m.genes == "FAM64A")]
cc.genes$g2m.genes[which(cc.genes$g2m.genes == "FAM64A")] <- "PIMREG"

cc.genes$g2m.genes[which(cc.genes$g2m.genes == "HN1")]
cc.genes$g2m.genes[which(cc.genes$g2m.genes == "HN1")] <- "JPT1"

T_obj <- CellCycleScoring(T_obj,
                          g2m.features = cc.genes$g2m.genes,
                          
                          s.features = cc.genes$s.genes)



se.obj <- T_obj 

  se.obj<- NormalizeData(se.obj, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  se.obj <- FindVariableFeatures(se.obj, selection.method = "vst", nfeatures = 2000)
  se.obj <- ScaleData(se.obj, 
                      vars.to.regress = c("S.Score", "G2M.Score"), 
                      features = rownames(se.obj))
  
  ElbowPlot(se.obj)
  ggsave(file=file.path(dir_path,'T', 'elbowplot.pdf'),width=10,height=10)
  se.obj <- RunPCA(se.obj, features = VariableFeatures(se.obj))
  
  ElbowPlot(se.obj)
  ggsave(file=file.path(dir_path,'T', 'elbowplot1.pdf'),width=10,height=10)
  
  se.obj <- FindNeighbors(se.obj, dims = 1:10)
  se.obj <- FindClusters(se.obj, resolution = 0.5)
  se.obj <- RunUMAP(se.obj, dims = 1:10)

al.list <- SplitObject(T_obj, split.by = "Patient_ID")
al.list <- lapply(X = al.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
ncells <- lapply(al.list, ncol)
al.list <- al.list[ncells > 100]
features <- SelectIntegrationFeatures(object.list = al.list)
al.list <- lapply(X = al.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = al.list, 
                                  reduction = "rpca", dims = 1:30)

al.integrated <- IntegrateData(anchorset = anchors, dims = 1:30) 

al.integrated <- ScaleData(al.integrated, verbose = FALSE)
al.integrated <- RunPCA(al.integrated, verbose = FALSE)

al.integrated <- RunUMAP(al.integrated, dims = 1:30) 


##########
pdf(file.path(dir_path, 
              paste0('umap.cluster.1.pdf')),7.5,6)
DimPlot(al.integrated, group.by = 'RNA_snn_res.0.5',label =T, label.size =5) 
dev.off()

markers <- c('CD4','CD8A','CD8B')
cyto.markers <- c('NKG7', 'CCL4', 'CST7','PRF1',
                  'GZMA','GZMB','IFNG','CCL3')

Exhau.markers <- c( 'TIGIT','LAG3', 'CTLA4')
#, 'FOXP3'
naive.markers <- c('CCR7','TCF7','LEF1','SELL')
use.markers <- c(markers, cyto.markers, Exhau.markers, naive.markers)

use.markers.1 <- c(use.markers,'FIR2DL3','FCGR3A')                                                                                                                            

DotPlot(al.integrated, features =use.markers,
        assay = 'RNA')
ggsave(file = file.path(dir_path,
                        paste0("DotPlot", 'use.markers.pdf')), 
       width = 13, height =10)


Idents(al.integrated) <- 'RNA_snn_res.0.5'
VlnPlot(al.integrated, features = use.markers,assay = 'RNA')
ggsave(file = file.path(dir_path,
                        paste0("VlnPlot", 'use.markers.pdf')), 
       width = 15, height =20)



# marker
naive.markers <- c('CCR7','TCF7','LEF1','SELL')
cyto.markers <- c('IL2','GNLY','GZMK',
                  'GZMA','PRF1','GZMB','IFNG','NKG7')
exhau.markers <- c('LAG3','TIGIT', 'PDCD1','HAVCR2','CTLA4')
treg.markers <- c('IL2RA','FOXP3','IKZF2')


use.markers <- c(markers, cyto.markers, exhau.markers, naive.markers,treg.markers)
DotPlot(al.integrated, features =use.markers,
        assay = 'RNA')
ggsave(file = file.path(dir_path,
                        paste0("DotPlot_ZEM_", 'use.markers.pdf')), 
       width = 20, height =10)


###===================================================================

## cell subtype##############################
VlnPlot(al.integrated, features = use.markers,assay = 'RNA')
ggsave(file = file.path(dir_path,
                        paste0("VlnPlot_ZZM_", 'use.markers.pdf')), 
       width = 15, height =20)


VlnPlot(al.integrated, features = c('CD3D', 'CD3E', 'TRAC', 'CD8A', 'CD8B'),assay = 'RNA')
ggsave(file = file.path(dir_path, 
                        paste0("VlnPlot_T_marker", 'use.markers.pdf')), 
       width = 15, height =20)


###===========================================



##annotate sub celltye========================
cellanno <- al.integrated@meta.data 
cellanno$cell.name <- rownames(cellanno)
cellanno <- cellanno %>% mutate(subtypes = 
                                  ifelse(RNA_snn_res.0.5 == 8, 'CD4-CD8-', 
                                         #ifelse(RNA_snn_res.0.5 ==1, 'CD8_Naive',
                                                ifelse(RNA_snn_res.0.5 %in% c(0,2,3,5,6), 'CD8_Effector',
                                                       ifelse(RNA_snn_res.0.5 ==10, 'CD8_Exhausted',
                                                              ifelse(RNA_snn_res.0.5 %in%  c(4,9), 'CD4_Naive', 
                                                                     ifelse(RNA_snn_res.0.5 ==7,'CD4_Exhausted','CD4_CD8_')
                                                              )))))
                                #)



if(identical(cellanno$cell.name, rownames(al.integrated@meta.data))){
  al.integrated@meta.data$subtypes <- cellanno$subtypes
}


FeaturePlot(al.integrated, features = use.markers)
ggsave(file = file.path(dir_path,
                        paste0('FeaturePlot', '.sub.markers.featureplot.pdf')), 
       width = 20, height =20)
##============================================



#tmp <- readRDS('/Volumes/G/Prj_SCLC/219/part3/T/al.integrated.rds')
#DotPlot(al.integrated, features =use.markers,
#        assay = 'RNA')
#ggsave(file = file.path(dir_path, 'part3', 'immune2_202205','T', 
#                       paste0("DotPlot_ZEM_1", 'use.markers.pdf')), 
#      width = 20, height =10)

#dim(tmp@meta.data)
dim(al.integrated@meta.data)

#'CD4_CD8','#FFD700', ; '#6B8E23','CD8_Naive',
T_colors <- c('#efafad','#c35691','#F4A460','#33A1C9','#2c9678','#b78d12')
names(T_colors) <- c('CD4_CD8_','CD8_Effector','CD8_Exhausted','CD4_Naive','CD4_Exhausted','CD4-CD8-')

#######################marker -T cells subtypes 

#Idents(al.integrated) <- 'subtypes'

markers <- c('CD4','CD8A','CD8B')
cyto.markers <- c('NKG7', 'CCL4', 'CST7','PRF1',
                  'GZMA','GZMB','IFNG','CCL3')

Exhau.markers <- c( 'TIGIT','LAG3', 'CTLA4')

naive.markers <- c('CCR7','TCF7','LEF1','SELL')
use.markers <- c(markers, cyto.markers, Exhau.markers, naive.markers)

Idents(al.integrated) <- 'RNA_snn_res.0.5'
avg.nec <- log1p(AverageExpression(al.integrated, verbose = FALSE)$RNA)


#avg.nec <- avg.nec[, colnames(avg.nec)!='Undetermined']

#ifelse(RNA_snn_res.0.5 ==1, 'CD8_Naive',


clust_subtypes <- data.frame(
  subtypes <- c( 'CD4-CD8-',
                 'CD4_CD8_', 
                 'CD4_Naive','CD4_Naive',
                 'CD4_Exhausted',
                 'CD8_Effector','CD8_Effector','CD8_Effector','CD8_Effector','CD8_Effector',
                 'CD8_Exhausted'
                 
  )
)       


#colnames(clust_subtypes) <- c('cluster','subtypes')
colnames(clust_subtypes) <- 'subtypes'
rownames(clust_subtypes) <- c('8',
                              '1',
                              '4', '9',
                              '7',
                              '0','2','3','5','6',
                              '10'
)

top_anno <- HeatmapAnnotation(
  df = clust_subtypes,
  #which ='column',
  col = list('subtypes' = T_colors))

pdf(file.path(dir_path, 
              paste0('sub.marker.pdf')),10,10)
pheatmap(avg.nec[use.markers, c('8',
                                '1',
                                '4', '9',
                                '7',
                                '0','2','3','5','6',
                                '10'
)],
#annotation_col = clust_subtypes,
top_annotation = top_anno,
color = colorRampPalette(c("#004e66", "#e1eef6", "#ff5f2e"))(200),
scale = 'row', 
cluster_cols =F , 
cluster_rows =F,
gaps_col = c(1,2,4,5,10),
gaps_row = c(3,11,14))

dev.off()


pdf(file.path(dir_path, 
              paste0('umap.subtypes_label.pdf')),7.5,6)
DimPlot(al.integrated, group.by ='subtypes', label =T, label.size =3, cols = T_colors)
dev.off()

pdf(file.path(dir_path, 
              paste0('umap.subtypes.pdf')),7.5,6)
DimPlot(al.integrated, group.by ='subtypes', label =F, label.size =3, cols = T_colors)
dev.off()


pdf(file.path(dir_path, 
              paste0('umap.Patient_ID.pdf')),7.5,6)
DimPlot(al.integrated, group.by ='Patient_ID', cols = patient_colors[c(5,6,7,8)])
dev.off()

pdf(file.path(dir_path, 
              paste0('umap.Treatment.pdf')),7.5,6)
DimPlot(al.integrated, group.by = 'Treatment', cols =treatment_colors)
dev.off()

pdf(file.path(dir_path,
              paste0('condition.pdf')),7.5,6)
DimPlot(al.integrated, group.by = 'condition', cols =condition_colors)
dev.off()




########################The fraction of T cell subtypes in each patient
table(al.integrated@meta.data$subtype)

static_cell <- function(type,PaID){
  num_subtypes = length(which(al.integrated@meta.data$subtypes == type & al.integrated@meta.data$Patient_ID == PaID))
  num_Patient_ID= length(which(al.integrated@meta.data$Patient_ID == PaID))
  #print(num_Effect)
  percet_subtypes = num_subtypes/num_Patient_ID
  static_percet_subtypes = sprintf("%0.4f", percet_subtypes)
  return(static_percet_subtypes)
  
}

P_inf_PID <- c()
P_inf_subtypes <- c()
P_inf_percent <- c()

Patient_ID_ID  <- c("P6","P5","P7","P8")
for (PID in c(Patient_ID_ID)){
  for (i in c( 'CD4_CD8_','CD4-CD8-','CD4_Naive', 'CD4_Exhausted','CD8_Effector', 'CD8_Exhausted' )){
    P_inf_PID <- c(P_inf_PID,PID)
    P_inf_subtypes <- c(P_inf_subtypes,i)
    P_inf_percent <- c(P_inf_percent,static_cell(i,PID))
  }
  
}

perc_subtypes  <- data.frame(Patient = P_inf_PID,
                             Percent = as.numeric((P_inf_percent)),
                             subtypes = P_inf_subtypes
)

#cl<- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")
#perc_subtypes$Patient1 <- factor(perc_subtypes$Patient,levels = c("P6","P5","P7","P8"))
#perc_subtypes$subtypes1 <- factor(perc_subtypes$subtypes,levels = c('CD8_Effector','CD8_Exhausted','CD4_CD8_','CD4-CD8-','CD4_Naive', 'CD4_Exhausted' ))


perc_subtypes$Patient1 <- factor(perc_subtypes$Patient,levels = c("P6","P5","P7","P8"))

perc_subtypes.1 <- perc_subtypes[order(perc_subtypes$Percent,decreasing=T),]

perc_subtypes$sample[perc_subtypes$Patient == "P6"] <- "S1"
perc_subtypes$sample[perc_subtypes$Patient == "P5"] <- "S2"
perc_subtypes$sample[perc_subtypes$Patient == "P7"] <- "S3"
perc_subtypes$sample[perc_subtypes$Patient == "P8"] <- "S4"

p<- ggplot(data = perc_subtypes,aes(sample,weight=Percent,fill=factor(subtypes)))+geom_bar()+theme_gray()+scale_fill_manual(values = T_colors)+
  theme(legend.title = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14))+ylab("Cellular fraction")+
  coord_cartesian(ylim = c(0,1.0)) +
  scale_x_discrete(labels = c("P6","P5","P7","P8"))

print(p)
pdf(file.path(dir_path,
              paste0('Percent.subtypes.pdf')),6,10)
print(p)
dev.off()

write.csv(perc_subtypes,file=file.path(dir_path,
                                       paste0('perc_subtypes.csv')))





#GSVA==========================================================
use.sub <- "CD8_Effector"
cluster1 <- subset(al.integrated, cells = rownames(
al.integrated@meta.data[al.integrated@meta.data$subtypes == use.sub, ]))
expr <- as.data.frame(cluster1@assays$RNA@data)

meta <- cluster1@meta.data[, "Treatment"]

expr=as.matrix(expr)
saveRDS(expr,'expr.RDS')


kegggmt <- read.gmt("/data/h.all.v7.2.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=15)


exprSet <- kegg2

group <- factor(meta,levels = c("Treatment_naive",
                                "Chemotherapy"),ordered = F)

design <- model.matrix(~group)

colnames(design) <- levels(group)

fit <- lmFit(exprSet,design)

fit2 <- eBayes(fit)

allDiff=topTable(fit2, adjust='fdr',coef=2,number=30,
                 p.value=0.05, lfc = 0.1)

allDiff=topTable(fit2, adjust='fdr',coef=2,number=30,
                 )

allDiff=topTable(fit2, adjust='fdr',coef=2,number=30,
                 p.value=0.05,lfc = 0.05)

#add pathway gene in allDiff

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
  # geom_text(aes(x = pathway, y = t, label = overlappedDEgene), 
  # 	#hjust = -0.5,
  # 	size = 2,
  # 	position = position_dodge(width = 1),
  # 	inherit.aes = TRUE)+
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
                       paste0(use.sub, '.gsva.path.pdf')), width =4, height =4)


allDiff$overlappedDEgene <- as.character(allDiff$overlappedDEgene)
write.csv(allDiff, file = file.path(dir_path, 
                                    paste0(use.sub, '.gsva.diff.csv')))

kegg3 <- as.data.frame(kegg2) 
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
                cellwidth = 0.06, cellheight = 13,
                fontsize = 10)
pdf(file.path(dir_path,  
              paste0(use.sub, ".gsva.heatmap.pdf")),width = 20,height = 6)
kk2
dev.off()

 

################### monocle
library(monocle)
CD8_obj <- subset(al.integrated, cells = rownames(
  al.integrated@meta.data[al.integrated@meta.data$subtypes %in% c("CD8_Exhausted","CD8_Effector"), ]))

expr_matrix <- as(as.matrix(CD8_obj@assays$RNA@counts), 'sparseMatrix')

p_data <- CD8_obj@meta.data 
f_data <- data.frame(gene_short_name = row.names(CD8_obj),row.names = row.names(CD8_obj))

f_data <- data.frame(gene_short_name = row.names(expr_matrix),row.names = row.names(expr_matrix))


pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

express_genes <- VariableFeatures(CD8_obj)
cds <- setOrderingFilter(cds, express_genes)



dir.create(file.path(dir_path,'T',"monocle1"))
setwd(file.path(dir_path,"T","monocle1"))

pdf('plot_ordering_gene-seurat.pdf')
plot_ordering_genes(cds)
dev.off()

diff <- differentialGeneTest(cds[ express_genes,],fullModelFormulaStr="~Treatment",cores=1) 
head(diff)

diff <- differentialGeneTest(cds[ express_genes,],fullModelFormulaStr="~subtypes",cores=1) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出1149个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
dim(deg)

write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)


ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()


cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')


cds <- orderCells(cds)
saveRDS(cds,"cds.CD8.rds")

pdf("train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()

pdf("train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="subtypes", color=T_colors,size=1,show_backbone=TRUE) +
  scale_color_manual(breaks = c("CD8_Naive","CD8_Effector", "CD8_Exhausted"), values=c("#6B8E23","#BC8F8F" , "#F4A460")) + theme(legend.position = "right")
dev.off()

pdf("train.monocle.treatment.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Treatment", color=treatment_colors,size=1,show_backbone=TRUE) +
  scale_color_manual(breaks = c(  "Treatment_naive","Chemotherapy" ), values=c("#00A087FF", "#3C5488FF" )) + theme(legend.position = "right")
dev.off()


pdf("train.monocle.patient.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Patient_ID", color=patient_colors,size=1,show_backbone=TRUE) +
  scale_color_manual(breaks = c("P5","P6","P7", "P8"), 
                     values=c("#91D1C2FF","#8491B4FF", "#F39B7FFF", "#3C5488FF"    )) + theme(legend.position = "right")
dev.off()


pdf("train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()




library(ggpubr)
df <- pData(cds) 
View(df)

pdf("Pseudotime-subtypes.pdf")
ggplot(df, aes(Pseudotime, colour = subtypes, fill=subtypes)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ scale_fill_manual(name = "", values = T_colors)+scale_color_manual(name = "", values = T_colors)
dev.off()



