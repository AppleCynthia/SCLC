

library(Seurat)
library(ggplot2)
library(dplyr)

##################################### processing data
source('/src/utilities.R')
clus_path <- '/data/part1'

merged.se.obj <- readRDS(file.path(clus_path, 
                                   'filter.merged.se.obj'))

dir_path <- '/data/part1'
setwd(dir_path)


merged.se.obj<- NormalizeData(merged.se.obj, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
merged.se.obj <- FindVariableFeatures(merged.se.obj, 
                                      selection.method = "vst", nfeatures = 2000)
merged.se.obj <- ScaleData(merged.se.obj)
gc()
merged.se.obj <- RunPCA(merged.se.obj, 
                        features = VariableFeatures(object = merged.se.obj))
gc()
ElbowPlot(merged.se.obj)
ggsave(file=file.path(clus_path, 'elbowplot.1.pdf'))
merged.se.obj <- FindNeighbors(merged.se.obj, dims = 1:25)
merged.se.obj <- FindClusters(merged.se.obj, resolution = 0.7)
merged.se.obj <- RunUMAP(merged.se.obj, dims = 1:25)
saveRDS(merged.se.obj,'merged.se.obj.RDS')



########################### umap
DimPlot(merged.se.obj, label=T,raster=FALSE) + NoLegend()
ggsave(file.path(dir_path,
                 '1.umap.cluster.pdf'),width=10,height=10)


use.markers <- c('CD3D', 'CD3E', 'TRAC', 'CD8A', 'CD8B', # T cell 
                 
                 'GNLY', 'NKG7', # nk cell
                 
                 'CD19', 'CD79A', # B cell 
                 
                 'LYZ', 'S100A8', # monocyte cell
                 'C1QA', 'MARCO', # macrophage cell
                 'PLD4','CD1C', # monocyte-derived dendritic cell 
                 'UCHL1', 'NCAM1', 'SYP', 'CHGA', #neuroendocine cell
                 'CLDN5',  'MMRN1',# endothelial cell
                 'LUM','ACTA2', 'COL1A1','DCN', #fibroblast cell
                 'MS4A2','RGS13',# Mast cell
                 'KRT5', 'KRT6A',# Basal 11
                 'SCGB3A2', 'SCGB3A1','SCGB1A1', # club cell
                 'AGER','CLIC5',# AT1 cell
                 'PGC','SFTPC','SFTPD',# AT2 cell
                 'FOXJ1','CCDC78',# ciliated cell
                 'PTPRC'#CD45 cell
)

for(i in c(use.markers)){
  p <- FeaturePlot(merged.se.obj, features = i,raster=FALSE,
                   #cols = c('lightgrey','#725e95'), order =TRUE, 
                   label=T) + NoLegend()
  ggsave(plot=p, file =file.path(dir_path, 
                                 '1.markers',paste0(i, '.pdf')),width=10,height=10)
}

############################ label cell type
new.meta <- merged.se.obj@meta.data
new.meta <- new.meta %>% mutate(
  all_CellType.1 = ifelse(RNA_snn_res.0.7 %in% c(3), 'NK',
                          ifelse(RNA_snn_res.0.7 %in% c(1,0,24,34), 'T',
                                 ifelse(RNA_snn_res.0.7 %in% c(23,27), 'B', #22,27,35
                                        ifelse(RNA_snn_res.0.7 %in% c(5,6,26,31,35), 'Myeloid',
                                               ifelse(RNA_snn_res.0.7 %in% c(11,14,25,28,30) , 'Epi', 
                                                      ifelse(RNA_snn_res.0.7 %in% c(4,17,19)  , 'Endo',  
                                                             ifelse(RNA_snn_res.0.7 %in% c(7,12,13), 'Fib',
                                                                    ifelse(RNA_snn_res.0.7 == 36,'Mast','NE')))))))))


merged.se.obj@meta.data$all_CellType.1 <- new.meta$all_CellType.1

new <- paste0(merged.se.obj@meta.data$all_CellType.1,'-',merged.se.obj@meta.data$seurat_clusters)
merged.se.obj@meta.data$all_CellType.1.seurat <- new


pdf("1.umap.all_CellType.1.seurat.pdf",10,10)
DimPlot(merged.se.obj, group.by ='all_CellType.1.seurat',raster=FALSE) 
dev.off()
                 

pdf("1.umap.all_CellType.1.seurat.1.pdf",10,10)
DimPlot(merged.se.obj, group.by ='all_CellType.1.seurat',label=T,raster=FALSE) 
dev.off()


pdf("1.umap.lib.pdf",10,10)
DimPlot(merged.se.obj, group.by ='lib',raster=FALSE) 
dev.off()

pdf("1.umap.condition.pdf",10,10)
DimPlot(merged.se.obj, group.by ='condition',cols=condition_colors,raster=FALSE) 
dev.off()

pdf("1.umap.Treatment.pdf",10,10)
DimPlot(merged.se.obj, group.by ='Treatment',cols=treatment_colors,raster=FALSE) 
dev.off()

pdf("1.umap.Patient_ID.pdf",10,10)
DimPlot(merged.se.obj, group.by ='Patient_ID',cols=patient_colors,raster=FALSE) 
dev.off()




all_CellType.1.seurat <- c("#ec7696","#2bae85",
"#ed5a65","#7cabb1","#fb8b05",
 "#82111f","#d1c2d3","#e2d849","#fa5d19","#f1908c","#681752","#314a43",
"#ee4863",
"#bf3553",
"#ec2c64","#66c18c","#ad9e5f","#fc8c23","#5d3131",
 "#951c48","#5bae23", "#a61b29","#0f95b0","#8b2671",  "#2d0c13","#fcd337","#495c69","#3b818c","#70887d","#835e1d","#826b48","#873d24",
"#de3f7c",
"#ee2c79","#fb9968","#815c94","#1a94bc")

names(all_CellType.1.seurat) <- c("B-23", "B-27" ,      
"Endo-17", "Endo-19", "Endo-4", 
"Epi-11","Epi-14","Epi-25","Epi-28", "Epi-30","Fib-12","Fib-13" ,   
 "Fib-7",   
 "Mast-36",  
"Myeloid-26","Myeloid-31","Myeloid-35","Myeloid-5","Myeloid-6", 
"NE-10","NE-15","NE-16","NE-18","NE-2","NE-20" , "NE-21" ,"NE-22","NE-29","NE-32","NE-33","NE-8","NE-9" ,     
 "NK-3",
"T-0","T-1","T-24","T-34" )

pdf("1.umap.all_CellType.1.seurat.2.pdf",10,10)
DimPlot(merged.se.obj, group.by ='all_CellType.1.seurat',cols = all_CellType.1.seurat,label=T,raster=FALSE) 
dev.off()


##############################  the expression of marker genes in the indicated cell types

use.markers <- c('CD3D', 'CD3E', 'TRAC', 'CD8A', 'CD8B', #tcell 
                 
                 'CD19', 'CD79A', #B cell  
                 
                 'LYZ', 'S100A8', #monocyte 
                 'C1QA', 'MARCO', #macrophage 
                 'PLD4','CD1C',#with lyz monocyte-derived dendritic cell mDC， 26 吧 
                 
                 'MS4A2','RGS13',#Mast 
                 
                 'GNLY', 'NKG7', #nk  
                 
                 'CLDN5',  'MMRN1',#endo 
                 
                 'LUM','ACTA2', 'COL1A1','DCN', #fibroblast
                 
                 "EPCAM",
                 
                 'UCHL1', 'NCAM1', 'SYP', 'CHGA' #neuroendocine(NE)
              
                 
)



p<- DotPlot(merged.se.obj, features = use.markers)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")


data.anno <- data.frame(
  features.plot = unique(p$data$features.plot),
  label = c("T","T","T","T","T",
            "B","B",
            "Myeloid","Myeloid","Myeloid","Myeloid","Myeloid","Myeloid",
            "Mast", "Mast",
            "NK",  "NK", 
            "Endo","Endo",
            "Fib", "Fib", "Fib", "Fib",
            "Epi",
            "NE","NE","NE","NE"
            
  )
)

data.usage <- p$data
df.plot <- plyr::join(data.usage,data.anno)

df.plot$id.1 <- factor(df.plot$id,levels = c(1,0,24,34,#T cell
                                             23,27,  # B cell
                                             5,6,26,31,35, # Myeloid cell
                                             36, # Mast
                                             3,# NK
                                             4,17,19, # Endo cell
                                             7,12,13,  # Fib cell
                                             11,14,25,28,30,# Epi cell
                                             2,8,9,10,15,16,18,20,21,22,29,
                                             32,33 # NE cell
)
)

df.plot$features.plot.1 <- factor(df.plot$features.plot,levels =
                                    c('CD3D', 'CD3E', 'TRAC', 'CD8A', 'CD8B', # 
                                      
                                      'CD19', 'CD79A', # B cell  
                                      
                                      'LYZ', 'S100A8', # monocyte 
                                      'C1QA', 'MARCO', # macrophage 
                                      'PLD4','CD1C',#  monocyte-derived dendritic cell mDC
                                      
                                      'MS4A2','RGS13',# Mast cell
                                      
                                      'GNLY', 'NKG7', # nk cell 
                                      
                                      'CLDN5',  'MMRN1',# endothelial cell 
                                      
                                      'LUM','ACTA2', 'COL1A1','DCN', # fibroblast  cell
                                      
                                      "EPCAM", # epithelial cell
                                      
                                      'UCHL1', 'NCAM1', 'SYP', 'CHGA' # neuroendocine cell
                                    )
)


df.plot$label.1 <- factor(df.plot$label,levels=
                            c("T", "B","Myeloid","Mast", "NK", "Endo","Fib", "Epi","NE"
                            ))

df.plot$id.celltype <-   NA 

df.plot <-  df.plot %>% mutate(
  celltype = ifelse(id %in% c(3), 'NK',
                    ifelse(id %in% c(0,1,24,34), 'T',
                           ifelse(id %in% c(23,27), 'B', 
                                  ifelse(id %in% c(5,6,26,31,35), 'Myeloid',
                                         ifelse(id %in% c(11,14,25,28,30) , 'Epi', 
                                                ifelse(id %in% c(4,17,19)  , 'Endo',  
                                                       ifelse(id %in% c(7,12,13), 'Fib',
                                                              ifelse(id == 36,'Mast','NE')))))))))




df.plot$id.celltype  <- paste0(df.plot$celltype,'-c',df.plot$id)
df.plot$id.celltype <- factor(df.plot$id.celltype,levels=
                                c("T-c0",  "T-c1" , "T-c24", "T-c34" ,
                                  "B-c23","B-c27",     
                                  "Myeloid-c5",  "Myeloid-c6" ,"Myeloid-c26","Myeloid-c31" ,"Myeloid-c35",
                                  "Mast-c36",
                                  "NK-c3" ,   
                                  "Endo-c4" ,"Endo-c17", "Endo-c19", 
                                  "Fib-c7","Fib-c12","Fib-c13" , 
                  
                                  "Epi-c11" , "Epi-c14"   ,"Epi-c25" ,"Epi-c28" ,"Epi-c30",
                                  #2,8,9,10,15,16,18,20,21,22, 20220929 补上了NE-c29.
                                  "NE-c2","NE-c8","NE-c9","NE-c10","NE-c15","NE-c16","NE-c18","NE-c20", "NE-c21","NE-c22" , "NE-c29","NE-c32","NE-c33"  
                                )  )

p <- ggplot(df.plot,aes(x=features.plot.1,y =  as.numeric(id.1),size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("Percent\nExpressed", range = c(0,6)) + 
  scale_color_gradientn(colours = viridis::viridis(20),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("Markers") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id.1)),labels = levels(df.plot$id.celltype),sec.axis = dup_axis())+ #复制 y轴 代替边框效果
  facet_grid(~label.1, scales="free_x",space = "free")+theme_classic() +
  theme(
    axis.text.x = element_text(size=12, angle=90, hjust=0.5, color="black",face="bold"),#x轴标签样式
    axis.text.y = element_text(size=12, color="black",face="bold"),
    axis.title.x = element_text(size=14,colour = 'black',vjust = -0.8,hjust = 0.5),#坐标轴标题
    
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    
    panel.spacing=unit(0, "mm"), 
    strip.text.x = element_text(size=15, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),
    strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
  )

immune_colors <- c("#e2c027","#ec7696","#a35c8f","#b2cf87","#66a9c9")
names(immune_colors) <- c("T","B","Myeloid","Mast","NK")

CellType_colors1 <- c("#91D1C2" ,
                      "#d22d79", 
                      "#806d9e","#D0978F")
names(CellType_colors1) <- c(  'Endo', 
                               'NE',"Epi",'Fib')
cols <- c(immune_colors, CellType_colors1)


g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <-  cols
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
plot(g)

pdf("/data/part1/1.marker.20220922.pdf",10,10)
plot(g)
dev.off()




###################### immune cell proportion
df <- merged.se.obj@meta.data %>% mutate(type =ifelse(CellType == 'Non_immune',
                                                      'Non_immune', 'Immune cells'))

df <- merged.se.obj@meta.data
ggplot(df, aes(x = LIB, fill = immune)) +
  geom_bar(position="fill") +
  scale_fill_manual(values = c('grey', '#5691ff')) +
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
        legend.position = 'top', 
        legend.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.x =element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(y = 'Fraction of Immune/Non_immune cells')+
  coord_flip()
ggsave(file.path(dir_path,
                 '1.fractionOFimmunce.pdf'), width = 4, height = 2)


DimPlot(merged.se.obj, group.by ='condition',raster=FALSE,cols=condition_colors) 
ggsave(file.path(dir_path,
                 '1.UMAP.condation.pdf'),width=10,height=10)







###########the proportion of cell types in tumor and paired normal tissue in treated and untreated patients

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(ggalluvial)
library(ggsci)
library(cowplot)
library(tastypie)

immune_colors <- c("#ec7696","#b2cf87","#a35c8f","#66a9c9","#e2c027")
names(immune_colors) <- c("B","Mast","Myeloid","NK","T")

CellType_colors1 <- c("#91D1C2" ,
                      "#d22d79", 
                      "#D0978F","#806d9e")
names(CellType_colors1) <- c(  'Endo', 
                               'Fib', 
                               'NE',"Epi")


all_CellType_colors <- c(immune_colors,CellType_colors1)




static_cell <- function(type,Trm){
  num_subtypes = length(which(merged.se.obj@meta.data$all_CellType.1 == type & merged.se.obj@meta.data$Treatment == Trm))
  num_Treatment_ID= length(which(merged.se.obj@meta.data$Treatment  == Trm ))
  percet_subtypes = num_subtypes/num_Treatment_ID
  static_percet_subtypes = sprintf("%0.4f", percet_subtypes)
  return(static_percet_subtypes)
  
}

P_inf_Treatment <- c()
P_inf_subtypes <- c()
P_inf_percent <- c()

P_Treatment <-  c("Treatment_naive","Chemotherapy") 

for (TD in c(P_Treatment)){
  for (i in c( "NE",  "Epi" , "Endo",  "Fib","T","B" , "Myeloid", "Mast","NK"          )){
    P_inf_Treatment <- c(P_inf_Treatment,TD)
    P_inf_subtypes <- c(P_inf_subtypes,i)
    P_inf_percent <- c(P_inf_percent,static_cell(i,TD))
  }
  
}

perc_subtypes  <- data.frame(Treatment = P_inf_Treatment,
                             Percent = as.numeric((P_inf_percent)),
                             subtypes = P_inf_subtypes
)


perc_subtypes$Treatment1 <- factor(perc_subtypes$Treatment,levels=c("Treatment_naive","Chemotherapy") )
perc_subtypes$subtypes1 <- factor(perc_subtypes$subtypes,levels=c("T","B", "Myeloid", "Mast","NK",  "Endo",  "Fib", "Epi" ,"NE" ) )

p<- ggplot(data = perc_subtypes,aes(Treatment1,weight=Percent,fill=factor(subtypes1)))+geom_bar(width=0.6)+theme_gray()+scale_fill_manual(values = all_CellType_colors)+
  theme(legend.title = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14))+ylab("Cellular fraction")+coord_cartesian(ylim = c(0,1.0)) +coord_flip()


pdf(file.path(dir_path,
              paste0('Percent.treatment.all_celltypes.pdf')),10,8)
print(p)
dev.off()


g=ggplot(perc_subtypes, aes(x = Treatment1,y=Percent,fill = factor(subtypes1), 
                  stratum = subtypes1, alluvium = subtypes1)) +
  geom_col(width = 0.4,color= NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0.5)   +#knot.pos可以使连线更直
  #geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 与geom_flow效果相似
  #scale_fill_manual(values = pal_npg()(9))+  # 与细胞亚型数量对应
  #scale_fill_manual(values = pal_npg()(9))+
  scale_fill_manual(values = all_CellType_colors)+
  theme_map()+
  theme(axis.text.x=element_text(size=20,vjust = 5),
        legend.position = 'none') +
  theme_classic() +
  labs(x="Treatment",y="Ratio of cell types")
  

pdf(file.path(dir_path,
              paste0('mulberry.Percent.treatment.all_celltypes.pdf')),6,8)
print(g)
dev.off()


new.meta.data <- merged.se.obj@meta.data[which(merged.se.obj@meta.data$condition == "Normal"),]
static_cell <- function(type,Trm){
  num_subtypes = length(which(new.meta.data$all_CellType.1 == type & new.meta.data$Treatment == Trm))
  num_Treatment_ID= length(which(new.meta.data$Treatment  == Trm ))
  percet_subtypes = num_subtypes/num_Treatment_ID
  static_percet_subtypes = sprintf("%0.4f", percet_subtypes)
  return(static_percet_subtypes)
  
}

P_inf_Treatment <- c()
P_inf_subtypes <- c()
P_inf_percent <- c()

P_Treatment <-  c("Treatment_naive","Chemotherapy") 

for (TD in c(P_Treatment)){
  for (i in c( "NE",  "Epi" , "Endo",  "Fib","T","B" , "Myeloid", "Mast","NK"          )){
    P_inf_Treatment <- c(P_inf_Treatment,TD)
    P_inf_subtypes <- c(P_inf_subtypes,i)
    P_inf_percent <- c(P_inf_percent,static_cell(i,TD))
  }
  
}

perc_subtypes  <- data.frame(Treatment = P_inf_Treatment,
                             Percent = as.numeric((P_inf_percent)),
                             subtypes = P_inf_subtypes
)

perc_subtypes$Treatment1 <- factor(perc_subtypes$Treatment,levels=c("Treatment_naive","Chemotherapy") )
perc_subtypes$subtypes1 <- factor(perc_subtypes$subtypes,levels=c("T","B", "Myeloid", "Mast","NK",  "Endo",  "Fib", "Epi" ,"NE" ) )

p<- ggplot(data = perc_subtypes,aes(Treatment1,weight=Percent,fill=factor(subtypes1)))+geom_bar(width=0.6)+theme_gray()+scale_fill_manual(values = all_CellType_colors)+
  theme(legend.title = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14))+ylab("Cellular fraction")+coord_cartesian(ylim = c(0,1.0)) +coord_flip()


pdf(file.path(dir_path,
              paste0('NormalTissue_Percent.treatment.all_celltypes.pdf')),10,8)
print(p)
dev.off()

g=ggplot(perc_subtypes, aes(x = Treatment1,y=Percent,fill = factor(subtypes1), 
                            stratum = subtypes1, alluvium = subtypes1)) +
  geom_col(width = 0.4,color= NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0.5)   +#knot.pos可以使连线更直
  #geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 与geom_flow效果相似
  #scale_fill_manual(values = pal_npg()(9))+  # 与细胞亚型数量对应
  #scale_fill_manual(values = pal_npg()(9))+
  scale_fill_manual(values = all_CellType_colors)+
  theme_map()+
  theme(axis.text.x=element_text(size=20,vjust = 5),
        legend.position = 'none') +
  theme_classic() +
  labs(x="Treatment",y="Ratio of cell types")


pdf(file.path(dir_path,
              paste0('NormalTissue_mulberry.Percent.treatment.all_celltypes.pdf')),6,8)
print(g)
dev.off()


new.meta.data <- merged.se.obj@meta.data[which(merged.se.obj@meta.data$condition == "Tumor"),]
static_cell <- function(type,Trm){
  num_subtypes = length(which(new.meta.data$all_CellType.1 == type & new.meta.data$Treatment == Trm))
  num_Treatment_ID= length(which(new.meta.data$Treatment  == Trm ))
  #print(num_Effect)
  percet_subtypes = num_subtypes/num_Treatment_ID
  static_percet_subtypes = sprintf("%0.4f", percet_subtypes)
  return(static_percet_subtypes)
  
}

P_inf_Treatment <- c()
P_inf_subtypes <- c()
P_inf_percent <- c()

P_Treatment <-  c("Treatment_naive","Chemotherapy") 

for (TD in c(P_Treatment)){
  for (i in c( "NE",  "Epi" , "Endo",  "Fib","T","B" , "Myeloid", "Mast","NK"          )){
    P_inf_Treatment <- c(P_inf_Treatment,TD)
    P_inf_subtypes <- c(P_inf_subtypes,i)
    P_inf_percent <- c(P_inf_percent,static_cell(i,TD))
  }
  
}

perc_subtypes  <- data.frame(Treatment = P_inf_Treatment,
                             Percent = as.numeric((P_inf_percent)),
                             subtypes = P_inf_subtypes
)


perc_subtypes$Treatment1 <- factor(perc_subtypes$Treatment,levels=c("Treatment_naive","Chemotherapy") )
perc_subtypes$subtypes1 <- factor(perc_subtypes$subtypes,levels=c("T","B", "Myeloid", "Mast","NK",  "Endo",  "Fib", "Epi" ,"NE" ) )

p<- ggplot(data = perc_subtypes,aes(Treatment1,weight=Percent,fill=factor(subtypes1)))+geom_bar(width=0.6)+theme_gray()+scale_fill_manual(values = all_CellType_colors)+
  theme(legend.title = element_blank(),panel.grid = element_blank(),axis.ticks.x = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14))+ylab("Cellular fraction")+coord_cartesian(ylim = c(0,1.0)) +coord_flip()


pdf(file.path(dir_path,
              paste0('TumorTissue_Percent.treatment.all_celltypes.pdf')),10,8)
print(p)
dev.off()



g=ggplot(perc_subtypes, aes(x = Treatment1,y=Percent,fill = factor(subtypes1), 
                            stratum = subtypes1, alluvium = subtypes1)) +
  geom_col(width = 0.4,color= NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0.5)   +
    scale_fill_manual(values = all_CellType_colors)+
  theme_map()+
  theme(axis.text.x=element_text(size=20,vjust = 5),
        legend.position = 'none') +
  theme_classic() +
  labs(x="Treatment",y="Ratio of cell types")


pdf(file.path(dir_path,
              paste0('TumorTissue_mulberry.Percent.treatment.all_celltypes.pdf')),6,8)
print(g)
dev.off()
