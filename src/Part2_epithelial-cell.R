
library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(snow)

source('/src/utilities.R')

dir_path <- '/data/part2'
setwd(file.path(dir_path,'Epi'))


merged.se.obj <- readRDS(file.path('/data/part1', '1.merged.se.obj.RDS'))



celltype <- c('Epi','NE')
obj <- subset(merged.se.obj, cells = 
                rownames(merged.se.obj@meta.data
                         [merged.se.obj@meta.data$all_CellType.1 %in% celltype &
                             merged.se.obj@meta.data$condition == 'Tumor' &
                             merged.se.obj@meta.data$Treatment == 'Treatment_naive', 
                         ]))





se.obj <- obj 

  se.obj<- NormalizeData(se.obj, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  se.obj <- ScaleData(se.obj)
  se.obj <- RunPCA(se.obj, features = VariableFeatures(object = se.obj))
  

  ElbowPlot(se.obj)
  ggsave(file=file.path(dir_path,'Epi-1', 'elbowplot.pdf'),width=10,height=10)
  se.obj <- FindNeighbors(se.obj, dims = 1:20)
  se.obj <- FindClusters(se.obj, resolution = 0.7)
  se.obj <- RunUMAP(se.obj, dims = 1:20)
 
use.markers <- c(
  'UCHL1', 'NCAM1', 'SYP', 'CHGA', # neuroendocine cell
  'KRT5', 'KRT6A',# Basal cell 
    'SCGB3A2', 'SCGB3A1','SCGB1A1', # club cell
  'AGER','CLIC5',# AT1 cell
  'PGC','SFTPC','SFTPD',# AT2 cell
  'FOXJ1','CCDC78'# ciliated cell
)

plot_markers <- function(x){
  FeaturePlot(se.obj, features = x,
               label=T) + NoLegend()
  ggsave(file =file.path(dir_path, 
                         'Epi-1', 'markers',paste0(x, '.pdf')),height=12,width=12)
}




dir.create(file.path(dir_path, 
                     'Epi-1', 'markers'))
for(i in c(use.markers)){
  plot_markers(i)
}

pdf(file.path(dir_path,   'Epi-1', paste0("all_CellType.1", '.umap.pdf')),8,6)
DimPlot(se.obj, group.by = 'all_CellType.1', label=T, label.size=5)
dev.off()



###########################


new.meta <- se.obj@meta.data
new.meta$cell <- rownames(new.meta)
new.meta <- new.meta %>%
  mutate(CellType = 
           #ifelse(seurat_clusters%in%c(2,17,10,12), 'Endo',
                  ifelse(seurat_clusters==19, 'Cil',
                         ifelse(seurat_clusters==25, 'AT2',
                                #ifelse(seurat_clusters%in%c(5,6,16,19,28), 'Fib',
                                       ifelse(seurat_clusters %in% c(6,12,23), 'Basal',
                                              #ifelse(seurat_clusters ==23, 'AT2',
                                                     ifelse(seurat_clusters==14, 'Club',
                                                            'NE')))))
#)))

#0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29

identical(new.meta$cell, rownames(se.obj@meta.data))
se.obj@meta.data$CellType.Epi <- new.meta$CellType

saveRDS(se.obj, 
        file = file.path(dir_path, 'Epi', 'Epi.rds'))

#######################******************plot
pdf(file.path(dir_path,'Epi', 'cluster.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, reduction ='umap', label = TRUE)
dev.off()


pdf(file.path(dir_path,'Epi', 'condition.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, reduction = 'umap', 
        group.by = 'condition', cols = condition_colors,
        ncol = 1, combine = FALSE)
dev.off()





pdf(file.path(dir_path,'Epi', 'Treatment.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, reduction = 'umap', 
        group.by = 'Treatment', cols = treatment_colors,
        ncol = 1, combine = FALSE)
dev.off()

pdf(file.path(dir_path,'Epi', 'P_ID.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, reduction = 'umap', 
        group.by = 'Patient_ID', 
        cols = patient_colors,
        ncol = 2, combine = FALSE)
dev.off()


CellType_colors <- c("#fcf0ac",  "#88B150", "#EA9236","#D0978F", "#247187")
names(CellType_colors) <- c( "Cil","Basal","AT1","NE","Club")


pdf(file.path(dir_path,'Epi', 'CellType.Epi.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, group.by = 'CellType.Epi', cols = CellType_colors,label =T,
        label.size = 6) 
dev.off()




pdf(file.path(dir_path,'Epi', 'CellType.Epi.nolabel.umap.pdf'), width = 8, height =7)
DimPlot(se.obj, group.by = 'CellType.Epi', cols = CellType_colors,label =F,
        label.size = 6) 
dev.off()


cluster_order <- c(
  0,1,3,4,7,9,11,14,15,18,22,25,26,29,
  6,12,23,
  14,
  25,
  21            
)




clusters_cell <- data.frame(Cluster= cluster_order, 
                            Celltype = as.factor(c(
                              rep('NE',14),
                              rep('Basal',3), 
                              rep('Endo',4), 
                              rep('Fib',5),
                              'Club',
                              'AT1',   
                              'AT2', 
                              'Cil'
                            )))

marker_cell <- data.frame(Markers = use.markers,
                          Celltype = as.factor(c(
                            rep('NE',4), rep('Basal',2), rep('Endo',2),
                            rep('Fib',4), 
                            rep('Club',3) ,rep('AT1', 2), rep('AT2' ,3), rep('Cil',2))))


p11 <- ggplot(marker_cell, aes(y = Markers, x = 'null')) +geom_tile(aes(fill = Celltype))+
  scale_fill_manual(values = CellType_colors) +
  scale_y_discrete(limits = rev(use.markers) ) +  
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
        legend.position = 'left', 
        legend.text=element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x =element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave(p11, file =file.path(dir_path, 'Part1', 'Marker.CellType.colors.pdf'), width = 3, height = 7)


p <- DotPlot(sclc, features = use.markers )

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(10), limits=c(-0.75,2.5))


df <- p$data
df$id <- as.character(df$id)

p1 <- ggplot(df, aes(y = features.plot,  x = as.character(id )))+
  geom_point(aes(col = avg.exp.scaled, size = pct.exp)) + sc + 
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
        legend.position = 'bottom', 
        legend.text=element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x =element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=12),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(limits = as.character(cluster_order), position = 'top') +
  scale_y_discrete(limits = rev(use.markers))

ggsave(p1, file = file.path(dir_path, 'Part1', 'use.markers.dotplot.pdf'),height=12,width=12)


p10 <- ggplot(marker_cell, aes(y = Markers, x = ' ')) +geom_tile(aes(fill = Celltype))+
  scale_fill_manual(values = CellType_colors) +
  scale_y_discrete(limits = rev(use.markers) ) +  
  theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
        legend.position = 'left', 
        legend.text=element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x =element_text(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=0),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "white"),
        axis.ticks = element_blank(),
  )


#patchwork
p101 <- p10 + p1 + plot_layout(widths = c(1, 12))
ggsave(p101, file = file.path(dir_path, 'Part1', 'Marker.CellType.dotplot.pdf'),height=8,width=12)










