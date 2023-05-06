#Utilities
#dependencies===================================================================================================
suppressPackageStartupMessages({library(Seurat)
library(patchwork)
library(ggrepel)
library(ggrastr)
library(ggsci)
library(dplyr)
library(SingleR)
library(pheatmap)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(GO.db)
library(topGO)
library(clusterProfiler) 
library(GSVA)
library(limma)
library(scales)
library(ComplexHeatmap)
library(ggpubr)
library(openxlsx)
})


###############colours
immune_colors <- c("#ec7696","#b2cf87","#a35c8f","#66a9c9","#e2c027")
names(immune_colors) <- c("B","Mast","Myeloid","NK","T")

CellType_colors1 <- c("#91D1C2" ,
                      "#d22d79", 
                      "#D0978F","#806d9e")
names(CellType_colors1) <- c(  'Endo', 
                               'Fib', 
                               'NE',"Epi")


all_CellType_colors <- c(immune_colors,CellType_colors1)




#colors used===================================================================================================
treatment_colors <- c('#3C5488FF','#00A087FF')
names(treatment_colors) <- c('Chemotherapy', 'Treatment_naive')

patient_colors <- c('#e9c994',
 rev(pal_npg()(10)[4:10] ), '#c0d8e7', '#b381a5')
names(patient_colors) <- paste0('P',c(1:10))



condition_colors <- c('#a5d196', '#f27744')
names(condition_colors) <- c('Normal', 'Tumor')


CellType_colors <- c("#D43F3A" ,"#91D1C2" ,"#fcf0ac" ,
	"#d22d79", "#88B150", 
  "#EA9236"  ,"#D0978F","#247187")


names(CellType_colors) <- c( 'AT2', 'Endo', 'Cil',
 'Fib', 'Basal',
  'AT1', 'NE', 'Club')

effect_colors <- c('#FF69B4','#EECD55','#00A087FF')
names(effect_colors) <- c('Resistance','No_Responce','Treatment_naive')

ANPY_colors <- c('#FF1493','#9ACD32','#B0E0E6','#FFD700','#808A87')
names(ANPY_colors) <- c("A","N","P","Y","Undeterimated")


#endothelium subtype marker#===========================
endo.tumor <- c('INSR', 'HSPG2', 'VWA1')
endo.tip <- c('RGCC', 'RAMP3', 'ADM', 'CXCR4','PGF','LXN')
endo.stalk <- c('ACKR1', 'SELP')
endo.epc <- c('TYROBP', 'C1QB')
endo.lymp <- c('LYVE1', 'PROX1')
endo.proliferating <- c('UHRF1','TYMS','MKI67', 'CENPM')
endo.sub <- c(endo.tumor, endo.tip, endo.stalk, endo.epc, endo.lymp)

endo.stalk <- c('ACKR1', 'AQP1', 'C1QTNF9', 'CD36', 
  'CSRP2', 'EHD4', 'FBLN5', 'HSPB1', 'LIGP1', 'IL6st', 'JAM2', 'LGALS3', 
  'LRG1', 'MEOX2', 'PLSCR2', 'SDPR', 'SELP', 'SPINT2', 'TGFBI', 'TGM2', 'TMEM176A',
'TMEM176B', 'TMEM252', 'TSPAN7', 'VWF')
#=======================================================

##normal NE marker=====================================
normal.NE <- c('KCNA6','SPTSSB','KIR2DL4','CABLES1','CYTL1','NID2','NCAM1',
                'JAKMIP2-AS1', 'LGALS9B','NCR1','PPP1R9A','XCL1','TNFRSF11A',
                'DLL1','CDHR1')

#=====================================


###Fib subtype=========================================
Fib_mf <- c('ACTA2', 'ACTG2', 'MYH11', 'MYLK', 'TAGLN')

#======================================================


NE.genes <- c('BEX1', 'ASCL1', 'INSM1', 'CHGA','TAGLN3', 'KIF5C', 'CRMP1',
				'SCG3', 'SYT4', 'RTN1', 'MYT1','SYP','KIF1A', 'TMSB15A', 'SYN1','SYT11',
				'RUNDC3A', 'TFF3','CHGB','FAM57B', 'SH3GL2','BSN','SEZ6',
				'TMSB15B', 'CELF3')
nonNE.genes <- c('RAB27B', 'TGFBR2', 'SLC16A5','S100A10', 'ITGB4','YAP1',
				'LGALS3', 'EPHA2','S100A16','PLAU','ABCC3','ARHGDIB',
				'CYR61','PTGES','CCND1','IFITM2','IFITM3','AHNAK','CAV2','TACSTD2',
				'TGFBI','EMP1','CAV1','ANXA1','MYOF')



NE.sig <- c(NE.genes, nonNE.genes)


ggtheme <- theme(plot.title = element_text(hjust = 0.5,face = 'bold', size = 25),
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




add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

    new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

   repelled.y <- function(d, d.select, k = repel.degree){
                                                                                                                                                                 
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- grid::segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)


  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions


  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                   grobs = new.flag,
                                   t = 4, 
                                   l = 4
  )

  
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label


  grid::grid.newpage()
  grid::grid.draw(heatmap)


  invisible(heatmap)
}






#
