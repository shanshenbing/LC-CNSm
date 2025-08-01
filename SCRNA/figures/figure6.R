library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(AUCell)



projid <- 'BMLM'
ppart <- 'otherStromalCell'
mainpath <- ''

setwd(mainpath)

source('./script/functions/plot.R')
source('./script/functions/singlecell.R')

#read common data
mycols <- readRDS(file = '/path_2_dir/mycols.rds')
namedColors <- readRDS(file = '/path_2_dir/namedColors.rds')
namedColorsMainAnno <- readRDS(file = '/path_2_dir/namedColors.mainAnno.rds' )
clin <- readRDS(file = '/path_2_dir/clin.rds')
#all meta data for all cells.
sr.fineAnno.meta <- readRDS(file = paste0('./seurat/data/detailAna/', projid, '.', ppart, '.fineAnno.addClin.AddPM.meta.rds'))
#pm <- read.table(file = './clin/pari-meta.txt', header = T)

srobj <- qs::qread(file = '/path_2_dir/otherStromalCell.sr.qs',nthreads = 30)


# 1.basic analysis --------------------------------------------------------
#part1 mapping colors of dimplot.
pdf(file = paste0(plotDir, '/', projid, '.',ppart,'.fineAnno.mappingColored.dimplot.pdf'),
    width = 4,height = 4)
DimPlot(object = srobj, cols =srobjCols, label = T, raster = T) + NoLegend()
dev.off()

plotRatio(mydf = sr.fineAnno.meta[sr.fineAnno.meta$mainCellType %in% 
                                    c("Endothelial_cells", 
                                      "Mural_vascular_cells",
                                      "Mesenchymal_cells/Fibroblasts"),],
          subInDf = 'tissue_type',
          otherClassifer = 'fineAnno',
          Col2CountRatio ='orig.ident' ,
          usingSub = c('LM','BM'),
          min.anum = 100,
          pjid = 'otherStromalCell',
          resdir = './plot/SCRNA/seurat/detailAna/otherStromalCell/cellratio',
          pH=3,pW = 2)

# 2. endo -----------------------------------------------------------------
srobj.mk.endo1 <- FindMarkers(srobj,
                              ident.1 = 'Endo1',
                              ident.2 = c('Endo2'),
                              only.pos = F,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
srobj.mk.endo1$SYMBOL <- rownames(srobj.mk.endo1)
#for using Findmarkers only.
srobj.mk.endo1$gene <- rownames(srobj.mk.endo1)
#srobj.mk.endo1$compar <- 'Endo1vsEndo2-3'
srobj.mk.endo1$compar <- 'Endo1vsEndo2'
srobj.mk.endo1$threshold <- factor(
  ifelse(
    srobj.mk.endo1$p_val_adj < 0.01 & abs(srobj.mk.endo1$avg_log2FC) >= 0.5, 
    ifelse(srobj.mk.endo1$avg_log2FC >= 0.5 ,'Up','Down'), 'NoSignifi'
  ), 
  levels=c('Up','Down','NoSignifi'))


srobj.cancer <- qs::qread('cancerCell.sr.qs')
srobj.meta.cancer <- srobj.cancer@meta.data
#using mki67 directly.
mkiratio <- FetchData(object = srobj.cancer, vars = 'MKI67')
mkiratio$cellid <- rownames(mkiratio) 
mkiratio <- left_join(x = mkiratio,y = srobj.meta.cancer, by='cellid')
mkiratio <- mkiratio %>% select(MKI67, orig.ident) %>% 
  group_by(orig.ident) %>% summarise(ratio=sum(MKI67 > 0)/length(MKI67))

mkiratio <- left_join(x = mkiratio,
                      y = clin[,c('sample_id', 'PM', 'tissue_type')],
                      by=c(orig.ident='sample_id'))

pdf(file = paste0(plotDir, '/', projid, '.MKI67.percent.loliplot.pdf'),
    width = 6,height = 3)
t1 <- mkiratio
t1$ratio <- signif(x = t1$ratio,digits = 2)
ggplot(t1, aes(x=orig.ident, y=ratio)) +
  geom_bar(stat = 'identity', width = 0.2, fill='#2b8cbe', alpha = 0.5) +
  geom_point(size=10, color='#2b8cbe') +
  geom_text(aes(label= ratio), vjust=0.8, color="white", size=3)+
  ylab('MKI67 Percent')+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank())
dev.off()
#cellratio vs cancer cell mki67 percent.
allR <- qs::qread('/BMLM.allCell.percent.qs')
#allR <- reshape2::dcast(data = allR, orig.ident~fineAnno,value.var = 'ratio',fill=0)
#allR <- column_to_rownames(.data = allR,var = 'orig.ident')
head(allR)
head(mkiratio)
mkiratio <- left_join(x = mkiratio,y = allR, by='orig.ident')
t1 <- reshape2::melt(data = mkiratio,
                     id.vars=c("orig.ident", "ratio", "PM", "tissue_type"),
                     value.name = 'cellPercent',
                     variable.name= 'cellType')
head(t1)
#endo ratio vs mki67 percent.
endoR <- sr.fineAnno.meta %>% 
  filter(fineAnno %in% c("Endo1", "Endo2", "Endo3", "Endo4")) %>% 
  dplyr::select(fineAnno, orig.ident) %>% 
  dplyr::group_by(fineAnno,orig.ident) %>% 
  dplyr::summarise(num=n()) %>% 
  dplyr::group_by(orig.ident) %>% 
  dplyr::mutate(ratio=num/sum(num)*100,
                anum=sum(num)) %>% 
  dplyr::filter(anum > 10)
endoR <- reshape2::dcast(data = endoR, orig.ident ~ fineAnno, value.var = 'ratio', fill=0 )
endoR <- reshape2::melt(data = endoR,
                        id.vars = c('orig.ident'),
                        value.name = 'ratio')
endoR <- dplyr::rename(.data = endoR, fineAnno='variable', cellPercent='ratio')

endoR <- left_join(x = endoR,y = mkiratio[,1:4], by='orig.ident')

pdf(file = paste0(plotDir, '/', projid, 
                  '.MKI67_percent-endoratioPercent.point.resize.pdf'),
    width = 12,height = 3)
ggplot(data = endoR,
       mapping = aes(x =cellPercent , y = ratio))+
  geom_point(aes(color=sequencedSampleType, size=2)) +
  geom_smooth(method = 'lm', se=F, color='black')+
  scale_color_aaas()+
  facet_wrap(facets = vars(fineAnno),ncol = 4)+
  stat_cor(label.y = 0.2,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 0.25)+
  xlab('Cell ratio') + ylab('MKI67 Percent') +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.5), face = 'bold'))

dev.off()




