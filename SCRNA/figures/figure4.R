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
ppart <- 'myeloidCell'
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


#functions
#prepare dir
dataDir <- paste0('./data/SCRNA/seurat/detailAna/', ppart)
plotDir <- paste0('./plot/SCRNA/seurat/detailAna/', ppart)
dir_mk(dataDir)
dir_mk(plotDir)

srobj <- read_Da(yourpath = '/path_2_dir/myeloidCell.sr.qs')
srobjCols <- namedColors[names(namedColors) %in% levels(srobj)]
# 1. basic metric  ----------------------------------------
pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.dimplot.pdf'),
    width = 4,
    height =4)
DimPlot(object = srobj,label = T, cols =srobjCols , raster = T) + NoLegend()
dev.off()

pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.splitBytissue_type.dimplot.pdf'),
    width = 12,
    height =4)
DimPlot(object = srobj,
        cols = srobjCols,
        split.by = 'tissue_type', 
        label = T,raster = T) + NoLegend()

dev.off()


pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.Macro_CXCL9.featureplot.pdf'),
    width = 3*2, height =3*2)
FeaturePlotNOaxis(yourobj = srobj,
                  genevec = c('CXCL9', 'CXCL10', 'CXCL11'),
                  pt.size=0.5,ncol = 2,
                  colVec =c('#abd3e0','white', 'red') ,
                  legendPos = c(0,0.3))
dev.off()







# 2. integration analysis -------------------------------------------------
#cell rank.
immR <- read_Da('./path_2_dir/Composition.Fil_min.anum100after_reshaped_arranged.tsv')
colnames(immR) <- c('tissue_type', 'orig.ident','fineAnno', 'ratio')
n_distinct(immR$orig.ident)
head(immR)

tgc <- Rmisc::summarySE(immR,measurevar="ratio", groupvars=c("tissue_type","fineAnno"))
tgc <- tgc %>% dplyr::group_by(tissue_type) %>% 
  dplyr::mutate(rankid=rank(-ratio))
tgc <- tgc %>% dplyr::arrange(tissue_type, desc(ratio))

view(tgc)

color_df <-as.data.frame(namedColors)
color_df <- color_df  %>%  
  rownames_to_column(var = 'fineAnno') %>% 
  dplyr::rename(point_col=namedColors)
tgc <- left_join(x = tgc,y = color_df)

pdf(file = paste0(plotDir, '/', projid, '.immune.rank.pdf'),
    width = 5,height = 3.5)
ggplot(tgc, aes(x=rankid, y=ratio,group=tissue_type)) + 
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), 
                width=.1, alpha= 0.3, color='grey') +
  geom_line(aes(color=tissue_type),size=0.5) +
  geom_point(size=2,color=tgc$point_col)+
  geom_label_repel(data = tgc %>% dplyr::group_by(tissue_type) %>% slice_head(n=3),
                   mapping = aes(label=fineAnno), 
                   max.overlaps = 10, color='#2b8cbe', fill = "white") +
  scale_color_lancet()+
  ylab('Ratio in immune cells')+xlab('Cell rank')+
  theme(legend.position = c(0.8,0.8))
dev.off()

t1 <- tgc
a1 <- distinct(sr.fineAnno.meta[, c('fineAnno', 'mainCellType')])
t1 <- left_join(x = t1, y = a1)
#sort as cell types.
t1 <- t1 %>% dplyr::arrange(tissue_type, fineAnno)
t3 <- reshape2::dcast(t1[,c('tissue_type','rankid', 'fineAnno', 'mainCellType')],
                      mainCellType+fineAnno ~ tissue_type, value.var = 'rankid')
#all cell.
temp_cor <- cor(t3[,c('BM', 'LM', 'primaryLungCancer')],method = "kendall")
col_fun <- circlize::colorRamp2(c(0.1, 0.5),c('#ffffb2', '#e31a1c'))

pdf(file = paste0(plotDir, '/', projid, '.TME_similarity_by_rank.all_imm.pdf'),
    width = 4,height = 4)
corrplot::corrplot(temp_cor, order = 'AOE', addCoef.col = 'black', 
                   tl.pos = 'd',
                   #cl.pos = 'n', 
                   col = rev(corrplot::COL2('PiYG')))
dev.off()
#TNK
temp_cor <- cor(t3[t3$mainCellType %in% c('T_cells/NK_cells','B_cells') ,c('BM', 'LM', 'primaryLungCancer')],method = "kendall")
col_fun <- circlize::colorRamp2(c(0.1, 0.7),c('#ffffb2', '#e31a1c'))
pdf(file = paste0(plotDir, '/', projid, '.TME_similarity_by_rank.lyphocytes.pdf'),
    width = 4,height = 4)
corrplot::corrplot(temp_cor, order = 'AOE', addCoef.col = 'black', 
                   tl.pos = 'd',
                   col = rev(corrplot::COL2('PiYG')))
dev.off()

#myeloid
temp_cor <- cor(t3[t3$mainCellType == 'Myeloid_cells',c('BM', 'LM', 'primaryLungCancer')],method = "kendall")
pdf(file = paste0(plotDir, '/', projid, '.TME_similarity_by_rank.myeloid_cells.pdf'),
    width = 4,height = 4)
corrplot::corrplot(temp_cor, order = 'AOE', addCoef.col = 'black', 
                   tl.pos = 'd',
                   #cl.pos = 'n', 
                   col = rev(corrplot::COL2('PiYG')))
dev.off()



immR <- read_Da('./path_2_dir/Composition.Fil_min.anum100after_reshaped_arranged.tsv')
colnames(immR) <- c('tissue_type', 'orig.ident','fineAnno', 'ratio')
n_distinct(immR$orig.ident)
head(immR)

allR <- immR

#cell correlation.
allR <- reshape2::dcast(data = allR, orig.ident~fineAnno,value.var = 'ratio',fill=0)

allR <- column_to_rownames(.data = allR,var = 'orig.ident')
color_df <-as.data.frame(namedColors)
color_df <- color_df  %>% rownames_to_column(var = 'fineAnno') %>% rename(point_col=namedColors)

#all col-----CXCL9.
using_vec <- c('Macro_SPP1/CCL3', 'Macro_SPP1', 'CD8Tprf', 'CD4T_CXCL13','Macro_CXCL9',
               'BCells1', 'pDC', 'cDC2')
rcc.cxcl9  <- allRpcc[rownames(allRpcc) %in% c( 'Macro_CXCL9'), ]
rcc.cxcl9  <- reshape2::melt(rcc.cxcl9)
colnames(rcc.cxcl9) <- c('PCC')
rcc.cxcl9$fineAnno  <- rownames(rcc.cxcl9) 
rcc.cxcl9 <- left_join(rcc.cxcl9, color_df, 'fineAnno')
rcc.cxcl9  <- rcc.cxcl9  %>% mutate(fineAnno=fct_reorder(fineAnno, PCC, .desc=T),
                                    lab_x_pos=ifelse(PCC > 0, 1.1, -0.1),
                                    lab_x_corlor= ifelse(fineAnno %in% using_vec, '#2b8cbe', 'grey25'),
                                    line_corlor=ifelse(fineAnno %in% using_vec, 'blue', 'grey'),
                                    point_col=ifelse(fineAnno %in% using_vec, point_col, 'grey'))
#remove self.
rcc.cxcl9 <- rcc.cxcl9[!(rcc.cxcl9$fineAnno %in% 'Macro_CXCL9'), ]


#all col-----CXCL9.
pdf(file = paste0( plotDir, '/',projid,'.',ppart, '.Macro_CXCL9.cor.pdf'),
    width = 7, height = 7)
ggplot(rcc.cxcl9, aes(x = fineAnno, y = PCC)) +
  #geom_segment(aes(xend = fineAnno, yend = 0), color = rcc.cxcl9$line_corlor, alpha=0.5) +  # 画线段
  geom_text(aes(y=0,label=fineAnno),
            hjust=rcc.cxcl9$lab_x_pos,
            color=rcc.cxcl9$lab_x_corlor,
            angle=90)+
  geom_point(color=rcc.cxcl9$point_col,size = 4)+
  geom_hline(mapping = aes(yintercept=0))+
  ylab('PCC')+ ggtitle('Macro_CXCL9')+
  scale_y_continuous(limits = c(-0.4, 0.8), n.breaks = c(6))+
  #scale_color_manual(values=namedColors[names(namedColors) %in% rcc.cxcl9$fineAnno])+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(arrow = 
                                     arrow(length = unit(0.3, 'cm'),
                                           ends = 'both',type = 'closed')),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position='none')
dev.off()



#all col-----CXCL13.
using_vec <- c('Macro_SPP1/CCL3', 'Macro_SPP1', 'CD4T_CXCL13','Macro_CXCL9',
               'PlasmaCells', 'CD8Tem_CCL', 'CD8Tef_GZMAL', 'BCells1')
rcc.cxcl9  <- allRpcc[rownames(allRpcc) %in% c( 'CD4T_CXCL13'), ]
rcc.cxcl9  <- reshape2::melt(rcc.cxcl9)
colnames(rcc.cxcl9) <- c('PCC')
rcc.cxcl9$fineAnno  <- rownames(rcc.cxcl9) 
rcc.cxcl9 <- left_join(rcc.cxcl9, color_df, 'fineAnno')
rcc.cxcl9  <- rcc.cxcl9  %>% mutate(fineAnno=fct_reorder(fineAnno, PCC, .desc=T),
                                    lab_x_pos=ifelse(PCC > 0, 1.1, -0.1),
                                    lab_x_corlor= ifelse(fineAnno %in% using_vec, '#2b8cbe', 'grey25'),
                                    line_corlor=ifelse(fineAnno %in% using_vec, 'blue', 'grey'),
                                    point_col=ifelse(fineAnno %in% using_vec, point_col, 'grey'))
#remove self.
rcc.cxcl9 <- rcc.cxcl9[!(rcc.cxcl9$fineAnno %in% 'CD4T_CXCL13'), ]


#all col-----CXCL13.
pdf(file = paste0( plotDir, '/',projid,'.',ppart, '.CD4T_CXCL13.cor.pdf'),
    width = 7, height = 7)
ggplot(rcc.cxcl9, aes(x = fineAnno, y = PCC)) +
  #geom_segment(aes(xend = fineAnno, yend = 0), color = rcc.cxcl9$line_corlor, alpha=0.5) +  # 画线段
  geom_text(aes(y=0,label=fineAnno),
            hjust=rcc.cxcl9$lab_x_pos,
            color=rcc.cxcl9$lab_x_corlor,
            angle=90)+
  geom_point(color=rcc.cxcl9$point_col,size = 4)+
  geom_hline(mapping = aes(yintercept=0))+
  ylab('PCC')+ggtitle('CD4T_CXCL13')+
  scale_y_continuous(limits = c(-0.4, 0.8), n.breaks = c(6))+
  #scale_color_manual(values=namedColors[names(namedColors) %in% rcc.cxcl9$fineAnno])+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(arrow = 
                                     arrow(length = unit(0.3, 'cm'),
                                           ends = 'both',type = 'closed')),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position='none')
dev.off()
