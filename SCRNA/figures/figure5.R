library(tidyverse)
library(RColorBrewer)
library(Seurat)
#library(pheatmap)
library(patchwork)
library(ggpubr)
#library(harmony)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(AUCell)
# library(SeuratData)
# library(SeuratDisk)




projid <- 'BMLM'
ppart <- 'cancerCell'
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

srobj <- qs::qread(file = '/path_2_dir/cancerCell.sr.qs',nthreads = 30)


#prepare dir
dataDir <- paste0('./data/SCRNA/seurat/detailAna/', ppart)
plotDir <- paste0('./plot/SCRNA/seurat/detailAna/', ppart)
dir_mk(mypath = c(plotDir, dataDir))


# 1. basic analysis -------------------------------------------------------
pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.splitBytissue_type.dimplot.pdf'),
    width = 12,height =4)
DimPlot(object = srobj,
        cols = srobjCols,
        split.by = 'tissue_type', 
        label = T,raster = T) + NoLegend()
dev.off()

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.labeled.nolegend.dimplot.pdf'),
    width = 4,height = 4)
DimPlot(object = srobj,label = T,cols = srobjCols, raster = T)+NoLegend()
dev.off()

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.groupByorig.ident.labeled.nolegend.dimplot.pdf'),
    width = 4,height = 4)
DimPlot(object = srobj,label = T,cols = mycols, group.by = 'orig.ident', raster = T) +NoLegend()
dev.off()

#for replot
pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.replotMarkerGene.featureplot.GreyBG.using_orig.pdf'),
    width = 6*3, height =6*3)
p1 <- FeaturePlot(object = srobj,
                  features =c('EPCAM', 'KRT18','SOX2','MYC', 'STAT3', 'EGFR'),
                  order = T,cols = c(c('lightblue','white', 'orangered')),
                  min.cutoff = 'q10',max.cutoff = 'q90',raster = F)
rasterize(p1, layers='Point', dpi=300)
dev.off()
