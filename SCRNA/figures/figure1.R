source('~/script/selfEnvfun.R')
library(tidyverse)
library(Seurat)

projid <- 'BMLM'
ppart <- 'allCell'
mainpath <- paste0('~/project/', projid)
setwd(mainpath)

source(file = paste0(mainpath, '/script/functions/plot.R'))
source(file = paste0(mainpath, '/script/functions/singlecell.R'))

#read common data
mycols <- readRDS(file = '/path_2_dir/mycols.rds')
namedColors <- readRDS(file = '/path_2_dir/namedColors.rds')
namedColorsMainAnno <- readRDS(file = '/path_2_dir/namedColors.mainAnno.rds' )
clin <- readRDS(file = '/path_2_dir/clin.rds')
#sr.fineAnno.meta <- readRDS(file = paste0('./seurat/data/detailAna/', projid, '.', ppart, '.fineAnno.addClin.AddPM.meta.rds'))
#pm <- read.table(file = './clin/pari-meta.txt', header = T)

srobj <- qs::qread(file = '/path_2_dir/allcell.sr.qs',nthreads = 30)


#prepair dir
plotDir <- './plot/SCRNA/seurat/detailAna/allCell'
dataDir <- './data/SCRNA/seurat/detailAna/allCell'
dir_mk(plotDir)
dir_mk(dataDir)

#for mainCell type
Idents(srobj) <- 'mainCellType'
pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.mainCellType.withlegend.dimplot.pdf'), height = 4,width = 7)
DimPlot(srobj, label = F, cols = namedColorsMainAnno)
dev.off()

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.mainCellType.splitByTech.withoutlegend.dimplot.pdf'), height = 4,width = 8)
DimPlot(srobj, label = F, split.by = 'tech', cols = namedColorsMainAnno)+NoLegend()
dev.off()

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.mainCellType.splitBytissue_type.withoutlegend.dimplot.pdf'), height = 4,width = 12)
DimPlot(srobj, label = F, split.by = 'tissue_type', cols = namedColorsMainAnno)+NoLegend()
dev.off()


Idents(srobj) <- 'fineAnno'
pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.fineAnno.withlegend.dimplot.pdf'), height = 4,width = 12)
DimPlot(srobj, label = F, cols = namedColors)
dev.off()

pdf(file = paste0(plotDir, '/', projid, '.', ppart, 'groupByOrig.ident.withlegend.dimplot.pdf'), height = 4,width = 6)
DimPlot(srobj, label = F,group.by = 'orig.ident', cols = mycols)
dev.off()
pdf(file = paste0(plotDir, '/', projid, '.', ppart, 'groupBytech.withlegend.dimplot.pdf'), height = 4,width = 5)
DimPlot(srobj, label = F,group.by = 'tech', cols = mycols)
dev.off()

#tissue types
pdf(file = paste0(plotDir, '/', projid, '.', ppart, 'groupBytissue_type.withoutlegend.dimplot.pdf'), height = 4,width = 16)
DimPlot(srobj, label = F, split.by = 'tissue_type', cols = mycols) + NoLegend()
dev.off()

#marker gene
mkdf <- list(
  Cancer_cells=c('EPCAM', 'KRT18','KRT19'),
  B_cells= c('JCHAIN','MZB1'),
  `T/NK_cells`=c('CD2','CD3D','NKG7'),
  Myeloid_cells=c('CD14','FCGR3A', 'CD68','LYZ','MS4A7'),
  Mural_vascular_cells=c('RGS5', 'ACTA2'),
  Endothelial_cells=c('CLDN5', 'VWF', 'PECAM1'),
  `Astrocytes/Oligodendrocytes`=c('CRYAB','UGT8','CNP','OLIG2','FA2H'),
  `mesenchymal_cells/Fibroblasts`=c('ISLR'),
  Neurons=c('GABRA2','GABRG1')
)

mkdf <- namedList2LongDf(mk_list = mkdf,
                         mklevels =c('Myeloid_cells', 'T/NK_cells','B_cells',
                                     'Cancer_cells', 'Astrocytes/Oligodendrocytes', 'Neurons',
                                     'Endothelial_cells', 'Mural_vascular_cells', 'mesenchymal_cells/Fibroblasts'))
mkdf <- mkdf %>% dplyr::arrange(CELLTYPE)
mkdf <- mkdf[mkdf$gene %in% rownames(srobj),]

Idents(srobj) <- 'mainCellType'
pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.markergene.20250120.dotplot.pdf'), height = 4,width = 10)
plotDot2(expr_mat = srobj, markerDf = mkdf, min_expr = 0,
         minPercent = 5)
dev.off()

plist <- vector(mode = 'list', length = nrow(mkdf))
for (i in 1:length(plist)) {
  #i <- 1
  plist[[i]] <- FeaturePlot(object = srobj,
                            features = mkdf$gene[i],
                            ncol = 1, order = T, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                            cols = c('navyblue', 'orangered'),pt.size = 0.5) +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          legend.position = c(0.8,0.25))
  
}

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.markergene.noumapNOLegend.featureplot.pdf'), 
    height = 3*(floor(nrow(mkdf)/4) + 1) ,width = 12)
print(patchwork::wrap_plots(plist, ncol = 4))
dev.off()
#cell raito for main cell types.

sr.fineAnno.meta <- srobj@meta.data


findPercent(mydf = sr.fineAnno.meta,
            subInDf = "orig.ident",
            otherClassifer = c("mainCellType"),
            colorTheme = NA,
            pjid = paste0(projid, '.noColor'),
            resdir = './plot/SCRNA/seurat/MetricCompar/',
            longLegend = 10,
            ncolLegened = 6,
            pH=4,
            mappingColor = T,
            namedCol = namedColors)

#stromal types are defined as cancer cells, immune cells, and other stromal cells.
plotRatio(mydf = sr.fineAnno.meta,
          subInDf = 'tissue_type',
          otherClassifer = 'stromalType',
          Col2CountRatio ='orig.ident' ,
          usingSub = c('LM','BM'),
          min.anum = 100,
          pjid = 'BMLM.stromalType.BMLM',
          resdir = './plot/SCRNA/seurat/detailAna/allCell/cellratio',
          pH=3,pW = 2)
#immune
plotRatio(mydf = sr.fineAnno.meta[sr.fineAnno.meta$mainCellType %in% c( "Myeloid_cells", "T_cells/NK_cells", "B_cells"),],
          subInDf = 'tissue_type',
          otherClassifer = 'mainCellType',
          Col2CountRatio ='orig.ident' ,
          usingSub = c('LM','BM'),
          min.anum = 100,
          pjid = 'BMLM.immune.BMLM',
          resdir = './plot/SCRNA/seurat/detailAna/allCell/cellratio',
          pH=3,pW = 2)
