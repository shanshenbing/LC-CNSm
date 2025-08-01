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
ppart <- 'TNKCell'
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

srobj <- qs::qread(file = '/path_2_dir/TNKCell.sr.qs',nthreads = 30)

#prepare dir
dataDir <- paste0('./data/SCRNA/seurat/detailAna/', ppart)
plotDir <- paste0('./plot/SCRNA/seurat/detailAna/', ppart)
dir_mk(dataDir)
dir_mk(plotDir)

srobjCols <- namedColors[names(namedColors) %in% levels(srobj)]


CD4.8CluVec <- c("CD4T_CXCL13", "CD4Tn", "CD4Tprf", "CD4Treg1", "CD4Treg2",
                 "CD4T", "CD4Tcm", "CD8Tef_GZMBL", "CD8Tem_CCL", "CD8Tem_GZMK",
                 "CD8Tef_GZMAH", "CD8Tef_GZMAL", "CD8Tprf", "Tab" )



#part1
# 1.  metric compare---------------------------------------------------------------------
pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.dimplot.pdf'),
    width = 4,
    height =4)
DimPlot(object = srobj,label = T, cols = srobjCols, raster = T) + NoLegend()
dev.off()

pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.splitBytissue_type.dimplot.pdf'),
    width = 12,
    height =4)
DimPlot(object = srobj,
        cols = srobjCols,
        split.by = 'tissue_type', 
        label = T,raster = T) + NoLegend()

dev.off()


# 2. scores ---------------------------------------------------------------
cells_rankings <- AUCell_buildRankings(srobj@assays$RNA@data)

#using Cellular_architecture_of_human_brain_metastases
cahBM.sig <- readxl::read_xlsx(path = './data/literature/Cellular_architecture_of_human_brain_metastases/data/marker_signature_genes.xlsx', sheet = 2,skip = 1)
head(cahBM.sig)

tem <- c("Pro_inflammatory","T Cell Exhaustion")

plotScore(resDataDir = paste0(dataDir, '/signatureScore'),
          resPlotDir = paste0(plotDir, '/signatureScore'),
          filePfx = paste0(projid, '.', ppart),
          fileSfx = 'cahBM.unsorted',
          sortPoints = F,
          pointsize = 2,
          legendPos = c(0, 0.7),
          pH=8, pW = 8)

#for gene promoting T cell proliferation.
using_gene <- c('PDCD1', 'LAG3', 'TIGIT','CTLA4', 
                'HAVCR2', 'LAYN')
sr.MeanExpr <- 
  AverageExpression(object = srobj, 
                    features = using_gene,
                    slot = 'data', assays = 'SCT')$SCT
sr.MeanExpr <- sr.MeanExpr %>% as.data.frame()

pdf(file = paste0(plotDir, '/', projid, '.',ppart, '.slected_EX_gene.meanExpr.manual_sorted.heatmap.pdf'),
    width = 5, height = 5)
ComplexHeatmap::pheatmap(mat =as.matrix(sr.MeanExpr),
                         scale = 'row',
                         cluster_rows = F, cluster_cols = F,
                         color = circlize::colorRamp2(c(-2.5, 0, 2.5), c("navy", "white", "firebrick3")),
                         row_title= NULL,
                         column_order=c(
                           "NK", "Tgd", "Tab", "CD4Tn", "CD4T", "CD4Tcm", "CD4Tprf",
                           "CD4T_CXCL13", "CD4Treg1",  "CD4Treg2", 
                           "CD8Tprf", "CD8Tem_GZMK", "CD8Tem_CCL",
                           "CD8Tef_GZMAL", "CD8Tef_GZMAH", "CD8Tef_GZMBL"
                         )
)
dev.off()


cytotoxicSig <- c('GZMA', 'GZMB','GNLY', 'PRF1', 'IFNG', 'NKG7', 'GZMK',
                  'TNFSF10','FASLG')
pdf(file = paste0(plotDir, '/', projid, '.',ppart, '.cytotoxicSig.vlnplot.pdf'),width = 7,height = 10)
StackedVlnPlot(srobj,
               cytotoxicSig,
               pt.size=0, cols=namedColors)
dev.off()



# 3. cell ratio -----------------------------------------------------------
t1 <- sr.fineAnno.meta[
  sr.fineAnno.meta$mainCellType %in% c('Myeloid_cells', 'T_cells/NK_cells', 'B_cells'),]
plotRatio(mydf = t1,
          subInDf = 'tissue_type',
          otherClassifer = 'fineAnno',
          Col2CountRatio = 'orig.ident',
          usingSub = c('BM', 'LM'),
          min.anum = 100,
          pjid = 'BMLM.fineAnno_all_immune',numCol = 8,
          resdir = './plot/cellratio',
          pH=3,pW = 2)
plotRatio(mydf = t1,
          subInDf = 'PM',
          otherClassifer = 'fineAnno',
          Col2CountRatio = 'orig.ident',
          usingSub = c('pri', 'meta'),
          min.anum = 100,
          colorTheme = scale_color_jco,
          fillTheme = scale_fill_jco,
          pjid = 'BMLM.fineAnno_all_immune.PM',numCol = 8,
          resdir = './plot/cellratio',
          pH=3,pW = 2)



# 4. pseudotime ----------------------------------------------------------------
srobj <- subset(srobj, fineAnno %in% c("CD8Tef_GZMBL", "CD8Tem_CCL", "CD8Tem_GZMK",
                                       "CD8Tef_GZMAH", "CD8Tef_GZMAL",
                                       "CD8Tprf"))
#srobjCols <- namedColors[names(namedColors) %in% levels(srobj)]
srobjCols <- namedColors[names(namedColors) %in% as.character(unique(srobj$fineAnno))]
#Store data in a cell_data_set object
exprDa <- GetAssayData(object = srobj,slot = 'counts', assay = 'RNA')
metaDa <- srobj@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(exprDa))
rownames(gene_annotation)<-rownames(exprDa)

cds <- new_cell_data_set(exprDa,
                         cell_metadata = metaDa,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds <- align_cds(cds, alignment_group = batchCol)
cds <- reduce_dimension(cds)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
cds <- choose_cells(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

qs::qsave(x = cds,nthreads = 10,
          file = paste0(dataDir, '/', projid, '.', ppart, '.cds.qs'))
all(names(pseudotime(cds)) == rownames(colData(cds)))
colData(cds)$pseudotime <- pseudotime(cds)
qs::qsave(x = colData(cds),nthreads = 10,
          file = paste0(dataDir, '/', projid, '.', ppart, '.cds.colData.qs'))

pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.graph.dimplot.pdf'),
    width = 4,
    height = 4)
plot_cells(cds,
           color_cells_by = "fineAnno",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.5,
           group_label_size = 5,
           rasterize = T) + scale_color_manual(values = srobjCols)
dev.off()
pdf(file = paste0(plotDir, '/', projid, '.', ppart, '.pesudotime.featureplot.pdf'),
    width = 5,
    height = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           graph_label_size=2,
           cell_size = 0.5,
           rasterize = T)
dev.off()