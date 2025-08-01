library(Seurat)
library(ggplot2)
library(tidyverse)
library(harmony)
library(EnhancedVolcano)

projid <- 'BMLM'
ppart <- 'BBB'
mainpath <- ''

setwd(mainpath)

source('./script/functions/singlecell.R')
source('./script/functions/plot.R')

#read common data
# mycols <- readRDS(file = '/path_2_dir/mycols.rds')
# namedColors <- readRDS(file = '/path_2_dir/namedColors.rds')
# namedColorsMainAnno <- readRDS(file = '/path_2_dir/namedColors.mainAnno.rds' )
# clin <- readRDS(file = '/path_2_dir/clin.rds')
# #all meta data for all cells.
# sr.fineAnno.meta <- readRDS(file = paste0('./seurat/data/detailAna/', projid, '.', ppart, '.fineAnno.addClin.AddPM.meta.rds'))
# #pm <- read.table(file = './clin/pari-meta.txt', header = T)




#prepare dir
#wd <- './literature/result/BBB'
dataDir <- './data/literature/result/BBB'
plotDir <- './plot/literature/result/BBB'
dir_mk(c(dataDir, plotDir))

#read bbb obj.
HBVA.bbb <- readRDS(file = './data/literature/HBVA/data/HBVA.bbb.rds')
cahBm.bbb <- readRDS(file = './data/literature/Cellular_architecture_of_human_brain_metastases/data/cahBM.bbb.raw.rds')
ostromCell <- read_Da('./data/SCRNA/seurat/detailAna/otherStromalCell/BMLM.otherStromalCell.fineAnno.addClin.qs')

levels(ostromCell)
self.bbb <- subset(x = ostromCell, idents = c("Endo1", "Endo2", "Endo3", "Endo4", 
                                              "Fibro1", "Fibro2", 
                                              "Peri1", "Peri2",
                                              "SMC1", "SMC2"))

self.bbb$fineAnno <- Idents(self.bbb)
self.bbb$Cell_Type <- Idents(self.bbb)

colnames(cahBm.bbb@meta.data)
colnames(HBVA.bbb@meta.data)
colnames(self.bbb@meta.data)


HBVA.bbb.meta <- HBVA.bbb@meta.data[, c( "orig.ident", "nCount_RNA", "nFeature_RNA", "Cell_Type")]
HBVA.bbb.meta$mainCellType <- plyr::mapvalues(x = HBVA.bbb.meta$Cell_Type,
                                              from = c("Veinous", "SMC", "Pericyte",  "Capillary", "Arterial", "P. Fibro", "M. Fibro"),
                                              to = c('EC', 'SMC', 'Peri', 'EC', 'EC', 'Fibro', 'Fibro'))
HBVA.bbb.meta$dataFrom <- 'HBVA'
HBVA.bbb@meta.data <- HBVA.bbb.meta

cahBm.bbb.meta <- cahBm.bbb@meta.data[,c( "orig.ident", "nCount_RNA", "nFeature_RNA", "Cell_Type")]
cahBm.bbb.meta$mainCellType <- plyr::mapvalues(x = cahBm.bbb.meta$Cell_Type,
                                               from = c("PC-3", "PC-1", "PC-2",  "vSMCs", "EC-3", "EC-1", "EC-2"),
                                               to = c('Peri','Peri','Peri','SMC', 'EC', 'EC', 'EC'))
cahBm.bbb.meta$dataFrom <- 'cahBm'
cahBm.bbb@meta.data <- cahBm.bbb.meta

self.bbb.meta <- self.bbb@meta.data[,c( "orig.ident", "nCount_RNA", "nFeature_RNA", "Cell_Type", 'tech')]
self.bbb.meta$mainCellType <- plyr::mapvalues(x = self.bbb.meta$Cell_Type,
                                              from = c( "Endo1",  "Endo2",  "Endo3",  "Endo4",  "Fibro1", "Fibro2", "Peri1", "Peri2", "SMC1", "SMC2"),
                                              to = c('EC','EC','EC','EC','Fibro', 'Fibro', 'Peri', 'Peri', 'SMC', 'SMC'))
self.bbb.meta$dataFrom <- paste0('BMLM_', self.bbb.meta$tech)
self.bbb.meta <- self.bbb.meta[c("orig.ident", "nCount_RNA", "nFeature_RNA", "Cell_Type", "mainCellType", "dataFrom")]
head(self.bbb.meta)
self.bbb@meta.data <- self.bbb.meta


#
b1 <- CreateSeuratObject(counts = HBVA.bbb@assays$RNA@counts, meta.data = HBVA.bbb@meta.data, min.cells = 10,min.features = 200)
b2 <- CreateSeuratObject(counts =cahBm.bbb@assays$RNA@counts, meta.data = cahBm.bbb@meta.data, min.cells = 10,min.features = 200)
b3 <- CreateSeuratObject(counts = self.bbb@assays$RNA@counts, meta.data = self.bbb@meta.data, min.cells = 10,min.features = 200)

bbb <- merge(x = b1, y = c(b2, b3))

srobj <- bbb

#
######################################################
srobj[['percent.mt']] <- PercentageFeatureSet(object = srobj, pattern = '^MT-')
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srobj <- CellCycleScoring(srobj,
                          s.features = s.genes,
                          g2m.features = g2m.genes,
                          set.ident = F)

srobj <- srobj %>% SCTransform(method = "glmGamPoi",
                               vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                               verbose = FALSE,
                               return.only.var.genes = F);
srobj <- RunPCA(object = srobj,features = VariableFeatures(srobj), verbose = F)
srobj <- srobj %>%  RunHarmony("dataFrom",
                               plot_convergence = TRUE,
                               assay.use = 'SCT')
#using 95% sdv
npcs <- which(cumsum(Stdev(object = srobj,reduction = 'harmony' ))/sum(Stdev(object = srobj,reduction = 'harmony' )) > 0.95)[1]
srobj <-  RunUMAP(srobj, reduction = "harmony", dims = 1:npcs, verbose = F)
#using redued dimention component to find negihbors.
srobj <- FindNeighbors(srobj, dims = 1:npcs, reduction = 'harmony')
srobj <- FindClusters(srobj, resolution = 0.5)

srobj$dataset <- plyr::mapvalues(x = srobj$dataFrom,
                                 from = c("HBVA","cahBm", "10x", "singleron"),
                                 to = c("HBVA","cahBm", "BMLM", "BMLM"))
srobj$TvsN <- paste(srobj$mainCellType, srobj$dataset, sep = '_')



Idents(srobj) <- 'Cell_Type'
levels(srobj) <- c("Arterial", "Capillary", "Veinous", "Endo1", "Endo2", "Endo3", "Endo4", "EC-1", "EC-2", "EC-3", 
                   "M. Fibro", "P. Fibro", "Fibro1", "Fibro2",
                   "Pericyte", "Peri1", "Peri2", "PC-1", "PC-2", "PC-3",
                   "SMC", "SMC1", "SMC2", "vSMCs")
srobj$Cell_Type <- 
  factor(x = srobj$Cell_Type,
         levels = c("Arterial", "Capillary", "Veinous", "Endo1", "Endo2", "Endo3", "Endo4", "EC-1", "EC-2", "EC-3", 
                    "M. Fibro", "P. Fibro", "Fibro1", "Fibro2",
                    "Pericyte", "Peri1", "Peri2", "PC-1", "PC-2", "PC-3",
                    "SMC", "SMC1", "SMC2", "vSMCs"))

pdf(file = paste0(plotDir, '/', projid, '.bbb.CellType.labelOnData.dimplot.pdf'), width = 5, height = 5)
DimPlot(object = srobj,label = T, cols = mycols,pt.size = 1.5,raster = T) + NoLegend()
dev.off()
pdf(file = paste0(plotDir, '/', projid, '.bbb.groupBydataset.legendRight.dimplot.pdf'), width = 6, height = 5)
DimPlot(object = srobj, cols = mycols,pt.size = 1.5,raster = T, group.by = 'dataset')
dev.off()


library('AUCell')
cells_rankings <- AUCell_buildRankings(srobj@assays$RNA@data)

geneSets <- readRDS(file = '/path_2_dir/PublicData/annotation/MSigDB/MSigDB.all.geneset.v7.5.1.list.rds')
#using KEGG_TIGHT_JUNCTION, KEGG_ABC_TRANSPORTERS, REACTOME_ION_CHANNEL_TRANSPORT only.
geneSets <- geneSets[names(geneSets) %in% c('REACTOME_ION_CHANNEL_TRANSPORT',
                                            'KEGG_TIGHT_JUNCTION',
                                            'KEGG_ABC_TRANSPORTERS',
                                            'REACTOME_METAL_ION_SLC_TRANSPORTERS')]

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

geneSet <- "KEGG_TIGHT_JUNCTION"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
srobj$TJ_Score  <- aucs
geneSet <- "KEGG_ABC_TRANSPORTERS"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
srobj$ABC_Score  <- aucs
geneSet <- "REACTOME_ION_CHANNEL_TRANSPORT"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
srobj$ION_Score  <- aucs
geneSet <- "REACTOME_METAL_ION_SLC_TRANSPORTERS"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
srobj$SLC_Score  <- aucs

#transporter plot dir.
tpPlotDir <- paste0(plotDir, '/tranporter')
dir_mk(mypath = tpPlotDir)
ECs <- c("Arterial", "Capillary", "Veinous", "Endo1", "Endo2", "Endo3", "Endo4", "EC-1", "EC-2", "EC-3")
pdf(file = paste0(tpPlotDir, '/', projid, '.bbb.transporter.EC.VlnPlot.pdf'),width = 5, height = 4)
VlnPlot(object = srobj, idents = ECs,
        features = c('ION_Score','ABC_Score', 'TJ_Score','SLC_Score'),
        cols = mycols, stack = T) + NoLegend()+
  theme(axis.title.x = element_blank())
dev.off()




#self and HBVA
bbb.hs <- subset(x = srobj, subset = dataset %in% c('HBVA', 'BMLM'))
bbb.hLM <- subset(x = bbb.hs, subset = tissue_type == 'BM', invert = T)
bbb.hBM <- subset(x = bbb.hs, subset = tissue_type == 'LM', invert = T)

job::job({
  DE.EC.LMvsN <- FindMarkers(object = bbb.hLM, ident.1 =c("EC_BMLM") , ident.2 = c("EC_HBVA"), logfc.threshold = 0)
  DE.Fibro.LMvsN <- FindMarkers(object = bbb.hLM, ident.1 =c("Fibro_BMLM") , ident.2 = c("Fibro_HBVA"), logfc.threshold = 0)
  DE.Peri.LMvsN <- FindMarkers(object = bbb.hLM, ident.1 =c("Peri_BMLM") , ident.2 = c("Peri_HBVA"), logfc.threshold = 0)
  DE.SMC.LMvsN <- FindMarkers(object = bbb.hLM, ident.1 =c("SMC_BMLM") , ident.2 = c("SMC_HBVA"), logfc.threshold = 0)
}, import = c(bbb.hLM))
#BMvsH
job::job({
  DE.EC.BMvsN <- FindMarkers(object = bbb.hBM, ident.1 =c("EC_BMLM") , ident.2 = c("EC_HBVA"), logfc.threshold = 0)
  DE.Fibro.BMvsN <- FindMarkers(object = bbb.hBM, ident.1 =c("Fibro_BMLM") , ident.2 = c("Fibro_HBVA"), logfc.threshold = 0)
  DE.Peri.BMvsN <- FindMarkers(object = bbb.hBM, ident.1 =c("Peri_BMLM") , ident.2 = c("Peri_HBVA"), logfc.threshold = 0)
  DE.SMC.BMvsN <- FindMarkers(object = bbb.hBM, ident.1 =c("SMC_BMLM") , ident.2 = c("SMC_HBVA"), logfc.threshold = 0)
}, import = c(bbb.hBM))


bbb.hc <- subset(x = srobj, subset = dataset %in% c('HBVA', 'cahBm'))
#pan bm vs H
job::job({
  Idents(srobj) <- 'TvsN'
  bbb.hc <- subset(x = srobj, subset = dataset %in% c('HBVA', 'cahBm'))
  DE.EC.panvsN <- FindMarkers(object = bbb.hc, ident.1 =c("EC_cahBm") , ident.2 = c("EC_HBVA"), logfc.threshold = 0)
  DE.Peri.panvsN <- FindMarkers(object = bbb.hc, ident.1 =c("Peri_cahBm") , ident.2 = c("Peri_HBVA"), logfc.threshold = 0)
  DE.SMC.panvsN <- FindMarkers(object = bbb.hc, ident.1 =c("SMC_cahBm") , ident.2 = c("SMC_HBVA"), logfc.threshold = 0)
}, import = c(srobj))

#DEddata: rbind all data: DE.EC.LMvsN,DE.Fibro.LMvsN,DE.Peri.LMvsN,DE.SMC.LMvsN, DE.EC.BMvsN, DE.Fibro.BMvsN, DE.Peri.BMvsN, DE.SMC.BMvsN, DE.EC.panvsN, DE.Peri.panvsN, DE.SMC.panvsN
#DEdata <- do.call('rbind', DEdata)
DEdata$SYMBOL <- DEdata$gene
DEdata$threshold <- factor(
  ifelse(
    DEdata$p_val_adj < 0.01 & abs(DEdata$avg_log2FC) >= 0.5, 
    ifelse(DEdata$avg_log2FC >= 0.5 ,'Up','Down'), 'NoSignifi'
  ), 
  levels=c('Up','Down','NoSignifi'))

DEG <- DEdata[DEdata$threshold %in% c('Up','Down'),]
#
targe_tb <- read.table(file = './data/literature/result/BBB/targets_and_families.tsv', sep = '\t',skip = 1, header = T)
targe_tb <- targe_tb %>% filter(Type == 'transporter' & !(HGNC.symbol == ""))
source('~/script/selfEnvfun.R')
targe_tb <- targe_tb[targe_tb$HGNC.symbol %in% DEG$SYMBOL,]

targe_tb$tpType <- ifelse(str_detect(string = targe_tb$HGNC.name, pattern = 'ATP binding cassette'),'ABC',
                          ifelse(str_detect(string = targe_tb$HGNC.name, pattern = 'ATPase'),'ATPase',
                                 ifelse(str_detect(string = targe_tb$HGNC.name, pattern = 'solute carrier family'),'SLC','Other')))

t1 <- FetchData(object = srobj, vars = unique(targe_tb$HGNC.symbol), slot = 'data') %>%
  mutate(label = (srobj$Cell_Type)) %>%
  group_by(label) %>%
  summarise(across(everything(),mean)) %>% 
  mutate(across(-label, ~ ((.x - mean(.x))/sd(.x)))) %>%
  column_to_rownames(var = 'label') %>% t()

nrow(t1)

t1 <- t1[unlist(targe_tb$HGNC.symbol), ]
all(rownames(t1) == targe_tb$HGNC.symbol)

col_fun <- colorRamp2(c(-2, 0, 2), c("navyblue", "white", "red"))
pdf(file = paste0(tpPlotDir, '/BMLM.bbb.DEGOnly.web.heatmap.nogenename.celltypeNotClustered.pdf') ,height = 5, width = 5)
Heatmap(matrix = t1,
        col = col_fun,
        cluster_columns = F,
        row_split = unlist(targe_tb$tpType),
        show_row_names = F,
        #column_km = 2,
        show_row_dend = F,
        show_column_dend = F,
        left_annotation = rowAnnotation(type=targe_tb$tpType,
                                        col = list(type =c('ABC' = '#e66101',
                                                           'ATPase' = '#fdb863',
                                                           'SLC' = '#b2abd2',
                                                           'Other' = '#5e3c99'))),
        row_gap = unit(c(1), "mm"),
        column_gap = unit(c(1), "mm"))
dev.off()



pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.chemokine.pdf'),
    width = round(nlevels(srobj)/5) + 3, height = 4 + round(10/5))
StackedVlnPlot(obj = srobj,
               features =c('ICAM1','ICAM2','VCAM1','PECAM1','SELE','SELP','CD34',
                           'CXCL2', 'CCL2', 'CXCL12',),
               cols = namedColors,
)
dev.off()

# #ECM gene
# pdf(file = paste0(plotDir, '/',projid,'.',ppart, '.ECMGene.pdf'),
#     width = round(nlevels(srobj)/5) + 3, height = 4 + round(10/5))
# StackedVlnPlot(obj = srobj,
#                features = c('COL3A1','COL4A1',
#                             'COL4A2', 'FN1',
#                             'OGN','LUM','ECM1',
#                             'SPARC','SPARCL1',
#                             'CTSW','CTSB','HMMR'),
#                cols = namedColors,
# )
# dev.off()
# 

ECs <- c("Arterial", "Capillary", "Veinous", "Endo1", "Endo2", "Endo3", "Endo4", "EC-1", "EC-2", "EC-3")
pdf(file = paste0(plotDir, '/',projid, '.selected_ABC_SLC.vlnplot.pdf'),
    width = 4,height = 4)
StackedVlnPlot(obj = srobj,
               features = c('ABCG2', 'ABCB1','SLC3A2', 'SLC2A1', 'SLC38A5'),
               cols=mycols,idents = ECs)
dev.off()



