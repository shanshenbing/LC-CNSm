library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(harmony)

projid <- ''
ppart <- 'all'
mainpath <- paste0('/path_2_dir/projects/', projid, '/')

setwd(mainpath)
#source('~/script/selfEnvfun.R')

#dir mk
rawAnaPlotDir <- paste0('./seurat/plot/rawAnalysis')
rawAnaDataDir <- paste0('./seurat/data/rawAnalysis')
dir_mk(mypath = rawAnaDataDir)
dir_mk(mypath = rawAnaPlotDir)


sr.merged <- readRDS(paste0(rawAnaDataDir, '/', projid, '.raw.sr.rds'))


png(paste0('./seurat/plot/RawQC/', projid, '.', ppart, '.rawQC.png'), type="cairo", width = 1600, height = 800)
#pdf(file = paste0('./seurat/plot/RawQC/', projid, '.', ppart, '.rawQC.pdf'), width = 8, height = 6)
VlnPlot(sr.merged, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt'), ncol = 3)
dev.off()

sr.merged <- subset(sr.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

sr.merged <- CellCycleScoring(sr.merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sr.merged <- sr.merged %>% RunHarmony("tech",
                                plot_convergence = TRUE,
                                assay.use = 'SCT')

pdf(file = paste0(rawAnaPlotDir, '/',projid,'.sr.merged.harmony.elboplotByharmony.pdf'), height = 4, width = 4)
ElbowPlot(object = sr.merged, ndims = 50, reduction = 'harmony')
dev.off()

sr.merged <-  RunUMAP(sr.merged, reduction = "harmony", dims = 1:40, verbose = F)
#using redued dimention component to find negihbors.
sr.merged <- FindNeighbors(sr.merged, dims = 1:40, reduction = 'harmony')
sr.merged <- FindClusters(sr.merged, resolution = 0.8)

pdf(file = paste0(rawAnaPlotDir,'/',projid, '.sr.merged.afterRunHarmony.UmapByTech.pdf'), height = 6, width = 6)
#DimPlot(object = t1.com, reduction = "harmony", pt.size = .1, group.by = "tech",cols = 'Set1')
DimPlot(object = sr.merged, pt.size = .1, group.by = "tech",cols = mycols)
DimPlot(sr.merged, cols = mycols, group.by = 'orig.ident', label = T, repel = T)
dev.off()

pdf(file = paste0(rawAnaPlotDir,'/',projid, '.sr.merged.afterRunHarmony.cluster.splitByTech.pdf'), height = 4, width = 8)
DimPlot(sr.merged, reduction = "umap", cols = mycols,split.by = 'tech', ncol = 2, label = T, repel = T) + NoLegend()
dev.off()

pdf(file = paste0(rawAnaPlotDir,'/',projid, '.sr.merged.afterRunHarmony.cluster.splitByOrig.ident.pdf'),height = (round(length(unique(sr.merged$orig.ident))/2)*4), width = 8)
DimPlot(sr.merged, reduction = "umap", cols = mycols,split.by = 'orig.ident', ncol = 2, label = T, repel = T) + NoLegend()
dev.off()

saveRDS(object = sr.merged, file = paste0(rawAnaDataDir, '/', projid, '.sr.merged.harmony.clustered.rds'))

