library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(clustree)
library(dendextend)

# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

# source('~/script/selfEnvfun.R')
# source('~/project/functions/singlecell.R')
# source('~/project/functions/plot.R')

# # #testing
main_path <- '/path_2_dir/LMST/'
setwd(main_path)

input_dir <- './data/ST/spaceranger/SAMPLE_ID/outs/'
sample_id <- ''
projid <- 'LMST'
ppart <- 'clustering'
mainAssay <- 'Spatial.008um'
output_data_dir <- '/path_2_dir/LMST/data/ST/seurat/raw'
output_plot_dir <- '/path_2_dir/LMST/plot/ST/seurat/raw'
um8_umi_thres <- 10


plotDir <- paste0(output_plot_dir,'/',ppart,'/',sample_id)
dataDir <- paste0(output_data_dir,'/',ppart,'/',sample_id)
dir_mk(plotDir)
dir_mk(dataDir)

# 1. basic analysi --------------------------------------------------------
localdir <- input_dir
#2 um is too large.c(2,8,16)
srobj <- Load10X_Spatial(data.dir = localdir,bin.size = c(8))
#orig.srobj <- srobj
# Setting default assay changes between 8um and 16um binning
Assays(srobj)
DefaultAssay(srobj) <- mainAssay

print('raw srobj')
srobj

# 创建过滤逻辑向量
DefaultAssay(srobj) <- 'Spatial.008um'
# srobj <- subset(srobj, 
#        subset =(!is.na(nCount_Spatial.008um) & nCount_Spatial.008um >= 10) | 
#          (!is.na(nCount_Spatial.016um) & nCount_Spatial.016um >= 20))
srobj <- subset(srobj,
                subset =(!is.na(nCount_Spatial.008um) & nCount_Spatial.008um >= um8_umi_thres))
print('after basic filtering ')
srobj
print(paste0('8 um bins:',sum(!is.na(srobj$nCount_Spatial.008um))))



srobj[["percent.mt"]] <- PercentageFeatureSet(srobj, pattern = "^MT-")
srobj[['percent.rb']] <- PercentageFeatureSet(srobj, pattern = "^RP[SL]")
# HB.genes <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'BHG1','BHG2','HBM', 'HBQ1', 'HBZ')
# #if gene names not in rownames of srobj will casue error
# HB.genes <- HB.genes[HB.genes %in% rownames(srobj)]
# srobj[['percent.HB']] <- PercentageFeatureSet(srobj, features = HB.genes)

color_list <- list()

Assays(srobj)
Reductions(srobj)


DefaultAssay(srobj) <- 'Spatial.008um'
#normalize
srobj <- NormalizeData(srobj)
#无监督聚类
# note that data is already normalized
srobj <- FindVariableFeatures(srobj)
srobj <- ScaleData(srobj)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srobj <- CellCycleScoring(srobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
qs::qsave(srobj,
          file = paste0(dataDir,'/',sample_id,'.sr.qs'),nthreads = 30)
#srobj <- read_Da(paste0(dataDir,'/',sample_id,'.sr.qs'))
qs::qsave(srobj@meta.data,
          file = paste0(dataDir,'/',sample_id,'.sr.meta.qs'),nthreads = 30)
# qs::qsave(color_list,
#           file = paste0(dataDir,'/',sample_id,'.color_list.qs'),nthreads = 30)


# 2.sketch data. ----------------------------------------------------------
print('#####################################')
print('#####################################starting sketching data##################################### ')
print('#####################################')

DefaultAssay(srobj) <- "Spatial.008um"
#sketchData --- scaled.data layer.
srobj <- SketchData(
  object = srobj,
  assay = "Spatial.008um",
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch.008um",
)
# switch analysis to sketched cells
DefaultAssay(srobj) <- "sketch.008um"
names(srobj@assays$sketch@layers)
names(srobj@assays$Spatial.008um@layers)

# perform clustering workflow
srobj <- FindVariableFeatures(srobj)

var_mt_genes <- VariableFeatures(srobj)
var_rb_genes <- VariableFeatures(srobj)
var_igh_genes <- VariableFeatures(srobj)

var_mt_genes <- var_mt_genes[str_detect(var_mt_genes,"^MT-")]
var_rb_genes <- var_rb_genes[str_detect(var_rb_genes, "^RP[SL]")]
var_igh_genes <- var_igh_genes[str_detect(var_igh_genes,"^IGH")]

length(var_mt_genes)
length(var_rb_genes)
length(var_igh_genes)

srobj <- NormalizeData(srobj)



srobj <- ScaleData(srobj,features = rownames(srobj))
srobj <- RunPCA(srobj, assay = "sketch.008um", 
                reduction.name = "pca.sketch.008um",features = rownames(srobj),
                npcs = 50)
#var(srobj)


pdf(file = paste0(plotDir,'/', projid,'.pca.sketch.008um.ElbowPlot.pdf'))
ElbowPlot(object = srobj,reduction = 'pca.sketch.008um',ndims = 50)
dev.off()

srobj <- FindNeighbors(srobj, assay = "sketch.008um", reduction = "pca.sketch.008um", dims = 1:10)
srobj <- RunUMAP(srobj, reduction = "pca.sketch.008um", reduction.name = "umap.sketch.008um", return.model = T, dims = 1:10)
srobj <- FindClusters(srobj, 
                      cluster.name = paste0("seurat_cluster.sketched.008um"),
                      resolution = 1)
Idents(srobj) <- 'seurat_cluster.sketched.008um'
color_list <- list()

cluster_color <- circlize::rand_color(n = nlevels(srobj),luminosity = 'bright')
names(cluster_color) <- levels(srobj)
color_list[['sketched.008um']] <- cluster_color

#cell umap.
p1 <- DimPlot(srobj, label = T,raster = T,cols = cluster_color,pt.size = 2) +
  coord_fixed()+
  ggtitle('sketched.008um') + theme(legend.position = "bottom")
pdf(file = paste0(plotDir,'/',projid,'.','sketched.008um','.cluster.dim.pdf'),width = 6,height = 6)
print(p1)
dev.off()

pdf(file = paste0(plotDir,'/',projid,'.','sketched.008um','.cluster.qc_feature.vln.pdf'),
    width = 8,height = 5)
VlnPlot(object = srobj,
        features = c("nCount_Spatial.008um","nFeature_Spatial.008um","percent.mt", "S.Score", "G2M.Score"),stack = T,flip = T)
dev.off()

#sketch marker
mk.sr <- FindAllMarkers(srobj, assay = 'sketch.008um', only.pos = TRUE)
qs::qsave(x = mk.sr,
          file = paste0(dataDir, '/',projid,'.',sample_id,'.sketch.008um.mk.sr.qs'))
write_tsv(x = mk.sr,
          file = paste0(dataDir, '/',projid,'.',sample_id,'.sketch.008um.mk.sr.tsv'))

# 3. projection sketch ----------------------------------------------------
srobj <- ProjectData(
  object  = srobj,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch.008um",
  sketched.assay = "sketch.008um",
  sketched.reduction = "pca.sketch.008um",
  umap.model = "umap.sketch.008um",
  dims = 1:10, #1:50 in tutoril. 1:10 using all genes to runPCA.
  refdata = list(seurat_cluster.projected.008um = "seurat_cluster.sketched.008um")
)
DefaultAssay(srobj) <- "Spatial.008um"

qs::qsave(srobj,
          file = paste0(dataDir,'/',sample_id,'.sr.qs'),nthreads = 30)
qs::qsave(srobj@meta.data,
          file = paste0(dataDir,'/',sample_id,'.sr.meta.qs'),nthreads = 30)
cluster_color <- c('#2a6f96','#100466','#9edf35',
  '#a832da','#ca2e1b','#ba246c',
  '#f4ee3b','#5cade1','#4b9ca1',
  '#2f1b96','#d0361b','#ca347b','#e83ab5')
names(cluster_color) <- levels(srobj$seurat_cluster.sketched.008um)
qs::qsave(cluster_color,
          file = paste0(dataDir,'/',sample_id,'.cluster_color.qs'),nthreads = 30)

#visulaizaiotn sketched.008um adn full dataset.
DefaultAssay(srobj) <- "sketch"
Idents(srobj) <- "seurat_cluster.sketched.008um"
p1 <- DimPlot(srobj, label = T,raster = T,cols = cluster_color) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")
# switch to full dataset
DefaultAssay(srobj) <- "Spatial.008um"
Idents(srobj) <- "seurat_cluster.projected.008um"
p2 <- DimPlot(srobj, label = T,raster = T,cols = cluster_color) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")
pdf(paste0(plotDir,'/', projid,'.sketched.008um.dim.pdf'),width = 12,height = 6)
print(p1 | p2)
dev.off()

#spatial umap.
p1 <- SpatialDimPlot(srobj, label = T, repel = T,
                     group.by = 'seurat_cluster.projected.008um',
                     pt.size.factor = 5,
                     label.size = 4,cols = cluster_color)
pdf(file = paste0(plotDir,'/',projid,'.','sketched.008um','.cluster.sdim.pdf'),
    width = 7,height = 7)
print(p1)
dev.off()


pdf(file = paste0(plotDir,'/',projid,'.sketched.008um.cluster.splity-cluster.sdim.pdf'),
    width = 3*4,height = ceiling(nlevels(srobj)/4)*3)
for (i in levels(srobj)) {
  sr_sub <- subset(srobj,idents = c(i))
  print(SpatialDimPlot(sr_sub, label = T, repel = T,
                       group.by = 'seurat_cluster.projected.008um',
                       pt.size.factor = 20,
                       label.size = 4,cols = cluster_color))
}
dev.off()

# Order clusters by similarity
mk.sr <- FindAllMarkers(srobj, assay = "Spatial.008um", only.pos = TRUE)
qs::qsave(x = mk.sr,
          file = paste0(dataDir, '/',projid,'.',sample_id,'.','Spatial.008um','.mk.sr.qs'))
write_tsv(x = mk.sr,
          file = paste0(dataDir, '/',projid,'.',sample_id,'.','Spatial.008um','.mk.sr.tsv'))

qs::qsave(srobj,
          file = paste0(dataDir,'/',sample_id,'.sr.qs'),nthreads = 30)
qs::qsave(srobj@meta.data,
          file = paste0(dataDir,'/',sample_id,'.sr.meta.qs'),nthreads = 30)

