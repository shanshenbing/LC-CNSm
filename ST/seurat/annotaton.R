library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scop)

# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

# source('~/script/selfEnvfun.R')
# source('~/project/functions/singlecell.R')
# source('~/project/functions/plot.R')



# # #testing
main_path <- '~/project/LMST/'
setwd(main_path)
# LMK05 011184-1 020021-1b
projid <- 'LMST'
ppart <- 'slef_clustering_anno'
sample_id <- 'your_sample_id'
mainAssay <- 'Spatial.008um'
input_sr_path <- './path_2_dir/your_sample_id.sr.qs'

plotDir <- paste0('./plot/ST/seurat/analysis','/',ppart,'/',sample_id)
dataDir <- paste0('./data/ST/seurat/analysis','/',ppart,'/',sample_id)
dir_mk(plotDir)
dir_mk(dataDir)

srobj <- read_Da(input_sr_path)
Reductions(srobj)
Assays(srobj)
DefaultAssay(srobj)
# 1. anno cluster---------------------------------------------------------------
sum(is.na(srobj$seurat_cluster.projected.008um))

new.cluster.ids <- c("0" ="Epi1-ADH1C", "1" ="Epi2-ELF3", "2" ="Fib1-LAMA4", "3" ="Fib2-LF", "4" ="Epi3-VEGFA", "5" ="Fib3-CTHRC1", "6" ="VF1-MCAM", "7" ="B1-IGHG1", "8" ="B2-IGHG3", "9" ="VF2-PLVAP", "10" ="Macro-CXCL9", "11" ="Macro-SPP1")
srobj <- RenameIdents(srobj, new.cluster.ids)
srobj$self_cluser_anno <- Idents(srobj)
#colors
self_cluster_color <- qs::qread('/path_2_dir/cluster_color.qs')

pdf(file = paste0(plotDir,'/', sample_id,'self_cluster.dim.pdf'),width = 5,height = 4)
DimPlot(srobj, label = TRUE, pt.size = 0.5,
        cols = self_cluster_color)+bsaxis()+saxis(aspect.ratio=1)
DimPlot(srobj, label = TRUE, 
        cols = self_cluster_color)+bsaxis()+saxis(aspect.ratio=1)

dev.off()
pdf(file = paste0(plotDir,'/', sample_id,'self_cluster.split.dim.pdf'),
    width = 3*3+1,height = ceiling(nlevels(srobj)/3)*3)
DimPlot(srobj,  
        split.by = 'self_cluser_anno',ncol = 3,
        cols = self_cluster_color)+bsaxis()+saxis(aspect.ratio=1)
dev.off()


srobj$self_cluser_anno <- factor(x =srobj$self_cluser_anno,
                                 levels = c("Epi1-ADH1C", "Epi2-ELF3", "Epi3-VEGFA",
                                            "Fib1-LAMA4", "Fib2-LF", "Fib3-CTHRC1",
                                            "B1-IGHG1","B2-IGHG3",
                                            "Macro-CXCL9", "Macro-SPP1",
                                            "VF1-MCAM","VF2-PLVAP"))
Idents(srobj) <- 'self_cluser_anno'
qs::qsave(x = srobj,
          file = paste0(dataDir,'/', projid,'.',sample_id,'.annoed.sr.qs'),nthreads = 30)
qs::qsave(x = srobj@meta.data,
          file = paste0(dataDir,'/', projid,'.',sample_id,'.annoed.sr.meta.qs'),nthreads = 30)
#export cell groups for cloupe.
tempobj <- srobj@meta.data[,'self_cluser_anno',drop=F]
tempobj <- tempobj %>% rownames_to_column(var = 'Barcode') %>% rename(Cluster=self_cluser_anno) %>% 
  relocate(Barcode,.before = Cluster)
write_csv(x = tempobj,
          file = paste0(dataDir, '/', projid,'.',sample_id, '.self_clustering.anno.csv'),num_threads = 10)

# t1 <- as.character(srobj$self_cluser_anno)
# t1 <- 
# 2.1 visualization-spatial loc. ------------------------------------------
Idents(srobj) <- 'self_cluser_anno'
#CXCL9-EPI
pdf(file = paste0(plotDir,'/',projid,'.Epi-VEGFA_CXCL9-SPP1.sdim.pdf'),
    width = 3*4,height = ceiling(nlevels(srobj)/4)*3)
p1 <- function(subobj=sr_sub){
  SpatialDimPlot(subobj,
                 pt.size.factor = 10,
                 label.size = 4,cols = self_cluster_color)  
}
sr_sub <- subset(srobj,idents = c("Epi1-ADH1C", "Macro-CXCL9", "Macro-SPP1"))
print(p1())
sr_sub <- subset(srobj,idents = c("Epi2-ELF3", "Macro-CXCL9", "Macro-SPP1"))
print(p1())
sr_sub <- subset(srobj,idents = c("Epi3-VEGFA", "Macro-CXCL9", "Macro-SPP1"))
print(p1())
dev.off()

pdf(file = paste0(plotDir,'/',projid,'.Epi_in_one-VEGFA_CXCL9-SPP1.sdim.pdf'),
    width = 3*4,height = ceiling(nlevels(srobj)/4)*3)
sr_sub <- subset(srobj,idents = c("Epi1-ADH1C","Epi2-ELF3", "Epi3-VEGFA", "Macro-CXCL9", "Macro-SPP1"))
tempcolrs <- c('#43a2ca','#43a2ca','#43a2ca','#d0361b','#ca347b')
names(tempcolrs) <- c("Epi1-ADH1C","Epi2-ELF3", "Epi3-VEGFA", "Macro-CXCL9", "Macro-SPP1")
SpatialDimPlot(sr_sub,
               pt.size.factor = 10,
               label.size = 4,cols = tempcolrs)  
dev.off()
#SPP1-EPI
#SPP1-CXCL9
#CXCL9-SPP1-EPI

pdf(file = paste0(plotDir,'/',projid,'.annoed_niche.splity-cluster.sdim.pdf'),
    width = 3*4,height = ceiling(nlevels(srobj)/4)*3)
for (i in levels(srobj)) {
  sr_sub <- subset(srobj,idents = c(i))
  print(SpatialDimPlot(sr_sub,
                       #group.by = '',
                       pt.size.factor = 10,
                       cols = self_cluster_color))
}
dev.off()


# 2.2 top genes -----------------------------------------------------------
mk.sr <- qs::qread(file = './/data/ST/seurat/raw/revised_clustering/020021-1b/LMST.020021-1b.Spatial.008um.mk.sr.qs')

sig_gene <- 
  list(
    Epi=c('EPCAM','KRT17','KRT19'),
    myeloid=c('CD68','LYZ'),
    SPP1=c('SPP1','CCL3','CCL3L1','CXCL8'),
    CXCL9=c('CXCL9','CXCL10','CXCL11'),
    T_cell=c('CD3D','CD8A','CXCL13','CD4'),
    effect=c('GZMA','GZMK','IL2','IFNG','TNF'),
    exhaustion=c('PDCD1','CD274','LAG3','TIGIT','HAVCR2','ENTPD1'),
    regulatory=c('FOXP3','IL2RA','TNFRSF4')
  )
mkdf <- namedList2LongDf(sig_gene)
pdf(file = paste0(plotDir,'/', projid,'.sig_gene.dot.pdf'),
    height = round(nlevels(srobj)/5) + 2,
    width = round(length(mkdf$gene)/5) +3)
# plotDot2(expr_mat = srobj,
#          markerDf = namedList2LongDf(sig_gene),
#          min_expr = 0,
#          minPercent = 0,
#          zscore.min = -2,
#          zscore.max = 2)
DotPlot(object = srobj,features = unlist(sig_gene))
dev.off()
