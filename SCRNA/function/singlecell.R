####################################################
#Sparse mat to common mat
#a function called as_matrix to convert sparse matrix to conmmon matrix.

#############################################
# Rcpp::sourceCpp(code='
# #include <Rcpp.h>
# using namespace Rcpp;
# 
# 
# // [[Rcpp::export]]
# IntegerMatrix asMatrix(NumericVector rp,
#                        NumericVector cp,
#                        NumericVector z,
#                        int nrows,
#                        int ncols){
# 
#   int k = z.size() ;
# 
#   IntegerMatrix  mat(nrows, ncols);
# 
#   for (int i = 0; i < k; i++){
#       mat(rp[i],cp[i]) = z[i];
#   }
# 
#   return mat;
# }
# ' )
# 
# 
# as_matrix <- function(mat){
#   
#   row_pos <- mat@i
#   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
#   
#   tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
#                   nrows =  mat@Dim[1], ncols = mat@Dim[2])
#   
#   row.names(tmp) <- mat@Dimnames[[1]]
#   colnames(tmp) <- mat@Dimnames[[2]]
#   return(tmp)
# }


#dealwith anaPart.
get_input_sr_path <- function(
    origpath='./data/SCRNA/seurat/rawAnalysis/subDaGBC-TLIVER.T_cells_NK_cells.srobj.addMainCellType.qs'){
  if(! grepl(pattern = 'Fil',x = anaPart,ignore.case = T)){
    #mypath <- './data/SCRNA/seurat/mainAnno/soupxRoundedSinglelet/subDa/BMLM.T_cells_NK_cells.srobj.addMainCellType.qs'
    #mypath <- './data/SCRNA/seurat/rawAnalysis/subDaGBC-TLIVER.Myeloid_cells.srobj.addMainCellType.qs'
    mypath <- origpath
  }else{
    strTial <- gsub(pattern = 'Fil', 
                    replacement = '',
                    x = anaPart,ignore.case = T)
    if(strTial == ''){
      newstr <- ''
    }else{
      strTial <- as.numeric(strTial)
      strTial <- strTial - 1
      if(strTial == 1){
        newstr <- 'Fil'
      }else{
        newstr <- paste0('Fil', strTial)
      }
    }
    mypath <- paste0('./data/SCRNA/seurat/fineAnno/',
                     ppart,'/', runPart, '/',newstr, '/', projid, '.', ppart,
                     '.afterSubset.sr.qs')
  }
  return(mypath)
}

getPfx <- function(){
  pfx <- ''
  if(exists('projid')){
    pfx <- projid
  }
  if(exists('ppart')){
    pfx <- paste0(pfx, '.', ppart)
  }
  return(pfx)
}

sr_basic <- function(sr=srobj,
                     norm_method = 'LogNorm',
                     batchCol,
                     usingAssay,
                     regressCol=c("nCount_RNA", "percent.mt", 
                                  "S.Score", "G2M.Score"),
                     sr_from_raw_counts=T){
  if(sr_from_raw_counts){
    sr[["percent.mt"]] <- PercentageFeatureSet(sr, pattern = "^MT-")
    sr[['percent.rb']] <- PercentageFeatureSet(sr, pattern = "^RP[SL]")
    #计算红细胞比例
    HB.genes <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'BHG1','BHG2','HBM', 'HBQ1', 'HBZ')
    #if gene names not in rownames of srobj will casue error
    HB.genes <- HB.genes[HB.genes %in% rownames(sr)]
    sr[['percent.HB']] <- PercentageFeatureSet(sr, features = HB.genes)
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    sr <- CellCycleScoring(sr, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  }
    
  
  if(norm_method == 'LogNorm'){
    #LogNorm
    sr <- NormalizeData(sr)
    sr <- FindVariableFeatures(sr, nfeatures = 2000)
    sr <- ScaleData(sr, vars.to.regress = regressCol)
  }else{
    if('nCount_RNA' %in% regressCol){
      regressCol <- regressCol[regressCol != 'nCount_RNA']}
    print("using SCT to scale.")
    print(paste0('the following vars are regressed:',
                 paste0(regressCol, collapse = ',')))
    sr <- SCTransform(sr,
                      method = "glmGamPoi",
                  vars.to.regress = regressCol,
                  verbose = T,
                  return.only.var.genes = T)
  }
  sr <- RunPCA(object = sr,features = VariableFeatures(sr), verbose = F)
  return(sr)
}


#v20240820
#sr to run scale and harmony.
library(harmony)
sr2harmony <- function(sr=srobj,
                       norm_method = 'LogNorm',
                       batchCol,
                       usingAssay,
                       regressCol=c("nCount_RNA", "percent.mt", 
                                    "S.Score", "G2M.Score"),
                       sr_from_raw_counts=T){
  sr <- sr_basic(sr = sr,
                 norm_method = norm_method,batchCol = batchCol,
                 usingAssay = usingAssay,sr_from_raw_counts=sr_from_raw_counts)
  print('completing basic seurat pipeline')
  # sr[["percent.mt"]] <- PercentageFeatureSet(sr, pattern = "^MT-")
  # sr[['percent.rb']] <- PercentageFeatureSet(sr, pattern = "^RP[SL]")
  # #计算红细胞比例
  # HB.genes <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'BHG1','BHG2','HBM', 'HBQ1', 'HBZ')
  # #if gene names not in rownames of srobj will casue error
  # HB.genes <- HB.genes[HB.genes %in% rownames(sr)]
  # sr[['percent.HB']] <- PercentageFeatureSet(sr, features = HB.genes)
  # s.genes <- cc.genes.updated.2019$s.genes
  # g2m.genes <- cc.genes.updated.2019$g2m.genes
  # sr <- CellCycleScoring(sr, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  # 
  # 
  # if(norm_method == 'LogNorm'){
  #   #LogNorm
  #   sr <- NormalizeData(sr)
  #   sr <- FindVariableFeatures(sr, nfeatures = 2000)
  #   #all.genes <- rownames(inobj)
  #   sr <- ScaleData(sr, vars.to.regress = regressCol)
  #   sr <- RunPCA(object = sr,
  #                features = VariableFeatures(sr),
  #                verbose = F)
  # }else{
  #   if('nCount_RNA' %in% regressCol){
  #     regressCol <- regressCol[regressCol != 'nCount_RNA']}
  #   print(paste0('the following vars are regressed:',
  #                paste0(regressCol, collapse = ',')))
  #   sr <- sr %>% 
  #     SCTransform(method = "glmGamPoi",
  #                 vars.to.regress = regressCol,
  #                 verbose = FALSE,
  #                 return.only.var.genes = F);
  #   sr <- RunPCA(object = sr,features = VariableFeatures(sr), verbose = F)
  # }
  
  if(!is.null(batchCol)){
    if(!exists(usingAssay)){
      if(norm_method == 'LogNorm'){
        usingAssay = 'RNA'
      }else{
        usingAssay = 'SCT'
      }
    }
    print('########################')
    print(paste0('using ', usingAssay, ' to run harmony!!!'))
    print('########################')
    #sct 2 harmony
    sr <- sr %>% 
      RunHarmony(batchCol,
                 plot_convergence = TRUE,
                 assay.use = usingAssay,
                 project.dim = F,
                 epsilon.cluster=-Inf,
                 epsilon.harmony=-Inf)
    sr@reductions$harmony@stdev <- apply(sr@reductions$harmony@cell.embeddings,2,sd)
  }
  return(sr)
}
pred_singleR <- function(sr,ref){
  #singleR prediction.
  #for human cell
  #auto celltype anno 
  sr.sceobj <- as.SingleCellExperiment(sr)
  #temdf <- as.data.frame(Idents(srobj))
  #colnames(temdf) <- 'ident'
  #temdf <- cbind(srobj@meta.data, temdf)
  #sr.sceobj@colData <- DataFrame(temdf)
  #main celltype prediction
  pred <- SingleR(test=sr.sceobj,ref=ref, labels=ref$label.main,
                  clusters = Idents(sr))
  table(pred$labels)
  tab <- table(Assigned=pred$pruned.labels,
               Cluster=sort(unique(Idents(sr))))
  if(nrow(tab) == 1){
    tab <- rbind(tab, rep(0,ncol(tab)))
  }
  p1 <- pheatmap::pheatmap(tab,
                           color=colorRampPalette(c("white", "blue"))(101),
                           cluster_cols = F,
                           cluster_rows = F)
  cleanDev()
  #fine celltype prediction
  predFiner <- SingleR(test=sr.sceobj, ref=ref, labels=ref$label.fine,
                       clusters = Idents(sr))
  tabFiner <- table(Assigned=predFiner$pruned.labels, 
                    Cluster=sort(unique(Idents(sr))))
  if(nrow(tabFiner) == 1){
    tabFiner <- rbind(tabFiner, rep(0,ncol(tabFiner)))
  }
  cleanDev()
  p2 <- pheatmap::pheatmap(tabFiner,
                           color=colorRampPalette(c("white", "blue"))(101),
                           cluster_cols = F,
                           cluster_rows = F)
  pred_labels <- data.frame(idents=levels(sr),
                            pred_main=pred$labels,
                            pred_fine=predFiner$labels)
  resl <- list()
  resl$plot <- list(p1,p2)
  names(resl$plot) <- c('pred_anno_main','pred_anno_fine')
  resl$data <- list(pred, predFiner, pred_labels)
  names(resl$data) <- c('pred_anno_main','pred_anno_fine', 'pred_labels')
  return(resl)
}




namedList2LongDf <- function(mk_list=bbrmk,
                             mklevels=NULL,
                             nameCol='CELLTYPE',
                             valueCol='gene'){
  #suppose mk_list is a named list, whose names are celltypes and element in list are genes.
  #mklevels can used to set cell type levels.
  mydf <- vector(mode = 'list',length = length(mk_list))
  for (i in 1:length(mk_list)) {
    df <- data.frame(CELLTYPE = names(mk_list)[i],
                     gene = mk_list[[i]])
    mydf[[i]] <- df
  }
  mydf <- do.call('rbind', mydf)
  if(is.null(mklevels)){
    mydf$CELLTYPE <- factor(x = mydf$CELLTYPE, levels = names(mk_list))
  }else{
    mydf$CELLTYPE <- factor(x = mydf$CELLTYPE, levels = mklevels)
  }
  
  colnames(mydf) <- c(nameCol, valueCol)
  return(mydf)
}

#the df's col named is a cell type, and row are genes
df2clearList <- function(yourdf){
  
  list2ClearList <- function(yourlist){
    #yourlist <- t1
    #remoe empty ele.
    yourlist <- yourlist[!sapply(yourlist, is.null)]
    #remove na or emtpy string.
    matBool <- sapply(yourlist,function(x){!(is.na(x)| (x == ''))})
    for (i in 1:length(yourlist)) {
      yourlist[[i]] <- yourlist[[i]][unlist(matBool[[i]])]
    }
    return(yourlist)
  }
  
  tems <- vector(mode = 'list', length = ncol(yourdf))
  col_s <- colnames(yourdf)
  for (i in 1:ncol(yourdf)) {
    temsi <- unique(unlist(yourdf[[col_s[i]]]))
    tems[[i]] <- temsi
  }
  names(tems) <- col_s
  tems <- list2ClearList(tems)
  return(tems)
}
list2ClearList <- function(yourlist){
  #yourlist <- geneSets[1:2]
  #remoe empty ele.
  yourlist <- yourlist[!sapply(yourlist, is.null)]
  #remove na or emtpy string.
  yourlist <- lapply(yourlist,function(x){x[!(is.na(x)| (x == ''))]})
  return(yourlist)
}


#receive a chr or factor vec
#return a sorted unique vec
ChrVec2uniqVec <- function(yourChrVec){
  if(all(is.factor(yourChrVec))){
    resVec <- levels(yourChrVec)
    resVec <- resVec[resVec %in% yourChrVec]
    return(resVec)
  }else{
    return(sort(unique(yourChrVec)))
  }
}

df2namedlist <- function(yourdf,
                         geneCol='gene',
                         nameCol='mp'){
  
  #df format.
  #  gene mp
  #g1 p1
  #g2 p1
  #g3 p1
  #g4 p2
  #g5 p2
  #g6 p2
  
  yourdf <- yourdf[,c(geneCol, nameCol)]
  colnames(yourdf) <- c('gene', 'name')
  nameVec <- ChrVec2uniqVec(yourChrVec = yourdf$name)
  
  resList <- vector(mode = 'list',length = length(nameVec))
  for (i in 1:length(nameVec)) {
    resList[[i]] <- unique(unlist(yourdf[yourdf$name == nameVec[i],'gene',drop=T]))
  }
  names(resList) <- nameVec
  return(resList)
}

#get cell marker top gene df from seurat fineAllMarker format res.
topMk <- function(mkDF,
                  topn=5,
                  groupbyCol='cluster',
                  orderbyCol='avg_log2FC',
                  returnAllCol=F){
  mkDF$gbyCol <- unlist(mkDF[,groupbyCol])
  mkDF$obyCol <- unlist(mkDF[,orderbyCol])
  #print(head(mkDF))
  if(!('gene' %in% colnames(mkDF))){
    print('gene col not in your mkDF! get gene name from row names!')
    mkDF$gene <- rownames(mkDF)
  }
  
  mkDF <- mkDF %>%
    dplyr::group_by(gbyCol) %>% 
    dplyr::slice_max(n= topn,
              order_by = obyCol) %>% 
    dplyr::mutate(CELLTYPE=gbyCol) %>% 
    dplyr::select(!c(obyCol, gbyCol))
  if(!returnAllCol){
    mkDF <- mkDF[,c('CELLTYPE','gene')]
  }
  return(mkDF)
}

runCOSG <- function(srobj=TNKCell, usingAssay='RNA',ngene=100){
  
  # Run COSG:
  marker_cosg <- COSG::cosg(
    srobj,
    groups='all',
    assay=usingAssay,
    slot='data',
    mu=1,
    n_genes_user=ngene)
  # Check the marker genes:
  #head(marker_cosg$names)
  #head(marker_cosg$scores)
  
  marker_cosg_df <- 
    left_join(x = 
                marker_cosg$names %>% mutate(geneRank=paste0('top', 1:nrow(marker_cosg$names))) %>% 
                pivot_longer(!geneRank, names_to = "cluster", values_to = "gene") %>% 
                dplyr::mutate(id=paste(cluster,geneRank,sep = '_')) %>% 
                dplyr::select(id, cluster, gene),
              y=marker_cosg$scores%>% mutate(geneRank=paste0('top', 1:nrow(marker_cosg$names))) %>% 
                pivot_longer(!geneRank, names_to = "cluster", values_to = "value") %>% 
                dplyr::mutate(id=paste(cluster,geneRank,sep = '_')) %>% 
                dplyr::select(id, value)
              ,
              by='id')
  res_list <- vector(mode = 'list', length = 2)
  res_list[[1]] <- marker_cosg
  res_list[[2]] <- marker_cosg_df
  names(res_list) <- c('cosg', 'cosgDf')
  return(res_list)
}

markerBasedClusterCor <- function(
    markerDF=srAllNCmk,
    seuratObj=srobj,
    ntop=5,
    usingClu='all',
    usingAssay='SCT',
    corMethod=c("pearson"),
    #FOR save res.
    writePlot=T,
    writeData=T,
    returnData=F,
    projectid=projid,
    suffix='TNKCell.allClu',
    resDir='./seurat/plot/detailAna/TNKCell/markerBaseedPCC'
    
){
  #usingClu -> 'all' or a vec in markerDF$cluster
  #usingAssay -> one of c('SCT', 'RNA')
  #corMethod -> one of c("pearson", "kendall", "spearman")
  if('all' %in% usingClu ){
    subDf <- markerDF %>% dplyr::group_by(cluster) %>%
      dplyr::slice_min(order_by = p_val_adj, n = ntop)
  }else{
    subDf <- markerDF %>%
      dplyr::filter(cluster %in% usingClu) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_min(order_by = p_val_adj, n = ntop)
  }
  
  
  if(usingAssay == 'SCT'){
    exprMat <- AverageExpression(object = seuratObj,
                                 features = unique(subDf$gene),
                                 slot = 'data', assays = 'SCT')$SCT
    
  }else{
    if(usingAssay == 'RNA'){
      exprMat <- AverageExpression(object = seuratObj,
                                   features = unique(subDf$gene),
                                   slot = 'data', assays = 'RNA')$RNA
      
    }else{
      print('error!!!!assay not found!')
      break
    }
  }
  if( !('all' %in% usingClu)){
    exprMat <- exprMat[,usingClu]
  }
  
  
  corMat <- cor(exprMat, method = corMethod)
  
  resDirPlot <- paste0(resDir, '/plot')
  resDirData <- paste0(resDir, '/data')
  dir_mk(mypath = resDirPlot)
  dir_mk(mypath = resDirData)
  
  plot1 <- pheatmap::pheatmap(corMat)
  
  if(writePlot){
    pdf(file = paste0(resDirPlot, '/', projectid, '.', suffix, '.', usingAssay, '.ntop', ntop , '.correlationBy',corMethod, '.pheatmap.pdf'),
        width = 5,height = 5)
    print(plot1)
    dev.off()
  }else{
    print(plot1)
  }
  
  if(returnData){
    return(list(exprMat=exprMat,
                corMat=corMat))
  }
  if(writeData){
    write.table(x = as.data.frame(exprMat),
                quote = F, sep = '\t', row.names = T, col.names = T,
                file = paste0(resDirData, '/', projectid, '.', suffix, '.', usingAssay, '.ntop', ntop , '.expr.','tsv') )
    write.table(x = as.data.frame(corMat),
                quote = F, sep = '\t', row.names = T, col.names = T,
                file = paste0(resDirData, '/', projectid, '.', suffix,'.', usingAssay, '.ntop', ntop , '.correlationBy',corMethod, '.','tsv'))
  }
  
  #corRes <- cor()
  
}

clusterCor <- function(
    usingFeatures=VariableFeatures(srobj),
    seuratObj=srobj,
    usingClu='all',
    usingAssay='SCT',
    corMethod=c("pearson"),
    #FOR save res.
    writePlot=T,
    writeData=T,
    returnData=F,
    projectid=projid,
    suffix='TNKCell.allClu',
    resDir='./seurat/plot/detailAna/TNKCell/clusterPCC'
    
){
  #usingClu -> 'all' or a vec in markerDF$cluster
  #usingAssay -> one of c('SCT', 'RNA')
  #corMethod -> one of c("pearson", "kendall", "spearman")
  usingFeatures <- usingFeatures[usingFeatures %in% rownames(seuratObj)]
  print(paste0('using featurs in srobj only and len of usingFeatures:',
               length(usingFeatures)))
  allClu <- sort(unique(Idents(seuratObj)))
  if(!('all' %in% usingClu)){
    allClu <- allClu[allClu %in% usingClu]
  }
  
  if(usingAssay == 'SCT'){
    exprMat <- AverageExpression(object = seuratObj,
                                 features = usingFeatures,
                                 slot = 'data', assays = 'SCT')$SCT
  }else{
    if(usingAssay == 'RNA'){
      exprMat <- AverageExpression(object = seuratObj,
                                   features = usingFeatures,
                                   slot = 'data', assays = 'RNA')$RNA
    }else{
      print('error!!!!assay not found!')
      break
    }
  }
  
  exprMat <- exprMat[,allClu]
  corMat <- cor(exprMat, method = corMethod)
  resDirPlot <- paste0(resDir, '/plot')
  resDirData <- paste0(resDir, '/data')
  dir_mk(mypath = resDirPlot)
  dir_mk(mypath = resDirData)
  plot1 <- pheatmap::pheatmap(corMat)
  if(writePlot){
    pdf(file = paste0(resDirPlot, '/', projectid, '.', suffix, '.', usingAssay, '.correlationBy',corMethod, '.pheatmap.pdf'),
        width = 5,height = 5)
    print(plot1)
    dev.off()
  }else{
    print(plot1)
  }
  
  if(returnData){
    return(list(exprMat=exprMat,
                corMat=corMat))
  }
  if(writeData){
    write.table(x = as.data.frame(exprMat),
                quote = F, sep = '\t', row.names = T, col.names = T,
                file = paste0(resDirData, '/', projectid, '.', suffix, '.', usingAssay, '.expr.','tsv') )
    write.table(x = as.data.frame(corMat),
                quote = F, sep = '\t', row.names = T, col.names = T,
                file = paste0(resDirData, '/', projectid, '.', suffix,'.', usingAssay, '.correlationBy',corMethod, '.','tsv'))
  }
}









#stacked violin plot
#self modify for more compact.
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    #theme_void() +
    scale_y_continuous(n.breaks =  3) +
    ylab(feature) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      #axis.line.x = element_blank(),
      legend.position = 'none',
      plot.margin = plot.margin,
      plot.title = element_blank(),
      axis.text.y = element_text(size=rel(1.5)),
      axis.title.y = element_text(hjust = 0.5,
                                  vjust = 0.5,
                                  angle = 0,
                                  size = rel(1.2),
                                  face = 'bold')
    )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90,face = 'bold',size = rel(1.2)))
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


#plot umap for features.
plot_measure_dim <- function (dataset, measures, split_by = NA, point.size = 1) 
{
  
  get_measure_data <- function (dataset, measures, return_df = TRUE) 
  {
    meta_vars <- colnames(dataset@meta.data)
    var_genes <- VariableFeatures(dataset)
    all_genes <- rownames(dataset)
    df <- dataset@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
      cbind(dataset@meta.data) %>% rownames_to_column(var = "barcode")
    idx <- c()
    for (i in seq_along(1:length(measures))) {
      if (measures[i] %in% meta_vars) {
        next
      }
      else if (measures[i] %in% var_genes) {
        df[[measures[i]]] <- as.numeric(unlist(FetchData(object = srobj, vars = measures[i], slot = 'scale.data')))
      }
      else if (measures[i] %in% all_genes) {
        df[[measures[i]]] <- as.numeric(unlist(FetchData(object = srobj, vars = measures[i], slot = 'scale.data')))
      }
      else {
        warning(paste("Measure name not found:", measures[i]), 
                immediate. = TRUE)
        idx <- append(idx, i)
      }
    }
    if (!is.null(idx)) 
      measures <- measures[-idx]
    if (return_df) {
      return(df)
    }
    else {
      return(list(df, measures))
    }
  }
  
  
  l <- get_measure_data(dataset = dataset, measures = measures, 
                        return_df = FALSE)
  
  split_by <- ifelse(split_by == "No Split", NA, split_by)
  
  df <- l[[1]]
  measures <- l[[2]]
  p <- list()
  for (i in seq_along(1:length(measures))) {
    qmin <- quantile(df[[measures[i]]], probs = 0.1)
    qmax <- quantile(df[[measures[i]]], probs = 0.9)
    p[[i]] <- ggplot(df) + geom_point(aes(x = .data$UMAP_1, 
                                          y = .data$UMAP_2, color = .data[[measures[i]]]), 
                                      size = point.size) +
      #scale_color_viridis_c(option = "A",
      #                      name = "",
      #                      direction = -1,
      #                      limits = c(quantile(.data[[measures[i]]], probs = 0.1),
      #                                 quantile(.data[[measures[i]]], probs = 0.9)),
      #                      oob = scales::squish) + 
      scale_color_gradient(low = 'navyblue', high = 'red',
                           limits = c(qmin, qmax),
                           oob = scales::squish) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
            panel.border = element_rect(colour = "black", fill = NA, 
                                        size = 1, linetype = 1), axis.line = element_blank()) + 
      labs(x = "UMAP_1", y = "UMAP_2") + ggtitle(measures[i]) + 
      if (!is.na(split_by)) {
        facet_wrap(as.formula(paste("~", split_by)))
      }
    else {
      theme(aspect.ratio = 1)
    }
  }
  patchwork::wrap_plots(p)
}


#old color: c('navyblue', 'orangered')
FeaturePlotNOaxis<- function(
  yourobj,
  genevec,
  ncol=4,
  legendPos = 'none',
  colVec=c('lightblue','white', 'firebrick'),
  sortPoints=T,
  ...){
  genevec <- unique(unlist(genevec))
  plist <- vector(mode = 'list', length = length(genevec))
  for (i in 1:length(plist)) {
    plist[[i]] <- FeaturePlot(object = yourobj,
                              features = genevec[i],
                              order = sortPoints, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                              cols = colVec, ...) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            legend.position = legendPos,
            legend.text = element_text(size = rel(1.2)),
            plot.title = element_text(size=rel(1.2)),
            legend.key.size = unit(0.8, "lines"))
  }
  names(plist) <- genevec
  #plist <- aplot::plot_list(gglist = plist,ncol = ncol)
  plist <- patchwork::wrap_plots(plist, ncol = ncol)
  return(plist)
}

#for self raster.
#because featureplot not sort cells when raster=T
library(ggrastr)
#old color: c('navyblue', 'orangered')
FeaturePlotNOaxis2 <- function(
    yourobj,
    genevec,
    ncol=4,
    legendPos = 'none',
    colVec=c('lightblue','white', 'firebrick'),
    sortPoints=T,
    ...){
  genevec <- unique(unlist(genevec))
  plist <- vector(mode = 'list', length = length(genevec))
  for (i in 1:length(plist)) {
    #cells.order.2 <- names(sort(yourobj[['SCT']]@data[genevec[i],], decreasing = T))
    pli <- FeaturePlot(object = yourobj,
                       features = genevec[i],
                       order = sortPoints,
                       #cells= cells.order.2,
                       raster =F, min.cutoff = 'q10', max.cutoff = 'q90',
                       cols = colVec, ...) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            legend.position = legendPos,
            legend.text = element_text(size = rel(1.2)),
            plot.title = element_text(size=rel(1.2)),
            legend.key.size = unit(0.8, "lines"))
    pli <- ggrastr::rasterize(pli, layers='Point', dpi=300)
    plist[[i]] <- pli
  }
  names(plist) <- genevec
  #plist <- aplot::plot_list(gglist = plist,ncol = ncol)
  plist <- patchwork::wrap_plots(plist, ncol = ncol)
  return(plist)
}



#define a function to write and plot DEG.
saveDEG <- function(DEres=lmg02P.Macro_APOE.mk,
                    label='',
                    resDataDir=csfcomparMyeDataDir,
                    resPlotDir=csfComparMyePlotDir,
                    projectid=projid,
                    suffix='LMG02'){
  dir_mk(resDataDir)
  dir_mk(resPlotDir)
  
  write.table(x = DEres,
              file = paste0(resDataDir, '/',projectid, '.', suffix, '.', label, '.DEG.txt'),
              quote = F, sep = '\t', row.names = T, col.names = T)
  saveRDS(object = DEres, file = paste0(resDataDir, '/',projectid, '.', suffix, '.', label, '.DEG.rds'))
  
  myplot <- EnhancedVolcano(DEres,
                            lab = rownames(DEres),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            xlab =  bquote(~Log[2]~ 'fold change'),
                            pCutoff = 10e-30,
                            #col=c('#a6bddb', '#a6bddb', '#a6bddb', '#d95f0e'),
                            colAlpha = 1,legendPosition = 'none',
                            drawConnectors = TRUE,
                            widthConnectors = 0.75,
                            title = label)
  pdf(file = paste0(resPlotDir, '/',projectid, '.', suffix, '.', label, '.withConnector.volcano.pdf'),width = 12, height = 8)
  print(myplot)
  dev.off()
  
  
  myplot <- EnhancedVolcano(DEres,
                            lab = rownames(DEres),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            pCutoff = 10e-30,
                            #col=c('#a6bddb', '#a6bddb', '#a6bddb', '#d95f0e'),
                            colAlpha = 1,legendPosition = 'none',
                            title = label)
  pdf(file = paste0(resPlotDir, '/',projectid, '.', suffix, '.', label, '.withoutConnector.volcano.pdf'),width = 12, height = 8)
  print(myplot)
  dev.off()
}

#
#tem <- c("T Cell Activation", "Cytotoxicity", "Interferon Type I", "Interferon Type II", "Pro_inflammatory", "T Cell Anergy", "T Cell Exhaustion", "Anti_inflammatory")
#define a function to run AUCcell and get mean(using quantile 0.1- 0.9) 
plotScore <- function(s=cahBM.sig[,tem,drop= F],
                      selected_sig=NULL,
                      yourobj = srobj,
                      cellranking=cells_rankings,
                      #celltypes = NULL,
                      returnPlot =F,
                      returnData =F,
                      writeData = T,
                      writePlot=T,
                      resDataDir = paste0(dataDir, '/signatureScore'),
                      resPlotDir = paste0(plotDir, '/signatureScore'),
                      filePfx='TNKCell',
                      fileSfx='hallmark',
                      pH=NULL,
                      pW=NULL,
                      plotNumCol = 4,
                      pointsize = 0.5,
                      legendPos = 'none',
                      seuratCols = namedColors,
                      simplifyPlot=TRUE,
                      sortPoints=T,
                      roundDA=F){
  #s your signature gene, can be list, whose element names are signature names or a df, whose colnames are sig names.
  #yourobj a prepared seurat obj.
  #cellranking -> a AUCcell obj.
  #celltypes -> a two col df,col1 is cellid, col2 is celltype
  #
  #when counting mean expr of sig
  #if a sig's gene num < 5, then using mean expr dircetly, otherise, using q10-q90 gene.
  
  #simplifyPlot=TRUE remove feature plot's axis and other plot ele.
  #sortPoints=T sort points in featureplot
  
  #pW and pH is used only in saving vln plot file.
  #if not given, will auto detected by len of pathways and levels of srobj.
  #each element in featureplot file is 3cm*3cm and the fig file size is saved by 
  #number of element and ncol.
  
  #using patchwork to wrap plots.
  
  #if roundDA=T using
  #meanDa <- round(x = meanDa,digits = 2)
  #aucdata <- round(x = aucdata, digits = 2)
  # to round result. NOTE!!! ROUND DATA WILL Cause vlnplot ugly.
  #testing
  
  # s <- geneSets
  # selected_sig=t1[1:2]
  # yourobj = srobj
  # cellranking=cells_rankings
  # #celltypes = NULL,
  # returnPlot =F
  # returnData =F
  # writeData = T
  # writePlot=T
  # resDataDir = paste0(dataDir, '/signatureScore')
  # resPlotDir = paste0(plotDir, '/signatureScore')
  # filePfx='TNKCell'
  # fileSfx='hallmark'
  # pH=NULL
  # pW=NULL
  # plotNumCol = 4
  # pointsize = 0.5
  # legendPos = 'none'
  # seuratCols = namedColors
  # simplifyPlot=TRUE
  # sortPoints=T
  # roundDA=F
  
  if(is.data.frame(s)){
    print('s is a df and covert s to a list')
    s <-df2clearList(yourdf = s)
  }else{
    if(is.list(s)){
      print('s is a list!')
      s <- list2ClearList(yourlist = s)
    }else{
      print('wrong s format!')
      break
    }
  }
  if(!(all(is.null(selected_sig)))){
    s <- s[names(s) %in% selected_sig]
  }
  if(all(is.null(s)|is.na(s))){
    print('s is empty after filter by selected_sig')
    break
  }
  print(head(names(s)))
  if(is.null(pH)){
    pH <- round(length(s)/2) + 3
  }
  if(is.null(pW)){
    pW <- round(nlevels(yourobj)/4)+3
  }
  #calcuateing AUC score
  cellsAUC <- AUCell_calcAUC(s, cellranking, aucMaxRank=nrow(cellranking)*0.1)
  ##set gene set of interest here for plotting
  aucdata <- vector(mode = 'list', length = length(s))
  for (i in 1:length(s)) {
    tem <- as.data.frame((getAUC(cellsAUC)[names(s)[i], ]))
    colnames(tem) <- names(s)[i]
    aucdata[[i]] <- tem
  }
  aucdata <- do.call('cbind', aucdata)
  #print(head(aucdata))
  
  #calcuatring mean expr
  #using activate assay's data.
  meanDa <- vector(mode = 'list', length = length(s))
  for (i in 1:length(s)) {
    i_names <- names(s)[i]
    i_genes <- s[[i]]
    adai <- FetchData(object = yourobj, vars = i_genes, slot = 'data')
    #print(head(adai))
    if(length(i_genes) >10){
      #using q1 - q9 to counting mean.
      gene_expr <- adai %>% as.tibble() %>%
        summarise(across(.cols = everything(), .fns = mean)) %>% unlist()
      q19 <- quantile(x = gene_expr, c(0.1,0.9))
      #print(q19)
      gene_expr <- gene_expr[((gene_expr>= q19[1]) & (gene_expr <= q19[2]))]
      #print(head(gene_expr))
      adai <- adai[,names(gene_expr)]
      #print(head(adai))
    }
    adai <- adai %>% t() %>% as.data.frame()
    adai <- apply(X = adai, MARGIN = 2, FUN = mean)
    adai <- as.data.frame(adai)
    meanDa[[i]] <- adai
  }
  meanDa <- do.call('cbind', meanDa)
  colnames(meanDa) <- names(s)
  if(roundDA){
    meanDa <- round(x = meanDa,digits = 2)
    aucdata <- round(x = aucdata, digits = 2)
  }
  #res <- vector(mode = 'list', length = 2)
  #res[[1]] <- aucdata
  #res[[2]] <- meanDa
  res <- list(aucdata, meanDa)
  names(res) <- c('AUC_data', 'mead_data')
  if(returnData){
    return(res)
  }
  orig.meta <- yourobj@meta.data
  meta.auc <- cbind(orig.meta, aucdata)
  meta.mean <- cbind(orig.meta, meanDa)
  pauc <- vector(mode = 'list', length = length(s))
  pmean <- vector(mode = 'list', length = length(s))
  
  yourobj@meta.data <- meta.auc
  pvln.auc <- StackedVlnPlot(obj = yourobj, features = names(s), cols = seuratCols)
  if(simplifyPlot){
    for (i in 1:length(s)) {
      #old color:'navyblue', 'orangered'
      pauc[[i]] <- FeaturePlot(object = yourobj,
                               features = c(names(s)[i]),
                               order = sortPoints, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                               cols = c('lightblue','white', 'firebrick'),pt.size = pointsize) +
        theme(legend.position = legendPos,
              legend.text = element_text(size = rel(1.2)),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.key.size = unit(0.8, "lines"))
    }
  }else{
    for (i in 1:length(s)) {
      pauc[[i]] <- FeaturePlot(object = yourobj,
                               features = c(names(s)[i]),
                               order = sortPoints, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                               cols = c('navyblue', 'orangered'),pt.size = pointsize)
    }
  }
  
  yourobj@meta.data <- meta.mean
  pvln.mean <- StackedVlnPlot(obj = yourobj, features = names(s), cols = seuratCols)
  if(simplifyPlot){
    for (i in 1:length(s)) {
      pmean[[i]] <- FeaturePlot(object = yourobj,
                                features = c(names(s)[i]),
                                order = sortPoints, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                                cols = c('navyblue', 'orangered'),pt.size = pointsize) +
        theme(legend.position = legendPos,
              legend.text = element_text(size = rel(1.2)),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.key.size = unit(0.8, "lines"))
    }
  }else{
    for (i in 1:length(s)) {
      pmean[[i]] <- FeaturePlot(object = yourobj,
                                features = c(names(s)[i]),
                                order = sortPoints, raster =T, min.cutoff = 'q10', max.cutoff = 'q90',
                                cols = c('navyblue', 'orangered'),pt.size = pointsize)
    }
  }
  
  names(pauc) <- names(s)
  names(pmean) <- names(s)
  
  if(returnPlot){
    resp <- list(pauc, pvln.auc, pmean, pvln.mean)
    names(resp) <- c('aucFeatureplot','aucVln', 'meanFeatureplot', 'meanVln')
    
    return(resp)
  }
  if(writePlot){
    dir_mk(mypath = resPlotDir)
    #write
    pdf(file = paste0(resPlotDir, '/', filePfx,'.', fileSfx, '.signatureScore.featureplot.pdf'),
        width = 3*plotNumCol,
        height =  6*(round(length(s)/plotNumCol)))
    #print(aplot::plot_list(gglist = c(pauc, pmean), ncol = plotNumCol))
    print(patchwork::wrap_plots(c(pauc, pmean),ncol = plotNumCol))
    dev.off()
    pdf(file = paste0(resPlotDir, '/', filePfx,'.', fileSfx, '.signatureScore.vlnplot.pdf'),
        width = pW,
        height = pH)
    print(pvln.auc)
    print(pvln.mean)
    dev.off()
  }
  if(writeData){
    #wreite data.
    dir_mk(mypath = resDataDir)
    write.table(x = aucdata, file = paste0(resDataDir, '/',filePfx, '.', fileSfx, '.auc.signatureScore.txt'),
                quote = T, sep = '\t', row.names = T, col.names = T)
    saveRDS(object = aucdata,
            file = paste0(resDataDir, '/',filePfx, '.', fileSfx, '.auc.signatureScore.rds'))
    write.table(x = meanDa, file = paste0(resDataDir, '/',filePfx, '.', fileSfx, '.mean.signatureScore.txt'),
                quote = T, sep = '\t', row.names = T, col.names = T)
    saveRDS(object = meanDa, file = paste0(resDataDir, '/',filePfx, '.', fileSfx, '.mean.signatureScore.rds'))
  }
}



#source('/picb/sysgenomics/shanshenbing/projects/LCBM/script/functions/plot.R')
#source('/picb/sysgenomics/shanshenbing/projects/LCBM/script/functions/singlecell.R')
#v2040819.
plotMainCellType <- function(obj=srobjSub,
                             projID=projid,
                             partID=ppart,
                             resDir=annoPlotDir,
                             mk=NULL,
                             isLabel=T,
                             plotFeature=T,
                             fineDot=T){
  #obj <- TNKCell
  #projID <- 'test'
  #partID='tnkcell'
  #resDir='./test'
  #mklist <- NULL
  #isLabel=T
  #plotFeature=F
  #fineDot=T
  #print(round(length(unique(unlist(mklist)))/5) +3)
  if(is.null(mk)){
    print('NOt given mklist and using marker list in this function!')
    mk <- list(
      Cancer_cells=c('EPCAM', 'KRT18','KRT19'),
      B_cells= c('MS4A1', 'CD19','JCHAIN', 'CD79A', 'MZB1', 'IGHA2','IGHG3','IGHM'),
      T_cells=c('TRAC','CD2','CD3D', 'CD3D','CD3E','CD3G'),
      CD8=c('CD8A', 'CD8B'),
      CD4=c('CD4'),
      NK_cells=c('NKG7','NCAM1', 'FCGR3A', 'KLRD1', 'KLRF1','KLRB1','GNLY'),
      #Macrophages=c('CD68'),
      #Monocytes=c('CD14','FCGR3A','LYZ','MS4A7'),
      Myeloid_cells=c('CD14','FCGR3A','CD68','MARCO','LYZ'),
      Neutrophils=c('CSF3R', 'FPR1', 'FCGR3B', 'NAMPT', 'MNDA'),
      conventional_DC_cells=c('FCER1A', 'CST3', 'FCER1A10'),
      plasmacytoid_DC_cells=c('IL3RA', 'CLEC4C'),
      Mast_cells=c('GATA2','MS4A2','KIT'),
      Platelet=c('PPBP', 'PF4'),
      Mural_vascular_cells=c('RGS5', 'ACTA2'),
      Endothelial_cells=c('RAMP2','FLT1','CLDN5', 'VWF', 'PECAM1'),
      `Astrocytes/Oligodendrocytes`=c('CRYAB','UGT8','CNP','OLIG2','FA2H'),
      #Oligodendrocytes=c('CLDN11','MOG','OLIG2','OLIG1'),
      Mesenchymal_cell=c('ISLR', 'CTHRC1'),
      Red_blood_cells=c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'BHG1','BHG2','HBM', 'HBQ1', 'HBZ'),
      Fibroblasts=c('COL1A2','COL1A1','THY1','DCN'),
      Hepatocyte=c( "TF", "ALB","APOC1","APOA1","FABP1","FGA", "FGB","HP")
    )
  }
  
  if(is.data.frame(mk)){
    mkdf <- mk
    mklist <- df2namedlist(yourdf = mk,geneCol = 'gene',nameCol = 'CELLTYPE')
  }else{
    mklist <- mk
    mkdf <- namedList2LongDf(mk_list = mklist,mklevels = names(mklist))
    mkdf <- mkdf %>% dplyr::arrange(CELLTYPE)
  }
  #prepare dir
  dir_mk <- function(mypath){
    if(!(dir.exists(paths = mypath))){
      dir.create(path = mypath, recursive = T)
    }
  }
  dir_mk(paste0(resDir, '/', partID))
  #feature plot.
  if(plotFeature){
    dir_mk(mypath = paste0(resDir, '/', partID, '/featurePlot'))
    for(i in 1:length(mklist)){
      gene_i <- mklist[[i]]
      gene_i <- gene_i[gene_i %in% rownames(obj)]
      if(length(gene_i) > 0){
        myfig <- FeaturePlot(object = obj,
                             features = gene_i,
                             label = isLabel,
                             ncol = 1,
                             raster = T) +
          plot_annotation(names(mklist)[i],
                          theme = theme(plot.title = element_text(size = 20, face = 'bold', colour = 'black')))
        
        celltypeid <- str_replace_all(string = names(mklist)[i], pattern = '[\\\\/]',replacement = '_')
        print(celltypeid)
        
        pdf(paste0(resDir, '/', partID, '/featurePlot/', projID, '.', partID, '.',celltypeid, '.FeaturePlot.pdf'),
            width = 4,
            height = 4*length(unique(mklist[[i]])))
        print(myfig)
        dev.off()
      }
    }
  }
  
  ##dot plot
  if(fineDot){
        myfig <- plotDot2(expr_mat = obj,
               markerDf = mkdf,
               min_expr = 0,
               minPercent = 0,
               zscore.min = -2,
               zscore.max = 2)
    pW <- round(nlevels(obj)/5) + 4
  }else{
    myfig <- DotPlot(object = obj, features = unique(unlist(mklist)),
                     cols = c('blue','red')) + RotatedAxis()
    pW <- round(nlevels(obj)/5)
  }
  
  pdf(file = paste0(resDir, '/', partID, '/', projID, '.', partID, '.mainType.markerGene.dotplot.pdf'),
      width = round(length(unique(unlist(mklist)))/5) +3,
      height = pW)
  print(myfig)
  dev.off()
  
}





genebasedCor <- function(srobj,
                         geneFile='',# a df with col:chr start end gene_symbol.
                         chromeSome=paste0('chr', 1:22),
                         sampleVec=c('LMG01','LMG01L'),
                         groupVar=NULL,
                         usingHighVariGene=F,
                         usingAssay='SCT',
                         corMethod='pearson',
                         resDir='./',
                         filePfx='test',
                         pW=4,
                         pH=4,
                         returnData=F,
                         writeData=T
){
  
  #srobj=srobj
  #geneFile='',# a df with col:chr start end gene_symbol.
  #chromeSome=paste0('chr', 1:22)
  #chromeSome='chr19'
  #sampleVec=c('LMG01','LMG01L')
  #groupVar=lmg01
  #corMethod='pearson'
  
  #filter geneFile
  colnames(geneFile) <- c('chr', 'start', 'end', 'gene')
  
  srobj <- subset(x = srobj,orig.ident %in% sampleVec)
  #using cell in groupVar.
  srobj <- srobj[,rownames(groupVar)]
  srobj$cellid <- rownames(srobj@meta.data)
  
  groupVar <- groupVar[,'group',drop=F]
  srobj <- AddMetaData(object = srobj,metadata = groupVar)
  srobj.meta <- srobj@meta.data
  Idents(srobj) <- 'group'
  
  geneFile <- geneFile[geneFile$chr %in% chromeSome,]
  
  if(usingHighVariGene){
    srobj <- FindVariableFeatures(object = srobj)
    geneFile <- geneFile[geneFile$gene %in% VariableFeatures(srobj),]
  }
  
  
  #fetch gene.
  if(usingAssay == 'SCT'){
    exprMat <- AverageExpression(object = srobj,
                                 features = unique(geneFile$gene),
                                 slot = 'data', assays = 'SCT')$SCT
  }else{
    exprMat <- AverageExpression(object = srobj,
                                 features = unique(geneFile$gene),
                                 slot = 'data', assays = 'RNA')$RNA
  }
  
  if(returnData){
    return(srobj)
  }
  
  corMat <- cor(exprMat, method = corMethod)
  
  ddir <- paste0(resDir, '/data')
  pdir <- paste0(resDir, '/plot')
  dir_mk(mypath = ddir)
  dir_mk(mypath = pdir)
  
  pdf(file = paste0(resDir, '/', filePfx, '.cor.heatmap.pdf'),
      width = pW,height = pH)
  pheatmap::pheatmap(mat = corMat)
  dev.off()
  
  
  if(writeData){
    saveRDS(object = srobj,
            file = paste0(ddir, '/', filePfx, '.sr.CNVCellOnly.rds'))
    
  }
  
  
}


#outlier point: Q3 + 1.5 * IQR、 Q1 - 1.5 * IQR
#extreme point: Q3 + 3 * IQR、Q1 - 3 * IQR
# Get cell embeddings and feature loadings
rmOutlier <- function(yourobj,
                      yourMeta,
                      clusterCol='fineAnno',
                      idinMetaCol='cellid',
                      reduc='umap',
                      n=3){
  #yourobj=TNKCell
  #yourMeta=sr.fineAnno.meta
  #clusterCol='fineAnno'
  #idinMetaCol='cellid'
  #reduc='umap'
  #n=3
  
  myembed <- Embeddings(object = yourobj, reduction  = reduc)
  cellidVec <- rownames(myembed)
  myembed <- as.data.frame(myembed)
  myembed$cellid <- cellidVec
  
  yourMeta$cellid <- unlist(yourMeta[,idinMetaCol])
  yourMeta$clusters <- unlist(yourMeta[,clusterCol])
  myembed <- left_join(myembed,yourMeta[,c('cellid', 'clusters')],by='cellid')
  
  getQuant <- function(x,num=n){
    q25 <- quantile(x, 0.25)
    q75 <- quantile(x, 0.75)
    iqr <- q75 -q25
    q_l <- q25 - num*iqr
    q_h <- q75 + num*iqr
    res=c(q_l,q_h)
  }
  #low quantile
  q_l <- myembed %>% group_by(clusters) %>% summarise(UMAP_1L=getQuant(x = UMAP_1)[1],
                                                      UMAP_2L=getQuant(x = UMAP_2)[1],)
  q_h <- myembed %>% group_by(clusters) %>% summarise(UMAP_1H=getQuant(x = UMAP_1)[2],
                                                      UMAP_2H=getQuant(x = UMAP_2)[2],)
  
  myembed <- left_join(myembed,x = q_l)
  myembed <- left_join(myembed,x = q_h)
  myembed <- myembed %>% group_by(clusters) %>% 
    filter((UMAP_1 > UMAP_1L) & (UMAP_1 < UMAP_1H) &
             (UMAP_2 > UMAP_2L) & (UMAP_2 < UMAP_2H))
  return(myembed)
}



#make UMAP1/2 arrow like axis.
bsaxis <- function(len=1.5){
  axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(len, "cm"))
  return(guides(x = axis, y = axis))
}
saxis <- function(...){
  tem_obj <- theme(axis.line = element_line(arrow = arrow(type = "closed",angle = '30',length = unit(0.25, 'cm'))),
                   axis.title = element_text(hjust = 0,size=rel(0.8)),
                   panel.border =element_blank(),
                   axis.ticks = element_blank(),axis.text = element_blank(),...)
  return(tem_obj)
}

var2Dimplot <- function(sr=srobj,
                        var_vec=c('idents','tissue_type'),
                        isLabel=T,
                        group_color_list=NULL,
                        pt_size='auto',
                        pw=5,ph=5,
                        plotDir=NULL,
                        writePlot=F){
  numCells <- ncol(sr)
  if(pt_size == 'auto'){
    #1k 3
    #1k-1w 2
    #1w-5w 1.5
    #>5w 1
    pt_size <- ifelse(numCells > 100000,1, 
                      ifelse(numCells > 10000, 1.5,
                             ifelse(numCells > 1000, 2,3)))
  }
  print(paste0('pt_size:', pt_size))
  #group_color_list=NULL a named list, MSUT MATCHING var_vec.
  pfx <- getPfx()
  # 检查plotDir参数是否为NULL
  if (is.null(plotDir)) {
    # 如果为NULL，则检查全局环境中是否存在plotDir变量
    if (exists("plotDir", envir = .GlobalEnv)) {
      # 如果存在，则将全局变量的值赋给plotDir
      plotDir <- get("plotDir", envir = .GlobalEnv)
    } else {
      # 如果不存在，则将plotDir设置为"./“
      plotDir <- "./"
    }
  }
  var_len <- length(var_vec)
  #res_l
  res_l <- list()
  sr@meta.data$idents <- Idents(sr)
  for(i in var_vec){
    vari_vec <- unique(sr@meta.data[,i])
    vari_len <- length(vari_vec)
    #check cols.
    rcols_vari <- group_color_list[[i]]
    if(is.null(rcols_vari)){
      rcols_vari <- circlize::rand_color(n = vari_len,luminosity = 'bright')
    }
    resi <- list()
    p1_group <- DimPlot(object = sr,raster = T,pt.size = pt_size,label = isLabel,
                        group.by = i,cols = rcols_vari)+
      bsaxis()+saxis(aspect.ratio = 1)
    pwi <- pw+ceiling(vari_len/10)
    
    resi$p <- p1_group
    resi$pw <- pwi
    resi$ph <- ph
    resi$file <- paste0(plotDir,'/', pfx,
                              '.groupBy.',i, '.dim.pdf')
    res_l[[i]] <- resi
  }
  #names(res_l) <- var_vec
  if(writePlot){
    dir_mk(plotDir)
    for (i in names(res_l)) {
      ploti <- res_l[[i]]
      pdf(file = ploti$file,
          width = ploti$pw,height = ploti$ph)
      print(ploti$p)
      dev.off()
      # if(exists(ploti[['split']])){
      #   pdf(file = ploti$split$file,
      #       width = ploti$split$pw,height = ploti$split$ph)
      #   print(ploti$split$plot)
      #   dev.off()
      }
  }else{
    return(res_l)
  }
}




ldf2stack <- function(ldf){
  #col1:x;col2:y
  ocol <- colnames(ldf)
  colnames(ldf) <- c('orig.ident', 'idents')
  dfi <- ldf %>%
    group_by(orig.ident, idents) %>% 
    summarise(num=n()) %>% 
    arrange(orig.ident, desc(idents)) %>% 
    group_by(orig.ident) %>% mutate(lab_y=cumsum(num))
  colnames(dfi)[1:2] <- ocol
  return(dfi)
}
sr2num <- function(sr=srobj,
                   by=NULL,
                   stack_col='orig.ident',
                   mappingColor=NULL,
                   plotLoli=T,
                   plotStack=T){
  if(is.data.frame(sr)){
    print('sr is a data.frame')
    sr.meta <- sr
    if(is.null(by)){
      print('by CANNOT BE NULL when sr is data.frame!!!')
      break
    }
  }else{
    sr.meta <- sr@meta.data
  }
  
  resl <- list()
  #for loli
  if(plotLoli){
    if(is.null(by)){
      #by <- 'idents'
      cell_num <- as.data.frame(table(Idents(sr)))
    }else{
      cell_num <- as.data.frame(table(sr[[by]]))
    }
    colnames(cell_num) <- c('idents', 'num')
    p_loli <- loli(df=cell_num,x='idents',y='num')
    resl$plot$loli <- p_loli
    resl$data$loli <- cell_num
  }
  #for stacked bar
  if(plotStack){
    if(is.null(by)){
      #by <- 'idents'
      num_by_sample_idents <- cbind(sr[[stack_col]], as.data.frame(Idents(sr)))
    }else{
      num_by_sample_idents <- sr.meta[,c(stack_col, by)]
    }
    colnames(num_by_sample_idents) <- c('orig.ident', 'idents')
    num_by_sample_idents <- ldf2stack(ldf = num_by_sample_idents)
    plist <- stackBar(df = num_by_sample_idents,mappingColor = mappingColor)
    p_stack <- plist$stack
    p_stack_lab <- plist$stack_lab
    # 
    # if(!is.null(mappingColor)){
    #   p_stack <- p_stack + scale_fill_manual(values = mappingColor)
    # }
    # p_stack_lab <- p_stack+ 
    #   geom_text(data = num_by_sample_idents,
    #             aes(x=orig.ident,y=lab_y, label= num), vjust=0.5, color="black", size=3)
    
    resl$plot$stack <- p_stack
    resl$plot$stack_lab <- p_stack_lab
    resl$data$stack <- num_by_sample_idents
  }
  return(resl)
}


readFineAnnoSr <- function(p_arg,a_arg){
  #project levels args.
  projid <- p_arg$projid
  #batchCol <- p_arg$batchCol
  #LogNorm or SCT.
  norm_method <- p_arg$norm_method
  
  #analysis part args.
  for (i in names(a_arg)) {
    assign(i,a_arg[[i]])
  }
  # ppart <- a_arg$ppart
  # runPart <- a_arg$runPart
  # anaPart <- a_arg$anaPart
  # using_res <- a_arg$using_res
  # dpfx <- a_arg$dpfx
  #projid.ppart
  pfx <- getPfx()
  dataDir <- paste0(mainpath, '/data/',dpfx,'/', ppart, '/', runPart, '/', anaPart)
  #plotDir <- paste0(mainpath, '/data/',dpfx,'/', ppart, '/', runPart, '/', anaPart)
  dir_mk(mypath = c(dataDir))
  input_sr_path <- paste0(dataDir, '/', pfx,'.sr.qs')
  print(paste0('using sr:',input_sr_path))
  sr <- read_Da(input_sr_path)
  ##using resi
  if(norm_method == 'SCT'){
    my_assay <- 'SCT'
  }else{
    my_assay <- 'RNA'
  }
  cluster_by_resn <- paste0(my_assay, '_snn_res.', using_res)
  Idents(sr) <- cluster_by_resn
  print(paste0('using resolution:',cluster_by_resn))
  return(sr)
}
check_cluster<- function(){
  
  #cell num qc
  if(exists('cluster_cols')){
    t1 <- sr2num(mappingColor = cluster_cols)
  }else{
    t1 <- sr2num()
  }
  
  write_tsv(x = t1$data$loli,
            file = paste0(dataDir, '/', pfx, '.cellnum.tsv'))
  write_tsv(x = t1$data$stack,
            file = paste0(dataDir, '/', pfx, '.cellnum_by_sample.tsv'))
  pdf(file = paste0(plotDir,  '/',
                    projid, '.', 'cellnumbers.pdf'),
      width = 4 + round(nlevels(srobj)/6),height = 6)
  print(t1$plot$loli)
  dev.off()
  pdf(file = paste0(plotDir,  '/',
                    projid, '.', 'cellnumbers_by_samples-main_cell_type.pdf'),
      width = 4 + round(n_distinct(srobj$orig.ident)/5),height = 6)
  print(t1$plot$stack)
  dev.off()
  pdf(file = paste0(plotDir,  '/',
                    projid, '.', 'cellnumbers_by_samples-main_cell_type.labeled.pdf'),
      width = 4 + round(n_distinct(srobj$orig.ident)/5),height = 8)
  print(t1$plot$stack_lab)
  dev.off()
  
  #metric qc
  pdf(paste0(plotDir, '/',pfx, '.rawQC.pdf'),
      height = 12,width = 3+ceiling(nlevels(srobj)/6))
  print(VlnPlot(srobj, 
                features = c("nFeature_RNA", "nCount_RNA", 
                             'percent.mt', 'percent.rb', 'percent.HB'),
                ncol = 1,log = T,pt.size=0))
  dev.off()
  
  
  #dimplot 
  using_var <- c('idents', batchCol, 'orig.ident')
  t2 <- var2Dimplot(var_vec = using_var)
  for (i in using_var) {
    ploti <- t2[[i]]
    pdf(file = ploti$file,width = ploti$pw,height =ploti$ph )
    print(ploti$p)
    dev.off()
  }
}

# 
# add_clin <- function(sr=srobj, clininfo=clin, byCol=c(orig.ident='sample_id')){
#   #orig.ident='sample_id'
#   #column MUST IN clin!!!
#   #only cols notin sr will used.
#   if(isS4(sr)){
#     sr.meta <- sr@meta.data
#   }else{
#     sr.meta <- sr
#   }
#   #orig cell id 
#   cellid_order <- rownames(sr.meta)
#   sr.meta$cellid <- rownames(sr.meta)
#   noCol <- colnames(clininfo)
#   noCol <- noCol[(noCol %in% colnames(sr.meta))]
#   noCol <- noCol[!(noCol %in%  byCol)]
#   clininfo <- clininfo[,!(noCol)]
#   #clininfo <- clininfo[,!(colnames(clininfo) %in% colnames(sr.meta))]
#   sr.meta <- left_join(x = sr.meta,
#                        y = clininfo,
#                        by=byCol)
#   rownames(sr.meta) <- sr.meta$cellid
#   if(all(rownames(sr.meta) == cellid_order)){
#     print('all columns in your meta is:')
#     print(colnames(sr.meta))
#     return(sr.meta)
#   }else{
#     print('all(rownames(sr.meta) == colnames(sr) output false!!!')
#     print('check your srobj and clin data!!!')
#     break
#   }
# }

#v20241118
col2diff <- function(df=srobj@meta.data,
                     c1='SCT_snn_res.1',c2='mainCellType',
                     decimal_places = 0){
  df=df[,c(c1,c2)]
  #print(head(df))
  colnames(df) <- c('c1', 'c2')
  df$c1 <- as.character(df$c1)
  df$c2 <- as.character(df$c2)
  m1 <- table(df$c1,df$c2)
  percent_table <- prop.table(m1, margin = 1) * 100  # 转换为百分比矩阵
  #percent_labels <- apply(percent_table, 2, function(x) sprintf("%.2f", x))
  # 动态生成格式字符串，保留指定的小数位数
  number_format <- paste0("%.", decimal_places, "f")
  percent_labels <- apply(percent_table, 2, function(x) sprintf(number_format, x))
  
  # 5. 绘制热图
  #by percents
  p1 <- pheatmap::pheatmap(
    percent_table,               # 数据矩阵（百分比）
    display_numbers = percent_labels,  # 在热图中标注百分比
    #number_format = "%.2f",       # 保留两位小数
    color = colorRampPalette(c("lightblue","white", "orangered"))(100), # 颜色渐变
    cluster_rows = F,         # 不进行行聚类
    cluster_cols = T,         # 不进行列聚类
    fontsize_number = 10,         # 标注数字的字体大小
    angle_col = 45,
    main = paste0("Correspondence between ", c1, ' and ', c2)
  )
  
  #by number
  # p2 <- pheatmap::pheatmap(
  #   m1,               # 数据矩阵（百分比）
  #   display_numbers = m1,  # 在热图中标注百分比
  #   #number_format = "%.2f",       # 保留两位小数
  #   color = colorRampPalette(c("lightblue","white", "orangered"))(100), # 颜色渐变
  #   cluster_rows = FALSE,         # 不进行行聚类
  #   cluster_cols = FALSE,         # 不进行列聚类
  #   fontsize_number = 10,         # 标注数字的字体大小
  #   angle_col = 45,
  #   #main = paste0("Correspondence between ", c1, ' and ', c2)
  # )
  # plist <- list(p1,p2)
  # names(plist) <- c('percent', 'number')
  # return(plist)
  return(p1)
}







# uniq_fineAnno_list <- list(
#   `T_cells/NK_cells`=c("CD4Tn_CCR7", "CD4T_CXCL13", "CD4Tcm_LMNA", "CD4Treg_LEF1", "CD4Treg_TNFRSF4",
#                        "CD8Tef_PRF1", "CD8Tex_CXCL13", "CD8Tprf_CXCL13", "CD8T_IFIT1", "CD8T_MT1X",
#                        "CD8Tem_GZMK", "CD8Tem_ZNF683", "Tmait", "Tgd", "NK_FCER1G",
#                        "NK_FCGR3A", "NK_KIT"),
#   Platelets=c("Platelets"),
#   Myeloid_cells=c("cDC1", "cDC2", "cDC3", "M_ACP5", "M_BAG3", "M_C1QA", "M_CCL20",
#                   "M_CCR2", "M_CD209", "M_CD5L", "M_CXCL9", "M_EREG", "M_KCNQ3",
#                   "M_LRRC23", "M_MT1X", "M_prf", 
#                   "Neu_FCGR3B","Neu_IFIT", "Neu_MME", "Neu_PI3"),
#   Mast_cells=c("Mast_cells"),
#   pDC=c("pDC"),
#   B_cells=c("Bn_IFIT3", "Bn_IGHD", "B_DZ_prf", "B_LZ_RGS13", "Bmem_ITGB1",
#             "Bmem_MAST4", "Plasma_CCR10", "Plasma_FNDC3B", "Plasma_PRDX4", "B_ACTG1",
#             "B_MT1X","B_TNF"),
#   Lym_prf=c("CD4Tn_CCR7", "CD4T_CXCL13", "CD4Tcm_LMNA", "CD4Treg_LEF1", "CD4Treg_TNFRSF4",
#             "CD8Tef_PRF1", "CD8Tex_CXCL13", "CD8Tprf_CXCL13", "CD8T_IFIT1", "CD8T_MT1X",
#             "CD8Tem_GZMK", "CD8Tem_ZNF683", "Tmait", "Tgd", "NK_FCER1G",
#             "NK_FCGR3A", "NK_KIT",
#             
#             "Bn_IFIT3", "Bn_IGHD", "B_DZ_prf", "B_LZ_RGS13", "Bmem_ITGB1",
#             "Bmem_MAST4", "Plasma_CCR10", "Plasma_FNDC3B", "Plasma_PRDX4", "B_ACTG1",
#             "B_MT1X","B_TNF"),
#   Endothelial_cells=c("Endo_ACKR1", "Endo_FCN2", "Endo_GNG11", "Endo_NOTCH4", "Endo_PROX1"),
#   Fibroblasts=c("Fibro_C3", "Fibro_COL11A1", "Fibro_COLEC11",
#                 "Fibro_CXCL10", "Fibro_ENTPD2", "Fibro_prf", "Fibro_SLC24A3","Fibro_WNT11"),
#   Mural_vascular_cells=c("Mural_ACTG2", "Mural_NOTCH3"),
#   Hepatocytes=c("Hepatocytes"),
#   Neuron_like_cells=c("Neuron_like_cells"),
#   Epithelial_cells=c("N_Epi_c01_DDT", "N_Epi_c02_PKHD1", "N_Epi_c03_ITPKC", "N_Epi_c04_REG1A",
#                      
#                      "CC_Epi_c01_ANO3", "CC_Epi_c02_TGM2", "CC_Epi_c03_RND1", 
#                      "CC_Epi_c04_S100A4","CC_Epi_c05_IGFBP2","CC_Epi_c06_HIST1H3J",
#                      "CC_Epi_c07_COL17A1", "CC_Epi_c08_MUC2",
#                      
#                      "A_Epi_c01_LEMD1", "A_Epi_c02_LDHB", "A_Epi_c03_ADGRF5",
#                      "A_Epi_c04_FXYD2", "A_Epi_c05_DMBT1", "A_Epi_c06_PGC",
#                      "A_Epi_c07_MKI67", "A_Epi_c08_NOS2", "A_Epi_c09_ZEB2",
#                      "A_Epi_c10_SPINK6", "A_Epi_c11_NEK10", "A_Epi_c12_ZBBX",
#                      "A_Epi_c13_WNT11","A_Epi_c14_FABP1",
#                      
#                      "Aw_Epi_c01_MOGAT1", "Aw_Epi_c02_ID3", "Aw_Epi_c03_CDC20", 
#                      "Aw_Epi_c04_CPE","Aw_Epi_c05_CXCL14","Aw_Epi_c06_SLC26A3",
#                      "Aw_Epi_c07_VGLL3", "Aw_Epi_c08_CALB1", "Aw_Epi_c09_KRT14",
#                      "Aw_Epi_c10_SULF1", "Aw_Epi_c11_DUSP27", "Aw_Epi_c12_CLCA1",
#                      "Aw_Epi_c13_KRT4",
#                      
#                      "O_Epi_c01_PFN1", "O_Epi_c02_TENM2", "O_Epi_c03_TM4SF4",
#                      "O_Epi_c04_CDC20", "O_Epi_c05_SPRR1B", "O_Epi_c06_CCDC190",
#                      "O_Epi_c07_VCAM1", "O_Ne_c01_NRXN3", "O_Ne_c02_ADGRB3",
#                      "O_Ne_c03_TRPM3", "O_Ne_c04_NPR3", "O_Ne_c05_HIST1H3F","O_Ne_c06_SCGB1D1"))
# all(unlist(uniq_fineAnno_list) %in% uniq_fineAnno)
ct_mismatch <- function(sr=srobj.aw,
                        anno_match_list=uniq_fineAnno_list
){
  sr.meta <- sr@meta.data
  mct_vec <- sort(unique(sr.meta$mainCellType))
  cell2rm <- list()
  for (mct in mct_vec) {
    meati <- sr.meta[sr.meta$mainCellType == mct, ]
    cell2rm[[mct]] <- meati[!(meati$fineAnno %in% anno_match_list[[mct]]),]
  }
  cell2rm <- do.call('rbind', cell2rm)
  print(paste0('Rm ', nrow(cell2rm),' cells'))
  print(table(cell2rm$mainCellType, cell2rm$fineAnno))
  #return(cell2rm)
  sr$c2r <- ifelse(sr$cellid %in% cell2rm$cellid, 'Y','N')
  print(paste0('Cell num in sr before filtering:  ', ncol(sr)))
  Idents(sr) <- 'c2r'
  sr <- subset(sr, c2r != 'Y')
  sr$c2r <- NULL
  Idents(sr) <- 'mainCellType'
  print(paste0('Cell num in sr after filtering:  ', ncol(sr)))
  return(sr)
}


ptsize <- function(sr){
  num <- ncol(sr)
  res <- case_when(
    num > 100000 ~ 1,    # 如果 x 是偶数，赋值 "Even"
    num > 20000  ~ 2,     # 如果 x 是奇数，赋值 "Odd"
    num > 10000  ~ 2.5,
    num > 5000   ~ 3,
    TRUE ~ 4
  )
  return(res)
}
