library(tidyverse)
library(ggplot2)



#mk <- data.frame(CELLTYPE= c('aSMC','aSMC','aSMC', 'aaSMC','aaSMC'),
#                 gene = c( "ACTA2", "TAGLN","SLIT3", "SLIT3", "CTNNA3"))
#plotDot2 is used for well annotated data's marker genes.
plotDot2 <- function(expr_mat = t1,
                     markerDf = mk,
                     min_expr = 0,
                     minPercent = 1,
                     zscore.min = NULL,
                     zscore.max = NULL,
                     returnData=FALSE,
                     returnPlot=TRUE){
  #expr_met -> a matrix or a seurat obj.
  #min_expr -> the threshold for counting expresion percent.
  #minPercent -> the min percent to plot a dot.
  #mk -> a df. a df with col 'CELLTYPE' 'gene'
  
  #if expr_meta is sr obj than using rna data in markerDf
  if(!(is_tibble(expr_mat)|is.data.frame(expr_mat)|is.matrix(expr_mat))){
    print('suppose expr_mat is a srobj')
    print('using Idents in expr_mat as levels')
    
    #t1 <- as.data.frame(t(as.matrix(expr_mat@assays$RNA@data[rownames(expr_mat) %in% unique(unlist(markerDf$gene)),])))
    t1 <- as.data.frame(FetchData(object = expr_mat, vars = unique(unlist(markerDf$gene))),row.names = Cells(expr_mat))
    #t1$label <- Idents(expr_mat)
    t1$label <- Idents(expr_mat)[names(Idents(expr_mat)) %in% rownames(t1)]
    t1$label <- factor(x = t1$label, levels = levels(expr_mat))
    expr_mat <- t1
  }
  
  #celltype levels
  if(all(is.factor(expr_mat$label))){
    ct_level <- levels(expr_mat$label)
  }else{
    ct_level <- sort(unique(unlist(expr_mat$label)))
  }
  expr_percent <- expr_mat %>%
    group_by(label) %>%
    summarise(across(everything(), .fns = function(x){sum(x > min_expr)/length(x)*100})) %>%
    column_to_rownames(var = 'label') %>% t()
  
  expr <- expr_mat
  
  sub_expr <- expr %>% 
    group_by(label) %>%
    summarise(across(everything(),mean))  %>%
    mutate(across(-label, ~ ((.x - mean(.x))/sd(.x)))) %>%
    column_to_rownames(var = 'label') %>% t()
  sub_expr <- as.data.frame(sub_expr)
  
  
  gene <- rownames(sub_expr)
  sub_expr <- sub_expr %>%
    dplyr::mutate(gene = gene) %>%
    pivot_longer(!gene, names_to = 'cell_type', values_to = 'zscore') %>%
    dplyr::mutate(id= paste0(gene, cell_type))
  
  expr_percent <- as.data.frame(expr_percent)
  expr_percent <- expr_percent %>%
    dplyr::mutate(gene = gene) %>%
    pivot_longer(!gene, names_to = 'cell_type', values_to = 'expr_percent') %>%
    dplyr::mutate(id= paste0(gene, cell_type))
  
  sub_expr <- left_join(x = sub_expr, y = expr_percent[,c('id', 'expr_percent')])
  sub_expr <- dplyr::select(.data = sub_expr, !id)
  
  
  mk <- markerDf
  if(all(is.factor(mk$CELLTYPE))){
    CT_level <- levels(mk$CELLTYPE)
    mk$CELLTYPE <- as.character(mk$CELLTYPE)
  }
  mk <- left_join(x = mk, y = sub_expr, by='gene')
  
  mk <- mk %>% dplyr::filter(expr_percent > minPercent)
  
  mk$cell_type <- factor(x = mk$cell_type, levels = ct_level)
  if(exists('CT_level')){
    mk$CELLTYPE <- factor(x = mk$CELLTYPE, levels = CT_level)
  }
  if(returnData){
    return(mk)
  }
  
  if(is.null(zscore.min)){
    zscore.min <- min(mk$zscore)
  }
  if(is.null(zscore.max)){
    zscore.max <- max(mk$zscore)
  }
  
  p1 <- mk %>% ggplot(aes(factor(gene, levels = unique(sub_expr$gene)), 
                          cell_type)) +
    geom_point(aes(size = expr_percent, colour = zscore)) +
    theme(strip.text.x = element_text(face = 'bold', colour = '#2b8cbe', size = rel(1.5)),
          axis.title = element_text(size = rel(1.5)),
          axis.text = element_text(size = rel(1.2)),
          legend.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1.2)), 
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black", angle = -90, vjust = 0.5,hjust = 0),
          panel.background = element_rect(colour = "black", fill = "white"), 
          panel.grid = element_line(colour = "grey", linetype = "dashed"), 
          panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)) + 
    facet_grid(. ~ CELLTYPE, scales = "free", space = "free") +
    scale_colour_distiller(palette = "RdYlBu",
                           limits = c(zscore.min, zscore.max),
                           oob = scales::squish) + 
    labs(x = "", y = "")
  
  if(returnPlot){
    return(p1) 
  }
}

#a fun to auto-classify celltype by expr.
roughCelltype <- function(srobj,
                          mker,
                          confidence=c(0.2,0.5),
                          minDE=0.2,
                          uncertainCellName='Uncertain'){
  #srobj = TNKCell
  #mker = mklist
  #confidence=c(0.2,0.5)
  #minDE=0.2
  #uncertainCellName='Uncertain'
  
  #define a func to get confidence.
  getConfi <- function(vec,confidence){
    return(ifelse(vec<0.2, 'low', ifelse(vec>0.5, 'high', 'middle')))
  }
  
  
  if(!(is.data.frame(mker)|is_tibble(mker))){
    mkdf <- namedList2LongDf(mk_list = mker)
  }
  mk <- plotDot2(expr_mat = srobj,
                 markerDf = mkdf,
                 min_expr = 0,
                 minPercent = 0,
                 zscore.min = NULL,
                 zscore.max = NULL,
                 returnData = T)
  #celltype levels
  if(all(is.factor(mk$cell_type))){
    ct_level <- levels(mk$cell_type)
  }else{
    ct_level <- sort(unique(unlist(mk$cell_type)))
  }
  if(all(is.factor(mk$CELLTYPE))){
    CT_level <- levels(mk$CELLTYPE)
  }else{
    CT_level <- sort(unique(unlist(mk$CELLTYPE)))
  }
  
  mk.mean <- mk %>% select(CELLTYPE, cell_type, zscore, expr_percent) %>% 
    group_by(CELLTYPE, cell_type) %>% summarise(across(.cols = everything(), mean))
  res_l <- vector(mode = 'list',length = length(ct_level))
  for (i in 1:length(ct_level)) {
    #i <- 2
    
    dfi <- mk.mean[mk.mean$cell_type == ct_level[i],]
    #dfi <- dfi %>% mutate(CELLTYPE=as.character(CELLTYPE), cell_type=as.character(cell_type))
    if(nrow(dfi) < 1){
      resDfi <- data.frame(CELLTYPE=uncertainCellName,
                           cell_type=ct_level[i],
                           confidence='High')
    }else{
      zs <- sort(unlist(dfi$zscore),decreasing = T)
      resDfi <- dfi[1,1:2]
      resDfi$confidence <- 'High'
      if(length(zs) < 2){
        resDfi$confidence <- getConfi(vec = dfi$expr_percent,confidence = confidence)
      }else{
        ValDelta <- (zs[1] -zs[2])/abs(zs[2])
        if(ValDelta>=minDE){
          idx <- which.max(dfi$zscore)[1]
          resDfi <- dfi[idx, 1:2]
          resDfi$confidence <- getConfi(vec = dfi$expr_percent[idx],confidence = confidence)
        }else{
          resDfi$CELLTYPE <- uncertainCellName
        }
      }
    }
    res_l[[i]] <- resDfi
  }
  res_l <- do.call('rbind', res_l)
  return(res_l)
}
