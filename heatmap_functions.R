getCutHeight <- function(mod,search.seq, n.clusters){
  
  cuth <-  sapply(search.seq,function(x){
    length(unique(cutree(mod,h=x)))
  })
  cuth <- cuth[cuth>=n.clusters]
  hld <- which.min(abs(n.clusters-cuth))
  cat("Cutting into ",cuth[hld]," clusters\n")
  cuth <- search.seq[hld]
  
  return(cuth)
}


makeHeatmap <- function(data,OIDs,color.vars=NULL,color.labels,
                        filename,out.width=20,out.height=20,
                        truncate.sd=4,scale=T,heatmap.color.scale=rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                                                         "#4393C3", "#2166AC", "#053061"))(1024)),
                        dist.func.r=function(x) as.dist(1-cor(t(x),method="pearson")),
                        dist.func.c=function(x) as.dist(1-cor(t(x),method="pearson")),
                        cut.search,
                        cut.n,
                        value="NPX",
                        return.cluster=F,
                        row.labels,row_cex,
                        add_gridlines=F,
                        return.row.cluster=F,
                        cut.n.row
                        ){
  
  if(!is.null(color.vars)){
    if(missing(color.labels)) color.labels <- color.vars
    names(color.labels) <- color.vars
  }
  
  data[["NPX"]] <- data[[value]]
  
  if(scale){
    
    tmp.wide <- data %>% 
      filter(OlinkID %in% OIDs ) %>%
      mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
      select(all_of(c("Assay_OlinkID","SampleID","NPX",color.vars))) %>% 
      group_by(Assay_OlinkID) %>% 
      mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>%
      mutate(NPX= ifelse(abs(NPX)>truncate.sd,truncate.sd*sign(NPX),NPX)) %>% ##Replace any values beyond truncate.sdsd with truncate.sd
      mutate(NPX= ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##Replace any na's with 0/mean
      ungroup() %>% 
      pivot_wider(names_from="Assay_OlinkID",values_from="NPX")
    
    
  } else if(!scale){
    tmp.wide <- data %>% 
      filter(OlinkID %in% OIDs ) %>%
      mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
      select(all_of(c("Assay_OlinkID","SampleID","NPX",color.vars))) %>% 
      group_by(Assay_OlinkID) %>% 
      mutate(NPX= ifelse(is.na(NPX)|!is.finite(NPX),median(NPX,na.rm=T),NPX)) %>% ##Replace any na's with 0/mean
      ungroup() %>% 
      pivot_wider(names_from="Assay_OlinkID",values_from="NPX")
    
    
  }
  
  if(!missing(color.vars)){
    
    for(i in color.vars) tmp.wide[[i]] <- factor(tmp.wide[[i]])
    color.pal <- c("#00C7E1FF", "#FE1F04FF", "#00559EFF", "#FFC700FF","#077183FF","#FF51B8FF","#27AE55FF","#6A27AEFF","#FF8C22FF","#AAB1B9FF","red","blue")
    
    pal.iter <- 0
    for(i in color.vars){
      tmp.wide[[paste0(i,".color")]] <- color.pal[(pal.iter+1):(pal.iter+length(levels(tmp.wide[[i]])))][as.numeric(tmp.wide[[i]])]
      pal.iter <- pal.iter + length(levels(tmp.wide[[i]]))
    }
    
    rowside.cols <- as.matrix(tmp.wide[,paste0(color.vars,".color")])
    colnames(rowside.cols) <- color.labels
    x <- as.matrix(tmp.wide %>% select(!all_of(c("SampleID",color.vars,paste0(color.vars,".color")))))
    rownames(x) <- tmp.wide$SampleID
    
    
    leg.labs <- NULL
    leg.cols <- NULL
    for(i in color.vars){
      tmp <- tmp.wide %>% select(all_of(c(i,paste0(i,".color")))) %>% unique() %>% arrange(!!rlang::ensym(i))
      leg.labs <- c(leg.labs,as.character(color.labels[i]),as.character(tmp[[i]]))
      leg.cols <- c(leg.cols,"white",tmp[[paste0(i,".color")]])
    }
    
    
  } else{
    x <- as.matrix(tmp.wide %>% select(!all_of(c("SampleID"))))
    rownames(x) <- tmp.wide$SampleID
    rowside.cols <- NULL
    
    leg.labs <- " "
    leg.cols <- "white"
    
  }
  
  
  if(!missing(cut.search) & !missing(cut.n)){
    tmp.mod <- hclust(dist.func.c((x)))
    cuth <- getCutHeight(tmp.mod,search.seq = cut.search,n.clusters=cut.n)
  } else{
    cuth <- NULL
  }

  
  min.col <- which.min(min(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  max.col <- which.min(max(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  
  
  CairoPDF(filename,
           width=out.width,
           height=out.height)
  
  
  args <- list(x=t(x),
               scale="none",
               method="complete",
               distfunR = dist.func.r,
               distfunC = dist.func.c,
               col = heatmap.color.scale[min.col:max.col],
               legendfun =  function() {showLegend(legend =leg.labs,
                                                   col=leg.cols,
                                                   cex=2.5,
                                                   lwd=5)})
  if(!missing(row_cex)){
    args$cexRow <- row_cex
  }
  if(!is.null(cuth)) args$ColSideCut <- cuth
  if(!is.null(rowside.cols)) args$ColSideColors = rowside.cols
  if(!missing(row.labels)){
    tmp.row.labels <- rep("",nrow(t(x)))
    tmp.row.labels[rownames(t(x)) %in% row.labels] <- rownames(t(x))[rownames(t(x)) %in% row.labels]
    args$labRow <- tmp.row.labels
  }
  if(add_gridlines){
    args$highlightCell = expand.grid(r=1:nrow(t(x)),c=1:ncol(t(x)),col="black",sz=1)
  }
  
  
  htmap <- do.call(heatmap3,args=args)

    ##Color bar
  plotfunctions::gradientLegend(valRange=range(x,na.rm = T),
                                color = heatmap.color.scale,
                                pos=c(0.125,.85,0.135,1.02), # xleft, ybottom, xright, ytop
                                fit.margin = F,cex=1.75)
  
  dev.off()
  
  #return col cluster only
  if(return.cluster){
    if(is.null(cuth)) stop("Must provide search sequence and desired number of clusters")
    out.clusters <- cutree(tmp.mod,h=cuth)
    cluster.list.col <- data.frame(SampleID = names(out.clusters),
                               Cluster=as.character(as.numeric(out.clusters))) %>% 
      left_join(data.frame(SampleID=colnames(t(x))[htmap$colInd],
                           Order_L2R = 1:length(htmap$colInd)),by="SampleID") %>% 
      arrange(Order_L2R) %>% 
      group_by(Cluster) %>% 
      mutate(med = median(Order_L2R)) %>% 
      ungroup() %>% 
      mutate(Cluster=as.character(as.numeric(as.factor(med)))) %>% 
      select(-med)
    
    
  }
  
  ##return row cluster only
  if(return.row.cluster){
    if(missing(cut.n.row)) stop("Must provide search sequence and desired number of clusters")
    
    
    tmp.mod <- hclust(dist.func.r(t(x)))
    # cuth <- getCutHeight(tmp.mod,search.seq = cut.search,n.clusters=cut.n)
    out.clusters <- cutree(tmp.mod,k=cut.n.row)
    
    
    cluster.list.row <- data.frame(Assay = names(out.clusters),
                               Cluster=as.character(as.numeric(out.clusters))) %>% 
      left_join(data.frame(Assay=colnames(t(t(x)))[htmap$rowInd],
                           Order_T2B = length(htmap$rowInd):1),by="Assay") %>% 
      arrange(Order_T2B) %>% 
      group_by(Cluster) %>% 
      mutate(med = median(Order_T2B)) %>% 
      ungroup() %>% 
      mutate(Cluster=as.character(as.numeric(as.factor(med)))) %>% 
      select(-med)

  }
  
  if(return.cluster & return.row.cluster){
    return(list(col.clusters = cluster.list.col,
                row.clusters = cluster.list.row))
  }
  if(return.cluster & !return.row.cluster){
    return(cluster.list.col)
  }
  if(!return.cluster & return.row.cluster){
    return(cluster.list.row)
  }
  
}



###### Assay correlation heatmaps -----------------

makeAssayCorrHeatmap <- function(data,OIDs,
                                 filename,out.width=20,out.height=20,
                                 truncate.sd=4,scale=T,heatmap.color.scale=rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                                                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                                                                  "#4393C3", "#2166AC", "#053061"))(1024)),
                                 cut.search,
                                 cut.n,
                                 value="NPX",
                                 method="pearson",
                                 return.cluster=F,
                                 row.labels,
                                 row_cex,
                                 add_gridlines=F){
  
  data[["NPX"]] <- data[[value]]
  
  
  tmp <- data %>% 
    filter(OlinkID %in% OIDs) %>% 
    mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
    select(Assay_OlinkID,SampleID,NPX) %>% 
    group_by(Assay_OlinkID) %>% 
    mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
    mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>%
    ungroup() %>% 
    pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
    select(-SampleID) %>% 
    cor(method=method)
  
  min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
  max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))
  
  if(!missing(cut.search) & !missing(cut.n)){
    tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
    cuth <- getCutHeight(tmp.mod,search.seq = cut.search,n.clusters=cut.n)
  } else{
    cuth <- NULL
  }
  
  
  CairoPDF(filename,
           width=out.width,
           height=out.height)
  
  args <- list(x=tmp,
           scale="none",
           method="complete",
           revC=T,
           distfun = function(x) as.dist(1 - (x)),
           col = heatmap.color.scale[min.col:max.col],
           legendfun =  function() {showLegend(legend = "",col="white")},
           margins=c(10,10))
  
  
  
  if(!is.null(cuth)) args$ColSideCut <- cuth
  if(!missing(row.labels)){
    tmp.row.labels <- rep(NA,nrow(tmp))
    names(tmp.row.labels) <- rownames(tmp)
    tmp.row.labels[row.names(tmp) %in% row.labels] <- rownames(tmp)[row.names(tmp) %in% row.labels]
    args$labRow <- tmp.row.labels
  }
  
  if(!missing(row_cex)){
    args$cexRow <- row_cex
  }
  if(add_gridlines){
    args$highlightCell = expand.grid(r=1:nrow(tmp),c=1:ncol(tmp),col="black",sz=1)
  }
  
  
  
  original.yaxp <- par("yaxp")
  par(yaxp=c(1,nrow(tmp),nrow(tmp)))
  
  htmap <- do.call(heatmap3,args=args)
  

  plotfunctions::gradientLegend(valRange=c(-1,1),
                                color = heatmap.color.scale,
                                pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                                fit.margin = F,cex=2)
  
  
  dev.off()
  
  par(yaxp=original.yaxp)

  if(return.cluster){
    if(is.null(cuth)) stop("Must provide search sequence and desired number of clusters")
    out.clusters <- cutree(tmp.mod,h=cuth)
    cluster.list <- data.frame(Assay = names(out.clusters),
                               Cluster=as.character(as.numeric(out.clusters))) %>% 
      left_join(data.frame(Assay=colnames(tmp)[htmap$colInd],
                           Order_L2R = 1:length(htmap$colInd)),by="Assay") %>% 
      arrange(Order_L2R) %>% 
      group_by(Cluster) %>% 
      mutate(med = median(Order_L2R)) %>% 
      ungroup() %>% 
      mutate(Cluster=as.character(as.numeric(as.factor(med)))) %>% 
      select(-med)
    

    return(cluster.list)
  }
  
  
  
}


