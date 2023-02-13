
compute_silhouette_rogue_scores <- function(seurat,res.from=.1, res.to=.8, by=.05,plot=T, plot.prefix=""){
  
  if(system.file(package="ROGUE") == ""){
    devtools::install_github("PaulingLiu/ROGUE")
  }
  
  library(ROGUE)
  library(cluster)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(Seurat)
  
  #
  # Silhouette
  #
  
  cat(paste0("Using assay: ",DefaultAssay(seurat),"\nUpdate your object before applying the function if you want to use another one."))
  
  
  expr <- as.matrix(GetAssayData(seurat))
  
  
  sil.results <- data.frame()
  
  for(resolution in seq(res.from,res.to,by = by)){
    
    # SILHOUETTE
    seurat<-FindClusters(seurat,resolution=resolution,verbose = F)
    dists <- dist(Embeddings(seurat,reduction = "umap"))
    combined.fac <- factor(paste0(seurat$seurat_clusters))
    clusters<-as.integer(seurat$seurat_clusters)
    sil <- silhouette(clusters, dist = dists)
   
    # ROGUE
    rogue.res <- rogue(expr, labels = clusters, samples = seurat$orig.ident, platform = "UMI", span = 1.6)
    
    if(plot){
      filename=paste0(plot.prefix,"silhouette.res=",resolution,".png")
      png(filename)
      plot(sil, border = NA)
      dev.off()
      
      rogue.bp <- rogue.boxplot(rogue.res)
      ggsave(plot = rogue.bp,filename = paste0(plot.prefix,"rogue.res=",resolution,".png"))
      cat("Plotted to: ",filename,"\n")
    }
    
    
    sil.results <- rbind(sil.results,c(resolution,mean(sil[,3]),length(levels(seurat$seurat_clusters)),median(reshape2::melt(rogue.res)[,2],na.rm = T)))
    
    
  }
  
  colnames(sil.results)<-c("Resolution","Average Si","nClusters", "Average ROGUE")
  
  plot_grid(
    ggplot(sil.results,aes(x=factor(Resolution),y=`Average Si`))+geom_bar(stat="identity")+theme_cowplot()+xlab("")+geom_bar(data=subset(sil.results, `Average Si`==max(`Average Si`)), aes(factor(Resolution),`Average Si`),
                                                                                                                             fill="red", stat="identity"),
    ggplot(sil.results,aes(x=factor(Resolution),y=`Average ROGUE`))+geom_bar(stat="identity")+theme_cowplot()+xlab("")+
      geom_bar(data=subset(sil.results, `Average ROGUE`==max(`Average ROGUE`)), aes(factor(Resolution),`Average ROGUE`),fill="red", stat="identity")+ coord_cartesian(ylim=c(min(sil.results$`Average ROGUE`)*.9,1)),
    ggplot(sil.results,aes(x=factor(Resolution),y=nClusters))+geom_bar(stat="identity")+theme_cowplot()+xlab("Resolution"),
    ncol=1)
  
  ggsave(paste0(plot.prefix,"silhouette_results.pdf"),w=8,h=6)
  cat("Combined results plotted to silhouette_results.pdf.\n")
  
  #
  # ROGUE
  #
  

  
  
  
  return(sil.results)
}
