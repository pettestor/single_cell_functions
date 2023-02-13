
compute_silhouette_rogue_scores <- function(seurat,res.from=.1, res.to=.8, by=.05,plot=T){
  
  if(system.file(package="ROGUE") == ""){
    devtools::install_github("PaulingLiu/ROGUE")
  }
  
  library(ROGUE)
  library(cluster)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(Seurat)
  
  
  if(plot){
    dir.create("silhouette_rogue_plots")
  }
  
  #
  # Silhouette
  #
  
  cat(paste0("Using assay: ",DefaultAssay(seurat),"\nUpdate your object before applying the function if you want to use another one.\n"))
  
  
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
    
    # PLOTTING 
    if(plot){
      filename.sil=paste0("silhouette_rogue_plots/silhouette.res=",resolution,".png")
      png(filename.sil)
      plot(sil, border = NA)
      dev.off()
      cat("Plotted to: ",filename.sil,"\n")
      
      filename.rogue = paste0("silhouette_rogue_plots/rogue.res=",resolution,".png")
      rogue.bp <- rogue.boxplot(rogue.res)
      ggsave(plot = rogue.bp,filename = filename.rogue)
      cat("Plotted to: ",filename.rogue,"\n")
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
  
  ggsave(paste0("silhouette_rogue_plots/silhouette_results.pdf"),w=8,h=6)
  cat("Combined results plotted to silhouette_rogue_plots/silhouette_results.pdf.\n")
  
  
  return(sil.results)
}
