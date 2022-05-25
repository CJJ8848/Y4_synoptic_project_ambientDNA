
unloadNamespace('SoupX')
unloadNamespace('Seurat')
unloadNamespace('sctransform')
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DropletUtils")
library(DropletUtils)
library(Seurat)
library(dplyr)
library(patchwork)
library(SoupX)
library(ggplot2)
library(scales)

####ambient DNA expression profile
##read in raw and filtered matrix manually
args = commandArgs(trailingOnly=TRUE)
print ("raw matrix:")
print (args[1])
print ("filtered matrix:")
print (args[2])

#tod <- Read10X("../synoptic_project_data/raw_feature_bc_matrix/")
tod <- Read10X(args[1])
#toc <- Read10X("../synoptic_project_data/filtered_feature_bc_matrix/")
toc <- Read10X(args[2])
sc <- CreateSeuratObject(toc)
sc <- SCTransform(sc, verbose = F)
sc <- RunPCA(sc, verbose = F)
sc <- RunUMAP(sc, dims = 1:20, verbose = F)
sc <- FindNeighbors(sc, dims = 1:30, verbose = F)
sc <- FindClusters(sc, verbose = T,resolution = 0.4)

##empty droplet selection
#read in raw data
fname <- file.path( "../synoptic_project_data/raw_feature_bc_matrix/")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)
sce.pbmc
#generate cell information
bcrank <- barcodeRanks(counts(sce.pbmc))
#only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
#plot the UMI counts of all cells (before applying droplet filtration)
pdf('fig1_a_umicount_before.pdf')
df<-data.frame(rank=bcrank$rank[uniq], total=bcrank$total[uniq])
ggplot(df,aes(rank, total))+xlab("Rank")+ylab("Total UMI count")+geom_point(shape=1)+scale_x_log10()+scale_y_log10()+
  geom_hline(yintercept = metadata(bcrank)$inflection, colour="green")+
  geom_hline(yintercept =metadata(bcrank)$knee, colour="blue")
dev.off()
#legend=c("Inflection", "Knee")

set.seed(100)
#apply the algorithm
e.out <- emptyDrops(counts(sce.pbmc))
summary(e.out$FDR <= 0.001)
#select out the droplet containing a cell (FDR <=0.001)
retained <- sce.pbmc[,which(e.out$FDR <= 0.001)]
#new cell info
bcrank <- barcodeRanks(counts(retained))
uniq <- !duplicated(bcrank$rank)
#plot the UMI counts of all cells (after applying droplet filtration)
pdf('fig1_b_umicount_after.pdf')
df<-data.frame(rank=bcrank$rank[uniq], total=bcrank$total[uniq])
ggplot(df,aes(rank, total))+xlab("Rank")+ylab("Total UMI count")+geom_point(shape=1)+scale_x_log10()+scale_y_log10()
dev.off()
is.cell <- e.out$FDR <= 0.001
sum(is.cell, na.rm=TRUE)
#table(Limited=e.out$Limited, Significant=is.cell)
#visualization of the empty droplets
pdf('fig2_emptyvscell.pdf')
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")
dev.off()
#soupChannel create
soup.channel = SoupChannel(tod,toc,calcSoupProfile=FALSE)
soup.channel$Barcode <- c(attributes(soup.channel[["nDropUMIs"]]))
soup.channel$Barcode<-unlist(soup.channel$Barcode)
#Retain table of droplets and generate the soup gene expression profile
estimateSoupmine = function(sc,retained,keepDroplets=FALSE){
  if(!is(sc,'SoupChannel'))
    stop("sc must be a SoupChannel object.")
  #Estimate the soup 
  #w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  #w = which(sc$nDropUMIs > 0 & sc$nDropUMIs < 100)
  w = which((soup.channel$Barcode%in%retained@colData@rownames))
  sc$soupProfile = data.frame(row.names=rownames(sc$tod),
                              est = rowSums(sc$tod[,w,drop=FALSE])/sum(sc$tod[,w]),
                              counts = rowSums(sc$tod[,w,drop=FALSE]))
  #Saves a lot of space if we can drop the droplets now we're done with them
  if(!keepDroplets)
    sc$tod=NULL
  return(sc)
}
soup.channel = estimateSoupmine(soup.channel,retained=retained,keepDroplets=TRUE)
meta <- sc@meta.data
umap <- sc@reductions$umap@cell.embeddings
soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel <- setDR(soup.channel, umap)
soup.channel <- autoEstCont(soup.channel)
#show the top20 soup genes
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
#adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

#plot the top100 soup profile
#fig3_soup_profile top 5 of top100
soup_profile<-data.frame(head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 5))
p<-ggplot(soup_profile,aes(x=reorder(rownames(soup_profile),-soup_profile$est),y=soup_profile$est))+ylim(0,0.26)+
  geom_bar(stat = 'identity')+xlab("")+ylab("")+theme(axis.text.x = element_text(angle=90, hjust=.5, vjust=.5,size=24))+
  geom_text(aes(label=percent(round(soup_profile$est,3))),position=position_dodge(.5),vjust=-.5,size=8)
p

pdf('fig3_soup_profile_1.pdf')
p
dev.off()
#fig3_soup_profile bottom 5 of top100
soup_profile2<-data.frame(head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 100))
#write out top100 ambientDNAprofile
write.table(soup_profile2,"top100_ambientDNAprofile.txt")
soup_profile2<-data.frame(head(soup_profile2[order(soup_profile2$est, decreasing = F), ], n = 5))
p2<-ggplot(soup_profile2,aes(x=reorder(rownames(soup_profile2),-soup_profile2$est),y=soup_profile2$est))+ylim(0,0.26)+
  geom_bar(stat = 'identity')+xlab("")+ylab("")+theme(axis.text.x = element_text(angle=90, hjust=.5, vjust=.5,size=20))+
  geom_text(aes(label=percent(round(soup_profile2$est,5))),position=position_dodge(.5),vjust=-.5,size=8)
p2
pdf('fig3_soup_profile_2.pdf')
p2
dev.off()

##SingleCellExperiment 
##dim: 32285 387298 
##4386 filtered
##383416 empty
##3882 not empty



