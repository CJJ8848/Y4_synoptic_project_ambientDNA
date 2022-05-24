unloadNamespace('SoupX')
unloadNamespace('Seurat')
unloadNamespace('sctransform')
library(Seurat)
library(dplyr)
library(patchwork)
library(SoupX)
####original soupx
#read in two matrix
tod <- Read10X("../synoptic_project_data/raw_feature_bc_matrix/")
toc <- Read10X("../synoptic_project_data/filtered_feature_bc_matrix/")
sc <- CreateSeuratObject(toc)
sc <- SCTransform(sc, verbose = F)
sc <- RunPCA(sc, verbose = F)
sc <- RunUMAP(sc, dims = 1:20, verbose = F)
sc <- FindNeighbors(sc, dims = 1:30, verbose = F)
sc <- FindClusters(sc, verbose = T,resolution = 0.4)
#create soupchannel
soup.channel  = SoupChannel(tod,toc)
meta    <- sc@meta.data
umap    <- sc@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
#run algorithm and generate soup gene expression profile
soup.channel  <- autoEstCont(soup.channel)
#show top 100 ambient DNA
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 100)
profie<-soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ]
#adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
#saveRDS(adj.matrix,"SoupX_corrected.rds")
