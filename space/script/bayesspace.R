timestart<-Sys.time()

gwd <- getwd()
#get full path of this script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
ttdir <- dirname(thisFile())
print(ttdir)
setwd(ttdir)

library(BayesSpace)
library(SingleCellExperiment)
library(argparser)
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
library(patchwork)
library(assertthat)
library(cowplot)
library(viridis)
library(gtools)

argv <- arg_parser("")
argv <- add_argument(argv,"--datapath", help="the root of spatial data",default=NULL)
argv <- add_argument(argv,"--platform",help="Spatial sequencing platform. Used to determine spot layout and neighborhood structure (Visium = hex, ST = square).",default="Visium")
argv <- add_argument(argv,"--pca",help="Number of principal components to compute. We suggest using the top 15 PCs in most cases",default="10")
argv <- add_argument(argv,"--enhance",help="nrep of enhance",default="200")
argv <- add_argument(argv,"--clustnum",help="n cluster",default='0')
argv <- add_argument(argv,"--hvg",help="Number of highly variable genes to run PCA upon",default="2000")
argv <- add_argument(argv,"--qs",help="The values of q to evaluate",default="2,10")
argv <- add_argument(argv,"--initmethod",help="ues 'mclust' or 'kmeans' to obtain initial cluster assignments",default="kmeans")
argv <- add_argument(argv,"--scale",help="Controls the amount of jittering. Small amounts of jittering are more likely to be accepted but result in exploring the space more slowly. We suggest take values around 3",default="3.5")
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--prefix", help="group name prefix")
argv <- parse_args(argv)

###color
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
clustcol <- c(color_protocol, clustcol)
source("/Public/Script/shouhou/bayesspace/script/plotpicture.r")        ##绘图数据+热图
source("/Public/Script/shouhou/bayesspace/script/bayesspace_plot.R")    ##Enhance画图
###Variables
datapath <- argv$datapath
method <- argv$initmethod
platform <- argv$platform
qs <- as.character(unlist(strsplit(argv$qs,split = ",")))
clustnum <- as.numeric(argv$clustnum)
prefix <- argv$prefix
outdir <- argv$outdir
dir.create(outdir)

###loading data#####
sce <- readVisium(datapath)
seu_data <- Read10X(data.dir=file.path(datapath, "filtered_feature_bc_matrix"))

#######pre-processing data###
set.seed(519)
sce <- spatialPreprocess(sce,platform = platform,n.PCs = as.numeric(argv$pca),n.HVGs = as.numeric(argv$hvg),
	log.normalize = TRUE, 
	skip.PCA = FALSE, 
	assay.type = "logcounts"
)

###selecting the number of clusters###
sce <- qTune(sce,qs = seq(qs[1],qs[2]),platform = platform,d = as.numeric(argv$pca),
	nrep = 10000,
	burn.in = 1000
)

if(clustnum == '0'){
  sce_loglik <- attr(sce,"q.logliks") %>% as.matrix()
  diffsce <- diff(sce_loglik,lag=1,differences=1) %>% as.data.frame()
  sce_ck <- c(abs(diffsce$loglik)) %>% as.data.frame()
  colnames(sce_ck)[1] <- c("k")
  sce_loglik <- as.data.frame(sce_loglik)
  sce_loglik$sort <- rownames(sce_loglik)
  sce_ck$sort <- rownames(sce_ck)
  sce_loglik <- merge(sce_loglik,sce_ck,by="sort")
  q <- sce_loglik[which.min(sce_loglik$k),"q"]
}else{
  q = clustnum
}

p_q <- qPlot(sce)
pdf(paste0(outdir,'/',prefix,'.elbow_cluster.pdf'),width = 7,height = 7)
print(p_q)
dev.off()
png(paste0(outdir,'/',prefix,'.elbow_cluster.png'),width = 7,height = 7,res = 300,units = 'in')
print(p_q)
dev.off()

fun_drawpicture1(sce)

#####spatial cluster
if(platform == 'Visium'){
        gamma <- 3
}else if(platform == 'ST'){
        gamma <- 2
}else {
        stop("please input platform")
}
set.seed(519)
sce_cluster <- spatialCluster(sce,q = q,use.dimred = "PCA",platform = platform,d = as.numeric(argv$pca),init.method = method,gamma = gamma,
	model = "normal",
	nrep = 10000,
	burn.in = 1000,
	mu0 = NULL,
	lambda0 = NULL,
	alpha = 1,
	beta = 0.01,
	save.chain = FALSE,
	chain.fname = NULL
)

#####predict cluster for the cell-cell of stlearn
mdata <- sce_cluster@colData
pre_cluster <- c()
for(i in 1:length(unique(mdata$spatial.cluster))){
	pre_cluster <- c(pre_cluster,paste0("c",i))
}
tmp <- as.data.frame(matrix(NA,nrow=nrow(mdata),ncol = length(unique(mdata$spatial.cluster))))
for(i in 1:length(unique(mdata$spatial.cluster))){
tmp[i] <- ifelse(mdata$spatial.cluster==i,1,0)
}
for(i in 1:length(unique(mdata$spatial.cluster))){
	colnames(tmp)[i] = pre_cluster[i]
}
tmp$celltype_prediction <- paste0("c",mdata$spatial.cluster)
row.names(tmp) <- rownames(sce_cluster@colData)
write.csv(tmp,file=paste0(outdir,"/",prefix,".normal","_celltype_prediction.csv"),quote = FALSE,row.names = T)


######Seurat data###
print("################### Star  Seurat #############################")
seu_data <- CreateSeuratObject(counts = seu_data, assay = "Spatial", project = prefix, min.cells = 5)
image <- Read10X_Image(image.dir = file.path(datapath, "spatial"), filter.matrix = TRUE)
image <- image[Cells(x = seu_data)]
DefaultAssay(object = image) <- "Spatial"
seu_data[["slice1"]] <- image
rm(image)
mito.genes <- grep(pattern = "^MT-", x=rownames(x=seu_data[["Spatial"]]@data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.mito"]] <- PercentageFeatureSet(seu_data, features = mito.genes)
Ribosomal <- grep(pattern = "^(RPL|RPS)", x = rownames(x = seu_data[["Spatial"]]@data), value = TRUE, ignore.case = TRUE)
seu_data[["percent.Ribo"]] <- PercentageFeatureSet(seu_data, features = Ribosomal)
mito.distrib <- function(sample, x){
  minors <- 0
  for(n in sample){
    if(n > x){
      minors <- minors + 1
    }
  }
  return (minors / length(sample))
}
percent.mito <- Matrix::colSums(seu_data[["Spatial"]]@counts[mito.genes,]) / Matrix::colSums(seu_data[["Spatial"]]@counts)
u <- as.data.frame(percent.mito)
mito <- u$percent.mito
p05 <- mito.distrib(mito, 0.1) * 100
p1 <- mito.distrib(mito, 0.2) * 100
p15 <- mito.distrib(mito, 0.3) * 100
p2 <- mito.distrib(mito, 0.5) * 100
percent_mito <- c(">10%", ">20%", ">30%", ">50%")
percent_spot_num <- c(p05, p1, p15, p2)
data.mito <- data.frame(mito_percent = percent_mito, spot_percent = percent_spot_num)
################
seu_data <- NormalizeData(seu_data, assay = "Spatial", verbose = FALSE)
all_genes <- rownames(seu_data@assays$Spatial@data)
seu_data <- ScaleData(seu_data, assay = "Spatial", verbose = FALSE, features = all_genes)
seu_data  <- FindVariableFeatures(seu_data, selection.method = "vst", nfeatures = 2000)
seu_data <- RunPCA(seu_data, assay = "Spatial", verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5)
seu_data <- FindNeighbors(seu_data, reduction = "pca", dims = 1:as.numeric(argv$pca))
seu_data <- FindClusters(seu_data, verbose = FALSE, resolution = 0.8)
seu_data <- RunUMAP(seu_data, reduction = "pca", dims = 1:as.numeric(argv$pca), verbose = FALSE)
#seu_data <- RunTSNE(seu_data, reduction = "pca", dims = 1:as.numeric(argv$pca))

####sce + seu###
bsclust <- data.frame(rownames(sce_cluster@colData))
bsclust$bscluster <- sce_cluster@colData$spatial.cluster
colnames(bsclust) <- c('barcode','bscluster')
spclust <- data.frame(rownames(seu_data@meta.data))
colnames(spclust) <- c('barcode')
spclust <- left_join(spclust,bsclust,by='barcode')
seu_data@meta.data$bscluster <- spclust$bscluster
seu_data <- SetIdent(seu_data,value=seu_data@meta.data$bscluster)
clust <- summary(seu_data@active.ident)
cluster_spot <- as.data.frame(clust)
spot_number <- sum(cluster_spot$clust)
print(spot_number)
pt_use <- 0.6
if(spot_number > 1000){
  pt_use <- 0.4
}
if(spot_number > 2500){
  pt_use <- 0.3
}
if(spot_number > 4000){
  pt_use <- 0.2
}
if(spot_number > 5500){
  pt_use <- 0.15
}
if(spot_number > 6500){
  pt_use <- 0.1
}


write.table(sce_cluster@colData,file = paste0(outdir,'/',prefix,'_normal_colData.xls'),sep='\t',quote = F,row.names=F)
saveRDS(sce_cluster,file = paste0(outdir,'/',prefix,'.BayesSpace_normal.rds'))
saveRDS(seu_data,file=paste0(outdir,'/',prefix,'_normal.rds'))

cluster_order = mixedsort(levels(seu_data))
color_order = clustcol[1:length(cluster_order)]
cols <- list()
for(i in 1:length(cluster_order)) {
  cols[cluster_order[i]] <- color_order[i]
}
p_cluster <- SpatialDimPlot(seu_data,group.by='bscluster', cols = cols) + labs(title = prefix) + theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(outdir,'/',prefix,'.spatial_cluster.pdf'),width = 7,height = 7)
print(p_cluster)
dev.off()
png(paste0(outdir,'/',prefix,'.spatial_cluster.png'),width = 7,height = 7,res = 300,units = 'in')
print(p_cluster)
dev.off()



#fun_drawpicture2(p_cluster)
#####寻找差异基因####
diffgene_dir <- file.path(outdir, "diffgene")
dir.create(diffgene_dir)
j <- length(levels(x = seu_data))
for (l in 1: j){
  cluster.markers <- FindMarkers(object = seu_data, ident.1 = l, min.pct = 0.1, logfc.threshold = 0.25)
  colnames(cluster.markers)[2] <- c("avg_logFC")
  cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE),]
  write.table(cluster.markers, file = paste0(diffgene_dir, '/', prefix, '_cluster', l, '_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = TRUE)
}
seu_data.marker <- FindAllMarkers(object = seu_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
colnames(seu_data.marker)[2] <- c("avg_logFC")
markergene <- seu_data.marker %>% group_by(cluster) %>% top_n(2, avg_logFC)
markergenetop10 <- seu_data.marker %>% group_by(cluster) %>% top_n(10, avg_logFC)

heatname <- paste0(prefix, ' top10 marker genes heatmap')
pdf(paste0(diffgene_dir, '/', prefix, '_diffgenetop10_DOHeatmapplot.pdf'))
DoHeatmap(seu_data, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 2)) + scale_fill_viridis()
dev.off()
png(paste0(diffgene_dir, '/', prefix, '_diffgenetop10_DOHeatmapplot.png'))
DoHeatmap(seu_data, features = markergenetop10$gene, size = 3, hjust = 0, angle = 0, slot = "scale.data", draw.lines = TRUE, raster=F, group.by= "ident", assay = "Spatial") + theme(axis.text.y = element_text(size = 2)) + scale_fill_viridis()
dev.off()
allmarker <- as.data.frame(markergene)
write.table(allmarker, file = paste0(diffgene_dir, '/', prefix, '_all_diffgenes.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
fun_drawpicture3(markergenetop10,seu_data)

#####enhance
sce_enhance <- spatialEnhance(sce,q = q,platform = platform,d = as.numeric(argv$pca),init.method = method,gamma = gamma,jitter_scale = as.numeric(argv$scale),
	use.dimred = "PCA",
	model = "normal",
	nrep = as.numeric(argv$enhance),
	burn.in = as.numeric(argv$enhance)/10,
	save.chain = FALSE,
	jitter_prior = 0.3
)
p_enhanced <- clusterplot(sce_enhance,palette=c(color_protocol, clustcol)) +
	theme_bw() +
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
	ylab("Row") +
	xlab("Column") +
	labs(fill = "BayesSpace\ncluster",title="Spatial clustering Enhanced")
pdf(paste0(outdir,'/',prefix,'.Enhanced_spatial_cluster.pdf'),width = 7,height = 7)
print(p_enhanced)
dev.off()
png(paste0(outdir,'/',prefix,'.Enhanced_spatial_cluster.png'),width = 7,height = 7,res = 300,units = 'in')
print(p_enhanced)
dev.off()

write.table(sce_enhance@colData,file = paste0(outdir,'/',prefix,'_enhanced_colData.xls'),sep='\t',quote = F,row.names=F)
saveRDS(sce_enhance,file = paste0(outdir,'/',prefix,'.BayesSpace_enhanced.rds'))

timeend <- Sys.time()
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Programe running time")
print(timeend-timestart)


file.copy("/Public/Script/shouhou/SCRIPT/BayesSpace/BayesSpace-README.pdf", outdir)

