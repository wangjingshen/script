# basic_plot

suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
Total_genus_threshold = 2

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--subcluster", help = "subcluster, default: all")
argv <- add_argument(argv, "--splitgroup", help = "splitgroup, default: T")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 01.basic_plot")
argv <- parse_args(argv)

rds <- argv$rds
subcluster <- ifelse(is.na(argv$subcluster), "all", argv$subcluster)
splitgroup <- ifelse(is.na(argv$splitgroup), "T", argv$splitgroup)
outdir <- ifelse(is.na(argv$outdir), "01.basic_plot", argv$outdir)


if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}


data_seurat <- readRDS(rds)
genus <- colnames(data_seurat@meta.data)[(which(colnames(data_seurat@meta.data) == "total_genus") + 1) : ncol(data_seurat@meta.data)]

# totol genus umi: cluster group --
options(repr.plot.height = 5, repr.plot.width =7)
data_seurat@meta.data %>% 
    group_by(cluster, group) %>%
    summarise(total_genus_umi = sum(total_genus)) %>%
    ggplot(aes(x=cluster, y=total_genus_umi, fill=group)) +
        geom_bar(color="black", stat = "identity", position = position_dodge()) +
        theme_classic() + 
        scale_fill_manual(values = color_protocol) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
ggsave(str_glue("{outdir}/barplot_total_genus_umi.png"), height = 5, width = 7)

# totol genus umi group --
options(repr.plot.height = 5, repr.plot.width =4)
data_seurat@meta.data %>% 
    group_by(group) %>%
    summarise(total_genus_umi = sum(total_genus)) %>%
    ggplot(aes(x=group, y=total_genus_umi, fill=group)) +
        geom_bar(color="black", stat = "identity", position = position_dodge(), width = 0.5) +
        theme_classic() + 
        scale_fill_manual(values = color_protocol) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
ggsave(str_glue("{outdir}/barplot_total_genus_umi_group.png"), height = 5, width = 4)


# DimPlot genus_detect
data_seurat$genus_detect <- factor(ifelse(data_seurat$total_genus > Total_genus_threshold, "Genus_detected", "No_genus_detected"), 
                                   levels = c("No_genus_detected","Genus_detected"))

DimPlot(data_seurat, group.by = "genus_detect", cols = c("lightgrey", "red"), pt.size = 0.2, order = T)
ggsave(str_glue("{outdir}/featureplot_genus_detect.png"), height = 5, width = 7)

DimPlot(data_seurat, group.by = "genus_detect", cols = c("lightgrey", "red"), pt.size = 0.2, order = T, split.by = "group")
ggsave(str_glue("{outdir}/featureplot_genus_detect_group.png"), height = 5, width = 13)


# sample cluster genus_count genus_umi
if(subcluster == "all"){
    subcluster = levels(data_seurat$cluster)
    com <- expand.grid(levels(data_seurat$cluster), levels(data_seurat$sample))
}else{
    subcluster = unlist(str_split(subcluster, ","))
    com <- expand.grid(subcluster, levels(data_seurat$sample))
}
print(subcluster)
# get df
df_stat <- do.call(rbind,lapply(unique(data_seurat$sample), function(sample){
    data <-  data_seurat@meta.data[ data_seurat$sample == sample,]
    
    df <- do.call(rbind, lapply(subcluster, function(cluster){
        #print(cluster)
        data_sub <- data[ data$cluster == cluster,]
        data.frame(cluster = cluster,
                   genus_count = sum(colSums(data_sub[,genus]) > 0),  # 种类
                   genus_umi = sum(data_sub[,genus]))  # umi
    }))
    df$sample = sample
    return(df)
}))

com_levels <- paste(com$Var2, com$Var1, sep = "_")

df_stat$sample_cluster <- factor(paste0(df_stat$sample, "_" ,df_stat$cluster), levels = com_levels)
width = 5+nrow(df_stat)/84*15

# plot
## genus_count
df_stat %>% ggplot(aes(x= sample_cluster, y = genus_count, fill=sample)) +
    geom_bar(stat = "identity",color="black") +
        scale_fill_manual(values = color_protocol) +
        theme_classic() +
        theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 90), legend.position = "none", axis.title.x = element_blank())
ggsave(str_glue("{outdir}/barplot_genus_count.pdf"), height = 7, width = width)
ggsave(str_glue("{outdir}/barplot_genus_count.png"), height = 7, width = width)

## genus_umi
df_stat %>% ggplot(aes(x= sample_cluster, y = genus_umi, fill=sample)) +
    geom_bar(stat = "identity",color="black") +
        scale_fill_manual(values = color_protocol) +
        theme_classic() +
        theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 90), legend.position = "none", axis.title.x = element_blank())
ggsave(str_glue("{outdir}/barplot_genus_umi.pdf"), height = 7, width = width)
ggsave(str_glue("{outdir}/barplot_genus_umi.png"), height = 7, width = width)

write.table(df_stat[,c(4,1:3)], str_glue("{outdir}/barplot_plot_data.tsv"), sep="\t", row.names=F, quote =F)  # del sample_cluster


## top10 genus
if(!dir.exists(str_glue("{outdir}/top10_genus_split/"))){
    dir.create(str_glue("{outdir}/top10_genus_split/"), recursive = T)
}

function_top10_genus_plot <- function(data, name){
    top10_genus <- names(sort(colSums(data@meta.data[, genus]>0), decreasing = T)[1:10])
    FeaturePlot(data, features = top10_genus, order = T)
    ggsave(str_glue("{outdir}/Featureplot_top10_genus_{name}.pdf"), height = 12, width = 18)
    ggsave(str_glue("{outdir}/Featureplot_top10_genus_{name}.png"), height = 12, width = 18)
    
    for (i in 1:length(top10_genus)){
        FeaturePlot(data, features = top10_genus[i], order = T)
        ggsave(str_glue("{outdir}/top10_genus_split/Featureplot_top10_genus_{name}_{i}.pdf"), height = 6, width = 7)
        ggsave(str_glue("{outdir}/top10_genus_split/Featureplot_top10_genus_{name}_{i}.png"), height = 6, width = 7)
    }
}

if(splitgroup == F){
    function_top10_genus_plot(data_seurat, "Main")
}else{
    out = lapply(unique(data_seurat$group), function(x){
        data_sub <- subset(data_seurat, subset = group == x)
        function_top10_genus_plot(data_sub, x)
    })
}

file.copy("/SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/doc/basic_plot/README.txt", str_glue("{outdir}/README.txt"))
cat("basic plot done. \n")