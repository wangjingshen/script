# transcriptome_diversity

suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
    library(vegan)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--split_group", help = "default:F")
argv <- add_argument(argv, "--analysis_genus", help = "analysis genus, default: total_genus")
argv <- add_argument(argv, "--nFeatures", help = "nFeatures, default:500")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 04.transcriptome_diversity")
argv <- parse_args(argv)

rds <- argv$rds
split_group <- ifelse(is.na(argv$split_group), "F", argv$split_group)
analysis_genus <- ifelse(is.na(argv$analysis_genus), "total_genus", argv$analysis_genus)
nFeatures <- ifelse(is.na(argv$nFeatures), 500, as.numeric(argv$nFeatures))
umi_threshold <- ifelse(analysis_genus == "total_genus", 2, 0)  # set umi threshold
outdir <- ifelse(is.na(argv$outdir), "04.transcriptome_diversity", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# readRDS
data_seurat <- readRDS(rds)
data_seurat$genus_detect <- factor(ifelse(data_seurat@meta.data[,analysis_genus] > umi_threshold, "Genus_detected", "No_genus_detected"), 
                                   levels = c("No_genus_detected", "Genus_detected"))

## shannon
hvg <- VariableFeatures(FindVariableFeatures(data_seurat, nfeatures = nFeatures))
norm <- t(as.matrix(data_seurat[["RNA"]]@data[hvg, ]))
#print(norm[1:3,1:3])
data_seurat$Shannon_diversity <- diversity(norm, index = "shannon")

df_shannon <- data.frame(Shannon_diversity = data_seurat$Shannon_diversity,
                         cluster = data_seurat$cluster,
                         group = data_seurat$group,
                         genus_detect = data_seurat$genus_detect)

# wilcox
if(split_group == "T"){
    df_shannon_stat <- ggpubr::compare_means(Shannon_diversity~genus_detect, df_shannon, method = "wilcox.test", group.by = c("cluster","group")) %>% dplyr::select(-.y.)
    df_shannon_stat$cluster2 = paste0(df_shannon_stat$cluster, "_", df_shannon_stat$group) # cluster_grop
}else{
    df_shannon_stat <- ggpubr::compare_means(Shannon_diversity~genus_detect, df_shannon, method = "wilcox.test", group.by = "cluster") %>% dplyr::select(-.y.)
    df_shannon_stat$cluster2 = df_shannon_stat$cluster   
}

function_get_mean <- function(df, name){
    if(split_group == "T"){
        df_mean <- aggregate(df$Shannon_diversity, by = list(df$cluster, df$group), FUN=mean)
        colnames(df_mean) <- c("cluster", "group", paste0("mean_", name))
        df_mean$cluster2 <- paste0(df_mean$cluster, "_", df_mean$group)
    }else{
        df_mean <- aggregate(df$Shannon_diversity, by = list(df$cluster), FUN=mean)
        colnames(df_mean) <- c("cluster", paste0("mean_", name))
        df_mean$cluster2 <- df_mean$cluster
    }
    return(df_mean[,c("cluster2", paste0("mean_", name))])
}

# cal mean
df_shannon_mean <- merge(function_get_mean(df_shannon %>% filter(genus_detect == "Genus_detected"), "Genus_detected"),
                         function_get_mean(df_shannon %>% filter(genus_detect == "No_genus_detected"), "No_genus_detected"), by = "cluster2")
df_shannon_stat <- merge(df_shannon_stat, df_shannon_mean, by = "cluster2")

# -log10p * 1(detect); -1(no_detect)
df_shannon_stat$p_plot = -log10(df_shannon_stat$p) * sign(df_shannon_stat$mean_Genus_detected - df_shannon_stat$mean_No_genus_detected)

# set color
df_color <- data.frame(cluster = levels(data_seurat$cluster),
                       color = color_protocol[1: length(levels(data_seurat$cluster))])

function_gplot <- function(df, g){
    if(split_group == "T"){
        df_gplot <- df %>%
            filter(group == g) %>%
            arrange(p_plot)
    }else{
        df_gplot <- df %>%
            arrange(p_plot)        
    }
    df_gplot$cluster <- factor(df_gplot$cluster, levels = df_gplot$cluster)
    p <- df_gplot %>%
        ggplot(aes(cluster, p_plot,fill=cluster)) + 
            geom_col(color="black", width=0.8)+
            scale_fill_manual(values = df_color$color[match(df_gplot$cluster, df_color$cluster)]) +  # 
            theme_classic() + 
            labs(x="Cluster",y="-log10(P)") + 
            ggtitle(paste0("group: ", g)) +
            theme(legend.position = "none") +
            coord_flip()
            #geom_hline(yintercept = -log10(0.05), color="red", linetype = 5)+  # dashed
    return(p)
}

if(split_group == "T"){
    options(repr.plot.height = 5, repr.plot.width = 12)
    do.call(`+`, lapply(unique(data_seurat$group), function_gplot, df = df_shannon_stat))
    ggsave(str_glue("{outdir}/shannon_diversity_barplot.png"), height = 5, width = 12)
    ggsave(str_glue("{outdir}/shannon_diversity_barplot.pdf"), height = 5, width = 12)
}else{
    options(repr.plot.height = 5, repr.plot.width = 7)
    function_gplot(df = df_shannon_stat, g = "all")
    ggsave(str_glue("{outdir}/shannon_diversity_barplot.png"), height = 5, width = 7)
    ggsave(str_glue("{outdir}/shannon_diversity_barplot.pdf"), height = 5, width = 7)
}

file.copy("/SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/doc/transcriptome_diversity/README.txt", str_glue("{outdir}/README.txt"))

cat("transcriptome diversity done. \n")