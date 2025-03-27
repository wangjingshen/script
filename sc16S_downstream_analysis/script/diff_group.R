# ref: Tumor microbiome links cellular programs and immunity in pancreatic cancer

suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--split_cluster", help = "default:F")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

rds <- argv$rds
split_cluster <- ifelse(is.na(argv$split_cluster), "F", argv$split_cluster)
outdir <- ifelse(is.na(argv$outdir), "02.diff_group", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

data_seurat <- readRDS(rds)

## 02 diff ----
genus <- colnames(data_seurat@meta.data)[(which(colnames(data_seurat@meta.data) == "total_genus") + 1) : ncol(data_seurat@meta.data)]

df_genus <- cbind.data.frame(data_seurat@meta.data[,genus], 
                             data_seurat@meta.data[,c("cluster","group")])


function_wilcox_group <- function(data, analysis_g, out){
    data$compare_g <- as.character(data$group)
    data$compare_g[ data$group != analysis_g] = "others"
    data$compare_g <- factor(data$compare_g, levels = c(analysis_g, "others"))  # set levels
    print(table(data$compare_g))

    genus_diff_stat <- do.call(rbind, lapply(genus, function(x){
        data$analysis_genus <- data@meta.data[,x]
        if(sum(data$analysis_genus)>0){  # filter all 0 genus
            stat_df <- ggpubr::compare_means(analysis_genus~compare_g, data@meta.data, method = "wilcox.test") %>% dplyr::select(-.y.)
            mean_df <- aggregate(data$analysis_genus, by = list(data$compare_g), FUN=mean)
            #print(mean_df)
            #      Group.1          x
            #1 ZCY2435298P 0.02742898
            #2      others 0.01848249
            for (i in 1:nrow(mean_df)){
                stat_df[,paste0("mean_", mean_df[i,1])] = mean_df[i,2]
            }
            stat_df <- as.data.frame(stat_df)
            #print(class(stat_df))
            colnames(stat_df)[ which(colnames(stat_df) == paste0("mean_", analysis_g))] = "mean1"
            colnames(stat_df)[ which(colnames(stat_df) == "mean_others")] = "mean2"
            stat_df$minus <- stat_df[, "mean1"] - stat_df[, "mean2"]  # minus
            stat_df$genus <- x
            return(stat_df[,c(ncol(stat_df), 1:(ncol(stat_df)-1))])
        }else{
            cat(paste0(x, ": all 0, skip.\n"))  # must use cat; if print, will merge to genus_diff_stat
        }
    }))
    genus_diff_stat <- arrange(genus_diff_stat, p.adj)
    write.table(genus_diff_stat, str_glue("{out}/{analysis_g}_wilcox_test_sigall.tsv"), sep = "\t", quote = F, row.names = F)
    return(genus_diff_stat)
    # p < 0.05
    #genus_diff_stat_sig <- genus_diff_stat[ genus_diff_stat$p <0.05,]
    #write.table(genus_diff_stat_sig, str_glue("{outdir}/{analysis_g}_wilcox_test_sig.tsv"), sep = "\t", quote = F, row.names = F)
}

if(split_cluster == "F"){
    res <- do.call(rbind, lapply(unique(data_seurat$group), function_wilcox_group, data = data_seurat, out = outdir))
    write.table(res, str_glue("{outdir}/wilcox_test_sigall.tsv"), sep = "\t", quote = F, row.names = F)

    file.copy("/SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/doc/diff_group/README.txt", str_glue("{outdir}/README.txt"))
}else{
    for (i in unique(data_seurat$cluster)){
        outpath = str_glue("{outdir}/cluster_{i}/")
        if(!dir.exists(outpath)){
            dir.create(outpath, recursive = T)
        }
        data_sub <- subset(data_seurat, subset = cluster == i)
        res <- do.call(rbind, lapply(unique(data_sub$group), function_wilcox_cluster, data = data_sub, out = outpath))
        write.table(res, str_glue("{outpath}/wilcox_test_sigall_cluster_{i}.tsv"), sep = "\t", quote = F, row.names = F)

        file.copy("/SGRNJ06/randd/USER/wangjingshen/script/sc16S_down/doc/diff_group/README.txt", str_glue("{outpath}/README.txt"))
    }
}

cat("diff done. \n")