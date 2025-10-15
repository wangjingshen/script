# diff cluster or group

suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
    library(circlize)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--mode", help = "mode, cluster or group, default: cluster")
argv <- add_argument(argv, "--split", help = "default:F")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 02.diff")
argv <- parse_args(argv)

rds <- argv$rds
mode <- ifelse(is.na(argv$mode), "cluster", argv$mode)
split <- ifelse(is.na(argv$split), "F", argv$split)
outdir <- ifelse(is.na(argv$outdir), "02.diff", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# readRDS
data_seurat <- readRDS(rds)
data_seurat$cluster <- as.character(data_seurat$cluster)
data_seurat$group <- as.character(data_seurat$group)
genus <- colnames(data_seurat@meta.data)[(which(colnames(data_seurat@meta.data) == "total_genus") + 1) : ncol(data_seurat@meta.data)]

function_wilcox <- function(data, analysis_thing, out){
    # analysis_thing
    # mode: cluster, clusterA
    # mode: group, groupA
    if(mode == "cluster"){
        data$compare_group <- as.character(data$cluster)
        if(length(unique(data$compare_group)) > 1){
            data$compare_group[ data$cluster != analysis_thing] = "others"
            data$compare_group <- factor(data$compare_group, levels = c(analysis_thing, "others"))  # set levels
        }else{
            print("just one compare_group, skip.")
            print(table(data$compare_group))
            return(NULL)
        }
    }
    if(mode == "group"){
        data$compare_group <- as.character(data$group)
        if(length(unique(data$compare_group)) > 1){
            data$compare_group[ data$group != analysis_thing] = "others"
            data$compare_group <- factor(data$compare_group, levels = c(analysis_thing, "others"))  # set levels
        }else{
            print("just one compare_group, skip.")
            print(table(data$compare_group))
            return(NULL)
        }
    }
    #print(table(data$compare_group))

    genus_diff_stat <- do.call(rbind, lapply(genus, function(x){
        data$analysis_genus <- data@meta.data[,x]
        if(sum(data$analysis_genus) > 0){
            stat_df <- ggpubr::compare_means(analysis_genus~compare_group, data@meta.data, method = "wilcox.test") %>% dplyr::select(-.y.)
            mean_df <- aggregate(data$analysis_genus, by = list(data$compare_group), FUN=mean)
            for (i in 1:nrow(mean_df)){
                stat_df[,paste0("mean_", mean_df[i,1])] = mean_df[i,2]
            }
            stat_df <- as.data.frame(stat_df)
            colnames(stat_df)[ which(colnames(stat_df) == paste0("mean_", analysis_thing))] = "mean1"
            colnames(stat_df)[ which(colnames(stat_df) == "mean_others")] = "mean2"
            stat_df$minus <- stat_df[, "mean1"] - stat_df[, "mean2"]  # minus
            stat_df$genus <- x
            return(stat_df[,c(ncol(stat_df), 1:(ncol(stat_df)-1))])
        }#else{
            #cat(paste0(x, ": all 0, skip.\n"))
        #}
    }))
    genus_diff_stat <- arrange(genus_diff_stat, p.adj)
    write.table(genus_diff_stat, str_glue("{out}/{analysis_thing}_wilcox_test_sigall.tsv"), sep = "\t", quote = F, row.names = F)
    return(genus_diff_stat)
}

if(mode == "group"){
    if(split == "F"){
        res <- do.call(rbind, lapply(unique(data_seurat$group), function_wilcox, data = data_seurat, out = outdir))
        write.table(res, str_glue("{outdir}/wilcox_test_sigall.tsv"), sep = "\t", quote = F, row.names = F)
    }else{
        for (i in unique(data_seurat$cluster)){
            outpath = str_glue("{outdir}/cluster_{i}/")
            if(!dir.exists(outpath)){
                dir.create(outpath, recursive = T)
            }
            print(paste0("analysis cluster: ", i))
            data_sub <- subset(data_seurat, subset = cluster == i)
            res <- do.call(rbind, lapply(unique(data_sub$group), function_wilcox, data = data_sub, out = outpath))
            write.table(res, str_glue("{outpath}/wilcox_test_sigall_cluster_{i}.tsv"), sep = "\t", quote = F, row.names = F)
        }
    }
}

if(mode == "cluster"){
    if(split == "F"){
        res <- do.call(rbind, lapply(unique(data_seurat$cluster), function_wilcox, data = data_seurat, out = outdir))
        write.table(res, str_glue("{outdir}/wilcox_test_sigall.tsv"), sep = "\t", quote = F, row.names = F)
    }else{
        for (i in unique(data_seurat$group)){
            outpath = str_glue("{outdir}/group_{i}/")
            if(!dir.exists(outpath)){
                dir.create(outpath, recursive = T)
            }
            print(paste0("analysis group: ", i))
            data_sub <- subset(data_seurat, subset = group == i)
            res <- do.call(rbind, lapply(unique(data_sub$cluster), function_wilcox, data = data_sub, out = outpath))
            write.table(res, str_glue("{outpath}/wilcox_test_sigall_group_{i}.tsv"), sep = "\t", quote = F, row.names = F)
        }
    }
}

file.copy("/SGRNJ06/randd/USER/wangjingshen/script_dev/sc16S_downstream_analysis/doc/diff.txt", str_glue("{outdir}/README.txt"))
cat("diff done. \n")