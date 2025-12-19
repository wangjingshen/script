suppressWarnings({
library(argparser)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(stringr)
library(DOSE)
library(ggplot2)
library(tibble)
library(enrichplot)
})
#source("/Public/Script/shouhou/SCRIPT/lib/scRNA/seuratIntegrate/seuratSingleron/seuratSingleron.R")

argv <- arg_parser('')
argv <- add_argument(argv, "--outdir", help="out dir", default='./')
argv <- add_argument(argv, "--prefix", help="prefix", default="")
argv <- add_argument(argv, "--difflist", help="diff list, two column")
argv <- add_argument(argv, "--species", help='the latin name of the species reference genome, usable: homo_sapiens, mus_musculus', default='homo_sapiens')
argv <- add_argument(argv, "--KEGG", help='run KEGG pathways enrich or not', default=T, type='logit')
argv <- add_argument(argv, "--GO", help='run GO pathways enrich or not', default=T, type='logit')
argv <- add_argument(argv, "--updown", help='the up- or down-regulated genes used to enrich, usable: up, down, updown, allsig, list', default='up')
argv <- add_argument(argv, "--FC_cutoff", 'the cutoff of fold_change to filter difference genes, must be positive', default=0.25,type='numeric')
argv <- add_argument(argv, '--Q_cutoff', help='the cutoff of pvals_adj to filter difference genes', default=1, type='numeric')
argv <- add_argument(argv, '--P_cutoff', help='the cutoff of p_value to filter difference genes', default=0.05, type='numeric')
argv <- add_argument(argv, '--selectpathway', help="the select pathway, split by ',', 非ID, Description")
argv <- add_argument(argv, '--enrichlist', help="enrich result list, tow column")
argv <- add_argument(argv, '--includeAll', help="logical, whether to fill the pathway", default = T, type='logit')
argv <- add_argument(argv, '--top', help="top pathway", default = 5, type='numeric')
argv <- add_argument(argv, '--ontology', help="which ontology to display", default='all')
argv <- parse_args(argv)

############################################################

############################################################


genome <- read.table("/SGRNJ03/pipline_test/chenkai/compareCluster/species_database_used.tsv", header = T, sep = '\t')
outdir <- argv$outdir
dir.create(outdir)

species <- argv$species
orgdb <- genome$orgdb_package[genome$species_type==species]
orgkg <- genome$organisms[genome$species_type==species]

diff_list <- argv$difflist
enrich_list <- argv$enrichlist

updown <- argv$updown

FC_cutoff <- as.numeric(argv$FC_cutoff)
Q_cutoff  <- as.numeric(argv$Q_cutoff)
P_cutoff  <- as.numeric(argv$P_cutoff)

prefix <- argv$prefix

selepath <- argv$selectpathway
if ( !is.na(selepath) ) { selepath <- unlist(strsplit(selepath, split = ',')) }

includeAll <- argv$includeAll

if (argv$ontology == 'all'){
    ontology <- c('BP','CC','MF')
} else {
    ontology <- unlist(strsplit(argv$ontology, split = ','))
}

top <- argv$top



# merge_data function
merge_data <- function(data_list, data_type='diff') {
    datalist <- read.table(data_list, header = F, sep = '\t')
    merged_result <- apply(datalist, 1, FUN=function(x) {
        data <- read.delim(x[length(x)], sep='\t', header=T)
        if (data_type == 'diff') {
            name2id <- bitr(data[, 1], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = orgdb, drop = F) %>% column_to_rownames(., var = 'SYMBOL')
            data <- cbind(data, name2id) %>% filter(., !is.na(ENTREZID))
        }
        data$group1 <- x[1]
        if (length(x) == 3) {data$group2 <- x[2]}
        return(data)
    }) %>% do.call(rbind, .)
    if (data_type == 'enrich') {merged_result$GeneRatio <- parse_ratio(merged_result$GeneRatio)}
    print(paste0('first_group: ', paste0(unique(merged_result$group1), collapse = ',')))
    if ('group2' %in% colnames(merged_result)) {print(paste0('second_group: ', paste0(unique(merged_result$group2), collapse = ',')))}
    return(merged_result)
}

# filter diffgene function
filter_diff <- function(data, updown) {
    type <- ifelse('names' %in% colnames(data), 'scanpy', 'seurat')
    switch(type,
        'scanpy' = {fc_col <- 'logfoldchanges'; pval_col <- 'pvals'; padj_col <- 'pvals_adj'},
        'seurat' = {fc_col <- 'avg_log2FC'; pval_col <- 'p_val'; padj_col <- 'p_val_adj'}
    )
    if (updown == 'up') {
        filtered_data <- filter(data, data[,fc_col] > FC_cutoff & data[,pval_col] <= P_cutoff & data[,padj_col] <= Q_cutoff )
    } else if (updown == 'down') {
        filtered_data <- filter(data, data[,fc_col] < -FC_cutoff & data[,pval_col] <= P_cutoff & data[,padj_col] <= Q_cutoff )
    } else if (updown == 'list') {
        filtered_data <- data
    }
    return(filtered_data)
}

# GO enrich function
GO_enrich <- function(data) {
    print("-------------------GO enrich---------------")
    if (! 'group2' %in% colnames(data)) {
        GO <- compareCluster(ENTREZID~group1, data=data, fun='enrichGO', OrgDb=orgdb, pvalueCutoff=1, pAdjustMethod="BH", qvalueCutoff=1, ont='ALL')
    } else {
        GO <- compareCluster(ENTREZID~group1+group2, data=data, fun='enrichGO', OrgDb=orgdb, pvalueCutoff=1, pAdjustMethod="BH", qvalueCutoff=1, ont='ALL')
    }
    GO <- pairwise_termsim(GO)				# 用于计算通路的相似性，用于simplify()去除冗余通路
    GO <- setReadable(GO, OrgDb = orgdb, keyType = 'ENTREZID')
    GO_result <- GO@compareClusterResult
    GO_result$GeneRatio <- parse_ratio(GO_result$GeneRatio)
    print('----------------------GO done-----------------')
    return(GO_result)
}

# KEGG enrich function
KEGG_enrich <- function(data) {
    print("--------------------KEGG enrich--------------------")
    if (! 'group2' %in% colnames(data)) {
        KEGG <- compareCluster(ENTREZID~group1, data=data, fun = 'enrichKEGG', organism = orgkg, pvalueCutoff=1, pAdjustMethod = "BH", qvalueCutoff = 1, use_internal_data=T)
    } else {
        KEGG <- compareCluster(ENTREZID~group1+group2, data=data, fun = 'enrichKEGG', organism = orgkg, pvalueCutoff=1, pAdjustMethod = "BH", qvalueCutoff = 1, use_internal_data=T)
    }
    KEGG <- setReadable(KEGG, OrgDb = orgdb, keyType = 'ENTREZID')
    KEGG_result <- KEGG@compareClusterResult
    KEGG_result$GeneRatio <- parse_ratio(KEGG_result$GeneRatio)
    print('-------------------KEGG done----------------------')
    return(KEGG_result)
}

# filter pathway function
filter_pathway <- function(enrich_result, selepath, top=5, includeAll=TRUE, updown=NULL, enrichtype) {
    # 筛选指定通路 或 按group1分组并按pvalue升序排列，每组取top
    if(all(!is.na(selepath))) {
        print('filter selected path')
        filtered_enrich <- filter(enrich_result, Description %in% selepath)
        if ( length(selepath[!selepath %in% enrich_result$Description]) ) {
            cat('The following pathways were not enriched:', paste0(selepath[!selepath %in% enrich_result$Description], collapse = ','), '\n')
        }
        if (length(intersect(selepath, enrich_result$Description)) == 0) {return(message('Not enriched for selected pathways, skip ', enrichtype))}
    } else {
        print(paste0('filter top', top, ' path'))
        if (enrichtype == 'GO') {enrich_result <- filter(enrich_result, ONTOLOGY %in% ontology)}
        filtered_enrich <- group_by(enrich_result, group1) %>% top_n(n=-top, wt=pvalue)
        if (includeAll) {
            selepath <- unique(filtered_enrich$Description)		
            filtered_enrich <- filter(enrich_result, Description %in% selepath)
        }
    }
    write.table(filtered_enrich, paste0(outdir, '/', prefix, '_', updown, '.', enrichtype, '_enrichment.xls'), col.names=T, row.names=F, sep='\t', quote=F)
    return(filtered_enrich)
}

# plot function
plot_fun <- function(data) {
    print('plot filtered enrich result')
    ONTOLOGY <- ifelse('ONTOLOGY' %in% colnames(data), TRUE, FALSE)
    group2 <- ifelse('group2' %in% colnames(data), TRUE, FALSE)

    # 无ONTOLOGY, 无第二分组
    p <- ggplot(data, aes(x = factor(group1, levels=unique(group1)), y = factor(Description, levels=rev(unique(Description))), color=p.adjust, size=GeneRatio)) + 
    geom_point()+
    scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + xlab(NULL) + DOSE::theme_dose(12) +
    scale_size_continuous(range=c(3, 8)) +
    theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color = "black")) +
    guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2))

    # 无ONTOLOGY, 有第二分组
    if (!ONTOLOGY & group2 ) {
        p <- p + facet_grid( ~ group2, scales = "free", space = "free")
    }

    # 有ONTOLOGY, 无第二分组 
    if (ONTOLOGY & !group2 ) {
        p <- p + facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
    }

    # 有ONTOLOGY, 有第二分组
    if (ONTOLOGY & group2) {
        p <- p + facet_grid(ONTOLOGY ~ group2, scales = "free", space = "free")
    }

    return(p)
}

# save figure function
save_fig <- function(p, filename, width, height) {
    pathway <- unique(p$data$Description) %>% as.vector(.)
    longest_pathway <- max(nchar(pathway),na.rm=T)/10
    xaxis_len <- max(nchar(unique(p$data$group1)))/10
    if ('group2' %in% colnames(p$data)) {xaxis_len <- xaxis_len * 2}
    width=4 + xaxis_len + longest_pathway
    height=4 + 0.23*length(pathway)

    pdf(paste0(outdir, '/', filename, '.pdf'), width=width, height=height)
    print(p)
    dev.off()

    png(paste0(outdir, '/', filename, '.png'), width=width, height=height, res = 600, units = 'in')
    print(p)
    dev.off()
}





# 读取基因list进行富集分析出图
if (!is.na(diff_list)) {
    diff_result <- merge_data(diff_list, data_type = 'diff')
    filtered_diff_list <- list()
    if (updown %in% c('updown', 'up')) {
        filtered_diff_list[['up']] <- filter_diff(diff_result, updown = 'up')
    } 
    if (updown %in% c('updown', 'down')) {
        filtered_diff_list[['down']] <- filter_diff(diff_result, updown = 'down')
    } 
    if (updown == 'list') {
        filtered_diff_list[['list']] <- diff_result
    }

    for (i in names(filtered_diff_list)) {
        if (argv$GO) {
            GO_result <- GO_enrich(filtered_diff_list[[i]])
            filtered_enrich <- filter_pathway(GO_result, selepath=selepath, top=top, includeAll=includeAll, updown = i, enrichtype = 'GO')
            if (nrow(filtered_enrich)!=0) {
                p <- plot_fun(filtered_enrich)
                save_fig(p, filename = paste0(prefix, '_', i, '.GO_dotplot'))
            }
        }
        if (argv$KEGG) {
            KEGG_result <- KEGG_enrich(filtered_diff_list[[i]])
            filtered_enrich <- filter_pathway(KEGG_result, selepath=selepath, top=top, includeAll=includeAll, updown = i, enrichtype = 'KEGG')
            if (!is.null(filtered_enrich)) {
                p <- plot_fun(filtered_enrich)
                save_fig(p, filename = paste0(prefix, '_', i, '.KEGG_dotplot'))
            }
        }
    }
}



# 读取富集分析结果出图
if (!is.na(enrich_list)) {
    enrich_result <- merge_data(enrich_list, data_type = 'enrich')
    enrichtype <- ifelse(grepl('GO',enrich_result$ID[1]), 'GO', 'KEGG')
    filtered_enrich <- filter_pathway(enrich_result, selepath=selepath, top=top, includeAll=includeAll, enrichtype=enrichtype, updown='enrich')
    p <- plot_fun(filtered_enrich)
    save_fig(p, filename=paste0(prefix, '.', enrichtype, '_dotplot'))
}