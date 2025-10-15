suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

dirnameScript <- function(){
    # get full directory path of current script located
    cmd = commandArgs(trailingOnly = FALSE)
    scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
    if (length(scriptName) > 0) {
        path = normalizePath(scriptName)
        dirname = dirname(path)
    } else {
        print('Warning: not a runtime environment, using current directory instead.')
        dirname = getwd()
    }
    return(dirname)
}
 
sourceScript <- function(x, dir='/SGRNJ06/randd/USER/wangjingshen/script_dev/') {
    # try source script from directory one by one: file itself, dir from argument, script directory
    scriptdir = dirnameScript()
    searchpath = c(x,
                   file.path(dir, x),
                   file.path(scriptdir, x))
    sourced = FALSE
    for (i in searchpath) {
        if (file.exists(i)) {source(i); print(i); sourced = TRUE; break}
    }
    if (!sourced) {stop(paste0('can not source: ', x))}
}
scriptdir = dirnameScript() #

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--data", help = "matrix or fragments")
argv <- add_argument(argv, "--percent_threshold", help = "percent threshold, Default: 0.8")
argv <- add_argument(argv, "--mode", help = "mode, cr or atac or rna")
argv <- add_argument(argv, "--name", help = "name")
argv <- parse_args(argv)

data <- argv$data
percent_threshold <- ifelse(is.na(argv$percent_threshold), 0.8, as.numeric(argv$percent_threshold))
mode <- argv$mode
name <- argv$name

if(!dir.exists(name)){
    dir.create(name, recursive = TRUE)
}

function_plot <- function(data_count, mode, name){
    if(mode == "cr_atac" | mode == "cr_rna" | mode == "atac"){
        hs_genes <- row.names(data_count)[grep("GRCh38", row.names(data_count))]
        mm_genes <- row.names(data_count)[grep("mm10", row.names(data_count))]
        hs_percent <- colSums(data_count[ hs_genes,])/ colSums(data_count)
        mm_percent <- colSums(data_count[ mm_genes,])/ colSums(data_count)
        if(mode == "cr_rna"){
            lab_name = "transcripts"
        }else{
            lab_name = "fragment"
        }
    }
    if(mode == "rna"){
        hs_genes <- intersect(row.names(data_count),read.table(str_glue("{scriptdir}/../data/hs_genes.tsv"))[,1])
        mm_genes <- intersect(row.names(data_count),read.table(str_glue("{scriptdir}/../data/mmu_genes.tsv"))[,1])
        hs_percent <- colSums(data_count[ hs_genes,])/ colSums(data_count)
        mm_percent <- colSums(data_count[ mm_genes,])/ colSums(data_count)
        lab_name = "transcripts"       
    }

    hs_bc <- colnames(data_count)[ hs_percent > percent_threshold & mm_percent < 1 - percent_threshold]
    mm_bc <- colnames(data_count)[ mm_percent > percent_threshold & hs_percent < 1 - percent_threshold]
    
    bc_df <- data.frame(barcode = colnames(data_count),
                        hs_umi = colSums(data_count[hs_genes,]),
                        mm_umi = colSums(data_count[mm_genes,]),
                        species = colnames(data_count))

    bc_df$species[ bc_df$barcode %in% hs_bc] ="human"
    bc_df$species[ bc_df$barcode %in% mm_bc] ="mouse"
    bc_df$species[ !bc_df$barcode %in% c(hs_bc, mm_bc)] ="doublets"
    print(table(bc_df$species))

    df_stat <- data.frame(
        hs_number = sum(bc_df$species == "human"),
        mm_number = sum(bc_df$species == "human"),
        doublets_number = sum(bc_df$species == "doublets"),
        median_umi_hs_in_human = median( bc_df$hs_umi[ bc_df$species == "human"]),
        median_umi_mm_in_human = median( bc_df$mm_umi[ bc_df$species == "human"]),
        median_umi_hs_in_mouse = median( bc_df$hs_umi[ bc_df$species == "mouse"]),
        median_umi_mm_in_mouse = median( bc_df$mm_umi[ bc_df$species == "mouse"]),

        mean_umi_hs_in_human = median( bc_df$hs_umi[ bc_df$species == "human"]),
        mean_umi_mm_in_human = median( bc_df$mm_umi[ bc_df$species == "human"]),
        mean_umi_hs_in_mouse = median( bc_df$hs_umi[ bc_df$species == "mouse"]),
        mean_umi_mm_in_mouse = median( bc_df$mm_umi[ bc_df$species == "mouse"])
    )
    write.table(t(df_stat), str_glue("{name}/{name}_stat.tsv"), col.names=F, sep="\t", quote=F)

    # 计算每个类别的点的数量
    category_counts <- bc_df %>%
        group_by(species) %>%
        summarise(count = n()) %>%
        ungroup()
    #print(category_counts)
    # plot
    options(repr.plot.height=5.5, repr.plot.width =7)
    #print(summary(bc_df$hs_umi))
    #print(summary(bc_df$mm_umi))
    ggplot(bc_df, aes(x=hs_umi,y=mm_umi,color=species)) +
        geom_point() +
        xlab(str_glue("Number of human {lab_name}")) +
        ylab(str_glue("Number of mouse {lab_name}")) +
        theme_classic() +
            scale_color_manual(
            values = c("human" = "red", "mouse" = "blue", "doublets" = "grey"),
            labels = paste(category_counts$species, " ( n=",category_counts$count, ")")) +
        theme(legend.title = element_text(size = 12),legend.text = element_text(size = 10))
    ggsave(str_glue("{name}/{lab_name}_doublets.pdf"), width = 7, height = 5.5)
    ggsave(str_glue("{name}/{lab_name}_doublets.png"), width = 7, height = 5.5)
}


if(mode =="cr"){
    data <- Read10X(data)
    function_plot(data$`Gene Expression`, str_glue("{mode}_rna"), name)
    function_plot(data$Peaks, str_glue("{mode}_atac"), name)
}

if(mode == "rna"){
    function_plot(Read10X(data), mode, name)
}

if(mode == "atac"){
    function_plot(Read10X_h5(data), mode, name)
}
