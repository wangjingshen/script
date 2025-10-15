suppressMessages(suppressWarnings({
    library(pixelatorR)
    library(SeuratObject)
    library(dplyr)
    library(stringr)
    library(here)
    library(tibble)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pxl", help = "pxl file, split by ,")
argv <- add_argument(argv, "--spname", help = "spname, split by ,")
argv <- add_argument(argv, "--nUMI_cutoff", help = "nUMI_cutoff, default: 10000")
argv <- add_argument(argv, "--isotype_fraction_cutoff", help = "isotype_fraction_cutoff, default: 0.001")
argv <- add_argument(argv, "--rm_batch", help="rm batch, yes or no, default: no")
argv <- add_argument(argv, "--ndims", help="ndims, default: 10")
argv <- add_argument(argv, "--resolution", help="resolution, default: 0.8")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

pxl <- unlist(strsplit(argv$pxl, split = ","))
spname <- unlist(strsplit(argv$spname, split = ","))
nUMI_cutoff <- ifelse(is.na(argv$nUMI_cutoff), 10000, as.numeric(argv$nUMI_cutoff))
isotype_fraction_cutoff <- ifelse(is.na(argv$isotype_fraction_cutoff), 0.001, as.numeric(argv$isotype_fraction_cutoff))
rm_batch <- ifelse(is.na(argv$rm_batch), "no", argv$rm_batch)
ndims <- ifelse(is.na(argv$ndims), 10, as.numeric(argv$ndims))
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# read data -- 
obj_list <- lapply(matrix_10X, ReadPNA_Seurat)
names(obj_list) <- spname

if(length(spname)==1){
    pg_data_combined <- obj_list[[1]]
    pg_data_combined <- RenameCells(pg_data_combined, add.cell.id = names(obj_list))
}else{
    pg_data_combined <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
}
pg_data_combined <- JoinLayers(pg_data_combined)
print(pg_data_combined)

# molecule rank plot
options(repr.plot.height = 5, repr.plot.width = 7)
MoleculeRankPlot(pg_data_combined, group_by = "condition")
ggsave(str_glue("{outdir}/molecule_rank_plot.png"), height = 5, width = 7)

# TauPlot
options(repr.plot.height = 5, repr.plot.width = 5 + 5*length(spname))
TauPlot(pg_data_combined, group_by = "condition")
ggsave(str_glue("{outdir}/TauPlot.png"), height = 5, width = 5 + 5*length(spname))

# isotype fraction
options(repr.plot.height = 5, repr.plot.width = 7)
pg_data_combined[[]] %>% 
    select(condition, isotype_fraction) %>% 
    ggplot(aes(condition, isotype_fraction, fill = condition)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_jitter(size = 0.4, alpha = 0.5) +
    theme_bw() +
    labs(y = "Fraction of isotype control counts") +
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = isotype_control_fraction_cutoff, linetype = "dashed")
ggsave(str_glue("{outdir}/isotype_fraction.png"), height = 5, width = 7)

# filter --
pg_data_combined <- pg_data_combined %>% 
    subset(n_umi >= nUMI_cutoff & tau_type == "normal" & isotype_fraction < isotype_fraction_cutoff)


# Normalize
isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")
# dsb
pg_data_combined <- pg_data_combined %>% 
    JoinLayers() %>% 
    Normalize(method = "dsb", isotype_controls = isotype_controls, assay = "PNA")


DefaultAssay(pg_data_combined) <- "PNA"
pg_data_combined <- pg_data_combined %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F)

if(rm_batch == "no"){
    pg_data_combined <- pg_data_combined %>%
        RunUMAP(dims = 1 : ndims, verbose = F) %>%
        FindNeighbors(dims = 1 : ndims, verbose = F) %>%
        FindClusters(random.seed = 1, resolution = resolution, verbose = F)
}else{
    pg_data_combined <- pg_data_combined %>%
        RunHarmony(group.by.vars = "condition", reduction = "pca", reduction.save = "harmony") %>%
        FindNeighbors(dims = 1 : 10, reduction = "harmony", verbose = F) %>%
        FindClusters(random.seed = 1, resolution = resolution, verbose = F)
}

