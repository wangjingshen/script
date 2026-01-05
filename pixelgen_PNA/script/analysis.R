suppressMessages(suppressWarnings({
    library(pixelatorR)
    library(Seurat)
    library(SeuratObject)
    library(tidyverse)
    library(harmony)
    library(here)
    library(tibble)
    library(argparser)
    library(logger)
    library(future)
    library(furrr)
}))
options(future.globals.maxSize = 8 * 1024 * 1024 * 1024)   # 8 GiB

isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pxl", help = "pxl file, split by ,")
argv <- add_argument(argv, "--sample", help = "sample, split by ,")
argv <- add_argument(argv, "--group", help = "group, split by ,")
argv <- add_argument(argv, "--nUMI_cutoff", help = "nUMI_cutoff, default: 10000")
argv <- add_argument(argv, "--isotype_fraction_cutoff", help = "isotype_fraction_cutoff, default: 0.001")
argv <- add_argument(argv, "--normalize_method", help="normalize_method, default: dsb")
argv <- add_argument(argv, "--rm_batch", help="rm batch, T or F, default: F")
argv <- add_argument(argv, "--rm_batch_var", help="rm batch var, default: sample")
argv <- add_argument(argv, "--ndims", help="ndims, default: 10")
argv <- add_argument(argv, "--resolution", help="resolution, default: 0.8")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

pxl <- unlist(strsplit(argv$pxl, split = ","))
sample <- unlist(strsplit(argv$sample, split = ","))
group <- unlist(strsplit(argv$group, split = ","))
nUMI_cutoff <- ifelse(is.na(argv$nUMI_cutoff), 10000, as.numeric(argv$nUMI_cutoff))
isotype_fraction_cutoff <- ifelse(is.na(argv$isotype_fraction_cutoff), 0.001, as.numeric(argv$isotype_fraction_cutoff))
normalize_method <- ifelse(is.na(argv$normalize_method), "dsb", argv$normalize_method)
rm_batch <- ifelse(is.na(argv$rm_batch), "F", argv$rm_batch)
rm_batch_var <- ifelse(is.na(argv$rm_batch_var), "sample", argv$rm_batch_var)
ndims <- ifelse(is.na(argv$ndims), 10, as.numeric(argv$ndims))
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

# Validate input lengths
if (!(length(pxl) == length(sample) && length(sample) == length(group))) {
    stop("Error: --pxl, --sample, and --group must have the same number of elements.")
}

if(!dir.exists(str_glue("{outdir}/qc/"))){
    dir.create(str_glue("{outdir}/qc/"), recursive = T)
}
if(!dir.exists(str_glue("{outdir}/cluster/"))){
    dir.create(str_glue("{outdir}/cluster/"), recursive = T)
}

# read data --
log_info("Reading PXL files...")
#message("Reading PXL files...")

function_read_sample <- function(i) {
    data <- ReadPNA_Seurat(pxl[i])
    data$sample <- sample[i]
    data$group  <- group[i]
    return(data)
}
obj_list <- future_map(1:length(pxl), function_read_sample, .options = furrr_options(seed = NULL))
names(obj_list) <- sample

pg_data_combined <- if (length(obj_list) == 1) {
    obj_list[[1]]
} else {
    merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
} %>% JoinLayers()

# qc --
p <- MoleculeRankPlot(pg_data_combined, group_by = "sample")
ggsave(str_glue("{outdir}/qc/molecule_rank_plot.png"), plot = p, height = 5, width = 7)

p <- TauPlot(pg_data_combined, group_by = "sample")
ggsave(str_glue("{outdir}/qc/TauPlot.png"), plot =p, height = 5, width = 5 + 5*length(sample))

p <- pg_data_combined[[]] %>% 
    select(sample, isotype_fraction) %>% 
    ggplot(aes(sample, isotype_fraction, fill = sample)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_jitter(size = 0.4, alpha = 0.5) +
    theme_bw() +
    labs(y = "Fraction of isotype control counts") +
    scale_y_continuous(labels = scales::percent)
ggsave(str_glue("{outdir}/qc/isotype_fraction.png"), plot = p, height = 5, width = 6 + 2*length(sample))

# filter --
pg_data_combined <- pg_data_combined %>% 
    subset(n_umi >= nUMI_cutoff & tau_type == "normal" & isotype_fraction < isotype_fraction_cutoff)

# Normalize --
if(normalize_method == "dsb"){
    pg_data_combined <- pg_data_combined %>% 
        JoinLayers() %>% 
        Normalize(method = "dsb", isotype_controls = isotype_controls, assay = "PNA")
}
if(normalize_method == "clr"){
    pg_data_combined[["clr_assay"]] <- pg_data_combined[["PNA"]]
    pg_data_combined <- pg_data_combined %>% 
        JoinLayers() %>% 
        Normalize(method = "clr", assay = "clr_assay")  
}

# cluster --
DefaultAssay(pg_data_combined) <- "PNA"
pg_data_combined <- pg_data_combined %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F)

if(rm_batch == "T"){
    pg_data_combined <- RunHarmony(pg_data_combined, group.by.vars = rm_batch_var, reduction = "pca", reduction.save = "harmony")
}

pg_data_combined <- pg_data_combined %>%
    RunUMAP(dims = 1 : ndims, reduction = if(rm_batch == "T") "harmony" else "pca", verbose = F) %>%
    FindNeighbors(dims = 1 : ndims, reduction = if(rm_batch == "T") "harmony" else "pca", verbose = F) %>%
    FindClusters(random.seed = 1, resolution = resolution, verbose = F)

# umap plot
function_save_plot <- function(plot_func, file_name, data, h = 5, w = NULL) {
    if (is.null(w)) w <- 7 + 5 * length(unique(data$sample))
    p <- plot_func(data)
    ggsave(file.path(outdir, file_name), plot = p, height = h, width = w)
    invisible(TRUE)
}

dim_plots <- list(
    sample = function(x) DimPlot(x, reduction = "umap", group.by = "sample") + ggtitle("sample") + coord_fixed() + theme(plot.title = element_text(hjust = 0.5)),
    group = function(x) DimPlot(x, reduction = "umap", group.by = "group") + ggtitle("group") + coord_fixed() + theme(plot.title = element_text(hjust = 0.5)),
    seurat_clusters = function(x) DimPlot(x, reduction = "umap", group.by = "seurat_clusters") + coord_fixed()
)
future_walk2(names(dim_plots), dim_plots,
             function(name, plot_obj) function_save_plot(plot_obj, str_glue("cluster/{name}.png"), data = pg_data_combined, h = 5, w = 8),
             .options = furrr_options(seed = NULL))

saveRDS(pg_data_combined, str_glue("{outdir}/data.rds"))
logger::log_success("Pipeline completed successfully.")