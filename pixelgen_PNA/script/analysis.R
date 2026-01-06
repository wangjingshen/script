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
    library(fs)
}))
options(future.globals.maxSize = 8 * 1024 * 1024 * 1024)   # 8 GiB

isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pxl", help = "pxl file, split by ,")
argv <- add_argument(argv, "--sample", help = "sample, split by ,")
argv <- add_argument(argv, "--group", help = "group, split by ,")
argv <- add_argument(argv, "--nUMI_cutoff", default = 10000, type = "integer", help = "nUMI_cutoff, default: 10000")
argv <- add_argument(argv, "--isotype_fraction_cutoff", default = 0.001, type = "double", help = "isotype_fraction_cutoff, default: 0.001")
argv <- add_argument(argv, "--normalize_method", default = "dsb", help="normalize_method, default: dsb")
argv <- add_argument(argv, "--rm_batch", default = "F", help="rm batch, T or F, default: F")
argv <- add_argument(argv, "--rm_batch_var", default ="sample", help="rm batch var, default: sample")
argv <- add_argument(argv, "--ndims", default = 10, type = "integer", help="ndims, default: 10")
argv <- add_argument(argv, "--resolution", default = 0.8, type = "double", help="resolution, default: 0.8")
argv <- add_argument(argv, "--spatial_vis", default = "F", help = "spatial_vis, default: F")
argv <- add_argument(argv, "--select_seurat_cluster", default = "0", help="select seurat cluster, default: 0")
argv <- add_argument(argv, "--select_component", default = "top1", help="select component, default: top1")
argv <- add_argument(argv, "--select_markers", default = "top1", help="select markers, default: top1")
argv <- add_argument(argv, "--outdir", default = "outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

pxl <- unlist(strsplit(argv$pxl, split = ","))
sample <- unlist(strsplit(argv$sample, split = ","))
group <- unlist(strsplit(argv$group, split = ","))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

# Validate input lengths
if (!(length(pxl) == length(sample) && length(sample) == length(group))) {
    stop("Error: --pxl, --sample, and --group must have the same number of elements.")
}

if(!dir_exists(outdir)){
    dir_create(outdir, "qc")
    dir_create(outdir, "cluster")
    #dir.create(str_glue("{outdir}/qc/"), recursive = T)
    #dir.create(str_glue("{outdir}/cluster/"), recursive = T)
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
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = argv$isotype_fraction_cutoff, linetype = "dashed")
ggsave(str_glue("{outdir}/qc/isotype_fraction.png"), plot = p, height = 5, width = 6 + 2*length(sample))

# filter --
pg_data_combined <- pg_data_combined %>% 
    subset(n_umi >= argv$nUMI_cutoff & tau_type == "normal" & isotype_fraction < argv$isotype_fraction_cutoff)

# Normalize --
if(argv$normalize_method == "dsb"){
    pg_data_combined <- pg_data_combined %>% 
        JoinLayers() %>% 
        Normalize(method = "dsb", isotype_controls = isotype_controls, assay = "PNA")
}
if(argv$normalize_method == "clr"){
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

if(argv$rm_batch == "T"){
    pg_data_combined <- RunHarmony(pg_data_combined, group.by.vars = argv$rm_batch_var, reduction = "pca", reduction.save = "harmony")
}

pg_data_combined <- pg_data_combined %>%
    RunUMAP(dims = 1 : argv$ndims, reduction = if(argv$rm_batch == "T") "harmony" else "pca", verbose = F) %>%
    FindNeighbors(dims = 1 : argv$ndims, reduction = if(argv$rm_batch == "T") "harmony" else "pca", verbose = F) %>%
    FindClusters(random.seed = 1, resolution = argv$resolution, verbose = F)

# umap plot
function_save_plot <- function(plot_func, file_name, data, h = 5, w = NULL) {
    if (is.null(w)) w <- 7 + 5 * length(unique(data$sample))
    p <- plot_func(data)
    ggsave(file.path(outdir, file_name), plot = p, height = h, width = w)
    invisible(TRUE)
}

dim_plots <- list(
    sample = function(x) DimPlot(x, reduction = "umap", group.by = "sample") + coord_fixed(),
    group = function(x) DimPlot(x, reduction = "umap", group.by = "group") + coord_fixed(),
    seurat_clusters = function(x) DimPlot(x, reduction = "umap", group.by = "seurat_clusters") + coord_fixed()
)
future_walk2(names(dim_plots), dim_plots,
             function(name, plot_obj) function_save_plot(plot_obj, str_glue("cluster/{name}.png"), data = pg_data_combined, h = 5, w = 8),
             .options = furrr_options(seed = NULL))

saveRDS(pg_data_combined, str_glue("{outdir}/data.rds"))

# deg
deg_clusters <- FindAllMarkers(pg_data_combined, test.use = "wilcox", assay = "PNA", slot = "data", 
                               logfc.threshold = 0, return.thresh = 1, mean.fxn = rowMeans, fc.name = "difference", verbose = F)
#write.table(deg_clusters[, c(7,1:6)], str_glue("{outdir}/cluster/deg_clusters.tsv"), sep="\t", row.names=F, quote =F)
write_tsv(deg_clusters, file.path(outdir, "cluster", "deg_clusters.tsv"))

if(argv$spatial_vis %in% c("T","True","TRUE")){
    if(!dir_exists(str_glue("{outdir}/spatial"))){
        dir_create(outdir, "spatial")
    }
    #if(!dir.exists(str_glue("{outdir}/spatial/"))){
    #    dir.create(str_glue("{outdir}/spatial/"), recursive = T)
    #}
    # Cell Visualization
    clustering_scores <- ProximityScores(pg_data_combined, meta_data_columns = c("seurat_clusters"), add_marker_counts = TRUE) %>% 
        filter(seurat_clusters == argv$select_seurat_cluster) %>% 
        filter(marker_1 == marker_2)

    select_component <- ifelse(argv$select_component == "top1",
                               clustering_scores$component[1],
                               unlist(strsplit(argv$select_component, split = ",")))

    select_markers <- ifelse(argv$select_markers == "top1", 
                             deg_clusters %>% filter(cluster == argv$select_seurat_cluster) %>% pull(gene) %>% head(1),
                             unlist(strsplit(argv$select_markers, split = ",")))
 
    pg_data_combined <- ComputeLayout(pg_data_combined, layout_method = "wpmds")
    pg_data_combined <- LoadCellGraphs(pg_data_combined, cells = select_component, add_layouts = TRUE, verbose = FALSE)

    with_progress({
        future_walk(select_markers,
                    function(m) {
                        p <- Plot2DGraph(pg_data_combined, cells = select_component, marker = m)
                        ggsave(file.path(outdir, "spatial", str_glue("{m}_2D.png")), p, height = 5, width = 6)
                    }, .options = furrr_options(seed = NULL))
    })

    with_progress({
        future_walk(select_markers,
                    function(m) {
                        cg <- CellGraphs(pg_data_combined)[[select_component]]
                        xyz <- cg@layout$wpmds_3d %>% mutate(node_val = cg@counts[, m])
                        render_rotating_layout(data = xyz, str_glue("{outdir}/spatial/{m}_3D.gif"), 
                                               max_degree = 180, frames = 200, show_first_frame = FALSE, 
                                               colors = c("lightgrey", "red"), width = 1e3, height = 1e3)
                    }, .options = furrr_options(seed = NULL))
    })
}

logger::log_success("Pipeline completed successfully.")