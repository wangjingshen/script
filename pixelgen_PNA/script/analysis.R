suppressMessages(suppressWarnings({
    library(pixelatorR)
    library(Seurat)
    library(SeuratObject)
    library(tidyverse)
    library(harmony)
    library(here)
    library(tibble)
    library(argparser)
}))

#
isotype_control_fraction_cutoff <- 0.001

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pxl", help = "pxl file, split by ,")
argv <- add_argument(argv, "--spname", help = "spname, split by ,")
argv <- add_argument(argv, "--gname", help = "gname, split by ,")
argv <- add_argument(argv, "--nUMI_cutoff", help = "nUMI_cutoff, default: 10000")
argv <- add_argument(argv, "--isotype_fraction_cutoff", help = "isotype_fraction_cutoff, default: 0.001")
argv <- add_argument(argv, "--normalize_method", help="normalize_method, default: dsb")
argv <- add_argument(argv, "--rm_batch", help="rm batch, yes or no, default: no")
argv <- add_argument(argv, "--rm_batch_var", help="rm batch var, default: spname")
argv <- add_argument(argv, "--ndims", help="ndims, default: 10")
argv <- add_argument(argv, "--resolution", help="resolution, default: 0.8")
argv <- add_argument(argv, "--species", help="species, default: human")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

pxl <- unlist(strsplit(argv$pxl, split = ","))
spname <- unlist(strsplit(argv$spname, split = ","))
gname <- unlist(strsplit(argv$gname, split = ","))
nUMI_cutoff <- ifelse(is.na(argv$nUMI_cutoff), 10000, as.numeric(argv$nUMI_cutoff))
isotype_fraction_cutoff <- ifelse(is.na(argv$isotype_fraction_cutoff), 0.001, as.numeric(argv$isotype_fraction_cutoff))
normalize_method <- ifelse(is.na(argv$normalize_method), "dsb", argv$normalize_method)
rm_batch <- ifelse(is.na(argv$rm_batch), "no", argv$rm_batch)
rm_batch_var <- ifelse(is.na(argv$rm_batch_var), "spname", argv$rm_batch_var)
ndims <- ifelse(is.na(argv$ndims), 10, as.numeric(argv$ndims))
resolution <- ifelse(is.na(argv$resolution), 0.8, as.numeric(argv$resolution))
species <- ifelse(is.na(argv$species), "human", as.numeric(argv$species))
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# read data -- 
obj_list <- lapply(1:length(pxl), function(x){
    data <- ReadPNA_Seurat(pxl[[x]])
    data$spname <- spname[[x]]
    data$gname <- gname[[x]]
    return(data)
})
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
p <- MoleculeRankPlot(pg_data_combined, group_by = "spname")
ggsave(str_glue("{outdir}/molecule_rank_plot.png"), plot = p, height = 5, width = 7)

# TauPlot
options(repr.plot.height = 5, repr.plot.width = 5 + 5*length(spname))
p <- TauPlot(pg_data_combined, group_by = "spname")
ggsave(str_glue("{outdir}/TauPlot.png"), height = 5, plot =p, width = 5 + 5*length(spname))

# isotype fraction
options(repr.plot.height = 5, repr.plot.width = 7)
p <- pg_data_combined[[]] %>% 
    select(spname, isotype_fraction) %>% 
    ggplot(aes(spname, isotype_fraction, fill = spname)) +
    geom_violin(draw_quantiles = 0.5) + 
    geom_jitter(size = 0.4, alpha = 0.5) +
    theme_bw() +
    labs(y = "Fraction of isotype control counts") +
    scale_y_continuous(labels = scales::percent)
    #geom_hline(yintercept = isotype_control_fraction_cutoff, linetype = "dashed")
ggsave(str_glue("{outdir}/isotype_fraction.png"), plot = p, height = 5, width = 7)

# filter --
pg_data_combined <- pg_data_combined %>% 
    subset(n_umi >= nUMI_cutoff & tau_type == "normal" & isotype_fraction < isotype_fraction_cutoff)


# Normalize
isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b")
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
        RunHarmony(group.by.vars = rm_batch_var, reduction = "pca", reduction.save = "harmony") %>%
        RunUMAP(reduction = "harmony", dims = 1:ndims, verbose = F) %>%
        FindNeighbors(dims = 1 : 10, reduction = "harmony", verbose = F) %>%
        FindClusters(random.seed = 1, resolution = resolution, verbose = F)
}

# singleR --
function_singleR <- function(data, species){
    if(species == "human"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/HumanPrimaryCellAtlasData.rds")
    }
    if(species == "mouse"){
        ref = readRDS("/SGRNJ06/randd/USER/wangjingshen/share/singleR_ref_rds/MouseRNAseqData.rds")
    }
    
    data_anno <- data@assays$PNA@layers$data
    row.names(data_anno) <- row.names(data)
    colnames(data_anno) <- colnames(data)
    
    anno <- SingleR(test = data_anno, ref = ref, 
                    clusters = unlist(data$seurat_clusters),
                    assay.type.test=1, labels = ref$label.main)
    cell_type_singleR <- plyr::mapvalues(x = data$seurat_clusters, from = row.names(anno), to = anno$labels)
    return(cell_type_singleR)
}
pg_data_combined$cluster <- function_singleR(pg_data_combined, "human")

# plot --
options(repr.plot.height = 6, repr.plot.width = 8)
p <- DimPlot(pg_data_combined, reduction = "umap", group.by = "spname") +
    coord_fixed()
ggsave(str_glue("{outdir}/spname.png"), plot = p, height = 5, width = 8)

p <- DimPlot(pg_data_combined, reduction = "umap", group.by = "gname") +
    coord_fixed()
ggsave(str_glue("{outdir}/gname.png"), plot = p, height = 5, width = 8)

p <- DimPlot(pg_data_combined, reduction = "umap", group.by = "seurat_clusters") +
    coord_fixed()
ggsave(str_glue("{outdir}/seurat_clusters.png"), plot = p, height = 5, width = 8)

p <- DimPlot(pg_data_combined, reduction = "umap", group.by = "clusters") +
    coord_fixed()
ggsave(str_glue("{outdir}/clusters.png"), plot = p, height = 5, width = 8)


# saveRDS
saveRDS(pg_data_combined, str_glue("{outdir}/data.rds"))