#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
    library(future)
    library(furrr)
    library(logger)
    library(glue)
}))



# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--spname", help = "sample name in rds, split by ,")
argv <- add_argument(argv, "--tsne_tag", help = "tsne_tag, split by ,")
argv <- add_argument(argv, "--tag_name", help = "tag name")
argv <- add_argument(argv, "--VlnPlot", help = "VlnPlot, default: F")
argv <- add_argument(argv, "--FeaturePlot", help = "FeaturePlot, default:F")
argv <- add_argument(argv, "--threads", help = "threads, default: 1")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

rds <- argv$rds
spname <- unlist(strsplit(argv$spname, split = ","))
tsne_tag <- unlist(strsplit(argv$tsne_tag, split = ","))

if (length(spname) != length(tsne_tag)) {
    log_error("The number of spname and tsne_tag is inconsistent ！")
    quit()
}

tag_name <- argv$tag_name
VlnPlot <- ifelse(is.na(argv$VlnPlot), "F", argv$VlnPlot)
FeaturePlot <- ifelse(is.na(argv$FeaturePlot), "F", argv$FeaturePlot)

threads <- ifelse(is.na(argv$threads), 1, as.numeric(argv$threads))
plan(multisession, workers = threads)   # use multi

outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
if(!dir.exists(outdir)){
    log_info(glue("Creating output directory: {outdir}."))
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# begin
log_info("Loading Seurat object...")
data_seurat <- readRDS(rds)
data_seurat$sample <- data_seurat$`Sample ID`
data_seurat$cluster <- data_seurat$annot_full

# VlnPlot --
if(VlnPlot %in% c("TRUE","True","T")){
    log_info("VlnPlot start...")
    function_test <- function(i){
        tsne_df <- input_df[i, 1]
        spname <- input_df[i, 2]
        tsne_df <- read.table(tsne_df, sep="\t",header = T, row.names = 1)
        row.names(tsne_df) <- paste0(spname, "_", row.names(tsne_df))
        return(tsne_df)
    }
    input_df <- data.frame(tsne_tag = tsne_tag, spname = spname)
    tag_df <- do.call(rbind, lapply(1:nrow(input_df), function_test))

    # Check consistency between Seurat object and tag file
    missing_cells <- setdiff(colnames(data_seurat), rownames(tag_df))
    if (length(missing_cells) > 0) {
        log_error(glue(
            "Inconsistency detected: {length(missing_cells)} cells in Seurat object are missing in tag file. ",
            "Example missing cells: {paste(head(missing_cells), collapse = ', ')}"
        ))
        quit(save = "no", status = 1)
    }

    tag_df <- tag_df[ colnames(data_seurat),]

    data_seurat$tag <- tag_df[colnames(data_seurat), tag_name] %>% as.numeric()
    data_seurat$log2_tag <- log2(data_seurat$tag + 1)
    data_seurat$cluster_sample <- paste0(data_seurat$cluster, "_", data_seurat$sample)

    # VlnPlot
    function_VlnPlot <- function(plot_var, cutoff, name){
        p <- VlnPlot(data_seurat, plot_var, group.by = "cluster_sample", pt.size = 0) + 
            theme(legend.position = "none", axis.title.x = element_blank(),plot.margin = margin(10, 10, 10, 50, "pt")) +
            ggtitle(tag_name) + 
            theme(plot.title = element_text(hjust = 0.5))
        if(!is.na(cutoff)){
            p <- p + ylim(c(0,10000))
        }
        ggsave(str_glue("{outdir}/VlnPlot_{tag_name}{name}.png"), plot = p, height = 8, width = 12)
        ggsave(str_glue("{outdir}/VlnPlot_{tag_name}{name}.pdf"), plot = p, height = 8, width = 12)
    }
    function_VlnPlot("tag", cutoff=NA, name="")
    function_VlnPlot("log2_tag", cutoff=NA, name="_log2")
    function_VlnPlot("tag", cutoff=10000, name="_cutoff_10k")

    log_success("Vlnplot completed successfully.")
}


# FeaturePlot --
if(FeaturePlot %in% c("TRUE","True","T")){
    log_info("FeaturePlot start...")
    function_plot <- function(data, tsne_tag, spname, outdir, tag_name){
        logger::log_info("Processing sample: {spname}")

        tag_df <- read.table(tsne_tag, sep="\t", header = T, row.names = 1)
        row.names(tag_df) <- paste0(spname, "_", row.names(tag_df))
        tryCatch({
            data_sub <- subset(data, subset = sample == spname)
        },error = function(e){
            log_error(glue("Sample is inconsistent with that in the RNA: {spname}"))
            quit()
        })
    
        tag_vec <- tag_df[colnames(data_sub), tag_name] %>% as.numeric()
        data_sub$tag <- if (max(tag_vec, na.rm = TRUE) > 100) log2(tag_vec + 1) else tag_vec

        jet256 <- colorRampPalette(c("#00008F","#0000FF","#0080FF","#00FFFF","#80FF80","#FFFF00","#FF8000","#FF0000","#800000"))(256)

        # get legend 
        max_lab <- round(max(data_sub$tag, na.rm = TRUE))
        min_lab <- round(min(data_sub$tag, na.rm = TRUE))
        mid_lab <- (max_lab + min_lab)/2

        p <- FeaturePlot(data_sub, "tag", reduction = "umap", cols = jet256)
        p <- p + scale_colour_gradientn(colours = jet256, 
                                        limits  = c(0, 250), 
                                        breaks  = c(0, 125, 250), 
                                        labels  = c(min_lab, mid_lab, max_lab)) + 
            guides(colour = guide_colourbar(draw.ulim = TRUE, frame.colour = "black")) +
            ggtitle(tag_name) +
            theme(plot.title = element_text(hjust = 0.5))

        ggsave(str_glue("{outdir}/{spname}_{tag_name}.png"), plot = p, height = 6, width = 7)
        ggsave(str_glue("{outdir}/{spname}_{tag_name}.pdf"), plot = p, height = 6, width = 7)
    
        log_info(glue("Successfully processed sample: {spname}."))
        return(invisible(TRUE))  # suppress output in lapply
    }

    input_df <- data.frame(tsne_tag = tsne_tag, spname = spname)
    input_df %>%
        mutate(res = future_map2(tsne_tag, spname, function_plot,
                                 data = data_seurat, outdir = outdir, tag_name = tag_name,
                                 .options = furrr_options(seed = NULL))) %>%
        invisible()

    log_success("Script completed successfully.")
}