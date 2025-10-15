suppressMessages({
    library(tidyverse)
    library(argparser)
})
options(scipen = 20)

argv <- arg_parser('')
argv <- add_argument(argv, "--virus_id_file", help = "barcode umi id tsv file, split by ','. The file can be get from get_virus_bam.py")
argv <- add_argument(argv, "--sample_name", help = "sample name, split by ','")
argv <- add_argument(argv, "--total_reads", help = "total reads count, split by ','")
argv <- add_argument(argv, "--downsample_dis", help = "downsample distance, 0.05 or 0.1, Default: 0.05")
argv <- add_argument(argv, "--outdir", help = "outdir to save fig. Default: './'")
argv <- parse_args(argv)
virus_id_file <- unlist(strsplit(argv$virus_id_file, split = ","))
total_reads <- unlist(strsplit(argv$total_reads, split = ","))
sample_name <- unlist(strsplit(argv$sample_name, split = ","))
downsample_dis <- ifelse(is.na(argv$downsample_dis), 0.05, as.numeric(argv$downsample_dis))
outdir <- ifelse(is.na(argv$outdir), './', argv$outdir)

# read data --
virus_id_list <- lapply(virus_id_file, function(x){
    data <- read.table(x, sep="\t")
    colnames(data) <- c("barcode_umi", "id")
    return(data)
})
names(virus_id_list) <- sample_name

# downsample --
function_downsample <- function(total_reads, percent, sample_name){
    data <- virus_id_list[[sample_name]]
    set.seed(100)
    downsample_id <- sample(as.numeric(total_reads), as.numeric(total_reads)*percent)
    downsample_data <- data[data$id %in% intersect(downsample_id, data$id),]
    return(length(unique(downsample_data$barcode_umi)))
}
downsample_res <- list()
for (i in 1:length(sample_name)){
    downsample_res[[sample_name[i]]]  <- sapply(seq(downsample_dis, 1, downsample_dis), function_downsample, total_reads = total_reads[i], sample_name = sample_name[i])
}
names(downsample_res) <- sample_name

plot_df <- do.call(rbind, lapply(sample_name, function(x){
    return(data.frame(percent = seq(downsample_dis, 1, downsample_dis),
                      umi_count = unlist(downsample_res[[x]]),
                      sample_name = x))
}))

# plot --
p <- plot_df %>%
    ggplot(aes(x = percent, y = umi_count, color = sample_name)) +
        geom_point(size=0.1, color= "white") +
        geom_smooth(se = F) +
        scale_y_continuous(name = "UMI count") +
        scale_x_continuous(name = "Sequencing depth") +
        theme_classic() + theme(legend.title = element_blank())

ggsave(str_glue('{outdir}/umi_sequencing_depth.pdf'), plot = p, height = 4, width = 6)
ggsave(str_glue('{outdir}/umi_sequencing_depth.png'), plot = p, height = 4, width = 6)
#write.table(plot_df, paste0(outdir, "/umi_sequencing_depth_plot_data.tsv"), sep = "\t", quote = F, row.names = F)