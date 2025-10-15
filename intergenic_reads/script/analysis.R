suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--snp_dir", help = "snp dir")
argv <- parse_args(argv)

snp_dir <- argv$snp_dir

df <- read.table(str_glue("{snp_dir}/UMI_merge_all.txt"), sep="\t", header = T)

df <- df %>% 
    group_by(sample_name,match_name,protein) %>%
    summarise(ncell = sum(ncell),
              total = sum(total)) %>%
    mutate(percent = paste0(round(100*ncell/total, 2)," %"))

write.table(df, str_glue("{snp_dir}/UMI_merge_all_sum.txt"), sep="\t", row.names=F, quote=F)