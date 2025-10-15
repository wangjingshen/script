suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--sample", help = "sample")
argv <- add_argument(argv, "--pathseq_tsv", help = "pathseq tsv")  # celescope_16S/Control-16S/02.pathseq/Control-16S_pathseq.tsv
argv <- add_argument(argv, "--downsample_tsv", help = "downsample tsv")  # celescope_16S/Control-16S/02.pathseq/Control-16S_downsample.tsv
argv <- add_argument(argv, "--tax_id", help = "tax_id")
argv <- add_argument(argv, "--tax_name", help = "tax_name")
argv <- add_argument(argv, "--match_dir", help = "match_dir")
argv <- add_argument(argv, "--outdir", help = "output dir")
argv <- parse_args(argv)

sample <- argv$sample
pathseq_tsv <- argv$pathseq_tsv
downsample_tsv <- argv$downsample_tsv
tax_id <- argv$tax_id
tax_name <- argv$tax_name
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)
match_dir <- argv$match_dir
match_bc <- str_glue("{match_dir}/outs/filtered/barcodes.tsv.gz")

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# get YP from pathseq
data_pathseq <- read.table(pathseq_tsv, sep="\t")
data_pathseq <- data_pathseq[ grep("^YP", data_pathseq$V2),]
data_pathseq$V2 <- gsub("YP:Z:","",data_pathseq$V2)

# 
tax_id <- read.table(tax_id, sep="\t")
function_get_tax_id <- function(x){
    list1 <- unlist(str_split(data_pathseq$V2[x], pattern = ",", simplify = F))
    return(length(intersect(list1, tax_id$V1)))
}
data_pathseq$n_tax_id <- unlist(lapply(1:nrow(data_pathseq), function_get_tax_id))
data_pathseq <- data_pathseq[ data_pathseq$n_tax_id > 0,]   # filter 0 tax_id

#  bc umi
data_downsample <- read.table(downsample_tsv, sep="\t")
data_downsample <- data_downsample[ match(data_pathseq$V1, data_downsample$V1), ]   # match reads_id
print(identical(data_pathseq$V1, data_downsample$V1))

data_downsample$barcode <- gsub("CR:Z:", "", data_downsample$V2)
data_downsample$V3 <- gsub("UR:Z:", "", data_downsample$V3)     # umi

# summarise
df_tax_id <- data_downsample %>% 
    group_by(barcode) %>%
    summarise(umi = length(unique(V3)),
              reads = length(V3))

# match bc
match_bc <- read.table(match_bc)[,1]

df_tax_id <- merge(data.frame(barcode = match_bc), df_tax_id,  by = "barcode", all.x = T)
df_tax_id [ is.na(df_tax_id)] <- 0
df_tax_id <- df_tax_id[ match(match_bc, df_tax_id$barcode),]

write.table(df_tax_id, str_glue("{outdir}/{sample}_{tax_name}.tsv"), row.names=F, sep="\t", quote=F)