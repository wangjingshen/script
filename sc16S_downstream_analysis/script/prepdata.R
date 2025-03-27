# ref: Tumor microbiome links cellular programs and immunity in pancreatic cancer

suppressWarnings(suppressMessages({
    library(Seurat)
    library(tidyverse)
    library(argparser)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")


# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--df_genus", help = "df_genus")
argv <- add_argument(argv, "--rna_spname", help = "rna sample name, split by ,")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 00.data")
argv <- parse_args(argv)

rds <- argv$rds
df_genus <- unlist(strsplit(argv$df_genus, split = ","))
rna_spname <- unlist(strsplit(argv$rna_spname, split = ","))
umi_threshold <- ifelse(is.na(argv$umi_threshold), 2, as.numeric(argv$umi_threshold))
outdir <- ifelse(is.na(argv$outdir), "00.data", argv$outdir)

print("df_genus: ")
print(df_genus)
print("rna sample: ")
print(rna_spname)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# read rds --
data_seurat <- readRDS(rds)
data_seurat$barcode <- colnames(data_seurat)
print(table(data_seurat$sample))

# read genus --
args_df <- data.frame(df_genus = df_genus,
                      rna_spname = rna_spname)

df_genus <- lapply(1:nrow(args_df), function(x){
    df <- t(read.table(args_df[x,1], sep="\t", row.names = 1, header = T))
    row.names(df) <- paste0(args_df[x,2], "_", row.names(df))   
    return(df)
})
all_genus <- Reduce(union, lapply(df_genus, colnames))

# fill all_genus for each df_genus
df_genus <- lapply(df_genus, function(x){
    df <- data.frame(genus = all_genus) %>%
        left_join(rownames_to_column(as.data.frame(t(x)), var = "genus") , by = c("genus" = "genus")) %>%
        t() %>%
        as.data.frame()
    colnames(df) <- df[1, ]
    df <- df[-1, ]
    df[is.na(df)] <- 0
    df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.numeric)  # character to numeric
    return(df)
})
df_genus <- do.call(rbind, df_genus)
#print(df_genus[1:2,1:2])
df_genus <- df_genus[ colnames(data_seurat),]
df_genus <- df_genus[, colSums(df_genus)>0]   # fix filter rds, some genus all 0
all_genus <- colnames(df_genus) # update all_genus

if(!identical(colnames(data_seurat), row.names(df_genus))){
    print("bc not identical, please check the input.")
    print(dim(data_seurat))
    print(dim(df_genus))
    print(colnames(data_seurat)[1:3])
    print(row.names(df_genus)[1:3])
    quit()
}else{
    print("bc identical.")
}

data_seurat$total_genus <- rowSums(df_genus)
data_seurat@meta.data <- cbind.data.frame(data_seurat@meta.data, df_genus)
saveRDS(data_seurat, str_glue('{outdir}/data_seurat.rds'))

#
# genus barcode count
#    g1 g2
# bc1
# bc2
genus_bc <- data.frame(genus = names(sort(colSums(df_genus > 0), decreasing = T)),
                       detect_barcode = as.numeric(sort(colSums(df_genus > 0), decreasing = T))) %>%
                       mutate(cell_number = ncol(data_seurat),
                              percent = paste(round(detect_barcode/cell_number*100, 2),"%"))
write.table(genus_bc, str_glue('{outdir}/genus_detect_barcode.tsv'), sep = "\t", quote = F, row.names = F)


# genus umi
genus_umi <- data.frame(genus = names(colSums(df_genus)),
                        umi = as.numeric(colSums(df_genus))) %>%
    arrange(desc(umi)) %>%
    mutate(total_umi = sum(umi),
           percent = paste(round(umi/total_umi*100, 2),"%"))
write.table(genus_umi, str_glue('{outdir}/genus_detect_umi.tsv'), sep = "\t", quote = F, row.names = F)