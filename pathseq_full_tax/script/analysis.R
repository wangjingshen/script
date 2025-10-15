suppressWarnings(suppressMessages({
    library(tidyverse)
    library(argparser)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--pathseq_score", help = "pathseq_score")  # "celescope_16S/Control-16S/02.pathseq/Control-16S_pathseq_score.txt"
argv <- add_argument(argv, "--df_genus", help = "df_genus")  # "/SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.lims/P25062301_16S/Control-16S/e5d06535-3b6a-4f9f-8f62-4c4af3accd6a/call-pathseq_cellranger3/execution/Control-16S_EmptyDrops_CR_raw_UMI_matrix.tsv.gz"
argv <- add_argument(argv, "--name", help = "name")
argv <- add_argument(argv, "--outdir", help = "output dir, default: outdir")
argv <- parse_args(argv)

pathseq_score <- argv$pathseq_score
df_genus <- argv$df_genus
name <- argv$name
outdir <- ifelse(is.na(argv$outdir), "outdir", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

#
pathseq_score <- data.table::fread(pathseq_score, sep="\t", header = T, check.names = F)[,1:2]   # get tax_id, tax_name
df_genus <- read.table(df_genus, sep="\t", header = T, row.names = 1)
print("df_genus:")
print(dim(df_genus))
#sort(row.names(df_genus))
pathseq_df <- pathseq_score[ grepl(paste0(row.names(df_genus),collapse = "|"), pathseq_score$taxonomy),]

pathseq_df <- lapply(pathseq_df$taxonomy, function(x){
    test <- unlist(str_split(x, pattern = "\\|", simplify = F))
    for(i in length(test):1){
        if(grepl("_", test[i])==FALSE){
            res = test[1:i]
            break
        }
    }
    return(res)
})

df_list <- lapply(unique(sapply(pathseq_df, length)), function(x){
   return(pathseq_df[lapply(pathseq_df, length)==x])
})

names(df_list) <- paste0("df_", unique(sapply(pathseq_df, length)))


function_gettax<- function(data){
    data <- as.data.frame(t(as.data.frame(data)))    
    data <- data[ !duplicated(data[,ncol(data)]),]
    
    row.names(data) <- 1:nrow(data)
    
    data <- data[data[,ncol(data)] %in% row.names(df_genus),]
    #print(data[1:2,])
    if(nrow(data)>0){   # rm undect genus(not in df_Genus)
        return(data)
    }
}

df_tax <- lapply(df_list, function_gettax)

function_get_trans <- function(df, x){
    res <- data.frame(genus = df[x,ncol(df)],
                      replace = paste0(df[x,], collapse ="|"))
    return(res)
}
df_trans <- do.call(rbind, lapply(df_tax, function(x){
    if(!is.null(x)){
        do.call(rbind,lapply(1:nrow(x), function_get_trans, df=x))
    }
}))

# process unclass
class_genus <- as.character(unlist(sapply(df_tax, function(x){
    x[,ncol(x)]
})))

unclass_genus <- setdiff(row.names(df_genus), class_genus) 
#unclass_genus

# process unclass_genus no (
if(length(grep("(",unclass_genus, fixed = T, invert = T)) > 0){
    unclass_genus1 <- unclass_genus[grep("(",unclass_genus, fixed = T, invert = T)]
    print("unclass_genus1:")
    print(unclass_genus1)

    for (i in unclass_genus1){
        tmp <- pathseq_score$taxonomy[ grep(i, pathseq_score$taxonomy)]
        tmp <- lapply(tmp, function(x){
            unlist(str_split(x, pattern ="\\|",simplify =F))
        })
        min_length <- min(sapply(tmp, length))
        tmp <- tmp[[which(sapply(tmp, length) == min_length)]]
        df <- data.frame(genus = i,
                   replace = paste0(tmp, collapse = "|"))
        df_trans <- rbind(df_trans, df)
    }
}

# with (
if(length(grep("(",unclass_genus, fixed = T)) > 0){
    unclass_genus2 <- (unclass_genus[grep("(",unclass_genus, fixed = T)])
    print("unclass_genus2:")
    print(unclass_genus2)

    for (i in unclass_genus2){
        tmp <- pathseq_score$taxonomy[ grep(i, pathseq_score$taxonomy, fixed = T)]
        tmp <- lapply(tmp, function(x){
            unlist(str_split(x, pattern ="\\|",simplify =F))
        })
        min_length <- min(sapply(tmp, length))
        tmp <- tmp[[which(sapply(tmp, length) == min_length)]]
        df <- data.frame(genus = i,
                         replace = paste0(tmp, collapse = "|"))
        df_trans <- rbind(df_trans, df)
    }
}

df_trans <- df_trans[ match(row.names(df_genus), df_trans$genus),]
print(paste0("check genus of df_trans and row.names of df_genus:", identical(df_trans$genus, row.names(df_genus))))
#write.table(df_trans, str_glue("{outdir}/{name}_test.tsv"), sep="\t", quote=F, row.names=F)
row.names(df_genus) <- df_trans$replace
write.table(df_genus, str_glue("{outdir}/{name}_raw_UMI_matrix.tsv"),sep="\t",quote=F)


# check df trans
#test = sapply(1:nrow(df_trans), function(x){
#    test = unlist(str_split(df_trans[x,2], pattern = "\\|"))
#    return(test[length(test)])
#})
#identical(test, df_trans$genus)