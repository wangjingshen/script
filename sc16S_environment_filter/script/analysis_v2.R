##
# support celescope

suppressWarnings(suppressMessages({
    library(argparser)
    library(tidyverse)
}))

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--fj_path", help = "fj path")
argv <- add_argument(argv, "--sample", help = "sample")
argv <- parse_args(argv)

fj_path <- argv$fj_path
sample_list <- unlist(strsplit(argv$sample, split = ","))

# DOI: 10.1186/s12915-014-0087-z
environment_bac <- c(
"Proteobacteria",
"Alpha-proteobacteria",
"Afipia","Aquabacterium","Asticcacaulis","Aurantimonas","Beijerinckia","Bosea","Bradyrhizobium","Brevundimonas",
"Caulobacter","Craurococcus","Devosia","Hoeflea","Mesorhizobium","Methylobacterium","Novosphingobium","Ochrobactrum",
"Paracoccus","Pedomicrobium","Phyllobacterium","Rhizobium","Roseomonas","Sphingobium","Sphingomonas","Sphingopyxis",

"Beta-proteobacteria",
"Acidovorax","Azoarcus","Azospira","Burkholderia","Comamonas","Cupriavidus","Curvibacter","Delftia","Duganella",
"Herbaspirillum","Janthinobacterium","Kingella","Leptothrix","Limnobacter","Massilia","Methylophilus","Methyloversatilis", 
"Oxalobacter","Pelomonas","Polaromonas","Ralstonia","Schlegelella","Sulfuritalea","Undibacterium","Variovorax",

"Gamma-proteobacteria",
"Acinetobacter","Enhydrobacter","Enterobacter","Escherichia","Nevskia","Pseudomonas","Pseudoxanthomonas","Psychrobacter",
"Stenotrophomonas","Xanthomonas",

"Actinobacteria",
"Aeromicrobium","Arthrobacter","Beutenbergia","Brevibacterium","Corynebacterium","Curtobacterium","Dietzia","Geodermatophilus",
"Janibacter","Kocuria","Microbacterium","Micrococcus","Microlunatus","Patulibacter","Propionibacterium","Rhodococcus","Tsukamurella",

"Firmicutes",
"Abiotrophia","Bacillus","Brevibacillus","Brochothrix","Facklamia","Paenibacillus","Streptococcus",

"Bacteroidetes",
"Chryseobacterium","Dyadobacter","Flavobacterium","Hydrotalea","Niastella","Olivibacter","Pedobacter","Wautersiella",

"Deinococcus Thermus",
"Deinococcus",

"Acidobacteria",
"Predominantly unclassified Acidobacteria Gp2 organisms"
)

function_filter <- function(fj_path, sample){
    # mkdir
    if(!dir.exists(str_glue('01.raw_mtx/{sample}'))){
        dir.create(str_glue('01.raw_mtx/{sample}'), recursive = TRUE)
    }
    if(!dir.exists(str_glue('02.filter_mtx/{sample}'))){
        dir.create(str_glue('02.filter_mtx/{sample}'), recursive = TRUE)
    }

    ## raw
    genus_df <- read.table(str_glue('{fj_path}/{sample}/outs/{sample}_raw_UMI_matrix.tsv.gz'), sep = "\t", row.names = 1, header = T)
    genus_df <- as.data.frame(t(genus_df))
    write.table(cbind.data.frame(barcode = row.names(genus_df),
                                 genus_df), 
               str_glue('01.raw_mtx/{sample}/{sample}_genus_umi.tsv'), sep = "\t", quote = F, row.names = F)
    
    # genus barcode count
    genus_bc <- data.frame(genus = names(sort(colSums(genus_df > 0), decreasing = T)),
                           detect_barcode = as.numeric(sort(colSums(genus_df > 0), decreasing = T))) %>%
                           mutate(cell_number = nrow(genus_df),
                                  percent = paste(round(detect_barcode/cell_number*100, 2), "%"))
    write.table(genus_bc, str_glue('01.raw_mtx/{sample}/{sample}_genus_detect_barcode.tsv'), sep = "\t", quote = F, row.names = F)

    # genus umi
    genus_umi <- data.frame(genus = names(colSums(genus_df)),
                            umi = as.numeric(colSums(genus_df))) %>%
        arrange(desc(umi)) %>%
        mutate(total_umi = sum(umi),
               percent = paste(round(umi/total_umi*100, 2), "%"))
    write.table(genus_umi, str_glue('01.raw_mtx/{sample}/{sample}_genus_detect_umi.tsv'), sep = "\t", quote = F, row.names = F)


    ## filter
    genus_df_filter <- genus_df[ ,setdiff(colnames(genus_df), environment_bac)]

    write.table(cbind.data.frame(barcode = row.names(genus_df_filter),
                                 genus_df_filter), 
               str_glue('02.filter_mtx/{sample}/{sample}_genus_umi.tsv'), sep = "\t", quote = F, row.names = F)

    # genus barcode count
    genus_bc <- data.frame(genus = names(sort(colSums(genus_df_filter > 0), decreasing = T)),
                           detect_barcode = as.numeric(sort(colSums(genus_df_filter > 0), decreasing = T))) %>%
                           mutate(cell_number = nrow(genus_df_filter),
                                  percent = paste(round(detect_barcode/cell_number*100, 2), "%"))
    write.table(genus_bc, str_glue('02.filter_mtx/{sample}/{sample}_genus_detect_barcode.tsv'), sep = "\t", quote = F, row.names = F)

    # genus umi
    genus_umi <- data.frame(genus = names(colSums(genus_df_filter)),
                            umi = as.numeric(colSums(genus_df_filter))) %>%
        arrange(desc(umi)) %>%
        mutate(total_umi = sum(umi),
               percent = paste(round(umi/total_umi*100, 2), "%"))
    write.table(genus_umi, str_glue('02.filter_mtx/{sample}/{sample}_genus_detect_umi.tsv'), sep = "\t", quote = F, row.names = F)
}

input_df <- data.frame(fj_path = fj_path,
                       sample = sample_list)

out <- lapply(1:nrow(input_df), function(x){
    function_filter(input_df[x,1], input_df[x,2])
})

print("Filter environment done.")