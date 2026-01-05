suppressWarnings(suppressMessages({
    library(Seurat)
    library(argparser)
    library(tidyverse)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")
color_protocol2 <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101","OrangeRed","SlateBlue","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink","Red","#4682B4","#FFDAB9",
"#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00",
"#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B",
"#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A",
"#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF",
"#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE",
"#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D",
"#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32",
"#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", help = "rds")
argv <- add_argument(argv, "--diff_df", help = "diff_df")
argv <- add_argument(argv, "--subset_cluster", help = "subset cluster")
argv <- add_argument(argv, "--subset_group", help = "subset group")
argv <- add_argument(argv, "--plot_var", help = "plot_var, default: group")
argv <- add_argument(argv, "--width", help = "width, default:5")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 02.diff_group")
argv <- parse_args(argv)

rds <- argv$rds
diff_df <- argv$diff_df
plot_var <- ifelse(is.na(argv$plot_var), "group", argv$plot_var)
width <- ifelse(is.na(argv$width), 5, as.numeric(argv$width))
outdir <- argv$outdir

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# 
data_seurat <- readRDS(rds)

genus <- colnames(data_seurat@meta.data)[(which(colnames(data_seurat@meta.data) == "total_genus") + 1) : ncol(data_seurat@meta.data)]
df_genus <- cbind.data.frame(data_seurat@meta.data[,genus], 
                             data_seurat@meta.data[,c("cluster","group")])
df_genus_plot <- gather(df_genus, key = "genus", value = 'umi', -cluster, -group)

# subset
if(!is.na(argv$subset_cluster)){
    df_genus_plot <- df_genus_plot[ df_genus_plot$cluster %in% unlist(strsplit(argv$subset_cluster, split = ",")),]
}
if(!is.na(argv$subset_group)){
    df_genus_plot <- df_genus_plot[ df_genus_plot$group %in% unlist(strsplit(argv$subset_group, split = ",")),]
    plot_var = "cluster"
}

function_diff_barplot <- function(sub_g, outdir){
    genus_filter <- genus_diff_stat[ genus_diff_stat$p < 0.05 & genus_diff_stat$group1 == sub_g & genus_diff_stat$minus > 0 ,]$genus
    if(length(genus_filter)>0){
        if(length(genus_filter) > 21){
            color_use = color_protocol2
        }else{
            color_use = color_protocol[-13]
        }
        # plot
        if(plot_var == "group"){
            p <- df_genus_plot  %>%
                filter(genus %in% genus_filter) %>%
                group_by(genus, group) %>%
                summarise(mean_umi = mean(umi)) %>%
                ggplot(aes(group, mean_umi, fill = genus))+
                    geom_bar(color="black", stat="identity", width = 0.5)+
                    scale_fill_manual(values = color_use)+
                    theme_classic()+
                    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        }
        if(plot_var == "cluster"){
            p <- df_genus_plot  %>%
                filter(genus %in% genus_filter) %>%
                group_by(genus, cluster) %>%
                summarise(mean_umi = mean(umi)) %>%
                ggplot(aes(cluster, mean_umi, fill = genus))+
                    geom_bar(color="black", stat="identity", width = 0.5)+
                    scale_fill_manual(values = color_use)+
                    theme_classic()+
                    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        }
        ggsave(str_glue("{outdir}/barplot_diff_genus_{sub_g}.png"), plot = p, width = width, height = 5)
        ggsave(str_glue("{outdir}/barplot_diff_genus_{sub_g}.pdf"), plot = p, width = width, height = 5)
        print(paste0(sub_g, " done."))
    }
}
genus_diff_stat <- read.table(diff_df, sep="\t", header = T)
if(plot_var == "group"){
    out <- lapply(unique(df_genus_plot$group), function_diff_barplot, outdir = outdir)
}
if(plot_var == "cluster"){
    out <- lapply(unique(df_genus_plot$cluster), function_diff_barplot, outdir = outdir)
}

print("diff barplot done.")