# circos each group

suppressMessages(suppressWarnings({
    library(Seurat)
    library(tidyverse)
    library(argparser)
    library(circlize)
}))
color_protocol <- c("#0067AA","#FF7F00","#00A23F","#FF1F1D","#A763AC","#B45B5D","#FF8AB6","#B6B800","#01C1CC","#85D5F8","#FFC981","#C8571B","#727272","#EFC800","#8A5626","#502E91","#59A4CE","#344B2B","#FBE29D","#FDD6E6","#849C8C","#F07C6F","#000101")

# args --
argv <- arg_parser('')
argv <- add_argument(argv, "--stat_df", help = "stat df")
argv <- add_argument(argv, "--obj", help = "clusterA or groupA or all")
argv <- add_argument(argv, "--top_n", help = "top_n genus, dafault: 10")
argv <- add_argument(argv, "--outdir", help = "output dir, default: 02.diff/circos/")
argv <- parse_args(argv)

stat_df <- argv$stat_df
obj <- argv$obj
top_n <- ifelse(is.na(argv$top_n), 10, as.numeric(argv$top_n))
outdir <- ifelse(is.na(argv$outdir), "02.diff/circos/", argv$outdir)

if(!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}

# read df
stat_df <- read.table(stat_df, sep="\t", header = T)

function_circos_plot <- function(analysis_obj){
    stat_df_sub <- stat_df[ stat_df$group1 == analysis_obj & stat_df$p < 0.05 & stat_df$minus > 0,]   # get sig + postive
    if(nrow(stat_df_sub) > 0){
        stat_df_sub <- arrange(stat_df_sub, desc(mean1)) # desc to get top_n
        if(nrow(stat_df_sub) > top_n){
            stat_df_sub <- stat_df_sub[1:top_n,]
        }
        stat_df_sub <- arrange(stat_df_sub, mean1) # arrange to plot
        df1 <- stat_df_sub[, c("genus","group1","mean1")]
        df2 <- stat_df_sub[, c("genus","group2","mean2")]
        colnames(df1) <- c("from","to","value")
        colnames(df2) <- c("from","to","value")
        df <- rbind(df1, df2)
        row.names(df) <- paste0(df$from, df$to)
        #return(df)
    
        # plot
        set.seed(123)
        pdf(str_glue("{outdir}/circos_{analysis_obj}.pdf"), width = 8, height = 8, bg='white')
        grid.col <- setNames(color_protocol[1:length(union(df$to, df$from))], union(df$to, df$from))   # fix cols
        circos.par(track.height = 0.3)
        chordDiagram(df, annotationTrack = "grid", preAllocateTracks = 1,
            annotationTrackHeight = c(0.08),
            grid.col = grid.col, 
            link.decreasing = FALSE)

        circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
            xlim = get.cell.meta.data("xlim")
            sector.name = get.cell.meta.data("sector.index")
            n<-length(sector.name)
            for (i in seq_len(n)) {
                circos.text(mean(xlim), 0, sector.name[i], adj = c(0, 0.5),
                            facing = "clockwise", niceFacing = TRUE, cex = 0.5)
            }
        }, bg.border = NA)
        dev.off()
    }else{
        cat(paste0(analysis_obj, " skip. \n"))
    }
}

if(obj == "all"){
    res <- lapply(unique(stat_df$group1), function_circos_plot)
}else{
    res <- function_circos_plot(obj)
}
file.copy("/SGRNJ06/randd/USER/wangjingshen/script_dev/sc16S_downstream_analysis/doc/circos.txt", str_glue("{outdir}/README.txt"))
cat("circos plot done. \n")