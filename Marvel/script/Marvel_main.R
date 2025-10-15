
suppressPackageStartupMessages({
    # Load MARVEL package
    library("MARVEL")

    # Load adjunct packages for selected MARVEL features
    # General data processing, plotting
    library("ggnewscale")
    library("ggrepel")
    library("reshape2")
    library("plyr")
    library("stringr")
    library("textclean")
    library("tidyverse")

    # Gene ontology analysis
    library("AnnotationDbi")
    library("clusterProfiler")
    library("org.Hs.eg.db")
    library("org.Mm.eg.db")

    # ad hoc gene candidate gene analysis
    library("gtools")

    # Visualising splice junction location
    library("GenomicRanges")
    library("IRanges")
    library("S4Vectors")
    library("wiggleplotr")

    # Load adjunct packages for this tutorial
    library("Matrix")
    library("data.table")
    library("ggplot2")
    library("gridExtra")

    #
    library("Seurat")
    library("dplyr")
    library("argparser")
})

argv <- arg_parser('')
argv <- add_argument(argv,"--gene_mtx",help="path to the gene count matrix 10x")
argv <- add_argument(argv,"--sj_mtx",help="path to sj count matrix")
argv <- add_argument(argv,"--gtf",help="path to gtf file",default=NA)
argv <- add_argument(argv,"--outdir",help="path of output files")
argv <- add_argument(argv,"--rerun_cluster", help="bool,TRUE/FALSE", default=FALSE)
argv <- add_argument(argv,"--tsne_coord",help="required if param.rerun_cluster=FALSE",default=NA)
argv <- add_argument(argv,"--prefix",help="prefix of output files",default = 'NA')
argv <- add_argument(argv,"--sj_out_tab",help='SJ.out.tab file produced by starsolo')
argv <- add_argument(argv,"--species",help='homo/homo_mus/mus,ignored if gtf supplied ')
argv <- add_argument(argv,"--assay",help='assay')
argv <- add_argument(argv,"--sample_file",help='批量运行多个样本，可以提供sample_file,文件包含2列，1-celescopev2的输出路径，2-prefix,prefix不能有重复')
argv <- parse_args(argv)

gtf.file = argv$gtf
outdir = argv$outdir
rerun_cluster = argv$rerun_cluster
species = argv$species
assay = argv$assay


if ( ! dir.exists(outdir)){dir.create (outdir)}

if(is.na(gtf.file)){
    if(species=='homo'){
        gtf.file='/SGRNJ06/randd/public/genome/rna/celescope_v2/hs/Homo_sapiens.GRCh38.99.gtf'
        print('Homo_sapiens.GRCh38.99 used for analyze')
    }else if (species=='homo_mus') {
        gtf.file = '/SGRNJ06/randd/public/genome/rna/celescope_v2/hs_mmu/homo_mus.gtf'
        print('celscope v2 homo_mus genome used to analyze')
    }else if(species == 'mus'){
        gtf.file = '/SGRNJ06/randd/public/genome/rna/celescope_v2/mmu/Mus_musculus.GRCm38.99.gtf'
        print('celscope v2 Mus_musculus.GRCm38.99 used to analyze')
    }else{
        stop('no gtf detected!')
    }
}
gtf <- as.data.frame(fread(gtf.file), 
                    sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gtf) = c('V1','V2','V3','V4','V5','V6','V7','V8','V9')
gtf$V1 <-  str_replace_all(gtf$V1,'_','-')

run_marvel<-function(gtf,
                    gene.count.matrix.dir,
                    rerun_cluster,
                    tsne_coord,
                    sj.count.matrix.dir,
                    outdir,
                    sj_out_tab
                    ){
    outpath = paste0(outdir,'/',prefix)
    if( ! dir.exists(outpath)){dir.create (outpath)}
    sobj = CreateSeuratObject(Read10X(gene.count.matrix.dir),min.cells = 1,min.features = 1)
    sobj.norm <- NormalizeData(sobj, normalization.method = "RC",scale.factor = 10000)

    if(rerun_cluster){
        nfeatures=2000
        Dims =20
        resolution=0.6
        sobj.umap <- NormalizeData(sobj, normalization.method = "LogNormalize",scale.factor = 10000)
        sobj.umap <- FindVariableFeatures(sobj.umap, selection.method = "vst", nfeatures = nfeatures, mean.cutoff = c(0.1, 8), 
                                    dispersion.cutoff = c(1, Inf),
                                    mean.function = ExpMean, dispersion.function = LogVMR)
        use.genes <- sobj.umap@assays$RNA@var.features
        sobj.umap <- ScaleData(sobj.umap, vars.to.regress = c("nCount_RNA","percent.mito"), features = use.genes) #no percent.mito

        sobj.umap <- RunPCA(object = sobj.umap, features = use.genes, do.print = FALSE)

        sobj.umap <- FindNeighbors(sobj.umap, dims = 1:Dims, force.recalc = TRUE, reduction = "pca")
        sobj.umap <- FindClusters(sobj.umap, resolution = resolution)

        sobj.umap <- RunTSNE(sobj.umap, dims = 1:Dims,check_duplicates = FALSE)
        tsne_coord = as.data.frame(sobj.umap@reductions$tsne@cell.embeddings)
        saveRDS(sobj.umap,paste0(outpath,'/geneMatrix.rds'))

    }else{
        tsne_coord = read.table(tsne_coord_file,header = TRUE,row.names=1)
    }
    df.gene.norm <- sobj.norm@assays$RNA@data
    df.gene.norm.pheno <- sobj.norm@meta.data
    df.gene.norm.pheno["cell.id"] <- Cells(sobj.norm)
    df.gene.norm.feature = data.frame(gene_short_name = rownames(sobj.norm))

    df.gene.count <- sobj.norm@assays$RNA@counts
    df.gene.count.pheno <- sobj.norm@meta.data
    df.gene.count.pheno["cell.id"] <- Cells(sobj.norm)
    df.gene.count.feature = data.frame(gene_short_name = rownames(sobj.norm))

    sobj.sj = CreateSeuratObject(Read10X(sj.count.matrix.dir,gene.column = 1),min.cells = 1,min.features = 1)
    sobj.sj = subset(sobj.sj,cells = Cells(sobj))
    df.sj.count = sobj.sj@assays$RNA@counts
    df.sj.count.pheno = data.frame(cell.id = colnames(df.sj.count))
    df.sj.count.feature = data.frame(coord.intron = rownames(df.sj.count))

    df.coord = data.frame(cell.id = rownames(tsne_coord),
                        x = tsne_coord$tSNE_1,
                        y = tsne_coord$tSNE_2
                        )
    rownames(df.coord) = df.coord$cell.id

    marvel <- CreateMarvelObject.10x(gene.norm.matrix=df.gene.norm,
                                    gene.norm.pheno=df.gene.norm.pheno,
                                    gene.norm.feature=df.gene.norm.feature,
                                    gene.count.matrix=df.gene.count,
                                    gene.count.pheno=df.gene.count.pheno,
                                    gene.count.feature=df.gene.count.feature,
                                    sj.count.matrix=df.sj.count,
                                    sj.count.pheno=df.sj.count.pheno,
                                    sj.count.feature=df.sj.count.feature,
                                    pca=df.coord,
                                    gtf=gtf
                                    )
    marvel <- AnnotateGenes.10x(MarvelObject=marvel)
    marvel <- AnnotateSJ.10x(MarvelObject=marvel)
    write.table(as.data.frame(table(marvel$sj.metadata$sj.type)),paste0(outpath,'/sj.type.stat.txt'),sep = ',',quote = FALSE,row.names = FALSE)

    marvel_flt <- ValidateSJ.10x(MarvelObject=marvel)
    marvel_flt <- FilterGenes.10x(MarvelObject=marvel)
    marvel_flt <- CheckAlignment.10x(MarvelObject=marvel)
    save(marvel_flt, file=paste(outpath,'/',prefix,".MARVEL.RData", sep=""))

    marvel = marvel_flt
    sj_read = read.table(sj_out_tab)

    print(paste0("Cell Number: ",as.character(nrow(marvel$sample.metadata))))
    print(paste0("SJ Number: ",as.character(nrow(marvel$sj.metadata))))

    colSums(marvel$sj.count.matrix>0)->cell_sj
    print(paste0("Mean SJs per cell: ",as.character(mean(cell_sj))))
    print(paste0("Median SJs per cell: ",as.character(median(cell_sj))))

    colSums(marvel$sj.count.matrix)->umi_sj
    print(paste0("Mean UMI per cell: ",as.character(mean(umi_sj))))
    print(paste0("Median UMI per cell: ",as.character(median(umi_sj))))

    sj_read%>%mutate(sj_name = paste0('chr',V1,':',V2,':',V3),sj_read = V7+V8)->sj_read_1
    sj_read_1$sj_name = str_replace_all(sj_read_1$sj_name,'_','-')
    sj_read_1%>%filter(sj_name %in% marvel$sj.metadata$coord.intron)->sj_read_flt
    print(paste0("Mean reads per SJ(all): ",as.character(mean(sj_read_flt$sj_read))))
    print(paste0("Median reads per SJ(all): ",as.character(median(sj_read_flt$sj_read))))


    rowSums(marvel$sj.count.matrix>0)->df
    cbind(sj_read_flt,data.frame(cellNum=df))->sj_read_flt
    sj_read_flt%>%filter(cellNum>0)%>%mutate(mean_read_sj = sj_read/cellNum)->sj_read_flt

    print(paste0("Mean reads per SJ(expressed in cells): ",as.character(mean(sj_read_flt$sj_read))))
    print(paste0("Median reads per SJ(expressed in cells): ",as.character(median(sj_read_flt$sj_read))))
    
    print(paste0("Mean reads per SJ per cell: ",as.character(mean(sj_read_flt$mean_read_sj))))
    print(paste0("Median reads per SJ per cell: ",as.character(median(sj_read_flt$mean_read_sj))))



    data.frame(`Cell numeber` = nrow(marvel$sample.metadata),
            `SJ number` = nrow(marvel$sj.metadata),
            `Mean SJs per cell` = mean(cell_sj),
            `Median SJs per cell` = median(cell_sj),
            `Mean UMIs per cell` = mean(umi_sj),
            `Median UMIs per cell` = median(umi_sj),
            `Mean reads per SJ` = mean(sj_read_flt$sj_read),
            `Median reads per SJ` = median(sj_read_flt$sj_read),
            `Mean reads per SJ per cell` = mean(sj_read_flt$mean_read_sj),
            `Median reads per SJ per cell` = median(sj_read_flt$mean_read_sj)
            )->stat_df
    data.frame(Categories = c('Cell numeber','SJ number','Mean SJs per cell','Median SJs per cell','Mean UMIs per cell','Median UMIs per cell','Mean reads per SJ','Median reads per SJ','Mean reads per SJ per cell','Median reads per SJ per cell'),
    Values = c(nrow(marvel$sample.metadata),nrow(marvel$sj.metadata),mean(cell_sj),median(cell_sj),mean(umi_sj),median(umi_sj),mean(sj_read_flt$sj_read),median(sj_read_flt$sj_read),mean(sj_read_flt$mean_read_sj),median(sj_read_flt$mean_read_sj)))->stat_df
    
    #stat_df = t(stat_df)
    print(stat_df)
    #colnames(stat_df) = c('Categories','Value')
    write.table(stat_df,paste0(outpath,'/',prefix,'.stat.txt'),sep = "\t",quote = FALSE,row.names =FALSE)
    write.table(stat_df,paste0(outpath,'/summary.txt'),sep = "\t",quote = FALSE,row.names = FALSE)

    sample_info = data.frame(`Sample ID` = prefix,`Species` = species,`Assay` = assay)
    write.table(sample_info,paste0(outpath,'/sample_info.txt'),sep = "\t",quote = FALSE,row.names = FALSE)
}

if(is.na(argv$sample_file)){
    gene.count.matrix.dir = argv$gene_mtx
    sj.count.matrix.dir = argv$sj_mtx
    tsne_coord_file = argv$tsne_coord
    prefix = argv$prefix
    sj_out_tab = argv$sj_out_tab

    run_marvel(gtf,gene.count.matrix.dir,rerun_cluster,tsne_coord,sj.count.matrix.dir,outdir,sj_out_tab)

}else{
    sample_list=read.table(argv$sample_file)
    sname=unlist(lapply(sample_list[,1],function(x){.y=str_split(x,'/')[[1]];.y=.y[.y!=""];return(.y[length(.y)])}))

    for(i in 1:length(sname)){
       
        gene.count.matrix.dir = paste0(sample_list[i,1],'/outs/filtered')
        sj.count.matrix.dir = paste0(outdir,'/SJ_mtx/',sample_list[i,2])
        dir.create (paste0(outdir,'/SJ_mtx'))
        dir.create (paste0(outdir,'/SJ_mtx/',sample_list[i,2]))
        solo_path = paste0(sample_list[i,1],'/01.starsolo/',sname[i],'_Solo.out')
        command=paste0("bash /SGRNJ06/randd/USER/fuxin/RD/FullLength/src/generate_sj_mtx.sh ",solo_path," ",sj.count.matrix.dir)
        tsne_coord_file = paste0(sample_list[i,1],'/outs/tsne_coord.tsv')
        prefix = sample_list[i,2]
        sj_out_tab = paste0(sample_list[i,1],'/01.starsolo/',sname[i],'_SJ.out.tab')

        system(command,
                intern = FALSE,ignore.stdout = FALSE,ignore.stderr = FALSE,
                wait = TRUE, input = NULL, show.output.on.console = TRUE,
                minimized = FALSE, invisible = TRUE, timeout = 0)
        run_marvel(gtf,gene.count.matrix.dir,rerun_cluster,tsne_coord,sj.count.matrix.dir,outdir,sj_out_tab)

        run_html_cmd = paste0("source activate target_env;python /SGRNJ06/randd/USER/wangjingshen/script_dev/report/run_html.py  --acc_out_dir ",outdir,"/",prefix,"  --out_dir  ",outdir,"/",prefix)
        system(run_html_cmd,
                intern = FALSE,ignore.stdout = FALSE,ignore.stderr = FALSE,
                wait = TRUE, input = NULL, show.output.on.console = TRUE,
                minimized = FALSE, invisible = TRUE, timeout = 0)

    }
}
