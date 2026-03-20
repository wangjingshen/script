source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/pathseq_full_tax/script/analysis.R \
    --pathseq_score /SGRNJ06/randd/PROJECT/P25062301/B1/r0_check/celescope_16S/Control-16S/02.pathseq/Control-16S_pathseq_score.txt \
    --df_genus /SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.lims/P25062301_16S/Control-16S/e5d06535-3b6a-4f9f-8f62-4c4af3accd6a/call-pathseq_cellranger3/execution/Control-16S_EmptyDrops_CR_raw_UMI_matrix.tsv.gz \
    --name Control \
    --outdir outdir