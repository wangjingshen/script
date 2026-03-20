source activate r4.1_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/HBV_UMI_depth/script/plot.R \
    --virus_id_file outdir/H0118_4_virus_id.tsv,outdir/hs0912_virus_id.tsv \
    --sample_name H0118_4,hs0912 \
    --total_reads 6313889,52407768 \
    --downsample_dis 0.05 \
    --outdir outdir
