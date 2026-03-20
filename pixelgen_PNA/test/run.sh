source activate pixelgenr_env

Rscript /SGRNJ06/randd/USER/wangjingshen/script/pixelgen_PNA/script/analysis.R \
    --pxl /SGRNJ06/randd/USER/wangjingshen/rd_project/2025/mpx/r1/P25112507/S01_pixelgen-singleron-P25112507.layout.pxl,/SGRNJ06/randd/USER/wangjingshen/rd_project/2025/mpx/r1/P25112507/S02_pixelgen-singleron-P25112507.layout.pxl \
    --sample S01,S02 \
    --group S01,S02 \
    --nUMI_cutoff 10000 \
    --isotype_fraction_cutoff 0.001 \
    --rm_batch no \
    --rm_batch_var spname \
    --ndims 10 \
    --resolution 0.8 \
    --outdir outdir \
    --spatial_vis T 