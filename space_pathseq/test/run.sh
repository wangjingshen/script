source activate r4.1_env

python /SGRNJ06/randd/USER/wangjingshen/script/space_pathseq/script/pipeline.py \
    --space_dir /SGRNJ06/randd/PROJECT/R25030501_Spatial_FFPE_tgx/20251127/Mint_FFPE_96_96_1119/ \
    --pathseq_df Mint_FFPE_96_96_1119/outs/Mint_FFPE_96_96_1119_raw_UMI_matrix.tsv.gz \
    --topn_genus 5\
    --outdir Mint_FFPE_96_96_1119 \
