source activate kb_env

python /SGRNJ06/randd/USER/wangjingshen/script/virus_positive_fq/script/VirusPositiveFastq.py \
    --bam /SGRNJ06/randd/USER/wangjingshen/script_dev/test_data/HBV/H0808_1_Hep3T3_HBV_FJ/04.star_virus/H0808_1_Hep3T3_HBV_FJ_virus_Aligned.out.bam \
    --filter_umi_file /SGRNJ06/randd/USER/wangjingshen/script_dev/test_data/HBV/H0808_1_Hep3T3_HBV_FJ/06.filter_virus/H0808_1_Hep3T3_HBV_FJ_filtered_UMI.csv \
    --filter_read_count_json /SGRNJ06/randd/USER/wangjingshen/script_dev/test_data/HBV/H0808_1_Hep3T3_HBV_FJ/06.filter_virus/H0808_1_Hep3T3_HBV_FJ_filtered_read_count.json \
    --outdir outdir \
    --sample H0808_1_Hep3T3_HBV_FJ \
    --threads 8
