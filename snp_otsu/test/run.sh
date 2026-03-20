source activate celescope3.0.0

python /SGRNJ06/randd/USER/wangjingshen/script/snp_otsu/script/filter_snp.py \
  --vcf /SGRNJ06/randd/USER/zhouyiqi/2024/scsnp/2024-10-22-NPM1/outs/snpeff_snpeff/BL21605_Target.ann.vcf \
  --outdir outdir \
  --sample BL21605_Target \
  --bcftools_filter "QUAL>=100" \
  --ref_threshold_method otsu \
  --alt_threshold_method otsu \
  --ref_min_support_read 2 \
  --alt_min_support_read 2 \
  --vaf 0.2 \
  --subparser_assay snp \
  --thread 10
