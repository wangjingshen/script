source activate kb_env

python /SGRNJ06/randd/USER/wangjingshen/script/intergenic_reads/script/pipeline.py \
    --mapfile mapfile \
    --genome_fa /SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/Genome/celescope_v2/mmu_ensembl_110/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --genome_gtf /SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/Genome/celescope_v2/mmu_ensembl_110/Mus_musculus.GRCm39.110.mkgtf.gtf \
    --gff /SGRNJ06/randd/USER/wangjingshen/rd_project/2025/RD24012902_FFPE/20250926_r9/ensembl_mus_musculus/v110/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20221007.gff \
    --minimum_overlap 0.5
