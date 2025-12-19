#!/usr/bin/env python

import pysam
import sys
import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

def snp_vaf(snp_path, name):
    # mkdir
    #os.makedirs(f"{name}", exist_ok=True) 
    Path(f"{name}").mkdir(parents=True, exist_ok=True)
    
    final_vcf = glob.glob(f'{snp_path}/09.analysis_snp/{name}_final.vcf')[0]
    with pysam.VariantFile(final_vcf) as vcf:
        header = vcf.header
        header.formats.add("VAF", 1, "Float", "Variant allele frequency")
        with pysam.VariantFile(f"{name}/{name}_final_vaf.vcf", "w", header=header) as out:
            for rec in vcf:
                for sam in rec.samples:
                    ad = rec.samples[sam]["AD"]
                    vaf = ad[1] / (ad[0] + ad[1]) if sum(ad) > 0 else None
                    rec.samples[sam]["VAF"] = vaf
                out.write(rec)
        
def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\s+")
    snp_path_list = df_mapfile['snp_path']
    sample_list = df_mapfile['sample']
    return snp_path_list, sample_list

def main():
    mapfile = sys.argv[1]
    snp_path_list, sample_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(snp_vaf, snp_path_list, sample_list)

if __name__ == '__main__':
    main()