#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob
import scanpy as sc
import numpy as np

def get_intergenic(genome_fa, genome_gtf, outdir):
    cmd1 = f'cp {genome_fa} genome.fa '
    cmd2 = f'samtools faidx genome.fa '
    cmd3 = f'cut -f1 genome.fa.fai > genome.order.txt '
    cmd4 = (
        r'''awk '$3=="gene"{print $1,$4-1,$5}' OFS="\t" '''
        f'{genome_gtf} | sort -k1,1 -k2,2n > gene.bed'
    )
    cmd5 = f'bedtools sort -i gene.bed -g genome.order.txt > gene.sorted.bed'
    cmd6 = f'bedtools complement -i gene.sorted.bed -g genome.fa.fai > {outdir}/intergenic.bed '
    cmd7 = f'rm genome.fa genome.fa.fai genome.order.txt gene.bed gene.sorted.bed'
        
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    subprocess.check_call(cmd3, shell=True)
    subprocess.check_call(cmd4, shell=True)
    subprocess.check_call(cmd5, shell=True)
    subprocess.check_call(cmd6, shell=True)
    subprocess.check_call(cmd7, shell=True)


#def parse_mapfile(mapfile):
#    df_mapfile = pd.read_csv(mapfile, sep='\s+', header=None)
#    genome_fa_list = df_mapfile.iloc[:,0].tolist()
#    genome_gtf_list = df_mapfile.iloc[:,1].tolist()
#    outdir_list = df_mapfile.iloc[:,2].tolist()
#    return genome_fa_list, genome_gtf_list, outdir_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_fa', help='genome_fa', required=True)
    parser.add_argument('--genome_gtf', help='genome_gtf', required=True)
    parser.add_argument('--outdir', help='outdir', required=True)

    args = parser.parse_args()
    print(args)
    get_intergenic(args.genome_fa, args.genome_gtf, args.outdir)
    #genome_fa_list, genome_gtf_list, outdir_list = parse_mapfile(args.mapfile)
    #with ProcessPoolExecutor(max_workers = 1) as executor:
    #    executor.map(get_intergenic, genome_fa_list, genome_gtf_list, outdir_list)


if __name__ == '__main__':
    main()