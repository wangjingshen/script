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


class Intergenic():
    def __init__(self, genome_fa, genome_gtf, gff, bam, outdir):
        '''
        '''
        self.genome_fa = genome_fa
        self.genome_gtf = genome_gtf
        self.gff = gff
        self.bam = bam
        self.outdir = outdir

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")


    def get_intergenic(self):
        cmd1 = f'cp {self.genome_fa} genome.fa '
        cmd2 = f'samtools faidx genome.fa '
        cmd3 = f'cut -f1 genome.fa.fai > genome.order.txt '
        cmd4 = (
            r'''awk '$3=="gene"{print $1,$4-1,$5}' OFS="\t" '''
            f'{self.genome_gtf} | sort -k1,1 -k2,2n > gene.bed'
        )
        cmd5 = f'bedtools sort -i gene.bed -g genome.order.txt > gene.sorted.bed'
        cmd6 = f'bedtools complement -i gene.sorted.bed -g genome.fa.fai > {self.outdir}/intergenic.bed '

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)
        subprocess.check_call(cmd4, shell=True)
        subprocess.check_call(cmd5, shell=True)
        subprocess.check_call(cmd6, shell=True)

    def get_promoter_enhancer(self):
        cmd1 = (
            f'cat {self.gff} | grep promoter | '
            r'''awk -F "\t" -v OFS="\t" '{print $1,$4,$5}' | '''
            f'sort -k1,1 -k2,2n > {self.outdir}/promoter.bed '
        )
        cmd2 = (
            f'cat {self.gff} | grep enhancer | '
            r'''awk -F "\t" -v OFS="\t" '{print $1,$4,$5}' | '''
            f'sort -k1,1 -k2,2n > {self.outdir}/enhancer.bed '
        )

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)

    def get_bam(self):
        cmd1 = f'bedtools intersect -a {self.bam} -b {self.outdir}/intergenic.bed -wa -u > {self.outdir}/intergenic.bam'
        cmd2 = f'bedtools intersect -a {self.bam} -b {self.outdir}/promoter.bed -wa -u > {self.outdir}/promoter.bam'
        cmd3 = f'bedtools intersect -a {self.bam} -b {self.outdir}/enhancer.bed -wa -u > {self.outdir}/enhancer.bam'
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)

    def rm_tmp(self):
        cmd1 = f'rm genome.fa genome.fa.fai genome.order.txt gene.bed gene.sorted.bed'
        subprocess.check_call(cmd1, shell=True)

    def run(self):
        self.get_intergenic()
        self.get_promoter_enhancer()
        self.get_bam()
        self.rm_tmp()


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep=r'\s+', header=0)
    genome_fa_list = df_mapfile.genome_fa
    genome_gtf_list = df_mapfile.genome_gtf
    gff_list = df_mapfile.gff
    bam_list = df_mapfile.bam
    outdir_list = df_mapfile.outdir
    
    return genome_fa_list, genome_gtf_list, gff_list, bam_list, outdir_list


def run_single(genome_fa, genome_gtf, gff, bam, outdir):
    runner = Intergenic(genome_fa, genome_gtf, gff, bam, outdir)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    genome_fa_list, genome_gtf_list, gff_list, bam_list, outdir_list = parse_mapfile(args.mapfile)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, genome_fa_list, genome_gtf_list, gff_list, bam_list, outdir_list)