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
import pysam

from gff2bed import get_promoter_enhancer
from intergenic_bed import get_intergenic

def get_reads_num(bam):
    '''
    #samtools view -F 0x900 common/promoter.bam | cut -f1 | sort -u | wc -l
    '''
    with pysam.AlignmentFile(bam, "rb") as bam:
        unique_ids = {read.query_name for read in bam if not (read.flag & 0x900)}
    return(len(unique_ids))


class Intergenic_reads():
    def __init__(self, genome_fa, genome_gtf, gff, bam, spname, minimum_overlap, step):
        '''
        '''
        self.genome_fa = genome_fa
        self.genome_gtf = genome_gtf
        self.gff = gff
        self.bam = bam
        self.spname = spname
        self.minimum_overlap = minimum_overlap
        self.step = step.strip().split(',')
        #print(self.gff)
        if not os.path.exists(self.spname):
            os.system(f"mkdir -p {self.spname}")


    def get_intergenic(self):
        get_intergenic(self.genome_fa, self.genome_gtf, self.spname)


    def get_promoter_enhancer(self):
        get_promoter_enhancer(self.gff, self.spname)


    def get_bam(self):
        cmd1 = f'bedtools intersect -a {self.bam} -b {self.spname}/intergenic.bed -wa -u -f {self.minimum_overlap} > {self.spname}/intergenic.bam'
        cmd2 = f'bedtools intersect -a {self.spname}/intergenic.bam -b {self.spname}/promoter.bed -wa -u -f {self.minimum_overlap} > {self.spname}/promoter.bam'
        cmd3 = f'bedtools intersect -a {self.spname}/intergenic.bam -b {self.spname}/enhancer.bed -wa -u -f {self.minimum_overlap} > {self.spname}/enhancer.bam'
        #cmd4 = f'bedtools intersect -a {self.spname}/promoter.bam -b {self.spname}/enhancer.bed -wa -u -f {self.minimum_overlap} > {self.spname}/promoter_enhancer.bam'
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)


    def get_stat_df(self):
        intergenic_reads = get_reads_num(f'{self.spname}/intergenic.bam')
        promoter_reads = get_reads_num(f'{self.spname}/promoter.bam')
        enhancer_reads = get_reads_num(f'{self.spname}/enhancer.bam')

        df = pd.DataFrame({
            'spname' : f'{self.spname}',
            'type': ['intergenic', 'promoter', 'enhancer', 'others'],
            'reads': [intergenic_reads, promoter_reads, enhancer_reads, intergenic_reads - (promoter_reads + enhancer_reads)]
            })
        df.to_csv(f"{self.spname}/intergenic_reads_number.tsv", index = False, sep="\t")

    def plot(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script_dev/intergenic_reads/script/plot.R '
            f'--df {self.spname}/intergenic_reads_number.tsv '
            f'--spname {self.spname}'
            )


    def run(self):
        if 'intergenic' in self.step:
            self.get_intergenic()
        if 'pr_en' in self.step:
            self.get_promoter_enhancer()
        if 'bam' in self.step:
            self.get_bam()
        if 'stat_df' in self.step:
            self.get_stat_df()
        if 'plot' in self.step:
            self.plot()

def parse_mapfile(mapfile, genome_fa, genome_gtf, gff, minimum_overlap, step):
    df_mapfile = pd.read_csv(mapfile, sep=r'\s+', header=0)
    df_mapfile['genome_fa'] = genome_fa
    df_mapfile['genome_gtf'] = genome_gtf
    df_mapfile['gff'] = gff
    df_mapfile['minimum_overlap'] = minimum_overlap
    
    if(step is None): # None run all steps
        step = 'intergenic,pr_en,bam,stat_df'
    df_mapfile['step'] = step

    bam_list = df_mapfile.bam
    spname_list = df_mapfile.spname
    genome_fa_list = df_mapfile.genome_fa
    genome_gtf_list = df_mapfile.genome_gtf
    gff_list = df_mapfile.gff
    minimum_overlap_list = df_mapfile.minimum_overlap
    step_list = df_mapfile.step

    return genome_fa_list, genome_gtf_list, gff_list, bam_list, spname_list, minimum_overlap_list, step_list


def run_single(genome_fa, genome_gtf, gff, bam, spname, minimum_overlap, step):
    runner = Intergenic_reads(genome_fa, genome_gtf, gff, bam, spname, minimum_overlap, step)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--genome_fa', help='mapfile', required=True)
    parser.add_argument('--genome_gtf', help='mapfile', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--minimum_overlap', default=0.5)
    parser.add_argument('--step')
    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    genome_fa_list, genome_gtf_list, gff_list, bam_list, spname_list, minimum_overlap_list, step_list = parse_mapfile(args.mapfile, args.genome_fa, args.genome_gtf, args.gff, args.minimum_overlap, args.step)
    #print(gff_list)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, genome_fa_list, genome_gtf_list, gff_list, bam_list, spname_list, minimum_overlap_list, step_list)