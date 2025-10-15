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

THREADS = 16
SEED = 100

class Diamond():
    def __init__(self, sample, n_reads):
        self.sample = sample
        self.n_reads = n_reads
        self.outdir = f'{self.sample}/diamond_outdir'
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def run(self):
        cmd1 = f'seqtk sample -s {SEED} {self.sample}/*star*/{self.sample}_Unmapped.out.mate1 {self.n_reads} | seqtk seq -a > {self.outdir}/{self.sample}_unmap_{self.n_reads}.fa'
        cmd2 = (
            f'/Public/Software/miniconda2/envs/old/bin/diamond blastx '
            f'--query {self.outdir}/{self.sample}_unmap_{self.n_reads}.fa '
            f'--out {self.outdir}/{self.sample}_blast.xml '
            f'--db /SGRNJ/Public/Database/nr/nr.dmnd '
            f'--outfmt 5 '
            f'--evalue 1e-5 '
            f'--max-target-seqs 1 '
            f'--threads {THREADS} '
        )
        cmd3 = (
            f'/SGRNJ/Public/Script/Pipeline/singlecell-rna-seq/scRNA_SEQ_V1/Module/python/get_xml_info '
            f'--xml {self.outdir}/{self.sample}_blast.xml '
            f'--outfile {self.outdir}/{self.sample}_contamination.txt '
        )
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header=0)
    sample_list = df_mapfile.loc[:,"sample"]
    n_reads_list = df_mapfile.loc[:,"n_reads"]
    return sample_list, n_reads_list

def run_single(sample, n_reads):
    runner = Diamond(sample, n_reads)
    runner.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    sample_list, n_reads_list = parse_mapfile(args.mapfile)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, sample_list, n_reads_list)