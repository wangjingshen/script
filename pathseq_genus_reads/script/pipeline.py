#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob

class Pathseq_genus_reads():
    def __init__(self, dir, sample, tax_id, tax_name, match_dir, outdir):
        self.sample = sample
        self.pathseq_bam = glob.glob(f"{dir}/{sample}/02.pathseq/{sample}_pathseq.bam")[0]
        self.downsample_bam = glob.glob(f"{dir}/{sample}/02.pathseq/{sample}_downsample.bam")[0]
        self.tax_id = tax_id
        self.tax_name = tax_name
        self.match_dir = match_dir
        self.outdir = outdir

        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def awk(self):
        cmd1 = (
            f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/samtools view {self.downsample_bam} | '
            '''awk -v FS="\t" -v OFS="\t" '{print $1,$17,$18}' >  '''
            f'{self.outdir}/downsample.tsv'
        )
        cmd2 = (
            f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/samtools view {self.pathseq_bam} | '
            '''awk -v FS="\t" -v OFS="\t" '{print $1,$16}' >  '''
            f'{self.outdir}/pathseq.tsv'
        )
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)

    def run_analysis(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script_dev/pathseq_genus_reads/script/analysis.R '
            f'--sample {self.sample} '
            f'--pathseq_tsv {self.outdir}/pathseq.tsv '
            f'--downsample_tsv {self.outdir}/downsample.tsv '
            f'--tax_id {self.tax_id} '
            f'--tax_name {self.tax_name} '
            f'--match_dir {self.match_dir} '
            f'--outdir {self.outdir} '
        )
        subprocess.check_call(cmd, shell=True)

    def rm_file(self):
        cmd = f'rm {self.outdir}/downsample.tsv {self.outdir}/pathseq.tsv'
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.awk()
        self.run_analysis()
        self.rm_file()

def parse_mapfile(mapfile, tax_id, tax_name):
    df_mapfile = pd.read_csv(mapfile, sep='\s+')
    df_mapfile['tax_id'] = tax_id
    df_mapfile['tax_name'] = tax_name

    dir_list = df_mapfile['dir']
    sample_list = df_mapfile['sample']
    tax_id_list = df_mapfile['tax_id']
    tax_name_list = df_mapfile['tax_name']
    match_dir_list = df_mapfile['match_dir']
    outdir_list = df_mapfile['outdir']
    return dir_list, sample_list, tax_id_list, tax_name_list, match_dir_list, outdir_list


def run_single(dir, sample, tax_id, tax_name, match_dir, outdir):
    runner = Pathseq_genus_reads(dir, sample, tax_id, tax_name, match_dir, outdir)
    runner.run()
    print(f"{sample} done.")


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--tax_id', help='tax_id', required=True)
    parser.add_argument('--tax_name', help='tax_name', required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    dir_list, sample_list, tax_id_list, tax_name_list, match_dir_list, outdir_list = parse_mapfile(args.mapfile, args.tax_id, args.tax_name)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, dir_list, sample_list, tax_id_list, tax_name_list, match_dir_list, outdir_list)